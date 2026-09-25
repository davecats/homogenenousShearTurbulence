# Next session: the kernels around the transposes

Copy the block at the end as the opening message of the next session.

State at the end of this one: `main` at the commit that adds this file,
pushed.  The multi-GPU path of the previous handoff is done as far as it
pays: the timer shows the transposes as their own phases, the alltoall
goes through NCCL (`make GPU=1 NCCL=1`, deck `transport`), four A100 are
3.3x one at 512^3 (81% parallel efficiency), and of Part 1 the reciprocal
pivots are in and the `buildrhs` loop order was tried and rejected
(FINDINGS.md, last section).  Not done: overlap of the alltoall with the
transforms (measured at 16% of the step, the ceiling of what it can hide)
and WP6.  Read `README.md`, then FINDINGS.md, then this file.

## Safety net, use it before and after every change

```bash
tests/run_tests.sh build-cpu 2       # 12 runs: transforms, solver kinds, Kelvin (5 decks),
tests/run_tests.sh build-gpu 2       #   pressure, Taylor-Green, forcing, conservation, Stokes
tests/regression.sh build-cpu 2      # 50 steps of three decks against tests/reference/*.fld at 1e-10
tests/regression.sh build-gpu 1
```

The CPU build must be made in a shell that has *not* sourced
`env/istm.sh` (that puts NVHPC's mpifort first on PATH).  An NCCL build
on the ISTM boxes (`make GPU=1 NCCL=1 BUILD=build-nccl`) can only be
exercised on istmcetus (two RTX A6000; NCCL refuses two ranks on one GPU,
`transport = 'auto'` falls back to MPI there):

```bash
ssh istmcetus 'cd ~/Codes/hst/homogenenousShearTurbulence && source env/istm.sh && tests/run_tests.sh build-nccl 2 && tests/regression.sh build-nccl 2'
```

On HoreKA (`~/hst` there is an rsync copy of the repository, not a
checkout; never rebuild it while a job on it is queued, put A/B variants
in `~/hst-exp`):

```bash
rsync -a --delete --exclude 'build-*' --exclude .git --exclude 'hst-*.out' ./ horeka:hst/
ssh horeka 'cd hst && source env/horeka.sh gpu && make GPU=1 GPU_ARCH=cc80 NCCL=1 -j8 && make GPU=1 GPU_ARCH=cc80 NCCL=1 test'
sbatch jobs/horeka_tests.slurm        # suite on 4 A100 (NCCL) + 4-GPU run against a CPU run
sbatch jobs/horeka_bench_all.slurm    # bench_64/256/512 on 1 and 4 A100 (the README table)
sbatch jobs/horeka_profile.slurm      # phase timer + nsys kernel summary on 1 A100
sbatch --ntasks-per-node=4 --gres=gpu:4 --export=ALL,NP=4 jobs/horeka_profile.slurm   # 4 GPUs, mpi and nccl
```

`dev_accelerated` runs one job per user at a time; a two-minute job waits
5-20 minutes behind the user's other dev jobs.  Timing on the RTX 3060 is
only good for A/B comparisons made back to back and only for FP64-bound
kernels (the solver); memory-access changes must be timed on the A100
(the `buildrhs` loop order looked neutral on the RTX 3060 and lost 80%
on the A100).

## Where the time goes now (FINDINGS.md, last section, has the tables)

Four A100 with NCCL, share of the step at 512^3: products + FFTs +
`buildrhs` 31%, pack/unpack 18%, alltoall 16%, implicit solves 13%,
recovery 7%, transform to physical 10%.  One A100 at 256^3, kernel time:
line solver 26%, `build_products` 13.5%, `buildrhs` 12.3%, the local
repack between pencil layouts 14.6%, cuFFT 22%.

## What is worth doing, in order

1. **Tiled transposes.**  `pack_*`, `unpack_*` and `repack_*_local` in
   `src/hst_mpi.f90` are plain index-swapping loops: one side coalesced,
   the other a stride of nzB or nxB complex numbers.  They run at 40-50%
   of the A100's bandwidth (2.0 ms for 1.2 GB at 256^3 on one GPU, 3 ms
   for 2.4 GB at 512^3 on four).  A 32x32 tile through team-local memory
   (`teams distribute` over tiles, `parallel do` inside, a team-private
   array; check with nsys that nvfortran puts it in shared memory) should
   give 2-3x on 15-18% of the step.  Measure on the A100 with
   `~/hst-exp` against `~/hst` (the A/B job of this session is described
   in FINDINGS.md; it was a one-off script, not in `jobs/`).
2. **The local repack on one GPU is pure overhead** of the pencil
   abstraction (14.6% at 256^3).  Either item 1, or transform the z lines
   in the x-pencil layout directly: cuFFT takes `istride = nxd+1`,
   `idist = 1` for a batch over ix, but the batch over y rows needs a
   second distance cuFFT's 1-D `many` plan does not have, so it is one
   call per y row and field (780 calls at 256^3), or a 2-D
   `cufftPlanMany` trick.  Item 1 first.
3. **`build_products` + `buildrhs` + `buildrhs_prepare`** are 30% of one
   A100.  `build_products` is a pure streaming kernel (read three real
   fields, write three) at 2.7 ms for 1.5 GB read + 0.75 GB written at
   256^3, i.e. 80% of bandwidth: nothing to gain in the kernel, only by
   not materialising `products` (write the products straight into the
   real view of `VVdp`?  it already does; the three-field batching
   requires the copy).  `buildrhs`: 15 stencil reads and 8 field accesses
   per element, L1-served, 2.4 ms; `buildrhs_prepare` builds `helm` and
   `biharm` from `der` and `k2` per element, 1.5 ms.  Fusing `prepare`
   into `buildrhs` saves one read+write of `rhs`/`oldrhs` (4 x 0.5 GB at
   256^3).  Small.
4. **Overlap** the alltoall with the transforms of the next batch: at most
   16% of the 4-GPU step, and it means splitting the three-field batches
   again (three alltoalls of one field each with double buffers, the
   channel's `CHANNEL_OVERLAPPING`).  With NCCL on the OpenMP stream the
   transfer is already asynchronous with respect to the host; the
   overlap would need a second stream for the transforms.  Only after 1.
5. **y decomposition (WP6)** when going beyond one node: DESIGN.md 7 (i)
   sketches the distributed cyclic Schur solve; the `ny0:nyN` slab
   indexing and the single `line_solve` entry point are in place.  The
   alltoall over two nodes goes through InfiniBand, where NCCL and MPI
   should be closer; measure with `--nodes=2` first.
6. **Memory per rank** at 512^3: 19.5 GB of transform buffers on one rank,
   4.9 GB on four; the line-solver workspace with all columns is 6 x 16 B
   x lines x ny (`line_chunk` bounds it); `nvidia-smi` in the bench job
   prints the total.

## Opening message for the next session

```
Repository ~/Codes/hst/homogenenousShearTurbulence (also ~/hst on HoreKA as
an rsync copy), a GPU/CPU DNS for homogeneous shear turbulence; read
README.md, FINDINGS.md and NEXT_SESSION.md.  Task: the kernels around the
transposes, NEXT_SESSION.md items 1-3 in that order (tiled transposes for
pack/unpack and the local repack, measured on the A100 with an A/B job;
then the product/buildrhs kernels), with tests/run_tests.sh and
tests/regression.sh green after every step and the HoreKA jobs for the
numbers; overlap (item 4) only if the alltoall is still the largest
single item after item 1.  Do not modify ~/Codes/hst/channel or
~/Codes/hst/hst-main.  Commit each step; push at the end.
```
