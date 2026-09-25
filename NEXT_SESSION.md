# Next session: the line solver, then the overlap

Copy the block at the end as the opening message of the next session.

State at the end of this one: `main` at the commit that adds this file,
pushed.  Of the previous handoff, items 1-3 are done as far as they pay:
the four transposing kernels are one tiled CUDA Fortran kernel through
shared memory (the two pack kernels were block copies already),
`build_products` reads each field once, and the `buildrhs` fusions were
costed at 2% of the step and left alone.  Item 4 (overlap of the
alltoall) was conditional on the alltoall being the largest item after
item 1, which it is not.  Per step on the A100 (FINDINGS.md, last
section): 256^3 one GPU 0.113 -> 0.098, four 0.045 -> 0.041; 512^3 one
GPU 0.957 -> 0.820, four 0.306 -> 0.268.  Read `README.md`, then
FINDINGS.md, then this file.

## Safety net, use it before and after every change

```bash
tests/run_tests.sh build-cpu 2       # 12 runs: transforms, solver kinds, Kelvin (5 decks),
tests/run_tests.sh build-gpu 2       #   pressure, Taylor-Green, forcing, conservation, Stokes
tests/regression.sh build-cpu 2      # 50 steps of three decks against tests/reference/*.fld at 1e-10
tests/regression.sh build-gpu 1
```

The CPU build must be made in a shell that has *not* sourced
`env/istm.sh` (that puts NVHPC's mpifort first on PATH), and the GPU
tests must be run in one that *has* (the system `mpirun` cannot launch
the NVHPC binary; the symptom is "run failed" for every deck).  An NCCL
build on the ISTM boxes (`make GPU=1 NCCL=1 BUILD=build-nccl`) can only
be exercised on istmcetus (two RTX A6000):

```bash
ssh istmcetus 'cd ~/Codes/hst/homogenenousShearTurbulence && source env/istm.sh && tests/run_tests.sh build-nccl 2 && tests/regression.sh build-nccl 2'
```

On HoreKA (`~/hst` there is an rsync copy of the repository, not a
checkout; never rebuild it while a job on it is queued; A/B variants go
in `~/hst-exp`, `~/hst-exp2`, ...):

```bash
rsync -a --delete --exclude 'build-*' --exclude .git --exclude 'hst-*.out' ./ horeka:hst/
ssh horeka 'cd hst && source env/horeka.sh gpu && make GPU=1 GPU_ARCH=cc80 NCCL=1 -j8 && make GPU=1 GPU_ARCH=cc80 NCCL=1 test'
sbatch jobs/horeka_tests.slurm        # suite on 4 A100 (NCCL) + 4-GPU run against a CPU run
sbatch jobs/horeka_bench_all.slurm    # bench_64/256/512 on 1 and 4 A100 (the README table)
sbatch jobs/horeka_profile.slurm      # phase timer + nsys kernel summary on 1 A100
sbatch --ntasks-per-node=4 --gres=gpu:4 --export=ALL,NP=4 jobs/horeka_profile.slurm   # 4 GPUs, mpi and nccl
sbatch --export=ALL,A_ROOT=$HOME/hst,B_ROOT=$HOME/hst-exp jobs/horeka_ab.slurm       # A/B: both builds, same node, timer + nsys
```

`dev_accelerated` runs one job per user at a time and queues at most
four; a two-minute job waits 5-40 minutes.  The A/B job is the way to
time a change: the same node, both builds back to back, and NCCL's
alltoall varies by 30% between runs, so never compare alltoall phases
across jobs.  The RTX 3060 is good for kernel-level A/Bs with nsys (`mpirun -np N
nsys profile -t cuda ...`, then `nsys stats --report cuda_gpu_kern_sum`;
`cuobjdump -res-usage` on the object file gives registers and shared
memory per kernel), and its bandwidth-bound kernels predicted the A100
well this time (transpose 1.6x there, 2.1x on the A100; products 1.4x
and 2.0x); the FP64-bound ones (solver, FFTs) do not.

## Where the time goes now (FINDINGS.md, last section)

One A100, bench_256, kernel time: line solver 27% (the sweep 2.05 ms per
call), cuFFT 25%, `buildrhs` 14%, tiled transpose 8%, `build_products`
8%, `buildrhs_prepare` 4.5%, `assemble_vvdz` 3.5%.  Four A100: the
products with their transforms and `buildrhs` 25-30%, the alltoall
18-23%, the implicit solves 13-17%, pack/unpack 12-13%.

## What is worth doing, in order

1. **The line solver is latency-bound, not bandwidth-bound.**  Its
   workspace is already interleaved (line index first, `U(il, iy, j)`),
   so the sweep's own reads and writes are coalesced; yet the sweep
   kernel (`line_solve` F1L126) moves 1.1 GB in 2.05 ms at 256^3, 35% of
   the A100, at 1% of its FP64 peak.  The facts to start from
   (`cuobjdump -res-usage build-gpu/hst_linsolve.o`): 136 registers per
   thread, so at most 3 blocks of 128 threads per SM (19% occupancy),
   and only (2nz+1)(nx+1) = 32640 threads at 256^3, i.e. 300 per SM on
   108 SMs, each walking a dependent chain of ny = 256 rows with one
   complex division and about fifteen dependent complex multiply-adds
   per row.  Nothing hides that latency.  Ways to more threads or
   shorter chains, to be measured with the A/B job: (a) both implicit
   systems (eta and d2v) in one launch, `kind` chosen per line, twice
   the threads for the same work; (b) fewer registers, e.g. the two
   border columns Y1, Y2 computed in a second, cheaper sweep (their
   right-hand sides are two stencil entries of rows 0 and 1, so they
   could also be tabulated once per substep since they depend on the
   mode only through k2 and the phase); (c) two threads per line, each
   half of the y range, joined by a 2x2 Schur step (the WP6 machinery of
   DESIGN.md 7 (i), which then serves both purposes).  `buildrhs` (14%)
   has a real access problem: its `VVdz` stencil reads are a plane apart
   between neighbouring threads (half of every sector wasted), and the
   other loop order was tried and lost; a fix means reading `VVdz` in
   its own layout and transposing on the way (a stencil-and-transpose
   tile), which is the CUDA Fortran pattern of `transpose_tiled`.
2. **Overlap of the alltoall with the transforms** on four GPUs (18-23%
   of the step): three alltoalls of one field each with double buffers
   (the channel's `CHANNEL_OVERLAPPING`), on a second stream for the
   transforms, since NCCL already runs asynchronously on the OpenMP
   stream.  Splits the three-field batches again.  Worth at most that
   share; measure with `--nodes=2` first if WP6 is the real goal.
3. **`buildrhs_prepare` into `buildrhs`**: costed at 2% (FINDINGS.md),
   only if 1(b) is done anyway.
4. **y decomposition (WP6)** when going beyond one node: DESIGN.md 7 (i).
5. **Memory per rank** at 512^3: 19.5 GB of transform buffers on one
   rank, 4.9 GB on four; the line-solver workspace with all columns is
   6 x 16 B x lines x ny (`line_chunk` bounds it).

Things learned this session that the next one should not relearn
(FINDINGS.md has the numbers): nvfortran 25.9's OpenMP cannot be made to
put a team-private tile in shared memory with good inner loops (`teams
distribute` + `parallel do` does the memory but not the loops, `teams
loop` the reverse, `allocate(omp_pteam_mem_alloc)` is unrecognised), so
a kernel that needs shared memory is CUDA Fortran, launched on the
OpenMP stream from `use_device_addr` pointers as `transpose_tiled` does;
and a Fortran array may not be called `tile` next to a parameter `TILE`.

## Opening message for the next session

```
Repository ~/Codes/hst/homogenenousShearTurbulence (also ~/hst on HoreKA as
an rsync copy), a GPU/CPU DNS for homogeneous shear turbulence; read
README.md, FINDINGS.md and NEXT_SESSION.md.  Task: NEXT_SESSION.md item 1,
the latency-bound line solver (start from the facts listed there, try
(a) both systems in one launch first, measure each idea on the A100 with
jobs/horeka_ab.slurm), then
item 2 (overlap of the alltoall on four GPUs) if it still shows as 18% or
more of the 4-GPU step; tests/run_tests.sh and tests/regression.sh green
after every step, HoreKA jobs for the numbers.  Do not modify
~/Codes/hst/channel or ~/Codes/hst/hst-main.  Commit each step; push at
the end.
```
