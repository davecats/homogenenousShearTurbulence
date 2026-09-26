# Next session: what is left on one node, or the second node

Copy the block at the end as the opening message of the next session.

State at the end of this one: `main` at the commit that adds this file,
pushed.  Both items of the previous handoff are done: the line solver is
two passes instead of five (half the bytes; the sweep was bandwidth-bound,
not latency-bound as the handoff said, see FINDINGS.md), with a register
cap of 168 for its kernel; and the alltoall overlaps the transforms
(per-field double-buffered transposes, NCCL on a second stream).  Per
step on the A100 (README table): 256^3 one GPU 0.099 -> 0.087, four
0.039 -> 0.034; 512^3 one GPU 0.819 -> 0.687, four 0.262 -> 0.227.  Read
`README.md`, then FINDINGS.md (the last two sections), then this file.

## Safety net, use it before and after every change

```bash
tests/run_tests.sh build-cpu 2       # 12 runs: transforms, solver kinds, Kelvin (5 decks),
tests/run_tests.sh build-gpu 2       #   pressure, Taylor-Green, forcing, conservation, Stokes
tests/regression.sh build-cpu 2      # 50 steps of three decks against tests/reference/*.fld at 1e-10
tests/regression.sh build-gpu 1
tests/regression.sh build-gpu 2      # the MPI transport on one GPU (two ranks)
```

The CPU build must be made in a shell that has *not* sourced
`env/istm.sh` (that puts NVHPC's mpifort first on PATH), and the GPU
tests must be run in one that *has* (the system `mpirun` cannot launch
the NVHPC binary; the symptom is "run failed" for every deck).  A change
to the Makefile's flags does not recompile the objects: `rm` the object
(or the build directory) first.  An NCCL build on the ISTM boxes
(`make GPU=1 NCCL=1 BUILD=build-nccl`) can only be exercised on istmcetus
(two RTX A6000; the second stream of the overlap is only used there and on
HoreKA):

```bash
ssh istmcetus 'cd ~/Codes/hst/homogenenousShearTurbulence && source env/istm.sh && tests/run_tests.sh build-nccl 2 && tests/regression.sh build-nccl 2'
```

On HoreKA (`~/hst` there is an rsync copy of the repository, not a
checkout; never rebuild it while a job on it is queued; A/B variants go
in `~/hst-exp`, `~/hst-exp2`, ... which are rsync copies too):

```bash
rsync -a --delete --exclude 'build-*' --exclude .git --exclude 'hst-*.out' ./ horeka:hst/
ssh horeka 'cd hst && source env/horeka.sh gpu && make GPU=1 GPU_ARCH=cc80 NCCL=1 -j8 && make GPU=1 GPU_ARCH=cc80 NCCL=1 test'
sbatch jobs/horeka_tests.slurm        # suite on 4 A100 (NCCL) + 4-GPU run against a CPU run
sbatch jobs/horeka_bench_all.slurm    # bench_64/256/512 on 1 and 4 A100 (the README table)
sbatch jobs/horeka_profile.slurm      # phase timer + nsys kernel summary on 1 A100
sbatch --ntasks-per-node=4 --gres=gpu:4 --export=ALL,NP=4 jobs/horeka_profile.slurm   # 4 GPUs, mpi and nccl
sbatch --export=ALL,A_ROOT=$HOME/hst,B_ROOT=$HOME/hst-exp jobs/horeka_ab.slurm       # A/B: both builds, same node, timer + nsys
sbatch --export=ALL,A_ROOT=$HOME/hst-exp,B_ROOT=$HOME/hst-exp2,C_ROOT=$HOME/hst-exp3 jobs/horeka_ab.slurm   # three-way
sbatch --export=ALL,HST_ROOT=$HOME/hst-exp,KERNEL=regex:buildrhs,SKIP=6,COUNT=2 jobs/horeka_ncu.slurm   # Nsight Compute counters
```

`dev_accelerated` runs one job per user at a time and queues at most
four; a two-minute job waits 5-40 minutes.  The A/B job is the way to
time a change: the same node, both builds back to back, and NCCL's
alltoall varies by 30% between runs, so never compare alltoall phases
across jobs.  **The hardware counters are open on the HoreKA compute
nodes** (`jobs/horeka_ncu.slurm`: DRAM bytes, achieved occupancy,
registers, stall reasons of a kernel family), not on the ISTM boxes
(`ERR_NVGPUCTRPERM`).  Use it before believing a byte count or a
"latency-bound" diagnosis: this session's handoff had the solver at 35%
of the bandwidth from an estimate, ncu said 82%.  The RTX 3060 is good
for kernel-level A/Bs with nsys (`mpirun -np N nsys profile -t cuda ...`,
then `nsys stats --report cuda_gpu_kern_sum`; `cuobjdump -res-usage` on
the object file gives registers and shared memory per kernel), but it
runs FP64 at 1/64 rate, so the solver and the FFTs are compute-bound
there and do not predict the A100 (the new solver kernel takes the same
time there as the three old ones).

## Where the time goes now (FINDINGS.md, last sections)

One A100, bench_256, kernel time: line solver 24% (1.86 ms per call for
1.59 GB, 55% of the bandwidth at 19% occupancy), cuFFT 28% (four
transforms), `buildrhs` 15%, tiled transpose 9%, `build_products` 8.5%,
`buildrhs_prepare` 5%, `assemble_vvdz` 4%.  Four A100: the timer phases
are "to physical: FFTs, transposes, CFL" 20% and "products: FFTs,
transposes, rhs" 49% (the transposes have no phase of their own any more,
only their exposed part counts), the implicit solves 15%, ghosts 11%;
the exposed alltoall is about 13% of the step.

## What is worth doing, in order

1. **`buildrhs` (15%, 2.45 ms per call at 256^3, about half the
   bandwidth).**  Its `VVdz` stencil reads are a plane apart between
   neighbouring threads (half of every 32-byte sector wasted); the other
   loop order was tried and lost (FINDINGS.md, session 3).  The fix is a
   stencil-and-transpose tile: read `VVdz` in its own layout (iz
   contiguous) into shared memory with the two ghost rows on each side,
   apply the three stencils along iy in the tile, write `rhs` in its
   layout (iy contiguous).  That is the CUDA Fortran pattern of
   `transpose_tiled` (DESIGN.md 7 allows CUDA Fortran in hst_mpi and
   hst_fft only; a kernel of hst_equations would need that rule
   extended, or the kernel placed in hst_mpi next to the tile kernel).
   `buildrhs_prepare` (5%) could then be folded into it (costed at 2% of
   the step, FINDINGS.md session 4).
2. **The line solver's remaining third (latency at 19% occupancy).**
   Fewer registers by construction: the six accumulators of the rows-0-
   and-1 recurrences could be a 2-vector recurrence with four; the ten
   border coefficients are computed per line after the sweep and could be
   hoisted; `coef` evaluates fifteen `der` products per row that a
   per-kind table `cf(iy, j, 0:2)` (coef = cf0 + k2 cf1 + k2^2 cf2,
   built once per call) would replace by three loads and two FMAs.  Or
   two rows of loads in flight in the backward sweep.  Measure each with
   the A/B job and `jobs/horeka_ncu.slurm` (occupancy, stall reasons);
   the sweep's upper bound is about 1.1 ms at the bandwidth, so at most
   another 7% of the 1-GPU step.
3. **Deeper overlap** (the alltoalls of the first product group behind
   the products and `buildrhs` of the second; six products in memory at
   once): at most the exposed 13% of the 4-GPU step, more likely half of
   it.  Only after 1 and 2.
4. **y decomposition (WP6)** when going beyond one node: DESIGN.md 7 (i).
   Measure `--nodes=2` with the present code first (8 ranks over
   InfiniBand: NCCL handles it, the alltoall grows by the inter-node
   share).
5. **Memory per rank** at 512^3: 19.5 GB of transform buffers on one
   rank, 4.9 GB on four; the line-solver workspace with all columns is
   5 x 16 B x lines x ny (`line_chunk` bounds it); the transpose buffers
   are now two pairs of one field each (two thirds of the old three-field
   pair).

Things learned this session that the next one should not relearn
(FINDINGS.md has the numbers): the cyclic solve needs no third pass
(rows 0 and 1 of the back substitution come from a forward recurrence);
a device routine cannot take an assumed-shape dummy for the field
(illegal address: pass explicit shape, and make the `line_solve` dummies
`contiguous`); a scalar `u1` clashes with an array `U1` (case); a
target-specific `FFLAGS +=` in the Makefile works for a per-file flag;
`-gpu=maxregcount:N` at 128 spills and loses; the cuFFT plans per field
cost nothing; NCCL on a second stream needs `cudaStreamNonBlocking` and
two events per buffer pair, and `MPI_Ialltoall` works on device buffers
under `use_device_addr` if the compute stream is synchronised first.

## Opening message for the next session

```
Repository ~/Codes/hst/homogenenousShearTurbulence (also ~/hst on HoreKA as
an rsync copy), a GPU/CPU DNS for homogeneous shear turbulence; read
README.md, FINDINGS.md (last two sections) and NEXT_SESSION.md.  Task:
NEXT_SESSION.md item 1, buildrhs as a stencil-and-transpose tile (measure
first with jobs/horeka_ncu.slurm, then A/B on the A100 with
jobs/horeka_ab.slurm), then item 2 (the line solver's occupancy) if the
ncu profile still shows it latency-bound; tests/run_tests.sh and
tests/regression.sh green after every step (CPU, GPU, and the NCCL build
on istmcetus), HoreKA jobs for the numbers.  Do not modify
~/Codes/hst/channel or ~/Codes/hst/hst-main.  Commit each step; push at
the end.
```
