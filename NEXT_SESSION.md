# Next session: the multi-GPU path

Copy the block at the end as the opening message of the next session.

State at the end of this one: `main` at the commit that adds this file,
pushed.  Everything in DESIGN.md is implemented except WP6 (y
decomposition) and the NCCL transport; the cleanup and the single-GPU
performance pass of the previous handoff are done (FINDINGS.md, last
section).  Read `README.md`, then FINDINGS.md, then this file.

## Safety net, use it before and after every change

```bash
tests/run_tests.sh build-cpu 2       # 12 runs: transforms, solver kinds, Kelvin (5 decks),
tests/run_tests.sh build-gpu 2       #   pressure, Taylor-Green, forcing, conservation, Stokes
tests/regression.sh build-cpu 2      # 50 steps of three decks against tests/reference/*.fld at 1e-10
tests/regression.sh build-gpu 1
```

The CPU build must be made in a shell that has *not* sourced
`env/istm.sh` (that puts NVHPC's mpifort first on PATH).  On HoreKA
(`~/hst` there is an rsync copy of the repository, not a checkout):

```bash
rsync -a --delete --exclude 'build-*' --exclude .git ./ horeka:hst/
ssh horeka 'cd hst && source env/horeka.sh gpu && make GPU=1 GPU_ARCH=cc80 -j8 && make GPU=1 GPU_ARCH=cc80 test'
sbatch jobs/horeka_tests.slurm        # suite on 4 A100 + 4-GPU run against a CPU run
sbatch jobs/horeka_bench_all.slurm    # bench_64/256/512 on 1 and 4 A100
sbatch --ntasks-per-node=4 --gres=gpu:4 --export=ALL,NP=4 jobs/horeka_profile.slurm   # phase timer on 4 GPUs
```

`dev_accelerated` runs one job per user at a time.  Timing on the RTX 3060
is only good for A/B comparisons made back to back: the card runs FP64 at
1/64 rate and is sometimes shared with other jobs.

## Where the time goes now (see the tables in FINDINGS.md and README)

On one A100 the line solver is still the largest single item, then the
product accumulation (`buildrhs`) and the transposes.  On four A100 the
alltoall (pack, MPI, unpack) is what limits the scaling; the phase timer
of `jobs/horeka_profile.slurm` with `NP=4` shows the split.

## Part 1: what is worth doing on one GPU

1. **`buildrhs` and `compute_pressure` read `VVdz` strided.**  The loop
   has `iy` innermost so that `rhs(iy, iz, ix)` is written contiguously,
   but the five-point stencil reads `VVdz(izd(iz)+1, ix-nx0+1, iy+j, p)`
   at a stride of a whole plane.  Either order the loop `(iy, ix, iz)`
   with `iz` innermost (contiguous reads of `VVdz`, strided writes of
   `rhs`/`oldrhs`, 4 accesses against 15 reads) or let the transpose
   deliver the z-pencil with `iy` fastest.  Measure with nsys first
   (`jobs/horeka_profile.slurm`): `buildrhs` was 12% of the A100 kernel
   time before batching.
2. **The line solver's complex divisions.**  Four per row (two
   multipliers, and the back substitution of each of the three right-hand
   sides divides by the pivot).  Storing the reciprocal of the pivot
   instead of the pivot turns three of them into multiplications.  This
   changes the rounding: run the regression and expect 1e-14, then
   `tests/regression.sh ... --update` is *not* needed (1e-10 tolerance).
3. **`buildrhs_prepare`** builds `helm` and `biharm` from `der` and `k2`
   per element; a small thing, 2% on the A100.

## Part 2: the multi-GPU path

1. **Measure**: `NP=4` profile of bench_256 and bench_512, and the same
   with `nsys` (`-t cuda,mpi` shows the alltoall waits).  Parallel
   efficiency of the original code was 51% at 512^3 on four A100.
2. **NCCL transport** for the alltoall on NVLink nodes.  The collective
   lives in `hst_mpi::alltoall` only; the channel code's
   `channel_comm_alltoall_complex` and `channel_nccl_p2p.c` are the
   drop-in (NCCL under `$NVHPC_ROOT/comm_libs/13.0/nccl`).  Worth it only
   if the timer says the MPI alltoall is a large fraction at 4 GPUs.
3. **Overlap**: with three fields per transpose there are three alltoalls
   per substep.  The channel code double-buffers so that the alltoall of
   one batch overlaps the transforms of the next (`CHANNEL_OVERLAPPING`);
   here that would mean splitting each three-field batch into three
   again, which is the launch count we just removed.  Only after 1 and 2.
4. **y decomposition (WP6)** when scaling beyond one node: DESIGN.md 7 (i)
   sketches the distributed cyclic Schur solve; the `ny0:nyN` slab
   indexing and the single `line_solve` entry point are in place.
5. **Memory per rank** at 512^3: three-field buffers cost 3x the transform
   memory (19.5 GB of buffers on one rank at 512^3, 4.9 GB on four); the
   line-solver workspace with all columns is 6 x 16 B x lines x ny.
   `nvidia-smi` inside the job prints the total.

## Opening message for the next session

```
Repository ~/Codes/hst/homogenenousShearTurbulence (also ~/hst on HoreKA as
an rsync copy), a GPU/CPU DNS for homogeneous shear turbulence; read
README.md, FINDINGS.md and NEXT_SESSION.md.  Task: the multi-GPU path of
NEXT_SESSION.md Part 2 (measure on 4 A100 first, then NCCL transport if
the alltoall dominates, then overlap only if measured to pay), with
tests/run_tests.sh and tests/regression.sh green after every step and
the HoreKA jobs for the numbers; Part 1 items if time remains.  Do not
modify ~/Codes/hst/channel or ~/Codes/hst/hst-main.  Commit each step;
push at the end.
```
