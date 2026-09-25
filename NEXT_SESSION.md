# Next session: simplification and cleanup, then performance

Copy the block at the end as the opening message of the next session.

State at the end of this one: `main` at the commit that adds this file,
pushed.  Everything in PLAN.md sections 1..10 is implemented except WP6
(y decomposition, NCCL).  Read `README.md`, then PLAN.md sections 8
(findings) and 9 (status), then this file.

## Safety net, use it before and after every change

```bash
tests/run_tests.sh build-cpu 2       # 12 runs: transforms, solver, Kelvin (5 decks),
tests/run_tests.sh build-gpu 2       #   pressure, Taylor-Green, forcing, conservation, Stokes
tests/regression.sh build-cpu 2      # 50 steps of three decks against tests/reference/*.fld at 1e-10
tests/regression.sh build-gpu 1
```

The regression references were produced by the code as it is now, on the
CPU; the GPU reproduces them to 1e-14.  `tests/regression.sh ... --update`
rewrites them and must only be run after a change that is *meant* to alter
the numbers (a different summation order does not count: 1e-10 leaves room
for that).  On HoreKA, `sbatch jobs/horeka_tests.slurm` runs the suite on
4 A100 and compares a 4-GPU run with a CPU run.

## Part 1: simplification and cleanup

Candidates, in the order I would take them.  None changes results;
`tests/regression.sh` proves it after each.

1. **One assembly kernel in `hst_linsolve`.**  `solve_component`, `apply_dy`
   and `unweight_d0` each spell the same chunk loop, the same wrap-phase
   assembly and the same gather/scatter, differing only in the row
   coefficients and the right-hand side.  A single `assemble_and_solve(kind,
   lambda, src, dst)` with a `select case (kind)` inside the kernel removes
   two of the three copies (~120 lines).
2. **The chunk loop itself.**  `line_chunk = 16` x columns per batch bounds
   the workspace `A, X, Y1, Y2`.  Measure whether a single batch of all
   `nxB` columns is affordable at 512^3 per rank (A is 5 x 16 B x ny x
   nlines: 10 GB at nx = 255, ny = 512, nz = 255 on one rank, 2.5 GB on
   four); if it is, the loop can go, or become a namelist parameter.
3. **`memrhs` as scratch.**  It is used as the assembly target in
   `buildrhs_prepare` and as scratch by `shift_unweighted`, `outstats`
   (component 1) and `write_pressure` (component 2).  Rename it `scratch`
   and document the four uses in `hst_params`, or split.
4. **Test boilerplate.**  Eight test programs repeat MPI init, device
   selection, deck reading and set-up.  A `tests/test_common.f90` with
   `test_start(deck)` / `test_finish(passed)` halves them.
5. **Statistics module.**  `outstats` now carries the mean flow, the seven
   second moments, the dissipation and the Stokes region sums in one
   kernel with a 13-entry reduction.  Consider one small kernel per group
   or keep it but name the slots (an `enum`-like parameter list instead of
   `glob(11)`).
6. **`hst_io`**: the two repack loops (write and read) and `field_write`
   share the ghost-row and component permutation; one helper each way.
7. **PLAN.md** has grown into a design document plus a lab notebook.  Split:
   `DESIGN.md` (sections 1..7, 10) and `FINDINGS.md` (8), keep a short
   status table in README.
8. Small things: the `no default(none)` workaround comment in
   `shear_shift` (check whether nvfortran 25.9 still needs it after the
   rewrite); `KIND_*` constants could become the `select case` above;
   `ystretch` and the Stokes parameters print only when active, the S2
   ones too, fine; `has_average` is used for the mean mode and for the
   pressure gauge, keep.

What not to touch: the file formats (CPL layout), the namelist names, the
`ny0:nyN` slab indexing (it is what makes WP6 an addition), the two
nvfortran rules in `src/physics/README.md` of the channel code
(no `use` in modules with `declare target` variables; no device procedure
called across modules).

## Part 2: performance

### Baselines (seconds per full step, three substeps, default settings)

| deck | grid | machine | s/step |
| --- | --- | --- | --- |
| examples/bench_64.in | 64 x 128 x 64 (default box) | RTX 3060, 1 rank | 0.178 |
| examples/bench_64.in | | istmio2 CPU, 1 / 4 ranks | 1.85 / 1.03 |
| examples/bench_64.in | | 1 x A100 | 0.055 |
| examples/bench_256.in | 256^3 | RTX 3060, 1 rank | 1.33 (2.05 with exact_shift) |
| examples/bench_256.in | | 1 x A100 / 4 x A100 | 0.33 / 0.145 |
| examples/bench_512.in | 512^3 | 4 x A100 | 0.99 |

The small grid is 32x fewer points than 256^3 but only 7x faster on the
RTX 3060 and 6x on the A100: it is bound by kernel launches and
synchronisations, not by arithmetic.  The large grids are bound by the transposes and the FFTs.

### First: measure, do not guess

Add a light per-phase timer before changing anything: `MPI_Wtime` around
the six phases of a substep (transform to physical, buildrhs_prepare,
products+buildrhs, shear_shift, the three solves in linsolve, ghosts and
recover), with a device synchronisation before each reading in the GPU
build, accumulated over the run and printed at the end when `timing =
.true.` in `&time_control`.  Twenty lines in `hst_equations`; no external
profiler needed on HoreKA.  Then `nsys profile` on istmio2 for the kernel
count and the gaps between kernels (NVHPC ships it under
`/opt/Nvidia/nvhpc/Linux_x86_64/25.9/profilers/Nsight_Systems`).

### Small problems (launch-bound)

Per substep today: about 18 kernels + 3 cuFFT calls for the forward
transforms, 6 x (product kernel + cuFFT + transpose kernels + cuFFT +
accumulate) for the products, 2 for `buildrhs_prepare`, 1 (or 4 x 6 with
`exact_shift`) for the shift, and for each of the three line solves
`(nx+1)/16` chunks x 3 kernels, plus 3 ghost fills and the u/w recovery.
Roughly 100 launches and 20 explicit `cudaDeviceSynchronize` calls per
substep.  In order of expected gain:

1. **Batch the six products.**  Form all six in one kernel into a
   `products(:, :, :, 6)` buffer, one cuFFT with batch x 6, one transpose
   (one alltoall with a 6x larger message instead of six), one FFT, one
   accumulate kernel that adds all six contributions.  Six launches and one
   collective instead of ~30 and six.  Memory: five more real buffers of
   the padded physical size (about 2 fields each).  This also helps the
   large grids (fewer, larger messages).
2. **Batch the three forward transforms** the same way (one buffer of
   three components, one cuFFT plan with batch x 3, one transpose).
3. **Remove the explicit `cudaDeviceSynchronize` around every cuFFT call**
   (inherited from the channel code).  They serialise the device.  The
   question to settle first is whether nvfortran's OpenMP offload and
   cuFFT run on the same stream; if not, `cufftSetStream` to the OpenMP
   stream (nvfortran offers `ompx_get_cuda_stream` in recent versions, to
   be checked) makes the syncs unnecessary.
4. **Fewer, larger line batches**: `line_chunk` = all columns when memory
   allows (item 2 of Part 1), so each solve is three launches, not
   3 x (nx+1)/16.
5. **Assemble inside the solve kernel**: each line thread builds its own
   five diagonals from `der`, `k2` and the wrap phase while it factors,
   instead of a separate assembly kernel writing `A` to global memory and
   the solve reading it back.  Halves the line-solve memory traffic and
   removes one launch per solve.  This also settles item 1 of Part 1 in a
   different way, so decide the two together.

### Large problems (transpose- and FFT-bound)

1. Read the timer: at 512^3 on 4 A100 the alltoall (pack, MPI, unpack)
   and the four cuFFT calls per component are expected to dominate.
2. **Overlap**: the channel code double-buffers so the alltoall of one
   component overlaps the transforms of the next (`CHANNEL_OVERLAPPING`).
   After batching (items 1 and 2 above) there is one big alltoall per
   direction per substep instead of nine, so overlap buys less; measure
   before re-adding the complexity.
3. **NCCL transport** for the alltoall on NVLink nodes: the alltoall lives
   in one routine (`hst_mpi::alltoall`), so this is a drop-in
   (channel's `channel_comm_alltoall_complex` and `channel_nccl_p2p.c`).
   Only worth it if the timer says the MPI alltoall is a large fraction at
   4 GPUs.
4. **y decomposition (WP6)** when scaling beyond one node: the all-to-all
   over all ranks and the `nproc | nx+1` constraint are the limits.  PLAN.md
   section 7 (i) sketches the distributed cyclic Schur solve.
5. Memory per rank at 512^3: check `nvidia-smi` in the job; the products
   batching adds buffers, the line workspace scales with `line_chunk`.

### CPU

The CPU build is `-O2` gfortran without OpenMP threads: one MPI rank per
core, the alltoall over all of them.  Two cheap things: `-O3
-march=native` (check the regression still holds at 1e-10; it should), and
`FFTW_MEASURE` instead of `FFTW_PATIENT` for shorter start-up.  Threads
(`-fopenmp` on the CPU build, which would turn the `target` directives into
host parallel loops) are a larger change; measure ranks-per-node first.

### Targets to aim for

- bench_64: from 0.18 to below 0.05 s/step on the RTX 3060, from 0.055 to
  about 0.02 on the A100.
- bench_256 on one A100: from 0.33 to about 0.2 s/step.
- bench_512 on four A100: parallel efficiency against one A100 (measure
  the 1-GPU number first; it fits in 40 GB).

## Opening message for the next session

```
Repository ~/Codes/hst/homogenenousShearTurbulence (also ~/hst on HoreKA), a
GPU/CPU DNS for homogeneous shear turbulence; read README.md, PLAN.md
sections 8-9 and NEXT_SESSION.md.  Two tasks, in this order:
(1) simplification and cleanup of the code following NEXT_SESSION.md Part 1,
    with tests/run_tests.sh and tests/regression.sh green after every step;
(2) performance on small and large problems, locally and on HoreKA,
    following Part 2: add the per-phase timer first, then the batched
    products and transforms, then the launch and synchronisation reductions;
    report the baselines table before and after.
Do not modify ~/Codes/hst/channel or ~/Codes/hst/hst-main.  Commit each
step; push at the end.
```
