# hst

Direct numerical simulation of homogeneous shear turbulence on CPUs and
NVIDIA GPUs, in about 2000 lines of plain Fortran.  Derived from the
[`channel`](https://github.com/davecats/channel) DNS code and from the CPL
HST code; the design and the reasoning behind every simplification are in
[DESIGN.md](DESIGN.md), the findings made on the way in
[FINDINGS.md](FINDINGS.md).

## What it solves

Incompressible Navier-Stokes fluctuations about the mean flow `U = S*y`,
optionally with an unsteady spanwise component `W = S2(t)*y`
(`S2 = A sin(2 pi (t - t0)/T)`, the `S2data.cpl` law) or a Stokes layer
(an oscillating spanwise profile at mid-box driven by its body force, the
`SLdata.cpl` model, on a grid clustered at mid-box with `ystretch`), in a
box that is periodic in `x` (streamwise) and `z` (spanwise) and
shear-periodic in `y`: Fourier in `x` and `z` with 3/2 dealiasing, compact
(sixth-order, five-point) finite differences in `y`, velocity-vorticity
(`v`, `eta`) formulation, RK3 (Rai-Moin) for the nonlinear terms,
Crank-Nicolson for the viscous ones, and the mean-shear advection
integrated exactly by a phase shift (Sekimoto, Dong & Jimenez).  The
phase shift is applied to the D0-weighted unknowns as in the CPL code, or,
with `exact_shift = .true.` in the deck, to the unweighted quantities,
which removes a second-order error at the cost of four extra line solves
per substep (FINDINGS.md).  The pressure is computed online at snapshot
times.

**Axis convention:** the shear direction is `y`, as in the channel code.
The CPL code calls that direction `z` (and its `w` is our `v`).

Stored modes are `0..nx` in `x` and `-nz..nz` in `z`; `ny` uniform points
span `ly` in `y`.  The number of MPI ranks must divide both `nx+1` and
`nzd` (`3*nz` rounded up to a power of two times at most one factor 3).
Each rank owns all of `y` (one GPU per rank, x-z pencils).

## Build

```bash
source env/istm.sh            # istmio2 / istmcetus / istmcorax
make                          # CPU: gfortran + MPI + FFTW  -> build-cpu/hst
make GPU=1                    # GPU: nvfortran + cuFFT       -> build-gpu/hst
make test                     # the test programs, same build directory
```

On HoreKA:

```bash
source env/horeka.sh gpu;  make GPU=1 GPU_ARCH=cc80     # A100 (cc90 for H100)
source env/horeka.sh cpu;  make
```

`GPU_ARCH` names the device (cc86 RTX 3060 / A6000, cc120 RTX 5090, cc80
A100, cc90 H100).  `BUILD=<dir>` keeps builds for different machines apart
on a shared home.  `FFTW_DIR` points the CPU build at FFTW.  No CMake, no
fypp, no Python in the build or the run.

## Run

From a directory containing `hst.in` (a Fortran namelist; the copy in this
repository is commented):

```bash
mpirun -np 4 /path/to/build-cpu/hst          # or build-gpu/hst, one rank per GPU
sbatch ~/hst/jobs/horeka_gpu.slurm           # HoreKA, 1 node, 4 GPUs
```

If `Dati.cart.out` exists it is read (with `time_from_restart = .true.`
the clock too); a field written by the CPL code works as well.  Otherwise
a seeded, divergence-free random field is generated (`&init`).

`timing = .true.` in `&time_control` prints the wall-clock time per phase
of the substep at the end of the run (`src/hst_timer.f90`).  `line_chunk`
in `&mesh` bounds the workspace of the line solver (x columns per batch;
0 = all columns on the GPU, 16 on the CPU, see `src/hst_linsolve.f90`).

## Output

All files are in the layout of the CPL code `hst-main`, so its
post-processing chain (`postprocess/`, `pressure_reconstruction/`) reads
them unchanged, and its fields can be used as restart files here.

- `Runtimedata`: one line per `dt_stat` with the CPL columns
  `time meanflowx meanflowy S S2 gamma_x gamma_y deltat cfl energy diss uw/2 vw/2`
  (integrals over the box height: `energy = ly/2 <u_i u_i>`,
  `diss = ly/2 <grad u : grad u>` without `nu`, computed here from the
  compact derivatives; CPL naming, so `uw/2` is `ly/2 <u v>` in our axes).
- `variances_runtime.dat`: `time uu vv ww uv` (CPL naming, `ly <..>`).
- `stokes_runtime.dat` (with a Stokes layer): `time energy_out energy_in
  diss_out diss_in`, region averages inside and outside the layer.
- `Dati.cart.out`: restart file, every `dt_save` and at the end.
- `fields/field<n>.fld`: velocity snapshots every `dt_field`;
  `p_fields/pField<n>.fld`: the pressure at the same times.  The two
  directories are created at start-up.

File layout (`src/hst_io.f90`): the CPL text header up to `Vfield=`, then
the complex128 array `(0..nx, -ny_cpl..ny_cpl, -1..nz_cpl+1)` of `(u, v, w)`
in C order with the four ghost rows, in CPL names (`ny_cpl` = our `nz`,
`nz_cpl` = our `ny + 1`, their `(v, w)` = our `(w, v)`).  Pressure files
are the same array without header and with one component.
`tests/compare_fields.py` reads both.

## Tests

```bash
tests/run_tests.sh build-cpu 2         # or build-gpu; second argument: ranks
```

- `test_roundtrip`: spectral-physical-spectral round trip and restart file,
  round-off on any rank count.
- `test_linsolve`: each system kind of the line solver (the two implicit
  systems, the Poisson equation, `d/dy`, the `D0` inverse) against the
  operator applied on the host through the ghost rows, round-off.
- `test_kelvin`: one Fourier mode in uniform shear against the closed-form
  Kelvin-mode solution with viscosity: 3.6e-4 at ny = 128 with the CPL
  treatment of the mean-shear advection (second order in dy, FINDINGS.md),
  4e-9 with `exact_shift = .true.`; also with constant and oscillating
  spanwise shear S2 (4e-9).
- `test_pressure`: Taylor-Green vortex against its exact pressure
  (1.9e-5 at ny = 64).
- `test_taylorgreen`: Taylor-Green vortices in two orientations decay
  exactly with the nonlinear terms on (7e-6 after 20 steps).
- `test_forcing`: the nonlinear forcing of `eta` against its closed form
  for `u = sin(ky)`, `w = sin(kx)`.
- `test_conservation`: energy input of each nonlinear product on a random
  field (`uu`, `vv`, `ww` vanish, the cross terms cancel).
- `test_stokes`: the mean profile of a body-force-driven Stokes layer
  against the analytic layer after the prescription window (2.3e-5).

The full solver agrees between CPU and GPU to 1e-13 after 50 steps, and
runs with different rank counts are bit-identical on the GPU.  Against the
CPL code `hst-main` from an identical start the box energy agrees to 5e-5
over 150 steps, and a run to S t = 100 gives S* = 6.3, -uv/q2 = 0.16 and
production = dissipation (FINDINGS.md).

## Performance

Seconds per full time step (three substeps), `examples/bench_*.in`, after
the performance pass of FINDINGS.md (the numbers in brackets are those of
the code before it):

| grid (dealiased) | 1 x A100 | 4 x A100 | 1 x RTX 3060 | istmio2 CPU, 4 ranks |
| --- | --- | --- | --- | --- |
| 64 x 128 x 64 | 0.019 (0.054) | | 0.17 (0.18) | 0.89 (1.03) |
| 256 x 256 x 256 | 0.115 (0.334) | 0.101 (0.140) | 1.25 (1.33) | |
| 512 x 512 x 512 | 0.97 (1.99) | 0.81 (0.98) | | |

The RTX 3060 runs double precision at 1/64 rate, so it gains little from
what helps the A100; `timing = .true.` prints where the time goes.  Four
A100 at 256^3 are limited by the alltoall (transposes are 85% of the
step), so a second node buys nothing before WP6.

## Status

| | |
| --- | --- |
| numerics, GPU, pressure, CPL files, S2, Stokes layer | done and validated (FINDINGS.md) |
| machines | istmio2, istmcetus, istmcorax (RTX 3060 / A6000 / RTX 5090), HoreKA (4 x A100 per node) |
| y decomposition, NCCL transport | not started (DESIGN.md 7, WP6) |
| safety net | `tests/run_tests.sh` (12 runs) and `tests/regression.sh` (three decks at 1e-10) on CPU and GPU |

## Layout

```
src/hst_params.f90       all state: mesh, parameters, clock, rank layout, fields
src/hst_input.f90        the namelist deck
src/hst_mpi.f90          x-z pencil decomposition, alltoall transpose, MPI-IO types
src/hst_fft.f90          FFTW / cuFFT (the only vendor-specific file besides hst_mpi)
src/hst_timer.f90        per-phase timer
src/hst_setup.f90        allocation and device mapping
src/hst_transforms.f90   spectral <-> physical, products, CFL
src/hst_initial.f90      seeded initial field
src/hst_io.f90           restart and snapshot files
src/hst_derivatives.f90  compact stencils, shear-periodic ghost rows
src/hst_linsolve.f90     cyclic pentadiagonal solves, one thread per mode, rows built on the fly
src/hst_equations.f90    the equations and one RK step
src/hst_stats.f90        Runtimedata
src/hst_pressure.f90     pressure at snapshot times
src/hst.f90              main program
tests/                   test programs on tests/test_common.f90, decks, runner, regression, field comparison
env/, jobs/              environment scripts and SLURM jobs
examples/                benchmark decks
```
