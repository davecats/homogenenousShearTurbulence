# hst

Direct numerical simulation of homogeneous shear turbulence on CPUs and
NVIDIA GPUs, in about 2000 lines of plain Fortran.  Derived from the
[`channel`](https://github.com/davecats/channel) DNS code and from the CPL
HST code; the design, the reasoning behind every simplification and the
findings made on the way are in [PLAN.md](PLAN.md).

## What it solves

Incompressible Navier-Stokes fluctuations about the mean flow `U = S*y` in a
box that is periodic in `x` (streamwise) and `z` (spanwise) and
shear-periodic in `y`: Fourier in `x` and `z` with 3/2 dealiasing, compact
(sixth-order, five-point) finite differences in `y`, velocity-vorticity
(`v`, `eta`) formulation, RK3 (Rai-Moin) for the nonlinear terms,
Crank-Nicolson for the viscous ones, and the mean-shear advection
integrated exactly by a phase shift (Sekimoto, Dong & Jimenez).  The
pressure is computed online at snapshot times.

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
the clock too); otherwise a seeded, divergence-free random field is
generated (`&init`).

## Output

- `Runtimedata`: one line per `dt_stat`:
  `time deltat cfl q2 eps uv uu vv ww`, box averages of the fluctuations
  (`q2 = <u_i u_i>`, `eps` the dissipation, then the Reynolds stresses).
- `Dati.cart.out`: restart file, every `dt_save` and at the end.
- `Dati.cart.<i>.out`: velocity snapshots every `dt_field`, and
  `Dati.cart.<i>.p.out` the pressure at the same time.

File format (`src/hst_io.f90`): 3 int32 `nx ny nz`, 6 float64
`alfa0 beta0 ly re S time`, then the complex128 array `(ny, 2nz+1, nx+1, 3)`
of the `(u, v, w)` modes in Fortran order (the pressure file has one
component).  `tests/compare_fields.py` reads it.

## Tests

```bash
tests/run_tests.sh build-cpu 2         # or build-gpu; second argument: ranks
```

- `test_roundtrip`: spectral-physical-spectral round trip and restart file,
  round-off on any rank count.
- `test_linsolve`: the cyclic pentadiagonal line solver against exact
  solutions, round-off.
- `test_kelvin`: one Fourier mode in uniform shear against the closed-form
  Kelvin-mode solution with viscosity (3.6e-4 at ny = 128; see PLAN.md 8
  for why this is second order in dy).
- `test_pressure`: Taylor-Green vortex against its exact pressure
  (1.9e-5 at ny = 64).
- `test_taylorgreen`: Taylor-Green vortices in two orientations decay
  exactly with the nonlinear terms on (7e-6 after 20 steps).
- `test_forcing`: the nonlinear forcing of `eta` against its closed form
  for `u = sin(ky)`, `w = sin(kx)`.
- `test_conservation`: energy input of each nonlinear product on a random
  field (`uu`, `vv`, `ww` vanish, the cross terms cancel).

The full solver agrees between CPU and GPU to 1e-13 after 50 steps, and
runs with different rank counts are bit-identical on the GPU.  Against the
CPL code `hst-main` from an identical start the box energy agrees to 5e-5
over 150 steps, and a run to S t = 100 gives S* = 6.3, -uv/q2 = 0.16 and
production = dissipation (PLAN.md, section 8).

## Performance

Seconds per full time step (three substeps), `examples/bench_*.in`:

| grid (dealiased) | 1 x A100 | 4 x A100 | 1 x RTX 3060 |
| --- | --- | --- | --- |
| 256 x 256 x 256 | 0.33 | 0.145 | 1.33 |
| 512 x 512 x 512 | | 0.99 | |

## Layout

```
src/hst_params.f90       all state: mesh, parameters, clock, rank layout, fields
src/hst_input.f90        the namelist deck
src/hst_mpi.f90          x-z pencil decomposition, alltoall transpose, MPI-IO types
src/hst_fft.f90          FFTW / cuFFT (the only vendor-specific file besides hst_mpi)
src/hst_setup.f90        allocation and device mapping
src/hst_transforms.f90   spectral <-> physical, products, CFL
src/hst_initial.f90      seeded initial field
src/hst_io.f90           restart and snapshot files
src/hst_derivatives.f90  compact stencils, shear-periodic ghost rows
src/hst_linsolve.f90     cyclic pentadiagonal solves, one thread per mode
src/hst_equations.f90    the equations and one RK step
src/hst_stats.f90        Runtimedata
src/hst_pressure.f90     pressure at snapshot times
src/hst.f90              main program
tests/                   test programs, decks, runner, field comparison
env/, jobs/              environment scripts and SLURM jobs
examples/                benchmark decks
```
