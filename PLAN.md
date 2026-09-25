# Plan: a GPU/CPU DNS code for homogeneous shear turbulence (HST)

Goal: take the Fortran `channel` DNS code (OpenMP-offload GPU solver for
turbulent channel flow) and turn it into a solver for homogeneous shear
turbulence, reproducing the numerics of the CPL code `hst-main`
(`scddns.cpl`, `S1data.cpl`, `linsolver_smw.cpl`).  The two source codes are
not modified; everything new lives in this repository.

Non-negotiable properties of the result:

- one binary, `hst`, built from ~12 plain `.f90` files with a Makefile
- runs on CPU (gfortran + MPI + FFTW) and on NVIDIA GPU (nvfortran + cuFFT)
- runs on istmio2 / istmcetus / istmcorax and on HoreKA with the same sources
- input is one Fortran namelist file, no Python, no runner scripts
- every kernel is readable Fortran: no fypp macros, no callback contexts,
  no autotuner; the first version has no NCCL and no y-decomposition, but
  is laid out so both can be added (section 7)
- pressure is available online, at snapshot times (3.7)
- vendor calls confined to two files, so an AMD backend is an addition
  (section 7)

---

## 1. What the two codes do, and what carries over

### 1.1 The CPL HST code (`hst-main`)

Velocity-vorticity formulation on a box that is Fourier in x (streamwise,
`alpha0`) and in the spanwise direction (`beta0`), and compact 4th-order
finite differences on a 5-point stencil in the shear direction, which the CPL
code calls z.  Mean flow `U = S z`.  Time integration: 3-stage RK (Rai-Moin)
for the nonlinear terms, Crank-Nicolson for viscous terms, and the mean-shear
advection `S z d/dx` integrated **analytically** (Sekimoto, Dong & Jimenez).

The three HST-specific ingredients (everything else is standard):

1. **Shear-periodic wrap.**  The z direction is periodic, but the field at
   `z + Lz` is the field at `x - S t Lz`.  In Fourier space this is a phase:
   every stencil point that wraps around the box is multiplied by
   `exp(i kx * gamma * (+-Lz))` with `gamma = S*time`.  This is the
   `EXP[I*alpha0*ix*S*time*(z(wrapped)-z(iz+k))]` factor that appears in
   every `D0/D1/D2/D4` and in the matrix corners.
2. **Integrating factor for the mean-shear advection.**  After each RK
   substep of length `ci*deltat`, the right-hand side and the stored explicit
   term are multiplied by `exp(-i kx S z ci deltat)` (the `~*EXP[-ialpha*S*z(iz)*ci*deltat]`
   lines in `buildrhs` and `linsolve`).  No grid remeshing, no tilted frame.
3. **Mean-shear tilting term** in the vertical-vorticity equation:
   `expl = ... - S*ibeta*D0(w)`.  The `D2w` (Orr-Sommerfeld) equation has no
   extra term because `U'' = 0`.

The linear systems per Fourier mode are pentadiagonal **plus six corner
entries** (cyclic band).  `linsolver_smw.cpl` solves them with a rank-6
Sherman-Morrison-Woodbury update.

Features of `hst-main` that are **out of scope** (stripped): unsteady
spanwise shear `S2`, Stokes layer, body force, `synthetic_init` and the CPL
post-processing tools.  The offline `pressure_reconstruction/` becomes an
online pressure computation (3.7).  Kept as a runtime switch
because it costs three lines and is needed for validation: `linear = .true.`
(nonlinear terms off).

### 1.2 The Fortran channel code (`channel`, branch `main_GPU`, 15k lines)

Same formulation (`v`, `eta`), same RK3/CN coefficients (`RK_rai` in
`channel_state.f90` equals `RK1/2/3_rai` in `scddnsdata.cpl`, substep lengths
8/15, 2/15, 1/3), same compact stencils (`setup_derivatives` in
`case_setup.fypp` is line-for-line the CPL `setup_derivatives`), same 3/2
dealiasing, same MPI x-z pencil transpose with an `alltoall`, same
`Runtimedata`/`Dati.cart.*.out` habits.

What carries over unchanged (after removing fypp and the `npy > 1` branches):

| channel file | reused as | notes |
| --- | --- | --- |
| `fft/ffts.fypp` (FFTW / cuFFT paths) | `hst_fft.f90` | drop HIP, drop byte_workspace: plain allocatables mapped once |
| `mpi/mpi_transpose.f90` (pack/unpack/alltoall, MPI-IO types, `init_MPI`) | `hst_mpi.f90` | drop y-slab, NCCL, HIP; `npy = 1` always |
| `numerics/channel_transforms.f90` | `hst_transforms.f90` | drop overlapping variant (one buffer) |
| `physics/channel_equations.fypp` `buildrhs_prepare`, `buildrhs`, `build_products`, `linsolve` recover-u-w | `hst_equations.f90` | macros written out; add the three HST ingredients |
| `numerics/stencil_coefficients.f90`, `setup_derivatives` | `hst_derivatives.f90` | drop all one-sided wall rows |
| `io/restart_io.f90` | `hst_io.f90` | header gains `S`, `ly`; keeps MPI-IO subarray write |
| `physics/initial_condition.f90` `uniform_from_key` | `hst_initial.f90` | seeded, decomposition-independent random field |
| `run/driver.f90` device selection + time loop | `hst.f90` | |
| `linsolve/y_line_solvers.fypp` custom penta kernel idea | `hst_linsolve.f90` | rewritten as a cyclic solver, see 3.4 |
| `io/pressure_output.fypp` (Poisson solve for p from the products) | `hst_pressure.f90` | shear-periodic Poisson with the cyclic solver; see 3.7 |

What is deleted outright: `y_schur_solver.f90` (958), `y_line_solvers.fypp`
(2512), `y_pipeline_nccl.f90`, `channel_nccl_p2p.c`, `mpi_autotune.f90`
(872), `convvelo*` (1650), `compact_component_solve.fypp` (487: callbacks
and wall closures), `wall_closure.fypph`, `channel_boundaries.fypph`,
`channel_bcs.f90`, `byte_workspace.fypp`, `env_options.f90`, `roctx.f90`,
`config.f90` (replaced by a namelist), all scalars (`nPhi`), flow-rate /
pressure-gradient correction, `post/`, `tests/*.py`, `postpro/`, cuSPARSE.
Section 7 records why each of the four larger deletions (y-decomposition,
NCCL, autotuner, cuSPARSE) was made and what it takes to bring it back.

Expected size of the result: about 2000 lines against 15000.

---

## 2. Conventions of the new code (decide once, here)

**Axes.**  Keep the channel layout, because every reused kernel is written
in it: `x` streamwise (Fourier, `alfa0`), `z` spanwise (Fourier, `beta0`),
`y` the compact-FD direction, which carries the shear: `U = S*y`.  So the CPL
`(x, y, z, u, v, w, zeta, D2w)` maps to our `(x, z, y, u, w, v, eta, d2v)`.
The README will say this in its first paragraph.

**Array layout.**  `V(iy, iz, ix, 1:3)`, complex, exactly as in the channel
code: `iy = -2:ny+1`, `iz = -nz:nz`, `ix = nx0:nxN` (this rank's x slab).
The unique grid points are `iy = 0 .. ny-1`, `y(iy) = iy*ly/ny`, uniform
(`ly = 2` by default, matching the CPL `zmax`).  The four ghost rows `-2, -1`
and `ny, ny+1` hold the shear-periodic images.  The stretched grid
(`htcoeff`) of the CPL code only served the Stokes layer and is dropped, but
`setup_derivatives` stays general so a non-uniform `y` costs nothing later.

**Unknowns.**  As in both codes: `V(:,:,:,2)` is `v`; `V(:,:,:,1)` is `eta =
i beta u - i alfa w` for every mode except `(0,0)`, where it packs the two
real mean profiles `(u_mean, w_mean)` as one complex number.  With `S` in the
base flow and no forcing, the `(0,0)` mode stays zero; it is kept for
symmetry and because the CPL code keeps it.

**Decomposition.**  `nproc` ranks split `nx+1` and `nzd` evenly (`npxz =
nproc`, `npy = 1`).  Each rank owns all of `y`, so every wall-normal solve is
rank-local, and the only communication per substep is the x-z `alltoall`
already in the channel code.  Constraint, inherited and checked at startup:
`nproc` divides `nx+1` and `nzd`.  The channel's slab indexing `ny0:nyN` is
kept in every array declaration and every kernel bound, and all line solves
go through one entry point, so that a distributed cyclic solve (7, item i)
can be added later without touching the physics kernels.

**GPU model.**  Unchanged from the channel code: OpenMP target offload with
`nvfortran -mp=gpu -cuda`, cuFFT through NVHPC's `cufft` module, all state
arrays mapped to the device once at startup (`!$omp target enter data`),
fields copied back only for I/O.  With gfortran the `!$omp target` lines are
comments and the same source runs on the CPU.  The two nvfortran rules from
`channel/src/physics/README.md` are kept: no `use` in modules that declare
`declare target` variables, and no device procedures called across module
boundaries (the cyclic solver is a `contains`-ed routine of the module whose
kernel calls it).

**Input.**  One namelist file `hst.in`:

```fortran
&mesh     nx=63, ny=128, nz=63, alfa0=2.0943951, beta0=6.2831853, ly=2.0 /
&physics  re=1000.0, s=1.0, linear=.false. /
&time     deltat=0.0, cflmax=1.0, t_max=100.0, nstep=1000000,
          dt_stat=0.01, dt_field=10.0, dt_save=10.0, time=0.0,
          time_from_restart=.false. /
&init     amplitude=1.0e-3, seed=1, kpeak=4.0 /
```

Fortran reads this in one statement per group; there is no parser to
maintain.  Unknown or missing keys are compile-time facts, not runtime
surprises.

**Output.**  `Runtimedata` (one line per `dt_stat`: `time, deltat, cfl,
2*tke, dissipation, <uv>, <uu>, <vv>, <ww>`, the same quantities as the CPL
`Runtimedata` + `variances_runtime.dat`), `Dati.cart.out` restart every
`dt_save`, `Dati.cart.<i>.out` snapshots every `dt_field`, all in the
channel's MPI-IO stream format with a 12-number header.

---

## 3. Numerics to implement (the parts that are actually new)

### 3.1 Shear-periodic ghost rows: `fill_ghosts(field)`

One kernel, run after every solve and before every transform, for each
`(iz, ix)` and component:

```
phase = exp(i * alfa0*ix * gamma * ly),   gamma = S*time   (mod lx/ly, to keep the argument bounded)
V(-1)   = V(ny-1) * conjg(phase)      V(-2)   = V(ny-2) * conjg(phase)
V(ny)   = V(0)    * phase             V(ny+1) = V(1)    * phase
```

(the sign follows `z(wrapped) - z(iz+k) = -+ly` in `S1data.cpl`; the Kelvin
mode test in 5.3 pins it down).  With the images in place, every explicit
stencil in `buildrhs_prepare` and `buildrhs` is the plain 5-point sum the
channel code already has, so those kernels port with the macros expanded and
nothing else changed.  The nonlinear products in the ghost rows come out
right automatically (a shifted field has a shifted product), which is why the
CPL code convolves rows `-2 .. nz+1` too.

### 3.2 Integrating factor: `shear_shift(dt_sub)`

After `buildrhs` has produced the RHS for `eta` and `d2v` (in `V(:,:,:,1:2)`)
and stored the explicit terms in `oldrhs`, multiply both by
`exp(-i*alfa0*ix*S*y(iy)*dt_sub)` with `dt_sub = 2/RK_rai(1,i)*deltat`, then
advance `time`.  This is the only place `y` appears in a phase.

### 3.3 Mean-shear tilting term

In `buildrhs_prepare`, the `eta` explicit part gets `- S*ibeta(iz)*D0(v)`
(accumulated through the same `rk_accum` bookkeeping as the nonlinear terms).
With `linear = .true.` the nonlinear `buildrhs` is skipped and the transforms
are only run for the CFL estimate.

### 3.4 Cyclic pentadiagonal solves: `hst_linsolve.f90`

Per substep and per mode there are three solves, exactly as in both codes:

| system | matrix rows `j = -2..2` | corner phase |
| --- | --- | --- |
| `d2v` | `lambda*(D2 - k2*D0) - ni*(D4 - 2 k2 D2 + k2^2 D0)` | yes |
| `eta` | `lambda*D0 - ni*(D2 - k2*D0)` | yes |
| `dv/dy` | `D0 x = D1 v` | yes |

The matrix is pentadiagonal on rows `0..ny-1` with the wrapped stencil
entries landing in the corners `(0,ny-2) (0,ny-1) (1,ny-1) (ny-2,0) (ny-1,0)
(ny-1,1)`, each multiplied by the wrap phase of 3.1.

Instead of the rank-6 Woodbury update of `linsolver_smw.cpl` (seven
pentadiagonal solves per line) use **bordering**: unknowns `x(ny-2)` and
`x(ny-1)` are the border.  Rows and columns `0..ny-3` form a *plain*
pentadiagonal matrix `P`, so

1. factor `P` once (in-place, the same forward elimination as the channel's
   custom penta kernel),
2. back-substitute three right-hand sides: `b`, and the two columns that
   couple to the border,
3. solve the 2x2 Schur complement for the border, then correct.

Three back-substitutions plus a 2x2 instead of seven solves and a 6x6, and no
dense LU on the device.  One GPU thread per `(iz, ix)` line, `ny` sequential
steps, complex arithmetic throughout (the corners make the matrix complex even
for `ix = 0`, where the phase is 1 and it degenerates to the real cyclic
case).  Assembly and solve are one kernel each; the assembly kernel is where
`lambda`, `ni`, `k2`, `der` and the phase meet, and it is 30 lines.

cuSPARSE is not used: `gpsvInterleavedBatch` cannot do cyclic systems, and
removing it removes a library, a workspace manager and an environment
variable.

### 3.5 Substep sequence (`hst.f90`)

```
do i = 1, 3
   fill_ghosts(V)                       ! images at the current gamma
   transform_to_physical()              ! V -> rVVdx, plus CFL on the last substep
   buildrhs_prepare(RK_rai(:,i))        ! implicit + carried explicit parts, tilting term
   transform_back_and_build_rhs(...)    ! six products -> spectral -> accumulate
   shear_shift(2/RK_rai(1,i)*deltat); time = time + 2/RK_rai(1,i)*deltat
   linsolve(RK_rai(1,i)/deltat)         ! assemble + cyclic solve for d2v, eta; dv/dy; recover u, w
end do
outstats()                              ! every dt_stat: reductions on device, one line to Runtimedata
```

This is `scddns.cpl`'s loop with `buildrhs -> time += -> linsolve -> pbc`
in the same order.

### 3.6 Statistics

`outstats` computes, on the device with `reduction(+:)`, the y-integrated
`sum |u|^2+|v|^2+|w|^2`, `Re(u conj(v))`, the three variances and the
dissipation `sum k2 |u|^2 + |du/dy|^2` (using `dv/dy` from the solve and
`du/dy, dw/dy` from a D1 stencil), the `ix = 0` plane counted once and all
others twice, then `MPI_Allreduce`.  This replaces both `getenergy()` and
`getvariances()` of `io.cpl`, minus the Stokes-layer branches.

### 3.7 Pressure, computed online

The pressure Poisson equation for the fluctuation `p` in a mean shear
`U = S y` is

```
(D2 - k2) p  =  -( ialfa^2 uu + ibeta^2 ww - D2 vv + 2 ialfa D1 uv + 2 ibeta D1 vw + 2 ialfa ibeta uw )  -  2 S ialfa v
```

with the same shear-periodic wrap as the velocity, so the left-hand side is
the `eta`-type cyclic matrix of 3.4 with `lambda = 0` and the corner phase of
3.1, solved by the same kernel.  The `(0,0)` mode is periodic and singular:
its solution is `p_00 = -vv_00 + const`, and the constant is fixed by zero
mean.  The six products are the ones `transform_back_and_build_rhs` already
forms; at an output time the code runs one extra product pass on the freshly
solved field (the products of the last substep belong to the field *before*
its solve), assembles the right-hand side with the `D1`/`D2` stencils on the
ghost-filled products, solves, and writes `Dati.cart.<i>.p.out` next to the
velocity snapshot in the same MPI-IO layout with one component.  This is the
channel's `pressure_output` reduced to a periodic box: no Neumann wall rows,
no `dp/dy` file (it is one stencil away in post-processing).  Cost: one
transform pass and one solve per snapshot, nothing per step.  The CPL
`pressure_reconstruction/prepare_pressure.cpl` is the offline equivalent.

### 3.8 Initial condition

Random solenoidal field: for each mode draw amplitude and phase from the
channel's decomposition-independent `uniform_from_key(seed, ...)`, scale by
`(k/kpeak)^2 exp(-(k/kpeak)^2)`, set `v` and `eta` directly (the solver
recovers `u, w` from continuity, so the field is divergence-free by
construction).  The CPL single-mode start (`initflow.cpl`) is a special case
(`kpeak` small, one mode) if ever needed for a side-by-side run.

---

## 4. Repository layout

```
homogenenousShearTurbulence/
  README.md              what it solves, conventions (2), build, run, file formats
  PLAN.md                this file
  Makefile               `make` (CPU, gfortran+FFTW)   `make GPU=1` (nvfortran+cuFFT)
  hst.in                 example deck (Sekimoto et al. 2016 box: Lx:Ly:Lz = 3:2:1, S=1)
  env/istm.sh            module-free: adds NVHPC 25.9 (/opt/Nvidia/nvhpc) and its MPI to PATH
  env/horeka.sh          module loads: toolkit/nvidia-hpc-sdk/25.3 (GPU) or compiler/gnu + openmpi + fftw (CPU)
  jobs/horeka_gpu.slurm  1 node, 4 GPUs, partition accelerated / accelerated-h100
  jobs/horeka_cpu.slurm  1 node, cpuonly
  src/
    hst.f90              program: initialize, time loop, finalize          (~150 lines)
    hst_params.f90       grid, wavenumbers, fields, S, ni, clock; use-free, declare target
    hst_input.f90        namelist read + echo
    hst_derivatives.f90  setup_derivatives (compact 5-point weights), fill_ghosts
    hst_fft.f90          FFTW / cuFFT plans and the four transforms
    hst_mpi.f90          init_MPI, pack/unpack/alltoall, MPI-IO datatypes
    hst_transforms.f90   transform_to_physical, build_products, transform_back_and_build_rhs, CFL
    hst_equations.f90    buildrhs_prepare (with tilting term), buildrhs, shear_shift, linsolve, recover u/w
    hst_linsolve.f90     assemble_* kernels and the cyclic pentadiagonal line solver
    hst_pressure.f90     Poisson right-hand side from the products, cyclic solve, p snapshot
    hst_io.f90           restart read/write, snapshots, outstats / Runtimedata
    hst_initial.f90      generate_initial_field, uniform_from_key
  tests/
    linsolve/            tiny program: cyclic solver vs dense complex LU, 1e-12
    kelvin/              hst.in with linear=.true., one mode; script compares with the closed form
    decay/               S=0 isotropic decay: CPU vs GPU vs 1-rank vs 4-rank
```

The Makefile has two blocks and an explicit source order; module
dependencies are resolved by that order, not by a generator.  CMake is not
needed and fypp is not needed (HoreKA has neither on the path).

---

## 5. Work packages, in order

Each package ends in something that runs; nothing in a later package is
needed to test an earlier one.

**WP0  Skeleton (repo, build, environment).**
Makefile, `env/*.sh`, `hst_params`, `hst_input`, `hst.f90` that reads the
deck, sets up MPI and the device, prints the configuration and exits.  Build
with gfortran here and with nvfortran 25.9 from `/opt/Nvidia/nvhpc`
(RTX 3060 on istmio2, sm_86; A6000 on istmcetus, sm_86; RTX 5090 on
istmcorax, sm_120) and on HoreKA with `toolkit/nvidia-hpc-sdk/25.3`.

**WP1  Infrastructure transplant (CPU first).**
`hst_fft`, `hst_mpi`, `hst_transforms`, `hst_io` restart, `hst_initial`:
copies of the channel modules with fypp expanded and `npy`, NCCL, HIP,
overlapping and workspace code removed.  Test: transform a field to physical
space and back, 1 and 4 ranks, error at round-off.

**WP2  HST numerics (CPU).**
`hst_derivatives` (stencils + `fill_ghosts`), `hst_linsolve` (assembly +
cyclic solver + `tests/linsolve`), `hst_equations` (ported `buildrhs*` with
macros written out, tilting term, `shear_shift`, `linsolve`).  First
end-to-end runs: `S = 0` decaying turbulence (energy monotone, dissipation
balance `d(tke)/dt = -eps` to discretisation accuracy), then `S = 1`.

**WP3  GPU.**
Add the `!$omp target` clauses to the WP2 kernels (they are the channel's
clauses with shorter shared lists), map arrays at startup, build with
nvfortran.  Test: GPU and CPU runs from the same seeded deck agree to a few
ULP after one step and to solver accuracy after 100; 1 GPU vs 2 GPUs
(istmcetus) bit-identical.

**WP3b  Pressure.**
`hst_pressure` (3.7), CPU then GPU.  Test: for a solenoidal field the
divergence of the momentum equation must vanish, so `sum_i D_i (RHS_i)` of the
velocity right-hand side equals the Poisson right-hand side to round-off;
and at `S = 0` the pressure of a Taylor-Green-like field has a closed form.

**WP4  HoreKA.**
Both builds, `jobs/*.slurm`, one node with 4 GPUs (`accelerated`) and one
with 4 H100 (`accelerated-h100`), timing per step for a
`nx=255, ny=512, nz=127`-class box, so the README can state the cost.

**WP5  Validation against physics and against `hst-main`.**
See section 6.  Then README, and a 40-line optional `postpro/read_field.py`
(numpy only) for people who want to look at a snapshot.

Optional, after WP5, if the CPL post-processing chain (`postprocess/`,
`pressure_reconstruction/`) is wanted on the new fields: a second snapshot
writer that emits the CPL `.fld` layout (text header, then the array
`(0..nx, -nz..nz, -1..ny+2)` of `(u, w, v)` with the CPL axis order).  About
60 lines of repacking on the host; it changes nothing in the solver.

---

## 6. Verification

1. **Solver unit test.**  Cyclic pentadiagonal solve against a dense complex
   LU for random matrices with the corner structure, `ny = 8 .. 64`.
2. **Round trip.**  Transform/transpose identity, 1 and 4 ranks.
3. **Kelvin mode (linear = .true.).**  A single Fourier mode in uniform shear
   has a closed-form solution (Kelvin 1887; Moffatt 1967): the wavevector
   tilts as `ky(t) = ky0 - S kx t` and the amplitudes follow an ODE with
   viscous decay `exp(-ni int k^2 dt)`.  This checks the wrap phase sign
   (3.1), the integrating factor (3.2), the tilting term (3.3) and the RK/CN
   coupling at once, to time-stepping accuracy.  The CPL code has the
   `linearNS` flag for the same purpose.
4. **Isotropic decay (S = 0).**  Energy budget `d(tke)/dt = -eps` closes to
   `O(deltat^2)`; CPU = GPU; rank count does not matter.
5. **Side by side with `hst-main`.**  Same box as `scddns.in` (`nx=96, ny=32,
   nz=192`, `S=1`, `Re=1000`), single-mode start in both codes, compare
   `Runtimedata` (energy, dissipation, `<uv>`) over the first 2-3 shear
   times; they should track to the difference in the initial random phase.
   Both compilers exist here (`~/bin/cpl`, `mpifort`).
6. **Long run.**  Sekimoto, Dong & Jimenez (2016) box `Lx:Ly:Lz = 3:2:1`
   (`alfa0 = 2pi/3, ly = 2, beta0 = 2pi`): compare `S* = S q^2/eps`,
   `<uv>/q^2` and the Reynolds-stress anisotropy with their table.

---

## 7. Decisions, with the reasoning and the way back

The four large deletions are simplicity calls, not necessities.  Each is
recorded with what it would take to reverse it, and the code is laid out so
that reversing it is an addition, not a rewrite.

- **(i) y-decomposition: `npy = 1` first, distributed cyclic solve as WP6.**
  The shear-periodic system *can* be solved distributed with the channel's
  y-Schur scheme: every y-slab eliminates its interior unknowns and leaves
  two interface unknowns per slab boundary; the reduced interface system is
  solved redundantly on all ranks after an allgather.  Periodicity adds one
  boundary (`npy` instead of `npy-1`) and makes the reduced system
  block-cyclic with the wrap phase in its corner blocks; that is a change to
  the reduced solve only.  Bordering (3.4) is the `npy = 1` case of the same
  algorithm.  It is deferred because it is the least validated part of the
  channel code (multi-rank GPU never ran on the ISTM boxes), because it is
  ~3500 lines with the y-slab transposes, and because for HST boxes the x-z
  split already reaches one rank per x column; `npy` matters when the global
  alltoall becomes the bottleneck at high GPU counts.  What is kept now so
  that WP6 is an addition: the `ny0:nyN` slab indexing in every array and
  kernel bound, and one line-solve entry point.
- **(ii) NCCL: MPI first, optional backend later.**  NCCL is only a transport
  for the alltoall (and the y pipeline).  CUDA-aware MPI (HPC-X on HoreKA,
  NVHPC's OpenMPI on the ISTM boxes) already does GPU-direct alltoall.  NCCL
  pays off on NVLink nodes (4x A100/H100) and costs ~600 lines, a C bridge
  and library discovery, and it is vendor-specific.  The alltoall lives in
  exactly one routine so the channel's `channel_comm_alltoall_complex` and
  `channel_nccl_p2p.c` can be dropped in as `comm = nccl` once HoreKA
  timings show what it buys.
- **(iii) Autotuner: dropped.**  It timed candidate `(npxz, npy)` splits,
  y-Schur pass hierarchies and pipelined-LU batch counts.  With `npy = 1`
  there is nothing to tune; if (i) is added, the split is one integer in the
  namelist, chosen from a hand-run scan and identical between runs.
- **(iv) cuSPARSE: dropped.**  `gpsvInterleavedBatch` solves plain
  pentadiagonal systems with one right-hand side and consumes the matrix; the
  cyclic system needs either seven solves of one matrix (Woodbury) or three
  right-hand sides against one factorisation (bordering), neither of which it
  offers without replicating the matrix in the batch.  Internally it is also
  one thread per line, the same as the channel's own custom kernel
  (`CHANNEL_YS_FORCE_CUSTOM_GPSV`), which is what 3.4 generalises.  The
  custom kernel is vendor-neutral.  The line solve is one routine, so a
  library call can replace it if profiling ever asks for it.
- **Pressure is computed online** (3.7), at snapshot times, with the same
  cyclic solver.
- **NVIDIA + CPU first, AMD later.**  The HIP/Cray branches of the channel
  code are not ported now, but the seams they need are kept: vendor calls
  (plan creation, transform execution, device synchronisation,
  `use_device_addr` around MPI) appear only in `hst_fft.f90` and
  `hst_mpi.f90`, behind `HAVE_CUDA` / `HAVE_FFTW` blocks with an empty
  `HAVE_HIP` slot; every kernel is OpenMP target only, no CUDA Fortran
  outside those two files.  Re-adding AMD is then the channel's existing
  hipfft/hipfort blocks in two files plus a Makefile block for `amdflang`
  or `ftn`.
- Namelist input instead of the INI parser; Makefile instead of CMake; no
  fypp.  Each removes a dependency that HoreKA does not provide by default.
- Uniform `y` grid, no `htcoeff`; the stencil code stays general.
- Cyclic solver by bordering, not by the rank-6 SMW of the CPL code.  Same
  answer, fewer operations, no dense solve on the device.
- The unsteady spanwise shear (`S2`), Stokes layer, body force and scalars are
  not ported.  `linear` is kept as a runtime flag for the Kelvin-mode test.

**WP6 (after WP5, on demand):** distributed cyclic y-Schur solve with
`npy > 1`; then the NCCL alltoall backend if timings justify it.

---

## 8. Findings during implementation

**Second-order error of the exact-advection step (WP2, Kelvin test).**
With shear off, the Kelvin-mode test agrees with the closed form to 1e-9.
With shear on, the error is 1.4e-3 at ny = 64 and 3.6e-4 at ny = 128 for a
tilt of ky from pi to 1.05 -- second order in dy, independent of the time
step and of viscosity.  The stencils themselves are sixth-order (modified
wavenumber error 4e-7 at ny = 64).  The cause is the exact-advection step
of 3.2.  The stored unknown is the D0-weighted Laplacian,
d2v(i) = sum_j d0_j (lap v)(y_{i+j}); for a mode exp(i ky y) that is the
Laplacian times the D0 symbol delta0(ky) = sum_j d0_j exp(i ky j dy)
= 1 - (1/6)(ky dy)^2 + ...  Multiplying the stored unknown by the node
phase exp(-i kx S y_i dt) tilts the mode to ky' but keeps the weight
delta0(ky) of the old wavenumber, and the solve at the new time divides by
delta0(ky').  Per substep the amplitude is off by delta0(ky)/delta0(ky');
the product over substeps telescopes to delta0(ky0)/delta0(ky(t)), so the
accumulated relative error is (1/6) dy^2 (ky0^2 - ky(t)^2): 1.43e-3 at
ny = 64 against 1.43e-3 measured, independent of dt and of nu, and
proportional to dy^2.  `S1data.cpl` shifts its D0-weighted right-hand
sides in the same way, so this is a property of the reference method,
kept deliberately.  The exact treatment would advect the unweighted
Laplacian: d2v_new(i) = sum_j d0_j (lap v)(y_{i+j}) exp(-i kx S y_{i+j} dt),
i.e. unweight with a D0 solve, apply the node phase, re-weight with D0, for
each shifted quantity (the two right-hand sides and the two carried
explicit terms): four extra line solves per substep, a local change to
shear_shift.  Applying the phase inside the stencil to v itself would be
wrong (it advects v instead of lap v and loses the Kelvin amplification).

This is implemented as the namelist switch `exact_shift` in `&physics`
(default `.false.`, i.e. the CPL method).  The prediction was checked on
four modes before the change (measured / predicted at ny = 64: 1.43e-3 /
1.43e-3, 3.56e-3 / 3.57e-3, 1.43e-3 / 1.43e-3, 5.69e-3 / 5.71e-3; the third
mode tilts through ky = 0 at twice the rate and gives the same error, as
only ky(t)^2 enters).  With `exact_shift = .true.` the Kelvin error is
7e-8 at ny = 64 and 4e-9 at ny = 128 (from 1.4e-3 and 3.6e-4), on CPU and
GPU alike (`tests/decks/kelvin_exact.in`).  Cost on the RTX 3060 for the
256^3 deck: 2.05 s/step against 1.33 s/step, i.e. the four extra line
solves add about 50% there; the A100 figure is to be measured.  On the
nonlinear side-by-side deck (ny = 191, section below) the two treatments
differ by 6e-5 in q2 after 150 steps, and the exact one is if anything
closer to the CPL run (1e-5 against 5e-5 at t = 0.3): at that resolution
the D0-weight error is already below the other differences between the
codes.

**Energy budget of isotropic decay (S = 0).**  On the deliberately coarse
16x32x16 deck the ratio -d(q2)/dt / (2 eps) stays within 5% of one after
the RK start-up transient, with eps from the compact derivatives.

---

## 9. Status (2026-09-25)

| package | state |
| --- | --- |
| WP0 skeleton, build, environments | done; builds on istmio2, istmcetus, istmcorax, HoreKA |
| WP1 transforms, transpose, restart I/O, initial field | done; round trip at round-off, 1-4 ranks, CPU and GPU |
| WP2 HST numerics on CPU | done; solver, Kelvin and decay checks |
| WP3 GPU | done; CPU = GPU to 1e-13, rank counts bit-identical, 4 x A100 on HoreKA |
| WP3b pressure | done; Taylor-Green test |
| WP4 HoreKA jobs and timings | done; 256^3 0.145 s/step and 512^3 0.99 s/step on 4 A100 (README) |
| WP5 validation against hst-main and Sekimoto et al. | done (section 8) |
| WP6 y decomposition, NCCL | not started |

Performance note: the small default deck (64x128x64) runs at 0.18 s/step
on the RTX 3060, only 7x faster than 256^3, so small grids are launch- and
latency-bound (many small kernels per substep, line batches of 16 x
columns).  Worth a pass later: larger line batches, fewer launches in
transform_to_physical, the per-line solver's memory traffic.

**Energy conservation of the nonlinear terms (WP5).**  `test_conservation`
takes one inviscid step from the random field with one product at a time:
uu, vv, ww alone change the energy by < 2e-5 per unit time relative, the
three cross products by 0.12, 0.14, 0.27 with sum 5e-3, and all six
together by 2e-5.  A viscous run at S = 0 from the same field satisfies
d(q2)/dt = -2 eps to 1e-4 at every step when q2 includes the (0,0) mode.
The first version of the initial field contained random mean profiles
(the (0,0) mode, 7% of the energy); their exchange with the fluctuations
made the fluctuation budget look 45% off and made the side-by-side start
differ from hst-main, which zeroes that mode at start-up.  The mode is now
left zero.  The nonlinear forcing of eta is checked separately against
its closed form (`test_forcing`, 3e-4 = O(dt)), and Taylor-Green vortices
in two orientations decay exactly with the nonlinear terms on
(`test_taylorgreen`, 7e-6).

**Side by side with `hst-main` (WP5).**  Both codes were started from the
same random field (ours, written by `tests/to_cpl_field.py` in the CPL
layout and read by `scddns` through `Vfield=`) on the `scddns.in` box
(nx = 96, 32 spanwise modes, 191 points over ly = 2, Re = 1000, S = 1),
fixed step 0.002, nonlinear.  Box-averaged energy and Reynolds stress:

| t | q2 hst | q2 CPL | rel. diff | uv hst | uv CPL | rel. diff |
| --- | --- | --- | --- | --- | --- | --- |
| 0.02 | 0.720601 | 0.720599 | 2e-6 | 0.033396 | 0.033394 | 3e-5 |
| 0.10 | 0.709871 | 0.709862 | 1e-5 | 0.030103 | 0.030096 | 2e-4 |
| 0.20 | 0.696121 | 0.696100 | 3e-5 | 0.024652 | 0.024636 | 6e-4 |
| 0.30 | 0.681737 | 0.681701 | 5e-5 | 0.016932 | 0.016909 | 1e-3 |

The remaining difference is at the level of the two codes' dealiasing
sizes and the CPL centred-difference dissipation; the CPL run needed
about 11 s/step on 4 CPU ranks against 0.66 s/step here on a shared
RTX 3060.

**CPL-compatible files (WP5).**  All output is now in the `hst-main`
layout: `Dati.cart.out` and `fields/field<n>.fld` with the CPL text header
and the C-ordered array with ghost rows, `p_fields/pField<n>.fld` headerless,
`Runtimedata` and `variances_runtime.dat` with the CPL columns (spanwShear
variant) as integrals over the box height.  Checked both ways on the
side-by-side box: `scddns` reads a file written here and reports the same
energy and Reynolds stress to all printed digits; a file written by
`scddns` is read here with its time and rewritten bit-identically; after
ten steps from the same field every `Runtimedata` and variance column
agrees with the CPL run to 1e-6 except the dissipation (3e-4, compact
against centred derivatives) and the CFL number (CPL subsamples it).
The earlier private format and `tests/to_cpl_field.py` are gone.

**Long sheared run (WP5).**  The default deck (box 3:2:1, 64x128x64 modes,
Re = 1000, S = 1, cflmax = 0.8) run to S t = 100 on the RTX 3060 (57324
steps, 0.18 s/step).  Averages over S t = 30..100 (701 samples):

| quantity | hst | literature |
| --- | --- | --- |
| production / dissipation, -S uv / eps | 1.006 | 1 in a statistically stationary box |
| S* = S q2 / eps | 6.3 | 5..7 (Rogers & Moin 1987; Sekimoto, Dong & Jimenez 2016) |
| -uv / q2 | 0.159 | 0.15 (Tavoularis & Karnik 1989) |
| b_uu, b_vv, b_ww | +0.10, -0.04, -0.06 | +0.2, -0.14, -0.06 at Re_lambda ~ 150..250 |
| Re_lambda | 33 | |

The anisotropy is weaker than the laboratory values, as expected at
Re_lambda = 33 in a small box.  Resolution of the default deck at these
statistics: dx/eta = 2.0, dy/eta = 1.0, dz/eta = 0.3, i.e. x is the coarse
direction; for a production run in this box choose nx about 3 nz (e.g.
nx = 191, nz = 63, ny = 128).  Snapshots and pressure files were written
every 20 time units.

---

## 10. Wished features (not yet implemented)

- **Unsteady spanwise shear `S2`** (`S2data.cpl`): a second, time-dependent
  mean-shear component `dW/dy = S2(t)` with its own shear-periodic
  displacement `gamma_y`.  Touches the wrap phase (it becomes
  `exp(-i (kx gamma_x + kz gamma_y))`), `shear_shift` (phase
  `exp(-i (kx S + kz S2) y dt)`), the tilting term (`+ S2 i alfa D0 v`),
  and the deck (`A`, `T`, `t0`).
- **Stokes layer and its body force** (`SLdata.cpl`): a prescribed
  oscillating spanwise profile `w(y, t)` localised at mid-box, imposed on
  the mean mode, plus the equivalent body force `fy(y, t)` on the mean
  `w` equation; needs the stretched grid (`htcoeff`) to resolve the layer,
  which the stencil code already supports.
- ~~CPL-compatible files~~: done, as the only format (section 8).
- **Pressure cadence**: the pressure is already computed only at snapshot
  times (`dt_field`); a separate `dt_pressure` would decouple the two.
- **y decomposition and NCCL transport** (section 7, WP6), and the
  small-grid performance pass (section 9).
