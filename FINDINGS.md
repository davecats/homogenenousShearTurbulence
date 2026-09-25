# Findings made while building and validating the code

Kept in the order they were made.  The design they refer to is in
[DESIGN.md](DESIGN.md); section numbers below are its sections.

## Second-order error of the exact-advection step, and `exact_shift`

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

## Energy budget of isotropic decay (S = 0)

On the deliberately coarse
16x32x16 deck the ratio -d(q2)/dt / (2 eps) stays within 5% of one after
the RK start-up transient, with eps from the compact derivatives.

---


## Performance of the small grids

The small default deck (64x128x64) runs at 0.18 s/step
on the RTX 3060, only 7x faster than 256^3, so small grids are launch- and
latency-bound (many small kernels per substep, line batches of 16 x
columns).  Worth a pass later: larger line batches, fewer launches in
transform_to_physical, the per-line solver's memory traffic.

## Energy conservation of the nonlinear terms (WP5)

`test_conservation`
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

## Side by side with `hst-main` (WP5)

Both codes were started from the
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

## CPL-compatible files (WP5)

All output is now in the `hst-main`
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

## Unsteady spanwise shear S2 (DESIGN.md section 10 wish, done)

`&physics`
takes `s2_amplitude`, `s2_period` (0 = constant) and `s2_start`, giving
the CPL `S2data.cpl` law `S2(t) = A sin(2 pi (t - t0)/T)` for `t >= t0`
(`A`, `T`, `t0_SL` there).  It enters in four places: the wrap phase of
the ghost rows and of the line-solve corners, now
`exp(-i (kx gamma_x + kz gamma_y))` with `gamma_y = int S2 dt` in closed
form; the exact-advection phase over a substep,
`exp(-i (kx S dt + kz int S2 dt) y)` (both variants of `shear_shift`);
the tilting term of eta, `(S2 i alfa - S i beta) D0 v`; and the rapid
term of the pressure, `-2 (S i alfa + S2 i beta) v`.  `S2` and `gamma_y`
are written to `Runtimedata` and to the field headers.  One deliberate
difference from `S2data.cpl`: there the carried explicit term of a
substep is shifted with the displacement of the *previous* substep
(`delta_gamma` is updated only after `buildrhs`), here with the current
one, as in `S1data.cpl`.

Checks: the Kelvin test generalises (`ky(t) = ky0 - S kx t - kz gamma_y`,
`eta` closed-form for constant S2, `v` for any S2 with the viscous
integral done numerically): 4.6e-9 for constant `S2 = 0.7` and 4.1e-9 for
`S2 = 0.7 sin(2 pi t/1.5)` with `exact_shift`; with the default advection
5.5e-5, again exactly the `(1/6) dy^2 (ky0^2 - ky(t)^2)` prediction with
the end points of the tilt.  Nonlinear side-by-side with `scddns`
(`A = 0.6, T = 0.5, t0 = 0`, 20 steps from the same field): `S2` and
`gamma_y` agree to 1e-15, energy to 2e-6, `uw/2` to 7e-5, `vw/2` to
2e-6 absolute.

## Stokes layer and stretched grid (DESIGN.md section 10 wish, done)

`&mesh`
takes `ystretch` (the `htcoeff` tanh clustering at mid-box of
`scddnsdata.cpl`; the stencil weights were already general, the row
spacing `dyl` now enters the CFL estimate and the statistics as
integration weights).  `&physics` takes `sl_amplitude`, `sl_period`,
`sl_delta`, `sl_start`, `sl_bodyforce` (default true, the CPL
`bodyforce`/`bf_dvw` flags) and `sl_ramp` (the `smoothStep` flag).
`hst_stokes.f90` holds the profile
`W = A exp(-s) cos(omega (t - t0) - s)`, `s = sqrt((y - ly/2)^2/delta^2 + 0.01)`,
and its body force `f = dW/dt - nu W''` in closed form (checked against
the `bodyF` expression of `SLdata.cpl` term by term).  With the body
force the mean `w` equation carries `f` and drops its Reynolds-stress
divergence (`bf_dvw`), so the mean profile is exactly `W` whatever the
turbulence does; the profile is prescribed outright during the first two
periods and on the two edge rows, as `apply_SL` does.  Without the body
force it is prescribed at every substep.  Two differences from the CPL
code, both deliberate: the force is added D0-weighted like every other
term (`SLdata.cpl` adds it raw, a second-order error in the applied
force), and the two-period window counts from `sl_start` rather than from
time zero.  `stokes_runtime.dat` carries the region averages of
`<u_i u_i>` and `<grad u : grad u>` inside and outside `|y - ly/2| < 8 delta`
(the `energy_in/out`, `diss_in/out` of the CPL `Runtimedata`).

Checks: `test_stokes` starts from no fluctuations and runs three
periods, one of them beyond the prescription window; the mean profile
then agrees with the analytic layer to 2.3e-5 (ny = 128, `ystretch = 3`,
`delta = 0.05`).  The error is spatial and comes from the `eps = 0.1`
smoothing of `|y|`, a feature of width `eps delta`: with `ystretch =
1.5` it is 2.6e-3, unchanged by halving the time step, and falls to
1.1e-4 with twice the points or twice `delta`.  The Kelvin test on the
stretched grid gives 7e-8.  No side-by-side with `hst-main` was possible:
its `StokesLayer` build does not compile as delivered (`SLdata` is used
after `S1data`, which needs its `fy`, and the S2 header bookkeeping is
undeclared on that path, and `updategamma` of the S2 path is called
unguarded); three scratch-copy fixes were tried and the build still
failed, so that path of `hst-main` is unmaintained.

## Long sheared run (WP5)

The default deck (box 3:2:1, 64x128x64 modes,
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


## Performance pass (cleanup and performance session, 2026-09-25)

Baselines of the code as it stood before this session, seconds per full
step (three substeps), `examples/bench_*.in`:

| deck | grid | RTX 3060, 1 rank | istmio2 CPU 1 / 4 ranks | A100 1 / 4 GPUs |
| --- | --- | --- | --- | --- |
| bench_64 | 64 x 128 x 64 | 0.180 | 1.86 / 1.03 | 0.054 / |
| bench_256 | 256^3 | 1.33 | | 0.334 / 0.140 |
| bench_512 | 512^3 | | | 1.99 / 0.98 |

**Where the time went.**  A per-phase timer (`timing = .true.`) and nsys
kernel summaries, on the RTX 3060 and on one A100, before any change:

- The line solver dominated on both machines: the solve kernel plus its
  assembly kernel were 30% of the GPU time on the RTX 3060 and 51% on the
  A100 (bench_256).  The assembly kernel alone (writing the five diagonals
  of every line and an `exp` per element) was 19% on the A100.
- On the RTX 3060 the four cuFFT calls were 40% of the time, on the A100
  12%: the GeForce card runs double precision at 1/64 of its single
  precision rate, so every FP64-heavy kernel (the FFTs, the complex
  divisions of the solver) is compute-bound there, while the A100 is
  bandwidth-bound.  Timings on the RTX 3060 are therefore not a guide to
  the A100 for the FFT and solver shares.
- `compute_cfl` read the physical field with y innermost, i.e. strided by
  a whole plane: 33 ms per call on the RTX 3060 (3% of the time) for a
  reduction over 900 MB.
- The 64^3-class deck was 83% kernel-busy on the RTX 3060 with 61 launches
  per substep, so launch latency was a smaller part than expected; on the
  A100 the same deck ran at 0.031 s/step after the solver batches were
  widened (below), 1.7x the original.
- The transposes (pack, alltoall, unpack) were 8% on one A100; the
  512^3 run on four A100 reached 51% parallel efficiency.

**What was done, in order, each with the regression at 1e-10 (most of
them bit-identical):**

1. *All x columns in one line-solver batch* on the GPU (`line_chunk`,
   default 16 before): 0.180 -> 0.154 s/step on the RTX 3060 for bench_64,
   0.054 -> 0.031 on the A100.  On the CPU the larger workspace falls out
   of the cache (1.86 -> 2.05 s/step), so the CPU default stays 16.
2. *The line solver generates its rows on the fly.*  One thread per line
   builds each row from `der`, `k2` and the wrap phase, eliminates with the
   two previous rows held in registers, and forward-substitutes the three
   right-hand sides of the bordering scheme in the same sweep; only the
   three upper diagonals are stored for the back substitution.  The
   assembly kernel and two of the five stored diagonals are gone; same
   operations in the same order, so the fields are bit-identical.  On the
   CPU: 0.60 -> 0.45 s/step in the solves of bench_64 (1 rank).
3. *`compute_cfl` with x innermost*: 33 -> ~3 ms on the RTX 3060.
4. *Three fields per transform and per transpose*: u, v, w together, the
   six products in two groups of three, the x transforms in place (the
   real buffers are pointer views of the complex ones, found on the device
   through the mapping of their buffer).  Launches per substep 100 -> ~35,
   collectives 9 -> 3.
5. *cuFFT on the OpenMP target stream* (`ompx_get_cuda_stream`): the
   eight `cudaDeviceSynchronize` per substep inherited from the channel
   code are gone; nsys shows every kernel on one stream.
6. *FFTW_MEASURE instead of FFTW_PATIENT*: planning of a 256^3 run on
   four CPU ranks 436 s -> 16 s at the same execution speed.  `-O3
   -march=native` was tried and gave nothing (2.10 against 2.04 s/step).

**Tried and dropped.**  Padding the x rows to 128 bytes or to a power of
two for cuFFT: the 3/2 length 384 = 128 x 3 uses cuFFT's "regular" kernels,
but 512 costs more than it gains (18.4 against 15.2 ms on the RTX 3060 for
the real transform, 6.7 against 5.1 ms for the complex one) and 400 or
392 change nothing.  In-place against out-of-place: no difference, so the
in-place layout is free.

The RTX 3060 was shared with another job during the second half of the
session; its numbers after item 3 are not comparable with the earlier
ones (a back-to-back A/B of two builds under the same load is).

**Result** (same decks and machines as the baseline table above; A100
numbers from `jobs/horeka_bench_all.slurm` on the same day):

| deck | RTX 3060, 1 rank | istmio2 CPU 1 / 4 ranks | A100 1 / 4 GPUs |
| --- | --- | --- | --- |
| bench_64 | 0.17 (shared card) | 1.54 / 0.89 | 0.0186 / |
| bench_256 | 1.25 (shared card) | | 0.115 / 0.101 |
| bench_512 | | | 0.970 / 0.809 |

The A100 is 2.9x faster on one GPU at 64^3-class and 256^3 and 2.1x at
512^3; four GPUs gained less (1.4x at 256^3, 1.2x at 512^3) because their
step is 85-90% transposes (pack, alltoall, unpack; `NP=4` profile), which
the batching made larger but not fewer bytes.  That is the next session's
subject (NEXT_SESSION.md).  One more memory item came out of the 512^3
run on one A100: with three-field batches cuFFT's per-plan work areas no
longer fit next to the fields (out of memory at 27 GB of fields and
buffers), so the four plans now share one work area
(`cufftSetAutoAllocation` off); the line-solver workspace with all
columns is 6.1 GB there (`line_chunk` bounds it).

---


## Multi-GPU pass (2026-09-25, session 3)

The handoff said four A100 were 85-90% transposes, inferred from the
phases the timer had then.  The first step was to make the timer say so
itself: pack/unpack and the alltoall are now phases of their own, charged
from inside `hst_mpi` (the transform and product phases keep the FFTs and
kernels around them).  Two things came out of that.

**The MPI alltoall was the step.**  On four A100 (`jobs/horeka_profile.slurm`,
`NP=4`, seconds per full step):

| deck | step | alltoall | share |
| --- | --- | --- | --- |
| bench_256 | 0.101 | 0.061 | 60 % |
| bench_512 | 0.750 | 0.534 | 68 % |

At 512^3 each GPU sends 0.92 GB per transpose to the other three, nine
times per step: 8.2 GB/step in each direction at 0.534 s is 15 GB/s per
GPU, a twentieth of what the NVLink of the node carries.  HPC-X's
MPI_Alltoall on device buffers is not using it well.

**NCCL transport** (`make GPU=1 NCCL=1`, deck `transport`, default
`auto`).  Not the channel code's C bridge: `hst_mpi.f90` declares the six
NCCL C prototypes it needs (90 lines) and issues one `ncclSend` and one
`ncclRecv` per peer inside a group (NCCL 2.25/2.27 on the two machines has
no alltoall), on the OpenMP target stream from `hst_fft::target_stream`,
so the unpack kernel queues behind the transfer and the host never waits
for it.  NVHPC's own Fortran `nccl` module was tried first and rejected:
its interfaces want CUDA Fortran `device` arrays, which OpenMP-mapped
arrays are not.  The NCCL communicator is created in `setup_decomposition`
from a unique id broadcast over MPI; `auto` falls back to MPI when ranks
share a GPU (two ranks on the RTX 3060), which NCCL refuses.  Validation:
the 12 tests with two ranks on the two RTX A6000 of istmcetus and the
regression against the references at 1e-13; the 4-GPU test job on HoreKA.

Same job, same day, the two transports back to back:

| deck | alltoall MPI | alltoall NCCL | step MPI | step NCCL | 1 A100 |
| --- | --- | --- | --- | --- | --- |
| bench_256 | 0.061 | 0.0096 | 0.101 | 0.048 | 0.115 |
| bench_512 | 0.534 | 0.060 | 0.750 | 0.298 | 0.970 |

The alltoall is 6-9x faster (137 GB/s per GPU per direction at 512^3),
the step 2.1-2.5x, and four A100 are now 3.3x one at 512^3 (81% parallel
efficiency; the original channel code reached 51% there) and 2.4x at
256^3.

**A timer artefact, found with nsys.**  The first split charged 17% to
"pack, unpack" on the A100 and 38% on the RTX 3060, but the nsys kernel
summary of the same run showed the pack/unpack kernels at 5% of the
kernel time and near bandwidth (250-380 us for 39 MB each on the RTX
3060).  The z transform before each transpose is a cuFFT call that
returns at once, and the pack's toc was the next synchronisation point,
so the transform's time landed in the pack.  `hst_transforms` now marks
its phase before each transpose.  The alltoall line was never affected
(the pack kernel is synchronous, so the transfer starts on a quiet
device).  Corrected split below.
