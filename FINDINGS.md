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

**Where the time goes now.**  Four A100 with NCCL, corrected timer, share
of the step (`jobs/horeka_profile.slurm`, `NP=4`):

| phase | bench_256 | bench_512 |
| --- | --- | --- |
| products, FFTs, buildrhs | 24 % | 31 % |
| transpose pack, unpack | 15 % | 18 % |
| transpose alltoall | 18 % | 16 % |
| implicit solves | 21 % | 13 % |
| ghosts, dv/dy, u and w | 12 % | 7 % |
| to physical: FFTs, CFL | 7 % | 10 % |

No single item dominates any more.  The pack/unpack kernels run at about
half the A100's bandwidth (2.4 GB in 3 ms at 512^3), an uncoalesced write
side; a tiled transpose would take them to 1/3 of that.  The alltoall at
16% is the upper bound of what overlapping it with the transforms of the
next batch (the channel's double buffering) could hide, at the price of
splitting the three-field batches again; not done, the pack kernels and
the products are the cheaper targets.

On one A100 (bench_256, nsys kernel summary of the same day): line solver
26% (three kernels), `build_products` 13.5%, `buildrhs` 12.3%, the local
repacks between the two pencil layouts 14.6% (`repack_xTOz_local` 2.0 ms
per call, 1.2 GB moved, i.e. 40% of bandwidth: the same uncoalesced
transpose as the pack kernels), cuFFT 22%.  On one rank the repack is pure
overhead of the pencil abstraction; cuFFT could transform the z lines in
the x-pencil layout with a strided plan, one call per y row, or the repack
could become a tiled transpose.

**Reciprocal pivots in the line solver** (handoff Part 1, item 2): the
forward sweep divided by the pivot twice per row and the back substitution
three times; `U(:, :, 0)` now stores 1/pivot and the sweeps multiply.  RTX
3060, bench_64, one rank, back to back: implicit solves 0.0286 -> 0.0186
s/step, the step 0.146 -> 0.131.  Regression 1e-14 (CPU) and 6e-14 (GPU)
against the stored references, which were not updated.

**Tried and dropped: `buildrhs` with `iz` innermost** (handoff Part 1,
item 1).  The loop has `iy` innermost so that `rhs` is written
contiguously while the fifteen stencil reads of `VVdz` are a plane apart
between neighbouring threads; the guess was that reading contiguously and
writing strided would win.  It loses, 2.44 -> 4.47 ms per call on the
A100 (bench_256, `~/hst-exp` against `~/hst` back to back): with `iy`
innermost the five-point windows of neighbouring threads overlap four
fifths, so the L1 cache serves most of the reads, and the write side is
the cheap one.  The loop stays as it was.  Same A/B, one A100: the
reciprocal pivots take the implicit solves from 0.0197 to 0.0186 s/step
(6%) and the step from 0.1157 to 0.1128.

**Result** of the session (`jobs/horeka_bench_all.slurm`, same day, all
changes in; the README table):

| deck | 1 A100 | 4 A100 | before the session, 1 / 4 |
| --- | --- | --- | --- |
| bench_64 | 0.0160 | | 0.0186 / |
| bench_256 | 0.113 | 0.0425 | 0.115 / 0.101 |
| bench_512 | 0.955 | 0.296 | 0.970 / 0.809 |

The single-GPU gain is the solver's reciprocal pivots; the 4-GPU gain is
NCCL.  What is left on four GPUs is spread over the products and their
transforms, the pack/unpack kernels, the alltoall and the solver, none of
them above a third of the step (NEXT_SESSION.md).

---


## The kernels around the transposes (2026-09-25, session 4)

The handoff named six pack/unpack/repack kernels as index-swapping loops
at half the bandwidth.  Reading them, only four are: the two `pack`
kernels copy whole blocks with the leading index kept (the alltoall
permutes blocks), and nsys on the RTX 3060 confirms it, 2.1 ms per call
for 0.61 GB moved, 81% of the card's bandwidth.  The change of leading
index happens in the two `unpack` kernels (4.7 ms for the same bytes,
36%) and, on one rank, in the two local repacks (6.0 ms for 1.23 GB,
57%; 40% on the A100).

**One tiled transpose kernel.**  The four became one routine,
`transpose_tiled`: `B(jb*(block-1) + j, i, plane) = A(i, j, plane, block)`
for every plane (a y row of one field) and block (a source rank), called
from the two transpose drivers with the layout of each use, and the
receive buffer viewed as `(nzB, nxB, ny+4, 3, npxz)`.  On the GPU a
thread block moves a 32x32 tile through shared memory (padded by one
against bank conflicts), reading `A` along `i` and writing `B` along
`j`, so both sides are coalesced.

**It had to be CUDA Fortran.**  Four OpenMP forms of the same tile were
built and timed on the RTX 3060 (bench_256, one rank, the local repack,
plain loop 6.0 ms per call):

| form | tile in | per call |
| --- | --- | --- |
| `teams distribute` + `parallel do` inside, 128 threads | shared memory (26 KB, the compiler says so) | 10.0 ms |
| the same with `thread_limit(256)` | shared memory | 16.4 ms |
| `teams loop` + `loop bind(parallel)` | global memory (0 bytes shared in the launch) | 7.6 ms |
| 4x4 blocks per thread, no tile | registers (compiler fused the loops away) | 13.4 ms |
| CUDA Fortran kernel, `shared` tile, 32x8 threads | shared memory (16.9 KB static) | **3.8 ms** |

The `teams distribute` form puts the team-private array in shared memory
but generates the inner `parallel do` loops badly (109 registers per
thread, no parallelisation report for them); the `loop` form generates
good loops but leaves the tile in global memory; the OpenMP 5 `allocate`
directive with `omp_pteam_mem_alloc`, which would fix that, is
"unrecognized" by nvfortran 25.9, as a clause and as a directive.  The
CUDA Fortran kernel (40 lines, `attributes(global)`, launched on the
OpenMP target stream from device pointers obtained through
`use_device_addr`) is the one kernel of the code that is not OpenMP; the
CPU path is the plain loop.  DESIGN.md 7 already confined CUDA Fortran to
`hst_mpi.f90` and `hst_fft.f90`.  Two pitfalls: Fortran is
case-insensitive, so a tile array named `tile` next to the parameter
`TILE` is one symbol ("vector expression used where scalar expression
required"); and an explicit-shape dummy inside a target region is mapped
implicitly, which works because the actual is already on the device.

RTX 3060, bench_256, per call: local repack 6.0 -> 3.8 ms (323 GB/s, 90%
of the card), unpack (two ranks) 4.7 -> 2.1 ms, the same as the pack
copy.  Regression bit-identical to the references (5e-14, the CPU/GPU
difference), 12 tests on CPU, GPU and the NCCL build.

**`build_products` in one pass** (handoff item 3).  The product index was
the outer collapsed loop, so each velocity field was streamed from memory
twice per call (six reads for three products; the RTX 3060 was at 355
GB/s, i.e. the reads were not served by L2).  A thread now reads u, v, w
once and writes its three products, which also reads better than the
`merge` index arithmetic it replaces.  Same operations in the same order,
bit-identical; 7.8 -> 5.7 ms per call on the RTX 3060.

**`buildrhs_prepare` and `buildrhs`, evaluated and left alone.**  On one
A100 at 256^3 they take 1.5 and 2 x 2.4 ms per substep, both at about
half the bandwidth (1.2 GB and 1.9 GB of sectors moved; the `VVdz`
stencil reads are a plane apart between neighbouring threads, so half of
every 32-byte sector is wasted, and the other loop order was tried last
session and lost).  Fusing `prepare` into the first `buildrhs` saves one
read and write of `rhs` and `oldrhs` (1.07 GB, about 0.7 ms) per substep,
2% of the step, and merging the two `buildrhs` calls would save the same
at the price of six products in memory at once; neither is worth the
loss of the "V-only part, then the products" structure.

**A100 numbers** (`jobs/horeka_ab.slurm`, new in this session: the same
benchmarks on `~/hst` and `~/hst-exp` back to back on one node).

Job 5163733, `~/hst` at the previous commit against `~/hst-exp` with the
tiled transpose, seconds per full step from the timer (phase "transpose
pack, unpack" and the total; on four GPUs with NCCL):

| deck, GPUs | pack/unpack before | after | step before | after |
| --- | --- | --- | --- | --- |
| bench_256, 1 | 0.0166 | 0.0085 | 0.1133 | 0.1052 |
| bench_256, 4 | 0.0075 | 0.0049 | 0.0454 | 0.0429 |
| bench_512, 1 | 0.1456 | 0.0655 | 0.9567 | 0.8766 |
| bench_512, 4 | 0.0549 | 0.0351 | 0.3056 | 0.3040 |

The nsys summary on one A100 (bench_256) has the local repack at 1.86
and 1.94 ms per call before and the tile kernel at 0.92 ms after (1.2 GB
moved: 1.3 TB/s, 85% of the A100's bandwidth); the kernel time of five
steps 0.596 -> 0.566 s.  One GPU gains 7-8%, four GPUs 5% at 256^3.  At
512^3 on four GPUs the pack/unpack phase lost 0.020 s but the alltoall
phase gained 0.018 s in this sample (0.0547 -> 0.0732; the previous
session measured 0.060), so the step barely moved.  The second job below repeated the tiled
build and found its alltoall at 0.0512 and its step at 0.2825, so that
sample was NCCL noise (the alltoall of the same build varies by 30%
between runs; the pack/unpack phase is reproducible to 1%).

Job 5163750, `~/hst-exp` (tiled) against `~/hst-exp2` (tiled and the
one-pass `build_products`), same layout:

| deck, GPUs | products phase before | after | step before | after |
| --- | --- | --- | --- | --- |
| bench_256, 1 | 0.0466 | 0.0390 | 0.1060 | 0.0981 |
| bench_256, 4 | 0.0126 | 0.0102 | 0.0436 | 0.0409 |
| bench_512, 1 | 0.3766 | 0.3202 | 0.8769 | 0.8201 |
| bench_512, 4 | 0.0936 | 0.0805 | 0.2825 | 0.2683 |

`build_products` 2.71 -> 1.35 ms per call on the A100 (nsys, bench_256),
i.e. the second read of each field was not served by L2 there either;
the kernel now moves 0.92 GB read + 0.92 GB written in 1.35 ms, 1.36
TB/s, 88% of the bandwidth.  Together the two changes take one A100 from
0.113 to 0.098 s/step at 256^3 (13%) and from 0.957 to 0.820 at 512^3
(14%), four A100 from 0.045 to 0.041 and from 0.306 to 0.268 (12%).

**Overlap of the alltoall** (handoff item 4) was conditional on the
alltoall being the largest single item after item 1.  It is not: on
four A100 the alltoall is 18% of the step at 512^3 and 23% at 256^3,
against 25-30% for the products with their transforms and `buildrhs`,
and on one GPU there is none.  Not done; it stays the next thing to do
for the multi-GPU step (NEXT_SESSION.md).

**Where the time goes now** (one A100, bench_256, nsys kernel summary of
the final build; four A100 from the timer of `jobs/horeka_profile.slurm`):
line solver 27% (three kernels, 2.05 ms for the sweep), `buildrhs` 14%,
cuFFT 25% (four transforms), tiled transpose 8%, `build_products` 8%,
`buildrhs_prepare` 4.5%, `assemble_vvdz` 3.5%; on four GPUs the products
phase 25-30%, the alltoall 18-23%, the implicit solves 13-17%, pack/unpack
12-13%.

**Result** of the session (`jobs/horeka_bench_all.slurm`, same day, all
changes in; the README table; the 4-GPU tests and the 4-GPU field
against the CPU run at 3e-14):

| deck | 1 A100 | 4 A100 | before the session, 1 / 4 |
| --- | --- | --- | --- |
| bench_64 | 0.0145 | | 0.0160 / |
| bench_256 | 0.0986 | 0.0385 | 0.113 / 0.0425 |
| bench_512 | 0.819 | 0.262 | 0.955 / 0.296 |

---

## The line solver (2026-09-26, session 5)

The handoff called the sweep kernel latency-bound: 136 registers per
thread, 19% occupancy, 1.1 GB moved in 2.05 ms at 256^3, "35% of the
A100".  The byte count was wrong.  Counting the three passes of the
old solve (forward: one read and six writes per row; backward: six reads
and three writes; correction: three reads and one write, 16 bytes each)
gives 320 bytes per row, 2.67 GB per call, and Nsight Compute on the A100
(`jobs/horeka_ncu.slurm`, new: the hardware counters are open on the
HoreKA compute nodes, not on the ISTM boxes) confirms it: 1.29 GB read +
1.39 GB written in 2.11 ms, **82% of the DRAM bandwidth** (the gather
kernel before it another 0.27 GB in 0.79 ms at 22%, the scatter 0.29 GB
in 0.41 ms).  The sweep was bandwidth-bound; the lever was bytes, not
threads.

**Two passes instead of five.**  The third pass existed because the
2x2 Schur complement for the border needs the back-substituted values
of the three right-hand sides at rows 0 and 1, the far end of the back
substitution.  But x(0) = e_0^T U^-1 x' is an inner product of the
forward-substituted right-hand side with the first row of the inverse of
the unit upper factor, and that row obeys a *forward* recurrence, w(i) =
-U1(i-1) w(i-1) - U2(i-2) w(i-2) (the same for row 1 with v(1) = 1).  Two
short recurrences and six accumulators in the forward sweep give the
four values the Schur complement needs (rows m-2, m-1 are the last two
rows of the sweep), so the border is known before any back
substitution, and one backward sweep of X - Y1 xb1 - Y2 xb2 is the
solution.  With the rows divided by their pivots as they are stored
(U1, U2, X, Y1, Y2 scaled; the elimination then needs no multipliers),
the solve is: forward sweep reading the right-hand side straight from
the field and writing five columns, backward sweep reading them and
writing the solution straight into the field, 192 bytes per row, one
kernel, no gather, no scatter, no correction pass.  The CPU path is the
same code.  Accuracy unchanged: `test_linsolve` at 2e-15 relative, the
regression decks at 7e-14 against the stored references (which were not
updated), the NCCL build the same.

A100, job 5164142 (`jobs/horeka_ab.slurm`, `~/hst` at the previous
commit against `~/hst-exp`):

| deck, GPUs | implicit solves before | after | step before | after |
| --- | --- | --- | --- | --- |
| bench_256, 1 | 0.0185 | 0.0131 | 0.0980 | 0.0918 |
| bench_256, 4 | 0.0070 | 0.0048 | 0.0413 | 0.0386 |
| bench_512, 1 | 0.1668 | 0.0763 | 0.8228 | 0.6929 |
| bench_512, 4 | 0.0382 | 0.0263 | 0.2691 | 0.2626 |

(The "ghosts, dv/dy, u and w" phase holds a fourth solve per substep and
went from 0.0943 to 0.0539 at 512^3.)  At 512^3 the solver is 2.2x
faster and the step 16%; at 256^3 only 1.3x, and the second Nsight
Compute run (job 5164154) says why: the new kernel moves 1.59 GB per
call (the 192 bytes per row) but at **33% of the bandwidth**, 3.1 ms
under the profiler's locked clocks, 2.39 ms in nsys.  It has 202
registers per thread on cc80 (178 on cc86), so two blocks of 128
threads per SM, 12.5% occupancy, and its 255 blocks make one full wave
and a partial one of 39 blocks; the stalls are waits on global loads.
Halving the bytes has turned a bandwidth-bound kernel into a
latency-bound one at the small grid; at 512^3 (1022 blocks, five waves)
the tail does not matter and the bandwidth is reached.

**Register cap.**  nvfortran gives the kernel 202 registers on cc80 (178
on cc86) for the ten complex values of the two previous rows, the two
recurrences, the six accumulators and the current row.  Capped with
`-gpu=maxregcount` (a per-file flag in the Makefile), three-way A/B on
the A100 (job 5164153):

| deck, GPUs | solves, 202 regs | 168 (3 blocks/SM) | 128 (4 blocks/SM) | step, 202 | 168 | 128 |
| --- | --- | --- | --- | --- | --- | --- |
| bench_256, 1 | 0.0131 | 0.0108 | 0.0117 | 0.0913 | 0.0869 | 0.0880 |
| bench_256, 4 | 0.0049 | 0.0048 | 0.0051 | 0.0362 | 0.0356 | 0.0360 |
| bench_512, 1 | 0.0762 | 0.0767 | 0.0874 | 0.689 | 0.686 | 0.700 |
| bench_512, 4 | 0.0261 | 0.0218 | 0.0243 | 0.2466 | 0.2388 | 0.2412 |

168 (no spills, `LOCAL:0`) takes the sweep kernel from 2.41 to 1.86 ms
per call at 256^3 (nsys) and helps wherever a GPU has about 255 blocks
(one GPU at 256^3, four at 512^3); 128 spills to local memory and loses
at 512^3.  The cap is in the Makefile for `hst_linsolve.o` only.  At
1.86 ms the kernel is at 0.85 TB/s, 55% of the bandwidth, so a third of
its time is still latency at 19% occupancy; what would take it further
is fewer registers by construction (the six accumulators could be four
if rows 0 and 1 were handled as one 2-vector recurrence, or the border
coefficients computed once per line outside the sweep), or two rows of
loads in flight in the backward sweep.  The solver is 12% of the step
now, so this is where the session stopped.

## Overlap of the alltoall with the transforms (2026-09-26, session 5)

Handoff item 2, conditional on the alltoall still being 18% or more of
the 4-GPU step: it was 22% at 256^3 and 19% at 512^3 after the solver
change.  The transforms and transposes now go field by field
(`hst_transforms`): the z transform of field m, then `transpose_start(m)`
(pack and start the alltoall), then the x side of field m-1
(`transpose_finish`: wait, unpack; padding; x transform), so the
alltoall of one field runs while the neighbours are transformed.  Two
send/receive buffer pairs alternate between fields; `start(m+2)` follows
`finish(m)`, which makes the double buffering safe (`hst_mpi`).  NCCL
runs on a second CUDA stream: `ev_packed` (recorded on the compute
stream after the pack) gates the transfer, `ev_done` (recorded on the
communication stream) gates the unpack, and the host never waits.  MPI
posts `MPI_Ialltoall` after `cudaStreamSynchronize` of the compute
stream and waits in `finish`; the CPU build takes the same path.  On one
rank `finish` does the tiled local transpose per field.  cuFFT and FFTW
plans are per field, which costs nothing (same kernels, a third of the
batch, three times: the 1-GPU step is unchanged to 0.5%), and the
results are bit-identical to the batched ones.  The timer's pack/unpack
and alltoall phases are gone: the alltoall has no boundary any more, only
its exposed part costs, and that shows up in the transform and product
phases; nsys shows the NCCL kernels on their own stream.

A100, job 5164159 (`~/hst-exp`, the two-pass solver without the register
cap, against `~/hst-exp4`, the same with the overlap):

| deck, GPUs | step before | after | of which transforms + transposes + products before | after |
| --- | --- | --- | --- | --- |
| bench_256, 1 | 0.0913 | 0.0922 | 0.0620 | 0.0628 |
| bench_256, 4 | 0.0388 | 0.0345 | 0.0282 (alltoall 0.0087) | 0.0238 |
| bench_512, 1 | 0.6907 | 0.6903 | 0.5051 | 0.5047 |
| bench_512, 4 | 0.2530 | 0.2350 | 0.1947 (alltoall 0.0493) | 0.1766 |

The 4-GPU step gains 11% at 256^3 and 7% at 512^3: about half of the
alltoall is hidden at 256^3, 37% at 512^3.  What stays exposed: within a
batch of three fields the first alltoall has nothing before it to hide
behind and the last only the unpack and x transform of field 2, so about
a third of the transfer time is structurally exposed; and the kernels it
overlaps are themselves bandwidth-bound, while NCCL's copies use the
same HBM (the timed phases grew by more than the hidden time would
predict).  A deeper pipeline across batches (the alltoalls of the first
product group behind the products and `buildrhs` of the second) would
address the first part at the price of six products in memory; not done.

**Result** of the session (`jobs/horeka_bench_all.slurm`, job 5164166,
all changes in: the two-pass solver with the register cap and the
overlap; the README table; the 4-GPU test suite and the 4-GPU field
against the CPU run at 2.8e-14, job 5164165):

| deck | 1 A100 | 4 A100 | before the session, 1 / 4 |
| --- | --- | --- | --- |
| bench_64 | 0.0134 | | 0.0145 / |
| bench_256 | 0.0874 | 0.0338 | 0.0986 / 0.0385 |
| bench_512 | 0.687 | 0.227 | 0.819 / 0.262 |

One A100: 11% at 256^3, 16% at 512^3 (the solver); four: 12% and 13%
(the solver and the overlap).  Kernel time on one A100 at 256^3 after
the session: line solver 24%, `buildrhs` 15%, cuFFT 28% (four
transforms), tiled transpose 9%, `build_products` 8.5%,
`buildrhs_prepare` 5%, `assemble_vvdz` 4%.

## The line solver's latency (2026-09-26, session 6)

Handoff item 1: the two-pass sweep kernel at 1.86 ms per call at 256^3
on one A100, 55% of the bandwidth at 19% occupancy, with the register
ideas of the handoff as the candidates.  Reading the kernel for them
turned up something better: **the interior matrix P is real.**  The
stencil coefficients are real, and the wrap phase enters only the two
border rows and, through the wrapped entries of rows 0 and 1, the two
border columns, where it is one common factor conjg(ph).  So the
forward sweep can factor P and forward-substitute the border columns in
real arithmetic (one real division per row instead of a complex one)
and keep complex arithmetic for the right-hand side alone; U1, U2 and
the border columns B1, B2 are stored as real and the phase is applied
to the border values in the back substitution.  The entries of rows
m-2, m-1 that point at the border are no longer moved out of P: they
stay in U1, U2 and act through the seed of the back substitution
(x(n-2), x(n-1)), which also removes the special cases from the sweep.
128 bytes per row instead of 192, a third of the complex flops, and the
kernel fits in 80 registers without spills (166 were used before).
`test_linsolve` at 2e-15, the regressions at 7e-14 (the references
unchanged); the CPU path is the same code.

Nsight Compute on the A100 (`jobs/horeka_ncu.slurm`, bench_256,
launches 10-15 of the kernel; ncu locks the SM clock at 1.09 GHz, so
its durations are 1.3x those of nsys):

| kernel | DRAM bytes per call | duration (ncu) | of the bandwidth | warp cycles per instruction |
| --- | --- | --- | --- | --- |
| two-pass complex, 168 registers (session 5) | 1.59 GB | 2.4 ms | 33% (55% in nsys) | |
| real factor (this session) | 1.08 GB | 1.51-1.59 ms | 44-47% | 20.5 |
| + backward sweep four rows at a time | 1.08 GB | 1.32-1.42 ms | 48-53% | 17.7 |
| + right-hand side of the next row prefetched | 1.08 GB | 1.28-1.41 ms | 48-54% | 16.7 |
| + stencil coefficients of the next row prefetched | 1.08 GB | 1.44 ms | 47% | 15.7 |

The SASS of the real-factor kernel explained the handoff's "two rows of
loads in flight": the compiler unrolls the backward sweep by four but
issues the six loads of a row only after the stores of the previous
one, so each row waits a full memory latency.  Loading four rows into
registers before any of them is used (a parameter `nb` in
`cyclic_penta_solve`, the rows left over one by one) is worth 11% of
the kernel.  The two prefetches after that gained nothing (the second
even costs instructions), so they are not in the code.  A diagnostic
build without the backward sweep put the split at 0.6 ms forward, 0.7
ms backward for the D2V and ETA kinds (ncu clocks), each at 65-70% of
the bandwidth, with L1 at 43%; the KIND_DY solve is 7% slower than the
others because its five-point stencil reads the strided field five
times per row (L1 at 78% in its forward sweep; a sliding window would
recover most of it, about 1% of the step, not done).  The kernel's
stalls are scoreboard waits on L1TEX: the strided reads of the field
(one line per thread per row, 32 sectors per warp request) and the
write-through of the workspace.  What is left is the structure itself,
one thread per line reading a y-innermost field; the next step would
not be local.

A100 timings (`jobs/horeka_ab.slurm`), A = session 5, B = real factor,
C = B + the four-row backward sweep:

| deck, GPUs | step A | B | C | implicit solves A | B | C | ghosts, dv/dy, u and w A | B | C |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| bench_256, 1 | 0.0877 | 0.0843 | 0.0830 | 0.0108 | 0.0084 | 0.0074 | 0.0070 | 0.0057 | 0.0051 |
| bench_256, 4 | 0.0339 | 0.0312 | 0.0302 | 0.0047 | 0.0030 | 0.0024 | 0.0033 | 0.0023 | 0.0020 |
| bench_512, 1 | 0.687 | 0.658 | 0.646 | 0.0770 | 0.0576 | 0.0503 | 0.0493 | 0.0401 | 0.0350 |
| bench_512, 4 | 0.2267 | 0.2185 | 0.2146 | 0.0221 | 0.0169 | 0.0147 | 0.0139 | 0.0111 | 0.0098 |

(job 5164287; the two phases hold the four line solves of a substep.)
The sweep kernel in nsys at 256^3: 1834 us per call, 1444, 1219; it is
1.5x faster than at the start of the session and 14% of the 1-GPU
kernel time instead of 19%.  The step gains 5% at 256^3 on one GPU, 11%
on four, 6% and 5% at 512^3.  The solver is now 15% of the step (the
two phases together) on one GPU and 14% on four.

Register cap: the real-factor kernel needs 80 registers, so the cap of
168 (three blocks per SM) could drop to 96 (five) or 128 (four).  At
256^3 on one GPU and 512^3 on four the grid is 255 blocks, all resident
at once, and the cap cannot matter; at 512^3 on one GPU (1022 blocks)
five blocks per SM make two waves at 95% instead of four at 79%.
Measured (job 5164306, C with the three caps):

| deck, GPUs | solves, cap 168 | 96 | 128 | step, 168 | 96 | 128 |
| --- | --- | --- | --- | --- | --- | --- |
| bench_256, 1 | 0.0073 | 0.0075 | 0.0073 | 0.0827 | 0.0829 | 0.0821 |
| bench_256, 4 | 0.0025 | 0.0025 | 0.0026 | 0.0305 | 0.0306 | 0.0301 |
| bench_512, 1 | 0.0503 | 0.0500 | 0.0492 | 0.646 | 0.648 | 0.646 |
| bench_512, 4 | 0.0146 | 0.0151 | 0.0148 | 0.2149 | 0.2166 | 0.2161 |

All within the run-to-run noise; the kernel at 256^3 in nsys is 1207 us
with 168, 1321 with 96 and 1252 with 128 (the freer schedule wins over
the occupancy, which the 255-block grid cannot use anyway).  The cap
stays at 168.

**Result** of the session (`jobs/horeka_bench_all.slurm`, job 5164344;
the RTX 3060 and CPU columns measured on istmio2):

| deck | 1 A100 | 4 A100 | 1 RTX 3060 | istmio2 CPU, 4 ranks | before the session |
| --- | --- | --- | --- | --- | --- |
| bench_64 | 0.0117 | | 0.109 | 0.70 | 0.0134 / - / 0.127 / 0.89 |
| bench_256 | 0.0817 | 0.0304 | 0.85 | | 0.0874 / 0.0338 / 0.98 / - |
| bench_512 | 0.646 | 0.216 | | | 0.687 / 0.227 |

The cards that run FP64 slowly gain the most (the real division and the
real elimination): 14% on the RTX 3060, 21% on the CPU.

## The second node (2026-09-27, session 6, handoff item 2)

`jobs/horeka_2node.slurm` (new): the same decks on 4 A100 of one node
and on 8 A100 of two nodes (hkn0424 and hkn0812, InfiniBand between
them), NCCL transport and, for reference, MPI; the session-5 build (the
solver difference is in the table above and does not touch the
transposes).  Job 5164343 on `accelerated` (the dev partition has three
nodes and never scheduled it):

| deck | 4 GPUs, 1 node | 8 GPUs, 2 nodes, nccl | 8 GPUs, mpi | of the 8-GPU nccl step: to physical / products |
| --- | --- | --- | --- | --- |
| bench_256 | 0.0345 | 0.0740 | 0.282 | 0.0213 / 0.0443 (4 GPUs: 0.0072 / 0.0170) |
| bench_512 | 0.2260 | 0.5208 | 2.489 | 0.1599 / 0.3333 (4 GPUs: 0.0517 / 0.1246) |

**Eight GPUs on two nodes are 2.1-2.3x slower than four on one.**  The
solve phases halve as they should (0.0083 to 0.0073 at 256^3, 0.0356 to
0.0206 at 512^3), so the whole loss is in the two phases that hold the
transposes: the alltoall of an x-z pencil decomposition sends 7/8 of
every field off the rank, and half of that leaves the node.  The
inter-node share per step at 256^3 is about 1.8 GB per node each way
(27 field alltoalls), and the growth of the two phases, 0.041 s, puts
the effective inter-node rate at about 45 GB/s per node: NCCL is on the
InfiniBand (a socket fallback would be ten times slower), at roughly
half of the node's four HDR links, and even at their full 100 GB/s the
inter-node part alone would be 18 ms of a 34 ms step.  Deeper overlap
cannot hide that (the transforms it would hide behind take 15 ms).
MPI's alltoall over the two nodes is 3.8x slower still than NCCL's.

Conclusion: with x-z pencils the code stops at one node, and the four
A100 of a node are the fastest configuration for these decks.  A second
node needs the y decomposition (WP6, DESIGN.md 7 (i)): each node keeps a
y slab, the x-z alltoalls stay inside the node on NVLink, and only the
line solver's border couplings and the ghost rows cross the node
boundary, a few KB per line instead of half the field.  That is the
only item left that adds real structure, and the next session should
ask before starting it.

**Can the transpose be made to perform better across nodes?** (asked
after the measurement above; job 5167650 on two nodes, hkn0433 and
hkn0514, with the extended `jobs/horeka_2node.slurm`.)  Three things
were tried, none helps, and the node's topology says why.

The topology (`nvidia-smi topo -m` in the job): each A100 node has one
InfiniBand adapter (`mlx5_0`, HDR, 25 GB/s each way) and one RoCE port
(`mlx5_2`), both on NUMA node 0 next to GPUs 0 and 1; GPUs 2 and 3 sit
on the other socket.  NCCL uses the InfiniBand adapter for every
connection (its log names net device 2 for all of them).  The bytes that
must cross: at 256^3 on 8 ranks a block is 3.2 MB, a rank sends four of
its seven blocks to the other node, 51 MB per node per field alltoall,
1.38 GB per node per step.  At 25 GB/s that is 55 ms; the 8-GPU step
exceeds the 4-GPU step by 42 ms, so NCCL is already at the wire rate
with part of it hidden behind the transforms.  Nothing below could
change that.

| bench_256, 20 steps (bench_512, 10 steps) | s/step |
| --- | --- |
| 4 GPUs, one node, nccl | 0.0297 (0.216) |
| 8 GPUs, two nodes, nccl, ranks unbound | 0.0717 (0.527 bound) |
| 8 GPUs, nccl, each rank bound to its GPU's NUMA domain (`jobs/bind_numa.sh`) | 0.0719 |
| 8 GPUs, two-level alltoall `nccl2`, unbound / bound | 0.1196 / 0.1189 (0.897) |
| 8 GPUs, mpi | 0.2896 |
| 8 GPUs, nccl, `NCCL_NET_GDR_LEVEL=SYS` (GPUDirect RDMA forced) | 0.1275 (0.950) |

1. *GPUDirect RDMA.*  `NCCL_DEBUG=INFO` shows the adapter GDR-capable
   ("GPU Direct RDMA (nvidia-peermem) enabled for HCA 0") but every
   connection "via NET/IBext_v9/2/Shared": the transfers go through
   host memory, because the GPUs are farther from the adapter (NODE and
   SYS distance) than NCCL's default `NCCL_NET_GDR_LEVEL` allows.
   Forcing it with `NCCL_NET_GDR_LEVEL=SYS` (job 5167650) makes every
   connection "GDRDMA/Shared" and the step 1.8x *slower* (0.1275 against
   0.0717 at 256^3, 0.950 against 0.527 at 512^3): the RDMA reads of
   GPUs 2 and 3 cross the socket interconnect, which is why NCCL's
   default keeps them off.  The default is right on this machine.
2. *NUMA binding.*  `jobs/bind_numa.sh` runs each rank under `numactl`
   on the cores and memory of its GPU's NUMA node (the wrapper takes
   the local rank, reads the GPU's PCI address and its `numa_node` in
   sysfs).  No effect: 0.0719 against 0.0717.  The host side of NCCL is
   not the limit.
3. *Two-level alltoall* (`transport = 'nccl2'`, commit c151db4, reverted
   in 492469b): level 1 inside the node, each rank sends to local rank
   l' its blocks for every rank with that local index, the ones bound
   for other nodes into a staging buffer of the local relay; level 2
   one message per rank and remote node, which lands in the receive
   buffer where those sources are contiguous, so the unpack is
   unchanged.  About 60 lines with a check that the ranks of a node are
   consecutive.  Correct (regression decks on 8 ranks over two nodes at
   7e-14) and 1.7x *slower*: the flat alltoall already saturates the
   link, its NVLink part overlaps its InfiniBand part, while the two
   levels are sequential and a single message per peer gets fewer NCCL
   channels than four concurrent ones.  Not kept in main.

So with x-z pencils the transpose cannot be made cheaper across nodes
on this machine: the volume is fixed by the decomposition and the wire
is already full.  Only the y decomposition (WP6) changes the volume.

## The y decomposition (2026-09-27, session 7, branch `multinode-y`)

The handout's plan, done in three steps.  The measurements of "The second
node" above say what it is for: with x-z pencils half of every alltoall
crosses the node's single InfiniBand adapter and eight A100 on two nodes
are 2.1-2.3x slower than four on one; only changing what crosses the
node can help, and y is the direction where the coupling is local.

**Phase 0, on `main` (commit 452597d): the physics files written for any
`npy`.**  The contract of DESIGN.md 7 (i): every y loop runs `ny0, nyN`,
ghost rows come only from `exchange_ghost_rows` in `hst_mpi` (on `main`
the wrap of the single slab), box integrals are partial sums plus an
allreduce (the pressure mean was the one place written differently),
the Stokes rows are touched by ownership, the initial field fills the
rank's rows by global index, and the file I/O writes the rows a rank
owns through a row-range MPI-IO view, the file's four ghost rows in two
more collective writes from the ranks that own their sources (a rank's
file rows are then one contiguous block, which is what a subarray view
can describe; reading takes the interior rows and `fill_ghosts` does the
rest).  `npy` is a deck parameter on `main` too (`&mesh`, must be 1
there), so that `hst_params` and `hst_input` are the same on both
branches.  Checked bit for bit against the handout commit on the three
regression decks: GPU one and two ranks exactly zero; CPU exactly zero
once both builds use `FFTW_ESTIMATE` (with the default `FFTW_MEASURE`
the CPU build differs from *itself* run to run by 5e-15, because the
timed plan choice picks different algorithms).

**Phases 1 and 2, on the branch (commit 2678b4a): two files differ from
`main`, `hst_mpi.f90` (+150 lines) and `hst_linsolve.f90` (+230);
`hst_io.f90` needed nothing beyond phase 0.**

- *Decomposition* (`hst_mpi`): `nproc = npxz*npy`, the `npxz` ranks of a
  slab consecutive (`ipy = iproc/npxz`), so `mpirun --map-by
  ppr:npxz:node` puts one slab on one node; `comm_xz` (the transposes,
  MPI or one NCCL communicator per slab, its id broadcast inside the
  slab) and `comm_y` (the y column).  `nyB = ny/npy >= 8`.
- *Ghost rows*: with `npy = 1` the wrap in place as before (so the
  single-slab run has no extra copies); otherwise a pack kernel of the
  two rows for each neighbour, two `MPI_Sendrecv` on the device buffers
  over `comm_y`, and an unpack kernel that applies the shear-periodic
  phase on the two slabs at the box edge.
- *The line solve*: the bordering of session 5 per slab.  The forward
  sweep is the old one over the slab's `nyB - 2` interior rows (the
  entries of its rows 0 and 1 that point at the slab below are the
  border columns, with the phase factored out only on the first slab);
  what changes is that the four rows the border rows reach (0, 1, m-2,
  m-1) are kept as affine functions of *two* borders, `b_{s-1}` and
  `b_s`, instead of one, and written as a record of 28 reals per line:
  four rows of (value, 4 real coefficients) and the right-hand side of
  the slab's two border rows (the handout counted 20 complex; the
  coefficients are real because the interior matrix is, and the border
  right-hand sides must travel too, since every rank assembles every
  block row).  An in-place `MPI_Allgather` over `comm_y` of the records
  (device buffers) between the two kernels; the second kernel assembles
  the `2 npy x 2 npy` block-cyclic system of its line from the `npy`
  records (the coupling of a slab to the slab above carries `ph` on the
  last slab, that to the slab below `conjg(ph)` on the first; the block
  columns are accumulated, since they coincide for `npy = 1` and `2`),
  solves it by Gaussian elimination without pivoting in thread-private
  storage (at most 16 x 16), and does the backward sweep seeded with
  `b_s` and the border columns times `b_{s-1}`.  With `npy = 1` this is
  the old 2 x 2 Schur complement in a different order of operations, so
  the branch at `npy = 1` agrees with `main` to round-off, not bit for
  bit (3e-14 on the regression decks after 50 steps, the same level as
  the CPU/GPU difference).
- *Validation*: the full suite and the three regression decks at 1e-10
  on the CPU with 2 ranks x 2 slabs, 2 pencils x 2 slabs and 1 x 4
  slabs (`nyB = 8` on `small`), on the RTX 3060 with 2 ranks x 2 slabs,
  on istmcetus (2 x A6000, NCCL build) with 1 x 2 and 2 x 1; the line
  solver's unit test at 1e-15 in every layout; `test_roundtrip`
  exercises the row-range I/O.  `tests/run_tests.sh` and
  `tests/regression.sh` take `npy` as a third argument.

**What the exchange costs where it could be measured.**  bench_256 on
the two A6000 of istmcetus (PCIe, CUDA-aware MPI, one slab per GPU):
of a 0.170 s step, the ghost rows take 2.1 ms and the allgathers 9.4 ms
of pure transfer (9 solves x 7.3 MB at 7 GB/s), 7% of the step; on
HoreKA the same bytes go over InfiniBand at 25 GB/s per node, shared by
the four ranks of the node.  The estimate for two A100 nodes at 256^3
is therefore about 3 ms of y traffic against 15 ms of compute and the
node-internal alltoall, i.e. a two-node step around 20 ms against 30 ms
on one node (1.5x), and 512^3 similarly; to be measured
(`jobs/horeka_2node.slurm` with `npy` in `CONFIGS`, `NPY_REG=2`).
If the allgather turns out exposed and large, the next levers are (a)
`ncclAllGather` on `comm_y` in place of MPI, (b) overlapping the
allgather of one `line_chunk` batch with the forward sweep of the next
(the batches exist; on the GPU the default is one batch), (c) halving
the record by caching the 16 matrix coefficients per line for the
system kinds whose matrix does not change between calls (`KIND_D0`,
`KIND_DY`; the implicit systems change with `lambda`).

**HoreKA, one A100 node (jobs 5167735 and 5167736, 2026-09-27
evening; the dev partition had nothing running for three hours before
they started).**  The suite with 2 pencils x 2 slabs (NCCL for the
alltoall inside each slab): all 12 tests pass, the 4-GPU `small` deck
agrees with the login node's CPU run at 3e-14.  The A/B of `main`
against the branch at `npy = 1` (`jobs/horeka_ab.slurm`, same node):

| s/step | main | branch, npy = 1 | solve phases main -> branch |
| --- | --- | --- | --- |
| bench_256, 1 A100 | 0.0827 | 0.0833 | implicit 7.3 -> 7.7 ms, ghosts/dv/dy 5.1 -> 5.6 |
| bench_256, 4 A100 | 0.0303 | 0.0304 | 2.45 -> 2.52, 1.96 -> 2.00 |
| bench_512, 1 A100 | 0.645 | 0.654 | 50.4 -> 55.9, 35.1 -> 38.7 |
| bench_512, 4 A100 | 0.215 | 0.217 | 14.6 -> 15.4, 9.8 -> 10.7 |

A wash at the step (0.5-1.3%), and the solve phases show what the
split costs with one slab: 6-11%, i.e. the record written per line, the
second launch and the assembly of a 2 x 2 system through the general
16 x 16 path in thread-private memory.  Not worth a special case; the
lever, if ever needed, is a `nslab == 1` branch in `penta_backward`
that reuses the old 2 x 2 Schur code (about 20 lines).

**nvfortran 25.9 and names in OpenMP clauses.**  Once `hst_mpi` (a CUDA
Fortran module) is visible in a file, even through a `use ..., only:`
chain, nvfortran rejects the names `kind` and `x` in the data-sharing
clauses of that file's target regions ("must appear in a SHARED or
PRIVATE clause", although they do).  On `main` the module-level `use
hst_mpi` in `hst_derivatives` triggered it in `hst_linsolve`; scoping
the `use` to `fill_ghosts_field` cured it.  On the branch `hst_linsolve`
uses `hst_mpi` itself, and there the dummy `kind` is called `sys` and
the workspace `X` is `XR`.

**Two A100 nodes (job 5167737, 2026-09-28 02:30 after 6.5 h in the
`accelerated` queue; hkn0701 and hkn0721).**  `jobs/horeka_2node.slurm`
with `npy` in `CONFIGS`, one slab per node (`--map-by ppr:4:node`, so
the four pencils of a slab share a node and every alltoall stays on
NVLink); the regression decks on 8 ranks with `npy = 2` pass at 5e-14
(5e-13 the Stokes deck, as everywhere on the GPU).

| s/step | 4 GPUs, 1 node, 4 x 1 | 4 GPUs, 1 node, 2 x 2 | 8 GPUs, 2 nodes, 4 x 2 | 8 GPUs, 2 nodes, 8 x 1 (x-z pencils) |
| --- | --- | --- | --- | --- |
| bench_256 | 0.0306 | 0.0370 | 0.0304 | 0.0714 |
| bench_512 | 0.2183 | 0.2423 | 0.1626 | 0.521 (session 6) |

So the decomposition does what it was built for on the transposes and
then loses most of it on its own exchange.  The two alltoall phases
halve as they should (256^3: 24.1 ms on one node, 13.5 on two; 512^3:
178 and 90 ms), and the solve phases, whose compute halves too, grow
instead (256^3: 4.5 to 15.9 ms; 512^3: 26 to 66 ms).  The timer's
"transfers only" line says where: the y exchange (the host's wait in
`MPI_Sendrecv` and `MPI_Allgather` on device buffers over `comm_y`,
after the device has finished the pack) is

| transfers per step | ghost rows | reduced systems | sum | of the step |
| --- | --- | --- | --- | --- |
| bench_256, 2 nodes | 7.0 ms | 6.7 ms | 13.7 ms | 45% |
| bench_512, 2 nodes | 24.9 ms | 29.2 ms | 54.1 ms | 33% |
| bench_256, one node as 2 x 2 (NVLink) | 0.95 ms | 1.4 ms | 2.3 ms | 6% |
| bench_512, one node as 2 x 2 | 1.5 ms | 2.9 ms | 4.4 ms | 2% |

Without it the two-node step would be 17 ms at 256^3 (1.8x one node)
and 109 ms at 512^3 (2.0x).  The bytes do not explain it.  Per rank a
ghost exchange sends and receives two rows of its 8160 lines at 256^3
(2 x 261 KB) or 32704 lines at 512^3 (2 x 1.05 MB), 12 times a step
(four fields per substep); a record is 28 reals per line, 1.8 MB at
256^3 and 7.3 MB at 512^3, allgathered 9 times a step (three solves per
substep; the D0 solve of `exact_shift` is off in the benchmarks).  Per
node that is 25 + 66 MB each way per step at 256^3 and 100 + 264 MB at
512^3, i.e. 3.6 ms and 14.5 ms at the 25 GB/s of the node's InfiniBand
adapter, against the 13.7 and 54 ms measured: HPC-X's MPI on device
buffers runs the y exchange at 5-7 GB/s per node, the same 4x below
the wire that "The second node" found for its alltoall (and the timer
also charges the wait for the partner slab, which the MPI path
serialises with the host).  Inside a node the same calls run at
NVLink rate (2.3 ms), so it is the inter-node path of MPI, not the
decomposition.

Two more numbers from the same job: one node as 2 pencils x 2 slabs is
21% slower than as 4 pencils at 256^3 and 11% at 512^3 (the y exchange
plus longer alltoall phases with two ranks per slab), so on one node the
pencils stay; and the 8 x 1 x-z-pencil run on two nodes reproduces the
session-6 number (0.0714 against 0.0717).

**NCCL on the y column (branch commit 54a3d2a; A/B job 5168032, the
branch's previous commit in `~/hst-y` against the variant in
`~/hst-exp`, same two nodes, `ROOTS` in `jobs/horeka_2node.slurm`).**
Since session 6 had already measured MPI's inter-node device path at a
quarter of NCCL's, the lever is the transport, not the volume: a second
NCCL communicator over `comm_y` (`use_nccl_y`, independent of the slab's
`use_nccl`, because one pencil per slab has no alltoall), the ghost
rows as two grouped send/receive pairs and the records through
`ncclAllGather` (in place, the slab's block as the send buffer), both
issued on the compute stream between the pack and the unpack, so the
host does not wait at all (with `timing` it does, before and after, to
time the transfer alone); MPI on the device buffers stays as the
fallback.  About 80 lines in `hst_mpi`, nothing elsewhere.  On the two
A6000 of istmcetus (PCIe, one node) the change is a wash, as it should
be: ghost rows 3.1 ms, records 8.3 ms (MPI 2.1 and 9.4).

| s/step, same two nodes (hkn0426, hkn0428) | branch, MPI y exchange | with NCCL on `comm_y` |
| --- | --- | --- |
| bench_256, 1 node, 4 x 1 | 0.0301 | 0.0306 |
| bench_256, 1 node, 2 x 2 | 0.0370 | 0.0375 |
| bench_256, 2 nodes, 4 x 2 | 0.0302 | **0.0235** |
| bench_256, 2 nodes, 8 x 1 | 0.0704 | 0.0703 |
| bench_512, 1 node, 4 x 1 | 0.2179 | 0.2169 |
| bench_512, 1 node, 2 x 2 | 0.2424 | 0.2440 |
| bench_512, 2 nodes, 4 x 2 | 0.1688 | **0.1265** |

The y exchange on two nodes goes from 5.9 + 7.4 ms (ghost rows +
records) to 5.2 + 5.5 ms at 256^3 and from 25.1 + 34.9 to 8.8 + 14.9 ms
at 512^3; everything on one node, `npy = 1` included, is a wash (the
NCCL communicator over the column is created only with `npy > 1`).  Two
A100 nodes are now 1.28x one node at 256^3 and 1.72x at 512^3, with the
MPI y exchange they were 1.0x and 1.3x.  Kept on the branch.

**What is left in the exchange.**  At 512^3 the records now move at 70%
of the wire (14.9 ms against 10.5) and the ghost rows at half (8.8
against 4.0); at 256^3 the 21 exchanges of a step cost about 0.5 ms
each whatever their size (12 ghost exchanges of 2 x 261 KB in 5.2 ms,
9 allgathers of 1.8 MB in 5.5 ms), i.e. a fixed cost per exchange of
two nodes meeting through NCCL's host-proxied InfiniBand path, plus
whatever skew the timer on rank 0 sees while the other slab arrives.
So 10.7 ms of the 23.5 ms step at 256^3 and 23.7 of 126.5 ms at 512^3
are still the y exchange, and what can be done about it, in the order
of expected gain per line of code, all to be measured with `ROOTS` in
`jobs/horeka_2node.slurm`:

1. One NCCL group for the ghost rows instead of two (NCCL takes several
   sends to the same peer in a group, in order): one kernel per
   exchange instead of two, worth up to half of the ghost-row time at
   256^3 if the cost is per kernel.  Five lines.
2. The allgather of one `line_chunk` batch overlapped with the forward
   sweep of the next (`chunk` is all the lines on the GPU today; two
   record buffers, NCCL on a second stream with events as the alltoall
   does): hides the bandwidth part of the records, about 10 ms of the
   14.9 at 512^3, nothing of the fixed cost, so little at 256^3.
3. A smaller record for the kinds whose matrix does not change between
   calls (the 16 coefficients cached per line; `KIND_DY` every substep,
   the implicit kinds too while `deltat` is fixed): fewer bytes, worth
   something only at 512^3.

The ghost rows' fixed cost is out of reach without computing the
interior rows while the exchange is in flight, which is structure in
the physics files that `main` does not want.

## The layout without user input (2026-09-28, session 9)

**The two branches stay separate and must run the same decks and field
files.**  `tests/crossbranch.sh <build A> <build B> [nranks] [npy A]
[npy B]` (on `main`) runs 25 steps of the `small` deck with one build and
the other 25 with the other, restarted from the first one's
`Dati.cart.out`, and compares the result with the 50-step reference at
1e-10.  Between `main` and the branch at `npy = 2`, in both directions,
the round trip agrees to 2e-14 on the CPU (4 ranks: 2 x 2 on the branch)
and 6e-14 on the two A6000 of istmcetus (1 pencil x 2 slabs through
NCCL on the y column), i.e. at the regression's own level: the file
written by either branch (interior rows plus the four ghost rows, each
slab writing its own on the branch) restarts the other exactly.

**`npy = 0`, the default, means "the code chooses".**  On `main` it is
one slab (one line in `setup_decomposition`); on the branch
`setup_decomposition` counts the ranks per node
(`MPI_Comm_split_type(MPI_COMM_TYPE_SHARED)`, once, also for
`setup_transport`) and takes one slab per node when every node holds the
same number of ranks, consecutive in `MPI_COMM_WORLD` (what `mpirun
--map-by ppr:N:node` and Slurm's block distribution give), the number
of nodes is at most 8 (`NPY_MAX`, the reduced system's size) and divides
`ny` with at least 8 rows per slab, and the pencils per node divide
`nx+1` and `nzd`; otherwise one slab, with a one-line notice when there
is more than one node.  An explicit `npy > 1` on ranks that are not
node-consecutive aborts with a message naming the two launch options,
since such a layout would put a slab's alltoall across the nodes.  About
30 lines in `hst_mpi`, nothing in the physics files or the job scripts.
Checked on one node (istmcetus, two GPUs: 2 pencils x 1 slab, the tests
and regressions at 1, 2 and 4 slabs unchanged on CPU, GPU and NCCL) and
on two real nodes with the CPU build across istmio2 and istmcetus
(system OpenMPI, 2 ranks each): block mapping and `npy = 0` give 2 x 2,
cyclic mapping (`--map-by node`) with `npy = 2` aborts, cyclic with
`npy = 0` prints the notice and runs 4 x 1, block with `npy = 2` runs
2 x 2.  The HoreKA two-node job with `npy = 0` in `CONFIGS` (job
5168929) is the check that 8 ranks on two A100 nodes come out 4 x 2 and
4 ranks 4 x 1.

**On two HoreKA A100 nodes (job 5168929, 2026-09-29 night; the branch
head in `~/hst-exp`, `npy = 0` in every deck).**  8 ranks with `--map-by
ppr:4:node` come out as "ranks = 8 (4 x-z pencils x 2 y slabs)", 4 ranks
as 4 x 1; the three regression decks on 8 ranks at `npy = 0` pass at
5.9e-14 (5.7e-13 the Stokes deck); the bench steps equal those with
`npy = 2` (0.0232 and 0.1268 s against 0.0233 and 0.1274 in job
5168845 the same night).  So one deck now runs unchanged on one GPU, on
the four GPUs of a node and on two nodes, on either branch.

## The fixed cost of the y exchange (2026-09-29, session 9, job 5168845)

The `small` deck (nx = nz = 15, ny = 32) on 8 ranks with two slabs has
exchanges that carry almost nothing (a ghost exchange 2 rows x 31 x 4
complex = 4 KB per direction per rank, a record gather 28 x 124 reals =
28 KB per rank), so its "transfers only" line is the per-exchange
latency plus the skew between the two slabs.  Same job, `~/hst-y`
(the branch before `npy = 0`), 50 steps:

| `small` | s/step | ghost rows / step | records / step | per exchange (12 ghosts, 9 gathers) |
| --- | --- | --- | --- | --- |
| 8 A100, 2 nodes, 4 x 2 | 0.00674 | 2.49 ms | 2.32 ms | 0.21 ms, 0.26 ms |
| 4 A100, 1 node, 2 x 2 (NVLink) | 0.00375 | 0.57 ms | 0.51 ms | 0.05 ms, 0.06 ms |
| 4 A100, 1 node, 4 x 1 | 0.00328 | - | - | - |

So an exchange between the nodes costs 0.2-0.26 ms before any bytes,
five times the NVLink value; at 256^3 the same job's exchanges cost
0.41 ms (ghost rows, 4.95 ms / 12) and 0.61 ms (records, 5.51 / 9), so
about half of the 10.5 ms of y exchange in the 23.3 ms step is this
fixed part and half scales with the size (bytes at the effective rate,
and skew that grows with the work), and at 512^3 (9.3 + 15.0 ms) the
fixed part is a fifth.  The fixed part is out of reach of the levers
below (fewer exchanges would need the two-field ghost exchange, a
change in the physics files); the size-dependent part is what (a) and
(c) address.

## Lever (a): one NCCL group for the ghost rows (job 5168929): dropped

The two send/receive pairs of `exchange_ghost_rows` (to the slab above
and to the slab below; with two slabs the same rank) in one
`ncclGroupStart/End` instead of two.  Same two nodes, back to back:

| 8 A100, 2 nodes, `npy = 0` | branch head | one group |
| --- | --- | --- |
| bench_256 s/step | 0.02323 | 0.02369 |
| its ghost rows per step | 4.88 ms | 5.13 ms |
| bench_512 s/step | 0.12683 | 0.12740 |
| its ghost rows per step | 8.92 ms | 10.93 ms |
| `small` s/step | 0.00779 | 0.00786 |

A wash at 256^3, worse at 512^3 (the one-node rows within 1.5%); the
cost of a ghost exchange is not per NCCL kernel.  Dropped (the
worktree and the branch `lever-a` deleted).

## Lever (c): the records' gather behind the next batch's sweep (job 5168954): kept

`allgather_y` became `allgather_y_start` and `allgather_y_wait`
(`hst_mpi`): the gather of one solver batch's records is issued on the
communication stream after an event recorded on the compute stream
(`ev_rec`), and the compute stream waits for its completion event
(`ev_gath`) only before that batch's backward sweep, so the forward
sweep of the next batch runs while the records travel and the host
never blocks.  A workspace's blocks are not contiguous over the slabs,
so the gather is one send and one receive per other slab in a group
(NCCL) or a shift loop of `MPI_Sendrecv` (MPI; a first version paired
each rank with `ipy + s` on both sides and deadlocked at four slabs).
`hst_linsolve` keeps two workspaces when `npy > 1` (the line index
offset by `nlines_max`, so the device routines are untouched) and
`line_solve` runs the batches as a two-deep pipeline (forward sweep of
batch b+1 before the backward sweep of b; with one workspace the loop
degenerates to forward, gather, backward).  The number of batches
stays `line_chunk`'s, and the default is two batches when each still
fills the GPU (`LINES_FULL` = 108 SMs x 128 threads, one thread per
line: 512^3 on four pencils has 32704 lines per rank, 256^3 8160), else
one.  The `timing` line now reports the *exposed* part of the gather.
Same two nodes, back to back, `lc<N>` = `line_chunk = N` (16 or 32 =
two batches, 16 at 512^3 four):

| 8 A100, 2 nodes, `npy = 0` | branch head, 1 batch | head, 2 batches | lever (c), 1 batch | lever (c), 2 batches | lever (c), 4 batches |
| --- | --- | --- | --- | --- | --- |
| bench_256 s/step | 0.02347 | 0.02574 | 0.02317 | 0.02356 | |
| records per step | 6.05 ms | 6.40 | 3.36 (exposed) | 1.59 (exposed) | |
| implicit solves per step | 4.24 ms | 5.86 | 3.83 | 4.13 | |
| bench_512 s/step | 0.12702 | 0.12775 | 0.12727 | **0.12181** | 0.12250 |
| records per step | 15.58 ms | 14.21 | 11.72 (exposed) | 4.61 (exposed) | 2.95 (exposed) |
| implicit solves per step | 15.97 ms | 16.49 | 15.89 | 12.46 | 13.16 |

At 512^3 two batches hide 11 of the 15.6 ms of records and the step
gains 4% (0.1270 to 0.1218 s; the 8-GPU 512^3 rows of five runs across
three jobs that night spread 0.1268-0.1274, so the gain is well outside
the noise); the one-node rows are a wash (0.2159 / 0.2157 s at 512^3,
0.0300 / 0.0305 at 256^3, within the 2% the 4-GPU 256^3 row moves
between jobs).  At 256^3 the halved batch costs more in the sweeps
(4080 lines, 32 blocks of 128 threads on 108 SMs) than the overlap
hides, so one batch stays there, and the default rule says so.  Two
A100 nodes are now 1.77x one node at 512^3 (1.70 in the same job
before) and 1.29x at 256^3.  Of the 121.8 ms step at 512^3 the y
exchange still shows 8.6 ms of ghost rows and 4.6 ms of exposed records
(11%); at 256^3 4.7 + 1.6 of 23.6 ms (27%), mostly the fixed cost
measured above.  Kept: merged into `multinode-y`.

## Lever (b): one ghost exchange for u and w (2026-09-29, session 10, job 5170040): kept

The recovery of u and w at the end of `linsolve` ended with
`fill_ghosts(1)` and `fill_ghosts(3)`, two ghost-row exchanges back to
back, each paying the inter-node fixed cost.  Now `fill_ghosts(c, c2)`
takes an optional second component and `exchange_ghost_rows(field,
shift_x, shift_z, field2)` an optional second field: on `main` the
shear-periodic wrap runs once per field (an internal subroutine, since
an absent optional must not appear in a target region: `present` is
tested on the host and the kernel launched per field); on the branch the
buffers `ghost_send`/`ghost_recv` became `(2 rows, z, x, field,
direction)`, pack and unpack are internal subroutines launched once per
present field, and one NCCL (or MPI) message per direction carries both
fields, the field index being inside the direction one.  `linsolve` ends
with `call fill_ghosts(1, 3)`: 9 instead of 12 ghost exchanges per step
(18 instead of 21 y exchanges).  The first change to a physics file made
for the parallel layer (DESIGN.md, departures), ten lines in
`hst_derivatives` and `hst_equations`; the same operations in the same
order, so the results are bit-identical: `cmp` of the `small` deck's
`Dati.cart.out` after 50 steps before and after, on `main` (RTX 3060, 1
and 2 ranks) and on the branch at 1 x 2 through MPI on the device
(istmio2) and through NCCL (istmcetus), plus the safety net on both
branches and the cross-branch round trips.  Same two A100 nodes, back
to back, `npy = 0` (job 5170040, `dev_accelerated`):

| A100, `npy = 0` | branch head | lever (b) |
| --- | --- | --- |
| bench_256, 4 GPU, s/step | 0.03069 | 0.03015 |
| bench_256, 8 GPU (2 nodes), s/step | 0.02307 | **0.02262** |
| its phase "ghosts, dv/dy, u and w" | 4.50 ms | 4.07 ms |
| its ghost rows / exposed records (transfers only) | 4.64 / 3.32 ms | 4.93 / 3.42 ms |
| bench_512, 4 GPU, s/step | 0.21622 | 0.21676 |
| bench_512, 8 GPU (2 nodes), s/step | 0.12178 | **0.12118** |
| its phase "ghosts, dv/dy, u and w" | 11.71 ms | 11.25 ms |
| its ghost rows / exposed records (transfers only) | 8.96 / 4.62 ms | 8.65 / 4.11 ms |
| `small`, 8 GPU (2 nodes), s/step | 0.00784 | **0.00703** |
| its ghost rows / exposed records (transfers only) | 3.08 / 1.88 ms | 2.24 / 1.68 ms |

The `small` deck shows the fixed cost saved directly: its ghost rows
drop by 27%, the predicted quarter (3 of 12 exchanges), 0.8 ms of its
7.8 ms step.  At 256^3 the 8-GPU step gains 2.0% (the phase that holds
the exchange 0.43 ms, three exchanges' fixed cost; the row's noise
between jobs is about 1%), at 512^3 0.5% (0.46 ms in the phase; noise
0.3%), the 4-GPU rows are a wash (-1.8% and +0.3%, within the 2% the
one-node 256^3 row moves between jobs).  The "transfers only" ghost-row
line at 256^3 went up while the phase went down: that line is the sum
of synchronised transfer times, i.e. transfers plus the skew between the
slabs, and with fewer synchronisation points the skew lands elsewhere;
the phase and the step are the numbers to believe.  Kept.  Two nodes are
now 1.33x one node at 256^3 and 1.79x at 512^3 (this job, where the
one-node rows are 0.03015 and 0.21676).

What is left in the y exchange after (a), (b), (c): at 512^3 8.65 ms of
ghost rows and 4.11 ms of exposed records of the 121.2 ms step (10.5%);
at 256^3 4.93 + 3.42 of 22.6 ms (the exposed records were 1.6 ms in job
5168954 and 3.3-3.4 ms in this one: that part moves between jobs).  The
cheap levers are used up: (d), smaller records for the system kinds
whose matrix does not change between calls, is bytes only (at most 3%
at 512^3), and the ghost rows' fixed cost can only be hidden by
computing interior rows while the exchange is in flight, which is a
change in the physics files' loop structure, not in the parallel layer.

## Scaling: 512^3 and 1024^3 on two and four A100 nodes (2026-09-29, session 11, jobs 5170167 and 5170168)

Measurement only, no code change: `jobs/horeka_2node.slurm` gained a
`BUILD` variable (the H100 build directory) and a `mem` marker that
samples `nvidia-smi` on the first node during a run and prints the peak
memory per GPU, and `examples/bench_1024.in` is `bench_512.in` with the
modes doubled.  Two-node job 5170167 (`dev_accelerated`, hkn0401-0402,
6 min for six runs and the regression decks; its twin on `accelerated`
cancelled), `npy = 0` for the 8-GPU rows:

| A100 40 GB, NCCL, s/step | 4 GPU (1 node) | 8 GPU (2 nodes) | ratio |
| --- | --- | --- | --- |
| bench_256 | 0.03033 | 0.02257 | 1.34x |
| bench_512 | 0.21613 | 0.12089 | 1.79x |
| bench_1024 | out of memory | 1.00537 | |
| bench_1024, peak memory per GPU | (36.7 GB when the allocation failed) | 39.3 GB | |
| bench_1024, phase "ghosts, dv/dy, u and w" | | 50.3 ms (5.0%) | |
| bench_1024, ghost rows / exposed records (transfers only) | | 21.0 / 3.8 ms (2.5%) | |
| bench_512, the same | | 8.6 / 4.5 ms (10.8%) | |

**1024^3 fits on two A100 nodes and nowhere smaller.**  The 8-rank run
peaks at 39.3 GB of the 40 GB (cuFFT work area 2322 MB, two line-solver
workspaces of 1562 MB, the rest the seven spectral components with ghost
rows, 1.08 GB each, and the physical-space and transpose buffers); the
4-rank run dies in `cuMemAlloc` at 36.7 GB while still allocating.  So
the one-node 1024^3 reference does not exist on the A100, and the deck's
first useful configuration is 8 GPUs; a 1024^3 production run has 1.7
GB of headroom per GPU, enough for the statistics and the I/O buffer
(`hst_io` allocates a host copy only).

**The step scales with the points, the exchange does not.**  1024^3 on 8
GPUs is 8.3x the 512^3 step on the same 8 GPUs for 8x the points (the
FFT phases 9.4x and 9.2x, the solves 4.9x, the ghost phase 4.4x): the
y exchange, 13.1 ms at 512^3 (10.8% of the step), is 24.9 ms at 1024^3
(2.5%): the ghost rows carry 4x the bytes (21.0 vs 8.6 ms, i.e. the
fixed cost is no longer what they cost), the exposed records the same
3.8 vs 4.5 ms (the gather hides behind a sweep that is now 4x longer).
So at the production size the two-node code is within 3% of what a
node-local exchange would give, and the levers left in the handoff of
session 10 (interior/boundary row splitting, lever (d)) are worth at
most that.  The 512^3 two-node ratio is unchanged since job 5170040
(1.79x; the phases within 0.5%).

**H100 (`accelerated-h100`, `GPU_ARCH=cc90`, `~/hst-y/build-h100`
built):** all 20 nodes of the partition (hkn0902-0922) are in the
reservation `hk2teal` until 2026-10-31 and the job (5170169) pends with
`ReqNodeNotAvail`; a submission with `--reservation=hk2teal` was
accepted by sbatch and pended the same way (cancelled).  Job 5170169 is
left in the queue; if the nodes come back it runs the two-node configs
with the H100 build into `~/hst-runs/scal-h2`.  Whether the 94 GB H100
holds 1024^3 on one node (estimated 75 GB per GPU) is part of what it
would measure.

**Four nodes (job 5170168, `accelerated`, hkn0701/0706/0707/0711, 7.8 min
in 15, `npy = 0` = 4 x-z pencils x 4 y slabs):**

| A100 40 GB, NCCL, s/step | 8 GPU (2 nodes) | 16 GPU (4 nodes) | ratio |
| --- | --- | --- | --- |
| bench_256 | 0.02257 (job 5170167) | 0.02221 | 1.02x |
| bench_512 | 0.12108 | 0.09534 | 1.27x |
| its FFT phases (to physical + products) | 91.0 ms | 47.2 ms | 1.93x |
| its "implicit solves" phase | 12.0 ms | 26.4 ms | |
| its ghost rows / exposed records (transfers only) | 8.2 / 3.9 ms | 10.2 / 30.3 ms | |
| bench_1024 | 1.00659 | 0.63023 | 1.60x |
| its FFT phases | 842 ms | 425 ms | 1.98x |
| its "implicit solves" phase | 59.7 ms | 107.3 ms | |
| its "ghosts, dv/dy, u and w" phase | 49.1 ms | 70.5 ms | |
| its ghost rows / exposed records (transfers only) | 20.2 / 4.6 ms | 23.9 / 114.2 ms | |
| bench_1024, peak memory per GPU | 39.3 GB | 20.8 GB | |

**The compute halves, the reduced systems do not.**  From two to four
nodes the FFT phases (transforms, node-local alltoalls, products) scale
perfectly, 1.93x and 1.98x, and the ghost rows cost about what they did
(the bytes per rank are the same, one more hop is not).  What breaks is
the gather of the reduced systems: the line solver's "exposed records"
go from 4.6 to 114 ms at 1024^3 (from 3.9 to 30 ms at 512^3), and the
"implicit solves" phase, whose compute should have halved, nearly
doubles.  Each slab gathers the records of all `npy` slabs over the y
column (`allgather_y`), so the bytes per rank double from `npy = 2` to 4
while the sweep they hide behind halves with `nyB`, and at four slabs
the gather crosses three node links per rank instead of one.  Without
the exposed records the 16-GPU steps would be 65 ms (1.86x) and 516 ms
(1.95x): the reduced-system gather is the whole loss.  At 256^3 the
step is flat (22.2 vs 22.6 ms): the exposed records (9.8 ms) and ghost
rows (4.8 ms) are two thirds of it.

So the picture for production: two nodes are the sweet spot per GPU
(1024^3 at 1.0 s/step, the exchange at 2.5%); four nodes buy 1.6x for
2x the GPUs at 1024^3, and the lever that would make them pay is not
the ghost rows any more but the records: (d) halving the bytes for the
system kinds whose matrix does not change (rated <= 3% at two nodes, it
is worth up to 18% of the four-node 1024^3 step), or making the gathered
bytes independent of `npy` (each rank solving the reduced systems of a
share of the lines after an alltoall of the records within the y
column, instead of every rank gathering all records and solving all
lines redundantly).  The regression decks of the job (8 ranks, `npy =
1`) agree with the references at 7e-14 as before.

## The merge of `multinode-y` into `main` (2026-09-29, session 11)

With the scaling numbers in (the y exchange is 2.5% of the 1024^3 step)
the user decided to merge: `xz-parallel` tags `main` before the merge
(x-z pencils only, 7422b43), `xyz-parallel` the merge commit (61c5ab1,
`git merge --no-ff multinode-y`, 648 insertions in `hst_mpi.f90` and
`hst_linsolve.f90`, nothing else).  `main` now runs `npy` slabs, `npy =
0` choosing one per node; `tests/crossbranch.sh` with the same build
twice is the restart test between `npy` values.  Safety net on the
merged `main`, all green: CPU (istmio2, gfortran) `run_tests.sh` at 2,
4 x 2 and 4 x 4 ranks, `regression.sh` at 2 and 4 x 2 (3e-14),
`crossbranch.sh build-cpu build-cpu 4 1 2` and `4 2 1` (2e-14); GPU
(RTX 3060) `run_tests.sh` at 2 and 2 x 2, `regression.sh` at 1, 2 and
2 x 2 (6e-14); NCCL (istmcetus, two A6000) `run_tests.sh` at 2 x 2 and
2, `regression.sh` at 2 x 2 and 2 (7e-14), the restart round trips 1 ->
2 and 2 -> 1 (6e-14).  HoreKA `~/hst` is the merged `main`, rebuilt
(`GPU_ARCH=cc80 NCCL=1`); `jobs/horeka_tests.slurm` with `NPY=2` on
four A100 (job 5170308): all passed, the `small` deck at 2 x 2 against
the login node's CPU run at 3e-14.  The branch `multinode-y` was then deleted
(locally and on origin) and the `-main` worktree removed; its HoreKA copies (`~/hst-y`, `~/hst-exp`,
`~/hst-exp2`) are now the code of `main` minus the docs.


## Production run (2026-09-30, session 12, handoff item 1)

**The H100 will not come.**  The user cancelled job 5170169: the
`accelerated-h100` partition stays in its reservation for good, HoreKA is
being migrated to a new machine.  `~/hst-y/build-h100` is dead weight.

**What a production run needs that the benchmarks did not.**  Two small
additions to `hst.f90` (safety net green on CPU, GPU and NCCL, the `npy`
round trips included): `wall_max` in `&time_control` ends the run after
that many seconds of wall-clock time, decided on rank 0 and broadcast
every step (the ranks' clocks differ, and a loop exit that one rank takes
alone is a deadlock), so that a SLURM segment writes its restart file and
the next job continues; and every snapshot and restart write prints its
duration.  `jobs/horeka_prod.slurm` runs a deck in segments: it copies
the deck to the run directory with `wall_max` = the job's time limit
minus `MARGIN` (600 s) and `time_from_restart = .true.`, and resubmits
itself (same partition, nodes and limit) while the log says the
wall-clock limit ended the run, up to `MAXSEG` segments.  The output of
the big run goes to a workspace, `/hkfs/work/workspace/scratch/xt8786-hst`
(`ws_allocate hst 60`, 2026-09-30).

**How long the 1024^3 restart file takes.**  Job 5170168's 8-rank run
(`~/hst-runs/scal-a4/2node-bench_1024-np8-npy0`) wrote its 25.8 GB
`Dati.cart.out` between the last `Runtimedata` line (04:21:02) and the
file's mtime (04:23:28): about 140 s, 185 MB/s, on the home file system
(not Lustre: `lfs getstripe` refuses it).  Each rank's subarray view is
made of 24 KB pieces (`ncomp * nyB` complexes), which is what collective
MPI-IO gets from that; a production run at that size with 10 velocity
snapshots, 10 pressure files and 10 restarts would spend about 1 h of a
30 h run writing.  Tolerable for the first run; the levers, if it
matters, are ROMIO's collective-buffering hints and the workspace file
system, both untested.

**The CPL post-processing on our files.**  `jobs/cpl_postprocess.sh`
writes `scddns.in` (CPL names: their `ny` = our `nz`, their `nz` = our
`ny + 1`) and `postpro.in` next to `hst.in`, builds `postpro.cpl` with
`mpicpl` from a copy of `hst-main` and runs it.  The copy needs one
patch: `hst-main/postprocess/convenience.cpl` and
`pressure_reconstruction/poisson.cpl` call `penta_smw_solve` with a fifth
argument `check_linsolve` that `linsolver_smw.cpl` no longer takes (the
CPL compiler says "function expected"), i.e. the post-processing of
`hst-main` is out of step with its solver; the script drops the argument
in its copy (`hst-main` itself is not touched).  Checked on the default
deck run to t = 2 from an energetic field (amplitude 0.003, kpeak 8: q2 =
0.12) with four snapshots on the RTX 3060:

| CPL result | against | agreement |
| --- | --- | --- |
| `rms.dat` (uu, vv, ww, uv per plane) of one field, plane average | the `variances_runtime.dat` line at that time | 1e-7 |
| the same over the four fields | the time average of our four lines | 0.2-1.4% (uu, vv, uv), exact for ww |
| `pField<n>.fld` recomputed by `prepare_pressure.cpl` (`PRESSURE=cpl`) | our online pressure, (0,0) mode excluded | 9e-5 relative |
| its (0,0) mode | ours | differs (ours is -<vv>_xz with zero mean; CPL's singular solve fixes the constant differently) |

The percent-level differences over several fields are the tool's
definition, not an error: `compute_re_stresses` subtracts the squared
*time-averaged* mean profile (the (0,0) mode, rms 0.02 here) from the
time-averaged second moments, while our lines subtract nothing (the mean
profile is part of the box energy) -- for w there is no mean profile and
the numbers are identical.  `uiuj.bin` and `mke.bin` (the budgets) are
produced from our pressure files; their contents were not checked against
anything.

**The deck (`examples/prod_re20000.in`).**  The size was left to the
user; this is the recommendation, chosen so that the run is resolved and
isotropic rather than as large as possible.  Resolution from the long
sheared run of WP5 (Re = 1000, Re_lambda 33, S* 6.3): eta = 0.0157
there (from Re_lambda and S* with q2 = 0.104, eps = 0.0165), i.e. the
default deck's resolution was dx/eta = dy/eta = 1.0 and dz/eta = 0.33
(the table of WP5 says dx/eta = 2 because it took nxd as the number of
x points; there are 2 nxd).  With the dissipation set by the box (eps =
S q2/S*, q2 about 0.1 whatever nu), eta = 0.0157 (1000/Re)^(3/4): Re =
20000 gives eta = 0.0017, and the grid 1536 x 1024 x 512 (nx = 511, ny =
1024, nz = 170, nzd = 512) has dx = dy = dz = 0.00195, dx/eta = 1.2, kmax
eta = 1.8 in x, 5 in z (the retained modes), Re_lambda about 33 sqrt(20)
= 150 (Sekimoto et al. reach 250).  That is one third of the points of
`bench_1024`, so about 13 GB per GPU on two A100 nodes; `nx = 3 nz` is
what WP5 asked for in this box, `nx = nz` (the benchmark decks) resolves
z three times finer than x for nothing.  CFL 1 with the fluctuations
alone in the CFL (u' about 0.2, peaks 3-4x): dt about 0.001, 80-100 k
steps to S t = 100; at a third of the 1024^3 step (0.33-0.45 s) 8-12 h on
two nodes.  Snapshots every 5 time units (8.6 GB velocity, 2.9 GB
pressure, 230 GB in all), restart every 5.  The initial field:
amplitude 0.003 at kpeak 8, i.e. q2 = 0.12 from the start (amplitude
0.05 gives q2 = 32 on any grid of this box, `bench_1024`'s log and the
64^3 run agree to all digits: the amplitude scales the potential, not the
rms), instead of the 1e-3 at kpeak 4 of the default deck (q2 = 4e-4,
S t = 30 of growth).

**The small production run (job 5171453, one A100 node, `dev_accelerated`,
`--time=00:15:00` so that `wall_max` = 300 s).**  The job script did its
part (the limit parsed, `wall_max` set, the run written to
`~/hst-runs/prod-small`), but the default deck as it stood went nowhere:
with `amplitude = 1e-3` at `kpeak = 4` the initial field has q2 = 3e-4
and decays (energy 2.9e-4 at t = 0, 1.2e-5 at t = 100) without ever
becoming turbulent, so the CFL step stays at 0.13 and S t = 100 takes
826 steps and 28 s.  The long sheared run of WP5 needed 57324 steps with
the "default deck": the initial-field module (the keyed random vector
potential) came later and its `amplitude` is not an rms, so the deck
that transitioned then does not now.  `hst.in` now starts from
`amplitude = 0.003` at `kpeak = 8` (q2 = 0.12; the same start went
turbulent in the local check above, dt = 0.005 at Re = 1000 on 64^3).
Measured on the way: 0.0075-0.0087 s/step for the 64 x 128 x 64 deck on
four A100 (`npy = 1`), and 1-2.3 s for each 51 MB velocity + 17 MB
pressure snapshot with its restart file, i.e. the collective MPI-IO of
these files costs a second of latency whatever the size.

**The small production run, second attempt (job 5171459, `prod-small-b`,
the deck with the energetic start).**  Turbulent from the start: 24523
steps to S t = 100 at 0.00719 s/step on four A100 (176 s; the step
settles at 0.0034-0.0039), which is less than the 300 s segment, so the
chain was not exercised by it (a third submission with `MARGIN=800`, i.e.
100 s segments, `prod-small-c`, does that).  Averages of `Runtimedata`
over S t = 30..100 (7001 samples) against the long sheared run of WP5:

| quantity | this run | WP5 (57324 steps, cflmax 0.8, RTX 3060) |
| --- | --- | --- |
| q2 | 0.132 | 0.104 |
| S* = S q2/eps | 6.09 | 6.3 |
| -uv/q2 | 0.166 | 0.159 |
| production / dissipation | 1.010 | 1.006 |
| Re_lambda | 37 | 33 |
| eta from eps = 0.0216 | 0.0147 (dx/eta = 1.06) | 0.0157 |

The CPL chain on HoreKA (`jobs/cpl_postprocess.sh ~/hst-runs/prod-small-b
2 4 10 1`, fields 4-10) needs `module load compiler/gnu/13` and then, as
a second command, `module load mpi/openmpi/5.0` (Lmod's hierarchy: in
one command the Intel-built OpenMPI stays loaded and gcc chokes on
`-xCORE-AVX2`), `~/.local/bin` on the PATH for `cpl`, and the `-L` paths
of `mpicc -show` in `LOADLIBES` (the script does that now: `mpicpl`
links with `-lmpi` alone).  `hst-main` is now on HoreKA too
(`~/hst-main`, an rsync copy).  Its `rms.dat` over the seven snapshots
gives q2 = 0.148 (the seven-sample average of a quantity whose
10000-sample average is 0.140) and the anisotropy b_uu = +0.13, b_vv =
-0.05, b_ww = -0.08 (our names; CPL's `vv` is our spanwise `ww`).

**The restart chain (jobs 5171488 and 5171583, `prod-small-c`,
`MARGIN=800`: 100 s segments).**  Segment 1 stopped at S t = 62.4 after
13507 steps ("wall-clock limit wall_max = 100 s reached"), wrote
`Dati.cart.out` and submitted segment 2, which waited 48 min in the
`dev_accelerated` queue, restarted from the file ("restarting from
Dati.cart.out") and ran to S t = 100 in 10199 steps: "the run reached
t_max: done", no third submission.  `Runtimedata` is continuous across
the two segments except that the restart repeats the line of its start
time (the initial `outstats` of every run; 10002 lines instead of
10001), harmless for averages and a visible marker of the restart.  The
chained run's averages over S t = 30..100 (S* = 6.17, -uv/q2 = 0.165,
P/eps = 1.018) agree with the unbroken run's (6.09, 0.166, 1.010) to the
sampling error: the two runs are different realizations, since the
statistics lines differ already at t = 0 in the last digits of `uw/2`
(the device reduction of the statistics is not bit-reproducible; the
fields of the regression are).  The CPL `rms.dat` of `prod-small-b`
over its seven snapshots gives -uv/q2 = 0.171 (its `uw` column is our
u v).

**Why the CPL pressure and ours differ by 1e-4 (the user asked).**  The
two codes discretize the same Poisson equation with the same compact
operators but different right-hand sides: CPL forms the nine velocity
gradients with the compact `D1`, multiplies them in physical space
(`2 (ux vy + ...)`), and applies `D0`; we form the six products and take
`D1` and `D2` of them (`D2 vv`, `2 i alfa D1 uv`, ...).  In the continuum
these are equal by continuity; discretely the product rule holds only
up to the truncation error, and the difference is that error.  Checked
on the default deck at t = 0.5 (one snapshot, 64 x ny x 64, the (0,0)
mode excluded), CPL's `prepare_pressure` against our online pressure:

| ny | max diff / max p | rms diff / rms p | diff/p in the kx = 0 plane by y-mode band m |
| --- | --- | --- | --- |
| 64 | 1.6e-3 | 4.4e-3 | 4e-4 (m < 4), 9e-4 (4-8), 1e-2 (8-16), 0.28 (16-32) |
| 128 | 9.2e-5 | 2.5e-4 | 2e-5, 5e-5, 7e-4 (8-16), 1e-2 (16-32), 0.45 (32-64) |
| 256 | 5.2e-6 | 1.4e-5 | 1e-6, 3e-6, 4e-5 (8-16), 7e-4 (16-32), 2e-2 (32-64), 1.2 (64-128) |

A factor 17.7 per halving of dy, i.e. fourth order, and the difference
sits in the highest y-modes: for a given band it also falls 16x per
halving.  The stencil says why fourth: the construction in
`setup_derivatives` makes `D4` exact on polynomials up to degree 8
(sixth order, relative error 2e-8 at k dy = 0.2) but `D1` and `D2` exact
only up to degree 4, fourth order (6.7e-6 at k dy = 0.2, 16x per
doubling of k), the classical five-point compact scheme of the CPL and
channel codes.  README.md called the scheme sixth-order; corrected.  So
the 1e-4 at ny = 128 is the y-truncation error of the pressure at that
deck's resolution (dy/eta = 1, k dy up to pi), not a defect of either
code, and the same number measures how far the pressure of the
production deck (dx/eta = 1.2) can be trusted in its dissipative range:
a few per cent at k dy = 1, as for any fourth-order DNS.

## The production run at Re = 20000 (2026-10-01 to 10-03, session 13, jobs 5171460 and 5174814)

`examples/prod_re20000.in` (1536 x 1024 x 512, Re = 20000, S = 1, box
3:2:1) to S t = 100 on two A100 nodes (8 ranks, 4 x-z pencils x 2 y
slabs, NCCL), output in `/hkfs/work/workspace/scratch/xt8786-hst/re20000`
(222 GB: 20 velocity snapshots of 8.6 GB, 20 pressure files of 2.9 GB,
the restart file).  Two segments of `jobs/horeka_prod.slurm`: job
5171460 waited 24 h in the `accelerated` queue, ran 11 h 51 min
(`wall_max` = 42600 s) to S t = 90.0 in 121119 steps, wrote its restart
file and submitted job 5174814, which waited 27 h and finished the last
10 time units in 10997 steps and 1 h 14 min.  132116 steps in all, 13.1
h of compute on 8 A100 (26 node-hours).

**Cost.**  0.344 s/step of compute, 0.352 with the I/O (the estimate was
0.33-0.45): snapshots (velocity + pressure, 11.5 GB) 19-29 s each,
restart files (8.6 GB) 9-18 s each in segment 1, i.e. 400-900 MB/s on
the workspace file system against 185 MB/s on the home file system for
the 1024^3 benchmark, and 1.6% of the segment; in segment 2 two restart
writes took 318 and 320 s (27 MB/s, the same code and file, other
nodes) and the last 9.7 s, so the workspace has its bad moments.
Memory: cuFFT 774 MB and two line-solver workspaces of 521 MB per rank;
the peak was not sampled (a third of the 1024^3 run's 39.3 GB, i.e.
about 13 GB).

**The time step could only shrink (bug, fixed in session 13).**
`compute_cfl` accumulated `cfl` with `reduction(max:cfl)` without ever
resetting it, so the CFL-chosen step was the smallest the run had ever
needed: 0.00132 from the initial field, 0.00074 from S t = 2 (the
transient's velocity peak) to the end of segment 1, and 0.00091 after
the restart, whose fresh `cfl` shows what the turbulent field itself
allows.  The small runs did the same (0.00507 -> 0.00341 in steps,
never up).  The step was therefore always safe, only short: the
production run spent about 20% more steps than CFL = 1 needs, and the
printed CFL column was 1.0 by construction (`cfl * deltat`).  Now `cfl =
0` at the top of `compute_cfl`; the fixed-step decks of the tests are
unaffected and the safety net is green.

**Statistics (`Runtimedata`, S t = 30..100, 7001 samples; `energy` =
<q2>, `diss` = <grad u : grad u> for ly = 2).**

| quantity | Re = 20000, 1536 x 1024 x 512 | Re = 1000, 64 x 128 x 64 (`prod-small-b`) | literature |
| --- | --- | --- | --- |
| q2 | 0.083 | 0.132 | |
| eps = diss/Re | 0.0112 | 0.0216 | |
| S* = S q2/eps | 7.4 | 6.1 | 5..7 (Rogers & Moin 1987), up to 8 at high Re_z (Sekimoto, Dong & Jimenez 2016) |
| -uv/q2 | 0.136 | 0.166 | 0.15 (Tavoularis & Karnik 1989) |
| production / dissipation | 1.005 | 1.010 | 1 |
| b_uu, b_vv, b_ww (our names: v vertical) | +0.14, -0.07, -0.07 | +0.13, -0.05, -0.08 | +0.2, -0.14, -0.06 |
| Re_lambda = (q2/3) sqrt(15/(nu eps)) | 143 | 37 | |
| eta = (nu^3/eps)^(1/4) | 0.00183 | 0.0147 | |
| dx/eta = dy/eta = dz/eta | 1.07 | 1.06 (dz/eta 0.35) | |

The resolution came out where the deck meant it to (eta = 0.0017 was
the estimate, 0.00183 the result): kmax eta = 2.0 in x, 6 in z, the
y grid at dy/eta = 1.07 with the compact scheme's fourth-order
truncation (the pressure section above).  The flow is stationary from
S t = 30 (the time series of q2 wanders between 0.05 and 0.11 over 10-20
time units, the usual bursting of a shear-periodic box), production
balances dissipation to 0.5%, and S* and -uv/q2 move the way the
literature says they do with Re_z (S* up, -uv/q2 down from the small
box).  Re_lambda = 143 is in Sekimoto et al.'s range (their largest
boxes reach 250).

**The chain as a workflow.**  Two things to keep: the segments cost
nothing but the queue (24 h and 27 h here, for 12 h and 1.25 h of
work), so one long segment beats several short ones on `accelerated`;
and the restart repeats one `Runtimedata` line (10002 lines for S t =
0..100 at dt_stat = 0.01).  The CPL post-processing of the 15 snapshots
from S t = 30 runs as `jobs/horeka_postprocess.slurm` on a `cpuonly`
node (job 5181558, 4 ranks, about 50 GB each).

## Post-processing in Fortran (2026-10-05, session 13; the user asked)

`src/postpro/postpro.f90`, `make postpro`, driven by the namelist
`postpro.in` (snapshots, deck, output directory, and switches: `mean`,
`stresses`, `spectra`, `budgets` = the list of components).  It is the
solver's modules with a loop over snapshots on top: `restart_read` reads
a snapshot (its time with it), `fill_ghosts` the shear-periodic rows,
`compute_pressure` the pressure (so `p_fields/` is not needed),
`line_solve(KIND_DY)` the compact d/dy of the three components, and the
nine gradient fields and the pressure go to physical space three at a
time through `V` and `transform_to_physical` (the physical field is
read off `rVVdx`, which the CFL already uses that way).  The plane
statistics are sums over this rank's modes or physical points and rows,
reduced over the ranks at the end, so any decomposition works.  The
profile derivatives (the production's mean slope, the transports, the
viscous diffusion) use the solver's own stencils: `D0 f' = D1 f` is a
periodic pentadiagonal system of size ny, stored dense and
LU-factorized once (second-order central differences, the first
attempt, left the derivative terms 1-5% off the CPL ones).  Output: text
files with a header, `mean.dat`, `stresses.dat`, `budget_<ij>.dat`,
`spectra_xz.dat` (y-averaged, kz folded), `spectra_x.dat`, `spectra_z.dat`.

Checks on the default deck's snapshot at t = 2 (64 x 128 x 64) against
the CPL chain on the same file: the stresses, the production, the
dissipation, the pressure-strain, both transports and the viscous
diffusion of all six components agree to 5e-11 (CPL writes its binaries
with that precision); the sums of `spectra_xz.dat` equal the y-averaged
stresses; two x-z pencils x two y slabs on four CPU ranks and the GPU
build give the same files to the printed eight digits; four snapshots
with `budgets = 'uu uv'` average as they should.  The CPL names differ
from ours (their v is our w): `jobs/cpl_postprocess.sh` keeps running the
CPL chain as the cross-check, now with CPL's `postpro.in` inside its
build directory so that it does not overwrite ours.  Not ported: the MKE
budget (zero in a box without mean pressure gradient), the pressure
decomposition, the pressure-strain spectra, the VTK output.
