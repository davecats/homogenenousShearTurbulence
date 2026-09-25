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
