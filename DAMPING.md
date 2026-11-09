# Large-scale damping in an enlarged box (branch `damping_from_main`)

Session 16 (2026-11-06).  Third route to a turbulence whose statistics do
not depend on the angle between the mean shear and the box, after the
hexagonal lattice (`feature_honeycomb_from_main`, HONEYCOMB.md: reduces
but does not remove the dependence) and the co-rotating frame
(`feature_rotating_shear`, ROTATING.md: a different flow, relaminarises
for fast swings).  Branched from `main`; nothing of the other two
branches is here.

## 1. Idea

In a periodic box the energy of homogeneous shear turbulence sits in the
lowest modes, (kx, kz) = (1, 0), (0, 1), (1, 1) of the 3:2:1 box: the box
periodicity is the only length scale, and three wavevectors sample the
(|k|, angle) plane very coarsely.  Turning the shear by theta changes
which |k| and which angles relative to the shear those three vectors
have, hence the dependence on theta.

The box is enlarged n times in x and z at the same resolution, so the
lattice of wavevectors is n times finer in every direction, and a body
force removes the modes that the small box did not have: those inside
the ellipse with the semi-axes 2 pi/Lx, 2 pi/Lz of the small box, defined
in the shear frame.  The modes on and outside the ellipse are the small
box's own modes and their n^2 - 1 new neighbours.  The turbulence keeps
the length scale set by the ellipse, and the ellipse, being a constraint
in spectral space and not a periodicity, can have any orientation and can
turn with the shear (which no lattice can, ROTATING.md section 1).  The
residual angle dependence should fall like 1/n.

A set of modes that is not a lattice (non-equispaced Fourier modes) is
not an alternative: the nonlinear term needs the set closed under
k1 + k2, so any such set is a Galerkin truncation on a finer lattice,
which is exactly this with a hard mask, at the cost of the fine box's
transforms in any case.

The box in y is unchanged (Ly = 2, the shear-periodic tilt as before): the
mask acts on the x-z columns (kx, kz) only.

## 2. Implementation

Deck group `&damping` (optional; absent = off):

    damp_kpar, damp_kperp   semi-axes along and across the shear (both > 0: on)
    damp_norm               exponent p of the norm r = ((|k_par|/kpar)^p + (|k_perp|/kperp)^p)^(1/p); 2 = ellipse
    damp_hard               .true.: the columns with r < 1 are zeroed every substep (hard mask)
    damp_rate               soft damping rate g inside (explicit; g deltat < 1)
    damp_edge               width in r of a smoothstep ramp of g to zero at r = 1 (0 = step)
    damp_mean               damp the (0,0) column (mean profiles) too (default no)
    damp_follow             the ellipse follows the shear direction atan2(S2(t), S) (default); else damp_angle (deg)

Columns on the boundary (r = 1 within 1e-6) are kept, so with the
semi-axes of a small box its own modes are the first ones kept.

* `hst_damping` (new): `gdamp(iz, ix)` (rate, or 1 for a masked column),
  rebuilt on the host and copied to the device whenever the direction
  changes (once for a constant S2; every substep for the swing);
  `damping_mask` drops the right-hand sides of the masked columns before
  the implicit solve, so they come out zero.
* `hst_equations`: soft damping as the explicit terms `-g d2v`, `-g eta`
  (mean column `-g (u, w)`) next to the tilting term; the hard mask after
  the nonlinear terms of every substep.
* `hst_stats`: `damping_runtime.dat` (time, direction in degrees, number
  of damped columns, energy in them, power 2 g <u.u> removed; zero for the
  hard mask, which removes no power: it is a Galerkin truncation and the
  resolved modes keep the ordinary budget production = dissipation +
  dE/dt).  `spectra/exz<i>.dat` with every field snapshot: the column
  energies E(kx, kz) averaged over the dt_stat samples since the previous
  snapshot (all runs, damped or not), to see where the energy sits.

## 3. Verification

`tests/test_kelvin` with three decks: `kelvin_damped.in` (the mode (1,1)
inside a soft ellipse, rate 2: the closed form times exp(-2 t), relative
error 1e-8), `kelvin_hard.in` (inside the hard mask: stays exactly zero),
`kelvin_outside.in` (outside: unchanged, 4e-9).  The whole suite (15
decks) passes on the CPU (1 and 2 ranks) and GPU builds.

## 4. Runs

`~/Codes/hst/damping_campaign`: `tools/make_runs.py` writes the decks of
the 3:2:1 box enlarged n = 2, 3 times (nx 96/144, nz 48/72, ny 95,
Re 1000, |S| = 1, direction theta, ellipse (2 pi/3, 2 pi) following the
shear, random start with amplitude 0.003/n) to be compared with the
honeycomb rect runs (n = 1).  Phase 1: hard mask, theta = 0 and 45.

Results so far (2026-11-07; averages after t = 400 in the frame of the
shear, errors from 50-unit blocks; t = length of the run; all 20 runs
complete, 2026-11-08):

    case     theta   t    q2      eps     S*     -<u_par v>/q2  b_par   b_v    b_perp  rho_par,perp  rho_v,perp
    n1 (small)  0  2040  0.108   0.0174  6.23    0.163          0.112  -0.049  -0.064   0.003        0.000
    n1 (small) 45  2040  0.096   0.0132  7.32    0.140          0.014  -0.123   0.110   0.322       -0.102
    n2 hard     0  1400  0.126   0.0202  6.25    0.161          0.075  -0.092   0.017  -0.001        0.000
    n2 hard    45  1400  0.082   0.0126  6.51    0.155          0.050  -0.125   0.075   0.042       -0.002
    n3 hard     0  1400  0.119   0.0189  6.31    0.159          0.066  -0.102   0.037  -0.004        0.000
    n3 hard    15  1400  0.124   0.0201  6.16    0.163          0.081  -0.105   0.023  -0.001        0.001
    n3 hard    30  1400  0.098   0.0156  6.26    0.160          0.068  -0.113   0.045   0.013        0.002
    n3 hard    45  1400  0.117   0.0188  6.22    0.162          0.088  -0.109   0.020   0.010       -0.002
    n4 hard     0  1400  0.121   0.0192  6.31    0.159          0.066  -0.104   0.038  -0.001        0.001
    n4 hard    45  1400  0.126   0.0204  6.16    0.163          0.081  -0.102   0.021  -0.012        0.003
    n2 hard    15  1400  0.122   0.0208  5.86    0.172          0.099  -0.100   0.001   0.006       -0.002
    n2 hard    30  1400  0.103   0.0169  6.08    0.166          0.072  -0.108   0.036  -0.025        0.004
    n4 hard    15  1400  0.121   0.0195  6.20    0.162          0.078  -0.106   0.029   0.029       -0.008
    n4 hard    30  1400  0.116   0.0183  6.35    0.158          0.069  -0.109   0.040   0.004        0.001
    n5 hard     0  1400  0.127   0.0200  6.33    0.158          0.070  -0.104   0.034   0.000        0.000
    n5 hard    15  1400  0.116   0.0187  6.22    0.161          0.074  -0.108   0.034  -0.003        0.000
    n5 hard    30  1400  0.107   0.0172  6.23    0.161          0.076  -0.112   0.036  -0.011        0.001
    n5 hard    45  1400  0.114   0.0180  6.32    0.158          0.070  -0.111   0.041  -0.005        0.002
    n10 hard    0  1400  0.116   0.0185  6.29    0.159          0.071  -0.109   0.037  -0.000       -0.000
    n3 rate 10  0  1400  0.108   0.0166  6.47    0.164          0.074  -0.106   0.032  -0.002        0.000
    n3 rate 10 15  1400  0.115   0.0182  6.31    0.168          0.091  -0.108   0.016  -0.005        0.002
    n3 rate 10 30  1400  0.090   0.0140  6.41    0.166          0.079  -0.117   0.038   0.012        0.000
    n3 rate 10 45  1400  0.106   0.0164  6.42    0.166          0.097  -0.114   0.016   0.003        0.002
    n3 rate 3   0  1400  0.096   0.0145  6.61    0.170          0.090  -0.109   0.020  -0.002        0.001
    n3 rate 3  15  1400  0.102   0.0157  6.50    0.174          0.102  -0.112   0.009  -0.005        0.002
    n3 rate 3  30  1400  0.075   0.0115  6.55    0.172          0.097  -0.122   0.025   0.009        0.004
    n3 rate 3  45  1400  0.089   0.0134  6.60    0.172          0.113  -0.120   0.007   0.009       -0.001
    n3 rate 1   0  1400  0.080   0.0118  6.76    0.174          0.107  -0.115   0.008  -0.003        0.001
    n3 rate 1  15  1400  0.083   0.0126  6.64    0.177          0.119  -0.117  -0.001  -0.009        0.004
    n3 rate 1  30  1400  0.061   0.0092  6.71    0.176          0.121  -0.129   0.008   0.010        0.005
    n3 rate 1  45  1400  0.070   0.0104  6.74    0.175          0.133  -0.128  -0.005   0.016       -0.005

Spread (max - min)/mean over the angles 0-45:

    case   angles   q2     eps    S*     -<u_par v>/q2   b_par  b_v   b_perp
    n1       7     0.18   0.29   0.16   0.15            3.0    1.0   5.9
    n2       4     0.41   0.46   0.11   0.10            0.66   0.32  2.3
    n3       4     0.23   0.25   0.02   0.02            0.30   0.10  0.80
    n4       4     0.08   0.11   0.03   0.03            0.20   0.06  0.58
    n5       4     0.17   0.15   0.02   0.02            0.08   0.07  0.20

(b_perp is 0.02-0.04 for n >= 3, so its relative spread means little; the
absolute spread is 0.02 at n = 3, 0.007 at n = 5.)

* The turbulence is sustained behind the hard mask (the masked columns
  are exactly zero; damping_runtime.dat).  The energy sits in the first
  ring outside the ellipse (1 <= r < 2: 56-58 %; 2 <= r < 3: 18-20 %; the
  (0,0) column 1-2 %) for every n and angle.
* The structure of the turbulence is angle independent from n = 3 on:
  S* = 6.2-6.3, -<u_par v>/q2 = 0.16, the anisotropy b_ij and the
  correlation coefficients are those of the small box at theta = 0 at
  every angle, within 1-3 %, and the cross correlations
  rho_par,perp, rho_v,perp of the turned small box (0.32, -0.10 at 45)
  are gone (|rho| < 0.015).  (With n = 2 the lattice, spacing half the
  ellipse axes, still decides which few modes sit just outside the
  ellipse: theta = 45 is a different flow, q2 0.082, S* 6.5.)
* The intensity level (q2, eps together, at constant S*) still depends
  on the angle and does not converge from n = 3 to n = 5 (spread of q2
  23, 8, 17 %; history.png: each run sits at its level from t = 0 on,
  the block errors are 1 %).  The level is a lattice-counting effect:
  the energy per column is largest just outside the ellipse, and the
  number of lattice columns in the innermost band r < 1.1 depends on
  the angle (counted from spectra/exz*.dat, columns of both signs):

        n   theta 0   15   30   45      q2 0 / 15 / 30 / 45
        2         4    4    0    0      0.126 0.122 0.103 0.082
        3        12   10    2    6      0.119 0.124 0.098 0.117
        4        16   14   12   16      0.121 0.121 0.116 0.126
        5        28   14   10   20      0.127 0.116 0.107 0.114
        10       70   66   68   70      0.116   -     -   (pending)

  q2 follows the count at every n.  At theta = 0 the lattice has modes
  exactly on the ellipse (the axis modes, and for n = 5 also (3,4) and
  (4,3): 12 columns with 8 % of the energy), at a generic angle none;
  the count in the band fluctuates by a factor 2-3 between angles, as
  the number of lattice points in a thin ring does, and its relative
  scatter falls like 1/n: at n = 10 the count is 66-70 at the four
  angles (6 % scatter), so the level spread should fall to a few per
  cent.  The structure, which is set by the whole first ring
  (1 <= r < 2 holds 57 % of the energy at every n and angle), is
  converged from n = 3 on.
* Soft damping (n = 3, rates g = 1, 3, 10 in units of S, step edge, four
  angles each; tools/soft.py, soft.png) is a worse mask than the hard
  one and does not reach it quickly.  Averages over t > 400, the same at
  all four angles within a few per cent:

      rate   E_in/q2    removed power g<u.u>_in/eps   q2 (theta 0)   S*
      1      2.6 %      0.17                          0.080          6.7
      3      0.6 %      0.12                          0.096          6.6
      10     0.09 %     0.057                         0.108          6.4
      hard   0          0                             0.119          6.3

  The hard mask is a Galerkin truncation: a triad with a zero mode
  contributes nothing to the other two, so the resolved modes keep the
  ordinary budget (production = dissipation) and the mask removes no
  power.  A finite rate leaves a small amplitude in the columns inside
  the ellipse, and the nonlinear transfer feeds them continuously, so the
  damping drains a sizeable fraction of the dissipation: the result is a
  weaker flow (q2 down 10-35 %) with a higher S* and stronger
  correlation (-<u_par v>/q2 0.164-0.176 against 0.16), which approaches
  the hard one only as the rate grows (the removed power falls roughly
  like 1/g; g deltat < 1 bounds the rate at about 100 with the deltat =
  0.0064 of these runs, so a rate near the hard limit is not reachable).
  The soft flows are as angle independent as the hard one (S* within
  2 %, |rho_par,perp| < 0.02) and have the same 25 % dip of the level at
  theta = 30, which is therefore a property of the n = 3 lattice and not
  of the mask.  Conclusion: the hard mask is the right choice; the soft
  rate stays available for the swinging shear should the hard mask prove
  too abrupt when the ellipse turns (phase 4).
* Phase 4 (2026-11-08/09): the n = 10 check, theta = 0 and 45 (nx 479,
  nz 240: four x-z pencils need nx + 1 divisible by four; one A100 node
  with four ranks per run, tools/horeka_node.slurm in self-resubmitting
  segments, jobs 5231128/5231129), to see whether the level spread falls
  as 1/n.  theta = 0 is done (242720 steps, 0.112 s/step, 7.5 h on one
  node; the processes hung on Lustre after the final restart write, so
  the job ran to its limit): q2 0.116, S* 6.29, -<u_par v>/q2 0.159,
  b_ij and rho_ij as at n = 3-5, the energy 58 % in 1 <= r < 2 and
  0.09 % inside the ellipse (161 masked columns).  theta = 45 pending.
* Tools: tools/analyse.py (angle_dependence.png with the correlation
  coefficients as the bottom row, results.csv), tools/spectrum.py
  (spectra.png in the shear frame), tools/history.py (history.png),
  tools/soft.py (soft.csv, soft.png: what the damping removes vs rate);
  analyse.py --only hard --out angle_dependence_hard.png for the
  hard-mask cases alone.
