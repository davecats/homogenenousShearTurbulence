# CPL post-processing of velocity snapshots

`post/cpl/` holds a post-processor written in CPL.  It reads the velocity snapshots
`fields/field<n>.fld` (and pressure snapshots `p_fields/pField<n>.fld`) saved by
`hst` during a run and computes statistics from them: means, Reynolds stresses,
Reynolds-stress budgets, velocity spectra, and optionally the pressure itself and
VTK files for ParaView.  Everything is computed per *phase* (see below), which is
how unsteady runs (oscillating spanwise shear `S2(t)`) are phase-averaged; with one
phase it is the plain average over all snapshots.

The statistics that `hst` itself writes while running (`Runtimedata`,
`variances_runtime.dat`) are not produced here; they are read by `post/python/`
and `post/cpl/meanandvar_utils/`.

## Reference system

The post-processor uses the same convention as the Fortran code `hst`:

| axis | direction | velocity | spectral / discrete |
|---|---|---|---|
| x | streamwise | `u` | Fourier modes `ix = 0..nx`, `kx = alpha0*ix` |
| y | vertical (shear direction) | `v` | `ny` mesh rows `iy = 0..ny-1`, compact finite differences |
| z | spanwise | `w` | Fourier modes `iz = -nz..nz`, `kz = beta0*iz` |

The mean flow is `U = S*y`, `W = S2(t)*y`.  `nx`, `ny`, `nz`, `ly` are the values of
the `&mesh` group of `hst.in`.  Rows `-2,-1,ny,ny+1` are shear-periodic ghost rows, as
in `hst`; the grid `y(iy)` (uniform, or clustered with `ystretch`) is that of `hst`.

All names follow this: `uv` is the streamwise-vertical correlation (the shear stress),
`uw` the streamwise-spanwise one, `vw` the vertical-spanwise one, `Euv` the
corresponding spectrum, and so on.

(The legacy CPL code called the spanwise axis `y` and the vertical one `z`.  The
velocity files still carry those legacy labels, because that is the file layout
`hst` writes; the translation is done in `src/hstfieldio.cpl` and nowhere else.)

## Build and run

Requirements: `cpl`/`mpicpl`, an MPI installation (`mpirun`), `python3`.

```bash
post/cpl/run_postprocess.sh [CASE_DIR]
```

It can be run from any folder.  The case-specific settings are at the top of the
script: `DATA_DIR` (or the first argument), `NF_MIN`, `NF_MAX` (snapshots
`field<NF_MIN>.fld .. field<NF_MAX>.fld`), `DN` (number of phases), `NPROCS`
(MPI ranks) and `BUILD` (set to 0 to skip recompiling).  `CASE_DIR` must contain
`hst.in` and the directories `fields/`, `p_fields/`, `statistics/` (and `vtks/` when
the VTK output is on).

The script converts `hst.in` to the flat token file the CPL program can read
(`src/convert_hstin.py`), writes `postpro.in`, and calls
`mpirun -np NPROCS postprocess postpro.in readinput.in`.  Both files are temporary.
The snapshots are distributed over the MPI ranks (each rank processes whole
fields), so the memory per rank is that of one field with its derivatives, and more
ranks than snapshots are useless.

## What is computed (compile-time flags)

The flags are in `src/flags.cpl` (`#define` = on, `! #define` = off):

| flag | what | output |
|---|---|---|
| `compute_pressure` | recomputes the pressure from the velocity and overwrites `p_fields/` | `p_fields/pField<n>.fld` |
| `compute_variances` | mean velocities and Reynolds stresses | `statistics/velstats.dat` |
| `compute_budget` | budgets of the Reynolds stresses | `statistics/budget_<ij>.dat` |
| `compute_spectra` | velocity spectra | `statistics/spectra*.{bin,dat}` |
| `convert_to_vtk` | velocity to VTK | `vtks/field<n>.vtk` |

Without `compute_pressure` the pressure files written by `hst` are used (they are
read even if no pressure-dependent statistic is requested, so `p_fields/` must
exist).  The pressure is only defined up to a constant; its mean `(0,0)` mode is
arbitrary and removed by the Reynolds decomposition.

### Phases and the Reynolds decomposition

Snapshot `nf` belongs to phase `phase_num = ((nf-1) mod dn) + 1`, so `dn` is the
number of snapshots per period of the forcing (`dn = 1`: no phase resolution).
The program makes three passes over the snapshots:

1. (optional) pressure;
2. means: the mean `<u>, <v>, <w>, <p>` of each phase is the average over the
   fields of the phase of the **y-average of the `(kx,kz) = (0,0)` mode**
   (statistically homogeneous flow, so the mean profile is averaged over y), plus the
   phase mean of `S2`;
3. statistics: each field is read again, the phase mean is subtracted from its
   `(0,0)` mode (`apply_rey_decomp`), and the fluctuation is accumulated in the
   statistics of its phase.

At the end each phase is divided by its number of snapshots.  Every phase needs at
least one snapshot.  All y-averages use the same quadrature (`yintegr`, a
three-point rule on pairs of rows, so `ny` must be even).

### Statistics files

* `velstats.dat`: one line per phase, `# phase U V W uu vv ww uv uw vw` (means of
  `u`, `v`, `w` and the six Reynolds stresses, all y-averaged).
* `budget_uu.dat`, `budget_vv.dat`, `budget_ww.dat`, `budget_uv.dat`,
  `budget_uw.dat`, `budget_vw.dat`: one line per phase, columns
  `# phase uiuj prod psdiss pstrain forcingprod`: the stress itself, and the terms of
  its budget: production by the mean shears `S` and (phase-mean) `S2`, pseudo-spectral
  (dealiased) dissipation, pressure-strain, and the production by the fluctuation of
  `S2` about its phase mean (`forcingprod`).
* `pField<n>.fld`: same array layout as the velocity but one component and no header.
* `vtks/field<n>.vtk`: binary legacy VTK, rectilinear grid, `2*nxd x ny x nzd`
  points (the 3/2-padded streamwise and spanwise grids), axes x streamwise, y
  vertical, z spanwise, velocity `(u,v,w)`.

## Velocity spectra

### Definition

For every mesh row `iy`, streamwise mode `ix` and spanwise mode `iz`, the fluctuating
Fourier coefficients `û_i(ix,iz,iy)` give

```
E_ij(iy, ix, |iz|) = sum over +iz and -iz of  c * Re[ û_i * conj(û_j) ],   c = 1 (ix = 0), 2 (ix > 0)
```

* **Components** `ij`, always in this order: `1:uu 2:vv 3:ww 4:uv 5:uw 6:vw`
  (`uv` is the streamwise-vertical co-spectrum, `uw` the streamwise-spanwise one).
  The cross-terms are co-spectra (real part only).
* **Factor `c`**: only the modes `ix >= 0` are stored (real fields); `c = 2` accounts
  for `-ix`.  This is the same weighting as in the variances.
* **Folding in `kz`**: the modes `+iz` and `-iz` are added, so `|iz| = 0..nz`.
* These are *mode energies*, not densities: the sum over all `ix` and `|iz|` of
  `E_ii` at a row is the variance at that row.  For a spectral density divide by
  `alpha0*beta0` (= `dkx*dkz`).
* The fluctuation is taken with the phase-wise Reynolds decomposition described
  above, so with `dn > 1` there is one spectrum per phase, averaged over the
  snapshots of that phase.  Only the `(0,0)` mode depends on the decomposition.
* The y-average of the spectra, `sum_iy w(iy) E(iy,...)` with the quadrature weights
  `w` of the code, summed over all `ix` and `|iz|`, gives exactly the corresponding
  entry of `velstats.dat`.

### Files

**`statistics/spectra.bin`: the 2D spectra of all phases, y-resolved.**

* Binary, raw 8-byte IEEE reals (`float64`), native (little-endian) byte order,
  **no header**.
* C order (last index fastest), i.e. an array of shape

  ```
  (dn, 6, ny, nx+1, nz+1)   indexed   (phase, component, iy, ix, |iz|)
  ```

  with `phase` = 0..dn-1 (phase 1 first), `component` as above, `iy` the mesh row
  (`y(iy)` as in `hst`), `ix` the streamwise mode (`kx = alpha0*ix`) and `|iz|` the
  absolute spanwise mode (`kz = beta0*|iz|`).
* Size in bytes: `8 * dn * 6 * ny * (nx+1) * (nz+1)`.  The file stores no sizes, so
  `nx, ny, nz` (from `hst.in`) and `dn` (the value used in `run_postprocess.sh`) are
  needed to read it.

**`statistics/spectra_kx.dat`, `statistics/spectra_kz.dat`**: text, y-averaged 1D
spectra, one block per phase, headers
`# phase ix kx Euu Evv Eww Euv Euw Evw` and `# phase iz kz Euu Evv Eww Euv Euw Evw`.
`E(kx)` is the 2D spectrum summed over `kz`, `E(kz)` (`kz >= 0`, `+-kz` folded) the
2D spectrum summed over `kx`.

### Reading the data in Python

Run from the repository root (or put it on `PYTHONPATH`), using the helpers of
`post/python`:

```python
import numpy as np
from post.python.config import load_config
from post.python.grid import make_grid

def y_weights(y, ly):
    """w(iy) such that sum_iy w(iy) f(iy) is the y-average of f (the quadrature of the CPL code)."""
    ny = len(y)
    ye = lambda i: y[i % ny] + ly * (i // ny)        # periodic extension of the mesh
    w = np.zeros(ny)
    for iy in range(0, ny, 2):                       # three-point rule on (iy-1, iy, iy+1)
        yp, ym = ye(iy + 1) - ye(iy), ye(iy - 1) - ye(iy)
        a1 = -ym / 3 + yp / 6 + yp * yp / (6 * ym)
        a3 = yp / 3 - ym / 6 - ym * ym / (6 * yp)
        w[(iy - 1) % ny] += a1
        w[iy] += yp - ym - a1 - a3
        w[(iy + 1) % ny] += a3
    return w / ly

case, dn = "/path/to/case", 2                        # dn as used in run_postprocess.sh
cfg = load_config(f"{case}/hst.in")
E = np.fromfile(f"{case}/statistics/spectra.bin", np.float64)
E = E.reshape(dn, 6, cfg.ny, cfg.nx + 1, cfg.nz + 1)         # (phase, component, iy, ix, |iz|)
kx = cfg.alfa0 * np.arange(cfg.nx + 1)
kz = cfg.beta0 * np.arange(cfg.nz + 1)
y, _ = make_grid(cfg)

Evv_plane = E[0, 1, 10]                                      # E_vv(kx, |kz|) of phase 1 at row iy = 10
E_avg = np.einsum("y,pcyxz->pcxz", y_weights(y, cfg.ly), E)  # y-averaged E(kx, |kz|)
E_kx = E_avg.sum(axis=3)                                     # y-averaged E(kx):  (phase, component, ix)
E_kz = E_avg.sum(axis=2)                                     # y-averaged E(|kz|): (phase, component, |iz|)
```

`E_kx` and `E_kz` reproduce `spectra_kx.dat` and `spectra_kz.dat` to round-off.
For a uniform mesh `w(iy)` is `2/(3 ny)` for odd and `4/(3 ny)` for even `iy`.
For a spectral density, divide by `cfg.alfa0 * cfg.beta0`.

### Memory

While running, every rank holds the spectra of all `dn` phases:
`8 * dn * 6 * ny * (nx+1) * (nz+1)` bytes (about 75 MB per phase for
`nx=191, ny=128, nz=63`), plus one phase as receive buffer.  The phases are reduced
over the ranks and appended to `spectra.bin` one at a time.

## Source files

| file | content |
|---|---|
| `postprocess.cpl` | main program: the passes over the snapshots, flag handling |
| `src/flags.cpl` | compile-time flags |
| `src/readinput.cpl` | reads `readinput.in` (mesh, physics, time control); `postpro.in` is read in `convenience_postprocess.cpl` |
| `src/convert_hstin.py` | converts `hst.in` (Fortran namelist) to the flat token file |
| `src/ydiscret.cpl` | y mesh, compact finite differences, shear-periodic shifts, y quadrature |
| `src/convenience_postprocess.cpl` | postprocessing input, MPI work distribution, field arrays, velocity gradient |
| `src/hstfieldio.cpl` | reading the velocity files (translation of the legacy layout) |
| `src/io_postprocess.cpl` | banner and `read_vel_field`/`read_press_field` |
| `src/linsolver_smw.cpl` | solver for the cyclic pentadiagonal systems |
| `src/poisson.cpl` | pressure from the velocity |
| `src/mod_rey_decomp.cpl` | phase means and the Reynolds decomposition |
| `src/mod_compute_variances.cpl` | `velstats.dat` |
| `src/mod_compute_budget.cpl` | budgets |
| `src/mod_compute_spectra.cpl` | spectra |
| `src/mod_fld2vtk.cpl` | VTK output |
| `meanandvar_utils/` | mean and uncertainty of a time series (`Runtimedata`), independent of the rest |

## Velocity file layout (as read by the post-processor)

Text header up to a line `Vfield=`, then a `complex128` array stored with the last
index fastest.  The file keeps the legacy CPL labels: the header reads `nx=`, `ny=`
(= our `nz`, spanwise modes), `nz=` (= our `ny + 1`), and the array is
`(0..nx, -nz..nz, rows)` of `(u, w, v)` (the spanwise component first) with
`rows = -2..ny+1` (`ny+4` rows including the ghost rows).  The scalars `time`, `S`,
`S2`, `gamma_x`, `gamma_y` (our spanwise image displacement) follow the header as raw
8-byte reals.  The file layout is written by `src/hst_io.f90` at the top level of
the repository, which is the reference.
