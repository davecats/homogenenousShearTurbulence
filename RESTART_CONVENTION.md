# Restart convention: old-CPL → new hst — verification

**Date:** 2026-09-28
**Files examined:** `run/Dati.cart.in.Re2000`, `src/hst_io.f90` (Fortran restart
reader), `post/python/io.py` (Python reader), `run/hst.in` (deck).

## 1. The two conventions

| | streamwise | spanwise | shear-wise |
|---|---|---|---|
| **Old CPL** (`Dati.cart.in.Re2000`) | x → u | **y → v** | **z → w** |
| **New hst** (this repo)             | x → u | **z → w** | **y → v** |

The change is a relabelling of two transverse directions and their velocity
components: the old **y/spanwise ↔ new z/spanwise** and old **z/shear ↔ new
y/shear**; and correspondingly the old spanwise velocity **v ↔ new w** and the
old shear velocity **w ↔ new v**.

## 2. What the on-disk format stores (convention-independent)

A CPL velocity file stores the array

```
ARRAY(0..nx, -ny_cpl..ny_cpl, -1..nz_cpl+1)   of (u, v_cpl, w_cpl)
```

in C order (last index fastest, three components innermost) with a text header
`nx=.. ny=.. nz=..`.  The header integers are the **CPL names**:

```
dim0 = streamwise,        size nx+1
dim1 = spanwise modes,    size 2*ny_cpl + 1      (ny_cpl = spanwise modes)
dim2 = shear rows,        size nz_cpl + 3        (nz_cpl = shear points + 1)
comp = (u, v_cpl, w_cpl)  where v_cpl = spanwise velocity, w_cpl = shear velocity
```

Crucially, the layout stores *physical axis types* (x, spanwise, shear), not the
old labels `y`/`z`.  So when a code switches from calling the spanwise direction
`y` to calling it `z`, **the array geometry does not change at all** — only the
names and the meaning of the two transverse velocity slots change.  That is why
an old-CPL file needs no transpose when re-read: file dim1 (spanwise) already
maps to the new `z` axis and file dim2 (shear) to the new `y` axis.

## 3. The mapping the readers apply

**Fortran** — `src/hst_io.f90::restart_read`:
 - Mesh compatibility: the reader accepts differing resolutions.  The deck
  tuple `(nx, nz, ny+1)` corresponds to file `(nx, ny_cpl, nz_cpl)`; common
  signed Fourier indices are copied and absent target modes remain zero.
  The physical box must be unchanged (in particular `ly`; `ly` is not stored
  in the legacy CPL header), and the Fourier fundamentals `alpha0` and
  `beta0` must match.
- `CPL_ORDER = [1, 3, 2]`; the copy loop
  `V(iy, iz, ix, CPL_ORDER) = buf(:, iy, iz, ix)` sends CPL component 1 (old
  spanwise v) to the new `w` slot and CPL component 2 (old shear w) to the new
  `v` slot.  `cpl_view` maps file row `r` → slot `r+2` and file dim1 (spanwise)
  onto the new `iz` index, file dim2 (shear rows) onto the new `iy` index.

**Python** — `post/python/io.py::read_field`:
- `a = raw.reshape((nx+1, 2*ny_cpl+1, nz_cpl+3, 3))`, transposes to
  `(row, span, ix, comp)`, keeps interior rows `2:2+ny_hst`.
- `vel = (u, v, w)` with `v = CPL comp2 (w_cpl, shear)`, `w = CPL comp1
  (v_cpl, spanwise)` — the same swap.

The Fortran reader additionally maps each target `y(iy)` to the nearest source
row on the source uniform grid.  A stretched target grid is therefore also
supported.  Both readers implement the component mapping required by section
1; the resampling operation is specific to the Fortran restart reader.

## 4. Verification result (run/Dati.cart.in.Re2000, deck run/hst.in)

The check `python3 tests/verify_restart_convention.py run/Dati.cart.in.Re2000
run/hst.in` prints:

```
mesh: deck (nx, nz, ny+1) = (95, 47, 192)   file (nx, ny_cpl, nz_cpl) = (95, 47, 192)
hst sizes: ny(shear)=191, nz(spanwise)=47, nx=95
physical field shape: (191, 95, 96, 3)
header: time=2032 S=1 S2=0 gamma_x=2032 gamma_y=0
round-trip (physical vs raw CPL): exact
<u_i u_i> (abs^2 sum): file 7.94032029  field 7.94032029  match=yes
PASSED: convention conversion verified (no change to restart bytes)
```

- **Mesh** matches → the Fortran reader accepts the file against the current
  deck.
- **Shape** `(191, 95, 96, 3)` = `(ny_shear=191, 2*nz_span+1=95, nx+1=96, 3)`.
- **Round-trip exact**: reversing the `v_cpl↔w_cpl` swap from the physical
  field reproduces the raw CPL component arrays to machine precision — the
  conversion is an exact permutation, with no data lost or reordered.
- **Scalars survive** `time, S, S2, gamma_x, gamma_y` (all are plain 8-byte
  reals; neither is a transverse velocity component, so no component
  adaptation applies — here `S2=gamma_y=0` because no unsteady shear is on).
- **Energy invariant**: the `v↔w` swap is a rotation in the transverse
  velocity plane, so `<u_i u_i>` is unchanged (7.94032029).

## 5. Conclusion

`Dati.cart.in.Re2000` is already stored in the convention-neutral CPL layout
and the shipped readers (`src/hst_io.f90::restart_read` and
`post/python/io.py::read_field`) map its axes and swap its velocity components
into the new hst convention correctly.  **No rewrite of the restart file is
required** — changing the file's bytes would make the reader double-apply the
swap and corrupt the field.  `restart_read` now also supports resolution
changes using the simple Fourier-copy/nearest-neighbour strategy described
above.  Increasing resolution zero-pads new Fourier modes; decreasing it
truncates modes outside the target range.

The verification test added to the repo (`tests/verify_restart_convention.py`)
is the executable record of these checks.
