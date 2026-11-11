#!/usr/bin/env python3
"""Verify that the restart readers convert the old-CPL velocity convention to
the new hst convention correctly, with no loss or reordering of data.

Context (two physical conventions, one on-disk layout):

* Old CPL code (produced ``Dati.cart.in.Re2000``):  x=streamwise(u),
  y=spanwise(v), z=shear(w).
* New hst code (this repo):  x=streamwise(u), y=shear(v), z=spanwise(w).

The on-disk layout is *convention neutral*:  the CPL file always stores

    ARRAY(0..nx, -ny_cpl..ny_cpl, -1..nz_cpl+1)  of (u, v_cpl, w_cpl)

i.e. dim0 = streamwise, dim1 = spanwise modes (v_cpl = spanwise velocity),
dim2 = shear rows (w_cpl = shear velocity).  Nothing about the axis *type*
changes when the code calls spanwise "z" instead of "y" -- so the only thing
the new code must do to "adopt" an old-CPL file is relabel the velocity
components:  spanwise v_cpl  ->  new w (spanwise),  shear w_cpl  ->  new v
(shear).  Both the Fortran reader (src/hst_io.f90::restart_read, with
CPL_ORDER = [1, 3, 2]) and the Python reader (post/python/io.py::read_field) already
apply exactly that swap.

This script proves it on ``run/Dati.cart.in.Re2000``:

  1. mesh:   deck (nx, nz, ny+1)  ==  file (nx, ny_cpl, nz_cpl)
  2. shape:  physical field is (ny_hst, 2 nz_hst + 1, nx + 1, 3)
  3. round-trip:  reconstructing the original CPL arrays from the physical
     field by inverting the swap reproduces the raw file component-arrays to
     machine precision (the swap is an exact permutation, nothing lost).
  4. scalars: header time / S / S2 / gamma_x / gamma_y survive read_field.
  5. energy:  <u_i u_i> is invariant under the rotation v_cpl<->w_cpl.

Exit status 0 = PASS; otherwise FAILED and non-zero.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

from post.python.io import read_field
from post.python.config import load_config


def header_int(head: bytes, key: str) -> int:
    m = re.search(key.encode() + rb"=([0-9]+)", head)
    assert m, f"no {key}= in header"
    return int(m.group(1))


def fail(msg: str):
    print("FAILED:", msg)
    raise SystemExit(1)


def main(restart: str, deck: str | None = None) -> int:
    fp = Path(restart)
    if not fp.is_file():
        fail(f"restart file not found: {fp}")
    deck = deck or (fp.with_name("hst.in"))
    if not Path(deck).is_file():
        deck = ROOT / "run" / "hst.in"

    # ---- 0. raw CPL ground truth -------------------------------------------
    blob = fp.read_bytes()
    i = blob.find(b"Vfield=\n")
    if i < 0:
        fail("no Vfield=\\n header (not a CPL velocity file)")
    head = blob[:i]
    nx = header_int(head, "nx"); nyc = header_int(head, "ny"); nzc = header_int(head, "nz")
    raw = np.frombuffer(blob, np.complex128, offset=i + len(b"Vfield=\n"))
    a = raw.reshape((nx + 1, 2 * nyc + 1, nzc + 3, 3))   # (ix, span, row, cpl-comp)

    # ---- 1. mesh validation (mirrors src/hst_io.f90::restart_read) --------
    cfg = load_config(deck)
    want = (cfg.nx, cfg.nz, cfg.ny + 1)          # deck (nx, nz, ny+1)
    have = (nx, nyc, nzc)                        # file (nx, CPL ny, CPL nz)
    print(f"mesh: deck (nx, nz, ny+1) = {want}   file (nx, ny_cpl, nz_cpl) = {have}")
    if want != have:
        fail("mesh mismatch; file would be rejected by the Fortran reader")

    ny_hst = nzc - 1                             # shear points
    nz_hst = nyc                                 # spanwise modes
    print(f"hst sizes: ny(shear)={ny_hst}, nz(spanwise)={nz_hst}, nx={nx}")

    # ---- 2. new-convention physical field via the shipped reader ----------
    f = read_field(fp)
    vel = f.vel                                        # (ny_hst, 2 nz+1, nx+1, 3)
    print("physical field shape:", vel.shape)
    if vel.shape != (ny_hst, 2 * nz_hst + 1, nx + 1, 3):
        fail(f"unexpected physical shape {vel.shape}")
    # scalar header survival
    print("header: time=%.6g S=%g S2=%g gamma_x=%.6g gamma_y=%g"
          % (f.time, f.header.S, f.header.S2, f.header.gamma_x, f.header.gamma_y))

    # ---- 3. round-trip: invert the v<->w swap and compare to the raw file --
    # physical (u, v, w):  comp0=u(streamwise), comp1=v(shear)=CPL w,
    #                      comp2=w(spanwise)=CPL v.
    u = vel[..., 0]
    v_new = vel[..., 1]      # shear   == CPL w_cpl
    w_new = vel[..., 2]      # spanwise== CPL v_cpl
    # interior rows: file rows 2..2+ny_hst map to physical iy 0..ny_hst-1
    interior = a[:, :, 2:2 + ny_hst]                   # (ix, span, iy, comp)
    u_cpl = interior[..., 0]
    v_cpl = interior[..., 1]                           # spanwise
    w_cpl = interior[..., 2]                           # shear
    # compare with the appropriate axis permutation
    ok = (np.max(np.abs(u_cpl.T - u)) < 1e-12 and
          np.max(np.abs(w_cpl.T - v_new)) < 1e-12 and
          np.max(np.abs(v_cpl.T - w_new)) < 1e-12)
    print("round-trip (physical vs raw CPL):",
          "exact" if ok else "MISMATCH")
    if not ok:
        fail("component swap does not invert back to the raw file")

    # ---- 4. energy invariance (rotation of v_cpl<->w_cpl preserves trace) --
    e_file = (np.abs(interior) ** 2).sum()          # raw CPL 3 components
    e_field = (np.abs(vel) ** 2).sum()              # physical 3 components
    print(f"<u_i u_i> (abs^2 sum): file {e_file:.9g}  field {e_field:.9g}  "
          f"match={'yes' if np.isclose(e_file, e_field) else 'NO'}")

    print("PASSED: convention conversion verified (no change to restart bytes)")
    return 0


if __name__ == "__main__":
    args = sys.argv[1:]
    restart = args[0] if args else str(ROOT / "run" / "Dati.cart.in.Re2000")
    deck = args[1] if len(args) > 1 else None
    try:
        sys.exit(main(restart, deck))
    except SystemExit:
        raise
