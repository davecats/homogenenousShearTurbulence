#!/usr/bin/env python3
"""Compare two hst field files (Dati.cart*.out).   usage: compare_fields.py A B

File layout (see src/hst_io.f90): 3 int32 (nx ny nz), 6 float64
(alfa0 beta0 ly re S time), then complex128 array (ny, 2nz+1, nx+1, 3) in
Fortran order.  Prints the largest absolute and relative difference.
"""
import sys
import numpy as np


def read(fname):
    with open(fname, 'rb') as f:
        nx, ny, nz = np.fromfile(f, np.int32, 3)
        hdr = np.fromfile(f, np.float64, 6)
        a = np.fromfile(f, np.complex128).reshape((3, nx + 1, 2*nz + 1, ny))
    return hdr, a


ha, a = read(sys.argv[1])
hb, b = read(sys.argv[2])
print(f"time {ha[5]:.6f} vs {hb[5]:.6f}")
d = np.abs(a - b).max()
m = np.abs(a).max()
print(f"max |A-B| = {d:.3e}   max |A| = {m:.3e}   relative {d/m:.3e}")
