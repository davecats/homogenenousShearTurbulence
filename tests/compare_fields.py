#!/usr/bin/env python3
"""Read and compare hst / CPL field files.

   compare_fields.py A [B]

Velocity files (Dati.cart.out, fields/field<n>.fld): CPL text header up to
"Vfield=\\n", then the array (0..nx, -ny_cpl..ny_cpl, -1..nz_cpl+1) of
(u, v, w) complex128 in C order, CPL naming (ny_cpl spanwise modes, nz_cpl-1
unique vertical points, their (v, w) = hst (w, v)).  Pressure files
(p_fields/pField<n>.fld) are the same array without header and with one
component; pass the sizes via the environment if needed (HST_NX, HST_NY,
HST_NZ in CPL naming) or a velocity file first, whose header sets them.
With one argument prints the header and the energy; with two, the largest
difference of the arrays (interior rows).
"""
import os
import sys
import numpy as np

sizes = {}


def read(fname):
    b = open(fname, 'rb').read()
    i = b.find(b'Vfield=\n')
    if i >= 0:
        head = b[:i].decode('latin-1')
        for key in ('nx', 'ny', 'nz'):
            j = head.find(key + '=')
            sizes[key] = int(head[j + len(key) + 1:].split()[0])
        t = np.frombuffer(b[b.find(b'time=\n') + 6:][:8], np.float64)[0]
        nx, ny, nz = sizes['nx'], sizes['ny'], sizes['nz']
        a = np.frombuffer(b[i + 8:], np.complex128).reshape((nx + 1, 2*ny + 1, nz + 3, 3))
        return t, a
    for key in ('nx', 'ny', 'nz'):
        if key not in sizes:
            sizes[key] = int(os.environ['HST_' + key.upper()])
    nx, ny, nz = sizes['nx'], sizes['ny'], sizes['nz']
    a = np.frombuffer(b, np.complex128).reshape((nx + 1, 2*ny + 1, nz + 3, 1))
    return None, a


def energy(a):
    w = np.where(np.arange(a.shape[0]) == 0, 1, 2)[:, None, None, None]
    return (np.abs(a[:, :, 2:-2, :])**2*w).sum()/(a.shape[2] - 4)


ta, a = read(sys.argv[1])
print(f"{sys.argv[1]}: CPL sizes nx={sizes['nx']} ny={sizes['ny']} nz={sizes['nz']}, time {ta}, <u_i u_i> = {energy(a):.9g}")
if len(sys.argv) > 2:
    tb, b = read(sys.argv[2])
    d = np.abs(a[:, :, 2:-2] - b[:, :, 2:-2]).max()
    m = np.abs(a).max()
    print(f"time {ta} vs {tb}:  max |A-B| = {d:.3e}   max |A| = {m:.3e}   relative {d/m:.3e}")
