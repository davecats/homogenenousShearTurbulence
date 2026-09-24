#!/usr/bin/env python3
"""Convert an hst field file into the CPL hst-main .fld layout.

   to_cpl_field.py Dati.cart.out template.fld out.fld

template.fld is any field written by the CPL code with the same nx, ny, nz
(its text header, up to and including the "Vfield=\\n" line, is copied).
Index and component mapping (CPL (x, y_spanwise, z_vertical, u, v, w) <->
hst (x, z_spanwise, y_vertical, u, w, v)):
   cpl[ix, iy_cpl + ny_cpl, iz_cpl + 1, :]  with  iz_cpl = 1 .. nz_cpl-1  <->  hst iy = iz_cpl - 1
   ny_cpl = hst nz,  nz_cpl = hst ny + 1;  cpl component (u, v, w) = hst (u, w, v).
The CPL ghost rows iz_cpl = -1, 0, nz, nz+1 are filled periodically (the CPL
code re-applies its periodic condition at start-up anyway).
"""
import sys
import numpy as np

src, template, out = sys.argv[1:4]
with open(src, 'rb') as f:
    nx, ny, nz = np.fromfile(f, np.int32, 3)
    hdr = np.fromfile(f, np.float64, 6)
    a = np.fromfile(f, np.complex128).reshape((3, nx + 1, 2*nz + 1, ny))   # C-order view of the Fortran array

b = open(template, 'rb').read()
i = b.find(b'Vfield=\n')
head = b[:i + 8]
ny_cpl, nz_cpl = nz, ny + 1
for tag, val in ((b'nx=', nx), (b'ny=', ny_cpl), (b'nz=', nz_cpl)):
    j = head.find(tag)
    assert int(head[j + len(tag):].split()[0]) == val, f"template {tag} does not match"

cpl = np.zeros((nx + 1, 2*ny_cpl + 1, nz_cpl + 3, 3), np.complex128)
# interior rows: cpl index iz_cpl+1 for iz_cpl = 1..nz_cpl-1  ->  2 .. ny+1
inner = cpl[:, :, 2:2 + ny, :]
inner[:, :, :, 0] = np.transpose(a[0], (0, 1, 2))   # u
inner[:, :, :, 1] = np.transpose(a[2], (0, 1, 2))   # cpl v (spanwise) = hst w
inner[:, :, :, 2] = np.transpose(a[1], (0, 1, 2))   # cpl w (vertical) = hst v
# periodic ghosts: iz_cpl = nz_cpl (index ny+2) = iz_cpl 1 (index 2), etc.
cpl[:, :, ny + 2, :] = cpl[:, :, 2, :]
cpl[:, :, ny + 3, :] = cpl[:, :, 3, :]
cpl[:, :, 1, :] = cpl[:, :, ny + 1, :]
cpl[:, :, 0, :] = cpl[:, :, ny, :]
with open(out, 'wb') as f:
    f.write(head)
    cpl.tofile(f)
print(f"wrote {out}: header {len(head)} bytes, array {cpl.nbytes} bytes, time in hst file {hdr[5]}")
