"""Extract physical-space slices from a velocity / pressure snapshot.

Usage::

    python3 -m post.python slices [deck] field1.fld [-z iz | -y j] [-c u|v|w|p] \\
        [-o out.png | --npz]

Reads a velocity field (or pressure) and produces a slice:

* ``-z iz``  an ``x-y`` plane at a fixed spanwise index ``iz`` (0..2*nz),
* ``-y j``   a ``y-z`` plane at a fixed vertical index ``j`` (0..ny-1).

The default is a mid-box ``x-y`` plane of streamwise velocity ``u``.  With
``--npz`` the slice is written as a numpy ``.npz`` with the physical
coordinates, so it can be reused by later tools without re-reading the binary.
"""

from __future__ import annotations

import dataclasses
import argparse
from pathlib import Path

import numpy as np

from ..io import read_field, read_pressure
from ..grid import make_grid, make_wavenumbers


def _override_cfg(cfg, **kw):
    """Return a shallow copy of ``cfg`` with the given attributes overridden.

    Used so tools can prefer the dimensions a binary field file carries in its
    header (the authoritative source) over the deck when the two disagree --
    e.g. when the deck defaulted and the snapshot is from a different run.
    ``dataclasses.replace`` would copy, but ``Config`` may grow; a manual
    shallow copy of the dataclass fields is enough.
    """
    return dataclasses.replace(cfg, **kw)


def run(cfg, argv) -> int:
    p = argparse.ArgumentParser(prog="slices", description="Slice a snapshot.")
    p.add_argument("field", nargs="?", default=None,
                   help="path to a fields/field<n>.fld (or a Dati.cart.out)")
    p.add_argument("-z", dest="iz", type=int, default=None,
                   help="spanwise index (0..2*nz) for an x-y plane")
    p.add_argument("-y", dest="jy", type=int, default=None,
                   help="vertical row (0..ny-1) for a y-z plane")
    p.add_argument("-c", "--component", default="u", choices=["u", "v", "w", "p"])
    p.add_argument("-o", "--out", default=None, help="output image (.png) or .npz")
    p.add_argument("--npz", action="store_true",
                   help="always write a .npz (several snapshots / later reuse)")
    args = p.parse_args(argv)

    if args.field is None:
        # pick the first snapshot available
        fields_dir = cfg.run_dir / "fields"
        snaps = sorted(fields_dir.glob("field*.fld"))
        if snaps:
            args.field = str(snaps[0])
        elif (cfg.run_dir / "Dati.cart.out").exists():
            args.field = str(cfg.run_dir / "Dati.cart.out")
        else:
            raise FileNotFoundError("no snapshot and no Dati.cart.out found")

    comp_i = {"u": 0, "v": 1, "w": 2}.get(args.component)
    if comp_i is None and args.component != "p":
        raise SystemExit("unknown component")

    path = Path(args.field)
    if comp_i is not None:
        fld = read_field(path)
        vel = fld.vel
        ny, nz2, nx1 = vel.shape[:3]
        # Velocity files carry their own authoritative dimensions / wavenumber
        # basis in the header, so build the grid from the header (overriding
        # any deck) rather than from the config -- the binary data always wins.
        h = fld.header
        eff = _override_cfg(cfg, nx=h.nx, ny=h.ny_hst, nz=h.nz_hst,
                            alfa0=h.alpha0, beta0=h.beta0,
                            ystretch=h.htcoeff if h.htcoeff > 0 else cfg.ystretch)
        y, dyl = make_grid(eff)
        kx, kz, k2, w = make_wavenumbers(eff)
        zcoord = np.arange(-eff.nz, eff.nz + 1) * eff.beta0
    else:
        # Pressure files are headerless: the dimensions come from the deck.
        fld = read_pressure(path, cfg=cfg)
        vel = fld.reshape(fld.shape + (1,))
        ny, nz2, nx1 = fld.shape
        kx, kz, k2, w = make_wavenumbers(cfg)
        y, dyl = make_grid(cfg)
        zcoord = np.arange(-cfg.nz, cfg.nz + 1) * cfg.beta0

    if args.iz is not None and args.jy is not None:
        raise SystemExit("choose at most one of -z or -y")
    if args.iz is None and args.jy is None:
        args.iz = nz2 // 2          # mid-box plane (default)

    if args.iz is not None:
        iz = args.iz
        data = vel[:, iz, :, comp_i if comp_i is not None else 0].T
        xs = kx
        ys = y
    else:
        jy = args.jy
        data = vel[jy, :, :, comp_i if comp_i is not None else 0].T
        xs = zcoord
        ys = y

    if args.npz or (args.out and args.out.endswith(".npz")):
        out = Path(args.out or (cfg.run_dir / f"slice_{args.component}.npz"))
        np.savez(out, x=xs, y=ys, data=data,
                 component=args.component, iz=args.iz, jy=args.jy,
                 time=getattr(fld, "time", 0.0))
        print(f"wrote {out}  ({data.shape})")
        return 0

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    X, Y = np.meshgrid(xs, ys, indexing="ij")
    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.pcolormesh(X, Y, data.real, shading="auto", cmap="RdBu_r")
    fig.colorbar(im, ax=ax, label=args.component)
    ax.set_xlabel("x (streamwise)" if args.iz is not None else "z (spanwise)")
    ax.set_ylabel("y (vertical)")
    plane = "x-y" if args.iz is not None else "y-z"
    ax.set_title(f"{path.name}  {plane}  {args.component}")
    fig.tight_layout()
    out = Path(args.out or (cfg.run_dir / f"slice_{args.component}.png"))
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")
    return 0
