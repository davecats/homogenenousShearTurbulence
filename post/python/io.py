"""Readers for the ``hst`` / CPL output files (single source of truth).

Two binary layouts, both little-endian native ``complex128`` (``C_DOUBLE_COMPLEX``)
as written with MPI-IO:

* **Velocity** (restart ``Dati.cart.out``, snapshots ``fields/field<n>.fld``):
  a CPL text header up to a lone ``Vfield=`` line, then the array
  ``ARRAY(0..nx, -ny_cpl..ny_cpl, -1..nz_cpl+1)`` of ``(u, v, w)`` with the
  last index fastest (C order) and the components innermost.  CPL naming:
  ``ny_cpl = nz`` (spanwise), ``nz_cpl = ny + 1``; their ``(v, w)`` is our
  ``(w, v)`` and their row ``iz_cpl`` is our ``iy = iz_cpl - 1``.

* **Pressure** (``p_fields/pField<n>.fld``): the same array without header
  and with a single component.

:func:`read_field` returns a :class:`Field` whose :attr:`vel` is the *physical*
array ``shape (ny, 2*nz+1, nx+1, 3)`` with components ``(u, v, w)`` (the CPL
``v↔w`` swap resolved) and the ghost rows dropped (interior rows only).
:func:`read_pressure` returns the same ``(iy, iz, ix)`` axes for the pressure.

Everything here is deliberately light: tools that only read the ASCII
statistics never load a snapshot, and tools that read many snapshots reuse the
parsed header.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np


def _find_header(binary: bytes):
    """Return ``(header_bytes, payload_offset, header_dict)`` for a velocity file."""
    i = binary.find(b"Vfield=\n")
    if i < 0:
        return b"", 0, {}
    head = binary[:i].decode("latin-1")
    hdr = {}
    for key in ("nx", "ny", "nz"):
        m = re.search(rf"{key}=([^\s\t\n]+)", head)
        if m:
            hdr[key] = int(m.group(1))
    for key in ("alpha0", "beta0", "htcoeff", "Re", "Pr", "deltat", "t_max",
                "dt_field", "dt_save", "t_field", "meanpx", "meanpy"):
        m = re.search(rf"{key}=([^\s\t\n]+)", head)
        if m:
            hdr[key] = float(m.group(1))
    return head.encode("latin-1"), i + len("Vfield=\n"), hdr


def _parse_header_raw(head_bytes: bytes) -> dict:
    """Extract the trailing raw-8-byte doubles (time, S, S2, gamma_x, gamma_y)."""
    # They live right before the final 'Vfield=\n', each on its own line set.
    out = {}
    keys = ("time", "S", "S2", "gamma_x", "gamma_y")
    # We locate each keyword as b'key=\n<8 raw bytes>\n'.
    text = head_bytes
    for key in keys:
        pat = key.encode() + b"=\n"
        j = text.rfind(pat)
        if j >= 0:
            k = j + len(pat)
            if k + 8 <= len(text):
                out[key] = np.frombuffer(text[k:k + 8], np.float64)[0]
    return out



@dataclass
class Header:
    """Parsed velocity-file header (CPL naming: ny = spanwise, nz = ny+1)."""

    nx: int
    ny: int
    nz: int
    alpha0: float
    beta0: float
    re: float = 0.0
    htcoeff: float = -1.0
    time: float = 0.0
    S: float = 0.0
    S2: float = 0.0
    gamma_x: float = 0.0
    gamma_y: float = 0.0
    t_field: float = 0.0

    @property
    def ny_hst(self) -> int:
        """hst vertical point count = nz(cpl) - 1."""
        return self.nz - 1

    @property
    def nz_hst(self) -> int:
        """hst spanwise mode count = ny(cpl)."""
        return self.ny


@dataclass
class Field:
    """A velocity snapshot: physical interior values plus the parsed header.

    ``vel`` has shape ``(ny_hst, 2*nz_hst+1, nx_hst+1, 3)`` indexed
    ``[iy][iz_index][ix][comp]``, comp ``0/1/2`` = physical ``u/v/w``.  Use
    :func:`post.grid.make_wavenumbers` to map ``iz_index -> iz = iz_index -
    nz_hst`` and ``[ix] -> kx``.
    """

    path: Path
    vel: np.ndarray            # (ny, 2*nz+1, nx+1, 3)
    header: Header
    time: float = 0.0


def read_field(path: str | Path) -> Field:
    """Read a velocity field file into a physical ``(u, v, w)`` array."""
    path = Path(path)
    b = path.read_bytes()
    head_bytes, off, hdr = _find_header(b)
    if off == 0:
        raise ValueError(
            f"{path} has no 'Vfield=' header (not a velocity field); "
            "use read_pressure for headerless pressure files."
        )
    nx = int(hdr["nx"]); nyc = int(hdr["ny"]); nzc = int(hdr["nz"])
    ny_hst = nzc - 1
    nz_hst = nyc

    raw = np.frombuffer(b, dtype=np.complex128, offset=off)
    a = raw.reshape((nx + 1, 2 * nyc + 1, nzc + 3, 3))   # (ix, iz_cpl, row, comp)

    # The stored array is (nx+1, 2*nyc+1, ny+4, 3) in CPL order, where the
    # four ghost rows are the lowest two and highest two positions.  The
    # interior physical rows iy = 0..ny_hst-1 occupy file rows 2..ny+1
    # (0-based), as the writer's cpl_view maps row r -> slot r+2 and
    # compare_fields.py reads them with a[:, :, 2:-2].
    p = np.transpose(a, (2, 1, 0, 3))     # (file_row, iz_cpl, ix, comp)
    interior = p[2:2 + ny_hst]

    # CPL (u, v, w) -> physical (u, v_vertical, w_spanwise): CPL v is the
    # spanwise (hst w), CPL w is the vertical (hst v).  We expose physical
    # (u, v, w) = component 0,1,2 = (CPL u, CPL w, CPL v).
    u = interior[..., 0]
    v = interior[..., 2]
    w = interior[..., 1]
    vel = np.stack([u, v, w], axis=-1)

    r = _parse_header_raw(head_bytes)
    h = Header(
        nx=nx, ny=nyc, nz=nzc,
        alpha0=float(hdr.get("alpha0", 0.0)), beta0=float(hdr.get("beta0", 0.0)),
        re=float(hdr.get("Re", 0.0)), htcoeff=float(hdr.get("htcoeff", -1.0)),
        t_field=float(hdr.get("t_field", 0.0)),
        time=float(r.get("time", 0.0)), S=float(r.get("S", 0.0)),
        S2=float(r.get("S2", 0.0)), gamma_x=float(r.get("gamma_x", 0.0)),
        gamma_y=float(r.get("gamma_y", 0.0)),
    )
    return Field(path=path, vel=vel, header=h, time=h.time)


def read_pressure(path: str | Path, cfg=None, ref_field: Field | None = None,
                  sizes: dict | None = None) -> np.ndarray:
    """Read a headerless pressure file (``p_fields/pField<n>.fld``).

    Returns the pressure with the same ``(iy, iz_index, ix)`` axes as
    :func:`read_field`.vel.  Sizes come from ``cfg`` (a run config), else a
    ``ref_field`` header, else a ``sizes`` dict keys ``nx/ny/nz`` in CPL
    naming, else the environment ``HST_NX/HST_NY/HST_NZ`` (CPL naming).
    """
    path = Path(path)
    b = path.read_bytes()
    _, off, _ = _find_header(b)   # pressure files: off == 0 (headerless)

    if cfg is not None:
        nx, ny_hst, nz_hst = int(cfg.nx), int(cfg.ny), int(cfg.nz)
    elif ref_field is not None:
        nx, ny_hst, nz_hst = ref_field.header.nx, ref_field.header.ny_hst, ref_field.header.nz_hst
    elif sizes:
        nx, ny_hst, nz_hst = int(sizes["nx"]), int(sizes["nz"]) - 1, int(sizes["ny"])
    else:
        nx = int(os.environ["HST_NX"])
        nyc = int(os.environ["HST_NY"]); nzc = int(os.environ["HST_NZ"])
        ny_hst, nz_hst = nzc - 1, nyc

    nyc = nz_hst                      # CPL spanwise count
    nzc = ny_hst + 1                  # CPL row count
    arr = np.frombuffer(b, dtype=np.complex128, offset=off)
    a = arr.reshape((nx + 1, 2 * nyc + 1, nzc + 3))     # (ix, iz_cpl, row)
    p = np.transpose(a, (2, 1, 0))     # (file_row, iz_cpl, ix)
    return p[2:2 + ny_hst]             # interior rows (file rows 2..ny+1)


# gfortran drops the "E" of three-digit exponents (1.2E-100 is written "1.2-100"),
# which happens when a quantity decays towards zero (e.g. a relaminarized run).
_MISSING_EXPONENT_E = re.compile(r"(?<=\d)([+-]\d{3})(?=\s|$)", re.M)


def read_text_repaired(path: str | Path) -> str:
    """Return the text of an ASCII statistics file with dropped exponent letters restored."""
    return _MISSING_EXPONENT_E.sub(r"E\1", Path(path).read_text())


def _loadtxt(path: str | Path):
    return np.atleast_2d(np.loadtxt(read_text_repaired(path).splitlines(), comments="#")).astype(float)


def find_run_file(run_dir: str | Path, name: str) -> Path:
    """A run-time statistics file of a run: ``<run_dir>/statistics/<name>`` (where
    ``hst`` writes it), else ``<run_dir>/<name>`` (runs made before the
    ``statistics/`` directory existed)."""
    run_dir = Path(run_dir)
    for candidate in (run_dir / "statistics" / name, run_dir / name):
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(
        f"could not find {name}; checked {run_dir / 'statistics' / name} and {run_dir / name}"
    )


def read_runtimedata(path: str | Path):
    """Read ``Runtimedata``: one line per ``dt_stat``, the 13 columns.

    Returns a structured numpy array with fields ``time, meanflowx,
    meanflowy, S, S2, gamma_x, gamma_y, deltat, cfl_dt, energy, diss, uv,
    vw``.  Naming follows the solver's (streamwise, shearwise, spanwise) =
    (u, v, w) convention.  All fluctuation quantities are y-averaged
    (box-mean) values built by Parseval's theorem along the spectral
    (streamwise, spanwise) directions:
    ``energy = <u_i u_i>/2`` (the TKE ``k``), ``diss = nu <∂_j u_i ∂_j u_i>``
    (kinematic viscosity ``nu = 1/re`` already folded in),
    ``uv = <u v>`` and ``vw = <v w>`` (the Reynolds shear stresses).
    ``meanflowx/meanflowy`` remain the box-height integrals of the mean
    streamwise / spanwise profiles (as in CPL).
    """
    names = ("time", "meanflowx", "meanflowy", "S", "S2",
             "gamma_x", "gamma_y", "deltat", "cfl_dt",
             "energy", "diss", "uv", "vw")
    data = _loadtxt(path)
    if data.shape[1] == len(names):
        return data.view(dtype=[(n, float) for n in names]).reshape(-1)
    return data


def read_variances(path: str | Path):
    """Read ``variances_runtime.dat``: ``time, uu, vv, ww, uw``.

    The four velocity second moments are y-averaged (box-mean) values:
    ``uu = <u u>``, ``vv = <v v>``, ``ww = <w w>`` (mean normal Reynolds
    stresses) and ``uw = <u w>`` (the streamwise-spanwise Reynolds shear
    stress; note this is ``<u w>``, not the ``<u v>`` logged in
    ``Runtimedata``).
    """
    names = ("time", "uu", "vv", "ww", "uw")
    data = _loadtxt(path)
    if data.shape[1] == len(names):
        return data.view(dtype=[(n, float) for n in names]).reshape(-1)
    return data


def read_stokes_runtime(path: str | Path):
    """Read ``stokes_runtime.dat``: ``time, energy_out, energy_in, diss_out, diss_in``."""
    names = ("time", "energy_out", "energy_in", "diss_out", "diss_in")
    data = _loadtxt(path)
    if data.shape[1] == len(names):
        return data.view(dtype=[(n, float) for n in names]).reshape(-1)
    return data
