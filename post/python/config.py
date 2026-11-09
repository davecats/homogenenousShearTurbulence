"""Resolve the run directory and parse the ``hst.in`` Fortran namelist.

The deck holds the mesh, box and physics metadata that the binary field files
do not always carry (the pressure files ``p_fields/pField<n>.fld`` are
headerless) and from which the physical grid and wavenumbers are rebuilt by
:mod:`post.grid`.

The files are written next to the deck (``fields/``, ``p_fields/``,
``Runtimedata``, ...), so given a deck path we also derive the run directory.

Fortran namelists are line-oriented ``name = value1, value2, ...`` with
``!`` comments; we only ever need the scalar keys used by ``hst``, so a small
tolerant reader is enough.  Values are parsed as Python scalars; Fortran
suffixes (``d0``, ``_true``) are stripped.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

# Keys we understand, per namelist group, with a python typecast.  All hst
# deck parameters are here so the config is a faithful record of a run even
# for keys no tool uses yet.
_GROUP_KEYS: dict[str, dict[str, str]] = {
    "mesh": {
        "nx": int, "ny": int, "nz": int,
        "alfa0": float, "beta0": float, "ly": float,
        "ystretch": float, "line_chunk": int,
        "transport": str, "npy": int,
    },
    "physics": {
        "re": float, "S": float, "linear": bool,
        "exact_shift": bool,
        "s2_mean": float, "s2_amplitude": float, "s2_period": float, "s2_start": float,
        "sl_amplitude": float, "sl_period": float, "sl_delta": float,
        "sl_start": float, "sl_bodyforce": bool, "sl_ramp": bool,
    },
    "time_control": {
        "deltat": float, "cflmax": float, "t_max": float, "nstep": int,
        "dt_stat": float, "dt_field": float, "dt_save": float,
        "time": float, "time_from_restart": bool, "timing": bool,
    },
    "init": {"amplitude": float, "seed": int, "kpeak": float},
}

# Defaults, kept pessimistically consistent with hst.in (deck values win).
_DEFAULTS: dict[str, dict[str, object]] = {
    "mesh": {"nx": 63, "ny": 128, "nz": 63,
             "alfa0": 2.0943951023931953, "beta0": 6.283185307179586,
             "ly": 2.0, "ystretch": 0.0, "line_chunk": 0,
             "transport": "auto", "npy": 1},
    "physics": {"re": 1000.0, "S": 1.0, "linear": False, "exact_shift": False,
                "s2_mean": 0.0, "s2_amplitude": 0.0, "s2_period": 0.0, "s2_start": 0.0,
                "sl_amplitude": 0.0, "sl_period": 1.0, "sl_delta": 0.02,
                "sl_start": 0.0, "sl_bodyforce": True, "sl_ramp": False},
    "time_control": {"deltat": 0.0, "cflmax": 1.0, "t_max": 100.0,
                     "nstep": 1000000, "dt_stat": 0.01, "dt_field": 10.0,
                     "dt_save": 10.0, "time": 0.0,
                     "time_from_restart": False, "timing": False},
    "init": {"amplitude": 1.0e-3, "seed": 1, "kpeak": 4.0},
}



def _strip_comment(text: str) -> str:
    # A ! inside an unquoted token starts a comment; hst decks quote strings
    # with ' so naive scan is fine here.
    out, quote = [], None
    for ch in text:
        if ch in ("'", '"'):
            if quote == ch:
                quote = None
            elif quote is None:
                quote = ch
        elif ch == "!" and quote is None:
            break
        out.append(ch)
    return "".join(out)


def _cast(value: str, kind: str):
    value = value.strip()
    if kind is float:
        # Fortran real suffixes: 0.0d0, 1.D0, 2.0_q8, or _rk from the deck
        value = re.sub(r"[dDqQ]\s*([+-]?\d+)\b", r"e\1", value)
        value = re.sub(r"_?\b[a-zA-Z_]\w*\s*$", "", value).strip()
    if kind is bool:
        return value.lower() in (".true.", "true", "t", "yes", "y", "1")
    if kind is int:
        # strip a trailing comma / any following 'name = val' on the line
        value = re.sub(r",.*$", "", value).strip()
        return int(value)
    if kind is float:
        value = re.sub(r",.*$", "", value).strip()
        return float(value)
    return value.strip().strip("'\"")


_LINE = re.compile(r"^[ \t]*([A-Za-z_][A-Za-z0-9_]*)[ \t]*=[ \t]*(.*)$")


def _split_assignments(line: str) -> list[tuple[str, str]]:
    """Split a namelist statement line like 'nx = 1, ny = 2, nz = 3' into pairs."""
    parts, buf, quote, depth = [], "", None, 0
    for ch in line:
        if ch in ("'", '"'):
            if quote == ch:
                quote = None
            elif quote is None:
                quote = ch
        elif ch in "([{" and quote is None:
            depth += 1
        elif ch in ")]}" and quote is None:
            depth -= 1
        if ch == "," and quote is None and depth == 0:
            parts.append(buf); buf = ""
        else:
            buf += ch
    parts.append(buf)
    out = []
    for part in parts:
        m = _LINE.match(part.strip())
        if m:
            out.append((m.group(1).lower(), m.group(2).strip()))
    return out


def _parse_namelist(text: str) -> dict[str, dict[str, object]]:
    """Return {group: {key: value}} from the whole deck text."""
    groups: dict[str, dict[str, object]] = {}
    current = None
    for raw in text.splitlines():
        raw = _strip_comment(raw)
        m = re.match(r"^\s*&([A-Za-z_][A-Za-z0-9_]*)\s*$", raw)
        if m:
            current = m.group(1).lower()
            groups.setdefault(current, {})
            continue
        if current is None:
            continue
        for key, val in _split_assignments(raw):
            if current in _GROUP_KEYS and key in _GROUP_KEYS[current]:
                groups[current][key] = _cast(val, _GROUP_KEYS[current][key])
    return groups


@dataclass
class Config:
    """Everything a tool needs to know about a run, resolved from a deck."""

    #: directory the deck lives in (where fields/, p_fields/, statistics/ are)
    run_dir: Path
    deck_path: Path
    # mesh
    nx: int = 63
    ny: int = 128
    nz: int = 63
    alfa0: float = 2.0943951023931953
    beta0: float = 6.283185307179586
    ly: float = 2.0
    ystretch: float = 0.0
    line_chunk: int = 0
    transport: str = "auto"
    npy: int = 1
    # physics
    re: float = 1000.0
    S: float = 1.0
    linear: bool = False
    # time control
    deltat: float = 0.0
    dt_stat: float = 0.01
    dt_field: float = 10.0
    dt_save: float = 10.0
    time: float = 0.0
    # everything else kept verbatim from the deck, for future tools
    extra: dict = field(default_factory=dict)

    @property
    def nu(self) -> float:
        """Kinematic viscosity 1/Re."""
        return 1.0 / self.re

    @property
    def lx(self) -> float:
        """Streamwise box length 2*pi/alfa0."""
        return 2.0 * 3.141592653589793 / self.alfa0

    @property
    def lz(self) -> float:
        """Spanwise box length 2*pi/beta0."""
        return 2.0 * 3.141592653589793 / self.beta0


def load_config(deck: str | Path | None = None) -> Config:
    """Build a :class:`Config` from an ``hst.in`` deck (or a run directory).

    ``deck`` defaults to ``hst.in`` in the current directory.  If it is a
    directory that contains ``hst.in``, that deck is used and the run
    directory is the directory itself.
    """
    if deck is None:
        path = Path("hst.in")
    else:
        path = Path(deck)
    if path.is_dir():
        path = path / "hst.in"
    if not path.is_file():
        raise FileNotFoundError(f"no deck found at {path} (pass a path or run dir)")

    text = path.read_text()
    groups = _parse_namelist(text)

    cfg = Config(run_dir=path.parent, deck_path=path)
    extra: dict[str, object] = {}
    for group, keys in _GROUP_KEYS.items():
        for key in keys:
            if group in groups and key in groups[group]:
                setattr(cfg, key, groups[group][key])
            else:
                setattr(cfg, key, _DEFAULTS[group][key])
        if group in groups:
            for key, value in groups[group].items():
                if key not in keys:
                    extra[f"{group}.{key}"] = value
    cfg.extra = extra
    return cfg

