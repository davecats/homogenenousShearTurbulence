"""Postprocessing tools for the `hst` homogeneous-shear-turbulence DNS.

Everything reads the hst / CPL output files through :mod:`post.python.io` (the
single source of truth for the binary layouts) and builds physical coordinates /
wavenumbers through :mod:`post.python.grid`.  Run tools from the command line
with::

    python3 -m post.python <tool> [deck] [args]

``<deck>`` is the run directory or an ``hst.in``-style namelist path (default
``hst.in`` in the current directory); it is parsed for the mesh, box and
physics metadata that the binary field headers do not always carry (the
pressure files are headerless, and the grid coordinates are reconstructed
from ``ny``, ``ly`` and ``ystretch``).

Adding a new tool is adding one file in :mod:`post.python.tools` and one line in
its ``TOOLS`` registry (see ``post/python/tools/__init__.py``).
"""

from .io import read_field, read_pressure, read_runtimedata, read_variances  # noqa: F401
from .grid import make_grid, make_wavenumbers  # noqa: F401
from .config import load_config  # noqa: F401

__version__ = "0.1.0"
