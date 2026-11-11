"""Reconstruct physical coordinates and wavenumbers for a run.

The binary field files carry no grid data, only the header values ``nx, ny,
nz, alpha0, beta0`` and the deck's ``ly`` and ``ystretch`` (as ``htcoeff``
when it is non-zero).  This module rebuilds everything a downstream tool needs
to place a field in physical or spectral space, matching the solver exactly:

* the vertical grid :math:`y`, uniform (:math:`i*ly/ny`) or clustered at
  mid-box by the two-sided tanh map ``ymap`` of ``hst_setup.f90``;
* the row spacing ``dyl`` used for box integrals; and
* the wavenumber arrays ``kx``, ``kz`` and the squared magnitude ``k2`` with
  the mode-counting weight the solver applies (``w = 1`` on the ``ix = 0``
  plane, ``w = 2`` elsewhere, the conjugate modes only stored once).

Axis/component convention matches :func:`post.io.read_field`, so ``u, v, w``
are the physical streamwise, vertical and spanwise fluctuations with
:math:`y` the vertical coordinate.
"""

from __future__ import annotations

import numpy as np

from .config import Config


def ymap(xi: np.ndarray, ystretch: float) -> np.ndarray:
    """xi in [0,1) -> [0,1); identical to hst_setup.f90::ymap."""
    xi = np.asarray(xi, dtype=float)
    if ystretch <= 1.0e-10:
        return xi
    out = np.empty_like(xi)
    low = xi <= 0.5
    out[low] = 0.5 * np.tanh(2.0 * ystretch * xi[low]) / np.tanh(ystretch)
    out[~low] = 0.5 * (2.0 + np.tanh(2.0 * ystretch * (xi[~low] - 1.0)) / np.tanh(ystretch))
    return out


def make_grid(cfg: Config) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(y, dyl)`` for the unique rows ``0..ny-1``.

    ``y`` is length ``ny`` (vertical coordinate over the box of height
    ``cfg.ly``, matches ``hst_setup``).  ``dyl`` is the row spacing used for
    box integrals: ``0.5*(y(iy+1) - y(iy-1))``.
    """
    ny = int(cfg.ny)
    iy = np.arange(ny)
    # y(iy) = ly*ymap(iy/ny) + ly*0  (the wrap term is zero for interior rows)
    y = cfg.ly * ymap(iy / ny, cfg.ystretch)
    # periodic image for the +1 row:
    y_next = cfg.ly * ymap(((iy + 1) % ny) / ny, cfg.ystretch) + cfg.ly * ((iy + 1) // ny)
    y_prev = cfg.ly * ymap(((iy - 1) % ny) / ny, cfg.ystretch) + cfg.ly * ((iy - 1) // ny)
    dyl = 0.5 * (y_next - y_prev)
    return y, dyl


def make_wavenumbers(cfg: Config):
    """Return ``(kx, kz, k2, w)`` arrays over the stored modes.

    ``kx`` has length ``nx+1`` for modes ``ix = 0..nx``, ``kz`` length
    ``2*nz+1`` for ``iz = -nz..nz`` (real ``kz = beta0*iz``; the solver's
    ``ialfa``/``ibeta`` are these times ``i`` and are recovered by tools that
    need the phase).  ``k2`` is shape ``(nx+1, 2*nz+1)`` in ``(ix, iz)``
    order.  ``w`` is the mode-counting weight of the solver's statistics
    (``(nx+1, 1)``), 1 on ``ix = 0`` and 2 elsewhere.
    """
    nx, nz = int(cfg.nx), int(cfg.nz)
    ix = np.arange(nx + 1)
    iz = np.arange(-nz, nz + 1)
    kx = cfg.alfa0 * ix
    kz = cfg.beta0 * iz
    k2 = kx[:, None] ** 2 + kz[None, :] ** 2
    w = np.where(ix == 0, 1.0, 2.0)[:, None]
    return kx, kz, k2, w
