"""Module to handle input parameters."""
from __future__ import annotations

from depmod._lib.config import _lib_Configuration


Config = _lib_Configuration


def init_config(
    gammadot: float,
    t_min: float = 0.,
    t_max: float | None = None,
    e_max: float | None = None,
    npts: int = 100,
    kpts: int = int | None,
    err_tresh: float = 1., # TODO: define err_tresh
) -> _lib_Configuration:
    
    gammadot = float(gammadot)
    npts = int(npts)
    t_min = float(t_min)

    if (t_max is not None) and (e_max is not None):
        raise ValueError("t_max and e_max can't be set at the same time")
    elif (t_max is not None) and (e_max is None):
        t_max = float(t_max)
        e_max = float(gammadot * t_max)
    elif (t_max is None) and (e_max is not None):
        e_max = float(e_max)
        t_max = float(e_max / gammadot)
    elif (t_max is None) and (e_max is None):
        raise ValueError("Either t_max or e_max must be set")

    if err_tresh is not None:
        k = int(((t_max - t_min) / (err_tresh * npts)) * gammadot ** 2)

    # overwrite iter
    if kpts is not None:
        k = kpts

    return _lib_Configuration(gammadot, t_min, t_max, npts, k)
