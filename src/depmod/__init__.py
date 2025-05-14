"""
Main module of the DEPMOD package.
"""
from __future__ import annotations

import numpy as np

from ._lib import read_lattice_from_file
from .deformation import DeformationPath, Deformation, Traction, Compression, PureShear

__version__ = "1.0.0"
__all__ = ["DeformationPath", "Deformation", "Traction", "Compression", "PureShear", ]


def read_lattice(file: str, format: str = "", compression: str = ""):
    lattice = np.zeros((3, 3), dtype=float)
    try:
        read_lattice_from_file(lattice, file, format, compression)
        return lattice
    except RuntimeError as err:
        raise RuntimeError("Fail to read the file") from err
