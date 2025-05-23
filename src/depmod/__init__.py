"""
Main module of the DEPMOD package.
"""

from __future__ import annotations

import numpy as np

from depmod._lib import read_lattice_from_file
from depmod.deformation import Compression, Deformation, DeformationPath, PureShear, Traction


def read_lattice(file: str, format: str = "", compression: str = ""):
    lattice = np.zeros((3, 3), dtype=float)
    try:
        read_lattice_from_file(lattice, str(file), format, compression)
        return lattice
    except RuntimeError as err:
        raise RuntimeError("Fail to read the file") from err


__version__ = "1.0.0"
__all__ = ["DeformationPath", "Deformation", "Traction", "Compression", "PureShear", "read_lattice_from_file"]
