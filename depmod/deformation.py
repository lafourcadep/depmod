from __future__ import annotations

import numpy as np

from depmod.maths import dirvec, unitvec, check_vector_orthogonality
from depmod.typing import ArrayLike

from depmod._lib.deformation import (
    _lib_Deformation,
    _lib_TractionCompression,
    _lib_TractionCompressionIsoV,
    _lib_PureShear
)

Deformation = _lib_Deformation

class TractionCompression:
    def __init__(self, comp=False, isoV=False):
        self.comp = comp
        self.isoV = isoV

    def from_axis(self, axis):
        trac = unitvec(axis)
        S = np.tensordot(trac, trac, axes=0).astype(float)
        if self.isoV:
            return _lib_TractionCompressionIsoV(S, self.comp)
        else:
            return _lib_TractionCompression(S, self.comp)

    def from_angles(self, theta: float, phi: float, degree: bool = True):
        if degree:
            theta = np.deg2rad(theta)
            phi = np.deg2rad(phi)
        return self.from_axis(dirvec(theta, phi))


class Traction(TractionCompression):
    def __init__(self, isoV: bool = False):
        TractionCompression.__init__(self, comp=False, isoV=isoV)


class Compression(TractionCompression):
    def __init__(self, isoV: bool = False):
        TractionCompression.__init__(self, comp=True, isoV=isoV)


class PureShear:
    def __init__(self, strict: bool = True):
        self.strict = True

    def from_axis(self, eta: ArrayLike, kappa: ArrayLike):
        eta = unitvec(eta)
        kappa = unitvec(kappa)
        
        if self.strict and not check_vector_orthogonality(eta, kappa):
            raise ValueError(f"eta: {eta} and kappa: {kappa} are not orthogonal.")

        S = np.tensordot(eta, kappa, axes=0).astype(float)
        return _lib_PureShear(S)

    def from_angles(self, theta: float, phi: float, degree=True):
        if degree:
            theta = np.deg2rad(theta)
            phi = np.deg2rad(phi)
        eta = dirvec(theta, phi)
        kappa = dirvec(theta + np.pi / 2, phi)

        return self.from_axis(eta, kappa)


# def MixedDeformation(dlist: list[Deformation]) -> Deformation:
#     if not all(isinstance(it, Deformation) for it in dlist):
#         raise TypeError("Not a deformation object.")
#     return _lib_MixedDeformation(dlist)
