from __future__ import annotations

import numpy as np

from ._lib import _Deformation, _DeformationPath
from .maths import check_vector_orthogonality, dirvec, spherical_basis, trim_array, unitvec


class Deformation(_Deformation):
    def __init__(self, S, cfac=1.0):
        _Deformation.__init__(self, S, cfac)


class TractionCompression(_Deformation):
    def __init__(self, S, comp: bool = False):
        cfac = 1.0 if not comp else -1.0
        _Deformation.__init__(self, S, cfac)

    @classmethod
    def from_axis(cls, eta, comp: bool = False, isoV: bool = False):
        eta = unitvec(eta)
        S = np.tensordot(eta, eta, axes=0).astype(float)

        if isoV:
            theta = np.arccos(eta[2])  # arccos(z)
            phi = np.arctan2(eta[1], eta[0])  # arctan2(y, x)
            m, n, p = spherical_basis(theta, phi, chi=0.0)
            assert np.allclose(m, eta)
            S -= 0.5 * (np.tensordot(n, n, axes=0) + np.tensordot(p, p, axes=0))

        return cls(trim_array(S, 1.0e-8), comp=comp)

    @classmethod
    def from_angles(cls, theta: float, phi: float, comp: bool = False, isoV: bool = False, degree: bool = False):
        if degree:
            theta, phi = np.deg2rad([theta, phi])
        return cls.from_axis(dirvec(theta, phi), comp=comp, isoV=isoV)


class Traction:
    @staticmethod
    def from_axis(eta, isoV: bool = False):
        return TractionCompression.from_axis(eta, False, isoV)

    @staticmethod
    def from_angles(theta: float, phi: float, isoV: bool = False, degree: bool = True):
        return TractionCompression.from_angles(theta, phi, False, isoV, degree)


class Compression(TractionCompression):
    @staticmethod
    def from_axis(eta, isoV: bool = False):
        return TractionCompression.from_axis(eta, True, isoV)

    @staticmethod
    def from_angles(theta: float, phi: float, isoV: bool = False, degree: bool = True):
        return TractionCompression.from_angles(theta, phi, True, isoV, degree)


class PureShear:
    @staticmethod
    def from_axis(eta, kappa, strict: bool = True, symmetric: bool = False):
        eta = unitvec(eta)
        kappa = unitvec(kappa)

        if strict and not check_vector_orthogonality(eta, kappa):
            raise ValueError(f"eta: {eta} and kappa: {kappa} are not orthogonal.")

        S = np.tensordot(eta, kappa, axes=0).astype(float)
        if symmetric:
            S += np.tensordot(kappa, eta, axes=0).astype(float)

        return Deformation(trim_array(S, 1.0e-8))

    @classmethod
    def from_angles(cls, theta: float, phi: float, chi: float = 0.0, symmetric: bool = False, degree: bool = True):
        if degree:
            theta, phi, chi = np.deg2rad([theta, phi, chi])
        eta, kappa, _ = spherical_basis(theta, phi, chi)
        return cls.from_axis(eta, kappa, symmetric=symmetric)


class DeformationPath(_DeformationPath):
    def __init__(self, deformation, strain_rate, tmax, tmin: float = 0.0, npts: int = 100, kpts: int = 100):
        _DeformationPath.__init__(self, deformation, strain_rate, tmin, tmax, npts, kpts)
