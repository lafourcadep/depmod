from __future__ import annotations

import numpy as np

from depmod._lib import _CustomExpressionStrainRate, _Deformation, _DeformationPath, _StrainRate
from depmod.maths import (check_vector_orthogonality, dirvec, spherical_basis,
                          trim_array, unitvec)

# Some typing definition
double = float
size_t = int
ArrayLike = np._typing.ArrayLike
NDArray = np._typing.NDArray


class Deformation(_Deformation):
    def __init__(self, S: NDArray, cfac: double = 1.0) -> None:
        _Deformation.__init__(self, S, cfac)


class TractionCompression(_Deformation):
    def __init__(self, S: NDArray, comp: bool = False) -> None:
        cfac = 1.0 if not comp else -1.0
        _Deformation.__init__(self, S, cfac)

    @classmethod
    def from_axis(cls, eta: ArrayLike, comp: bool = False, isoV: bool = False) -> TractionCompression:
        eta = unitvec(eta)
        S = np.tensordot(eta, eta, axes=0).astype(float)

        if isoV:
            theta = np.arccos(eta[2])  # arccos(z)
            phi = np.arctan2(eta[1], eta[0])  # arctan2(y, x)
            m, n, p = spherical_basis(theta, phi, chi=0.0)
            if not np.allclose(m, eta):
                raise ValueError("Something went wrong in computing the spherical basis")
            S -= 0.5 * (np.tensordot(n, n, axes=0) + np.tensordot(p, p, axes=0))

        return cls(trim_array(S, 1.0e-8), comp=comp)

    @classmethod
    def from_angles(
        cls, theta: double, phi: double, comp: bool = False, isoV: bool = False, degree: bool = True
    ) -> TractionCompression:
        if degree:
            theta, phi = np.deg2rad([theta, phi])
        return cls.from_axis(dirvec(theta, phi), comp=comp, isoV=isoV)


class Traction:
    @staticmethod
    def from_axis(eta, isoV: bool = False) -> TractionCompression:
        return TractionCompression.from_axis(eta, False, isoV)

    @staticmethod
    def from_angles(theta: double, phi: double, isoV: bool = False, degree: bool = True) -> TractionCompression:
        return TractionCompression.from_angles(theta, phi, False, isoV, degree)


class Compression(TractionCompression):
    @staticmethod
    def from_axis(eta: ArrayLike, isoV: bool = False) -> TractionCompression:
        return TractionCompression.from_axis(eta, True, isoV)

    @staticmethod
    def from_angles(theta: double, phi: double, isoV: bool = False, degree: bool = True) -> TractionCompression:
        return TractionCompression.from_angles(theta, phi, True, isoV, degree)


class PureShear:
    @staticmethod
    def from_axis(eta: ArrayLike, kappa: ArrayLike, strict: bool = True, symmetric: bool = False) -> Deformation:
        eta = unitvec(eta)
        kappa = unitvec(kappa)

        if strict and not check_vector_orthogonality(eta, kappa):
            raise ValueError(f"eta: {eta} and kappa: {kappa} are not orthogonal.")

        S = np.tensordot(eta, kappa, axes=0).astype(float)
        if symmetric:
            S += np.tensordot(kappa, eta, axes=0).astype(float)

        return Deformation(trim_array(S, 1.0e-8))

    @classmethod
    def from_angles(
        cls, theta: double, phi: double, chi: double = 0.0, symmetric: bool = False, degree: bool = True
    ) -> Deformation:
        if degree:
            theta, phi, chi = np.deg2rad([theta, phi, chi])
        eta, kappa, _ = spherical_basis(theta, phi, chi)
        return cls.from_axis(eta, kappa, symmetric=symmetric)


class DeformationPath(_DeformationPath):
    def __init__(
        self,
        deformation: Deformation,
        strain_rate: double | str,
        tmax: double,
        tmin: double = 0.0,
        npts: size_t = 100,
        kpts: size_t = 100,
        kw: dict | None = None,
    ) -> None:
        strain_rate_obj = None
        kw = kw or {}

        if isinstance(strain_rate, (double, int)):
            strain_rate_obj = _StrainRate(strain_rate)
        elif isinstance(strain_rate, str):
            expr_str = strain_rate
            for k, v in kw.items():
                expr_str = expr_str.replace(k, str(v))

            strain_rate_obj = _CustomExpressionStrainRate(expr_str)

        if not strain_rate_obj:
            raise RuntimeError("Invalid strain_rate...")

        _DeformationPath.__init__(self, deformation, strain_rate_obj, tmin, tmax, npts, kpts)
