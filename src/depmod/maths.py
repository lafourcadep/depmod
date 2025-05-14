from __future__ import annotations

import numpy as np
import scipy.optimize

ArrayLike = np._typing.ArrayLike
ndarray = np.ndarray


def check_vector_orthogonality(u: ndarray, v: ndarray, tol: float = 1.0e-8):
    """Check if two vector are orthogonal.

    Two vectors x, y in R are orthogonal if x·y = 0

    Parameters
    ----------
    u, v: np.ndarray
        2 vector
    tol: float, optional, default=1.0e-8

    Returns
    -------
    bool
    """
    return np.dot(u, v) < tol


def trim_array(v: ndarray, tol: float = 1.0e-8) -> ndarray:
    w = v.copy()
    w[np.where(np.abs(w) <= tol)] = 0.0
    return w


def trim_array_inplace(v: ndarray, tol: float = 1.0e-8) -> None:
    v[np.where(np.abs(v) <= tol)] = 0.0


def unitvec(v: ndarray, tol: float = 1.0e-8) -> ndarray:
    norml2 = np.linalg.norm(v)
    if norml2 > tol:
        return v / norml2
    return v


def spherical_basis(theta, phi, chi):
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)

    cos_phi = np.cos(phi)
    sin_phi = np.sin(phi)

    cos_chi = np.cos(chi)
    sin_chi = np.sin(chi)

    u_rho = np.array([sin_theta * cos_phi, sin_theta * sin_phi, cos_theta], dtype=float)
    e_theta = np.array([cos_theta * cos_phi, cos_theta * sin_phi, -sin_theta], dtype=float)
    e_phi = np.array([-sin_phi, cos_phi, 0.0], dtype=float)

    u_theta = cos_chi * e_theta + sin_chi * e_phi
    u_phi = cos_chi * e_phi - sin_chi * e_theta

    return (u_rho, u_theta, u_phi)


def dirvec(theta: float, phi: float) -> ndarray:
    """Return the direction vector defined by the angle teta and phi.

    Parameters
    ----------
    theta, phi: float
        Values in radian for theta and phi angles.

    Returns
    -------
    np.ndarray
        The direction vector

    """
    return spherical_basis(theta, phi, chi=0.0)[0]

def fit_n_order_polynomial_with_fixed_offset(
        x: ndarray,
        y: ndarray,
        n: int = 3,
        offset: float = 0.,
        include_zero: bool = False
) -> ndarray:
    """Fit coefficients of arbitrary N-order polynomial with fixing the offset.
    
    Parameters
    ----------
    x, y: ArrayLike
        The data to be fitted

    Returns
    -------
    np.ndarray
        The normalized vector
    """
    def _cost_function(coeffs, x, y):
        return y - np.polynomial.Polynomial((offset, *coeffs))(x)

    c0 = np.polyfit(x, y, deg=n)[:n][::-1] # coefficient are returned highest order first...
    result = scipy.optimize.least_squares(_cost_function, c0, args=(x, y))
    if include_zero:
        return np.array((0., *result.x), dtype=float)
    return np.array(result.x, dtype=float)
