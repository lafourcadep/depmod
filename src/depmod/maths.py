import numpy as np

import scipy.optimize

from depmod.typing import ArrayLike, ndarray


def trim_array_inplace(v: ndarray, tol: float = 1.0e-8) -> None:
    v[np.where(np.abs(v) <= tol)] = 0.0


def trim_array(v: ndarray, tol: float = 1.0e-8) -> ndarray:
    w = v.copy()
    w[np.where(np.abs(w) <= tol)] = 0.0
    return w


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
    return np.array(
        [
            np.sin(theta) * np.cos(phi),
            np.sin(theta) * np.sin(phi),
            np.cos(theta)
        ],
        dtype=float
    )


def unitvec(v: ndarray) -> ndarray:
    """Normalize the input vector.

    Parameters
    ----------
    v: np.ndarray
        The input vector
    tol: float, optional, default=1.e-8


    Returns
    -------
    np.ndarray
        The normalized vector

    """
    tol = 1.0e-8
    norml2 = np.linalg.norm(v)
    if norml2 > tol:
        return v / norml2
    return v


def  check_vector_orthogonality(u: ndarray, v: ndarray, tol: float = 1.0e-8):
    """Check if two vector are orthogonal.

    Two vectors x, y in R are orthogonal if x·y = 0
    """
    return np.dot(u, v) < tol


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


def rotation_matrix_from_vectors(v1: ndarray, v2: ndarray) -> ndarray:
    """Find the rotation matrix Q, that aligns v1 to v2.

    Parameters
    ----------
    v1, v2: np.ndarray
        A 3D "source" and "destination" vectors

    Returns
    -------
    np.ndarray
        A (3x3) rotation matrix which aligns v1 to v2

    """
    a = unitvec(v1)
    b = unitvec(v2)

    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    K = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]], dtype=np.float64)
    Q = np.eye(3) + K + np.dot(K, K) * ((1 - c) / (s**2))

    return Q


def build_rotation_matrix_around_axis(axis: ndarray, angle: float) -> ndarray:
    """Build a rotation matrix Q, that rotate around a given axis.

    Parameters
    ----------
    # TODO

    Returns
    -------
    # TODO

    """
    angle = np.deg2rad(angle)  # convert angle to rad
    a = unitvec(axis)
    K = np.array([[0.0, a[2], -a[1]], [-a[2], 0.0, a[0]], [a[1], -a[0], 0.0]], dtype=float)
    Q = np.eye(3) - K * np.sin(angle) + np.dot(K, K) * (1.0 - np.cos(angle))

    return Q


def align_to_lammps_convention(
    H: ndarray, zaxis: ndarray, tol: float = 1.0e-5
) -> tuple[ndarray, ndarray]:
    """Align a lattice to respect the LAMMPS convention.

    For more details see: https://docs.lammps.org/Howto_triclinic.html

    Parameters
    ----------
    # TODO

    Returns
    -------
    # TODO

    """
    # 1st step : bring back a and b in (x,y) plane by aligning a^b with z-axis
    H1 = np.zeros((3, 3), dtype=float)
    a = H[:, 0]
    b = H[:, 1]
    axb = unitvec(np.cross(a, b))
    if np.linalg.norm(axb - zaxis) < tol:
        H1[:, :] = H[:, :]
    else:
        Q1 = rotation_matrix_from_vectors(axb, zaxis)
        H1[:, :] = np.dot(Q1, H)

    # 2nd step: bring back a along x axis through a rotation around z axis
    xdir = 180.0 * np.arctan(H1[1, 0] / H1[0, 0]) / np.pi
    Q2 = build_rotation_matrix_around_axis(zaxis, -1.0 * xdir)
    # H2 = np.dot(Q2, H1)

    # H2 = Q2 * Q1 * H
    Q = np.dot(Q2, Q1)
    H2 = np.dot(Q, H)

    return H2, Q
