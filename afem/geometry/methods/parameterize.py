from __future__ import division

import OCC.BSplCLib as CLib
from numpy import array, diff, float64, sqrt, sum, zeros
from numpy.linalg import norm

__all__ = []


def uniform(n, a=0., b=1.):
    """
    Generate a uniform parameters.

    :param n: Number of parameters.
    :param float a: Beginning domain if not 0.
    :param float b: Ending domain if not 1.

    :return: Uniformly spaced parameters between [a, b].
    :rtype: ndarray
    """
    u = zeros(n, dtype=float64)
    u[0] = a
    u[-1] = b
    for i in range(1, n - 1):
        u[i] = a + i * (b - a) / n
    return u


def chord_length(pnts, a=0., b=1.):
    """
    Generate parameters using chord length method.

    :param pnts: List or array of ordered points.
    :type pnts: Points or array_like
    :param float a: Beginning domain if not 0.
    :param float b: Ending domain if not 1.

    :return: Parameters between [a, b].
    :rtype: ndarray
    """
    pnts = array(pnts, dtype=float64)
    n = len(pnts)
    dtotal = sum(norm(diff(pnts, axis=0), axis=1))
    u = zeros(n, dtype=float64)
    u[0] = a
    u[-1] = b
    if dtotal <= 0.:
        return u
    for i in range(1, n - 1):
        di = norm(pnts[i] - pnts[i - 1]) / dtotal
        u[i] = u[i - 1] + di
    return u


def centripetal(pnts, a=0., b=1.):
    """
    Generate parameters using centripetal method.

    :param pnts: List or array of ordered points.
    :type pnts: Points or array_like
    :param float a: Beginning domain if not 0.
    :param float b: Ending domain if not 1.

    :return: Parameters between [a, b].
    :rtype: ndarray
    """
    pnts = array(pnts, dtype=float64)
    n = len(pnts)
    dtotal = sum(norm(diff(pnts, axis=0), axis=1))
    u = zeros(n, dtype=float64)
    u[0] = a
    u[-1] = b
    if dtotal <= 0.:
        return u
    for i in range(1, n - 1):
        di = sqrt((norm(pnts[i] - pnts[i - 1]) / dtotal))
        u[i] = u[i - 1] + di
    return u


def reparameterize_knots(u1, u2, tcol_knots):
    """
    Reparameterize the knot values between *u1* and *u2*
    """
    CLib.bsplclib_Reparametrize(u1, u2, tcol_knots)
