from __future__ import division

import OCC.BSplCLib as CLib
from OCC.TColStd import TColStd_Array1OfReal
from numpy import hstack, zeros

from ...utils.tcol import to_np_from_tcolstd_array1_real

__all__ = []


def local_to_global_param(a, b, *args):
    if args[0] is None:
        return None
    local_u = args
    is_list = True
    if len(local_u) == 1:
        is_list = False
    global_u = []
    for ui in local_u:
        if ui < 0.:
            ui = 0.
        if ui > 1.:
            ui = 1.
        global_u.append(a + ui * (b - a))
    if is_list:
        return global_u
    else:
        return global_u[0]


def global_to_local_param(a, b, *args):
    if args[0] is None:
        return None
    global_u = args
    is_list = True
    if len(global_u) == 1:
        is_list = False
    local_u = []
    for ui in global_u:
        u = (ui - a) / (b - a)
        if u < 0.:
            u = 0.
        if u > 1.:
            u = 1.
        local_u.append(u)
    if is_list:
        return local_u
    else:
        return local_u[0]


def homogenize_array1d(cp, w):
    _w = w.reshape(-1, 1)
    return hstack((cp * _w, _w))


def homogenize_array2d(cp, w):
    n, m, _ = cp.shape
    cpw = zeros((n, m, 4), dtype=float)
    for i in range(n):
        for j in range(m):
            cpw[i, j, :3] = cp[i, j] * w[i, j]
            cpw[i, j, 3] = w[i, j]
    return cpw


def dehomogenize_array1d(cpw):
    w = cpw[:, -1]
    cp = cpw[:, :-1] / cpw[:, -1].reshape(-1, 1)
    return cp, w


def dehomogenize_array2d(cpw):
    cp = cpw[:, :, :3]
    w = cpw[:, :, -1]
    n, m = cp.shape[0:2]
    for i in range(n):
        for j in range(m):
            cp[i, j] /= w[i, j]
    return cp, w


def build_knot_vector_from_occ(tcol_knots, tcol_mult, p, is_periodic):
    """
    Build knot sequence from OCC data.
    """
    n = CLib.bsplclib_KnotSequenceLength(tcol_mult, p, is_periodic)
    tcol_knots_seq = TColStd_Array1OfReal(1, n)
    CLib.bsplclib_KnotSequence(tcol_knots, tcol_mult, p, is_periodic,
                               tcol_knots_seq)
    return to_np_from_tcolstd_array1_real(tcol_knots_seq)
