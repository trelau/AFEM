from OCC.Extrema import Extrema_ExtPC
from numpy import sqrt

__all__ = []


def distance_point_to_curve(point, curve):
    """
    Find the minimum distance between a point and a curve.
    """
    # OCC extrema.
    ext_pc = Extrema_ExtPC(point, curve.adaptor)
    if not ext_pc.IsDone():
        return None

    # Find the minimum result.
    n_ext = ext_pc.NbExt()
    for i in range(1, n_ext + 1):
        if ext_pc.IsMin(i):
            d = ext_pc.SquareDistance(i)
            return sqrt(d)

    return None


def curve_nearest_point(point, curves):
    """
    Find the curve nearest to the point.
    """
    ncrvs = len(curves)
    if ncrvs == 0:
        return None
    if ncrvs == 1:
        return curves[0]

    cmin = curves[0]
    dmin = distance_point_to_curve(point, cmin)
    for c in curves[1:]:
        di = distance_point_to_curve(point, c)
        if di < dmin:
            dmin = di
            cmin = c

    return cmin
