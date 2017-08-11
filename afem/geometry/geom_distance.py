# TODO Complet geometry distance tools
from math import sqrt

from OCC.Extrema import Extrema_ExtPC

from afem.geometry.geom_check import CheckGeom

__all__ = ["DistancePointToCurve"]


class DistancePointToCurve(object):
    """
    Calculate the extrema between a point and a curve.

    :param point_like pnt: The point.
    :param afem.geometry.geom_entities.Curve crv: The curve.
    :param float tol: The tolerance.

    :raise RuntimeError: If the extrema algorithm fails.

    Usage:

    >>> from afem.geometry import *
    >>> c = NurbsCurveByPoints([(0. ,0., 0.), (10., 0., 0.)]).curve
    >>> p = Point(5., 1., 0.)
    >>> dist = DistancePointToCurve(p, c)
    >>> dist.nsol
    3
    """

    def __init__(self, pnt, crv, tol=1.0e-10):
        tool = Extrema_ExtPC(pnt, crv.adaptor, tol)

        if not tool.IsDone():
            msg = 'Extrema between point and curve failed.'
            raise RuntimeError(msg)

        self._nsol = tool.NbExt()
        results = []
        for i in range(1, self._nsol + 1):
            di = sqrt(tool.SquareDistance(i))
            ext_pnt = tool.Point(i)
            gp_pnt = ext_pnt.Value()
            ui = ext_pnt.Parameter()
            pi = CheckGeom.to_point(gp_pnt)
            results.append((di, ui, pi))

        results.sort(key=lambda tup: tup[0])

        # TODO Remove duplicate points
        self._dmin = results[0][0]
        self._dist = [row[0] for row in results]
        self._prms = [row[1] for row in results]
        self._pnts = [row[2] for row in results]

    @property
    def nsol(self):
        """
        :return: The number of solutions.
        :rtype: int
        """
        return self._nsol

    @property
    def dmin(self):
        """
        :return: The minimum distance.
        :rtype: float
        """
        return self._dmin

    @property
    def distances(self):
        """
        :return: Sorted distances.
        :rtype: list[float]
        """
        return self._dist

    @property
    def parameters(self):
        """
        :return: Sorted parameters on curve.
        :rtype: list[float]
        """
        return self._prms

    @property
    def points(self):
        """
        :return: Sorted points on curve.
        :rtype: list[afem.geometry.geom_entities.Point]
        """
        return self._pnts


class DistancePointToSurface(object):
    pass


class DistanceCurveToCurve(object):
    pass


class DistanceCurveToSurface(object):
    pass
