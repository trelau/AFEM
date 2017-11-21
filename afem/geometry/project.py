from OCCT.Geom import Geom_Line
from OCCT.GeomAPI import (GeomAPI_ExtremaCurveCurve,
                          GeomAPI_ExtremaCurveSurface,
                          GeomAPI_ProjectPointOnCurve,
                          GeomAPI_ProjectPointOnSurf)
from OCCT.GeomProjLib import GeomProjLib

from afem.geometry.check import CheckGeom
from afem.geometry.entities import Curve

__all__ = ["PointProjector", "ProjectPointToCurve",
           "ProjectPointToSurface", "CurveProjector", "ProjectCurveToPlane",
           "ProjectCurveToSurface"]


class PointProjector(object):
    """
    Base class for point projections.
    """

    def __init__(self):
        self._npts = 0
        self._results = []

    @property
    def npts(self):
        """
        :return: Number of projection results.
        :rtype: int
        """
        return len(self._results)

    @property
    def success(self):
        """
        :return: *True* if successful, *False* if not.
        :rtype: bool
        """
        return self.npts > 0

    @property
    def points(self):
        """
        :return: Projected points.
        :rtype: list[afem.geometry.entities.Point]
        """
        if self.npts <= 0:
            return []
        return [results[0] for results in self._results]

    @property
    def parameters(self):
        """
        :return: Parameters of projected points.
        :rtype: list[float]
        """
        if self.npts <= 0:
            return []
        return [results[1] for results in self._results]

    @property
    def nearest_point(self):
        """
        :return: Nearest projection result to original point.
        :rtype: afem.geometry.entities.Point
        """
        return self._results[0][0]

    @property
    def nearest_param(self):
        """
        :return: Parameter of nearest point.
        :rtype: float
        """
        return self._results[0][1]

    @property
    def dmin(self):
        """
        :return: Minimum distance of all projection results.
        :rtype: float
        """
        return self._results[0][2]

    def point(self, indx=1):
        """
        Return the point result by index.

        :param int indx: Index for point.

        :return: Projected point.
        :rtype: afem.geometry.entities.Point
        """
        return self._results[indx - 1][0]

    def parameter(self, indx=1):
        """
        Return the parameter result by index.

        :param int indx: Index for parameter.

        :return: Parameter of point. For a curve projection a single float *u*
            will be returned. For a surface projection a tuple containing the
            *u* and *v* parameters will be returned (u, v).
        :rtype: float or tuple(float, float)
        """
        return self._results[indx - 1][1]

    def distance(self, indx=1):
        """
        Return the projection distance by index.

        :param int indx: Index for distance.

        :return: Projection distance between original point and projection
            result.
        :rtype: float
        """
        return self._results[indx - 1][2]


class ProjectPointToCurve(PointProjector):
    """
    Project a point to a curve.

    :param point_like pnt: Point to project.
    :param afem.geometry.entities.Curve crv: Curve to project to.
    :param array_like direction: Direction of projection. If *None* then a
        normal projection will be performed. By providing a direction the
        tool actually performs a line-curve intersection. This is generally
        not recommended but provided by request.
    :param bool update: Option to update the point's location to match the
        nearest point.

    For more information see GeomAPI_ProjectPointOnCurve_.

    .. _GeomAPI_ProjectPointOnCurve: https://www.opencascade.com/doc/occt-7.2.0/refman/html/class_geom_a_p_i___project_point_on_curve.html

    Usage:

    >>> from afem.geometry import *
    >>> p0 = Point()
    >>> v = Direction(1., 0., 0.)
    >>> line = LineByVector(p0, v).line
    >>> p = Point(5., 5., 0.)
    >>> proj = ProjectPointToCurve(p, line)
    >>> assert proj.success
    >>> proj.npts
    1
    >>> proj.nearest_point
    Point(5.000, 0.000, 0.000)
    >>> proj.nearest_param
    5.0
    >>> proj.dmin
    5.0
    """

    def __init__(self, pnt, crv, direction=None, update=False):
        super(ProjectPointToCurve, self).__init__()

        # Perform
        pnt = CheckGeom.to_point(pnt)
        direction = CheckGeom.to_direction(direction)
        self._results = []

        if not direction:
            # OCC projection.
            proj = GeomAPI_ProjectPointOnCurve(pnt, crv.handle)
            npts = proj.NbPoints()
            for i in range(1, npts + 1):
                ui = proj.Parameter(i)
                di = proj.Distance(i)
                pi = crv.eval(ui)
                self._results.append([pi, ui, di])
        else:
            # Use minimum distance between line and curve to project point
            # along a direction.
            geom_line = Geom_Line(pnt, direction)
            extrema = GeomAPI_ExtremaCurveCurve(crv.handle,
                                                geom_line)
            npts = extrema.NbExtrema()
            for i in range(1, npts + 1):
                ui, _ = extrema.Parameters(i, 0., 0.)
                di = extrema.Distance(i)
                pi = crv.eval(ui)
                self._results.append([pi, ui, di])

        # Sort by distance and return.
        if self._results:
            self._results.sort(key=lambda lst: lst[2])

        if update:
            pnt.set_xyz(self.nearest_point)


class ProjectPointToSurface(PointProjector):
    """
    Project a point to a surface.

    :param point_like pnt: Point to project.
    :param afem.geometry.entities.Surface srf: Surface to project to.
    :param array_like direction: Direction of projection. If *None* then a
        normal projection will be performed. By providing a direction the
        tool actually performs a line-surface intersection. This is generally
        not recommended but provided by request.
    :param bool update: Option to update the point's location to match the
        nearest point.

    For more information see GeomAPI_ProjectPointOnSurf_.

    .. _GeomAPI_ProjectPointOnSurf: https://www.opencascade.com/doc/occt-7.2.0/refman/html/class_geom_a_p_i___project_point_on_surf.html

    Usage:

    >>> from afem.geometry import *
    >>> p0 = Point()
    >>> n = Direction(0., 0., 1.)
    >>> pln = PlaneByNormal(p0, n).plane
    >>> p = Point(1., 1., 1.)
    >>> proj = ProjectPointToSurface(p, pln)
    >>> assert proj.success
    >>> proj.npts
    1
    >>> proj.nearest_point
    Point(1.000, 1.000, 0.000)
    >>> proj.nearest_param
    (1.0, 1.0)
    >>> proj.dmin
    1.0
    """

    def __init__(self, pnt, srf, direction=None, update=False):
        super(ProjectPointToSurface, self).__init__()

        # Perform
        pnt = CheckGeom.to_point(pnt)
        direction = CheckGeom.to_direction(direction)
        self._results = []

        if not direction:
            # OCC projection.
            proj = GeomAPI_ProjectPointOnSurf(pnt, srf.handle)
            npts = proj.NbPoints()
            for i in range(1, npts + 1):
                ui, vi = proj.Parameters(i, 0., 0.)
                di = proj.Distance(i)
                pi = srf.eval(ui, vi)
                self._results.append([pi, (ui, vi), di])
        else:
            # Use minimum distance between line and surface to project point
            # along a direction.
            geom_line = Geom_Line(pnt, direction)
            extrema = GeomAPI_ExtremaCurveSurface(geom_line,
                                                  srf.handle)
            npts = extrema.NbExtrema()
            for i in range(1, npts + 1):
                _, ui, vi = extrema.Parameters(i, 0., 0., 0.)
                di = extrema.Distance(i)
                pi = srf.eval(ui, vi)
                self._results.append([pi, (ui, vi), di])

        # Sort by distance and return.
        if self._results:
            self._results.sort(key=lambda lst: lst[2])

        if update:
            pnt.set_xyz(self.nearest_point)


class CurveProjector(object):
    """
    Base class for curve projections.
    """

    def __init__(self):
        self._crv = None

    @property
    def success(self):
        """
        :return: *True* if successful, *False* if not.
        :rtype: bool
        """
        return self._crv is not None

    @property
    def curve(self):
        """
        :return: The projected curve.
        :rtype: afem.geometry.entities.Curve
        """
        return self._crv


class ProjectCurveToPlane(CurveProjector):
    """
    Project a curve to a plane along a direction.

    :param afem.geometry.entities.Curve crv: Curve to project.
    :param afem.geometry.entities.Plane pln: Plane to project to.
    :param array_like direction: Direction of projection. If *None* is
        provided, then the curve is projected normal to the plane.

    :raise RuntimeError: If the OCC method fails to project the curve to the
        plane.

    For more information see GeomProjLib_ProjectOnPlane_.

    .. _GeomProjLib_ProjectOnPlane: https://www.opencascade.com/doc/occt-7.2.0/refman/html/class_geom_proj_lib.html

    Usage:

    >>> from afem.geometry import *
    >>> qp = [Point(), Point(5., 5., 1.), Point(10., 5., 1.)]
    >>> c = NurbsCurveByInterp(qp).curve
    >>> pln = PlaneByNormal(Point(), Direction(0., 0., 1.)).plane
    >>> proj = ProjectCurveToPlane(c, pln, [0., 0., 1.])
    >>> assert proj.success
    >>> cnew = proj.curve
    >>> cnew.p1
    Point(0.000, 0.000, 0.000)
    >>> cnew.p2
    Point(10.000, 5.000, 0.000)
    """

    def __init__(self, crv, pln, direction=None, keep_param=True):
        super(ProjectCurveToPlane, self).__init__()

        direction = CheckGeom.to_direction(direction)
        if not CheckGeom.is_direction(direction):
            direction = pln.handle.Pln().Axis().Direction()

        # OCC projection
        hcrv = GeomProjLib.ProjectOnPlane_(crv.handle, pln.handle, direction,
                                           keep_param)

        self._crv = Curve(hcrv)


class ProjectCurveToSurface(CurveProjector):
    """
    Project a curve to a surface. Only normal projections are supported.

    :param afem.geometry.entities.Curve crv: Curve to project.
    :param afem.geometry.entities.Surface srf: Surface to project to.

    :raise RuntimeError: If the OCC method fails to project the curve to the
        plane.

    For more information see GeomProjLib_Project_.

    .. _GeomProjLib_Project: https://www.opencascade.com/doc/occt-7.2.0/refman/html/class_geom_proj_lib.html

    Usage:

    >>> from afem.geometry import *
    >>> c = NurbsCurveByPoints([(0., 5., 6.), (10., 5., 6.)]).curve
    >>> c1 = NurbsCurveByPoints([(0., 0., 0.), (10., 0., 0.)]).curve
    >>> c2 = NurbsCurveByPoints([(0., 5., 5.), (10., 5., 5.)]).curve
    >>> c3 = NurbsCurveByPoints([(0., 10., 0.), (10., 10., 0.)]).curve
    >>> s = NurbsSurfaceByApprox([c1, c2, c3]).surface
    >>> proj = ProjectCurveToSurface(c, s)
    >>> assert proj.success
    >>> cproj = proj.curve
    >>> cproj.eval(0.5)
    Point(5.000, 5.000, 5.000)
    """

    def __init__(self, crv, srf):
        super(ProjectCurveToSurface, self).__init__()

        # OCC projection
        hcrv = GeomProjLib.Project_(crv.handle, srf.handle)
        self._crv = Curve(hcrv)


if __name__ == "__main__":
    import doctest

    doctest.testmod()
