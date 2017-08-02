from OCC import GeomProjLib
from OCC.Geom import Geom_Line
from OCC.GeomAPI import (GeomAPI_ExtremaCurveCurve,
                         GeomAPI_ExtremaCurveSurface,
                         GeomAPI_ProjectPointOnCurve,
                         GeomAPI_ProjectPointOnSurf)
from OCC.GeomAbs import GeomAbs_BSplineCurve, GeomAbs_BezierCurve, GeomAbs_Line
from OCC.GeomAdaptor import GeomAdaptor_Curve

from afem.geometry.check import CheckGeom
from afem.geometry.methods.create import (create_line_from_occ,
                                          create_nurbs_curve_from_occ)

__all__ = ["ProjectGeom", "PointProjector", "ProjectPointToCurve",
           "ProjectPointToSurface", "CurveProjector", "ProjectCurveToPlane",
           "ProjectCurveToSurface"]


class ProjectGeom(object):
    """
    Project geometry.
    """

    @staticmethod
    def _update(point, proj, update=False):
        """
        Return projection or update point location.
        """
        if proj.success and update:
            point.set_xyz(proj.nearest_point.xyz)
            return point
        return proj

    @staticmethod
    def point_to_geom(point, geom, update=False, direction=None):
        """
        Project a point to a curve or surface.

        :param point: Point to project.
        :type point: :class:`.Point`
        :param geom: Curve or surface entity to project point to.
        :param bool update: Option to update the location of the *point*
            rather than returning a projection object. If *True*, then the
            nearest point will be used to update the point location and a
            boolean will be returned indicating a successful *True* or
            unsuccessful *False* projection.
        :param direction: Direction for point projection. Normal projection
            by default.
        :type direction: :class:`.Direction`

        :return: Projection object depending on *geom* type.
        """
        point = CheckGeom.to_point(point)

        # Project point to curve.
        if CheckGeom.is_curve_like(geom):
            proj = ProjectPointToCurve(point, geom, direction)
            return ProjectGeom._update(point, proj, update)

        # Project point to surface.
        if CheckGeom.is_surface_like(geom):
            proj = ProjectPointToSurface(point, geom, direction)
            return ProjectGeom._update(point, proj, update)

        # # Project point to plane.
        # if CheckGeom.is_plane(geom):
        #     proj = ProjectPointToPlane(point, geom)
        #     return ProjectGeom._update(point, proj, update)
        return None

    @staticmethod
    def invert(point, geom):
        """
        Return the parameters of the nearest orthogonal projection to the
        geometry.

        :param point: Point to project.
        :type point: :class:`.Point`
        :param geom: Curve or surface entity to find parameters.

        :return: Parameter(s) on the geometry of the nearest orthogonal
            projection. For a curve a single float *u* is returned,
            for a surface a tuple (u, v) is returned.
        :rtype: float or tuple

        .. note::
            *None* is returned if no points are found or the method fails. A
            tuple (*None*, *None*) is returned for a surface.
        """
        point = CheckGeom.to_point(point)

        if not CheckGeom.is_point(point):
            if CheckGeom.is_curve_like(geom):
                return None
            return None, None

        if CheckGeom.is_curve_like(geom):
            proj = ProjectPointToCurve(point, geom)
            if not proj.success:
                return None
            return proj.nearest_param

        if CheckGeom.is_surface_like(geom):
            proj = ProjectPointToSurface(point, geom)
            if not proj.success:
                return None, None
            return proj.nearest_param

        if CheckGeom.is_curve_like(geom):
            return None
        return None, None

    @staticmethod
    def curve_to_plane(curve, plane, v=None, keep_param=True):
        """
        Project a curve to a plane along a direction.

        :param curve:
        :param plane:
        :param v:
        :param keep_param:

        :return:
        """
        if CheckGeom.is_curve(curve) and CheckGeom.is_plane(plane):
            return ProjectCurveToPlane(curve, plane, v, keep_param)

    @staticmethod
    def curve_to_surface(curve, surface):
        """
        Project a curve to a surface.

        :param curve:
        :param surface:

        :return:
        """
        if CheckGeom.is_curve(curve) and CheckGeom.is_surface(surface):
            return ProjectCurveToSurface(curve, surface)


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
        :return: *True* if at least one projection was found, *False* if not.
        :rtype: bool
        """
        if self.npts > 0:
            return True
        return False

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
            return None
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
        if indx > self.npts:
            return self._results[-1][0]
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
        if indx > self.npts:
            return self._results[-1][1]
        return self._results[indx - 1][1]

    def distance(self, indx=1):
        """
        Return the projection distance by index.

        :param int indx: Index for distance.

        :return: Projection distance between original point and projection
            result.
        :rtype: float
        """
        if indx > self.npts:
            return self._results[-1][2]
        return self._results[indx - 1][2]


class ProjectPointToCurve(PointProjector):
    """
    Project a point to a curve.

    :param point_like pnt: Point to project.
    :param curve_like crv: Curve to project to.
    :param array_like direction: Direction of projection. If *None* then a
        normal projection will be performed. By providing a direction the
        tool actually performs a line-curve intersection. This is generally
        not recommended but provided by request.

    Usage:

    >>> from afem.geometry import Direction, Line, Point, ProjectPointToCurve
    >>> p0 = Point()
    >>> v = Direction(1., 0., 0.)
    >>> line = Line(p0, v)
    >>> p = Point(5., 5., 0.)
    >>> proj = ProjectPointToCurve(p, line)
    >>> assert proj.success
    >>> proj.npts
    1
    >>> proj.nearest_point
    Point(5.0, 0.0, 0.0)
    >>> proj.nearest_param
    5.0
    >>> proj.dmin
    5.0
    """

    def __init__(self, pnt, crv, direction=None):
        super(ProjectPointToCurve, self).__init__()
        pnt = CheckGeom.to_point(pnt)
        direction = CheckGeom.to_direction(direction)
        self._perform(pnt, crv, direction)

    def _perform(self, pnt, crv, direction):
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
            extrema = GeomAPI_ExtremaCurveCurve(crv.GetHandle(),
                                                geom_line.GetHandle())
            npts = extrema.NbExtrema()
            for i in range(1, npts + 1):
                ui, _ = extrema.Parameters(i)
                di = extrema.Distance(i)
                pi = crv.eval(ui)
                self._results.append([pi, ui, di])

        # Sort by distance and return.
        if self._results:
            self._results.sort(key=lambda lst: lst[2])


class ProjectPointToSurface(PointProjector):
    """
    Project a point to a surface.

    :param point_like pnt: Point to project.
    :param surface_like srf: Surface to project to.
    :param array_like direction: Direction of projection. If *None* then a
        normal projection will be performed. By providing a direction the
        tool actually performs a line-surface intersection. This is generally
        not recommended but provided by request.

    Usage:

    >>> from afem.geometry import *
    >>> p0 = Point()
    >>> n = Direction(0., 0., 1.)
    >>> pln = Plane(p0, n)
    >>> p = Point(1., 1., 1.)
    >>> proj = ProjectPointToSurface(p, pln)
    >>> assert proj.success
    >>> proj.npts
    1
    >>> proj.nearest_point
    Point(1.0, 1.0, 0.0)
    >>> proj.nearest_param
    (1.0, 1.0)
    >>> proj.dmin
    1.0
    """

    def __init__(self, pnt, srf, direction=None):
        super(ProjectPointToSurface, self).__init__()
        pnt = CheckGeom.to_point(pnt)
        direction = CheckGeom.to_direction(direction)
        self._perform(pnt, srf, direction)

    def _perform(self, point, surface, direction):
        self._results = []
        if not direction:
            # OCC projection.
            proj = GeomAPI_ProjectPointOnSurf(point, surface.handle)
            npts = proj.NbPoints()
            for i in range(1, npts + 1):
                ui, vi = proj.Parameters(i)
                di = proj.Distance(i)
                pi = surface.eval(ui, vi)
                self._results.append([pi, (ui, vi), di])
        else:
            # Use minimum distance between line and surface to project point
            # along a direction.
            geom_line = Geom_Line(point, direction)
            extrema = GeomAPI_ExtremaCurveSurface(geom_line.GetHandle(),
                                                  surface.GetHandle())
            npts = extrema.NbExtrema()
            for i in range(1, npts + 1):
                _, ui, vi = extrema.Parameters(i)
                di = extrema.Distance(i)
                pi = surface.eval(ui, vi)
                self._results.append([pi, (ui, vi), di])

        # Sort by distance and return.
        if self._results:
            self._results.sort(key=lambda lst: lst[2])


class CurveProjector(object):
    """
    Base class for curve projections.
    """

    def __init__(self):
        self._crv = None

    @property
    def success(self):
        """
        :return: *True* if a projected curve exists, *False* if not.
        :rtype: bool
        """
        return CheckGeom.is_curve(self._crv)

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

    :param curve_like crv: Curve to project.
    :param afem.geometry.entities.Plane pln: Plane to project to.
    :param array_like direction: Direction of projection. If *None* is
        provided, then the curve is projected normal to the plane.

    Usage:

    >>> from afem.geometry import *
    >>> qp = [Point(), Point(5., 5., 1.), Point(10., 5., 1.)]
    >>> c = NurbsCurveByInterp(qp).curve
    >>> pln = Plane(Point(), Direction(0., 0., 1.))
    >>> proj = ProjectCurveToPlane(c, pln, [0., 0., 1.])
    >>> assert proj.success
    >>> cnew = proj.curve
    >>> cnew.p1
    Point(0.0, 0.0, 0.0)
    >>> cnew.p2
    Point(10.0, 5.0, 0.0)
    """

    def __init__(self, crv, pln, direction=None, keep_param=True):
        super(ProjectCurveToPlane, self).__init__()
        self._perform(crv, pln, direction, keep_param)

    def _perform(self, curve, plane, v, keep_param):
        v = CheckGeom.to_direction(v)
        if not CheckGeom.is_direction(v):
            v = plane.Pln().Axis().Direction()

        # OCC projection.
        hcrv = GeomProjLib.geomprojlib_ProjectOnPlane(curve.handle,
                                                      plane.handle,
                                                      v, keep_param)
        if hcrv.IsNull():
            return None

        # TODO Support other curve types.
        crv = None
        adp_crv = GeomAdaptor_Curve(hcrv)
        if adp_crv.GetType() == GeomAbs_Line:
            lin = create_line_from_occ(adp_crv.Line())
            crv = create_line_from_occ(lin)

        if adp_crv.GetType() in [GeomAbs_BezierCurve,
                                 GeomAbs_BSplineCurve]:
            occ_crv = adp_crv.BSpline().GetObject()
            crv = create_nurbs_curve_from_occ(occ_crv)

        self._crv = crv


class ProjectCurveToSurface(CurveProjector):
    """
    Project a curve to a surface.

    :param curve_like crv: Curve to project.
    :param surface_like srf: Surface to project to.

    Usage:

    """

    # TODO Add usage docstring
    def __init__(self, crv, srf):
        super(ProjectCurveToSurface, self).__init__()
        self._perform(crv, srf)

    def _perform(self, curve, surface):
        # OCC projection. Catch error in case curve is outside the surface
        # boundaries.
        try:
            hcrv = GeomProjLib.geomprojlib_Project(curve.handle,
                                                   surface.handle)
            if hcrv.IsNull():
                return None
        except RuntimeError:
            return None

        # TODO Support other curve types.
        crv = None
        adp_crv = GeomAdaptor_Curve(hcrv)
        if adp_crv.GetType() == GeomAbs_Line:
            lin = create_line_from_occ(adp_crv.Line())
            crv = create_line_from_occ(lin)

        if adp_crv.GetType() in [GeomAbs_BezierCurve, GeomAbs_BSplineCurve]:
            occ_crv = adp_crv.BSpline().GetObject()
            crv = create_nurbs_curve_from_occ(occ_crv)

        self._crv = crv


if __name__ == "__main__":
    import doctest

    doctest.testmod()
