from .checker import CheckGeom
from .methods.project import project_curve_to_surface, \
    project_point_to_curve, project_point_to_surface, project_curve_to_plane


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
    def point_to_geom(point, geom, update=False):
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

        :return: Projection object depending on *geom* type.
        """
        point = CheckGeom.to_point(point)

        # Project point to curve.
        if CheckGeom.is_curve_like(geom):
            proj = ProjectPointToCurve(point, geom)
            return ProjectGeom._update(point, proj, update)

        # Project point to surface.
        if CheckGeom.is_surface_like(geom):
            proj = ProjectPointToSurface(point, geom)
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
    Base for handling point projections to curves and surfaces.
    """

    def __init__(self):
        self._npts = 0
        self._results = []

    @property
    def npts(self):
        return len(self._results)

    @property
    def success(self):
        if self.npts > 0:
            return True
        return False

    @property
    def points(self):
        if self.npts <= 0:
            return None
        return [results[0] for results in self._results]

    @property
    def parameters(self):
        if self.npts <= 0:
            return None
        return [results[1] for results in self._results]

    @property
    def nearest_point(self):
        return self._results[0][0]

    @property
    def nearest_param(self):
        return self._results[0][1]

    @property
    def dmin(self):
        return self._results[0][2]

    def point(self, indx=1):
        """
        Return the point result by index.

        :param int indx: Index for point selection.

        :return: Projection point at index.
        :rtype: :class:`.Point`
        """
        if indx > self.npts:
            return self._results[-1][0]
        return self._results[indx - 1][0]

    def parameter(self, indx=1):
        """
        Return the parameter result by index.

        :param int indx: Index for parameter selection.

        :return: Projection parameter at index. For a curve projection a
            single float *u* will be returned. For a surface projection result
            a tuple containing the *u* and *v* parameters will be returned
            (u,v).
        :rtype: float or tuple
        """
        if indx > self.npts:
            return self._results[-1][1]
        return self._results[indx - 1][1]

    def distance(self, indx=1):
        """
        Return the projection distance by index.

        :param int indx: Index for distance selection.

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
    """

    def __init__(self, point, curve):
        super(ProjectPointToCurve, self).__init__()
        if CheckGeom.is_point(point) and CheckGeom.is_curve_like(curve):
            self._perform(point, curve)

    def _perform(self, point, curve):
        """
        Perform the projection.
        """
        results = project_point_to_curve(point, curve)
        self._results = results


class ProjectPointToSurface(PointProjector):
    """
    Project a point to a surface.
    """

    def __init__(self, point, surface):
        super(ProjectPointToSurface, self).__init__()
        point = CheckGeom.to_point(point)
        if CheckGeom.is_point(point) and CheckGeom.is_surface_like(surface):
            self._perform(point, surface)

    def _perform(self, point, surface):
        """
        Perform the point to surface projection.
        """
        results = project_point_to_surface(point, surface)
        self._results = results


class CurveProjector(object):
    """
    Base class for curve projections.
    """

    def __init__(self):
        self._crv = None

    @property
    def success(self):
        return CheckGeom.is_curve(self._crv)

    @property
    def result(self):
        return self._crv


class ProjectCurveToPlane(CurveProjector):
    """
    Project a curve to a plane along a direction.
    """

    def __init__(self, curve, plane, v=None, keep_param=True):
        super(ProjectCurveToPlane, self).__init__()
        if CheckGeom.is_curve(curve) and CheckGeom.is_plane(plane):
            self._perform(curve, plane, v, keep_param)

    def _perform(self, curve, plane, v, keep_param):
        v = CheckGeom.to_direction(v)
        if not CheckGeom.is_direction(v):
            v = plane.Pln().Axis().Direction()
        self._crv = project_curve_to_plane(curve, plane, v, keep_param)


class ProjectCurveToSurface(CurveProjector):
    """
    Project a curve to a surface.
    """

    def __init__(self, curve, surface):
        super(ProjectCurveToSurface, self).__init__()
        self._crv = None
        if CheckGeom.is_curve(curve) and CheckGeom.is_surface(surface):
            self._perform(curve, surface)

    def _perform(self, curve, surface):
        self._crv = project_curve_to_surface(curve, surface)
