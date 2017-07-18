from .checker import CheckGeom
from .curves import Line
from .methods.create import create_crv_by_approx_pnts, \
    create_crv_by_interp_pnts, create_isocurve, create_line_from_occ, \
    create_nurbs_curve, create_nurbs_curve_from_occ, \
    create_nurbs_curve_from_occ2d, create_nurbs_surface, \
    create_nurbs_surface_from_occ, create_plane_by_axes, \
    create_plane_by_fit_points, create_plane_by_points, create_plane_on_curve, \
    create_planes_along_curve, create_planes_between_planes, \
    create_point_from_other, create_points_along_curve, \
    create_srf_by_approx_crvs, create_srf_by_interp_crvs
from .points import Point, Point2D
from .surfaces import Plane
from .vectors import Direction, Vector


class CreateGeom(object):
    """
    Geometry creator.
    """

    @staticmethod
    def copy_geom(geom):
        """
        Return a new copy of the geometry.

        :param geom:
        :return:
        """
        if CheckGeom.is_point_like(geom):
            geom = CheckGeom.to_point(geom)
            return geom.copy()

        if CheckGeom.is_surface(geom):
            return create_nurbs_surface_from_occ(geom)

        if CheckGeom.is_curve(geom):
            return create_nurbs_curve_from_occ(geom)

        if CheckGeom.is_line(geom):
            return create_line_from_occ(geom)

        if CheckGeom.is_curve2d(geom):
            return create_nurbs_curve_from_occ2d(geom)

        if CheckGeom.is_plane(geom):
            return Plane(geom.Pln())

        return None

    @staticmethod
    def point(xyz=(0., 0., 0.)):
        """
        Create a Point.

        :param xyz:

        :return:
        """
        return Point(xyz[0], xyz[1], xyz[2])

    @staticmethod
    def point2d(xy=(0., 0.)):
        """
        Create a Point.

        :param xy:

        :return:
        """
        return Point2D(xy[0], xy[1])

    @staticmethod
    def point_by_xyz(x=0., y=0., z=0.):
        """
        Create a Point by x-, y-, and z-location.

        :param float x: x-location of point.
        :param float y: y-location of point.
        :param float z: z-location of point.

        :return: Point at (x, y, z).
        :rtype: :class:`.Point`
        """
        return Point(x, y, z)

    @staticmethod
    def nearest_point(p, pnts):
        """
        Find the point nearest to a given point.

        :param p:
        :param pnts:

        :return:
        """
        p = CheckGeom.to_point(p)
        pnts = [CheckGeom.to_point(pi) for pi in pnts]
        dmin = p.distance(pnts[0])
        pmin = pnts[0]
        for pi in pnts[1:]:
            di = p.distance(pi)
            if di < dmin:
                dmin = di
                pmin = pi
        return pmin

    @staticmethod
    def direction_by_xyz(vx=1., vy=0., vz=0.):
        """
        Create a unit vector by components.

        :param vx:
        :param vy:
        :param vz:

        :return:
        """
        return Direction(vx, vy, vz)

    @staticmethod
    def plane_by_normal(origin=(0., 0., 0.), vnorm=(1., 0., 0.)):
        """
        Create a plane by an origin and normal vector.

        :param origin:
        :param vnorm:

        :return:
        """
        p0 = CheckGeom.to_point(origin)
        vn = CheckGeom.to_direction(vnorm)
        return Plane(p0, vn)

    @staticmethod
    def plane_by_axes(origin=(0., 0., 0.), axes='yz'):
        """
        Create a plane by origin and basic axes.

        :param origin:
        :param axes:

        :return:
        """
        origin = CheckGeom.to_point(origin)
        if not origin:
            return None

        return create_plane_by_axes(origin, axes)

    @staticmethod
    def plane_by_points(p1, p2, p3):
        """
        Create a plane using three points.

        :param p1:
        :param p2:
        :param p3:
        :return:
        """
        p1 = CheckGeom.to_point(p1)
        p2 = CheckGeom.to_point(p2)
        p3 = CheckGeom.to_point(p3)
        return create_plane_by_points(p1, p2, p3)

    @staticmethod
    def fit_plane(pnts, tol=1.0e-7):
        """
        Fit a plane to points.

        :param pnts:
        :param float tol:
        :return:
        """
        _pnts = []
        for p in pnts:
            pi = CheckGeom.to_point(p)
            if pi:
                _pnts.append(pi)
        if not _pnts:
            return None
        return create_plane_by_fit_points(_pnts, tol)

    @staticmethod
    def plane_on_curve(curve, u=None, dx=None, pnt=None, pref=None):
        """
        Create a plane on a curve.

        :param curve:
        :param u:
        :param dx
        :param pnt:
        :param pref:

        :return:
        """
        pnt = CheckGeom.to_point(pnt)
        return create_plane_on_curve(curve, u, dx, pnt, pref)

    @staticmethod
    def planes_along_curve(curve, maxd=None, npts=None, pref=None, u1=None,
                           u2=None, s1=None, s2=None):
        """
        Create planes along a curve.

        :param curve:
        :param maxd:
        :param npts:
        :param pref:
        :param u1:
        :param u2:
        :param float s1: Distance from *u1* for first plane.
        :param float s2: Distance from *u2* for last plane.

        :return:
        """
        return create_planes_along_curve(curve, maxd, npts, pref, u1, u2,
                                         s1, s2)

    @staticmethod
    def planes_between_planes(pln1, pln2, maxd=None, nplns=None,
                              s1=None, s2=None):
        """
        Create planes between two other planes.

        :param pln1:
        :param pln2:
        :param maxd:
        :param nplns:
        :param float s1:
        :param float s2:

        :return:
        """
        return create_planes_between_planes(pln1, pln2, maxd, nplns, s1, s2)

    @staticmethod
    def line_by_vector(origin, v):
        """
        Create a line by an origin an a vector.

        :param origin:
        :param v:
        :return:
        """
        origin = CheckGeom.to_point(origin)
        v = CheckGeom.to_direction(v)
        if not origin or not v:
            return None
        return Line(origin, v)

    @staticmethod
    def line_by_points(p1, p2):
        """
        Create a line between two points.

        :param p1:
        :param p2:

        :return:
        """
        p1 = CheckGeom.to_point(p1)
        p2 = CheckGeom.to_point(p2)
        if not p1 or not p2:
            return None
        v = CheckGeom.to_direction(p2.xyz - p1.xyz)
        if not v:
            return None
        return Line(p1, v)

    @staticmethod
    def nurbs_curve(cp, knots, mult, p, weights=None, is_periodic=False):
        """
        Create a NURBS curve from data.

        :param cp:
        :param knots:
        :param mult:
        :param p:
        :param weights:
        :param is_periodic:

        :return:
        """
        return create_nurbs_curve(cp, knots, mult, p, weights, is_periodic)

    @staticmethod
    def interp_points_by_curve(qp, is_periodic=False, tol=None):
        """
        Interpolate points with cubic NURBS curve.

        :param qp:
        :param is_periodic:
        :param tol:

        :return:
        """
        return create_crv_by_interp_pnts(qp, is_periodic, tol)

    @staticmethod
    def approx_points_by_curve(qp, dmin=3, dmax=8, continuity='C2',
                               tol=1.0e-3):
        """
        Fit a NURBS curve to an array of points.

        :param qp:
        :param dmin:
        :param dmax:
        :param continuity:
        :param tol:

        :return:
        """
        return create_crv_by_approx_pnts(qp, dmin, dmax, continuity, tol)

    @staticmethod
    def linear_curve(qp):
        """
        Create a linear curve through the points.

        :param qp:

        :return:
        """
        return CreateGeom.approx_points_by_curve(qp, 1, 1, 'C0', 0.)

    @staticmethod
    def approx_curves_by_surface(curves, mind=3, maxd=8, tol3d=1.0e-3,
                                 tol2d=1.0e-6, niter=5, method='chord',
                                 continuity='C2'):
        """
        Create a surface by approximating a surface through the curves.

        :param curves:
        :param int mind:
        :param int maxd:
        :param float tol3d:
        :param float tol2d:
        :param int niter:
        :param str method:
        :param str continuity:

        :return:
        """
        _crvs = [crv for crv in curves if CheckGeom.is_curve(crv)]
        if not _crvs:
            return None
        return create_srf_by_approx_crvs(_crvs, mind, maxd, tol3d, tol2d,
                                         niter, method, continuity)

    @staticmethod
    def interp_curves_by_surface(curves, q=3, method='chord'):
        """
        Create a NURBS surface by interpolating curves.

        :param curves:
        :param q:
        :param method:

        :return:
        """
        _crvs = [c for c in curves if CheckGeom.is_curve(c)]
        if not _crvs:
            return None
        return create_srf_by_interp_crvs(_crvs, q, method)

    @staticmethod
    def linear_surface(curves):
        """
        Create a surface that is linear passing through the curves.

        :param curves:

        :return:
        """
        return CreateGeom.interp_curves_by_surface(curves, 1, 'chord')

    @staticmethod
    def point_from_other(curve, dx, u0):
        """
        Create a point at a curvilinear distance from a parameter.

        :param curve:
        :param dx:
        :param u0:

        :return:
        """
        return create_point_from_other(curve, dx, u0)

    @staticmethod
    def points_along_curve(curve, maxd=None, npts=None, u1=None, u2=None,
                           s1=None, s2=None):
        """
        Generate points along a curve.

        :param curve: Curve used to generate points.
        :type curve: :class:`.BezierCurve` or :class:`.NurbsCurve`
        :param float maxd: Maximum allowed point spacing based on arc-length of
            curve. This number will be adjusted to keep equally spaced points.
        :param int npts: Number of desired points along the curve. This will be
            treated as a minimum number if both *maxd* and *npts* are provided.
        :param float u1: Starting parameter.
        :param float u2: Ending parameter.
        :param float s1: Distance from *u1* for first point.
        :param float s2: Distance from *u2* for last point.

        :return: Class containing point along curve results.
        :rtype: :class:`.CreatedPoints`
        """
        if u1 is None:
            u1 = curve.u1
        if u2 is None:
            u2 = curve.u2
        results = create_points_along_curve(curve, maxd, npts, u1, u2, s1, s2)
        return CreatedPoints(results[1], results[2])

    @staticmethod
    def nurbs_surface(cp, uknots, vknots, umult, vmult, p, q, weights=None,
                      is_u_periodic=False, is_v_periodic=False):
        """
        Create a NURBS surface.

        :param cp:
        :param uknots:
        :param vknots:
        :param umult:
        :param vmult:
        :param p:
        :param q:
        :param weights: Not yet implemented.
        :param is_u_periodic:
        :param is_v_periodic:

        :return:
        """
        return create_nurbs_surface(cp, uknots, vknots, umult, vmult, p, q,
                                    weights, is_u_periodic, is_v_periodic)

    @staticmethod
    def isocurve(surface, u=None, v=None):
        """
        Create an isocurve from the surface.

        :param surface:
        :param u:
        :param v:

        :return:
        """
        if u is None and v is None:
            return None

        if not CheckGeom.is_surface_like(surface):
            return None

        return create_isocurve(surface, u, v)

    @staticmethod
    def vector(vxyz=(0., 0., 1.)):
        """
        Create a vector.

        :param vxyz:

        :return:
        """
        return Vector(vxyz[0], vxyz[1], vxyz[2])

    @staticmethod
    def vector_by_points(p1, p2):
        """
        Create a vector between two points.

        :param p1:
        :param p2:

        :return:
        """
        p1 = CheckGeom.to_point(p1)
        p2 = CheckGeom.to_point(p2)
        return Vector(p1, p2)


class CreatedPoints(object):
    """
    Class to handle results of point creation methods.
    """

    def __init__(self, pnts, params):
        self._pnts = pnts
        self._params = params

    @property
    def npts(self):
        return len(self._pnts)

    @property
    def pnts(self):
        return self._pnts

    @property
    def pnt(self):
        try:
            return self._pnts[0]
        except IndexError:
            return None

    @property
    def params(self):
        return self._params

    @property
    def param(self):
        try:
            return self._params[0]
        except IndexError:
            return None

    @property
    def zip_results(self):
        return zip(self._pnts, self._params)
