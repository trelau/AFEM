from math import ceil, radians

from OCC.Approx import (Approx_Centripetal, Approx_ChordLength,
                        Approx_IsoParametric)
from OCC.GCPnts import GCPnts_AbscissaPoint, GCPnts_UniformAbscissa
from OCC.Geom import Geom_Curve
from OCC.GeomAPI import GeomAPI_Interpolate, GeomAPI_PointsToBSpline
from OCC.GeomAbs import (GeomAbs_BSplineCurve, GeomAbs_C0, GeomAbs_C1,
                         GeomAbs_C2, GeomAbs_C3, GeomAbs_G1, GeomAbs_G2,
                         GeomAbs_Line)
from OCC.GeomAdaptor import GeomAdaptor_Curve
from OCC.TColStd import TColStd_Array1OfInteger, TColStd_Array1OfReal
from OCC.TColgp import TColgp_Array1OfPnt
from OCC.gp import gp_Quaternion, gp_Trsf

from afem.geometry.check import CheckGeom
from afem.geometry.entities import *
from afem.geometry.methods.create import create_crv_by_approx_pnts, \
    create_crv_by_interp_pnts, create_isocurve, create_line_from_occ, \
    create_nurbs_curve, create_nurbs_curve_from_occ, \
    create_nurbs_curve_from_occ2d, create_nurbs_surface, \
    create_nurbs_surface_from_occ, create_plane_by_axes, \
    create_plane_by_fit_points, create_plane_by_points, create_plane_on_curve, \
    create_planes_along_curve, create_planes_between_planes, \
    create_point_from_other, create_points_along_curve, \
    create_srf_by_approx_crvs, create_srf_by_interp_crvs
from afem.occ.utils import (to_tcolgp_array1_pnt, to_tcolgp_harray1_pnt,
                            to_tcolstd_array1_integer, to_tcolstd_array1_real)

__all__ = ["CreateGeom", "CreatedPoints", "PointByXYZ", "PointByArray",
           "PointsAlongCurveByNumber",
           "PointFromParameter", "PointsAlongCurveByDistance",
           "DirectionByXYZ", "DirectionByPoints",
           "DirectionByArray", "VectorByXYZ", "VectorByArray",
           "VectorByPoints", "LineByVector", "LineByPoints",
           "NurbsCurveByData",
           "NurbsCurveByInterp", "NurbsCurveByApprox", "NurbsCurveByPoints",
           "CurveByUIso",
           "CurveByVIso", "PlaneByNormal", "PlaneByAxes", "PlaneByPoints",
           "PlaneByFit", "PlanesAlongCurve", "PlanesBetweenPlanes",
           "SurfaceByData", "SurfaceByInterp", "SurfaceByFit"]

_occ_continuity = {'C0': GeomAbs_C0,
                   'G1': GeomAbs_G1,
                   'C1': GeomAbs_C1,
                   'G2': GeomAbs_G2,
                   'C2': GeomAbs_C2,
                   'C3': GeomAbs_C3}

_occ_parm_type = {'u': Approx_IsoParametric,
                  'uniform': Approx_IsoParametric,
                  'centripetal': Approx_Centripetal,
                  'c': Approx_ChordLength,
                  'chord': Approx_ChordLength,
                  'chord length': Approx_ChordLength}


def _create_nurbs_curve_from_occ(crv):
    """
    Create a NURBS curve from an OCC curve.
    """
    # Gather OCC data.
    tcol_poles = TColgp_Array1OfPnt(1, crv.NbPoles())
    crv.Poles(tcol_poles)
    tcol_weights = TColStd_Array1OfReal(1, crv.NbPoles())
    crv.Weights(tcol_weights)
    tcol_knots = TColStd_Array1OfReal(1, crv.NbKnots())
    crv.Knots(tcol_knots)
    tcol_mult = TColStd_Array1OfInteger(1, crv.NbKnots())
    crv.Multiplicities(tcol_mult)
    p = crv.Degree()
    is_periodic = crv.IsPeriodic()
    return NurbsCurve(tcol_poles, tcol_weights, tcol_knots, tcol_mult, p,
                      is_periodic)


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
    def plane_by_orientation(origin=(0., 0., 0.), axes='xz',
                             alpha=0., beta=0., gamma=0.):
        """
        Create a plane by orientation angles.

        :param origin:
        :param axes:
        :param alpha:
        :param beta:
        :param gamma:

        :return:
        """
        pln = CreateGeom.plane_by_axes((0., 0., 0.), axes)
        if not pln:
            return None

        # Build a quaternion for rotation angles
        r = gp_Quaternion()
        r.SetEulerAngles(2, radians(alpha), radians(beta), radians(gamma))

        # Build transformation matrix and rotate about global origin and
        # translate to desired plane origin.
        tf = gp_Trsf()
        v = CreateGeom.vector(origin)
        tf.SetTransformation(r, v)
        pln.Transform(tf)

        return pln

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


class GeomBuilder(object):
    """
    Base class for geometry builder.
    """

    def __init__(self):
        self._success = False
        self._performed = False
        self._results = {}

    @property
    def success(self):
        """
        Status of builder.

        :return: *True* if successful, *False* if not.
        :rtype: bool
        """
        if not self._performed:
            self.build()
        return self._success

    def _set_results(self, key, value):
        self._results[key] = value

    def _get_results(self, key):
        if not self._performed:
            self.build()
        try:
            return self._results[key]
        except KeyError:
            return None

    def build(self):
        """
        Build geometry. Should be overridden in subclasses.

        :return: None.
        :rtype: None
        """
        pass


class PointByXYZ(GeomBuilder):
    """
    Create a point by x, y, and z location.

    :param float x: x-location.
    :param float y: y-location.
    :param float z: z-location.

    Usage:

    >>> from afem.geometry import PointByXYZ
    >>> PointByXYZ(1., 2., 3.).point
    Point(1.0, 2.0, 3.0)
    """

    def __init__(self, x=0., y=0., z=0.):
        super(PointByXYZ, self).__init__()
        self._x = float(x)
        self._y = float(y)
        self._z = float(z)

    def build(self):
        self._performed = True

        p = Point(self._x, self._y, self._z)
        self._set_results('p', p)
        self._success = isinstance(p, Point)

    @property
    def point(self):
        """
        :return: The created point.
        :rtype: afem.geometry.entities.Point
        """
        return self._get_results('p')


class PointByArray(PointByXYZ):
    """
    Create a point from an array-like object.

    :param array_like xyz: Array-like object describing point location.

    Usage:

    >>> from afem.geometry import PointByArray
    >>> PointByArray([1., 2., 3.]).point
    Point(1.0, 2.0, 3.0)
    """

    def __init__(self, xyz=(0., 0., 0.)):
        assert len(xyz) == 3, "Invalid array size in PointByArray"
        x, y, z = xyz
        super(PointByArray, self).__init__(x, y, z)


class PointFromParameter(GeomBuilder):
    """
    Create a point along a curve at a specified distance from a parameter.

    :param curve_like c: The curve.
    :param float u0: The initial parameter.
    :param float ds: The distance along the curve from the given parameter.
    :param float tol: Tolerance.

    Usage:

    >>> from afem.geometry import LineByPoints, Point, PointFromParameter
    >>> p1 = Point()
    >>> p2 = Point(10., 0., 0.)
    >>> line = LineByPoints(p1, p2).line
    >>> builder = PointFromParameter(line, 0., 1.)
    >>> assert builder.success
    >>> builder.point
    Point(1.0, 0.0, 0.0)
    >>> builder.parameter
    1.0
    """

    def __init__(self, c, u0, ds, tol=1.0e-7):
        super(PointFromParameter, self).__init__()
        self._c = c
        self._ds = ds
        self._u0 = u0
        self._tol = tol

    def build(self):
        self._performed = True
        if not CheckGeom.is_curve_like(self._c):
            return None

        adp_curve = GeomAdaptor_Curve(self._c.handle)

        ap = GCPnts_AbscissaPoint(self._tol, adp_curve, self._ds, self._u0)
        if not ap.IsDone():
            self._success = False
            return None

        u = ap.Parameter()
        self._set_results('u', u)
        p = self._c.eval(u)
        self._set_results('p', p)
        self._success = isinstance(p, Point)

    @property
    def point(self):
        """
        :return: The created point.
        :rtype: afem.geometry.entities.Point
        """
        return self._get_results('p')

    @property
    def parameter(self):
        """
        :return: The parameter on the curve.
        :rtype: float
        """
        return self._get_results('u')


class PointsAlongCurveByNumber(GeomBuilder):
    """
    Create a specified number of points along a curve. The points will be
    equidistant.

    :param curve_like c: The curve.
    :param int n: Number of points to create (*n* > 0).
    :param float u1: The parameter of the first point (default=c.u1).
    :param float u2: The parameter of the last point (default=c.u2).
    :param float d1: An offset distance for the first point.This is typically
        a positive number indicating a distance from *u1* towards *u2*.
    :param float d2: An offset distance for the last point. This is typically
        a negative number indicating a distance from *u2* towards *u1*.
    :param float tol: Tolerance.

    Usage:

    >>> from afem.geometry import LineByPoints, Point, PointsAlongCurveByNumber
    >>> p1 = Point()
    >>> p2 = Point(10., 0., 0.)
    >>> line = LineByPoints(p1, p2).line
    >>> builder = PointsAlongCurveByNumber(line, 3, 0., 10.)
    >>> assert builder.success
    >>> builder.npts
    3
    >>> builder.points
    [Point(0.0, 0.0, 0.0), Point(5.0, 0.0, 0.0), Point(10.0, 0.0, 0.0)]
    >>> builder.parameters
    [0.0, 5.0, 10.0]
    """

    def __init__(self, c, n, u1=None, u2=None, d1=None, d2=None, tol=1.0e-7):
        super(PointsAlongCurveByNumber, self).__init__()
        self._c = c
        self._n = int(n)
        self._u1 = u1
        self._u2 = u2
        self._d1 = d1
        self._d2 = d2
        self._tol = tol

    def build(self):
        self._performed = True
        if not CheckGeom.is_curve_like(self._c) or self._n <= 0:
            return None

        adp_crv = GeomAdaptor_Curve(self._c.handle)

        # Set u1 and u2
        if self._u1 is None:
            self._u1 = adp_crv.FirstParameter()
        if self._u2 is None:
            self._u2 = adp_crv.LastParameter()

        # Adjust u1 and u2 if d1 or d2 != 0
        if self._d1 is not None:
            ap = GCPnts_AbscissaPoint(self._tol, adp_crv, self._d1, self._u1)
            if ap.IsDone():
                self._u1 = ap.Parameter()
        if self._d1 is not None:
            ap = GCPnts_AbscissaPoint(self._tol, adp_crv, self._d2, self._u2)
            if ap.IsDone():
                self._u2 = ap.Parameter()

        # Create uniform abscissa
        ua = GCPnts_UniformAbscissa(adp_crv, self._n, self._u1, self._u2,
                                    self._tol)
        if not ua.IsDone():
            return None

        # Gather results
        npts = ua.NbPoints()
        pnts = []
        prms = []
        for i in range(1, npts + 1):
            u = ua.Parameter(i)
            p = Point()
            adp_crv.D0(u, p)
            pnts.append(p)
            prms.append(u)
        self._set_results('npts', npts)
        self._set_results('prms', prms)
        self._set_results('pnts', pnts)

        # Point spacing
        if npts > 1:
            ds = pnts[0].distance(pnts[1])
            self._set_results('ds', ds)

        self._success = True

    @property
    def npts(self):
        """
        :return: The number of points.
        :rtype: int
        """
        return self._get_results('npts')

    @property
    def points(self):
        """
        :return: The points.
        :rtype: list[afem.geometry.entities.Point]
        """
        return self._get_results('pnts')

    @property
    def parameters(self):
        """
        :return: The parameters.
        :rtype: list[float]
        """
        return self._get_results('prms')

    @property
    def spacing(self):
        """
        :return: The spacing between the first and second points if there
            are more than one point. Otherwise *None*.
        :rtype: float or None
        """
        return self._get_results('ds')


class PointsAlongCurveByDistance(GeomBuilder):
    """
    Create points along a curve by distance between points. The points will
    be equidistant.

    :param curve_like c: The curve.
    :param float maxd: The maximum allowed spacing between points. The
        actual spacing will be adjusted to not to exceed this value.
    :param float u1: The parameter of the first point (default=c.u1).
    :param float u2: The parameter of the last point (default=c.u2).
    :param float d1: An offset distance for the first point.
    :param float d2: An offset distance for the last point.
    :param int nmin: Minimum number of points to create.
    :param float tol: Tolerance.

    Usage:

    >>> from afem.geometry import LineByPoints, Point, PointsAlongCurveByDistance
    >>> p1 = Point()
    >>> p2 = Point(10., 0., 0.)
    >>> line = LineByPoints(p1, p2).line
    >>> builder = PointsAlongCurveByDistance(line, 5., 0., 10.)
    >>> assert builder.success
    >>> builder.npts
    3
    >>> builder.points
    [Point(0.0, 0.0, 0.0), Point(5.0, 0.0, 0.0), Point(10.0, 0.0, 0.0)]
    >>> builder.parameters
    [0.0, 5.0, 10.0]
    """

    def __init__(self, c, maxd, u1=None, u2=None, d1=None, d2=None, nmin=0,
                 tol=1.0e-7):
        super(PointsAlongCurveByDistance, self).__init__()
        self._c = c
        self._maxd = float(maxd)
        self._u1 = u1
        self._u2 = u2
        self._d1 = d1
        self._d2 = d2
        self._nmin = int(nmin)
        self._tol = tol

    def build(self):
        self._performed = True
        if not CheckGeom.is_curve_like(self._c) or self._maxd <= 0:
            return None

        adp_crv = GeomAdaptor_Curve(self._c.handle)

        # Set u1 and u2
        if self._u1 is None:
            self._u1 = adp_crv.FirstParameter()
        if self._u2 is None:
            self._u2 = adp_crv.LastParameter()

        # Adjust u1 and u2 if d1 or d2 != 0
        if self._d1 is not None:
            ap = GCPnts_AbscissaPoint(self._tol, adp_crv, self._d1, self._u1)
            if ap.IsDone():
                self._u1 = ap.Parameter()
        if self._d1 is not None:
            ap = GCPnts_AbscissaPoint(self._tol, adp_crv, self._d2, self._u2)
            if ap.IsDone():
                self._u2 = ap.Parameter()

        # Determine number of points
        arc_length = GCPnts_AbscissaPoint.Length(adp_crv, self._u1,
                                                 self._u2, self._tol)
        n = ceil(arc_length / self._maxd) + 1
        if n < self._nmin:
            n = self._nmin

        # Create uniform abscissa
        ua = GCPnts_UniformAbscissa(adp_crv, n, self._u1, self._u2, self._tol)
        if not ua.IsDone():
            return None

        # Gather results
        npts = ua.NbPoints()
        pnts = []
        prms = []
        for i in range(1, npts + 1):
            u = ua.Parameter(i)
            p = Point()
            adp_crv.D0(u, p)
            pnts.append(p)
            prms.append(u)
        self._set_results('npts', npts)
        self._set_results('prms', prms)
        self._set_results('pnts', pnts)

        # Point spacing
        if npts > 1:
            ds = pnts[0].distance(pnts[1])
            self._set_results('ds', ds)

        self._success = True

    @property
    def npts(self):
        """
        :return: The number of points.
        :rtype: int
        """
        return self._get_results('npts')

    @property
    def points(self):
        """
        :return: The points.
        :rtype: list[afem.geometry.entities.Point]
        """
        return self._get_results('pnts')

    @property
    def parameters(self):
        """
        :return: The parameters.
        :rtype: list[float]
        """
        return self._get_results('prms')

    @property
    def spacing(self):
        """
        :return: The spacing between the first and second points if there
            are more than one point. Otherwise *None*.
        :rtype: float or None
        """
        return self._get_results('ds')


class DirectionByXYZ(GeomBuilder):
    """
    Create a direction (i.e., unit vector) by x-, y-, and z-components.

    :param float x: x-component.
    :param float y: y-component.
    :param float z: z-component.

    Usage:

    >>> from afem.geometry import DirectionByXYZ
    >>> DirectionByXYZ(10., 0., 0.).direction
    Direction(1.0, 0.0, 0.0)
    """

    def __init__(self, x=0., y=0., z=0.):
        super(DirectionByXYZ, self).__init__()
        self._x = float(x)
        self._y = float(y)
        self._z = float(z)

    def build(self):
        self._performed = True

        d = Direction(self._x, self._y, self._z)
        self._set_results('d', d)
        self._success = isinstance(d, Direction)

    @property
    def direction(self):
        """
        :return: The direction.
        :rtype: afem.geometry.entities.Direction
        """
        return self._get_results('d')


class DirectionByArray(DirectionByXYZ):
    """
    Create a direction (i.e., unit vector) from an array-like object.

    :param array_like xyz: Array-like object defining xyz-components.

    Usage:

    >>> from afem.geometry import DirectionByArray
    >>> DirectionByArray([10., 0., 0.]).direction
    Direction(1.0, 0.0, 0.0)
    """

    def __init__(self, xyz=(1., 0., 0.)):
        assert len(xyz) == 3, "Invalid array size in DirectionByArray"
        x, y, z = xyz
        super(DirectionByArray, self).__init__(x, y, z)


class DirectionByPoints(DirectionByArray):
    """
    Create a direction (i.e., unit vector) between two points.

    :param point_like p1: The first point.
    :param point_like p2: The last point.

    Usage:

    >>> from afem.geometry import DirectionByPoints, Point
    >>> p1 = Point()
    >>> p2 = Point(10., 0., 0.)
    >>> DirectionByPoints(p1, p2).direction
    Direction(1.0, 0.0, 0.0)
    """

    def __init__(self, p1, p2):
        p1 = CheckGeom.to_point(p1)
        p2 = CheckGeom.to_point(p2)
        super(DirectionByPoints, self).__init__(p2 - p1)


class VectorByXYZ(GeomBuilder):
    """
    Create a vector by x-, y-, and z-components.

    :param float x: x-component.
    :param float y: y-component.
    :param float z: z-component.

    Usage:

    >>> from afem.geometry import VectorByXYZ
    >>> VectorByXYZ(1., 2., 3.).vector
    Vector(1.0, 2.0, 3.0)
    """

    def __init__(self, x=0., y=0., z=0.):
        super(VectorByXYZ, self).__init__()
        self._x = float(x)
        self._y = float(y)
        self._z = float(z)

    def build(self):
        self._performed = True

        v = Vector(self._x, self._y, self._z)
        self._set_results('v', v)
        self._success = isinstance(v, Direction)

    @property
    def vector(self):
        """
        :return: The vector.
        :rtype: afem.geometry.entities.Vector
        """
        return self._get_results('v')


class VectorByArray(VectorByXYZ):
    """
    Create a vector from an array-like object.

    :param array_like xyz: Array-like object defining xyz-components.

    Usage:

    >>> from afem.geometry import VectorByArray
    >>> VectorByArray([1., 2., 3.]).vector
    Vector(1.0, 2.0, 3.0)
    """

    def __init__(self, xyz=(1., 0., 0.)):
        assert len(xyz) == 3, "Invalid array size in VectorByArray"
        x, y, z = xyz
        super(VectorByArray, self).__init__(x, y, z)


class VectorByPoints(GeomBuilder):
    """
    Create a vecotr between two points.

    :param point_like p1: The first point.
    :param point_like p2: The last point.

    >>> from afem.geometry import Point, VectorByPoints
    >>> p1 = Point()
    >>> p2 = Point(1., 2., 3.)
    >>> VectorByPoints(p1, p2).vector
    Vector(1.0, 2.0, 3.0)
    """

    def __init__(self, p1, p2):
        super(VectorByPoints, self).__init__()
        self._p1 = CheckGeom.to_point(p1)
        self._p2 = CheckGeom.to_point(p2)

    def build(self):
        self._performed = True

        if None in [self._p1, self._p2]:
            return None

        v = Vector(self._p1, self._p2)
        self._set_results('v', v)
        self._success = isinstance(v, Vector)

    @property
    def vector(self):
        """
        :return: The vector.
        :rtype: afem.geometry.entities.Vector
        """
        return self._get_results('v')


class LineByVector(GeomBuilder):
    """
    Create a line by an origin and a vector.

    :param point_like p: Origin of line.
    :param vector_like v: Direction of line.

    Usage:

    >>> from afem.geometry import Point, Vector, LineByVector
    >>> p = Point()
    >>> v = Vector(1., 0., 0.)
    >>> builder = LineByVector(p, v)
    >>> assert builder.success
    >>> line = builder.line
    """

    def __init__(self, p, v):
        super(LineByVector, self).__init__()
        self._p = CheckGeom.to_point(p)
        self._d = CheckGeom.to_direction(v)
        assert isinstance(self._p, Point), "Invalid point in LineByVector"
        assert isinstance(self._d, Direction), "Invalid vector in LineByVector"

    def build(self):
        self._performed = True

        lin = Line(self._p, self._d)
        self._set_results('lin', lin)
        self._success = isinstance(lin, Line)

    @property
    def line(self):
        """
        :return: The line.
        :rtype: afem.geometry.entities.Line
        """
        return self._get_results('lin')


class LineByPoints(GeomBuilder):
    """
    Create a line through two points.

    :param point_like p1: The first point.
    :param point_like p2: The last point.

    Usage:

    >>> from afem.geometry import LineByPoints, Point
    >>> p1 = Point()
    >>> p2 = Point(10., 0. ,0.)
    >>> builder = LineByPoints(p1, p2)
    >>> assert builder.success
    >>> line = builder.line
    """

    def __init__(self, p1, p2):
        super(LineByPoints, self).__init__()
        self._p1 = CheckGeom.to_point(p1)
        self._p2 = CheckGeom.to_point(p2)
        assert isinstance(self._p1, Point), ("Invalid first point in "
                                             "LineByPoints")
        assert isinstance(self._p2, Point), ("Invalid last point in "
                                             "LineByPoints")

    def build(self):
        self._performed = True

        d = DirectionByArray(self._p2 - self._p1).direction

        lin = Line(self._p1, d)
        self._set_results('lin', lin)
        self._success = isinstance(lin, Line)

    @property
    def line(self):
        """
        :return: The line.
        :rtype: afem.geometry.entities.Line
        """
        return self._get_results('lin')


class NurbsCurveByData(GeomBuilder):
    """
    Create a NURBS curve by data.

    :param list[point_like] cp: Control points of curve.
    :param list[float] knots: Knot vector of curve.
    :param list[int] mult: Multiplicities of curve knot vector.
    :param list[float] weights: Weights of control points.
    :param bool is_periodic: Flag for curve periodicity.

    Usage:

    >>> from afem.geometry import NurbsCurveByData, Point
    >>> cp = [Point(), Point(10., 0., 0.)]
    >>> knots = [0., 1.]
    >>> mult = [2, 2]
    >>> p = 1
    >>> builder = NurbsCurveByData(cp, knots, mult, p)
    >>> assert builder.success
    >>> c = builder.curve
    >>> c.knots
    array([ 0.,  1.])
    >>> c.mult
    array([2, 2])
    >>> c.eval(0.5)
    Point(5.0, 0.0, 0.0)
    """

    def __init__(self, cp, knots, mult, p, weights=None, is_periodic=False):
        super(NurbsCurveByData, self).__init__()
        self._cp = cp
        self._knots = knots
        self._mult = mult
        self._p = int(p)
        self._weights = weights
        self._is_periodic = is_periodic

    def build(self):
        self._performed = True

        tcol_cp = to_tcolgp_array1_pnt(self._cp)
        tcol_knots = to_tcolstd_array1_real(self._knots)
        tcol_mult = to_tcolstd_array1_integer(self._mult)
        if self._weights is None:
            self._weights = [1.] * tcol_cp.Length()
        tcol_weights = to_tcolstd_array1_real(self._weights)

        c = NurbsCurve(tcol_cp, tcol_weights, tcol_knots, tcol_mult, self._p,
                       self._is_periodic)
        self._set_results('c', c)
        self._success = isinstance(c, NurbsCurve)

    @property
    def curve(self):
        """
        :return: The NURBS curve.
        :rtype: afem.geometry.entities.NurbsCurve
        """
        return self._get_results('c')


class NurbsCurveByInterp(GeomBuilder):
    """
    Create a cubic curve by interpolating points.

    :param list[point_like] qp: List of points to interpolate.
    :param bool is_periodic: Flag for curve periodicity. If *True* the curve
        will be periodic and closed.
    :param vector_like v1: Tangent to match at first point.
    :param vector_like v2: Tangent to match at last point.
    :param float tol: Tolerance used to check for coincident points and the
        magnitude of end vectors.

    For more information see GeomAPI_Interpolate_.

    .. _GeomAPI_Interpolate: https://www.opencascade.com/doc/occt-7.1.0/refman/html/class_geom_a_p_i___interpolate.html

    Usage:

    >>> from afem.geometry import NurbsCurveByInterp, Point
    >>> qp = [Point(), Point(5., 5., 0.), Point(10., 0., 0.)]
    >>> builder = NurbsCurveByInterp(qp)
    >>> assert builder.success
    >>> c = builder.curve
    >>> c.p
    2
    >>> c.u1
    0.0
    >>> c.u2
    14.142135623730951
    >>> c.eval(5.)
    Point(3.5355339059327373, 4.571067811865475, 0.0)
    """

    def __init__(self, qp, is_periodic=False, v1=None, v2=None, tol=1.0e-7):
        super(NurbsCurveByInterp, self).__init__()
        self._qp = qp
        self._v1 = v1
        self._v2 = v2
        self._is_periodic = is_periodic
        self._tol = tol

    def build(self):
        self._performed = True

        tcol_hpnts = to_tcolgp_harray1_pnt(self._qp)
        # TODO Remove use of GetHandle
        interp = GeomAPI_Interpolate(tcol_hpnts.GetHandle(),
                                     self._is_periodic, self._tol)

        if None not in [self._v1, self._v2]:
            v1 = CheckGeom.to_vector(self._v1)
            v2 = CheckGeom.to_vector(self._v2)
            if v1 and v2:
                interp.Load(v1, v2)

        interp.Perform()
        if not interp.IsDone():
            return None

        # TODO Remove use of GetObject
        occ_crv = interp.Curve().GetObject()
        c = _create_nurbs_curve_from_occ(occ_crv)
        self._set_results('c', c)
        self._success = isinstance(c, NurbsCurve)

    @property
    def curve(self):
        """
        :return: The NURBS curve.
        :rtype: afem.geometry.entities.NurbsCurve
        """
        return self._get_results('c')


class NurbsCurveByApprox(GeomBuilder):
    """
    Create a NURBS curve by approximating points.

    :param list[point_like] qp: List of points to approximate.
    :param int dmin: Minimum degree.
    :param int dmax: Maximum degree.
    :param str continuity: Desired continuity of curve ('C0', 'G1', 'C1',
        'G2', 'C2', 'C3').
    :param str parm_type: Parametrization type ('uniform', 'chord',
        'centripetal').
    :param float tol: The tolerance used for approximation. The distance
        from the points to the resulting curve should be lower than *tol*.

    For more information see GeomAPI_PointsToBSpline_.

    .. _GeomAPI_PointsToBSpline: https://www.opencascade.com/doc/occt-7.1.0/refman/html/class_geom_a_p_i___points_to_b_spline.html

    Usage:

    >>> from afem.geometry import NurbsCurveByApprox, Point
    >>> qp = [Point(), Point(5., 5., 0.), Point(10., 0., 0.)]
    >>> builder = NurbsCurveByApprox(qp)
    >>> assert builder.success
    >>> c = builder.curve
    >>> c.p
    3
    >>> c.u1
    0.0
    >>> c.u2
    1.0
    >>> c.eval(0.5)
    Point(5.0, 5.0, 0.0)
    """

    def __init__(self, qp, dmin=3, dmax=8, continuity='C2',
                 parm_type='chord', tol=1.0e-3):
        super(NurbsCurveByApprox, self).__init__()
        self._qp = qp
        self._dmin = int(dmin)
        self._dmax = int(dmax)
        self._cont = continuity
        self._parm_type = parm_type
        self._tol = float(tol)

    def build(self):
        self._performed = True

        tcol_pnts = to_tcolgp_array1_pnt(self._qp)

        try:
            cont = _occ_continuity[self._cont.upper()]
        except (KeyError, AttributeError):
            cont = GeomAbs_C2

        try:
            parm_type = _occ_parm_type[self._parm_type.lower()]
        except (KeyError, AttributeError):
            parm_type = Approx_ChordLength

        fit = GeomAPI_PointsToBSpline(tcol_pnts, parm_type, self._dmin,
                                      self._dmax, cont, self._tol)
        if not fit.IsDone():
            return None

        # TODO Remove use of GetObject
        occ_crv = fit.Curve().GetObject()
        c = _create_nurbs_curve_from_occ(occ_crv)
        self._set_results('c', c)
        self._success = isinstance(c, NurbsCurve)

    @property
    def curve(self):
        """
        :return: The NURBS curve.
        :rtype: afem.geometry.entities.NurbsCurve
        """
        return self._get_results('c')


class NurbsCurveByPoints(NurbsCurveByApprox):
    """
    Create a linear curve (i.e., a polyline) between points.

    :param list[point_like] qp: List of points.

    Usage:

    >>> from afem.geometry import NurbsCurveByPoints, Point
    >>> qp = [Point(), Point(5., 5., 0.), Point(10., 0., 0.)]
    >>> builder = NurbsCurveByPoints(qp)
    >>> assert builder.success
    >>> c = builder.curve
    >>> c.p
    1
    >>> c.u1
    0.0
    >>> c.u2
    1.0
    >>> c.eval(0.5)
    Point(5.0, 5.0, 0.0)
    """

    def __init__(self, qp):
        super(NurbsCurveByPoints, self).__init__(qp, 1, 1, 'C0')


class CurveByUIso(GeomBuilder):
    """
    Create an isocurve from a surface at a constant u-parameter.

    :param surface_like s: The surface.
    :param float u: The parameter.

    The following curve types are created for a given surface:

    * Plane -> Line
    * NurbsSurface -> NurbsCurve

    Usage:

    >>> from afem.geometry import CurveByUIso, Direction, Plane, Point
    >>> p0 = Point()
    >>> vn = Direction(0., 0., 1.)
    >>> pln = Plane(p0, vn)
    >>> builder = CurveByUIso(pln, 1.)
    >>> assert builder.success
    >>> builder.is_line
    True
    >>> builder.is_nurbs
    False
    >>> line = builder.curve
    >>> line.eval(0.)
    Point(1.0, 0.0, 0.0)
    """

    def __init__(self, s, u):
        super(CurveByUIso, self).__init__()
        self._s = s
        self._u = float(u)

    def build(self):
        self._performed = True

        if not CheckGeom.is_surface_like(self._s):
            return None

        # TODO Handle other curve types
        hcrv = self._s.UIso(self._u)
        adp_curve = GeomAdaptor_Curve(hcrv)
        if adp_curve.GetType() == GeomAbs_Line:
            gp_lin = adp_curve.Line()
            c = Line(gp_lin)
        elif adp_curve.GetType() == GeomAbs_BSplineCurve:
            occ_crv = adp_curve.BSpline().GetObject()
            c = _create_nurbs_curve_from_occ(occ_crv)
        else:
            return None

        self._set_results('c', c)
        self._success = isinstance(c, Geom_Curve)

    @property
    def is_line(self):
        """
        :return: *True* if the isocurve is a line, *False* if not.
        :rtype: bool
        """
        return isinstance(self.curve, Line)

    @property
    def is_nurbs(self):
        """
        :return: *True* if the isocurve is a NURBS curve, *False* if not.
        :rtype: bool
        """
        return isinstance(self.curve, NurbsCurve)

    @property
    def curve(self):
        """
        :return: The curve.
        :rtype: afem.geometry.entities.Line or afem.geometry.entities.NurbsCurve
        """
        return self._get_results('c')


class CurveByVIso(GeomBuilder):
    """
    Create an isocurve from a surface at a constant v-parameter.

    :param surface_like s: The surface.
    :param float v: The parameter.

    The following curve types are created for a given surface:

    * Plane -> Line
    * NurbsSurface -> NurbsCurve

    Usage:

    >>> from afem.geometry import CurveByVIso, Direction, Plane, Point
    >>> p0 = Point()
    >>> vn = Direction(0., 0., 1.)
    >>> pln = Plane(p0, vn)
    >>> builder = CurveByVIso(pln, 1.)
    >>> assert builder.success
    >>> builder.is_line
    True
    >>> builder.is_nurbs
    False
    >>> line = builder.curve
    >>> line.eval(0.)
    Point(0.0, 1.0, 0.0)
    """

    def __init__(self, s, v):
        super(CurveByVIso, self).__init__()
        self._s = s
        self._v = float(v)

    def build(self):
        self._performed = True

        if not CheckGeom.is_surface_like(self._s):
            return None

        # TODO Handle other curve types
        hcrv = self._s.VIso(self._v)
        adp_curve = GeomAdaptor_Curve(hcrv)
        if adp_curve.GetType() == GeomAbs_Line:
            gp_lin = adp_curve.Line()
            c = Line(gp_lin)
        elif adp_curve.GetType() == GeomAbs_BSplineCurve:
            occ_crv = adp_curve.BSpline().GetObject()
            c = _create_nurbs_curve_from_occ(occ_crv)
        else:
            return None

        self._set_results('c', c)
        self._success = isinstance(c, Geom_Curve)

    @property
    def is_line(self):
        """
        :return: *True* if the isocurve is a line, *False* if not.
        :rtype: bool
        """
        return isinstance(self.curve, Line)

    @property
    def is_nurbs(self):
        """
        :return: *True* if the isocurve is a NURBS curve, *False* if not.
        :rtype: bool
        """
        return isinstance(self.curve, NurbsCurve)

    @property
    def curve(self):
        """
        :return: The curve.
        :rtype: afem.geometry.entities.Line or afem.geometry.entities.NurbsCurve
        """
        return self._get_results('c')


class PlaneByNormal(GeomBuilder):
    """
    Create a plane by an origin and a normal vector.
    """
    pass


class PlaneByAxes(GeomBuilder):
    """
    Create a plane by an origin and basic axes.
    """
    pass


class PlaneByPoints(GeomBuilder):
    """
    Create a plane by three points.
    """
    pass


class PlaneByFit(GeomBuilder):
    """
    Create a plane by fitting points.
    """
    pass


class PlanesAlongCurve(GeomBuilder):
    """
    Create planes along a curve.
    """
    pass


class PlanesBetweenPlanes(GeomBuilder):
    """
    Create planes between two other planes.
    """
    pass


class SurfaceByData(GeomBuilder):
    """
    Create a surface from data.
    """

    def __init__(self):
        super(SurfaceByData, self).__init__()


class SurfaceByInterp(GeomBuilder):
    """
    Create a surface by interpolating curves.
    """
    pass


class SurfaceByFit(GeomBuilder):
    """
    Create a surface by fitting curves.
    """
    pass


if __name__ == "__main__":
    import doctest

    doctest.testmod()
