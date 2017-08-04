from math import ceil, radians

import OCC.BSplCLib as CLib
from OCC.Approx import (Approx_Centripetal, Approx_ChordLength,
                        Approx_IsoParametric)
from OCC.GCPnts import GCPnts_AbscissaPoint, GCPnts_UniformAbscissa
from OCC.Geom import Geom_Curve
from OCC.GeomAPI import (GeomAPI_IntCS, GeomAPI_Interpolate,
                         GeomAPI_PointsToBSpline)
from OCC.GeomAbs import (GeomAbs_BSplineCurve, GeomAbs_C0, GeomAbs_C1,
                         GeomAbs_C2, GeomAbs_C3, GeomAbs_G1, GeomAbs_G2,
                         GeomAbs_Line)
from OCC.GeomAdaptor import GeomAdaptor_Curve
from OCC.GeomFill import (GeomFill_AppSurf, GeomFill_Line,
                          GeomFill_SectionGenerator)
from OCC.GeomPlate import GeomPlate_BuildAveragePlane
from OCC.TColStd import TColStd_Array1OfInteger, TColStd_Array1OfReal
from OCC.TColgp import TColgp_Array1OfPnt
from OCC.gp import gp_Pln, gp_Quaternion, gp_Trsf
from numpy import array, cross, mean, ones, zeros
from numpy.linalg import norm
from scipy.linalg import lu_factor, lu_solve

from afem.geometry.check import CheckGeom
from afem.geometry.entities import *
from afem.geometry.methods.create import create_crv_by_approx_pnts, \
    create_crv_by_interp_pnts, create_isocurve, create_line_from_occ, \
    create_nurbs_curve, \
    create_nurbs_curve_from_occ2d, create_nurbs_surface, \
    create_nurbs_surface_from_occ, create_plane_by_axes, \
    create_plane_by_fit_points, create_plane_by_points, create_plane_on_curve, \
    create_planes_along_curve, create_planes_between_planes, \
    create_point_from_other, create_points_along_curve, \
    create_srf_by_approx_crvs, create_srf_by_interp_crvs
from afem.geometry.utils import (basis_funs, centripetal_parameters,
                                 chord_parameters, dehomogenize_array2d,
                                 find_span, homogenize_array1d,
                                 uniform_parameters)
from afem.occ.utils import (to_np_from_tcolgp_array1_pnt,
                            to_np_from_tcolstd_array1_real,
                            to_tcolgp_array1_pnt, to_tcolgp_array2_pnt,
                            to_tcolgp_harray1_pnt,
                            to_tcolstd_array1_integer,
                            to_tcolstd_array1_real,
                            to_tcolstd_array2_real)

__all__ = ["CreateGeom", "CreatedPoints", "PointByXYZ", "PointByArray",
           "PointFromParameter", "PointsAlongCurveByNumber",
           "PointsAlongCurveByDistance", "DirectionByXYZ", "DirectionByArray",
           "DirectionByPoints", "VectorByXYZ", "VectorByArray",
           "VectorByPoints", "LineByVector", "LineByPoints",
           "NurbsCurveByData", "NurbsCurveByInterp", "NurbsCurveByApprox",
           "NurbsCurveByPoints", "CurveByUIso", "CurveByVIso",
           "PlaneByNormal", "PlaneByAxes", "PlaneByPoints", "PlaneByApprox",
           "PlanesAlongCurveByNumber", "PlanesAlongCurveByDistance",
           "PlanesBetweenPlanesByNumber", "PlanesBetweenPlanesByDistance",
           "NurbsSurfaceByData", "NurbsSurfaceByInterp",
           "NurbsSurfaceByApprox"]

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


def create_nurbs_curve_from_occ(crv):
    """
    Create a NURBS curve from an OCC Geom_BSplineCurve.

    :param OCC.Geom.Geom_BSplineCurve crv: The OCC curve.

    :return: The AFEM NurbsCurve.
    :rtype: afem.geometry.entities.NurbsCurve
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
        self._results = {}

    def _set_results(self, key, value):
        self._results[key] = value

    def _get_results(self, key):
        try:
            return self._results[key]
        except KeyError:
            return None


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

        # Build
        p = Point(x, y, z)
        self._set_results('p', p)
        self._success = True

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

    :raise ValueError: If `len(xyz) != 3`.

    Usage:

    >>> from afem.geometry import PointByArray
    >>> PointByArray([1., 2., 3.]).point
    Point(1.0, 2.0, 3.0)
    """

    def __init__(self, xyz=(0., 0., 0.)):
        if len(xyz) != 3:
            msg = 'Invalid array size.'
            raise ValueError(msg)
        x, y, z = xyz
        super(PointByArray, self).__init__(x, y, z)


class PointFromParameter(GeomBuilder):
    """
    Create a point along a curve at a specified distance from a parameter.

    :param curve_like c: The curve.
    :param float u0: The initial parameter.
    :param float ds: The distance along the curve from the given parameter.
    :param float tol: Tolerance.

    :raise RuntimeError: If OCC method fails.

    Usage:

    >>> from afem.geometry import LineByPoints, Point, PointFromParameter
    >>> p1 = Point()
    >>> p2 = Point(10., 0., 0.)
    >>> line = LineByPoints(p1, p2).line
    >>> builder = PointFromParameter(line, 0., 1.)
    >>> builder.point
    Point(1.0, 0.0, 0.0)
    >>> builder.parameter
    1.0
    """

    def __init__(self, c, u0, ds, tol=1.0e-7):
        super(PointFromParameter, self).__init__()

        # Build
        adp_curve = GeomAdaptor_Curve(c.handle)

        ap = GCPnts_AbscissaPoint(tol, adp_curve, ds, u0)
        if not ap.IsDone():
            msg = "GCPnts_AbscissaPoint failed."
            raise RuntimeError(msg)

        u = ap.Parameter()
        self._set_results('u', u)
        p = c.eval(u)
        self._set_results('p', p)
        self._success = True

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
    :param float d1: An offset distance for the first point. This is typically
        a positive number indicating a distance from *u1* towards *u2*.
    :param float d2: An offset distance for the last point. This is typically
        a negative number indicating a distance from *u2* towards *u1*.
    :param float tol: Tolerance.

    :raise RuntimeError: If OCC method fails.
    :raise RuntimeWarning: If OCC methods fails to update *u1* if *d1* is
        provided.
    :raise RuntimeWarning: If OCC methods fails to update *u2* if *d2* is
        provided.

    For more information see GCPnts_UniformAbscissa_.

    .. _GCPnts_UniformAbscissa: https://www.opencascade.com/doc/occt-7.1.0/refman/html/class_g_c_pnts___uniform_abscissa.html

    Usage:

    >>> from afem.geometry import LineByPoints, Point, PointsAlongCurveByNumber
    >>> p1 = Point()
    >>> p2 = Point(10., 0., 0.)
    >>> line = LineByPoints(p1, p2).line
    >>> builder = PointsAlongCurveByNumber(line, 3, 0., 10.)
    >>> builder.npts
    3
    >>> builder.points
    [Point(0.0, 0.0, 0.0), Point(5.0, 0.0, 0.0), Point(10.0, 0.0, 0.0)]
    >>> builder.parameters
    [0.0, 5.0, 10.0]
    """

    def __init__(self, c, n, u1=None, u2=None, d1=None, d2=None, tol=1.0e-7):
        super(PointsAlongCurveByNumber, self).__init__()
        n = int(n)

        # Build
        adp_crv = GeomAdaptor_Curve(c.handle)

        # Set u1 and u2
        if u1 is None:
            u1 = adp_crv.FirstParameter()
        if u2 is None:
            u2 = adp_crv.LastParameter()

        # Adjust u1 and u2 if d1 or d2 != 0
        if d1 is not None:
            ap = GCPnts_AbscissaPoint(tol, adp_crv, d1, u1)
            if ap.IsDone():
                u1 = ap.Parameter()
            else:
                msg = "GCPnts_AbscissaPoint failed for u1."
                raise RuntimeWarning(msg)
        if d1 is not None:
            ap = GCPnts_AbscissaPoint(tol, adp_crv, d2, u2)
            if ap.IsDone():
                u2 = ap.Parameter()
            else:
                msg = "GCPnts_AbscissaPoint failed for u2."
                raise RuntimeWarning(msg)

        # Create uniform abscissa
        ua = GCPnts_UniformAbscissa(adp_crv, n, u1, u2, tol)
        if not ua.IsDone():
            msg = "GCPnts_UniformAbscissa failed."
            raise RuntimeError(msg)

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
    be equidistant. This method calculates the number of points given the
    curve length and then uses :class:`.PointsAlongCurveByNumber`.

    :param curve_like c: The curve.
    :param float maxd: The maximum allowed spacing between points. The
        actual spacing will be adjusted to not to exceed this value.
    :param float u1: The parameter of the first point (default=c.u1).
    :param float u2: The parameter of the last point (default=c.u2).
    :param float d1: An offset distance for the first point. This is typically
        a positive number indicating a distance from *u1* towards *u2*.
    :param float d2: An offset distance for the last point. This is typically
        a negative number indicating a distance from *u2* towards *u1*.
    :param int nmin: Minimum number of points to create.
    :param float tol: Tolerance.

    :raise RuntimeError: If OCC method fails.
    :raise RuntimeWarning: If OCC methods fails to update *u1* if *d1* is
        provided.
    :raise RuntimeWarning: If OCC methods fails to update *u2* if *d2* is
        provided.

    Usage:

    >>> from afem.geometry import LineByPoints, Point, PointsAlongCurveByDistance
    >>> p1 = Point()
    >>> p2 = Point(10., 0., 0.)
    >>> line = LineByPoints(p1, p2).line
    >>> builder = PointsAlongCurveByDistance(line, 5., 0., 10.)
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
        maxd = float(maxd)
        nmin = int(nmin)

        # Build
        adp_crv = GeomAdaptor_Curve(c.handle)

        # Set u1 and u2
        if u1 is None:
            u1 = adp_crv.FirstParameter()
        if u2 is None:
            u2 = adp_crv.LastParameter()

        # Adjust u1 and u2 if d1 or d2 != 0
        if d1 is not None:
            ap = GCPnts_AbscissaPoint(tol, adp_crv, d1, u1)
            if ap.IsDone():
                u1 = ap.Parameter()
            else:
                msg = "GCPnts_AbscissaPoint failed for u1."
                raise RuntimeWarning(msg)
        if d1 is not None:
            ap = GCPnts_AbscissaPoint(tol, adp_crv, d2, u2)
            if ap.IsDone():
                u2 = ap.Parameter()
            else:
                msg = "GCPnts_AbscissaPoint failed for u2."
                raise RuntimeWarning(msg)

        # Determine number of points
        arc_length = GCPnts_AbscissaPoint.Length(adp_crv, u1, u2, tol)
        n = ceil(arc_length / maxd) + 1
        if n < nmin:
            n = nmin

        # Create uniform abscissa
        ua = GCPnts_UniformAbscissa(adp_crv, n, u1, u2, tol)
        if not ua.IsDone():
            msg = "GCPnts_UniformAbscissa failed."
            raise RuntimeError(msg)

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

        # Build
        d = Direction(x, y, z)
        self._set_results('d', d)
        self._success = True

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

    :raise ValueError: If `len(xyz) != 3`.

    Usage:

    >>> from afem.geometry import DirectionByArray
    >>> DirectionByArray([10., 0., 0.]).direction
    Direction(1.0, 0.0, 0.0)
    """

    def __init__(self, xyz=(1., 0., 0.)):
        if len(xyz) != 3:
            msg = 'Invalid array size.'
            raise ValueError(msg)
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

        # Build
        v = Vector(x, y, z)
        self._set_results('v', v)
        self._success = True

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

    :raise ValueError: If `len(xyz) != 3`.

    Usage:

    >>> from afem.geometry import VectorByArray
    >>> VectorByArray([1., 2., 3.]).vector
    Vector(1.0, 2.0, 3.0)
    """

    def __init__(self, xyz=(1., 0., 0.)):
        if len(xyz) != 3:
            msg = 'Invalid array size.'
            raise ValueError(msg)
        x, y, z = xyz
        super(VectorByArray, self).__init__(x, y, z)


class VectorByPoints(GeomBuilder):
    """
    Create a vecotr between two points.

    :param point_like p1: The first point.
    :param point_like p2: The last point.

    :raise TypeError: If *p1* or *p2* cannot be converted to a :class:`.Point`.

    >>> from afem.geometry import Point, VectorByPoints
    >>> p1 = Point()
    >>> p2 = Point(1., 2., 3.)
    >>> VectorByPoints(p1, p2).vector
    Vector(1.0, 2.0, 3.0)
    """

    def __init__(self, p1, p2):
        super(VectorByPoints, self).__init__()

        # Build
        p1 = CheckGeom.to_point(p1)
        p2 = CheckGeom.to_point(p2)
        if None in [p1, p2]:
            msg = 'Invalid points provided.'
            raise TypeError(msg)

        v = Vector(p1, p2)
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
    :param vector_like d: Direction of line.

    :raise TypeError: If *p* cannot be converted to a :class:`.Point`
    :raise TypeError: If *d* cannot be converted to a :class:`.Direction`

    Usage:

    >>> from afem.geometry import Point, LineByVector
    >>> p = Point()
    >>> d = Direction(1., 0., 0.)
    >>> builder = LineByVector(p, d)
    >>> line = builder.line
    """

    def __init__(self, p, d):
        super(LineByVector, self).__init__()

        # Build
        p = CheckGeom.to_point(p)
        d = CheckGeom.to_direction(d)
        if not isinstance(p, Point):
            msg = "Invalid point."
            raise TypeError(msg)
        if not isinstance(d, Direction):
            msg = "Invalid direction."
            raise TypeError(msg)

        lin = Line(p, d)
        self._set_results('lin', lin)
        self._success = True

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

    :raise TypeError: If *p1* or *p2* cannot be converted to a :class:`.Point`.

    Usage:

    >>> from afem.geometry import LineByPoints, Point
    >>> p1 = Point()
    >>> p2 = Point(10., 0. ,0.)
    >>> builder = LineByPoints(p1, p2)
    >>> line = builder.line
    """

    def __init__(self, p1, p2):
        super(LineByPoints, self).__init__()

        # Build
        p1 = CheckGeom.to_point(p1)
        p2 = CheckGeom.to_point(p2)
        if not isinstance(p1, Point):
            msg = "Invalid first point."
            raise TypeError(msg)
        if not isinstance(p2, Point):
            msg = "Invalid second point."
            raise TypeError(msg)

        d = DirectionByArray(p2 - p1).direction

        lin = Line(p1, d)
        self._set_results('lin', lin)
        self._success = True

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

    :param list[point_like] cp: Control points.
    :param list[float] knots: Knot vector.
    :param list[int] mult: Multiplicities of knot vector.
    :param int p: Degree.
    :param list[float] weights: Weights of control points.
    :param bool is_periodic: Flag for periodicity.

    Usage:

    >>> from afem.geometry import NurbsCurveByData, Point
    >>> cp = [Point(), Point(10., 0., 0.)]
    >>> knots = [0., 1.]
    >>> mult = [2, 2]
    >>> p = 1
    >>> builder = NurbsCurveByData(cp, knots, mult, p)
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
        p = int(p)

        # Build
        tcol_cp = to_tcolgp_array1_pnt(cp)
        tcol_knots = to_tcolstd_array1_real(knots)
        tcol_mult = to_tcolstd_array1_integer(mult)
        if weights is None:
            weights = [1.] * tcol_cp.Length()
        tcol_weights = to_tcolstd_array1_real(weights)

        c = NurbsCurve(tcol_cp, tcol_weights, tcol_knots, tcol_mult, p,
                       is_periodic)

        self._set_results('c', c)
        self._success = True

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

    :raise RuntimeError: If OCC method fails.

    For more information see GeomAPI_Interpolate_.

    .. _GeomAPI_Interpolate: https://www.opencascade.com/doc/occt-7.1.0/refman/html/class_geom_a_p_i___interpolate.html

    Usage:

    >>> from afem.geometry import NurbsCurveByInterp, Point
    >>> qp = [Point(), Point(5., 5., 0.), Point(10., 0., 0.)]
    >>> builder = NurbsCurveByInterp(qp)
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

        # Build
        tcol_hpnts = to_tcolgp_harray1_pnt(qp)
        # TODO Remove use of GetHandle
        interp = GeomAPI_Interpolate(tcol_hpnts.GetHandle(),
                                     is_periodic, tol)

        if None not in [v1, v2]:
            v1 = CheckGeom.to_vector(v1)
            v2 = CheckGeom.to_vector(v2)
            if v1 and v2:
                interp.Load(v1, v2)

        interp.Perform()
        if not interp.IsDone():
            msg = "GeomAPI_Interpolate failed."
            raise RuntimeError(msg)

        # TODO Remove use of GetObject
        occ_crv = interp.Curve().GetObject()
        c = create_nurbs_curve_from_occ(occ_crv)
        self._set_results('c', c)
        self._success = True

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

    :raise RuntimeError: If OCC method fails to interpolate the points with
        a curve.

    For more information see GeomAPI_PointsToBSpline_.

    .. _GeomAPI_PointsToBSpline: https://www.opencascade.com/doc/occt-7.1.0/refman/html/class_geom_a_p_i___points_to_b_spline.html

    Usage:

    >>> from afem.geometry import NurbsCurveByApprox, Point
    >>> qp = [Point(), Point(5., 5., 0.), Point(10., 0., 0.)]
    >>> builder = NurbsCurveByApprox(qp)
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

        # Build
        dmin = int(dmin)
        dmax = int(dmax)
        tcol_pnts = to_tcolgp_array1_pnt(qp)

        try:
            cont = _occ_continuity[continuity.upper()]
        except (KeyError, AttributeError):
            cont = GeomAbs_C2

        try:
            parm_type = _occ_parm_type[parm_type.lower()]
        except (KeyError, AttributeError):
            parm_type = Approx_ChordLength

        fit = GeomAPI_PointsToBSpline(tcol_pnts, parm_type, dmin,
                                      dmax, cont, tol)
        if not fit.IsDone():
            msg = "GeomAPI_PointsToBSpline failed."
            raise RuntimeError(msg)

        # TODO Remove use of GetObject
        occ_crv = fit.Curve().GetObject()
        c = create_nurbs_curve_from_occ(occ_crv)
        self._set_results('c', c)
        self._success = True

    @property
    def curve(self):
        """
        :return: The NURBS curve.
        :rtype: afem.geometry.entities.NurbsCurve
        """
        return self._get_results('c')


class NurbsCurveByPoints(NurbsCurveByApprox):
    """
    Create a linear curve (i.e., a polyline) between points. This method uses
    :class:`.NurbsCurveByApprox` to fit a linear curve.

    :param list[point_like] qp: List of points.

    Usage:

    >>> from afem.geometry import NurbsCurveByPoints, Point
    >>> qp = [Point(), Point(5., 5., 0.), Point(10., 0., 0.)]
    >>> builder = NurbsCurveByPoints(qp)
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

    :raise RuntimeError: If an unsupported curve type is created.

    The following curve types are created for a given surface:

    * Plane -> Line
    * NurbsSurface -> NurbsCurve

    Usage:

    >>> from afem.geometry import CurveByUIso, Direction, Plane, Point
    >>> p0 = Point()
    >>> vn = Direction(0., 0., 1.)
    >>> pln = Plane(p0, vn)
    >>> builder = CurveByUIso(pln, 1.)
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

        # Build
        # TODO Handle other curve types
        hcrv = s.UIso(u)
        adp_curve = GeomAdaptor_Curve(hcrv)
        if adp_curve.GetType() == GeomAbs_Line:
            gp_lin = adp_curve.Line()
            c = Line(gp_lin)
        elif adp_curve.GetType() == GeomAbs_BSplineCurve:
            occ_crv = adp_curve.BSpline().GetObject()
            c = create_nurbs_curve_from_occ(occ_crv)
        else:
            msg = 'Curve type not yet supported.'
            raise RuntimeError(msg)

        self._set_results('c', c)
        self._success = True

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

    :raise RuntimeError: If an unsupported curve type is created.

    The following curve types are created for a given surface:

    * Plane -> Line
    * NurbsSurface -> NurbsCurve

    Usage:

    >>> from afem.geometry import CurveByVIso, Direction, Plane, Point
    >>> p0 = Point()
    >>> vn = Direction(0., 0., 1.)
    >>> pln = Plane(p0, vn)
    >>> builder = CurveByVIso(pln, 1.)
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

        # Build
        # TODO Handle other curve types
        hcrv = s.VIso(v)
        adp_curve = GeomAdaptor_Curve(hcrv)
        if adp_curve.GetType() == GeomAbs_Line:
            gp_lin = adp_curve.Line()
            c = Line(gp_lin)
        elif adp_curve.GetType() == GeomAbs_BSplineCurve:
            occ_crv = adp_curve.BSpline().GetObject()
            c = create_nurbs_curve_from_occ(occ_crv)
        else:
            msg = 'Curve type not yet supported.'
            raise RuntimeError(msg)

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

    :param point_like origin: The origin.
    :param vector_like vnorm: The normal vector.

    :raise TypeError: If *origin* cannot be converted to :class:`.Point`.
    :raise TypeError: If *vnorm* cannot be converted to :class:`.Direction`.

    Usage:

    >>> from afem.geometry import PlaneByNormal
    >>> builder = PlaneByNormal((0., 0., 0.), (0., 0., 1.))
    >>> pln = builder.plane
    >>> pln.eval(1., 1.)
    Point(1.0, 1.0, 0.0)
    """

    def __init__(self, origin=(0., 0., 0.), vnorm=(0., 0., 1.)):
        super(PlaneByNormal, self).__init__()

        # Build
        p0 = CheckGeom.to_point(origin)
        vn = CheckGeom.to_direction(vnorm)
        if not isinstance(p0, Point):
            msg = 'Invalid point.'
            raise TypeError(msg)
        if not isinstance(vn, Direction):
            msg = 'Invalid normal vector.'
            raise TypeError(msg)

        pln = Plane(p0, vn)
        self._set_results('pln', pln)
        self._success = True

    @property
    def plane(self):
        """
        :return: The plane.
        :rtype: afem.geometry.entities.Plane
        """
        return self._get_results('pln')


class PlaneByAxes(GeomBuilder):
    """
    Create a plane by an origin and basic axes.

    :param point_like origin: The origin.
    :param vector_like axes: The axes ('xy', 'xz', 'yz').

    :raise TypeError: If *origin* cannot be converted to :class:`.Point`.
    :raise Value: If *axes* is not a supported option.

    Usage:

    >>> from afem.geometry import PlaneByAxes
    >>> builder = PlaneByAxes((0., 0., 0.), 'xz')
    >>> pln = builder.plane
    >>> pln.eval(1., 1.)
    Point(1.0, 0.0, 1.0)
    """

    def __init__(self, origin=(0., 0., 0.), axes='xz'):
        super(PlaneByAxes, self).__init__()

        # Build
        origin = CheckGeom.to_point(origin)
        axes = axes.lower()
        if not isinstance(origin, Point):
            msg = 'Invalid point.'
            raise TypeError(msg)
        if axes not in ['xy', 'yx', 'xz', 'zx', 'yz', 'zy']:
            msg = 'Invalid axes.'
            raise ValueError(msg)

        if axes in ['xy', 'yx']:
            pln = Plane(origin, Direction(0., 0., 1))
        elif axes in ['yz', 'zy']:
            pln = Plane(origin, Direction(1., 0., 0))
        elif axes in ['xz', 'zx']:
            pln = Plane(origin, Direction(0., 1., 0))
        else:
            msg = 'Unknown axes.'
            raise RuntimeError(msg)

        self._set_results('pln', pln)
        self._success = True

    @property
    def plane(self):
        """
        :return: The plane.
        :rtype: afem.geometry.entities.Plane
        """
        return self._get_results('pln')


class PlaneByPoints(GeomBuilder):
    """
    Create a plane by three points. The points must not be collinear. The
    plane will be defined by (*p2* - *p1*) x (*p3* - *p1*).

    :param point_like p1: First point. This point will be used as the origin.
    :param point_like p2: Second point.
    :param point_like p3: Third point.

    :raise TypeError: if *p1*, *p2*, or *p3* cannot be converted to a
        :class:`.Point`.
    :raise ValueError: If the points are collinear.

    Usage:

    >>> from afem.geometry import PlaneByPoints, Point
    >>> p1 = Point()
    >>> p2 = Point(1., 0., 0.)
    >>> p3 = Point(0., 1., 0.)
    >>> builder = PlaneByPoints(p1, p2, p3)
    >>> pln = builder.plane
    >>> pln.eval(1., 1.)
    Point(1.0, 1.0, 0.0)
    """

    def __init__(self, p1, p2, p3):
        super(PlaneByPoints, self).__init__()

        # Build
        p1 = CheckGeom.to_point(p1)
        if not isinstance(p1, Point):
            msg = 'Invalid point 1.'
            raise TypeError(msg)
        p2 = CheckGeom.to_point(p2)
        if not isinstance(p2, Point):
            msg = 'Invalid point 2.'
            raise TypeError(msg)
        p3 = CheckGeom.to_point(p3)
        if not isinstance(p3, Point):
            msg = 'Invalid point 3.'
            raise TypeError(msg)

        vx = p2 - p1
        v31 = p3 - p1
        vn = cross(vx, v31)
        if norm(vn) <= 1.0e-12:
            msg = 'Points are collinear.'
            raise ValueError(msg)

        n = Direction(*vn)
        vx = Direction(*vx)
        ax = Axis3(p1, n, vx)
        pln = Plane(gp_Pln(ax))

        self._set_results('pln', pln)
        self._success = True

    @property
    def plane(self):
        """
        :return: The plane.
        :rtype: afem.geometry.entities.Plane
        """
        return self._get_results('pln')


class PlaneByApprox(GeomBuilder):
    """
    Create a plane by fitting points. The points must not be collinear.

    :param list[point_like] pnts: Points to fit plane. Should not be collinear.
    :param float tol: Tolerance used to check for collinear points.

    :raise ValueError: If the number of points is less than three.
    :raise RuntimeError: If a plane cannot be fit to the points.

    For more information see GeomPlate_BuildAveragePlane_.

    .. _GeomPlate_BuildAveragePlane: https://www.opencascade.com/doc/occt-7.1.0/refman/html/class_geom_plate___build_average_plane.html

    Usage:

    >>> from afem.geometry import PlaneByApprox, Point
    >>> pnts = [Point(), Point(1., 0., 0.), Point(0., 1., 0.)]
    >>> builder = PlaneByApprox(pnts)
    >>> pln = builder.plane
    """

    def __init__(self, pnts, tol=1.0e-7):
        super(PlaneByApprox, self).__init__()

        # Build
        _pnts = [CheckGeom.to_point(p) for p in pnts if
                 CheckGeom.is_point_like(p)]
        if len(_pnts) < 3:
            msg = "Need at least three points to fit a plane."
            raise ValueError(msg)

        tcol_pnts = to_tcolgp_harray1_pnt(pnts)
        avg_pln = GeomPlate_BuildAveragePlane(tcol_pnts.GetHandle(),
                                              tcol_pnts.Length(),
                                              tol, 1, 1)
        if not avg_pln.IsPlane():
            msg = ('Failed to create a plane by fitting points. Check for '
                   'collinear points.')
            raise RuntimeError(msg)

        # Move to centroid.
        gp_pln = avg_pln.Plane().GetObject().Pln()
        pcg = Point(*mean(pnts, axis=0))
        gp_pln.SetLocation(pcg)
        pln = Plane(gp_pln)

        self._set_results('pln', pln)
        self._success = True

    @property
    def plane(self):
        """
        :return: The plane.
        :rtype: afem.geometry.entities.Plane
        """
        return self._get_results('pln')


class PlanesAlongCurveByNumber(GeomBuilder):
    """
    Create planes along a curve using a specified number. The origin of the
    planes will be equidistant along the curve.

    :param curve_like c: The curve.
    :param int n: Number of planes to create (*n* > 0).
    :param afem.geometry.entities.Plane ref_pln: The normal of this plane
        will be used to define the normal of all planes along the curve. If
        no plane is provided, then the first derivative of the curve will
        define the plane normal.
    :param float u1: The parameter of the first plane (default=c.u1).
    :param float u2: The parameter of the last plane (default=c.u2).
    :param float d1: An offset distance for the first plane. This is typically
        a positive number indicating a distance from *u1* towards *u2*.
    :param float d2: An offset distance for the last plane. This is typically
        a negative number indicating a distance from *u2* towards *u1*.
    :param float tol: Tolerance.

    :raise RuntimeError: If :class:`.PointsAlongCurveByNumber` fails to
        generate points along the curve.

    Usage:

    >>> from afem.geometry import LineByPoints, PlanesAlongCurveByNumber
    >>> line = LineByPoints((0., 0., 0.), (10., 0., 0.)).line
    >>> builder = PlanesAlongCurveByNumber(line, 3, u1=0., u2=10.)
    >>> builder.nplanes
    3
    >>> builder.parameters
    [0.0, 5.0, 10.0]
    >>> builder.spacing
    5.0
    """

    def __init__(self, c, n, ref_pln=None, u1=None, u2=None, d1=None,
                 d2=None, tol=1.0e-7):
        super(PlanesAlongCurveByNumber, self).__init__()

        # Build
        pnt_builder = PointsAlongCurveByNumber(c, n, u1, u2, d1, d2, tol)
        if pnt_builder.npts == 0:
            msg = ('Failed to generate points along the curve for creating '
                   'planes along a curve by number.')
            raise RuntimeError(msg)
        npts = pnt_builder.npts
        pnts = pnt_builder.points
        prms = pnt_builder.parameters
        spacing = pnt_builder.spacing

        plns = []
        if isinstance(ref_pln, Plane):
            gp_pln = ref_pln.Pln()
            ax1 = gp_pln.Axis()
            dn = ax1.Direction()
            for p in pnts:
                pln = Plane(p, dn)
                plns.append(pln)
        else:
            for i in range(npts):
                p = pnts[i]
                u = prms[i]
                vn = c.deriv(u, 1)
                dn = Direction(vn)
                pln = Plane(p, dn)
                plns.append(pln)

        self._set_results('plns', plns)
        self._set_results('nplns', npts)
        self._set_results('prms', prms)
        self._set_results('ds', spacing)
        self._success = True

    @property
    def nplanes(self):
        """
        :return: The number of planes.
        :rtype: int
        """
        return self._get_results('nplns')

    @property
    def planes(self):
        """
        :return: The planes.
        :rtype: list[afem.geometry.entities.Plane]
        """
        return self._get_results('plns')

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


class PlanesAlongCurveByDistance(GeomBuilder):
    """
    Create planes along a curve by distance between points. The origin of the
    planes will be equidistant along the curve. This method calculates the
    number of points given the curve length and then uses
    :class:`.PlanesAlongCurveByNumber`.

    :param curve_like c: The curve.
    :param float maxd: The maximum allowed spacing between planes. The
        actual spacing will be adjusted to not to exceed this value.
    :param afem.geometry.entities.Plane ref_pln: The normal of this plane
        will be used to define the normal of all planes along the curve. If
        no plane is provided, then the first derivative of the curve will
        define the plane normal.
    :param float u1: The parameter of the first plane (default=c.u1).
    :param float u2: The parameter of the last plane (default=c.u2).
    :param float d1: An offset distance for the first plane. This is typically
        a positive number indicating a distance from *u1* towards *u2*.
    :param float d2: An offset distance for the last plane. This is typically
        a negative number indicating a distance from *u2* towards *u1*.
    :param int nmin: Minimum number of planes to create.
    :param float tol: Tolerance.

    :raise RuntimeError: If :class:`.PointsAlongCurveByDistance` fails to
        generate points along the curve.

    Usage:

    >>> from afem.geometry import LineByPoints, Point, PlanesAlongCurveByDistance
    >>> p1 = Point()
    >>> p2 = Point(10., 0., 0.)
    >>> line = LineByPoints(p1, p2).line
    >>> builder = PlanesAlongCurveByDistance(line, 5., u1=0., u2=10.)
    >>> builder.nplanes
    3
    >>> builder.parameters
    [0.0, 5.0, 10.0]
    """

    def __init__(self, c, maxd, ref_pln=None, u1=None, u2=None, d1=None,
                 d2=None, nmin=0, tol=1.0e-7):
        super(PlanesAlongCurveByDistance, self).__init__()

        # Build
        pnt_builder = PointsAlongCurveByDistance(c, maxd, u1, u2, d1, d2,
                                                 nmin, tol)
        if pnt_builder.npts == 0:
            msg = ('Failed to generate points along the curve for creating '
                   'planes along a curve by distance.')
            raise RuntimeError(msg)
        npts = pnt_builder.npts
        pnts = pnt_builder.points
        prms = pnt_builder.parameters
        spacing = pnt_builder.spacing

        plns = []
        if isinstance(ref_pln, Plane):
            gp_pln = ref_pln.Pln()
            ax1 = gp_pln.Axis()
            dn = ax1.Direction()
            for p in pnts:
                pln = Plane(p, dn)
                plns.append(pln)
        else:
            for i in range(npts):
                p = pnts[i]
                u = prms[i]
                vn = c.deriv(u, 1)
                dn = Direction(vn)
                pln = Plane(p, dn)
                plns.append(pln)

        self._set_results('plns', plns)
        self._set_results('nplns', npts)
        self._set_results('prms', prms)
        self._set_results('ds', spacing)
        self._success = True

    @property
    def nplanes(self):
        """
        :return: The number of planes.
        :rtype: int
        """
        return self._get_results('nplns')

    @property
    def planes(self):
        """
        :return: The planes.
        :rtype: list[afem.geometry.entities.Plane]
        """
        return self._get_results('plns')

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


class PlanesBetweenPlanesByNumber(GeomBuilder):
    """
    Create planes between two other planes. This method will create a line
    normal to the first plane at its origin and then intersect
    that line with the second plane. Planes are then created along this line
    using :class:`.PlanesAlongCurveByNumber`.

    :param afem.geometry.entities.Plane pln1: The first plane. This will be
        the reference plane to define the orientation of the new planes.
    :param afem.geometry.entities.Plane pln2: The last plane.
    :param int n: The number of planes.
    :param float d1: An offset distance for the first plane. This is typically
        a positive number indicating a distance from *u1* towards *u2*.
    :param float d2: An offset distance for the last plane. This is typically
        a negative number indicating a distance from *u2* towards *u1*.

    :raise RuntimeError: If the line extending from the normal of the first
        plane cannot be intersected with the second plane.

    Usage:

    >>> from afem.geometry import PlaneByNormal, PlanesBetweenPlanesByNumber
    >>> pln1 = PlaneByNormal((0., 0., 0.), (1., 0., 0.)).plane
    >>> pln2 = PlaneByNormal((10., 0., 0.), (1., 0., 0.)).plane
    >>> builder = PlanesBetweenPlanesByNumber(pln1, pln2, 3)
    >>> builder.nplanes
    3
    >>> builder.spacing
    5.0
    >>> pln1 = builder.planes[0]
    >>> pln1.eval(0., 0.)
    Point(0.0, 0.0, 0.0)
    >>> pln2 = builder.planes[1]
    >>> pln2.eval(0., 0.)
    Point(5.0, 0.0, 0.0)
    >>> pln3 = builder.planes[2]
    >>> pln3.eval(0., 0.)
    Point(10.0, 0.0, 0.0)
    """

    def __init__(self, pln1, pln2, n, d1=None, d2=None):
        super(PlanesBetweenPlanesByNumber, self).__init__()

        # Build
        p1 = pln1.eval(0., 0.)
        vn = pln1.norm(0., 0.)
        line = LineByVector(p1, vn).line
        csi = GeomAPI_IntCS(line.handle, pln2.handle)
        if csi.NbPoints() == 0:
            msg = ('Failed to intersect the second plane to create planes '
                   'between them.')
            raise RuntimeError(msg)
        p2 = csi.Point(1)

        line = NurbsCurveByPoints([p1, p2]).curve
        builder = PlanesAlongCurveByNumber(line, n, pln1, d1=d1, d2=d2)

        self._set_results('plns', builder.planes)
        self._set_results('nplns', builder.nplanes)
        self._set_results('ds', builder.spacing)
        self._success = True

    @property
    def nplanes(self):
        """
        :return: The number of planes.
        :rtype: int
        """
        return self._get_results('nplns')

    @property
    def planes(self):
        """
        :return: The planes.
        :rtype: list[afem.geometry.entities.Plane]
        """
        return self._get_results('plns')

    @property
    def spacing(self):
        """
        :return: The spacing between the first and second points if there
            are more than one point. Otherwise *None*.
        :rtype: float or None
        """
        return self._get_results('ds')


class PlanesBetweenPlanesByDistance(GeomBuilder):
    """
    Create planes between two other planes by distance. This method will
    create a line normal to the first plane at its origin and then intersect
    that line with the second plane. Planes are then created along this line
    using :class:`.PlanesAlongCurveByNumber`.

    :param afem.geometry.entities.Plane pln1: The first plane. This will be
        the reference plane to define the orientation of the new planes.
    :param afem.geometry.entities.Plane pln2: The last plane.
    :param float maxd: The maximum allowed spacing between planes. The
        actual spacing will be adjusted to not to exceed this value.
    :param float d1: An offset distance for the first plane. This is typically
        a positive number indicating a distance from *u1* towards *u2*.
    :param float d2: An offset distance for the last plane. This is typically
        a negative number indicating a distance from *u2* towards *u1*.
    :param int nmin: Minimum number of planes to create.

    :raise RuntimeError: If the line extending from the normal of the first
        plane cannot be intersected with the second plane.

    Usage:

    >>> from afem.geometry import PlaneByNormal, PlanesBetweenPlanesByDistance
    >>> pln1 = PlaneByNormal((0., 0., 0.), (1., 0., 0.)).plane
    >>> pln2 = PlaneByNormal((10., 0., 0.), (1., 0., 0.)).plane
    >>> builder = PlanesBetweenPlanesByDistance(pln1, pln2, 5.)
    >>> builder.nplanes
    3
    >>> builder.spacing
    5.0
    >>> pln1 = builder.planes[0]
    >>> pln1.eval(0., 0.)
    Point(0.0, 0.0, 0.0)
    >>> pln2 = builder.planes[1]
    >>> pln2.eval(0., 0.)
    Point(5.0, 0.0, 0.0)
    >>> pln3 = builder.planes[2]
    >>> pln3.eval(0., 0.)
    Point(10.0, 0.0, 0.0)
    """

    def __init__(self, pln1, pln2, maxd, d1=None, d2=None, nmin=0):
        super(PlanesBetweenPlanesByDistance, self).__init__()

        # Build
        p1 = pln1.eval(0., 0.)
        vn = pln1.norm(0., 0.)
        line = LineByVector(p1, vn).line
        csi = GeomAPI_IntCS(line.handle, pln2.handle)
        if csi.NbPoints() == 0:
            msg = ('Failed to intersect the second plane to create planes '
                   'between them.')
            raise RuntimeError(msg)
        p2 = csi.Point(1)

        line = NurbsCurveByPoints([p1, p2]).curve
        builder = PlanesAlongCurveByDistance(line, maxd, pln1, d1=d1, d2=d2,
                                             nmin=nmin)

        self._set_results('plns', builder.planes)
        self._set_results('nplns', builder.nplanes)
        self._set_results('ds', builder.spacing)
        self._success = True

    @property
    def nplanes(self):
        """
        :return: The number of planes.
        :rtype: int
        """
        return self._get_results('nplns')

    @property
    def planes(self):
        """
        :return: The planes.
        :rtype: list[afem.geometry.entities.Plane]
        """
        return self._get_results('plns')

    @property
    def spacing(self):
        """
        :return: The spacing between the first and second points if there
            are more than one point. Otherwise *None*.
        :rtype: float or None
        """
        return self._get_results('ds')


class NurbsSurfaceByData(GeomBuilder):
    """
    Create a NURBS surface by data.

    :param list[list[point_like]] cp: Two-dimensional list of control points.
    :param list[float] uknots: Knot vector for u-direction.
    :param list[float] vknots: Knot vector for v-direction.
    :param list[int] umult: Multiplicities of knot vector in u-direction.
    :param list[int] vmult: Multiplicities of knot vector in v-direction.
    :param int p: Degree in u-direction.
    :param int q: Degree in v-direction.
    :param list[list[float]] weights: Two-dimensional list of control
        point weights.
    :param bool is_u_periodic: Flag for periodicity in u-direction.
    :param bool is_v_periodic: Flag for periodicity in v-direction.

    Usage:

    >>> from afem.geometry import NurbsSurfaceByData, Point
    >>> cp = [[Point(), Point(10., 0., 0.)], [Point(0., 10., 0.), Point(10., 10., 0.)]]
    >>> uknots = [0., 1.]
    >>> vknots = [0., 1.]
    >>> umult = [2, 2]
    >>> vmult = [2, 2]
    >>> p = 1
    >>> q = 1
    >>> builder = NurbsSurfaceByData(cp, uknots, vknots, umult, vmult, p, q)
    >>> s = builder.surface
    >>> s.eval(0.5, 0.5)
    Point(5.0, 5.0, 0.0)
    """

    def __init__(self, cp, uknots, vknots, umult, vmult, p, q, weights=None,
                 is_u_periodic=False, is_v_periodic=False):
        super(NurbsSurfaceByData, self).__init__()

        # Build
        tcol_cp = to_tcolgp_array2_pnt(cp)
        tcol_uknots = to_tcolstd_array1_real(uknots)
        tcol_umult = to_tcolstd_array1_integer(umult)
        tcol_vknots = to_tcolstd_array1_real(vknots)
        tcol_vmult = to_tcolstd_array1_integer(vmult)
        p, q = int(p), int(q)
        if weights is None:
            weights = ones((tcol_cp.ColLength(), tcol_cp.RowLength()))
        tcol_weights = to_tcolstd_array2_real(weights)

        s = NurbsSurface(tcol_cp, tcol_uknots, tcol_vknots,
                         tcol_umult, tcol_vmult, p, q, is_u_periodic,
                         is_v_periodic)

        # Set the weights since using in construction causes an error.
        for i in range(1, tcol_weights.ColLength() + 1):
            for j in range(1, tcol_weights.RowLength() + 1):
                s.SetWeight(i, j, tcol_weights.Value(i, j))

        self._set_results('s', s)
        self._success = True

    @property
    def surface(self):
        """
        :return: The NURBS surface.
        :rtype: afem.geometry.entities.NurbsSurface
        """
        return self._get_results('s')


class NurbsSurfaceByInterp(GeomBuilder):
    """
    Create a surface by interpolating curves. This method was developed from
    scratch using "The NURBS Book" since OpenCASCADE did not support
    interpolating surfaces with *q* = 1.

    :param list[curve_like] crvs: List of curves to interpolate.
    :param int q: Degree. The parameter will be adjusted if the number of
        curves provided does not support the desired degree.
    :param str parm_type: Parametrization type ('uniform', 'chord',
        'centripetal').
    :param float tol2d: 2-D tolerance.

    Usage:

    >>> from afem.geometry import *
    >>> c1 = NurbsCurveByPoints([(0., 0., 0.), (10., 0., 0.)]).curve
    >>> c2 = NurbsCurveByPoints([(0., 5., 5.), (10., 5., 5.)]).curve
    >>> c3 = NurbsCurveByPoints([(0., 10., 0.), (10., 10., 0.)]).curve
    >>> builder = NurbsSurfaceByInterp([c1, c2, c3], 2)
    >>> assert builder.surface
    >>> s = builder.surface
    >>> s.eval(0.5, 0.5)
    Point(5.0, 5.0, 5.0)
    """

    def __init__(self, crvs, q=3, parm_type='chord', tol2d=1.0e-9):
        super(NurbsSurfaceByInterp, self).__init__()

        # Build
        q = int(q)

        ncrvs = len(crvs)
        if ncrvs - 1 < q:
            q = ncrvs - 1

        # Make curves compatible using section generator.
        sec_gen = GeomFill_SectionGenerator()
        for c in crvs:
            sec_gen.AddCurve(c.handle)
        sec_gen.Perform(tol2d)

        # Use old method since OCC struggles to interpolate curves of low
        # order. Gather all data for old method.
        p = sec_gen.Degree()
        is_u_periodic = sec_gen.IsPeriodic()
        tcol_uknots = TColStd_Array1OfReal(1, sec_gen.NbKnots())
        tcol_umult = TColStd_Array1OfInteger(1, sec_gen.NbKnots())
        sec_gen.KnotsAndMults(tcol_uknots, tcol_umult)
        temp = []
        for i in range(1, ncrvs + 1):
            tcol_poles = TColgp_Array1OfPnt(1, sec_gen.NbPoles())
            sec_gen.Poles(i, tcol_poles)
            tcol_weights = TColStd_Array1OfReal(1, sec_gen.NbPoles())
            sec_gen.Weights(i, tcol_weights)
            cp = to_np_from_tcolgp_array1_pnt(tcol_poles)
            w = to_np_from_tcolstd_array1_real(tcol_weights)
            cpw = homogenize_array1d(cp, w)
            temp.append(cpw)

        # Compute v-direction parameters between [0, 1].
        # Find parameters between each curve by averaging each segment.
        temp = array(temp, dtype=float)
        pnts_matrix = temp.transpose((1, 0, 2))
        n = sec_gen.NbPoles() - 1
        m = ncrvs - 1
        v_matrix = zeros((n + 1, m + 1), dtype=float)
        for i in range(0, n + 1):
            if parm_type.lower() in ['u', 'uniform']:
                vknots = uniform_parameters(pnts_matrix[i, :], 0., 1.)
            elif parm_type.lower() in ['ch', 'chord']:
                vknots = chord_parameters(pnts_matrix[i, :], 0., 1.)
            else:
                vknots = centripetal_parameters(pnts_matrix[i, :], 0., 1.)
            v_matrix[i] = vknots
        # Average each column.
        vknots = mean(v_matrix, axis=0, dtype=float)
        vknots[0] = 0.0
        vknots[-1] = 1.0

        s = m + q + 1
        vk = zeros(s + 1, dtype=float)
        vk[s - q:] = 1.0
        for j in range(1, m - q + 1):
            temp = 0.
            for i in range(j, j + q):
                temp += vknots[i]
            vk[j + q] = 1.0 / q * temp
        # Compute OCC vknots and vmult.
        tcol_vknot_seq = to_tcolstd_array1_real(vk)
        nv = CLib.bsplclib_KnotsLength(tcol_vknot_seq, False)
        tcol_vknots = TColStd_Array1OfReal(1, nv)
        tcol_vmult = TColStd_Array1OfInteger(1, nv)
        CLib.bsplclib_Knots(tcol_vknot_seq, tcol_vknots, tcol_vmult, False)

        # Perform n + 1 interpolations in v-direction to generate surface
        # control points.
        cpw = zeros((n + 1, m + 1, 4), dtype=float)
        for i in range(0, n + 1):
            qp = pnts_matrix[i]
            # Set up system of linear equations.
            a = zeros((m + 1, m + 1), dtype=float)
            for j in range(0, m + 1):
                # span1, _ = CLib.bsplclib_LocateParameter(pnt, tcol_vknots,
                #                                          tcol_vmult,
                #                                         vknots[j], False)
                # Evaluate basis function. Figure out how to use math_Matrix
                # and use OCC in future.
                # mat = occ_math.math_Matrix(1, q, 1, q)
                # CLib.bsplclib_EvalBsplineBasis(span, q, q, tcol_vknot_seq,
                #                                vknots[j], mat)
                span = find_span(m, q, vknots[j], vk)
                a[j, span - q: span + 1] = basis_funs(span, vknots[j], q, vk)
            # Solve for [a][cp] = [qp] using LU decomposition.
            lu, piv = lu_factor(a, overwrite_a=True, check_finite=False)
            cpw[i, :] = lu_solve((lu, piv), qp, trans=0, overwrite_b=True,
                                 check_finite=True)

        # Create surface.
        cp, w = dehomogenize_array2d(cpw)
        tcol_poles = to_tcolgp_array2_pnt(cp)
        tcol_weights = to_tcolstd_array2_real(w)
        s = NurbsSurface(tcol_poles, tcol_weights, tcol_uknots, tcol_vknots,
                         tcol_umult, tcol_vmult, p, q, is_u_periodic, False)

        self._set_results('s', s)
        self._success = True

    @property
    def surface(self):
        """
        :return: The NURBS surface.
        :rtype: afem.geometry.entities.NurbsSurface
        """
        return self._get_results('s')


class NurbsSurfaceByApprox(GeomBuilder):
    """
    Create a NURBS surface by approximating curves.

    :param list[curve_like] crvs: List of curves.
    :param int dmin: Minimum degree.
    :param int dmax: Maximum degree.
    :param float tol3d: 3-D tolerance.
    :param float tol2d: 2-D tolerance.
    :param int niter: Number of iterations.
    :param str continuity: Desired continuity of curve ('C0', 'G1', 'C1',
        'G2', 'C2', 'C3').
    :param str parm_type: Parametrization type ('uniform', 'chord',
        'centripetal').

    :raise RuntimeError: If OCC method fails to approximate the curves with a
        surface.

    For more information see GeomFill_AppSurf_.

    .. _GeomFill_AppSurf: https://www.opencascade.com/doc/occt-7.1.0/refman/html/class_geom_fill___app_surf.html

    Usage:

    >>> from afem.geometry import *
    >>> c1 = NurbsCurveByPoints([(0., 0., 0.), (10., 0., 0.)]).curve
    >>> c2 = NurbsCurveByPoints([(0., 5., 5.), (10., 5., 5.)]).curve
    >>> c3 = NurbsCurveByPoints([(0., 10., 0.), (10., 10., 0.)]).curve
    >>> builder = NurbsSurfaceByApprox([c1, c2, c3])
    >>> assert builder.surface
    >>> builder.told3d_reached
    0.0
    >>> builder.told2d_reached
    0.0
    >>> s = builder.surface
    >>> s.eval(0.5, 0.5)
    Point(5.0, 5.0, 5.0)
    """

    def __init__(self, crvs, dmin=3, dmax=8, tol3d=1.0e-3, tol2d=1.0e-6,
                 niter=5, continuity='C2', parm_type='chord'):
        super(NurbsSurfaceByApprox, self).__init__()

        # Build
        dmin = int(dmin)
        dmax = int(dmax)

        # Initialize approximation tool
        app_tool = GeomFill_AppSurf(dmin, dmax, tol3d, tol2d, niter)

        # Set parametrization type
        try:
            parm_type = _occ_parm_type[parm_type.lower()]
        except (KeyError, AttributeError):
            parm_type = Approx_ChordLength
        app_tool.SetParType(parm_type)

        # Set continuity
        try:
            cont = _occ_continuity[continuity.upper()]
        except (KeyError, AttributeError):
            cont = GeomAbs_C2
        app_tool.SetContinuity(cont)

        # Use section generator to make all curves compatible
        sec_gen = GeomFill_SectionGenerator()
        for c in crvs:
            sec_gen.AddCurve(c.handle)
        sec_gen.Perform(tol3d)

        # Create line tool
        line_tool = GeomFill_Line(len(crvs))

        # Perform the approximation
        app_tool.Perform(line_tool.GetHandle(), sec_gen, False)
        if not app_tool.IsDone():
            msg = 'OCC method failed.'
            raise RuntimeError(msg)

        # Create a surface
        tcol_poles = app_tool.SurfPoles()
        tcol_weights = app_tool.SurfWeights()
        tcol_uknots = app_tool.SurfUKnots()
        tcol_vknots = app_tool.SurfVKnots()
        tcol_umult = app_tool.SurfUMults()
        tcol_vmult = app_tool.SurfVMults()
        p = app_tool.UDegree()
        q = app_tool.VDegree()
        is_u_periodic = sec_gen.IsPeriodic()
        is_v_periodic = False
        s = NurbsSurface(tcol_poles, tcol_weights, tcol_uknots, tcol_vknots,
                         tcol_umult, tcol_vmult, p, q, is_u_periodic,
                         is_v_periodic)

        tol3d_reached, tol2d_reached = app_tool.TolReached()
        self._set_results('s', s)
        self._set_results('tol3d_reached', tol3d_reached)
        self._set_results('tol2d_reached', tol2d_reached)
        self._success = True

    @property
    def surface(self):
        """
        :return: The NURBS surface.
        :rtype: afem.geometry.entities.NurbsSurface
        """
        return self._get_results('s')

    @property
    def told3d_reached(self):
        """
        :return: 3-D tolerance achieved.
        :rtype: float
        """
        return self._get_results('tol3d_reached')

    @property
    def told2d_reached(self):
        """
        :return: 2-D tolerance achieved.
        :rtype: float
        """
        return self._get_results('tol2d_reached')


if __name__ == "__main__":
    import doctest

    doctest.testmod()
