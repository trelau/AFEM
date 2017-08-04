from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeEdge
from OCC.GeomAPI import GeomAPI_IntCS
from OCC.GeomAbs import GeomAbs_BSplineCurve, GeomAbs_BezierCurve, GeomAbs_Line
from OCC.GeomAdaptor import GeomAdaptor_Curve
from OCC.GeomInt import GeomInt_IntSS
from OCC.IntTools import IntTools_EdgeEdge
from OCC.ShapeFix import ShapeFix_ShapeTolerance
from OCC.TopAbs import TopAbs_VERTEX
from numpy import float64, inf, mean, zeros
from scipy.spatial import KDTree

from afem.geometry.check import CheckGeom
from afem.geometry.create import create_nurbs_curve_from_occ
from afem.geometry.entities import Line, Point
from afem.geometry.methods.distance import curve_nearest_point

__all__ = ["IntersectGeom", "CurveIntersector", "IntersectCurveCurve",
           "IntersectCurveSurface", "SurfaceIntersector",
           "IntersectSurfaceSurface", "IntersectError"]


class IntersectGeom(object):
    """
    Intersect geometry.
    """

    @staticmethod
    def perform(geom1, geom2, itol=None):
        """
        Perform the intersection of two geometries.

        :param geom1: Geometry entity 1.
        :param geom2: Geometry entity 2.
        :param float itol: Intersection tolerance.

        :return: Intersector object depending on type of *geom1* and *geom2*.
            *None* if returned if entities are not supported.
        """
        # Curve intersection with...
        if CheckGeom.is_curve_like(geom1):
            # Curve
            if CheckGeom.is_curve_like(geom2):
                return IntersectCurveCurve(geom1, geom2, itol)
            # Surface/plane
            if CheckGeom.is_surface_like(geom2):
                return IntersectCurveSurface(geom1, geom2)

        # Surface intersection with...
        if CheckGeom.is_surface_like(geom1):
            # Curve
            if CheckGeom.is_curve_like(geom2):
                return IntersectCurveSurface(geom2, geom1)
            # Surface/plane
            if CheckGeom.is_surface_like(geom2):
                return IntersectSurfaceSurface(geom1, geom2, itol)

        # Return error if combination not supported.
        return IntersectError()


class IntersectError(object):
    """
    Class for handling intersection errors.
    """

    @property
    def success(self):
        return False

    @property
    def npts(self):
        return 0

    @property
    def ncrvs(self):
        return 0

    @property
    def points(self):
        return []

    @property
    def icrvs(self):
        return []


class CurveIntersector(object):
    """
    Base class for handling curve intersection methods and results.
    """

    def __init__(self, c1, c2):
        self._c1 = c1
        self._c2 = c2
        self._npts = 0
        self._results = []
        self._kdt = None

    def _set_results(self, npts, results):
        """
        Set curve intersection results.
        """
        if npts > 0:
            self._npts = npts
            # Replace ndarrays with point instances and build kd-tree.
            data = zeros((npts, 3), dtype=float64)
            for i in range(npts):
                data[i] = results[i][1]
            self._results = results
            self._kdt = KDTree(data)

    @property
    def npts(self):
        """
        :return: Number of intersection points.
        :rtype: int
        """
        return self._npts

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
        :return: List of intersection points.
        :rtype: list[afem.geometry.entities.Point]
        """
        if self._npts <= 0:
            return []
        return [results[1] for results in self._results]

    @property
    def parameters(self):
        """
        :return: List of intersection parameters. For a curve-curve
            intersection this will be a list of tuples containing the
            parameters of each curve [(u1, u2), (u1, u2), ...]. For a
            curve-surface intersection this will be a list of tuples
            containing the parameters for the surface and then the curve
            [(u, v, t), (u, v, t), ...].
        :rtype: list[tuple(float)]
        """
        if self._npts <= 0:
            return []
        return [results[0] for results in self._results]

    def point(self, indx=1):
        """
        Return the point result by index.

        :param int indx: Index for point.

        :return: Intersection point.
        :rtype: afem.geometry.entities.Point
        """
        return self._results[indx - 1][1]

    def query_point(self, p0, distance_upper_bound=inf):
        """
        Find the intersection result nearest to the provided point.

        :param point_like p0: Point to search from.
        :param float distance_upper_bound: Return only results within this
            distance.

        :return: Distance to nearest intersection result and its index (d, i).
            Returns (None, None) if no results are available.
        :rtype: tuple
        """
        if not self.success:
            return None, None
        p0 = CheckGeom.to_point(p0)
        d, i = self._kdt.query(p0.xyz, 1, 0., 2, distance_upper_bound)
        return d, i


class IntersectCurveCurve(CurveIntersector):
    """
    Curve-curve intersection. This method converts the curves to edges and
    performs the intersections that way. This proved to be more robust than
    OpenCASCADE's native curve-curve intersection tool.

    :param curve_like crv1: The first curve.
    :param curve_like crv2: The second curve.
    :param float itol: The intersection tolerance.

    For more information see IntTools_EdgeEdge_.

    .. _IntTools_EdgeEdge: https://www.opencascade.com/doc/occt-7.1.0/refman/html/class_int_tools___edge_edge.html

    Usage:

    >>> from afem.geometry import *
    >>> c1 = NurbsCurveByPoints([(0., 0., 0.), (10., 0., 0.)]).curve
    >>> c2 = NurbsCurveByPoints([(5., 0., 0.), (5., 5., 0.)]).curve
    >>> cci = IntersectCurveCurve(c1, c2)
    >>> assert cci.success
    >>> cci.npts
    1
    >>> cci.point(1)
    Point(4.999999999999886, 2.5000361106880673e-08, 0.0)
    """

    def __init__(self, crv1, crv2, itol=1.0e-7):
        super(IntersectCurveCurve, self).__init__(crv1, crv2)

        # Perform
        # TODO Reevaluate if IntTools_EdgeEdge is needed.

        # Build edges from curve.
        e1 = BRepBuilderAPI_MakeEdge(crv1.handle).Edge()
        e2 = BRepBuilderAPI_MakeEdge(crv2.handle).Edge()

        # Set tolerance to be half intersection tolerance.
        shp_tol = ShapeFix_ShapeTolerance()
        tol = itol / 2.
        shp_tol.SetTolerance(e1, tol)
        shp_tol.SetTolerance(e2, tol)

        # Perform edge-edge intersection
        cci = IntTools_EdgeEdge(e1, e2)
        cci.Perform()

        # Gather results of point intersection only.
        results = []
        common_parts = cci.CommonParts()
        for i in range(1, common_parts.Length() + 1):
            common_part = common_parts.Value(i)
            if not common_part.Type() == TopAbs_VERTEX:
                continue
            u1 = common_part.VertexParameter1()
            u2 = common_part.VertexParameter2()
            p1 = crv1.eval(u1)
            p2 = crv2.eval(u2)
            pi = mean([p1, p2], axis=0)
            pi = Point(*pi)
            results.append([(u1, u2), pi])

        npts = len(results)
        self._set_results(npts, results)


class IntersectCurveSurface(CurveIntersector):
    """
    Curve-surface intersection.

    :param curve_like crv: The curve.
    :param surface_like srf: The surface.

    For more information see GeomAPI_IntCS_.

    .. _GeomAPI_IntCS: https://www.opencascade.com/doc/occt-7.1.0/refman/html/class_geom_a_p_i___int_c_s.html

    Usage:

    >>> from afem.geometry import *
    >>> c = NurbsCurveByPoints([(5., 5., 10.), (5., 5., -10.)]).curve
    >>> c1 = NurbsCurveByPoints([(0., 0., 0.), (10., 0., 0.)]).curve
    >>> c2 = NurbsCurveByPoints([(0., 5., 5.), (10., 5., 5.)]).curve
    >>> c3 = NurbsCurveByPoints([(0., 10., 0.), (10., 10., 0.)]).curve
    >>> s = NurbsSurfaceByApprox([c1, c2, c3]).surface
    >>> csi = IntersectCurveSurface(c, s)
    >>> assert csi.success
    >>> csi.npts
    1
    >>> csi.point(1)
    Point(5.0, 5.0, 5.0)
    """

    def __init__(self, crv, srf):
        super(IntersectCurveSurface, self).__init__(crv, srf)

        # Perform

        # OCC intersection.
        csi = GeomAPI_IntCS(crv.handle, srf.handle)
        results = []
        for i in range(1, csi.NbPoints() + 1):
            u, v, t = csi.Parameters(i)
            pc = crv.eval(t)
            ps = srf.eval(u, v)
            pi = mean([pc, ps], axis=0)
            pi = Point(*pi)
            results.append([(t, u, v), pi])

        npts = len(results)
        self._set_results(npts, results)

    @property
    def curve_parameters(self):
        return [results[0][0] for results in self._results]

    @property
    def surface_parameters(self):
        return [results[0][1:] for results in self._results]


class SurfaceIntersector(object):
    """
    Base class for handling surface intersection methods and results.
    """

    def __init__(self):
        self._crvs = []

    @property
    def ncrvs(self):
        """
        :return: Number of intersection curves.
        :rtype: int
        """
        return len(self._crvs)

    @property
    def success(self):
        """
        :return: *True* if successful, *False* if not.
        :rtype: bool
        """
        return self.ncrvs > 0

    @property
    def curves(self):
        """
        :return: The intersection curves.
        :rtype: list[afem.geometry.entities.Curve]
        """
        return self._crvs

    def curve(self, indx=1):
        """
        Generate an intersection curve.

        :param int indx: Index of intersection curve.

        :return: Intersection curve.
        :rtype: afem.geometry.entities.Curve
        """
        return self._crvs[indx - 1]

    def curve_nearest_point(self, pnt):
        """
        Get the index of the intersection curve that is nearest to the given
        reference point.

        :param point_like pnt: Reference point.

        :return: Index of curve nearest point.
        :rtype: int
        """
        return curve_nearest_point(pnt, self._crvs)


class IntersectSurfaceSurface(SurfaceIntersector):
    """
    Surface-surface intersection.

    :param surface_like srf1: The first surface.
    :param surface_like srf2: The second surface.
    :param float itol: Intersection tolerance.

    For more information see GeomAPI_IntSS_.

    .. _GeomAPI_IntSS: https://www.opencascade.com/doc/occt-7.1.0/refman/html/class_geom_a_p_i___int_s_s.html

    Usage:

    >>> from afem.geometry import *
    >>> c1 = NurbsCurveByPoints([(0., 0., 0.), (10., 0., 0.)]).curve
    >>> c2 = NurbsCurveByPoints([(0., 5., 5.), (10., 5., 5.)]).curve
    >>> c3 = NurbsCurveByPoints([(0., 10., 0.), (10., 10., 0.)]).curve
    >>> s = NurbsSurfaceByApprox([c1, c2, c3]).surface
    >>> pln = PlaneByNormal((5., 5., 0.), (1., 0., 0.)).plane
    >>> ssi = IntersectSurfaceSurface(s, pln)
    >>> assert ssi.success
    >>> ssi.ncrvs
    1
    >>> c = ssi.curve(1)
    >>> c.eval(0.5)
    Point(5.0, 4.999983574755282, 4.99999990730418)
    """

    def __init__(self, srf1, srf2, itol=1.0e-7):
        super(IntersectSurfaceSurface, self).__init__()

        # Perform

        # OCC intersect
        ssi = GeomInt_IntSS(srf1.handle, srf2.handle, itol,
                            True, False, False)

        # Build curves
        ncrvs = ssi.NbLines()
        crvs = []
        for i in range(1, ncrvs + 1):
            hcrv = ssi.Line(i)
            adp_crv = GeomAdaptor_Curve(hcrv)
            if adp_crv.GetType() == GeomAbs_Line:
                gp_lin = adp_crv.Line()
                crv = Line(gp_lin)
                crvs.append(crv)
            elif adp_crv.GetType() in [GeomAbs_BezierCurve,
                                       GeomAbs_BSplineCurve]:
                crv = adp_crv.BSpline().GetObject()
                crv = create_nurbs_curve_from_occ(crv)
                crvs.append(crv)
            else:
                msg = 'Unsupported line created in IntersectSurfaceSurface.'
                raise RuntimeError(msg)

        self._crvs = crvs


if __name__ == "__main__":
    import doctest

    doctest.testmod()
