from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeEdge
from OCC.GeomAPI import GeomAPI_IntCS
from OCC.GeomAbs import GeomAbs_BSplineCurve, GeomAbs_BezierCurve, GeomAbs_Line
from OCC.GeomAdaptor import GeomAdaptor_Curve
from OCC.GeomInt import GeomInt_IntSS
from OCC.IntTools import IntTools_EdgeEdge
from OCC.ShapeFix import ShapeFix_ShapeTolerance
from OCC.TopAbs import TopAbs_VERTEX
from numpy import mean

from afem.geometry.entities import Line, Point
from .create import create_nurbs_curve_from_occ
from ...config import Settings

__all__ = []


def intersect_curve_curve(curve1, curve2, itol=None):
    """
    Find the intersection points of two curves.
    """
    if itol is None:
        itol = Settings.gtol

    # Build edges from curve.
    e1 = BRepBuilderAPI_MakeEdge(curve1.handle).Edge()
    e2 = BRepBuilderAPI_MakeEdge(curve2.handle).Edge()

    # Set tolerance to be half intersection tolerance.
    shp_tol = ShapeFix_ShapeTolerance()
    tol = itol / 2.
    shp_tol.SetTolerance(e1, tol)
    shp_tol.SetTolerance(e2, tol)

    # Perform edge-edge intersection
    cci = IntTools_EdgeEdge(e1, e2)
    cci.Perform()
    if not cci.IsDone():
        return 0, []

    # Gather results of point intersection only.
    results = []
    common_parts = cci.CommonParts()
    for i in range(1, common_parts.Length() + 1):
        common_part = common_parts.Value(i)
        if not common_part.Type() == TopAbs_VERTEX:
            continue
        u1 = common_part.VertexParameter1()
        u2 = common_part.VertexParameter2()
        p1 = curve1.eval(u1)
        p2 = curve2.eval(u2)
        pi = mean([p1, p2], axis=0)
        pi = Point(*pi)
        results.append([(u1, u2), pi])

    npts = len(results)
    return npts, results


def intersect_curve_surface(curve, surface):
    """
    Find the intersection points between a curve and a surface.
    """

    # OCC intersection.
    csi = GeomAPI_IntCS(curve.handle, surface.handle)
    if not csi.IsDone():
        return 0, []

    # Gather results of point intersections.
    results = []
    for i in range(1, csi.NbPoints() + 1):
        u, v, t = csi.Parameters(i)
        pc = curve.eval(t)
        ps = surface.eval(u, v)
        pi = mean([pc, ps], axis=0)
        pi = Point(*pi)
        results.append([(t, u, v), pi])

    npts = len(results)
    return npts, results


def intersect_surface_surface(surface1, surface2, itol=None):
    """
    Intersect two surfaces.
    """
    if itol is None:
        itol = Settings.gtol

    # OCC intersect.
    ssi = GeomInt_IntSS(surface1.handle, surface2.handle, itol,
                        True, False, False)
    if not ssi.IsDone():
        return []

    # Build curves.
    ncrvs = ssi.NbLines()
    crvs = []
    for i in range(1, ncrvs + 1):
        hcrv = ssi.Line(i)
        adp_crv = GeomAdaptor_Curve(hcrv)
        if adp_crv.GetType() == GeomAbs_Line:
            gp_lin = adp_crv.Line()
            crv = Line(gp_lin)
            crvs.append(crv)
        elif adp_crv.GetType() in [GeomAbs_BezierCurve, GeomAbs_BSplineCurve]:
            crv = adp_crv.BSpline().GetObject()
            crv = create_nurbs_curve_from_occ(crv)
            crvs.append(crv)

    return crvs
