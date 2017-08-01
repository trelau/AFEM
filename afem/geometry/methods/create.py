import OCC.BSplCLib as CLib
from OCC.Approx import Approx_Centripetal, Approx_ChordLength, \
    Approx_IsoParametric
from OCC.GCPnts import GCPnts_AbscissaPoint, GCPnts_UniformAbscissa
from OCC.GeomAPI import GeomAPI_Interpolate, GeomAPI_PointsToBSpline, \
    GeomAPI_ProjectPointOnCurve, GeomAPI_ProjectPointOnSurf
from OCC.GeomAbs import GeomAbs_C0, GeomAbs_C1, GeomAbs_C2, GeomAbs_C3, \
    GeomAbs_G1, GeomAbs_G2
from OCC.GeomAdaptor import GeomAdaptor_Curve
from OCC.GeomFill import GeomFill_AppSurf, GeomFill_Line, \
    GeomFill_SectionGenerator
from OCC.TColStd import TColStd_Array1OfInteger, TColStd_Array1OfReal, \
    TColStd_Array2OfReal
from OCC.TColgp import TColgp_Array1OfPnt, TColgp_Array1OfPnt2d, \
    TColgp_Array2OfPnt
from OCC.gp import gp_Dir, gp_Pln, gp_Pnt
from numpy import array, cross, mean, ones, zeros
from numpy.linalg import norm
from scipy.linalg import lu_factor, lu_solve

from ..utils import dehomogenize_array2d, homogenize_array1d
from ...occ.utils import to_np_from_tcolgp_array1_pnt, \
    to_np_from_tcolstd_array1_real, to_tcolgp_array1_pnt, \
    to_tcolgp_array1_pnt2d, to_tcolgp_array2_pnt, to_tcolgp_harray1_pnt, \
    to_tcolstd_array1_integer, to_tcolstd_array1_real, \
    to_tcolstd_array2_real
from .evaluate import basis_funs, find_span
from .parameterize import centripetal, chord_length, uniform
from afem.geometry.entities import *
from ...config import Settings

__all__ = []

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


def create_null_line():
    """
    Create Line for downcasting.
    """
    origin = gp_Pnt()
    d = gp_Dir(1., 0., 0.)
    return Line(origin, d)


def create_line_from_occ(lin):
    """
    Create Line from OCC data.
    """
    return Line(lin.Lin())


def create_nurbs_curve(cp, knots, mult, p, weights=None, is_periodic=False):
    """
    Create a NURBS curve from data.
    """
    # Convert data.
    tcol_cp = to_tcolgp_array1_pnt(cp)
    tcol_knots = to_tcolstd_array1_real(knots)
    tcol_mult = to_tcolstd_array1_integer(mult)
    p = int(p)
    if weights is None:
        weights = [1.] * tcol_cp.Length()
    tcol_weights = to_tcolstd_array1_real(weights)

    # Create the curve.
    c = NurbsCurve(tcol_cp, tcol_weights, tcol_knots, tcol_mult, p,
                   is_periodic)
    return c


def create_nurbs_curve2d(cp, knots, mult, p, weights=None, is_periodic=False):
    """
    Create a 2-D NURBS curve from data.
    """
    # Convert data.
    tcol_cp = to_tcolgp_array1_pnt2d(cp)
    tcol_knots = to_tcolstd_array1_real(knots)
    tcol_mult = to_tcolstd_array1_integer(mult)
    p = int(p)
    if weights is None:
        weights = [1.] * tcol_cp.Length()
    tcol_weights = to_tcolstd_array1_real(weights)

    # Create the curve.
    c = NurbsCurve2D(tcol_cp, tcol_weights, tcol_knots, tcol_mult, p,
                     is_periodic)
    return c


def create_nurbs_curve_from_occ(crv):
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

    c = NurbsCurve(tcol_poles, tcol_weights, tcol_knots, tcol_mult, p,
                   is_periodic)
    return c


def create_nurbs_curve_from_occ2d(crv2d):
    """
    Create a 2-D NURBS curve from an OCC curve.
    """
    # Gather OCC data.
    tcol_poles = TColgp_Array1OfPnt2d(1, crv2d.NbPoles())
    crv2d.Poles(tcol_poles)
    tcol_weights = TColStd_Array1OfReal(1, crv2d.NbPoles())
    crv2d.Weights(tcol_weights)
    tcol_knots = TColStd_Array1OfReal(1, crv2d.NbKnots())
    crv2d.Knots(tcol_knots)
    tcol_mult = TColStd_Array1OfInteger(1, crv2d.NbKnots())
    crv2d.Multiplicities(tcol_mult)
    p = crv2d.Degree()
    is_periodic = crv2d.IsPeriodic()

    c = NurbsCurve2D(tcol_poles, tcol_weights, tcol_knots, tcol_mult, p,
                     is_periodic)
    return c


def create_nurbs_surface(cp, uknots, vknots, umult, vmult, p, q, weights=None,
                         is_u_periodic=False, is_v_periodic=False):
    """
    Create a NURBS surface from data.
    """
    # Convert data.
    tcol_cp = to_tcolgp_array2_pnt(cp)
    tcol_uknots = to_tcolstd_array1_real(uknots)
    tcol_umult = to_tcolstd_array1_integer(umult)
    tcol_vknots = to_tcolstd_array1_real(vknots)
    tcol_vmult = to_tcolstd_array1_integer(vmult)
    p, q = int(p), int(q)
    if weights is None:
        weights = ones((tcol_cp.ColLength(), tcol_cp.RowLength()))
    tcol_weights = to_tcolstd_array2_real(weights)

    # Create the surface.
    s = NurbsSurface(tcol_cp, tcol_uknots, tcol_vknots,
                     tcol_umult, tcol_vmult, p, q, is_u_periodic,
                     is_v_periodic)
    # Set the weights since using in construction causes an error.
    for i in range(1, tcol_weights.ColLength() + 1):
        for j in range(1, tcol_weights.RowLength() + 1):
            s.SetWeight(i, j, tcol_weights.Value(i, j))

    return s


def create_nurbs_surface_from_occ(srf):
    """
    Create a NURBS surface from an OCC surface.
    """
    # Gather OCC data.
    tcol_poles = TColgp_Array2OfPnt(1, srf.NbUPoles(), 1, srf.NbVPoles())
    srf.Poles(tcol_poles)
    tcol_weights = TColStd_Array2OfReal(1, srf.NbUPoles(), 1, srf.NbVPoles())
    srf.Weights(tcol_weights)
    tcol_uknots = TColStd_Array1OfReal(1, srf.NbUKnots())
    srf.UKnots(tcol_uknots)
    tcol_vknots = TColStd_Array1OfReal(1, srf.NbVKnots())
    srf.VKnots(tcol_vknots)
    tcol_umult = TColStd_Array1OfInteger(1, srf.NbUKnots())
    srf.UMultiplicities(tcol_umult)
    tcol_vmult = TColStd_Array1OfInteger(1, srf.NbVKnots())
    srf.VMultiplicities(tcol_vmult)
    p = srf.UDegree()
    q = srf.VDegree()
    is_u_periodic = srf.IsUPeriodic()
    is_v_periodic = srf.IsVPeriodic()

    s = NurbsSurface(tcol_poles, tcol_weights, tcol_uknots, tcol_vknots,
                     tcol_umult, tcol_vmult, p, q, is_u_periodic,
                     is_v_periodic)
    return s


def create_point_from_other(curve, dx, u0):
    """
    Create a point on the curve at a specified distance from a parameter.
    """
    pac = GCPnts_AbscissaPoint(Settings.gtol, curve.adaptor, dx, u0)
    if not pac.IsDone():
        return None, None

    u = pac.Parameter()
    return u, curve.eval(u)


def create_points_along_curve(curve, maxd, npts, u1, u2, s1=None, s2=None):
    """
    Create equally spaced points along the curve.
    """
    # Check parameters.
    if u1 is None:
        u1 = curve.u1
    if u2 is None:
        u2 = curve.u2
    if u1 > u2:
        u1, u2 = u2, u1

    # Adjust u1 and u2 if s1 and s2 are not None.
    if s1 is not None:
        u1, _ = create_point_from_other(curve, s1, u1)
    if s2 is not None:
        u2, _ = create_point_from_other(curve, s2, u2)

    if maxd is not None:
        # Adjust step size if necessary.
        arc_len = curve.arc_length(u1, u2)
        nb_pts = int(arc_len / maxd) + 1
    else:
        nb_pts = int(npts)

    # Minimum number of points if maxd and npts are provided.
    if maxd is not None and npts is not None:
        if nb_pts < npts:
            nb_pts = int(npts)

    # OCC uniform abscissa.
    occ_pnts = GCPnts_UniformAbscissa(curve.adaptor, nb_pts, u1, u2,
                                      Settings.gtol)
    if not occ_pnts.IsDone():
        return 0, [], []

    # Gather results.
    npts = occ_pnts.NbPoints()
    points = []
    params = []
    for i in range(1, npts + 1):
        u = occ_pnts.Parameter(i)
        p = curve.eval(u)
        params.append(u)
        points.append(p)

    return npts, points, params


def create_crv_by_interp_pnts(qp, is_periodic, tol):
    """
    Interpolate the points with degree 3 curve.
    """
    # Convert points to OCC harray.
    tcol_hpnts = to_tcolgp_harray1_pnt(qp)
    if not tcol_hpnts:
        return None

    if tol is None:
        tol = Settings.gtol

    # Perform interpolation.
    interp = GeomAPI_Interpolate(tcol_hpnts.GetHandle(), is_periodic, tol)
    interp.Perform()
    if not interp.IsDone():
        return None

    # Create curve.
    occ_crv = interp.Curve().GetObject()

    # Convert to AFEM.
    return create_nurbs_curve_from_occ(occ_crv)


def create_crv_by_approx_pnts(qp, dmin, dmax, continuity, tol):
    """
    Fit a NURBS curve to an array of points.
    """
    # Convert points to OCC array.
    tcol_pnts = to_tcolgp_array1_pnt(qp)
    if not tcol_pnts:
        return None

    # Get OCC continuity.
    try:
        cont = _occ_continuity[continuity.upper()]
    except (KeyError, AttributeError):
        cont = GeomAbs_C2

    # Fit the points.
    fit = GeomAPI_PointsToBSpline(tcol_pnts, dmin, dmax, cont, tol)
    if not fit.IsDone():
        return None

    # Get the curve.
    occ_crv = fit.Curve().GetObject()

    # Convert to AFEM.
    return create_nurbs_curve_from_occ(occ_crv)


def create_isocurve(surface, u=None, v=None):
    """
    Create an isocurve from the surface.
    """
    if u is not None:
        hcrv = surface.UIso(u)
        adp_crv = GeomAdaptor_Curve(hcrv)
        occ_crv = adp_crv.BSpline().GetObject()
        return create_nurbs_curve_from_occ(occ_crv)
    if v is not None:
        hcrv = surface.VIso(v)
        adp_crv = GeomAdaptor_Curve(hcrv)
        occ_crv = adp_crv.BSpline().GetObject()
        return create_nurbs_curve_from_occ(occ_crv)
    return None


def create_plane_by_axes(origin, axes):
    """
    Create a plane define by an origin and axes.
    """
    if not isinstance(axes, str):
        return None
    if axes.lower() not in ['xy', 'yx', 'xz', 'zx', 'yz', 'zy']:
        return None

    if axes.lower() in ['xy', 'yx']:
        n = Direction(0., 0., 1.)
        vx = Direction(1., 0., 0.)
        ax = Axis3(origin, n, vx)
        pln = gp_Pln(ax)
        return Plane(pln)

    if axes.lower() in ['yz', 'zy']:
        n = Direction(1., 0., 0.)
        vx = Direction(0., 1., 0.)
        ax = Axis3(origin, n, vx)
        pln = gp_Pln(ax)
        return Plane(pln)

    if axes.lower() in ['xz', 'zx']:
        n = Direction(0., 1., 0.)
        vx = Direction(1., 0., 0.)
        ax = Axis3(origin, n, vx)
        pln = gp_Pln(ax)
        return Plane(pln)

    return None


def create_plane_by_points(p1, p2, p3):
    """
    Create a plane by points.
    """
    vx = p2.xyz - p1.xyz
    v31 = p3.xyz - p1.xyz
    vn = cross(vx, v31)
    if norm(vn) <= 1.0e-12:
        return None

    n = Direction(*vn)
    vx = Direction(*vx)
    ax = Axis3(p1, n, vx)
    pln = gp_Pln(ax)
    return Plane(pln)


def create_plane_on_curve(curve, u=None, dx=None, pnt=None, pref=None):
    """
    Create a plane at a point on the curve.
    """
    if u is None:
        u = curve.u1

    if isinstance(pnt, Point):
        proj = GeomAPI_ProjectPointOnCurve(pnt, curve.handle)
        if proj.NbPoints() < 0:
            u = proj.LowerDistanceParameter()

    if dx is not None:
        u, p0 = create_point_from_other(curve, dx, u)
    else:
        p0 = curve.eval(u)

    if isinstance(pref, Plane):
        gp_pln = pref.Pln()
        ax1 = gp_pln.Axis()
        dn = ax1.Direction()
        return Plane(p0, dn)

    vn = curve.deriv(u, 1)
    dn = Direction(vn)
    return Plane(p0, dn)


def create_planes_along_curve(curve, maxd, npts, pref, u1, u2,
                              s1=None, s2=None):
    """
    Create planes along a curve.
    """
    # Create points along the curve.
    npts, pnts, prms = create_points_along_curve(curve, maxd, npts, u1, u2,
                                                 s1, s2)
    if npts == 0:
        return []

    # Create a plane at each point using the reference plane or the curve
    # tangent.
    planes = []
    if isinstance(pref, Plane):
        gp_pln = pref.Pln()
        ax1 = gp_pln.Axis()
        dn = ax1.Direction()
        for p in pnts:
            pln = Plane(p, dn)
            planes.append(pln)
    else:
        for i in range(npts):
            p = pnts[i]
            u = prms[i]
            vn = curve.deriv(u, 1)
            dn = Direction(vn)
            pln = Plane(p, dn)
            planes.append(pln)

    return planes


def create_planes_between_planes(pln1, pln2, maxd, nplns, s1=None, s2=None):
    """
    Create planes between two other planes.
    """
    p1 = pln1.eval()
    proj = GeomAPI_ProjectPointOnSurf(p1, pln2.handle)
    u, v = proj.Parameters(1)
    p2 = pln2.eval(u, v)
    line = create_crv_by_interp_pnts([p1, p2], False, Settings.gtol)
    planes = create_planes_along_curve(line, maxd, nplns, None, None, None,
                                       s1, s2)
    if not planes:
        return []
    if s1 is None:
        planes.pop(0)
    if s2 is None:
        planes.pop()
    return planes


def create_plane_by_fit_points(pnts, tol=None):
    """
    Fit a plane to points.
    """
    # Convert points to array.
    pnts = array(pnts, dtype=float)
    if pnts.shape[0] < 3:
        return None

    # Put int tcol array and build plane.
    tcol_pnts = to_tcolgp_harray1_pnt(pnts)
    from OCC.GeomPlate import GeomPlate_BuildAveragePlane

    avg_pln = GeomPlate_BuildAveragePlane(tcol_pnts.GetHandle(),
                                          tcol_pnts.Length(), Settings.gtol,
                                          1, 1)
    if not avg_pln.IsPlane():
        return None

    gp_pln = avg_pln.Plane().GetObject().Pln()

    # Move to centroid.
    # Calculate average to use as the plane origin.
    pcg = Point(*mean(pnts, axis=0))
    gp_pln.SetLocation(pcg)

    # Create the plane.
    pln = Plane(gp_pln)

    # Check that all points are within tolerance to the plane.
    if not tol:
        return pln

    for p in pnts:
        p = gp_Pnt(*p)
        di = gp_pln.Distance(p)
        if di > tol:
            return None
    return pln


def create_srf_by_approx_crvs(curves, mind, maxd, tol3d, tol2d, niter,
                              method='chord', continuity='C2'):
    """
    Interpolate curves.
    """
    # Initialize approximation tool.
    app_tool = GeomFill_AppSurf(mind, maxd, tol3d, tol2d, niter)

    # Set parametrization type.
    try:
        parm_type = _occ_parm_type[method.lower()]
    except (KeyError, AttributeError):
        parm_type = Approx_ChordLength
    app_tool.SetParType(parm_type)

    # Set continuity.
    try:
        cont = _occ_continuity[continuity.upper()]
    except (KeyError, AttributeError):
        cont = GeomAbs_C2
    app_tool.SetContinuity(cont)

    # Use section generator to make all curves compatible.
    sec_gen = GeomFill_SectionGenerator()
    for c in curves:
        sec_gen.AddCurve(c.handle)
    sec_gen.Perform(Settings.ptol)

    # Create line tool.
    line_tool = GeomFill_Line(len(curves))

    # Perform the approximation.
    app_tool.Perform(line_tool.GetHandle(), sec_gen, False)
    if not app_tool.IsDone():
        return None

    # Create a surface.
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
    return s


def create_srf_by_interp_crvs(curves, q=3, method='chord'):
    """
    Create surface by interpolating curves.
    """
    ncrvs = len(curves)
    if ncrvs - 1 < q:
        q = ncrvs - 1

    # Make curves compatible using section generator.
    sec_gen = GeomFill_SectionGenerator()
    for c in curves:
        sec_gen.AddCurve(c.handle)
    sec_gen.Perform(Settings.ptol)

    # Use old method since OCC struggles to interpolate curves of low order.
    # Gather all data for old method.
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
        if method.lower() in ['u', 'uniform']:
            vknots = uniform(pnts_matrix[i, :], 0., 1.)
        elif method.lower() in ['ch', 'chord']:
            vknots = chord_length(pnts_matrix[i, :], 0., 1.)
        else:
            vknots = centripetal(pnts_matrix[i, :], 0., 1.)
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

    # Perform n + 1 interpolations in v-direction to generate surface control
    # points.
    cpw = zeros((n + 1, m + 1, 4), dtype=float)
    for i in range(0, n + 1):
        qp = pnts_matrix[i]
        # Set up system of linear equations.
        a = zeros((m + 1, m + 1), dtype=float)
        for j in range(0, m + 1):
            # span1, _ = CLib.bsplclib_LocateParameter(pnt, tcol_vknots,
            #                                          tcol_vmult,
            #                                         vknots[j], False)
            # Evaluate basis function. Figure out how to use math_Matrix and
            # use OCC in future.
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

    return s
