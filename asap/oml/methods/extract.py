from OCC.BRep import BRep_Tool
from OCC.BRepAdaptor import BRepAdaptor_Curve
from OCC.BRepAlgoAPI import BRepAlgoAPI_Section
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.ShapeAnalysis import ShapeAnalysis_FreeBounds_ConnectEdgesToWires
from OCC.TopAbs import TopAbs_EDGE
from OCC.TopExp import TopExp_Explorer
from OCC.TopTools import Handle_TopTools_HSequenceOfShape, \
    TopTools_HSequenceOfShape
from OCC.TopoDS import topods_Edge, topods_Wire

from ...geometry import CheckGeom, CreateGeom, ProjectGeom
from ...geometry.methods.create import create_nurbs_curve_from_occ
from ...geometry.methods.distance import curve_nearest_point
from ...topology import ShapeTools

_brep_tool = BRep_Tool()


def extract_wing_plane(wing, uv1, uv2):
    """
    Create a planar surface between the wing parameters.
    """
    if CheckGeom.is_point_like(uv1):
        uv1 = ProjectGeom.invert(uv1, wing.sref)
    if CheckGeom.is_point_like(uv2):
        uv2 = ProjectGeom.invert(uv2, wing.sref)
    p0 = wing.eval(uv1[0], uv1[1])
    p1 = wing.eval(uv2[0], uv2[1])
    vn = wing.norm(uv1[0], uv1[1])
    p2 = p0.xyz + vn.ijk
    return CreateGeom.plane_by_points(p0, p1, p2)


def extract_wing_ref_curve(wing, uv1, uv2, rshape):
    """
    Extract a curve on the wing reference surface.
    """
    # Wing reference surface.
    sref = wing.sref
    sref_shape = wing.sref_shape
    if not CheckGeom.is_surface(sref):
        return None

    # Generate points and/or parameters.
    p1, p2 = None, None
    if not CheckGeom.is_point_like(uv1):
        p1 = sref.eval(*uv1)
    elif CheckGeom.is_point_like(uv1):
        p1 = uv1
        uv1 = ProjectGeom.invert(p1, wing.sref)
        if uv1[0] is None:
            return None
    if not CheckGeom.is_point_like(uv2):
        p2 = sref.eval(*uv2)
    elif CheckGeom.is_point_like(uv2):
        p2 = uv2
        uv2 = ProjectGeom.invert(p2, wing.sref)
        if uv2[0] is None:
            return None

    # Generate intersecting shape.
    if CheckGeom.is_surface_like(rshape):
        rshape = BRepBuilderAPI_MakeFace(rshape.handle, 0.).Face()
    elif rshape is None:
        pln = extract_wing_plane(wing, uv1, uv2)
        rshape = BRepBuilderAPI_MakeFace(pln.handle, 0.).Face()

    # Generate section edges with BOP Section using sref shape.
    bop = BRepAlgoAPI_Section(rshape, sref_shape)
    if bop.ErrorStatus() != 0:
        return None
    bop.RefineEdges()
    shape = bop.Shape()

    # Gather all the edges and build wire(s).
    top_exp = TopExp_Explorer(shape, TopAbs_EDGE)
    edges = TopTools_HSequenceOfShape()
    tol = []
    while top_exp.More():
        ei = topods_Edge(top_exp.Current())
        edges.Append(ei)
        top_exp.Next()
        tol.append(ShapeTools.get_tolerance(ei, 1))
    hwires = Handle_TopTools_HSequenceOfShape()
    ShapeAnalysis_FreeBounds_ConnectEdgesToWires(edges.GetHandle(),
                                                 max(tol), False, hwires)
    crvs = []
    wires_obj = hwires.GetObject()
    for i in range(1, wires_obj.Length() + 1):
        wire = topods_Wire(wires_obj.Value(i))
        edge = ShapeTools.concatenate_wire(wire)
        adp_crv = BRepAdaptor_Curve(edge)
        occ_crv = adp_crv.BSpline().GetObject()
        crv = create_nurbs_curve_from_occ(occ_crv)
        crvs.append(crv)

    if len(crvs) == 1:
        crv = crvs[0]
    else:
        # Find curve nearest p1.
        crv = curve_nearest_point(p1, crvs)

    # Segment curve between p1 and p2 and reverse if necessary.
    u1 = ProjectGeom.invert(p1, crv)
    u2 = ProjectGeom.invert(p2, crv)
    if u1 > u2:
        crv.reverse()
        u1, u2 = crv.reversed_u(u1), crv.reversed_u(u2)
    crv.segment(u1, u2)
    return crv
