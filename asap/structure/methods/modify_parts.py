from OCC.BRepAdaptor import BRepAdaptor_Surface
from OCC.BRepAlgoAPI import BRepAlgoAPI_Cut
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, \
    BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeWire
from OCC.BRepClass3d import BRepClass3d_SolidClassifier
from OCC.BRepExtrema import BRepExtrema_DistShapeShape, BRepExtrema_IsInFace
from OCC.BRepGProp import brepgprop
from OCC.BRepPrimAPI import BRepPrimAPI_MakePrism
from OCC.Extrema import Extrema_ExtFlag_MIN
from OCC.GC import GC_MakeCircle
from OCC.GProp import GProp_GProps
from OCC.ShapeBuild import ShapeBuild_ReShape
from OCC.TopAbs import TopAbs_IN
from OCC.gp import gp_Dir, gp_Pnt, gp_Vec

from .reshape_parts import reshape_surface_parts
from ...config import Settings
from ...geometry import CreateGeom
from ...topology import ShapeTools


def discard_faces_by_distance(part, shape, dmin=None):
    """
    Discard faces based on their distance to the given shape.
    """
    faces = part.faces
    if not faces:
        return False

    if dmin is None:
        dmin = Settings.gtol

    # For each face, check the distance to the shape.
    reshape = ShapeBuild_ReShape()
    modified = False
    for f in faces:
        dss = BRepExtrema_DistShapeShape(f, shape, Extrema_ExtFlag_MIN)
        if not dss.IsDone():
            continue
        if dss.Value() <= dmin:
            reshape.Remove(f)
            modified = True

    if not modified:
        return False

    # Remove the faces.
    new_shape = reshape.Apply(part)
    part.set_shape(new_shape)
    return True


def discard_faces_by_solid(part, solid, tol=None):
    """
    Discard faces of a part using a solid.
    """
    faces = part.faces
    if not faces:
        return False

    if tol is None:
        tol = Settings.gtol

    # For each face, check if face is inside the solid.
    reshape = ShapeBuild_ReShape()
    modified = False
    sprops = GProp_GProps()
    classify = BRepClass3d_SolidClassifier(solid)
    for f in faces:
        brepgprop.SurfaceProperties(f, sprops, Settings.gtol)
        cg = sprops.CentreOfMass()
        classify.Perform(cg, tol)
        if classify.State() == TopAbs_IN:
            reshape.Remove(f)
            modified = True

    if not modified:
        return False

    # Remove the faces.
    new_shape = reshape.Apply(part)
    part.set_shape(new_shape)
    return True


def discard_wing_part_faces(part):
    """
    Try to automatically discard wing part faces.
    """
    if not part.cref:
        return False
    cref = part.cref

    # Create vectors at each end of the reference curve pointing "out" of
    # the part.
    u1, u2 = cref.u1, cref.u2
    v1 = cref.deriv(u1, 1)
    v2 = cref.deriv(u2, 1)
    # Reverse v1 so it's "out" of the part.
    v1.reverse()

    # Create planes at each end.
    p1 = cref.eval(u1)
    p2 = cref.eval(u2)
    pln1 = CreateGeom.plane_by_normal(p1, v1.xyz)
    pln2 = CreateGeom.plane_by_normal(p2, v2.xyz)

    # Translate points to define half space.
    pref1 = p1 + 100. * v1.xyz
    pref2 = p2 + 100. * v2.xyz
    hs1 = ShapeTools.make_halfspace(pln1, pref1)
    hs2 = ShapeTools.make_halfspace(pln2, pref2)

    # Discard
    status1 = discard_faces_by_solid(part, hs1)
    status2 = discard_faces_by_solid(part, hs2)

    return status1 or status2


def cut_wing_part_with_circle(part, dx, r):
    """
    Cut a wing part with a circle.
    """
    # Calculate a point at dx from u1 on the reference curve of the part.
    u, pc = CreateGeom.point_from_other(part.cref, dx, part.cref.u1)
    if None in [u, pc]:
        return False

    # Project point to the part.
    vc = ShapeTools.to_vertex(pc)
    dss = BRepExtrema_DistShapeShape(part, vc, Extrema_ExtFlag_MIN)
    if not dss.IsDone():
        return False

    # Get the face the point is on and the normal.
    if dss.SupportTypeShape1(1) != BRepExtrema_IsInFace:
        return False
    face = dss.SupportOnShape1(1)
    u, v = dss.ParOnFaceS1(1)
    face = ShapeTools.to_face(face)
    adp_srf = BRepAdaptor_Surface(face)
    dmin = dss.Value()

    # Create a circle
    pc = gp_Pnt()
    du, dv = gp_Vec(), gp_Vec()
    adp_srf.D1(u, v, pc, du, dv)
    vn = du.Crossed(dv)
    vn.Normalize()
    # For now scale the normal since cutting with an infinite prism is
    # causing errors.
    vn.Scale(2. * dmin + 1.)
    dn = gp_Dir(du.Crossed(dv))
    circle = GC_MakeCircle(pc, dn, r).Value()
    edge = BRepBuilderAPI_MakeEdge(circle).Edge()
    wire = BRepBuilderAPI_MakeWire(edge).Wire()
    face = BRepBuilderAPI_MakeFace(wire, True).Face()

    # Make a prism.
    circle = BRepPrimAPI_MakePrism(face, vn).Shape()

    # Cut the circle from the part.
    bop = BRepAlgoAPI_Cut(part, circle)
    if bop.ErrorStatus() != 0:
        return False

    # Replace modified face(s) of result into original shapes.
    status = reshape_surface_parts(bop, [part])

    return status


def cut_part(part, cutter):
    """
    Cut a part with a shape.
    """
    if part.IsNull() or cutter.IsNull():
        return False

    bop = BRepAlgoAPI_Cut(part, cutter)
    if bop.ErrorStatus() != 0:
        return False

    # Replace modified face(s) of result into original shapes.
    return reshape_surface_parts(bop, [part])
