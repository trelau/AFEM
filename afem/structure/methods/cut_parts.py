from OCC.BRepAdaptor import BRepAdaptor_Surface
from OCC.BRepAlgoAPI import BRepAlgoAPI_Cut
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, \
    BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeWire
from OCC.BRepExtrema import BRepExtrema_DistShapeShape, BRepExtrema_IsInFace
from OCC.BRepPrimAPI import BRepPrimAPI_MakePrism
from OCC.Extrema import Extrema_ExtFlag_MIN
from OCC.GC import GC_MakeCircle
from OCC.gp import gp_Dir, gp_Pnt, gp_Vec

from .reshape_parts import reshape_parts
from ...geometry import CreateGeom
from ...topology import ShapeTools


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
    status = reshape_parts(bop, [part])

    return status


def cut_part(part, cutter):
    """
    Cut a part.
    """
    if part.IsNull() or cutter.IsNull():
        return False

    bop = BRepAlgoAPI_Cut(part, cutter)
    if bop.ErrorStatus() != 0:
        return False

    # Replace modified shape(s).
    return reshape_parts(bop, [part])


def cut_surface_parts(parts, cutter):
    """
    Cut surface parts.
    """
    if len(parts) == 0 or cutter.IsNull():
        return False

    # Put parts into compound.
    parts = list(parts)
    compound = ShapeTools.make_compound(parts)
    # Cut the parts. Putting them in a compound avoids cutting them with each
    # other.
    bop = BRepAlgoAPI_Cut(compound, cutter)
    if bop.ErrorStatus() != 0:
        return False

    # Replace modified face(s) of result into original shapes.
    return reshape_parts(bop, parts)
