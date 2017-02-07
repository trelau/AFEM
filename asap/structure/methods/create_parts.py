from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeWire
from OCC.BRepOffsetAPI import BRepOffsetAPI_MakeOffset
from OCC.BRepAlgo import brepalgo_ConcatenateWireC0

from ..bulkhead import Bulkhead
from ..floor import Floor
from ..frame import Frame
from ..rib import Rib
from ..spar import Spar
from ..surface_part import SurfacePart
from ...geometry import CheckGeom
from ...oml import CheckOML
from ...topology import ShapeTools


def create_surface_part(name, rshape, *bodies):
    """
    Create a surface frame.
    """
    rshape = ShapeTools.to_shape(rshape)
    if not rshape:
        return None

    part = SurfacePart(name, rshape)

    for b in bodies:
        part.form(b)

    return part


def create_wing_part_by_params(etype, name, wing, u1, v1, u2, v2, rshape,
                               build):
    """
    Create a spar by parameters.
    """
    if not CheckOML.is_wing(wing):
        return None

    # If the reference surface is None, use a plane normal to the wing
    # reference surface at (u1, v1). If it is surface-like, convert it to a
    # face.
    rshape = ShapeTools.to_shape(rshape)
    if rshape is None:
        pln = wing.extract_plane((u1, v1), (u2, v2))
        if not pln:
            return None
        rshape = ShapeTools.to_face(pln)
    if not rshape:
        return None

    # Create the wing frame.
    if etype in ['spar']:
        wing_part = Spar(name, wing, rshape)
    else:
        wing_part = Rib(name, wing, rshape)

    # Set reference curve.
    cref = wing.extract_curve((u1, v1), (u2, v2), rshape)
    if cref:
        wing_part.set_cref(cref)

    # Form with wing.
    wing_part.form(wing)

    if build:
        wing_part.build()

    return wing_part


def create_wing_part_by_points(etype, name, wing, p1, p2, rshape, build):
    """
    Create a spar between points.
    """
    if not CheckOML.is_wing(wing):
        return None

    p1 = CheckGeom.to_point(p1)
    p2 = CheckGeom.to_point(p2)
    if not CheckGeom.is_point(p1) or not CheckGeom.is_point(p2):
        return None

    u1, v1 = wing.invert_point(p1)
    u2, v2 = wing.invert_point(p2)
    if None in [u1, v1, u2, v2]:
        return None

    return create_wing_part_by_params(etype, name, wing, u1, v1, u2, v2,
                                      rshape, build)


def create_frame_by_sref(name, fuselage, rshape, h):
    """
    Create a frame using a reference shape.
    """
    rshape = ShapeTools.to_face(rshape)
    if not rshape:
        return None

    if not CheckOML.is_fuselage(fuselage):
        return None

    # Find initial face using BOP Common.
    faces = ShapeTools.bcommon(fuselage, rshape, 'face')
    if not faces:
        return None

    f = faces[0]
    offset = BRepOffsetAPI_MakeOffset(f)
    offset.Perform(-h)
    if not offset.IsDone():
        return None
    w = ShapeTools.to_wire(offset.Shape())
    if not w:
        return None
    # Force concatenation of wire to avoid small edges.
    e = brepalgo_ConcatenateWireC0(w)
    w = BRepBuilderAPI_MakeWire(e).Wire()
    builder = BRepBuilderAPI_MakeFace(f)
    builder.Add(w)
    face = builder.Face()
    # TODO Offset wires are causing issues for STEP export.

    # Make the face a shell.
    shell = ShapeTools.to_shell(face)

    # Create the frame.
    frame = Frame(name, fuselage, rshape)

    # Set Frame shape to shell.
    frame.set_shape(shell)
    return frame


def create_bulkhead_by_sref(name, fuselage, rshape, build):
    """
    Create a bulkhead using a reference shape.
    """
    if not CheckOML.is_fuselage(fuselage):
        return None

    rshape = ShapeTools.to_face(rshape)
    if not rshape:
        return None

    # Create bulkhead.
    bulkhead = Bulkhead(name, fuselage, rshape)

    # Form with fuselage.
    bulkhead.form(fuselage)

    if build:
        bulkhead.build()

    return bulkhead


def create_floor_by_sref(name, fuselage, rshape, build):
    """
    Create a floor using a reference shape.
    """
    if not CheckOML.is_fuselage(fuselage):
        return None

    rshape = ShapeTools.to_face(rshape)
    if not rshape:
        return None

    # Create bulkhead.
    floor = Floor(name, fuselage, rshape)

    # Form with fuselage.
    floor.form(fuselage)

    if build:
        floor.build()

    return floor
