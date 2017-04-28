from OCC.BRep import BRep_Tool
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeWire
from OCC.BRepOffsetAPI import BRepOffsetAPI_MakeOffset

from ..bulkhead import Bulkhead
from ..curve_part import CurvePart
from ..floor import Floor
from ..frame import Frame
from ..rib import Rib
from ..skin import Skin
from ..spar import Spar
from ..surface_part import SurfacePart
from ...geometry import CheckGeom, CreateGeom, IntersectGeom
from ...oml import CheckOML
from ...topology import ShapeTools
from ...utils import pairwise

_brep_tool = BRep_Tool()


def create_curve_part(label, curve_shape):
    """
    Create a curve part.
    """
    shape = ShapeTools.to_shape(curve_shape)
    if not shape:
        return None

    return CurvePart(label, shape)


def create_curve_part_by_section(label, shape1, shape2):
    """
    Create a curve part from the section between two shapes. 
    """
    shape1, shape2 = ShapeTools.to_shape(shape1), ShapeTools.to_shape(shape2)
    if not shape1 or not shape2:
        return None

    # Find intersection between two shapes.
    section = ShapeTools.bsection(shape1, shape2)
    if not section:
        return None

    # Create wires from section.
    wires = ShapeTools.wires_from_shape(section)
    if not wires:
        return None

    # Create shape for curve part.
    if len(wires) == 1:
        shape = wires[0]
    else:
        shape = ShapeTools.make_compound(wires)

    return CurvePart(label, shape)


def create_surface_part(label, surface_shape, *bodies):
    """
    Create a surface part.
    """
    surface_shape = ShapeTools.to_shape(surface_shape)
    if not surface_shape:
        return None

    part = SurfacePart(label, surface_shape)

    for b in bodies:
        part.form(b)

    return part


def create_wing_part_by_params(etype, label, wing, u1, v1, u2, v2,
                               surface_shape, build):
    """
    Create a spar by parameters.
    """
    if not CheckOML.is_wing(wing):
        return None

    # If the reference surface is None, use a plane normal to the wing
    # reference surface at (u1, v1). If it is surface-like, convert it to a
    # face.
    surface_shape = ShapeTools.to_shape(surface_shape)
    if surface_shape is None:
        pln = wing.extract_plane((u1, v1), (u2, v2))
        if not pln:
            return None
        surface_shape = ShapeTools.to_face(pln)
    if not surface_shape:
        return None

    # Create the wing part.
    if etype in ['spar']:
        wing_part = Spar(label, surface_shape)
    else:
        wing_part = Rib(label, surface_shape)

    # Set reference curve.
    cref = wing.extract_curve((u1, v1), (u2, v2), surface_shape)
    if cref:
        wing_part.set_cref(cref)

    # Form with wing.
    wing_part.form(wing)

    if build:
        wing_part.build()

    return wing_part


def create_wing_part_by_points(etype, label, wing, p1, p2, surface_shape,
                               build):
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

    return create_wing_part_by_params(etype, label, wing, u1, v1, u2, v2,
                                      surface_shape, build)


def create_wing_part_by_sref(etype, label, wing, surface_shape, build):
    """
    Create wing part by reference surface.
    """
    if not CheckOML.is_wing(wing):
        return None

    surface_shape = ShapeTools.to_shape(surface_shape)
    if not surface_shape:
        return None

    # Intersect the wing reference surface with the reference shape.
    edges = ShapeTools.bsection(surface_shape, wing.sref, 'edge')
    if not edges:
        return None

    # Build wires
    tol = ShapeTools.get_tolerance(surface_shape)
    wires = ShapeTools.connect_edges(edges, tol)
    if not wires:
        return None

    # Use only one wire and concatenate and create reference curve.
    w = ShapeTools.longest_wire(wires)
    e = ShapeTools.concatenate_wire(w)
    cref = ShapeTools.curve_of_edge(e)

    # Orient the cref such that its first point is closest to the corner of
    # the wing.
    umin, vmin = wing.sref.u1, wing.sref.v1
    p0 = wing.eval(umin, vmin)
    p1 = cref.eval(cref.u1)
    p2 = cref.eval(cref.u2)
    d1 = p0.distance(p1)
    d2 = p0.distance(p2)
    if d2 < d1:
        cref.reverse()

    # Create the wing part.
    if etype in ['spar']:
        wing_part = Spar(label, surface_shape)
    else:
        wing_part = Rib(label, surface_shape)
    wing_part.set_cref(cref)
    wing_part.form(wing)

    if build:
        wing_part.build()

    return wing_part


def create_wing_part_between_geom(etype, label, wing, geom1, geom2,
                                  surface_shape, build):
    """
    Create a wing part between geometry.
    """
    if not CheckOML.is_wing(wing):
        return None

    surface_shape = ShapeTools.to_shape(surface_shape)
    if not surface_shape:
        return None

    # Intersect the wing reference surface with the reference shape.
    edges = ShapeTools.bsection(surface_shape, wing.sref, 'edge')
    if not edges:
        return None

    # Build wires
    tol = ShapeTools.get_tolerance(surface_shape)
    wires = ShapeTools.connect_edges(edges, tol)
    if not wires:
        return None

    # Use only one wire and concatenate and create reference curve.
    w = ShapeTools.longest_wire(wires)
    e = ShapeTools.concatenate_wire(w)
    cref = ShapeTools.curve_of_edge(e)

    # Intersect geometry and create by points.
    tol = ShapeTools.get_tolerance(w)
    ci = IntersectGeom.perform(cref, geom1, tol)
    if not ci.success:
        return None
    p1 = ci.point(1)

    ci = IntersectGeom.perform(cref, geom2, tol)
    if not ci.success:
        return None
    p2 = ci.point(1)

    return create_wing_part_by_points(etype, label, wing, p1, p2,
                                      surface_shape, build)


def create_frame_by_sref(label, fuselage, surface_shape, h):
    """
    Create a frame using a reference shape.
    """
    surface_shape = ShapeTools.to_face(surface_shape)
    if not surface_shape:
        return None

    if not CheckOML.is_fuselage(fuselage):
        return None

    # Find initial face using BOP Common.
    faces = ShapeTools.bcommon(fuselage, surface_shape, 'face')
    if not faces:
        return None

    outer_face = faces[0]
    offset = BRepOffsetAPI_MakeOffset()
    outer_wire = ShapeTools.outer_wire(outer_face)
    offset.AddWire(outer_wire)
    # Negative offset so it trims inner part of face.
    offset.Perform(-h)
    if not offset.IsDone():
        return None
    w = ShapeTools.to_wire(offset.Shape())
    if not w:
        return None
    # Force concatenation of wire to avoid small edges.
    e = ShapeTools.concatenate_wire(w)
    w = BRepBuilderAPI_MakeWire(e).Wire()
    # Create an inner face to cut from outer face.
    builder = BRepBuilderAPI_MakeFace(w, True)
    if not builder.IsDone():
        return None
    inner_face = builder.Face()

    shape = ShapeTools.bcut(outer_face, inner_face)
    if not shape:
        return None

    # Create the frame.
    frame = Frame(label, surface_shape)

    # Set Frame shape.
    frame.set_shape(shape)
    return frame


def create_bulkhead_by_sref(label, fuselage, surface_shape, build):
    """
    Create a bulkhead using a reference shape.
    """
    if not CheckOML.is_fuselage(fuselage):
        return None

    surface_shape = ShapeTools.to_face(surface_shape)
    if not surface_shape:
        return None

    # Create bulkhead.
    bulkhead = Bulkhead(label, surface_shape)

    # Form with fuselage.
    bulkhead.form(fuselage)

    if build:
        bulkhead.build()

    return bulkhead


def create_floor_by_sref(label, fuselage, surface_shape, build):
    """
    Create a floor using a reference shape.
    """
    if not CheckOML.is_fuselage(fuselage):
        return None

    surface_shape = ShapeTools.to_face(surface_shape)
    if not surface_shape:
        return None

    # Create bulkhead.
    floor = Floor(label, surface_shape)

    # Form with fuselage.
    floor.form(fuselage)

    if build:
        floor.build()

    return floor


def create_skin_from_body(label, body, copy=True):
    """
    Create skin from outer shell of body.
    """
    if not CheckOML.is_body(body):
        return None

    outer_shell = body.shell
    if copy:
        outer_shell = ShapeTools.copy_shape(outer_shell, False)
    skin = Skin(label, outer_shell)
    skin.set_shape(outer_shell)

    return skin


def create_skin_from_solid(label, solid, copy=True):
    """
    Create skin from outer shell of solid. 
    """
    solid = ShapeTools.to_solid(solid)
    if not solid:
        return None

    outer_shell = ShapeTools.outer_shell(solid)
    if copy:
        outer_shell = ShapeTools.copy_shape(outer_shell, False)
    skin = Skin(label, outer_shell)
    skin.set_shape(outer_shell)

    return skin


def create_frames_between_planes(label, fuselage, planes, h, maxd=None,
                                 nplns=None, indx=1):
    """
    Create frames evenly spaced between planes.
    """
    if not CheckOML.is_fuselage(fuselage):
        return []

    # Generate planes between consecutive planes.
    plns = []
    for pln1, pln2 in pairwise(planes):
        plns += CreateGeom.planes_between_planes(pln1, pln2, maxd, nplns)
    if not plns:
        return []

    # Create frames.
    frames = []
    for pln in plns:
        flabel = ' '.join([label, str(indx)])
        frame = create_frame_by_sref(flabel, fuselage, pln, h)
        if not frame:
            continue
        frames.append(frame)
        indx += 1

    return frames


def create_frames_at_shapes(label, fuselage, shapes, h, indx=1):
    """
    Create frames at shapes.
    """
    if not CheckOML.is_fuselage(fuselage):
        return []

    frames = []
    for shape in shapes:
        shape = ShapeTools.to_shape(shape)
        flabel = ' '.join([label, str(indx)])
        frame = create_frame_by_sref(flabel, fuselage, shape, h)
        if not frame:
            continue
        frames.append(frame)
        indx += 1

    return frames


def create_wing_parts_between_planes(etype, label, wing, planes, geom1, geom2,
                                     maxd=None, nplns=None, indx=1):
    """
    Create wing part evenly spaced between planes.
    """
    if not CheckOML.is_wing(wing):
        return []

    # Generate planes between consecutive planes.
    plns = []
    for pln1, pln2 in pairwise(planes):
        plns += CreateGeom.planes_between_planes(pln1, pln2, maxd, nplns)
    if not plns:
        return []

    # Create wing parts.
    parts = []
    for pln in plns:
        rlabel = ' '.join([label, str(indx)])
        part = create_wing_part_between_geom(etype, rlabel, wing, geom1,
                                             geom2, pln, True)
        if not part:
            continue
        parts.append(part)
        indx += 1

    return parts


def create_wing_parts_along_curve(etype, label, wing, curve, geom1, geom2,
                                  maxd=None, npts=None, ref_pln=None,
                                  u1=None, u2=None, s1=None, s2=None, indx=1):
    """
    Create wing parts along a curve.
    """
    if not CheckOML.is_wing(wing) or not CheckGeom.is_curve_like(curve):
        return []

    # Generate planes along the curve.
    plns = CreateGeom.planes_along_curve(curve, maxd, npts, ref_pln, u1, u2,
                                         s1, s2)

    # Create wing parts.
    parts = []
    for pln in plns:
        rlabel = ' '.join([label, str(indx)])
        part = create_wing_part_between_geom(etype, rlabel, wing, geom1,
                                             geom2, pln, True)
        if not part:
            continue
        parts.append(part)
        indx += 1

    return parts

# TODO Create wing parts between shapes.
