from math import radians, tan

from OCC.BOPAlgo import BOPAlgo_BOP, BOPAlgo_FUSE
from OCC.BRep import BRep_Tool
from OCC.BRepAdaptor import BRepAdaptor_CompCurve
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeWire
from OCC.BRepOffsetAPI import BRepOffsetAPI_MakeOffset
from OCC.GCPnts import GCPnts_AbscissaPoint

from ..bulkhead import Bulkhead
from ..curve_part import CurvePart
from ..floor import Floor
from ..frame import Frame
from ..rib import Rib
from ..skin import Skin
from ..spar import Spar
from ..stiffeners import Stiffener1D, Stiffener2D
from ..surface_part import SurfacePart
from ...geometry import CheckGeom, CreateGeom, IntersectGeom
from ...oml import CheckOML
from ...topology import ShapeTools
from ...utils import pairwise

_brep_tool = BRep_Tool()


def _form_with_solid(shape, solid):
    """
    Form common shapes.
    """
    # Perform BOP Common and get the resulting faces.
    faces = ShapeTools.bcommon(shape, solid, 'face')
    if not faces:
        faces = _hack_form_planar(shape, solid)
    if not faces:
        return None

    # Put all faces in a compound
    compound = ShapeTools.make_compound(faces)
    return compound


def _hack_form_planar(shape, solid):
    """
    Hopefully a temporary hack for a failing BOP with a planar shape.
    """
    edges = ShapeTools.bsection(shape, solid, 'edge')
    if not edges:
        return []

    wires = ShapeTools.connect_edges(edges)
    if not wires:
        return []

    faces = []
    for w in wires:
        if not w.Closed():
            continue
        builder = BRepBuilderAPI_MakeFace(w, True)
        if builder.IsDone():
            faces.append(builder.Face())
    return faces


def _build_surface_part(surface_shape, bodies):
    """
    Build surface part shape.
    """
    # Find common shapes.
    shapes = []
    for body in bodies:
        shape = _form_with_solid(surface_shape, body)
        if not shape:
            continue
        shapes.append(shape)
    if len(shapes) == 0:
        return None

    # Fuse shapes.
    if len(shapes) > 1:
        # Use BOP Fuse algorithm to fuse together all common shapes.
        bop = BOPAlgo_BOP()
        bop.SetOperation(BOPAlgo_FUSE)
        bop.AddArgument(shapes[0])
        for shape in shapes[1:]:
            bop.AddTool(shape)
        bop.Perform()
        if bop.ErrorStatus() != 0:
            return None
        shape = bop.Shape()
    else:
        shape = shapes[0]

    return shape


def create_curve_part(label, curve_shape):
    """
    Create a curve part.
    """
    cref = None
    if CheckGeom.is_curve_like(curve_shape):
        cref = curve_shape

    curve_shape = ShapeTools.to_shape(curve_shape)
    if not curve_shape:
        return None

    if not cref:
        cref = ShapeTools.curve_of_shape(curve_shape)

    return CurvePart(label, curve_shape, cref)


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
        wire = wires[0]
    else:
        wire = ShapeTools.longest_wire(wires)

    return create_curve_part(label, wire)


def create_surface_part(label, surface_shape, bodies=()):
    """
    Create a surface part.
    """
    sref = None
    if CheckGeom.is_surface_like(surface_shape):
        sref = surface_shape

    surface_shape = ShapeTools.to_shape(surface_shape)
    if not surface_shape:
        return None

    if not sref:
        sref = ShapeTools.surface_of_shape(surface_shape)

    # Return shape is no bodies are provided.
    if not bodies:
        return SurfacePart(label, surface_shape, sref=sref)

    # Build part shape.
    shape = _build_surface_part(surface_shape, bodies)
    if not shape:
        return None

    return SurfacePart(label, shape, sref=sref)


def create_wing_part_by_params(etype, label, wing, u1, v1, u2, v2,
                               surface_shape, bodies=()):
    """
    Create a wing part by parameters.
    """
    if not CheckOML.is_wing(wing):
        return None

    sref = None
    if CheckGeom.is_surface_like(surface_shape):
        sref = surface_shape

    # If the reference surface is None, use a plane normal to the wing
    # reference surface at (u1, v1). If it is surface-like, convert it to a
    # face.
    surface_shape = ShapeTools.to_shape(surface_shape)
    if surface_shape is None:
        pln = wing.extract_plane((u1, v1), (u2, v2))
        sref = pln
        if not pln:
            return None
        surface_shape = ShapeTools.to_face(pln)
    if not surface_shape:
        return None

    # Build reference surface.
    if not sref:
        sref = ShapeTools.surface_of_shape(surface_shape)

    # Build reference curve.
    cref = wing.extract_curve((u1, v1), (u2, v2), surface_shape)

    # Build part shape.
    shape = _build_surface_part(surface_shape, [wing] + list(bodies))
    if not shape:
        return None

    # Create the wing part.
    if etype in ['spar']:
        wing_part = Spar(label, shape, cref, sref)
    else:
        wing_part = Rib(label, shape, cref, sref)

    return wing_part


def create_wing_part_by_points(etype, label, wing, p1, p2, surface_shape,
                               bodies=()):
    """
    Create a wing part between points.
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
                                      surface_shape, bodies)


def create_wing_part_by_sref(etype, label, wing, surface_shape, bodies=()):
    """
    Create wing part by reference surface.
    """
    if not CheckOML.is_wing(wing):
        return None

    sref = None
    if CheckGeom.is_surface_like(surface_shape):
        sref = surface_shape

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

    # Use only one wire and create reference curve.
    w = ShapeTools.longest_wire(wires)
    e = ShapeTools.concatenate_wire(w)
    cref = ShapeTools.curve_of_edge(e)

    # Orient cref such that its first point is closest to the corner of
    # the wing.
    umin, vmin = wing.sref.u1, wing.sref.v1
    p0 = wing.eval(umin, vmin)
    p1 = cref.eval(cref.u1)
    p2 = cref.eval(cref.u2)
    d1 = p0.distance(p1)
    d2 = p0.distance(p2)
    if d2 < d1:
        cref.reverse()

    # Build reference surface.
    if not sref:
        sref = ShapeTools.surface_of_shape(surface_shape)

    # Build part shape.
    shape = _build_surface_part(surface_shape, [wing] + list(bodies))
    if not shape:
        return None

    # Create the wing part.
    if etype in ['spar']:
        wing_part = Spar(label, shape, cref, sref)
    else:
        wing_part = Rib(label, shape, cref, sref)

    return wing_part


def create_wing_part_between_geom(etype, label, wing, geom1, geom2,
                                  surface_shape, bodies=()):
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
                                      surface_shape, bodies)


def create_frame_by_sref(label, fuselage, surface_shape, h):
    """
    Create a frame using a reference shape.
    """
    if not CheckOML.is_fuselage(fuselage):
        return None

    sref = None
    if CheckGeom.is_surface_like(surface_shape):
        sref = surface_shape

    surface_shape = ShapeTools.to_shape(surface_shape)
    if not surface_shape:
        return None

    # Find initial face using BOP Common.
    faces = ShapeTools.bcommon(fuselage, surface_shape, 'face')
    if not faces:
        return None

    outer_face = ShapeTools.largest_face(faces)
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

    # Build reference surface.
    if not sref:
        sref = ShapeTools.surface_of_shape(surface_shape)

    # Build reference curve.
    # TODO Build reference curve of frame using outer wire.
    cref = ShapeTools.curve_of_wire(outer_wire)

    # Create the frame.
    frame = Frame(label, shape, cref, sref)

    return frame


def create_bulkhead_by_sref(label, fuselage, surface_shape, bodies=()):
    """
    Create a bulkhead using a reference shape.
    """
    if not CheckOML.is_fuselage(fuselage):
        return None

    sref = None
    if CheckGeom.is_surface_like(surface_shape):
        sref = surface_shape

    surface_shape = ShapeTools.to_shape(surface_shape)
    if not surface_shape:
        return None

    # Build part shape.
    shape = _build_surface_part(surface_shape, [fuselage] + list(bodies))
    if not shape:
        return None

    if not sref:
        sref = ShapeTools.surface_of_shape(surface_shape)

    # Create bulkhead.
    return Bulkhead(label, shape, sref=sref)


def create_floor_by_sref(label, fuselage, surface_shape, bodies=()):
    """
    Create a floor using a reference shape.
    """
    if not CheckOML.is_fuselage(fuselage):
        return None

    sref = None
    if CheckGeom.is_surface_like(surface_shape):
        sref = surface_shape

    surface_shape = ShapeTools.to_face(surface_shape)
    if not surface_shape:
        return None

    # Build part shape.
    shape = _build_surface_part(surface_shape, [fuselage] + list(bodies))
    if not shape:
        return None

    if not sref:
        sref = ShapeTools.surface_of_shape(surface_shape)

    # Create floor.
    return Floor(label, shape, sref=sref)


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
                                             geom2, pln)
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
                                             geom2, pln)
        if not part:
            continue
        parts.append(part)
        indx += 1

    return parts

    # TODO Create wing parts between shapes.


def create_stiffener1d(surface_part, stiffener, label):
    """
    Create a 1-D stiffener on a surface part. 
    """
    if not isinstance(surface_part, SurfacePart):
        return None

    cref = None
    if CheckGeom.is_curve_like(stiffener):
        cref = stiffener

    # If stiffener is already a Stiffener1D part, split the two parts.
    # Otherwise, find the intersection and create the stiffener.
    if not isinstance(stiffener, Stiffener1D):
        shape = ShapeTools.to_shape(stiffener)
        if not shape:
            return None
        edges = ShapeTools.bsection(surface_part, shape, 'edge')
        if not edges:
            return None
        wires = ShapeTools.connect_edges(edges)
        if not wires:
            return None
        if len(wires) == 1:
            curve_shape = wires[0]
            w = curve_shape
        else:
            curve_shape = ShapeTools.make_compound(wires)
            w = ShapeTools.longest_wire(wires)

        # Build reference curve from longest wire.
        if not cref:
            cref = ShapeTools.curve_of_wire(w)

        # Build stiffener
        stiffener = Stiffener1D(label, curve_shape, cref)

    # Split the parts.
    surface_part.split(stiffener)

    # Add to surface part.
    surface_part.add_subpart(stiffener.label, stiffener)

    return stiffener


def create_stiffener2d_by_section(label, surface_part, surface_shape, h,
                                  runout_angle):
    """
    Create a 2-D stiffener on a surface part by intersection. 
    """
    if not isinstance(surface_part, SurfacePart):
        return None

    if h <= 0. or runout_angle <= 0:
        return None

    sref = None
    if CheckGeom.is_surface_like(surface_shape):
        sref = surface_shape

    surface_shape = ShapeTools.to_shape(surface_shape)
    if not surface_shape:
        return None

    # Build spine
    edges = ShapeTools.bsection(surface_part, surface_shape, 'edge')
    if not edges:
        return None
    wires = ShapeTools.connect_edges(edges)
    if not wires:
        return None
    if len(wires) == 1:
        spine = wires[0]
    else:
        spine = ShapeTools.longest_wire(wires)

    # Build reference curve from spine.
    cref = ShapeTools.curve_of_wire(spine)

    # Calculate dx along spine to get proper run-out angle
    dx = h / tan(radians(runout_angle))

    # Find points at dx from each end.
    adp_crv = BRepAdaptor_CompCurve(spine)
    abs_pnt = GCPnts_AbscissaPoint(adp_crv, dx, adp_crv.FirstParameter())
    u = abs_pnt.Parameter()
    p1 = adp_crv.Value(u)
    abs_pnt = GCPnts_AbscissaPoint(adp_crv, -dx, adp_crv.LastParameter())
    u = abs_pnt.Parameter()
    p2 = adp_crv.Value(u)

    # Create profiles
    profile1 = ShapeTools.shape_normal_profile(p1, surface_part, h)
    profile2 = ShapeTools.shape_normal_profile(p2, surface_part, h)
    if not profile1 or not profile2:
        return None

    # Make pipe shell
    shape = ShapeTools.make_pipe_shell(spine, [profile1, profile2],
                                       surface_shape, True)
    if not shape:
        return None

    if not sref:
        sref = ShapeTools.surface_of_shape(shape)

    # TODO Trim shape with runout angles.

    # Create stiffener.
    stiffener = Stiffener2D(label, shape, cref, sref)

    # Fuse the surface part and the stiffener.
    surface_part.fuse(stiffener)
    # bop = ShapeTools.bfuse(surface_part, stiffener, 'builder')
    # if not bop.IsDone():
    #     return None
    # surface_part.set_shape(bop.Shape())
    #
    # # Reshape stiffener part.
    # stiffener.reshape(bop)

    # Add as subpart.
    # surface_part.add_subpart(stiffener.label, stiffener)
    # TODO How to handle stiffener 2-d as subpart.

    return stiffener
