from OCC.BRep import BRep_Tool
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeWire
from OCC.BRepOffsetAPI import BRepOffsetAPI_MakeOffset
from OCC.GeomAdaptor import GeomAdaptor_Curve

from ..bulkhead import Bulkhead
from ..floor import Floor
from ..frame import Frame
from ..rib import Rib
from ..skin import Skin
from ..spar import Spar
from ..surface_part import SurfacePart
from ...geometry import CheckGeom, CreateGeom, IntersectGeom
from ...geometry.methods.create import create_nurbs_curve_from_occ
from ...oml import CheckOML
from ...topology import ShapeTools
from ...utils import pairwise

_brep_tool = BRep_Tool()


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

    # Create the wing part.
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


def create_wing_part_by_sref(etype, name, wing, rshape, build):
    """
    Create wing part by reference surface.
    """
    if not CheckOML.is_wing(wing):
        return None

    rshape = ShapeTools.to_shape(rshape)
    if not rshape:
        return None

    # Intersect the wing reference surface with the reference shape.
    edges = ShapeTools.bsection(rshape, wing.sref, 'edge')
    if not edges:
        return None

    # Build wires
    tol = ShapeTools.get_tolerance(rshape)
    wires = ShapeTools.connect_edges(edges, tol)
    if not wires:
        return None

    # Use only one wire and concatenate and create reference curve.
    w = wires[0]
    e = ShapeTools.concatenate_wire(w)
    crv_data = _brep_tool.Curve(e)
    hcrv = crv_data[0]
    adp_crv = GeomAdaptor_Curve(hcrv)
    occ_crv = adp_crv.BSpline().GetObject()
    cref = create_nurbs_curve_from_occ(occ_crv)

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
        wing_part = Spar(name, wing, rshape)
    else:
        wing_part = Rib(name, wing, rshape)
    wing_part.set_cref(cref)
    wing_part.form(wing)

    if build:
        wing_part.build()

    return wing_part


def create_wing_part_between_geom(etype, name, wing, geom1, geom2, rshape,
                                  build):
    """
    Create a wing part between geometry.
    """
    if not CheckOML.is_wing(wing):
        return None

    rshape = ShapeTools.to_shape(rshape)
    if not rshape:
        return None

    # Intersect the wing reference surface with the reference shape.
    edges = ShapeTools.bsection(rshape, wing.sref, 'edge')
    if not edges:
        return None

    # Build wires
    tol = ShapeTools.get_tolerance(rshape)
    wires = ShapeTools.connect_edges(edges, tol)
    if not wires:
        return None

    # Use only one wire and concatenate and create reference curve.
    w = wires[0]
    e = ShapeTools.concatenate_wire(w)
    crv_data = _brep_tool.Curve(e)
    hcrv = crv_data[0]
    adp_crv = GeomAdaptor_Curve(hcrv)
    occ_crv = adp_crv.BSpline().GetObject()
    cref = create_nurbs_curve_from_occ(occ_crv)

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

    return create_wing_part_by_points(etype, name, wing, p1, p2, rshape, build)


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
    e = ShapeTools.concatenate_wire(w)
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


def create_skin_from_body(name, body):
    """
    Create skin from outer shell of body.
    """
    if not CheckOML.is_body(body):
        return None

    outer_shell = body.shell
    skin = Skin(name, outer_shell)
    skin.set_shape(outer_shell)

    return skin


def create_frames_between_planes(name, fuselage, planes, h, maxd=None,
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
        fname = ' '.join([name, str(indx)])
        frame = create_frame_by_sref(fname, fuselage, pln, h)
        if not frame:
            continue
        frames.append(frame)
        indx += 1

    return frames


def create_frames_at_shapes(name, fuselage, shapes, h, indx=1):
    """
    Create frames at shapes.
    """
    if not CheckOML.is_fuselage(fuselage):
        return []

    frames = []
    for shape in shapes:
        shape = ShapeTools.to_shape(shape)
        fname = ' '.join([name, str(indx)])
        frame = create_frame_by_sref(fname, fuselage, shape, h)
        if not frame:
            continue
        frames.append(frame)
        indx += 1

    return frames


def create_wing_parts_between_planes(etype, name, wing, planes, geom1, geom2,
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
        rname = ' '.join([name, str(indx)])
        part = create_wing_part_between_geom(etype, rname, wing, geom1,
                                             geom2, pln, True)
        if not part:
            continue
        parts.append(part)
        indx += 1

    return parts


def create_wing_parts_along_curve(etype, name, wing, curve, geom1, geom2,
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
        rname = ' '.join([name, str(indx)])
        part = create_wing_part_between_geom(etype, rname, wing, geom1,
                                             geom2, pln, True)
        if not part:
            continue
        parts.append(part)
        indx += 1

    return parts
