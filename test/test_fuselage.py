from __future__ import print_function

from asap.geometry import CreateGeom
from asap.graphics import Viewer
from asap.io import ImportVSP
from asap.structure import CreatePart, AssemblyData, PartTools
from asap.topology import ShapeTools

# Import model
fn = './test_io/777-200LR_mod_vsp350_sref.stp'
ImportVSP.step_file(fn)

# Get components.
wing = ImportVSP.get_body('Wing')
fuselage = ImportVSP.get_body('Fuselage')
other_wing = ImportVSP.get_body('Wing.2')

AssemblyData.create_assy('fuselage assy')

# Fwd bulkhead at 25% wing chord
p0 = wing.eval(0.25, 0.)
pln = CreateGeom.plane_by_axes(p0, 'yz')
fwd_bh = CreatePart.bulkhead.by_sref('fwd bh', fuselage, pln)

# Rear bulkhead at 65% wing chord
p0 = wing.eval(0.65, 0.)
pln = CreateGeom.plane_by_axes(p0, 'yz')
rear_bh = CreatePart.bulkhead.by_sref('rear bh', fuselage, pln)

# Aft bulkhead at 85% wing chord
p0 = wing.eval(0.85, 0.)
pln = CreateGeom.plane_by_axes(p0, 'yz')
aft_bh = CreatePart.bulkhead.by_sref('aft bh', fuselage, pln)

# Floor
pln = CreateGeom.plane_by_axes([0, 0, -24], 'xy')
floor = CreatePart.floor.by_sref('floor', fuselage, pln)

# Skin
fskin = CreatePart.skin.from_body('fuselage skin', fuselage)

# Cutting solids.
p0 = wing.eval(0., 0.)
fwd_pln = CreateGeom.plane_by_axes(p0, 'yz')
fwd_cut = ShapeTools.box_from_plane(fwd_pln, 1e6, 1e6, -1e6)

p0 = wing.eval(1., 0.)
aft_pln = CreateGeom.plane_by_axes(p0, 'yz')
aft_cut = ShapeTools.box_from_plane(aft_pln, 1e6, 1e6, 1e6)

above_floor = ShapeTools.make_halfspace(floor, [0, 0, 1e6])
below_floor = ShapeTools.make_halfspace(floor, [0, 0, -1e6])

xz_pln = CreateGeom.plane_by_axes([0., 0., 0.], 'xz')
left_cut = ShapeTools.box_from_plane(xz_pln, 1e6, 1e6, -1e6)

# Frames
plns = [fwd_pln, fwd_bh.sref, rear_bh.sref, aft_bh.sref, aft_pln]
frames = CreatePart.frame.between_planes('frame', fuselage, plns, 4., 24.)

# Frames at bulkheads
plns = [fwd_bh.sref, rear_bh.sref, aft_bh.sref]
indx = len(frames) + 1
bh_frames = CreatePart.frame.at_shapes('frame', fuselage, plns, 4., indx)
PartTools.cut_parts(bh_frames, below_floor)
frames += bh_frames

# Cut and discard to get half model
PartTools.cut_parts([fskin, floor], fwd_cut)
PartTools.cut_parts([fskin, floor], aft_cut)
PartTools.cut_parts([fwd_bh, rear_bh, aft_bh], above_floor)

# Cut a piece of the wing to cut fuselage faces inside it.
cut1 = ShapeTools.box_from_plane(fwd_bh.sref, 1e6, 1e6, -1e6)
cut2 = ShapeTools.box_from_plane(aft_bh.sref, 1e6, 1e6, 1e6)
joined_wing = wing.fuse(other_wing)
cutter = joined_wing.bop_algo([cut1, cut2], 'cut')
PartTools.cut_parts([fskin, fwd_bh, rear_bh, aft_bh, floor] + frames, cutter)

# # Join
PartTools.fuse_parts([fwd_bh, rear_bh, aft_bh, floor, fskin] + frames)

# Cut out windows
# edges = ShapeTools.bsection(floor, xz_pln, 'edge')
# wire = ShapeTools.connect_edges(edges)[0]
#
# pnts = ShapeTools.points_along_shape(wire, 24., s1=12., s2=-12.)
# v = CreateGeom.vector([0, 0, 48])
# dn = CreateGeom.direction_by_xyz(0, 1, 0)
# vp1 = CreateGeom.vector([0, 1000, 0])
# vp2 = CreateGeom.vector([0, -1000, 0])
# from OCC.GC import GC_MakeCircle
# from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeWire, BRepBuilderAPI_MakeFace
# windows = []
# for p in pnts:
#     # Translate to centroid
#     p.translate(v)
#
#     # Create circle
#     circle = GC_MakeCircle(p, dn, 8).Value().GetObject()
#
#     # Make face from circle.
#     edge = ShapeTools.to_edge(circle)
#     wire = BRepBuilderAPI_MakeWire(edge).Wire()
#     wire = ShapeTools.fix_shape(wire)
#     face = BRepBuilderAPI_MakeFace(wire, True).Face()
#     face = ShapeTools.fix_shape(face)
#
#     # Make prism
#     prism1 = ShapeTools.make_prism(face, [0, 200, 0])
#     prism2 = ShapeTools.make_prism(face, [0, -200, 0])
#     prism = ShapeTools.bfuse(prism1, prism2)
#
#     # Viewer.display(prism)
#     windows.append(prism)
#
# # Cut the skin.
# windows = ShapeTools.make_compound(windows)
# fskin.cut(windows)

# Viewing
fskin.set_transparency(0.5)
floor.set_transparency(0.5)
Viewer.add_items(fskin, floor, fwd_bh, rear_bh, aft_bh, *frames)
Viewer.show()

AssemblyData.mesh_assy(maxh=4., quad_dominated=False)
Viewer.add_meshes(*AssemblyData.get_parts())
Viewer.show_mesh()
