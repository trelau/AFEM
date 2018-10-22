import time

from afem.config import Settings
from afem.geometry import *
from afem.graphics import Viewer
from afem.oml import Body
from afem.structure import *
from afem.topology import *

Settings.log_to_console()

# Import geometry and get bodies
fn = r'../../models/777-200LR.xbf'
bodies = Body.load_bodies(fn)
wing = bodies['Wing']
fuselage = bodies['Fuselage']
htail = bodies['Htail']

wing.set_transparency(0.5)
fuselage.set_transparency(0.5)

# WING
wing_group = GroupAPI.create_group('wing group')
rib_group = GroupAPI.create_group('rib group', wing_group)

root = RibByParameters('root', 0.15, 0.05, 0.70, 0.05, wing).part
tip = RibByParameters('tip', 0.15, 0.2, 0.70, 0.2, wing).part

spar_group = GroupAPI.create_group('spar group', wing_group)

fspar = SparByPoints('fspar', root.cref.p1, tip.cref.p1, wing).part
rspar = SparByPoints('rspar', root.cref.p2, tip.cref.p2, wing).part

rib_group.activate()
RibsAlongCurveByNumber('rib', rspar.cref, 5, fspar.shape, rspar.shape, wing,
                       d1=30., d2=-50)

# Join groups and discard
FuseGroups([rib_group, spar_group])
DiscardByCref(wing_group.get_parts())

# Wing skin
wing_group.activate()
skin = SkinByBody('wing skin', wing).part
parts = rib_group.get_parts() + spar_group.get_parts()
skin.fuse(*parts)
skin.discard_by_dmin(wing.sref_shape, 1.)
skin.fix()
skin.set_transparency(0.5)

# FUSELAGE
fuse_group = GroupAPI.create_group('fuse group')

# Fuselage skin
skin = SkinByBody('fuselage skin', fuselage).part
skin.set_transparency(0.5)

# Discard half
pln = PlaneByAxes().plane
solid = HalfspaceBySurface(pln, (0, -1, 0)).solid
skin.discard_by_solid(solid)

# Frame
pln = PlaneByAxes(root.point_from_parameter(0.5, is_rel=True), 'yz').plane
frame = FrameByPlane('frame', pln, fuselage, 3.).part

# Join fuselage and frame
skin.fuse(frame)

# HTAIL
htail_group = GroupAPI.create_group('htail group')

root = RibByParameters('root', 0.15, 0.05, 0.70, 0.05, htail).part
tip = RibByParameters('tip', 0.15, 0.5, 0.70, 0.5, htail).part

fspar = SparByPoints('fspar', root.cref.p1, tip.cref.p1, htail).part
rspar = SparByPoints('rspar', root.cref.p2, tip.cref.p2, htail).part

pln = PlaneByAxes().plane
RibsAlongCurveByNumber('rib', rspar.cref, 5, fspar.shape, rspar.shape, htail,
                       d1=12, d2=-12, ref_pln=pln)

parts = htail_group.get_parts()
FuseSurfacePartsByCref(parts)
DiscardByCref(parts)

skin = SkinByBody('htail skin', htail).part
skin.fuse(*parts)
skin.discard_by_dmin(htail.sref_shape, 1.)
skin.fix()
skin.set_transparency(0.5)

# JOIN
start = time.time()
print('Performing assembly join...')
bop = FuseGroups([wing_group, fuse_group, htail_group])
print('Joining complete in ', time.time() - start, ' seconds.')

shape = bop.shape
tool = ExploreFreeEdges(shape)

parts = GroupAPI.get_master().get_parts()

# Mesh
the_mesh = MeshVehicle(4.)
mesh_start = time.time()
print('Computing mesh...')
status = the_mesh.compute()
if not status:
    print('Failed to compute mesh')
else:
    print('Meshing complete in ', time.time() - mesh_start, ' seconds.')

gui = Viewer()
gui.add(*tool.free_edges)
gui.add(the_mesh)
gui.start()
