from afem.config import Settings
from afem.exchange import brep
from afem.geometry import *
from afem.graphics import Viewer
from afem.smesh import *
from afem.oml import Body
from afem.structure import *
from afem.topology import *

Settings.set_units('in')
Settings.log_to_console()

# IMPORT SOLIDS ---------------------------------------------------------------
fn1 = r'..\models\uCRM\fuselage.brep'
fn2 = r'..\models\uCRM\lhs_wing.brep'
fn3 = r'..\models\uCRM\rhs_wing.brep'

shape1 = brep.read_brep(fn1)
shape2 = brep.read_brep(fn2)
shape3 = brep.read_brep(fn3)

fuselage = Body(shape1, 'fuselage')
lhs_wing = Body(shape2, 'lhs wing')
rhs_wing = Body(shape3, 'rhs wing')

fuselage.set_transparency(0.75)
fuselage.set_color(0.5, 0.5, 0.5)
lhs_wing.set_transparency(0.75)
lhs_wing.set_color(0.5, 0.5, 0.5)
rhs_wing.set_transparency(0.75)
rhs_wing.set_color(0.5, 0.5, 0.5)

# BUILD REFERENCE SURFACES ----------------------------------------------------
p1 = [908.101552, 0., 174.148325]
p2 = [1438.059406, 0., 112.513068]
c1 = NurbsCurveByPoints([p1, p2]).curve

p1 = [992.832371, 120.251998, 175.831087]
p2 = [1458.666073, 120.263028, 140.264227]
c2 = NurbsCurveByPoints([p1, p2]).curve

p1 = [1225.471431, 428.235526, 177.690824]
p2 = [1511.118796, 428.460574, 169.450682]
c3 = NurbsCurveByPoints([p1, p2]).curve

p1 = [1778.118761, 1159.253946, 182.662558]
p2 = [1885.469001, 1159.590941, 181.982083]
c4 = NurbsCurveByPoints([p1, p2]).curve

sref = NurbsSurfaceByInterp([c1, c2, c3, c4], 1).surface
rhs_wing.set_sref(sref)

p1 = [908.101552, 0., 174.148325]
p2 = [1438.059406, 0., 112.513068]
c1 = NurbsCurveByPoints([p1, p2]).curve

p1 = [992.832371, -120.251998, 175.831087]
p2 = [1458.666073, -120.263028, 140.264227]
c2 = NurbsCurveByPoints([p1, p2]).curve

p1 = [1225.471431, -428.235526, 177.690824]
p2 = [1511.118796, -428.460574, 169.450682]
c3 = NurbsCurveByPoints([p1, p2]).curve

p1 = [1778.118761, -1159.253946, 182.662558]
p2 = [1885.469001, -1159.590941, 181.982083]
c4 = NurbsCurveByPoints([p1, p2]).curve

sref = NurbsSurfaceByInterp([c1, c2, c3, c4], 1).surface
lhs_wing.set_sref(sref)

# WING STRUCTURE --------------------------------------------------------------
xz_pln = PlaneByAxes().plane
xz_face = FaceBySurface(xz_pln).face

# Use the knots of the reference surface to determine parameters of wing
# stations.
vknots = rhs_wing.sref.vknots

# Root rib at the first wing station using xz-plane. Restrict to 15 and 75%
# chord.
pln = PlaneByAxes(rhs_wing.sref.eval(0., vknots[1]), 'xz').plane
root_rib = RibByShape('root rib', pln, rhs_wing).part
p1 = root_rib.point_from_parameter(0.15, is_rel=True)
p2 = root_rib.point_from_parameter(0.75, is_rel=True)
root_rib.set_p1(p1)
root_rib.set_p2(p2)

# Tip rib
v2 = vknots[-1]
tip_rib = RibByParameters('tip rib', 0.15, v2, 0.75, v2, rhs_wing).part

# Front center spar
pln = PlaneByAxes(root_rib.cref.p1, 'yz').plane
fc_spar = SparBetweenShapes('fc spar', xz_face, root_rib.shape, rhs_wing,
                            pln).part

# Rear center spar
pln = PlaneByAxes(root_rib.cref.p2, 'yz').plane
rc_spar = SparBetweenShapes('rc spar', xz_face, root_rib.shape, rhs_wing,
                            pln).part

# Front spar is between the second and last wing stations. Use intersection
# of center spar and root rib for orientation.
pln = PlaneByIntersectingShapes(root_rib.shape, fc_spar.shape,
                                tip_rib.cref.p1).plane
fspar = SparByPoints('fspar', root_rib.cref.p1, tip_rib.cref.p1, rhs_wing,
                     pln).part

# Rear spar 1 is between second and third wing stations. Use intersection
# of center spar and root rib for orientation.
v2 = vknots[2]
p2 = rhs_wing.sref.eval(0.75, v2)
pln = PlaneByIntersectingShapes(root_rib.shape, rc_spar.shape, p2).plane
rspar1 = SparByPoints('rspar 1', root_rib.cref.p2, p2, rhs_wing, pln).part

# Rear spar 2 is between the third and last wing stations
rspar2 = SparByPoints('rspar 2', rspar1.cref.p2, tip_rib.cref.p2,
                      rhs_wing).part

parts = GroupAPI.get_parts()

# Fuse together boundary of wing box
FuseSurfacePartsByCref(parts)
DiscardByCref(parts)

# Distribute ribs along the front spar at 25"
aft_shape = CompoundByShapes([root_rib.shape, rspar1.shape,
                              rspar2.shape]).compound

# Start/stopping point for inboard ribs to avoid rib at kink
uref = fspar.cref.invert(rspar1.cref.p1)
builder1 = RibsAlongCurveByDistance('rib', fspar.cref, 25., fspar.shape,
                                    aft_shape, rhs_wing, d1=25., d2=-12.5,
                                    u2=uref)

uref_kink = fspar.cref.invert(rspar1.cref.p2)
builder2 = RibsAlongCurveByDistance('rib', fspar.cref, 25., fspar.shape,
                                    aft_shape, rhs_wing, d1=12.5, d2=-12.5,
                                    u1=uref, u2=uref_kink,
                                    first_index=builder1.next_index)

builder3 = RibsAlongCurveByDistance('rib', fspar.cref, 25., fspar.shape,
                                    aft_shape, rhs_wing, d1=12.5, d2=-12.,
                                    u1=uref_kink,
                                    first_index=builder2.next_index)

ribs = builder1.parts + builder2.parts + builder3.parts

# Fuse the ribs with their interfacing structure
FuseSurfaceParts(ribs, [root_rib, fspar, rspar1, rspar2])
DiscardByCref(ribs)

# Center wing ribs
ribs = RibsBetweenPlanesByNumber('center rib', xz_pln, root_rib.plane, 4,
                                 fc_spar.shape, rc_spar.shape, rhs_wing).parts
FuseSurfaceParts(ribs, [fc_spar, rc_spar])
DiscardByCref(ribs)

internal_parts = GroupAPI.get_parts()

# Skin
skin = SkinByBody('skin', rhs_wing).part
skin.fuse(*internal_parts)
skin.discard_by_dmin(rhs_wing.sref_shape, 1.)
skin.fix()
skin.set_transparency(0.5)

parts = GroupAPI.get_parts()

# Mesh
print('Computing mesh...')
the_mesh = MeshVehicle(4)

# Structured quad mesh if applicable
quad_hyp = QuadrangleHypo2D(the_mesh.gen)
quad_alg = QuadrangleAlgo2D(the_mesh.gen)
for part in internal_parts:
    the_mesh.set_quadrangle_2d(part.shape)
    # for face in part.shape.faces:
    #     if quad_alg.is_applicable(face, True):
    #         the_mesh.add_controls([quad_alg, quad_hyp], face)

the_mesh.compute()

# View
gui = Viewer()
gui.add(GroupAPI.get_master())
gui.start()
gui.clear()
gui.add(the_mesh)
gui.start()
