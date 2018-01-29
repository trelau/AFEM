from afem.config import Settings
from afem.exchange import ImportVSP
from afem.graphics import Viewer
from afem.smesh import *
from afem.structure import *
from afem.topology import *

Settings.log_to_console()

# v = Viewer()

# Set units to inch.
Settings.set_units('in')

# Import model
fn = r'../models/simple_wing.stp'
ImportVSP.step_file(fn)
wing = ImportVSP.get_body('WingGeom')

# Build structure
wingbox = AssemblyAPI.create_assy('wing box')
fspar = SparByParameters('front spar', 0.15, 0., 0.15, 1., wing).spar
rspar = SparByParameters('rear spar', 0.70, 0., 0.70, 1., wing).spar
RibByPoints('root rib', fspar.p1, rspar.p1, wing)
RibByPoints('tip rib', fspar.p2, rspar.p2, wing)
RibsAlongCurveByDistance('rib', rspar.cref, 30, fspar.shape, rspar.shape,
                         wing, d1=30, d2=-30)
internal_parts = wingbox.get_parts()
skin = SkinByBody('skin', wing).skin
cref = wing.sref.u_iso(0.5)
skin.discard_by_dmin(cref, 1.0)

FuseSurfaceParts([skin], internal_parts)

# Mesh
the_shape = wingbox.prepare_shape_to_mesh()
print('Creating mesh')
the_gen = MeshGen()
the_mesh = MeshGen.create_mesh(the_gen)
the_mesh.shape_to_mesh(the_shape)

# 1-d
hyp1d = LocalLength1D(the_gen, 4.)
alg1d = Regular1D(the_gen)
edges = ExploreShape.get_edges(the_shape)
cmp = CompoundByShapes(edges).compound
the_mesh.add_hypothesis(hyp1d, cmp)
the_mesh.add_hypothesis(alg1d, cmp)

# Unstructured quad-dominant
hyp2d = NetgenSimple2D(the_gen, 4.)
alg2d = NetgenAlgo2D(the_gen)
the_mesh.add_hypothesis(hyp2d, the_shape)
the_mesh.add_hypothesis(alg2d, the_shape)

# Apply mapped quadrangle to internal structure
# mapped_hyp = QuadrangleHypo2D(the_gen)
# mapped_alg = QuadrangleAlgo2D(the_gen)
# for part_ in internal_parts + [skin]:
#     for face in part_.faces:
#         if mapped_alg.is_applicable(face, True):
#             the_mesh.add_hypothesis(mapped_hyp, face)
#             the_mesh.add_hypothesis(mapped_alg, face)

the_gen.compute(the_mesh, the_shape)

# View
# skin.set_transparency(0.5)
v = Viewer()
v.add(wingbox)
v.start()
v.add(the_mesh)
v.start()
