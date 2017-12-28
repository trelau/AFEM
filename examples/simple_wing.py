from afem.config import Settings
from afem.exchange import ImportVSP
from afem.graphics import Viewer
from afem.smesh import MeshGen, NetgenSimple2D, NetgenAlgo2D, MaxLength1D, \
    Regular1D
from afem.structure import *

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
RibsAlongCurveByDistance('rib', rspar.cref, 60, fspar, rspar, wing, d1=30,
                         d2=-30)
internal_parts = wingbox.get_parts()
skin = SkinByBody('skin', wing).skin
cref = wing.isocurve(u=0.5)
skin.discard_by_dmin(cref, 1.0)

FuseSurfaceParts([skin], internal_parts)

# Mesh
the_shape = wingbox.prepare_shape_to_mesh()
print('Creating mesh')
the_gen = MeshGen()
the_mesh = MeshGen.create_mesh(the_gen)
the_mesh.shape_to_mesh(the_shape)

# 1-d
hyp1d = MaxLength1D(the_gen, 4.)
alg1d = Regular1D(the_gen)
the_mesh.add_hypothesis(hyp1d, the_shape)
the_mesh.add_hypothesis(alg1d, the_shape)

# Unstructured quad-dominant
hyp2d = NetgenSimple2D(the_gen, 4.)
alg2d = NetgenAlgo2D(the_gen)
the_mesh.add_hypothesis(hyp2d, the_shape)
the_mesh.add_hypothesis(alg2d, the_shape)

# Apply mapped quadrangle to internal structure
# mapped_hyp = QuadrangleHypo2D(the_gen)
# mapped_alg = QuadrangleAlgo2D(the_gen)
# for part_ in internal_parts + [skin]:
#     if mapped_alg.is_applicable(part_, True):
#         the_mesh.add_hypothesis(mapped_hyp, part_)
#         the_mesh.add_hypothesis(mapped_alg, part_)
#     else:
#         print('Not applicable: {}'.format(part_.label))

the_gen.compute(the_mesh, the_shape)

# View
v = Viewer()
v.display_assy(wingbox)
v.start()
v.display_mesh(the_mesh.object)
v.start()
