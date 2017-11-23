from afem.config import Settings
from afem.fem import MeshAPI
from afem.graphics.display import display_shape
from afem.io import ImportVSP
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
RibsAlongCurveByDistance('rib', rspar.cref, 30, fspar, rspar, wing, d1=30,
                         d2=-30)
internal_parts = wingbox.get_parts()
skin = SkinByBody('skin', wing).skin
FuseSurfaceParts([skin], internal_parts)

# Mesh
the_shape = wingbox.prepare_shape_to_mesh()
print('Creating mesh')
the_mesh = MeshAPI.create_mesh('the mesh', the_shape)

# Unstructured quad-dominant
MeshAPI.hypotheses.create_netgen_simple_2d('netgen', 4.)
MeshAPI.hypotheses.create_netgen_algo_2d('netgen algo')
MeshAPI.add_hypothesis('netgen')
MeshAPI.add_hypothesis('netgen algo')

MeshAPI.hypotheses.create_max_length_1d('max length', 4.)
MeshAPI.hypotheses.create_regular_1d('algo 1d')
MeshAPI.add_hypothesis('max length')
MeshAPI.add_hypothesis('algo 1d')

# Apply mapped quadrangle to internal structure
mapped_hyp = MeshAPI.hypotheses.create_quadrangle_parameters('quad hyp')
mapped_algo = MeshAPI.hypotheses.create_quadrangle_aglo('quad algo')
for part_ in internal_parts + [skin]:
    if mapped_algo.is_applicable(part_, True):
        MeshAPI.add_hypothesis(mapped_hyp, part_)
        MeshAPI.add_hypothesis(mapped_algo, part_)
    else:
        print('Not applicable: {}'.format(part_.label))

MeshAPI.compute_mesh()

display_shape(None, the_mesh.handle)
