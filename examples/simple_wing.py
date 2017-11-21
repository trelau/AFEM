from afem.config import Settings
from afem.fem import MeshAPI
# from afem.graphics import Viewer
from afem.io import ImportVSP
from afem.structure import *

Settings.log_to_console(True)

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
RibsAlongCurveByDistance('rib', rspar.cref, 762, fspar, rspar, wing, d1=762, d2=-762)
internal_parts = wingbox.get_parts()
skin = SkinByBody('skin', wing).skin
FuseSurfaceParts([skin], internal_parts)

skin.set_transparency(0.75)
# v.add('wing box')
# v.set_display_shapes()
# v.show()

# Mesh
the_shape = wingbox.prepare_shape_to_mesh()
the_mesh = MeshAPI.create_mesh('the mesh', the_shape)

# Unstructured quad-dominant
MeshAPI.hypotheses.create_netgen_simple_2d('netgen', 101.6)
MeshAPI.hypotheses.create_netgen_algo_2d('netgen algo')
MeshAPI.add_hypothesis('netgen')
MeshAPI.add_hypothesis('netgen algo')
MeshAPI.compute_mesh()

# v.add_meshes(MeshAPI.get_active())
# v.set_display_shapes()
# v.show()

from afem.graphics.display import display_shape

display_shape(None, the_mesh.handle)
