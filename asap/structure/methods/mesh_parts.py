from __future__ import print_function

from ...fem import MeshShape


def mesh_surface_part(part, maxh, quad_dominated=True):
    """
    Mesh a surface part.
    """
    print('Meshing part: ', part.name)
    for f in part.faces:
        MeshShape.perform(f, maxh, quad_dominated=quad_dominated)
    return True
