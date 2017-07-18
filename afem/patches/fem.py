from OCC.SMESH import SMESH_Mesh
from OCC.TopoDS import TopoDS_Edge, TopoDS_Face

from ..fem.elements import Elm1D, Elm2D
from ..fem.meshes import MeshData


def _elements_on_edge(self):
    """
    Get the elements on an edge in the active mesh.
    """
    mesh = MeshData.get_active()
    if not mesh:
        return []
    mesh = mesh.smesh_obj
    if not isinstance(mesh, SMESH_Mesh):
        return []

    submesh = mesh.GetSubMesh(self)
    if submesh.IsEmpty():
        return []

    submesh_ds = submesh.GetSubMeshDS()
    if not submesh_ds:
        return []

    elm_iter = submesh_ds.GetElements()
    elm_set = set()
    while elm_iter.more():
        elm = Elm1D(elm_iter.next())
        elm_set.add(elm)
    return elm_set


def _elements_on_face(self):
    """
    Get the elements on a face in the active mesh.
    """
    mesh = MeshData.get_active()
    if not mesh:
        return []
    mesh = mesh.smesh_obj
    if not isinstance(mesh, SMESH_Mesh):
        return []

    submesh = mesh.GetSubMesh(self)
    if submesh.IsEmpty():
        return []

    submesh_ds = submesh.GetSubMeshDS()
    if not submesh_ds:
        return []

    elm_iter = submesh_ds.GetElements()
    elm_set = set()
    while elm_iter.more():
        elm = Elm2D(elm_iter.next())
        elm_set.add(elm)
    return elm_set


TopoDS_Edge.elements = property(_elements_on_edge)
TopoDS_Face.elements = property(_elements_on_face)


# TODO Assign hypotheses
