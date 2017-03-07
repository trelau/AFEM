from OCC.MeshVS import MeshVS_BP_Mesh, MeshVS_DA_DisplayNodes, \
    MeshVS_DA_EdgeColor, MeshVS_Drawer, MeshVS_Mesh, MeshVS_MeshPrsBuilder
from OCC.Quantity import Quantity_Color
from OCC.SMESH import SMESH_MeshVSLink


def display_smesh(display, mesh):
    """
    Display an SMESH generated mesh.
    """
    ds = SMESH_MeshVSLink(mesh.smesh_obj)
    mesh_vs = MeshVS_Mesh(True)
    prs_builder = MeshVS_MeshPrsBuilder(mesh_vs.GetHandle(), 1,
                                        ds.GetHandle(), 0, MeshVS_BP_Mesh)
    mesh_vs.SetDataSource(ds.GetHandle())
    mesh_vs.AddBuilder(prs_builder.GetHandle(), True)
    mesh_vs_drawer = mesh_vs.GetDrawer().GetObject()
    assert isinstance(mesh_vs_drawer, MeshVS_Drawer)
    mesh_vs_drawer.SetBoolean(MeshVS_DA_DisplayNodes, False)
    black = Quantity_Color(0., 0., 0., 0)
    mesh_vs_drawer.SetColor(MeshVS_DA_EdgeColor, black)
    context = display.Context
    context.Display(mesh_vs.GetHandle())
    context.Deactivate(mesh_vs.GetHandle())
    return True
