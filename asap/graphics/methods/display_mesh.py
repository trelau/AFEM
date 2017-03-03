from OCC.MeshVS import MeshVS_BP_Mesh, MeshVS_DA_DisplayNodes, MeshVS_Mesh, \
    MeshVS_MeshPrsBuilder
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
    mesh_vs.GetDrawer().GetObject().SetBoolean(MeshVS_DA_DisplayNodes, False)
    context = display.Context
    context.Display(mesh_vs.GetHandle())
    context.Deactivate(mesh_vs.GetHandle())
    return True
