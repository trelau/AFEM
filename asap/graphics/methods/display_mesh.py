from OCC.AIS import AIS_Triangulation
from OCC.Graphic3d import Graphic3d_AspectFillArea3d
from OCC.Poly import Poly_Array1OfTriangle, Poly_Triangle, Poly_Triangulation
from OCC.Prs3d import Prs3d_ShadingAspect
from OCC.Quantity import Quantity_Color

from ...utils.tcol import to_tcolgp_array1_pnt


def display_part_mesh(display, part):
    """
    Display the mesh from a part.
    """
    from ...fem.data import MeshData
    for f in part.faces:
        mesh = MeshData.mesh_from_shape(f)
        if not mesh:
            continue
        display_trimesh(display, mesh, part.color)


def display_trimesh(display, mesh, color=None):
    """
    Display a triangular mesh.
    """
    # For now, just build arrays with duplicate nodes...
    elm_lst = mesh.elements
    pnt_lst = []
    tri_lst = []
    indx = 1
    for e in elm_lst:
        if not e.is_tri:
            continue
        for n in e.nodes:
            pnt_lst.append(n.pnt)
        tri_lst.append([indx, indx + 1, indx + 2])
        indx += 3

    if not tri_lst:
        return False

    # Create OCC data.
    tcol_pnts = to_tcolgp_array1_pnt(pnt_lst)
    ntri = len(tri_lst)
    arr_tri = Poly_Array1OfTriangle(1, ntri)
    indx = 1
    for i1, i2, i3 in tri_lst:
        tri = Poly_Triangle(i1, i2, i3)
        arr_tri.SetValue(indx, tri)
        indx += 1
    poly_tri = Poly_Triangulation(tcol_pnts, arr_tri)

    # Make graphics classes to render the elements and edges.
    gr = Graphic3d_AspectFillArea3d()
    gr.SetEdgeColor(Quantity_Color(0, 0, 0, 0))
    if color is None:
        gr.SetInteriorColor(mesh.color)
    else:
        gr.SetInteriorColor(color)
    gr.SetInteriorStyle(3)
    gr.SetEdgeOn()

    prs = Prs3d_ShadingAspect()
    prs.SetAspect(gr.GetHandle())

    ais_tri = AIS_Triangulation(poly_tri.GetHandle())
    ais_tri.Attributes().GetObject().SetShadingAspect(prs.GetHandle())

    display.SetModeWireFrame()
    display.Context.Display(ais_tri.GetHandle(), False)
    return True
