from OCC.Display.SimpleGui import init_display

from asap.geometry import CreateGeom
from asap.utils.tcol import to_tcolgp_array1_pnt

display, start_display, add_menu, add_function_to_menu = init_display()

from OCC.AIS import AIS_Triangulation
from OCC.Poly import Poly_Triangle, Poly_Array1OfTriangle, \
    Poly_Triangulation
from OCC.Graphic3d import Graphic3d_AspectFillArea3d
from OCC.Prs3d import Prs3d_ShadingAspect
from OCC.Quantity import Quantity_Color

p1 = CreateGeom.point_by_xyz(1, 1, 0)
p2 = CreateGeom.point_by_xyz(0, 0, 0)
p3 = CreateGeom.point_by_xyz(2, 0, 0)
p4 = CreateGeom.point_by_xyz(1, -1, 0)

tcol_pnts = to_tcolgp_array1_pnt([p1, p2, p3, p4])

tri1 = Poly_Triangle(1, 2, 3)
tri2 = Poly_Triangle(2, 4, 3)

tri_arr = Poly_Array1OfTriangle(1, 2)
tri_arr.SetValue(1, tri1)
tri_arr.SetValue(2, tri2)

triangulation = Poly_Triangulation(tcol_pnts, tri_arr)

g3d = Graphic3d_AspectFillArea3d()
g3d.SetEdgeColor(Quantity_Color(0, 0, 0, 0))
g3d.SetInteriorColor(Quantity_Color(1, 0, 0, 0))
g3d.SetInteriorStyle(3)
g3d.SetEdgeOn()

prs = Prs3d_ShadingAspect()
prs.SetAspect(g3d.GetHandle())

ais_tri = AIS_Triangulation(triangulation.GetHandle())

ais_tri.Attributes().GetObject().SetShadingAspect(prs.GetHandle())

display.SetModeWireFrame()

display.Context.Display(ais_tri.GetHandle(), False)

display.FitAll()
display.Repaint()

start_display()
