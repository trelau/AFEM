from OCC.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCC.Display.SimpleGui import init_display
from OCC.Geom import Geom_Geometry
from OCC.MeshVS import (MeshVS_BP_Mesh, MeshVS_DA_DisplayNodes,
                        MeshVS_DA_EdgeColor, MeshVS_Drawer, MeshVS_Mesh,
                        MeshVS_MeshPrsBuilder)
from OCC.Quantity import Quantity_Color, Quantity_TOC_RGB
from OCC.SMESH import SMESH_MeshVSLink
from OCC.TopoDS import TopoDS_Shape
from OCC.gce import gce_MakeMirror
from numpy.random import rand

__all__ = ["Viewer", "ViewableItem"]


class ViewableItem(object):
    """
    Class for items that can be viewed.
    """

    def __init__(self):
        r, g, b = rand(1, 3)[0]
        self.color = Quantity_Color(r, g, b, Quantity_TOC_RGB)
        self.transparency = 0.
        self.mirror = None

    def set_color(self, r, g, b):
        """
        Set color (0. <= r, g, b <= 1.).

        :param float r: Red.
        :param float g: Green.
        :param float b: Blue.
        """
        if r > 1.:
            r /= 255.
        if g > 1.:
            g /= 255.
        if b > 1.:
            b /= 255.
        self.color = Quantity_Color(r, g, b, Quantity_TOC_RGB)

    def set_transparency(self, transparency):
        """
        Set the opacity for graphics.

        :param float transparency: Level of transparency (0 to 1).
        """
        if transparency < 0.:
            transparency = 0.
        elif transparency > 1.:
            transparency = 1.
        self.transparency = transparency

    def set_mirror(self, pln):
        """
        Set a plane to mirror the item.

        :param pln:
        """
        self.mirror = pln

    def get_mirrored(self):
        """
        Get the mirrored item.

        :return:
        """
        if not self.mirror:
            return None
        trsf = gce_MakeMirror(self.mirror.Pln()).Value()
        builder = BRepBuilderAPI_Transform(self, trsf, True)
        if not builder.IsDone():
            return None
        return builder.Shape()


class Viewer(object):
    """
    View objects using PythonOCC simple GUI.
    """
    _items = []
    _meshes = []

    @classmethod
    def clear(cls):
        """
        Clear entities from viewer.
        """
        cls._items = []
        cls._meshes = []
        cls._entities = []

    @classmethod
    def show(cls, clear=True, view='iso'):
        """
        Show the viewer.
        """
        disp, start_display, _, _ = init_display()
        for item in cls._items:
            disp.DisplayShape(item, color=item.color,
                              transparency=item.transparency)
            if item.mirror:
                mirrored = item.get_mirrored()
                disp.DisplayShape(mirrored, color=item.color,
                                  transparency=item.transparency)
        for mesh in cls._meshes:
            _display_smesh(disp, mesh)

        view = view.lower()
        if view in ['i', 'iso', 'isometric']:
            disp.View_Iso()
        elif view in ['t', 'top']:
            disp.View_Top()
        elif view in ['b', 'bottom']:
            disp.View_Bottom()
        elif view in ['l', 'left']:
            disp.View_Left()
        elif view in ['r', 'right']:
            disp.View_Right()
        elif view in ['f', 'front']:
            disp.View_Front()
        elif view in ['rear']:
            disp.View_Rear()

        disp.FitAll()
        disp.Repaint()
        start_display()

        if clear:
            cls.clear()

    @classmethod
    def add(cls, *items):
        """
        Add items to the viewer.

        :param items:

        :return:
        """
        for item in items:
            if isinstance(item, (ViewableItem, TopoDS_Shape, Geom_Geometry)):
                cls._items.append(item)

    @classmethod
    def add_meshes(cls, *meshes):
        """
        Add meshes to the viewer.

        :param meshes:
        :return:
        """
        for mesh in meshes:
            cls._meshes.append(mesh)

    @classmethod
    def capture(cls, filename='capture.png', clear=True, view='iso'):
        """
        Capture a screenshot from the viewer.
        """
        disp, start_display, _, _ = init_display()
        for item in cls._items:
            disp.DisplayShape(item, color=item.color,
                              transparency=item.transparency)
        for mesh in cls._meshes:
            _display_smesh(disp, mesh)

        view = view.lower()
        if view in ['i', 'iso', 'isometric']:
            disp.View_Iso()
        elif view in ['t', 'top']:
            disp.View_Top()
        elif view in ['b', 'bottom']:
            disp.View_Bottom()
        elif view in ['l', 'left']:
            disp.View_Left()
        elif view in ['r', 'right']:
            disp.View_Right()
        elif view in ['f', 'front']:
            disp.View_Front()
        elif view in ['rear']:
            disp.View_Rear()

        disp.FitAll()
        disp.Repaint()

        disp.View.Dump(filename)

        if clear:
            cls.clear()


def _display_smesh(display, mesh):
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
