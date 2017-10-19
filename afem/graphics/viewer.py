from OCC.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCC.Display.SimpleGui import init_display
from OCC.Graphic3d import Graphic3d_NOM_DEFAULT
from OCC.V3d import V3d_ColorScale
from OCC.MeshVS import (MeshVS_BP_Mesh, MeshVS_DA_DisplayNodes,
                        MeshVS_DA_EdgeColor, MeshVS_Drawer, MeshVS_Mesh,
                        MeshVS_MeshPrsBuilder)
from OCC.Quantity import Quantity_Color, Quantity_TOC_RGB
from OCC.SMESH import SMESH_MeshVSLink
from OCC.TopoDS import TopoDS_Shape
from OCC.gce import gce_MakeMirror
from numpy.random import rand
import ctypes
import gc

__all__ = ["Viewer", "ViewableItem"]


class ViewableItem(object):
    """
    Base class for items that can be viewed.

    :var OCC.Quantity.Quantity_Color color: The OCC color quantity. The
        color is set randomly during initialization.
    :var float transparency: The transparency level.
    :var afem.geometry.entities.Plane mirror: The plane to mirror the object
        about. If provided then object will be mirrored about the plane for
        visualization purposes only.
    """

    def __init__(self):
        r, g, b = rand(1, 3)[0]
        self.color = Quantity_Color(r, g, b, Quantity_TOC_RGB)
        self.transparency = 0.
        self.mirror = None
        return

    def set_color(self, r, g, b):
        """
        Set color (0. <= r, g, b <= 1.).

        :param float r: Red.
        :param float g: Green.
        :param float b: Blue.

        :return: None.
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

        :return: None.
        """
        if transparency < 0.:
            transparency = 0.
        elif transparency > 1.:
            transparency = 1.
        self.transparency = transparency

    def set_mirror(self, pln):
        """
        Set a plane to mirror the item.

        :param afem.geometry.entities.Plane pln: The plane.

        :return: None.
        """
        self.mirror = pln

    def get_mirrored(self):
        """
        Get the mirrored shape.

        :return: The mirrored shape.
        :rtype: OCC.TopoDS.TopoDS_Shape
        """
        if not self.mirror:
            return None

        trsf = gce_MakeMirror(self.mirror.object.Pln()).Value()
        builder = BRepBuilderAPI_Transform(self, trsf, True)
        if not builder.IsDone():
            return None
        return builder.Shape()


class Viewer:
    def __init__(
            self, shapes=list(), meshes=list(), view='iso', white_bg=False, fit=True, xy1=None, xy2=None,
            draw_edges=True
    ):
        """
        Initializes a viewer object

        :param shapes: Item(s) to view.
        :type shapes: list of afem.graphics.viewer.ViewableItem or TopoDS.TopoDS_Shape
        :param list of afem.fem.meshes.Mesh meshes: Mesh(es) to view.
        :param str view: The viewing angle ('iso', 'top', 'bottom', 'left',
            'right', 'front', 'bottom').
        :param bool white_bg: Option to make the background white.
        :param bool fit: Option to fit contents to screen.
        :param array_like xy1: Lower corner for zooming.
        :param array_like xy2: Upper corner for zooming.
        :param bool draw_edges: Display edges as black lines if True
        """
        self._items = shapes
        self._meshes = meshes
        self._draw_edges = draw_edges
        view_scale = 0.9
        screen_res = (int(view_scale * ctypes.windll.user32.GetSystemMetrics(0)),
                      int(view_scale * ctypes.windll.user32.GetSystemMetrics(1)))
        self._disp, self._start_display, add_menu, add_function_to_menu = init_display(size=screen_res)
        self._win = add_menu.__closure__[0].cell_contents
        self.set_display_shapes()
        self.change_view(view)

        # self._disp.View.ColorScale()
        # self._disp.View.ColorScaleDisplay()
        self.draw_edges = draw_edges
        self._disp.default_drawer.SetFaceBoundaryDraw(self.draw_edges)

        if white_bg:
            self._disp.set_bg_gradient_color(255, 255, 255, 255, 255, 255)
            self._disp.Repaint()

        if xy1 is not None and xy2 is not None:
            self._disp.View.FitAll(xy1[0], xy1[1], xy2[0], xy2[1])
            fit = False

        if fit:
            self._disp.FitAll()

        self._disp.Repaint()
        self._win.close()
        return

    @property
    def draw_edges(self):
        return self._draw_edges

    @draw_edges.setter
    def draw_edges(self, tf):
        self._draw_edges = tf

    def clear_all(self):
        """
        Clear entities from viewer and from object memory
        :return: None.
        """
        self._items = []
        self._meshes = []
        self._disp.EraseAll()
        return

    def close(self):
        """
        Close the viewer window
        :return: None
        """
        self._win.close()
        return

    def delete_viewer(self):
        """
        Clear the viewer canvas object from memory
        :return: None
        """
        self._win.canva.close()
        self._win.canva.deleteLater()
        self._win.menu_bar.close()
        self._win.menu_bar.deleteLater()
        self._win = None
        self._disp = None
        self._start_display = None
        gc.collect()
        return

    def change_view(self, v):
        """
        Change the viewing angle of the current viewer.
        :param str v: The viewing angle ('iso', 'top', 'bottom', 'left',
            'right', 'front', 'bottom').
        :return:
        """
        v = v.lower()
        if v in ['i', 'iso', 'isometric']:
            self._disp.View_Iso()
        elif v in ['t', 'top']:
            self._disp.View_Top()
        elif v in ['b', 'bottom']:
            self._disp.View_Bottom()
        elif v in ['l', 'left']:
            self._disp.View_Front()
        elif v in ['r', 'right']:
            self._disp.View_Rear()
        elif v in ['f', 'front']:
            self._disp.View_Left()
        elif v in ['rear']:
            self._disp.View_Right()
        else:
            raise ValueError
        return

    def show(self):
        """
        Display the current viewer for interactive use
        :return: None
        """
        self._disp.Repaint()
        self._win.show()
        self._start_display()
        return

    def add(self, *items):
        """
        Add viewable items.

        :param items: Item(s) to view.
        :type items: afem.graphics.viewer.ViewableItem or TopoDS.TopoDS_Shape

        :return: None.
        """
        for item in items:
            if isinstance(item, (ViewableItem, TopoDS_Shape)):
                self._items.append(item)
        return

    def set_display_shapes(self):
        """
        Draws but does not display the current _items list and _meshes list in the viewer
        :return: None
        """
        for item in self._items:
            try:
                self._disp.DisplayShape(
                    item, color=item.color,
                    transparency=item.transparency,
                    material=Graphic3d_NOM_DEFAULT
                )
            except TypeError:
                # Hack for some geometry items
                self._disp.DisplayShape(
                    item.object,
                    color=item.color,
                    transparency=item.transparency,
                    material=Graphic3d_NOM_DEFAULT
                )
            if item.mirror:
                mirrored = item.get_mirrored()
                self._disp.DisplayShape(
                    mirrored,
                    color=item.color,
                    transparency=item.transparency,
                    material=Graphic3d_NOM_DEFAULT
                )
        for mesh in self._meshes:
            _display_smesh(self._disp, mesh)
        self._disp.Repaint()
        return

    def add_meshes(self, *meshes):
        """
        Add meshes to the viewer.

        :param afem.fem.meshes.Mesh meshes: The mesh(es) to add.

        :return: None.
        """
        for mesh in meshes:
            self._meshes.append(mesh)
        return

    def capture(self, filename='capture.png', view='iso', fit=True):
        """
        Take a screen capture of the current viewer configuration
        :param str filename: name of the file to be captured
        :param str view: The viewing angle ('iso', 'top', 'bottom', 'left',
            'right', 'front', 'bottom').
        :param bool fit: Option to fit contents to screen.
        :return: None
        """
        self._win.show()
        self.change_view(view)
        if fit:
            self._disp.FitAll()
        self._disp.Repaint()
        self._disp.View.Dump(filename)
        self._win.close()
        return


def _display_smesh(display, mesh):
    """
    Display an SMESH generated mesh.
    """
    ds = SMESH_MeshVSLink(mesh.object)
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
