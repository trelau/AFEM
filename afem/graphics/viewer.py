import os

from OCCT.AIS import AIS_InteractiveContext, AIS_Shaded, AIS_Shape
from OCCT.Aspect import Aspect_DisplayConnection, Aspect_TOTP_LEFT_LOWER
from OCCT.BRepBuilderAPI import BRepBuilderAPI_Transform, \
    BRepBuilderAPI_MakeVertex, BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeFace
from OCCT.Graphic3d import Graphic3d_MaterialAspect
from OCCT.MeshVS import (MeshVS_DA_DisplayNodes, MeshVS_DA_EdgeColor,
                         MeshVS_Mesh, MeshVS_MeshPrsBuilder)
from OCCT.OpenGl import OpenGl_GraphicDriver
from OCCT.Quantity import (Quantity_TOC_RGB, Quantity_NOC_WHITE,
                           Quantity_Color,
                           Quantity_NOC_BLACK)
from OCCT.TopoDS import TopoDS_Shape
from OCCT.SMESH import SMESH_MeshVSLink
from OCCT.V3d import V3d_Viewer, V3d_TypeOfOrientation
from OCCT.WNT import WNT_Window
from OCCT.gce import gce_MakeMirror
from PySide import QtCore
from PySide.QtGui import QWidget, QApplication, QPalette, QIcon
from numpy.random import rand

__all__ = ["Viewer", "ViewableItem"]

# Icon location
_icon = os.path.dirname(__file__) + '/resources/main.png'


class ViewableItem(object):
    """
    Base class for items that can be viewed.

    :var OCCT.Quantity.Quantity_Color color: The OCC color quantity. The
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
        :rtype: OCCT.TopoDS.TopoDS_Shape
        """
        if not self.mirror:
            return None

        trsf = gce_MakeMirror(self.mirror.handle.Pln()).Value()
        builder = BRepBuilderAPI_Transform(self, trsf, True)
        if not builder.IsDone():
            return None
        return builder.Shape()


class Viewer(QWidget):
    """
    Simple class for viewing items.

    :param int width: Window width.
    :param int height: Window height.

    :var OCCT.V3d.V3d_Viewer my_viewer: The viewer.
    :var OCCT.V3d.V3d_View my_view: The view.
    :var OCCT.AIS.AIS_InteractiveContext: The context.
    :var OCCT.Prs3d.Prs3d_Drawer: The default drawer.
    """

    _app = QApplication.instance()
    if _app is None:
        _app = QApplication([])

    def __init__(self, width=800, height=600):
        super(Viewer, self).__init__()

        # Qt settings
        self.setBackgroundRole(QPalette.NoRole)
        self.setMouseTracking(True)
        self.resize(width, height)
        self.setWindowTitle('AFEM')
        icon = QIcon(_icon)
        self.setWindowIcon(icon)

        # Display connection
        display_connection = Aspect_DisplayConnection()

        # Graphics driver
        graphics_driver = OpenGl_GraphicDriver(display_connection)

        # Window handle
        window_handle = self.winId()

        # Windows window
        wind = WNT_Window(window_handle)

        # Create viewer and view
        self.my_viewer = V3d_Viewer(graphics_driver)
        self.my_view = self.my_viewer.CreateView()
        self.my_view.SetWindow(wind)

        # Map window
        wind.Map()

        # AIS interactive context
        self.my_context = AIS_InteractiveContext(self.my_viewer)

        # Some default settings
        self.my_context.SetAutomaticHilight(True)
        self._white = Quantity_Color(Quantity_NOC_WHITE)
        self._black = Quantity_Color(Quantity_NOC_BLACK)
        self.my_viewer.SetDefaultLights()
        self.my_viewer.SetLightOn()
        self.my_view.SetBackgroundColor(Quantity_TOC_RGB, 0.5, 0.5, 0.5)
        self.my_context.SetDisplayMode(AIS_Shaded, True)
        self.my_drawer = self.my_context.DefaultDrawer()
        self.my_drawer.SetFaceBoundaryDraw(True)
        self.my_view.TriedronDisplay(Aspect_TOTP_LEFT_LOWER, self._black, 0.08)

        # Keyboard map
        self._keys = {
            ord('F'): self.my_view.FitAll
        }

        # Values for mouse movement
        self._x0, self._y0 = 0., 0.

    def start(self, fit=True):
        """
        Start the application.

        :param bool fit: Option to fit the contents before starting.

        :return: None.
        """
        if fit:
            self.fit()
        self._app.exec_()

    def resizeEvent(self, *args, **kwargs):
        self.my_view.MustBeResized()

    def keyPressEvent(self, e):
        try:
            self._keys[e.key()]()
        except KeyError:
            pass

    def wheelEvent(self, e):
        if e.delta() > 0:
            zoom = 1.5
        else:
            zoom = 0.75
        self.my_view.SetZoom(zoom)

    def mousePressEvent(self, e):
        pos = e.pos()
        x, y = pos.x(), pos.y()
        self._x0, self._y0 = x, y
        self.my_view.StartRotation(x, y)

    def mouseDoubleClickEvent(self, e):
        pos = e.pos()
        x, y = pos.x(), pos.y()
        self._x0, self._y0 = x, y
        self.my_view.StartRotation(x, y)

    def mouseMoveEvent(self, e):
        pos = e.pos()
        x, y = pos.x(), pos.y()
        button = e.buttons()

        # Rotate
        if button == QtCore.Qt.LeftButton:
            self.my_view.Rotation(x, y)
        # Pan
        elif button in [QtCore.Qt.MidButton, QtCore.Qt.RightButton]:
            dx, dy = x - self._x0, y - self._y0
            self._x0, self._y0 = x, y
            self.my_view.Pan(dx, -dy)

    def display(self, ais_shape, update=True):
        """
        Display an AIS_Shape.

        :param OCCT.AIS.AIS_Shape ais_shape: The AIS shape.
        :param bool update: Option to update the viewer.

        :return: None.
        """
        self.my_context.Display(ais_shape, update)

    def display_body(self, body):
        """
        Display a body.

        :param afem.oml.entities.Body body: The body.

        :return: The AIS_Shape created for the body.
        :rtype: OCCT.AIS.AIS_Shape
        """
        return self.display_shape(body, body.color, body.transparency)

    def display_shape(self, shape, rgb=None, transparency=None, material=None):
        """
        Display a shape.

        :param OCCT.TopoDS.TopoDS_Shape shape: The shape.
        :param rgb: The RGB color (r, g, b).
        :type param: collections.Sequence[float] or OCCT.Quantity.Quantity_Color
        :param float transparency: The transparency (0 to 1).
        :param OCCT.Graphic3d.Graphic3d_NameOfMaterial material: The material.

        :return: The AIS_Shape created for the part.
        :rtype: OCCT.AIS.AIS_Shape
        """
        ais_shape = AIS_Shape(shape)

        if isinstance(rgb, (tuple, list)):
            r, g, b = rgb
            color = Quantity_Color(r, g, b, Quantity_TOC_RGB)
            ais_shape.SetColor(color)
        if isinstance(rgb, Quantity_Color):
            ais_shape.SetColor(rgb)

        if transparency is not None:
            ais_shape.SetTransparency(transparency)

        if material is not None:
            ma = Graphic3d_MaterialAspect(material)
            ais_shape.SetMaterial(ma)

        self.my_context.Display(ais_shape, True)
        return ais_shape

    def display_geom(self, geom):
        """
        Display a geometric entity.

        :param afem.geometry.entities.Geometry geom: The geometry. Must be
            either a Point, Curve, or Surface.

        :return: The AIS_Shape created for the geometry. Returns *None* if the
            entity cannot be converted to a shape.
        :rtype: OCCT.AIS.AIS_Shape or None
        """
        from afem.geometry.entities import Point, Curve, Surface

        if isinstance(geom, Point):
            shape = BRepBuilderAPI_MakeVertex(geom).Vertex()
        elif isinstance(geom, Curve):
            shape = BRepBuilderAPI_MakeEdge(geom.handle).Edge()
        elif isinstance(geom, Surface):
            shape = BRepBuilderAPI_MakeFace(geom.handle, 1.0e-7).Face()
        else:
            return None

        return self.display_shape(shape, geom.color, geom.transparency)

    def display_mesh(self, mesh, mode=1):
        """
        Display a mesh.

        :param OCCT.SMESH_SMESH_Mesh mesh: The mesh.
        :param int mode: Display mode for mesh elements (1=wireframe, 2=solid).

        :return: The MeshVS_Mesh created for the mesh.
        :rtype: OCCT.MeshVS.MeshVS_Mesh
        """
        vs_link = SMESH_MeshVSLink(mesh)
        mesh_vs = MeshVS_Mesh()
        mesh_vs.SetDataSource(vs_link)
        prs_builder = MeshVS_MeshPrsBuilder(mesh_vs)
        mesh_vs.AddBuilder(prs_builder)
        mesh_vs_drawer = mesh_vs.GetDrawer()
        mesh_vs_drawer.SetBoolean(MeshVS_DA_DisplayNodes, False)
        mesh_vs_drawer.SetColor(MeshVS_DA_EdgeColor, self._black)
        mesh_vs.SetDisplayMode(mode)
        self.my_context.Display(mesh_vs, True)
        return mesh_vs

    def display_part(self, part):
        """
        Display a part.

        :param afem.structure.entities.Part part: The part.

        :return: The AIS_Shape created for the part.
        :rtype: OCCT.AIS.AIS_Shape
        """
        return self.display_shape(part, part.color, part.transparency)

    def display_parts(self, parts):
        """
        Display a sequence of parts.

        :param collections.Sequence[afem.structure.entities.Part] parts: The
            parts.

        :return: None.
        """
        for part in parts:
            self.display_part(part)

    def display_assy(self, assy, include_subassy=True):
        """
        Display all parts of an assembly.

        :param afem.structure.assembly.Assembly assy: The assembly.
        :param bool include_subassy: Option to recursively include parts
            from any sub-assemblies.

        :return: None.
        """
        for part in assy.get_parts(include_subassy):
            self.display_part(part)

    def add(self, *items):
        """
        Add items to be displayed.

        :param items: The items.

        :return: None.
        """
        from afem.geometry.entities import Geometry
        from afem.oml.entities import Body
        from afem.structure.assembly import Assembly
        from afem.structure.entities import Part

        for item in items:
            if isinstance(item, Part):
                self.display_part(item)
            elif isinstance(item, Assembly):
                self.display_assy(item)
            elif isinstance(item, Geometry):
                self.display_geom(item)
            elif isinstance(item, Body):
                self.display_body(item)
            elif isinstance(item, TopoDS_Shape):
                self.display_shape(item)

    def fit(self):
        """
        Fit the contents.

        :return: None.
        """
        self.my_view.FitAll()

    def set_bg_color(self, r, g, b):
        """
        Set the background color.

        :param float r: The r-value.
        :param float g: The g-value.
        :param float b: The b-value.

        :return: None.
        """
        self.my_view.SetBackgroundColor(Quantity_TOC_RGB, r, g, b)

    def set_white_background(self):
        """
        Set the background color to white.

        :return: None.
        """
        self.my_view.SetBackgroundColor(Quantity_TOC_RGB, 1., 1., 1.)

    def view_iso(self):
        """
        Isometric view.

        :return: None.
        """
        self.my_view.SetProj(V3d_TypeOfOrientation.V3d_XposYnegZpos)

    def view_top(self):
        """
        Top view.

        :return: None.
        """
        self.my_view.SetProj(V3d_TypeOfOrientation.V3d_Zpos)

    def view_bottom(self):
        """
        Bottom view.

        :return: None.
        """
        self.my_view.SetProj(V3d_TypeOfOrientation.V3d_Zneg)

    def view_front(self):
        """
        Front view.

        :return: None.
        """
        self.my_view.SetProj(V3d_TypeOfOrientation.V3d_Xneg)

    def view_rear(self):
        """
        Rear view.

        :return: None.
        """
        self.my_view.SetProj(V3d_TypeOfOrientation.V3d_Xpos)

    def view_left(self):
        """
        Left view.

        :return: None.
        """
        self.my_view.SetProj(V3d_TypeOfOrientation.V3d_Yneg)

    def view_right(self):
        """
        Right view.

        :return: None.
        """
        self.my_view.SetProj(V3d_TypeOfOrientation.V3d_Ypos)

    def capture(self, fn):
        """
        Capture the screen contents and save to a file. The type of file will be
        determined by the extension.

        :param str fn: The filename.

        :return: None.
        """
        self.my_view.Dump(fn)

    def remove_all(self):
        """
        Remove all items from the context.

        :return: None.
        """
        self.my_context.RemoveAll(True)
