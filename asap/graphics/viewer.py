from OCC.Display.SimpleGui import init_display
from OCC.Quantity import Quantity_Color, Quantity_TOC_RGB
from numpy.random import rand

from .methods.display_mesh import display_smesh


class ViewableItem(object):
    """
    Class for items that can be viewed.
    """

    def __init__(self):
        r, g, b = rand(1, 3)[0]
        self._color = Quantity_Color(r, g, b, Quantity_TOC_RGB)
        self._transparency = 0.

    @property
    def color(self):
        return self._color

    def set_color(self, r, g, b):
        """
        Set color (0. <= r, g, b <= 1.).

        :param float r: Red.
        :param float g: Green.
        :param float b: Blue.
        """
        self._color = Quantity_Color(r, g, b, Quantity_TOC_RGB)

    @property
    def transparency(self):
        return self._transparency

    def set_transparency(self, transparency):
        """
        Set the opacity for graphics.

        :param float transparency: Level of transparency (0 to 1).
        """
        if transparency < 0.:
            transparency = 0.
        elif transparency > 1.:
            transparency = 1.
        self._transparency = transparency


class Viewer(object):
    """
    View objects using PythonOCC simple GUI.
    """
    _items = []
    _meshes = []
    _entities = []

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
        for entity, color, transparency in cls._entities:
            disp.DisplayShape(entity, color=color, transparency=transparency)
        for mesh in cls._meshes:
            display_smesh(disp, mesh)

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
    def add_items(cls, *items):
        """
        Add viewable items to the viewer.

        :param items:

        :return:
        """
        for item in items:
            if isinstance(item, ViewableItem):
                cls._items.append(item)

    @classmethod
    def add_entity(cls, entity, color=None, transparency=0.):
        """
        Add an entity to the viewer.
        """
        if isinstance(color, (tuple, list)):
            color = Quantity_Color(color[0], color[1], color[2],
                                   Quantity_TOC_RGB)
        elif color in ['rand', 'random']:
            r, g, b = rand(1, 3)[0]
            color = Quantity_Color(r, g, b, Quantity_TOC_RGB)
        if not isinstance(color, Quantity_Color):
            color = None
        cls._entities.append([entity, color, transparency])

    @classmethod
    def add(cls, *items):
        """
        Add items to the viewer.

        :param items:

        :return:
        """
        for item in items:
            if isinstance(item, ViewableItem):
                cls.add_items(item)
            else:
                cls.add_entity(item)

    @classmethod
    def add_meshes(cls, *meshes):
        """
        Add meshes to the viewer.

        :param meshes:
        :return:
        """
        for mesh in meshes:
            cls._meshes.append(mesh)
