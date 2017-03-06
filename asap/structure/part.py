from __future__ import print_function

from OCC.BRepCheck import BRepCheck_Analyzer
from OCC.ShapeAnalysis import ShapeAnalysis_ShapeTolerance
from OCC.ShapeFix import ShapeFix_Shape
from OCC.TopoDS import TopoDS_Shape

from .assembly import AssemblyData
from ..graphics.viewer import ViewableItem


class Part(TopoDS_Shape, ViewableItem):
    """
    Base class for all parts.
    """

    def __init__(self, name):
        super(Part, self).__init__()
        ViewableItem.__init__(self)
        self._name = name
        self._metadata = {}
        AssemblyData.add_parts(None, self)
        print('Creating part: ', name)

    @property
    def name(self):
        return self._name

    @property
    def tol(self):
        shp_tol = ShapeAnalysis_ShapeTolerance()
        shp_tol.AddTolerance(self)
        return shp_tol.GlobalTolerance(0)

    @property
    def metadata(self):
        return self._metadata

    def add_metadata(self, key, value):
        """
        Add metadata to the body.

        :param key:
        :param value:

        :return:
        """
        self._metadata[key] = value

    def get_metadata(self, key):
        """
        Get metadata.

        :param key:

        :return:
        """
        try:
            return self._metadata[key]
        except KeyError:
            return None

    def set_shape(self, shape):
        """
        Set the shape of the part.

        :param shape:

        :return:
        """
        if not isinstance(shape, TopoDS_Shape):
            return False

        self.TShape(shape.TShape())
        self.Location(shape.Location())
        self.Orientation(shape.Orientation())
        return True

    def check(self):
        """
        Check the shape of the part.

        :return:
        """
        check = BRepCheck_Analyzer(self, True)
        return check.IsValid()

    def fix(self):
        """
        Attempt to fix the shape of the part.

        :return:
        """
        fix = ShapeFix_Shape(self)
        fix.Perform()
        shape = fix.Shape()
        self.set_shape(shape)
        return self.check()
