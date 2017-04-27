from OCC.TopoDS import TopoDS_Shape

from .assembly import AssemblyData
from .methods.cut_parts import cut_part
from .methods.split_parts import split_part
from ..graphics.viewer import ViewableItem
from ..topology import ShapeTools


class Part(TopoDS_Shape, ViewableItem):
    """
    Base class for all parts.
    """
    _indx = 1
    _all = {}

    def __init__(self, label, add_to_assy=True):
        super(Part, self).__init__()
        ViewableItem.__init__(self)
        self._label = label
        self._metadata = {}
        self._subparts = {}
        self._id = Part._indx
        Part._all[self._id] = self
        Part._indx += 1
        # Store in active assembly.
        if add_to_assy:
            AssemblyData.add_parts(None, self)
        print('Creating part: ', label)

    @property
    def is_null(self):
        return self.IsNull()

    @property
    def label(self):
        return self._label

    @property
    def id(self):
        return self._id

    @property
    def tol(self):
        return ShapeTools.get_tolerance(self, 0)

    @property
    def max_tol(self):
        return ShapeTools.get_tolerance(self, 1)

    @property
    def min_tol(self):
        return ShapeTools.get_tolerance(self, 2)

    @property
    def metadata(self):
        return self._metadata

    @property
    def subparts(self):
        return self._subparts.values()

    def nullify(self):
        """
        Destroy reference to underlying shape.
        
        :return: 
        """
        self.Nullify()
        return True

    def add_metadata(self, key, value):
        """
        Add metadata to the part.

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

    def get_subpart(self, label):
        """
        Get a sub-part.
        
        :param label:
         
        :return: 
        """
        try:
            return self._subparts[label]
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
        return ShapeTools.is_valid(self)

    def fix(self):
        """
        Attempt to fix the shape of the part.

        :return:
        """
        new_shape = ShapeTools.fix_shape(self)
        if not new_shape:
            return False
        return self.set_shape(new_shape)

    def cut(self, cutter):
        """
        Cut the part with a shape.

        :param cutter:

        :return:
        """
        cutter = ShapeTools.to_shape(cutter)
        if not cutter:
            return False
        return cut_part(self, cutter)

    def split(self, splitter, split_both=True):
        """
        Split the part with another part or shape.
        
        :param splitter: 
        :param split_both:
         
        :return: 
        """
        return split_part(self, splitter, split_both)
