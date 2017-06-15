from OCC.Geom import Geom_Curve, Geom_Surface
from OCC.TopoDS import TopoDS_Shape

from .assembly import AssemblyData
from .methods.cut_parts import cut_part
from .methods.reshape_parts import reshape_parts
from .methods.split_parts import split_part
from ..geometry import ProjectGeom
from ..graphics.viewer import ViewableItem
from ..topology import ShapeTools


class Part(TopoDS_Shape, ViewableItem):
    """
    Part base.
    """
    _indx = 1
    _all = {}

    def __init__(self, label, shape, cref=None, sref=None,
                 add_to_assy=True):
        # Initialize
        super(Part, self).__init__()
        ViewableItem.__init__(self)
        self._cref, self._sref = None, None
        self._metadata = {}
        self._subparts = {}

        # Set ID
        self._id = Part._indx
        Part._all[self._id] = self
        Part._indx += 1

        # Set attributes
        self._label = label
        shape = ShapeTools.to_shape(shape)
        assert isinstance(shape, TopoDS_Shape), "Invalid shape used in Part."
        self.set_shape(shape)
        self.set_cref(cref)
        self.set_sref(sref)

        # Store in active assembly
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

    @property
    def cref(self):
        return self._cref

    @property
    def sref(self):
        return self._sref

    @property
    def has_cref(self):
        return isinstance(self._cref, Geom_Curve)

    @property
    def has_sref(self):
        return isinstance(self._sref, Geom_Surface)

    @property
    def p1(self):
        try:
            return self._cref.eval(self._cref.u1)
        except AttributeError:
            return None

    @property
    def p2(self):
        try:
            return self._cref.eval(self._cref.u2)
        except AttributeError:
            return None

    @property
    def edges(self):
        return ShapeTools.get_edges(self)

    @property
    def nedges(self):
        return len(self.edges)

    @property
    def faces(self):
        return ShapeTools.get_faces(self)

    @property
    def nfaces(self):
        return len(self.faces)

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

    def add_subpart(self, label, subpart):
        """
        Add a subpart to the part.
        
        :param label: 
        :param subpart:
         
        :return: 
        """
        if not isinstance(subpart, Part):
            return False
        try:
            self._subparts[label] = subpart
            return True
        except KeyError:
            return False

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

    def reshape(self, tool):
        """
        Reshape the part with a tool.
        
        :param tool: 
        :return: 
        """
        return reshape_parts(tool, [self])

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

    def set_cref(self, cref):
        """
        Set the part reference curve.

        :param cref:

        :return: 
        """
        if not isinstance(cref, Geom_Curve):
            return False
        self._cref = cref
        return True

    def set_sref(self, sref):
        """
        Set the part reference surface.
        
        :param sref:
         
        :return: 
        """
        if not isinstance(sref, Geom_Surface):
            return False
        self._sref = sref
        return True

    def local_to_global_u(self, u):
        """
        Convert local parameter from 0 <= u <= 1 to u1 <= u <= u2.

        :param u:

        :return:
        """
        try:
            return self._cref.local_to_global_param(u)
        except AttributeError:
            return None

    def eval_cref(self, u):
        """
        Evaluate point on reference curve.

        :param u:
        :return:
        """
        try:
            return self._cref.eval(u)
        except AttributeError:
            return None

    def eval_sref(self, u, v):
        """
        Evaluate point on reference surface.

        :param u:
        :param v:
        :return:
        """
        try:
            return self._sref.eval(u, v)
        except AttributeError:
            return None

    def eval_dx(self, dx, u0=None, is_local=False):
        """
        Evaluate point on reference curve at a distance from a parameter.

        :param float dx:
        :param u0:
        :param bool is_local:
        
        :return:
        """
        try:
            return self._cref.eval_dx(dx, u0, is_local)
        except AttributeError:
            return None

    def spaced_points(self, dx, s1=None, s2=None, u1=None, u2=None,
                      shape1=None, shape2=None):
        """
        Create points along the reference curve.

        :param dx:
        :param float s1:
        :param float s2:
        :param float u1:
        :param float u2:
        :param shape1:
        :param shape2:

        :return:
        """
        if not self.has_cref:
            return []

        if isinstance(dx, int):
            npts = dx
            dx = None
        else:
            npts = None
        return ShapeTools.points_along_shape(self.cref, dx, npts, u1, u2,
                                             s1, s2, shape1, shape2)

    def point_to_cref(self, pnt, direction=None):
        """
        Project a point to reference curve.
        """
        if not self.has_cref:
            return False

        ProjectGeom.point_to_geom(pnt, self.cref, True, direction)
        return True

    def points_to_cref(self, pnts, direction=None):
        """
        Project points to reference curve.
        """
        if not self.has_cref:
            return False

        for p in pnts:
            ProjectGeom.point_to_geom(p, self.cref, True, direction)
        return True

    def point_to_sref(self, pnt, direction=None):
        """
        Project a point to reference surface.
        """
        if not self.has_sref:
            return False

        ProjectGeom.point_to_geom(pnt, self.sref, True, direction)
        return True

    def points_to_sref(self, pnts, direction=None):
        """
        Project points to reference surface.
        """
        if not self.has_sref:
            return False

        for p in pnts:
            ProjectGeom.point_to_geom(p, self.sref, True, direction)
        return True

    def invert_cref(self, pnt):
        """
        Invert the point on the reference curve.

        :param pnt:

        :return:
        """
        if not self.has_cref:
            return False
        return ProjectGeom.invert(pnt, self.cref)

    def invert_sref(self, pnt):
        """
        Invert the point on the reference surface.

        :param pnt:

        :return:
        """
        if not self.has_sref:
            return False
        return ProjectGeom.invert(pnt, self.sref)

    def distance(self, other):
        """
        Find the minimum distance between the part and other shape.
        
        :param other: Other part or shape.
        
        :return: The minimum distance and the points of minimum distance on 
            each shape (dmin, p1, p2).
        :rtype: tuple
        """
        return ShapeTools.min_distance(self, other)
