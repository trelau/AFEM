from OCC.Geom import Geom_Curve, Geom_Surface
from OCC.SMESH import SMESH_Mesh
from OCC.TopAbs import TopAbs_COMPSOLID, TopAbs_SOLID
from OCC.TopoDS import TopoDS_Shape

from afem.fem.elements import Elm2D
from afem.fem.meshes import MeshData
from afem.geometry import CreateGeom, ProjectGeom
from afem.graphics.viewer import ViewableItem
from afem.structure.assembly import AssemblyData
from afem.structure.methods.cut_parts import cut_part, \
    cut_wing_part_with_circle
from afem.structure.methods.explore_parts import get_shared_edges, \
    get_shared_nodes
from afem.structure.methods.fuse_parts import fuse_surface_part
from afem.structure.methods.merge_parts import merge_surface_part
from afem.structure.methods.modify_parts import discard_faces_by_distance, \
    discard_faces_by_solid, discard_wing_part_faces, unify_surface_part
from afem.structure.methods.reshape_parts import reshape_parts
from afem.structure.methods.sew_parts import sew_surface_parts
from afem.structure.methods.split_parts import split_part
from afem.topology import ShapeTools

__all__ = ["Part", "CurvePart", "Beam", "SurfacePart", "WingPart", "Spar",
           "Rib", "FuselagePart", "Bulkhead", "Floor", "Frame", "Skin",
           "Stiffener1D", "Stiffener2D", "Stringer"]


class Part(TopoDS_Shape, ViewableItem):
    """
    Part base.
    """
    _indx = 1
    _all = {}

    def __init__(self, label, shape, cref=None, sref=None, assy=None):
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

        # Add to assembly
        assy = AssemblyData.get_assy(assy)
        if assy:
            AssemblyData.add_parts(assy, self)
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


class CurvePart(Part):
    """
    Base class for curve-based parts.
    """

    def __init__(self, label, shape, cref=None, assy=None):
        super(CurvePart, self).__init__(label, shape, cref, None, assy)

    @property
    def reshapes(self):
        return self.edges

    @property
    def length(self):
        return ShapeTools.shape_length(self)


class Beam(CurvePart):
    """
    Beam.
    """

    def __init__(self, label, shape, cref=None, assy=None):
        super(Beam, self).__init__(label, shape, cref, assy)


class SurfacePart(Part):
    """
    Base class for surface-based parts.
    """

    def __init__(self, label, shape, cref=None, sref=None, assy=None):
        super(SurfacePart, self).__init__(label, shape, cref, sref, assy)

    @property
    def reshapes(self):
        return self.faces

    @property
    def area(self):
        return ShapeTools.shape_area(self)

    @property
    def stiffeners(self):
        return [part for part in self.subparts]

    @property
    def elements(self):
        smesh_mesh = MeshData.get_mesh().smesh_obj
        if not isinstance(smesh_mesh, SMESH_Mesh):
            return []
        compound = ShapeTools.get_faces(self, True)
        submesh = smesh_mesh.GetSubMesh(compound)
        if submesh.IsEmpty():
            return []
        submesh_ds = submesh.GetSubMeshDS()
        if not submesh_ds:
            return []
        elm_iter = submesh_ds.GetElements()
        elm_set = set()
        while elm_iter.more():
            elm = Elm2D(elm_iter.next())
            elm_set.add(elm)
        return elm_set

    @property
    def nodes(self):
        node_set = set()
        for e in self.elements:
            for n in e.nodes:
                node_set.add(n)
        return node_set

    def fuse(self, *other_parts):
        """
        Fuse with other parts.

        :param other_parts:

        :return:
        """
        _other_parts = []
        for part in other_parts:
            if isinstance(part, TopoDS_Shape) and not part.IsNull():
                _other_parts.append(part)
        if not _other_parts:
            return False
        return fuse_surface_part(self, *_other_parts)

    def sew(self, *other_parts):
        """
        Sew with other parts.

        :param other_parts:

        :return:
        """
        _other_parts = []
        for part in other_parts:
            if isinstance(part, TopoDS_Shape) and not part.IsNull():
                _other_parts.append(part)
        if not _other_parts:
            return False
        return sew_surface_parts([self] + _other_parts)

    def merge(self, other, unify=False):
        """
        Merge other part or shape with this one.

        :param other:
        :param bool unify:

        :return:
        """
        return merge_surface_part(self, other, unify)

    def unify(self, edges=True, faces=True, concat_bsplines=False):
        """
        Attempt to unify the same domains of the part shape.

        :param edges:
        :param faces:
        :param concat_bsplines:

        :return:
        """
        return unify_surface_part(self, edges, faces, concat_bsplines)

    def discard(self, shape, tol=None):
        """
        Discard faces of the part.

        :param shape:
        :param tol:

        :return:
        """
        shape = ShapeTools.to_shape(shape)
        if not shape:
            return False
        if shape.ShapeType() in [TopAbs_SOLID, TopAbs_COMPSOLID]:
            return discard_faces_by_solid(self, shape, tol)
        return discard_faces_by_distance(self, shape, tol)

    def shared_edges(self, other_part):
        """
        Get edges shared between the two parts.

        :param other_part:

        :return:
        """
        return get_shared_edges(self, other_part)

    def shared_nodes(self, other_part):
        """
        Get nodes shared between the two parts.

        :param other_part:

        :return:
        """
        return get_shared_nodes(self, other_part)


class WingPart(SurfacePart):
    """
    Base class for wing parts.
    """

    def __init__(self, label, shape, cref=None, sref=None, assy=None):
        super(WingPart, self).__init__(label, shape, cref, sref, assy)

    def discard(self, shape=None, tol=None):
        """
        Discard faces of the part that are inside the solid or use automated
        if *solid* is *None*.

        :param shape:
        :param tol:

        :return:
        """
        if not shape:
            # Automatic method for wing parts.
            return discard_wing_part_faces(self)
        else:
            # Call original method.
            return super(WingPart, self).discard(shape, tol)

    def cut_hole(self, dx, r):
        """
        Cut a circular hole in the part (in development).

        :param dx:
        :param r:

        :return:
        """
        return cut_wing_part_with_circle(self, dx, r)

    def get_plane(self, u=None, dx=None, pnt=None, ref_pln=None):
        """
        Get a plane along the reference curve.

        :param u:
        :param dx:
        :param pnt:
        :param ref_pln:

        :return:
        """
        return CreateGeom.plane_on_curve(self.cref, u, dx, pnt, ref_pln)


class Spar(WingPart):
    """
    Wing spar.
    """

    def __init__(self, label, shape, cref=None, sref=None, assy=None):
        super(Spar, self).__init__(label, shape, cref, sref, assy)


class Rib(WingPart):
    """
    Wing rib.
    """

    def __init__(self, label, shape, cref=None, sref=None, assy=None):
        super(Rib, self).__init__(label, shape, cref, sref, assy)


class FuselagePart(SurfacePart):
    """
    Base class for fuselage parts.
    """

    def __init__(self, label, shape, cref=None, sref=None, assy=None):
        super(FuselagePart, self).__init__(label, shape, cref, sref, assy)


class Bulkhead(FuselagePart):
    """
    Bulkhead.
    """

    def __init__(self, label, shape, cref=None, sref=None, assy=None):
        super(Bulkhead, self).__init__(label, shape, cref, sref, assy)


class Floor(FuselagePart):
    """
    Floor.
    """

    def __init__(self, label, shape, cref=None, sref=None, assy=None):
        super(Floor, self).__init__(label, shape, cref, sref, assy)


class Frame(FuselagePart):
    """
    Frame.
    """

    def __init__(self, name, shape, cref=None, sref=None, assy=None):
        super(Frame, self).__init__(name, shape, cref, sref, assy)


class Skin(SurfacePart):
    """
    Skin part.
    """

    def __init__(self, name, shape, cref=None, sref=None, assy=None):
        super(Skin, self).__init__(name, shape, cref, sref, assy)


class Stiffener1D(CurvePart):
    """
    1-D stiffener for surface parts.
    """

    def __init__(self, label, shape, cref=None, assy=False):
        super(Stiffener1D, self).__init__(label, shape, cref, assy)


class Stiffener2D(SurfacePart):
    """
    2-D stiffener for surface parts.
    """

    def __init__(self, label, shape, cref=None, sref=None, assy=False):
        super(Stiffener2D, self).__init__(label, shape, cref, sref, assy)


class Stringer(SurfacePart):
    """
    Discretely modeled stringer.
    """

    def __init__(self, label, shape, cref=None, sref=None, assy=None):
        super(Stringer, self).__init__(label, shape, cref, sref, assy)
