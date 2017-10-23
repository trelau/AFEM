from warnings import warn

from OCC.BOPAlgo import BOPAlgo_MakerVolume
from OCC.BRepAlgoAPI import (BRepAlgoAPI_Common, BRepAlgoAPI_Cut,
                             BRepAlgoAPI_Fuse, BRepAlgoAPI_Section)
from OCC.BRepFeat import BRepFeat_MakeCylindricalHole, BRepFeat_SplitShape
from OCC.GEOMAlgo import GEOMAlgo_Splitter
from OCC.TopoDS import TopoDS_Face, topods_Solid

from afem.geometry.check import CheckGeom
from afem.occ.utils import (to_lst_from_toptools_listofshape,
                            to_toptools_listofshape)
from afem.topology.check import CheckShape
from afem.topology.explore import ExploreShape

__all__ = ["BopAlgo", "FuseShapes", "CutShapes", "CommonShapes",
           "IntersectShapes", "SplitShapes", "VolumesFromShapes",
           "CutCylindricalHole", "LocalSplit"]


class BopCore(object):
    """
    Core class for Boolean operations and enabling attributes and methods for
    rebuilding shapes.
    """

    def __init__(self):
        self._bop = None

    def build(self):
        """
        Build the results.

        :return: None.
        """
        if isinstance(self._bop, (GEOMAlgo_Splitter, BOPAlgo_MakerVolume)):
            self._bop.Perform()
        else:
            self._bop.Build()

    @property
    def is_done(self):
        """
        :return: *True* if operation is done, *False* if not.
        :rtype: bool
        """
        if isinstance(self._bop, (GEOMAlgo_Splitter, BOPAlgo_MakerVolume,
                                  BRepFeat_MakeCylindricalHole)):
            return self._bop.ErrorStatus() == 0
        return self._bop.IsDone()

    @property
    def shape(self):
        """
        :return: The resulting shape.
        :rtype: OCC.TopoDS.TopoDS_Shape
        """
        return self._bop.Shape()

    def modified(self, shape):
        """
        Return a list of shapes modified from the given shape.

        :param OCC.TopoDS.TopoDS_Shape shape: The shape.

        :return: List of modified shapes.
        :rtype: list[OCC.TopoDS.TopoDS_Shape]
        """
        return to_lst_from_toptools_listofshape(self._bop.Modified(shape))

    def generated(self, shape):
        """
        Return a list of shapes generated from the given shape.

        :param OCC.TopoDS.TopoDS_Shape shape: The shape.

        :return: List of generated shapes.
        :rtype: list[OCC.TopoDS.TopoDS_Shape]
        """
        return to_lst_from_toptools_listofshape(self._bop.Generated(shape))

    def is_deleted(self, shape):
        """
        Check to see if shape is deleted.

        :param OCC.TopoDS.TopoDS_Shape shape: The shape.

        :return: *True* if deleted, *False* if not.
        :rtype: bool
        """
        return self._bop.IsDeleted(shape)


class BopAlgo(BopCore):
    """
    Base class for Boolean operations.

    :param shape1: The first shape.
    :type shape1: OCC.TopoDS.TopoDS_Shape or None
    :param shape2: The second shape.
    :type shape2: OCC.TopoDS.TopoDS_Shape or None
    :param bool parallel: Option to run in parallel mode.
    :param float fuzzy_val: Fuzzy tolerance value.
    :param bop: The OpenCASCADE class for the Boolean operation.

    .. note::

        If *shape1* or *shape2* is *None* then the user is expected to manually
        set the arguments and tools and build the result.
    """

    def __init__(self, shape1, shape2, parallel, fuzzy_val, bop):
        super(BopAlgo, self).__init__()

        if CheckShape.is_shape(shape1) and CheckShape.is_shape(shape2):
            self._bop = bop(shape1, shape2)
        else:
            self._bop = bop()

        if parallel:
            self._bop.SetRunParallel(True)
        else:
            self._bop.SetRunParallel(False)

        if fuzzy_val is not None:
            self._bop.SetFuzzyValue(fuzzy_val)

    def set_args(self, shapes):
        """
        Set the arguments.

        :param list[OCC.TopoDS.TopoDS_Shape] shapes: The arguments.

        :return: None.
        """
        if isinstance(self._bop, (GEOMAlgo_Splitter, BOPAlgo_MakerVolume)):
            for shape in shapes:
                self._bop.AddArgument(shape)
            return None
        args = to_toptools_listofshape(shapes)
        self._bop.SetArguments(args)

    def set_tools(self, shapes):
        """
        Set the tools.

        :param list[OCC.TopoDS.TopoDS_Shape] shapes: The tools.

        :return: None.
        """
        if isinstance(self._bop, BOPAlgo_MakerVolume):
            warn('Setting tools not available. Doing nothing.', RuntimeWarning)
            return None

        if isinstance(self._bop, GEOMAlgo_Splitter):
            for shape in shapes:
                self._bop.AddTool(shape)
            return None

        tools = to_toptools_listofshape(shapes)
        self._bop.SetTools(tools)

    @property
    def vertices(self):
        """
        :return: The vertices of the resulting shape.
        :rtype: list[OCC.TopoDS.TopoDS_Vertex]
        """
        return ExploreShape.get_vertices(self.shape)

    @property
    def edges(self):
        """
        :return: The edges of the resulting shape.
        :rtype: list[OCC.TopoDS.TopoDS_Edge]
        """
        return ExploreShape.get_edges(self.shape)

    def refine_edges(self):
        """
        Fuse C1 edges.

        :return: None.
        """
        if isinstance(self._bop, (GEOMAlgo_Splitter, BOPAlgo_MakerVolume,
                                  BRepFeat_MakeCylindricalHole)):
            warn('Refining edges not available. Doing nothing.',
                 RuntimeWarning)
        else:
            self._bop.RefineEdges()

    @property
    def fuse_edges(self):
        """
        :return: The result flag of edge refining.
        :rtype: bool
        """
        if isinstance(self._bop, (GEOMAlgo_Splitter, BOPAlgo_MakerVolume,
                                  BRepFeat_MakeCylindricalHole)):
            return False
        else:
            return self._bop.FuseEdges()

    @property
    def section_edges(self):
        """
        :return: A list of section edges as a result of intersection between
            the shapes.
        :rtype: list[OCC.TopoDS.TopoDS_Edge]
        """
        if isinstance(self._bop, (GEOMAlgo_Splitter, BOPAlgo_MakerVolume,
                                  BRepFeat_MakeCylindricalHole)):
            warn('Getting section edges not available. Returning an empty '
                 'list.', RuntimeWarning)
            return []
        else:
            return to_lst_from_toptools_listofshape(self._bop.SectionEdges())

    @property
    def has_modified(self):
        """
        :return: *True* if there is at least one modified shape.
        :rtype: bool
        """
        return self._bop.HasModified()

    @property
    def has_generated(self):
        """
        :return: *True* if there is at least one generated shape.
        :rtype: bool
        """
        return self._bop.HasGenerated()

    @property
    def has_deleted(self):
        """
        :return: *True* if there is at least one deleted shape.
        :rtype: bool
        """
        return self._bop.HasDeleted()


class FuseShapes(BopAlgo):
    """
    Boolean fuse operation.

    :param shape1: The first shape.
    :type shape1: OCC.TopoDS.TopoDS_Shape or None
    :param shape2: The second shape.
    :type shape2: OCC.TopoDS.TopoDS_Shape or None
    :param bool parallel: Option to run in parallel mode.
    :param float fuzzy_val: Fuzzy tolerance value.

    .. note::

        If *shape1* or *shape2* is *None* then the user is expected to manually
        set the arguments and tools and build the result.

    For more information see BRepAlgoAPI_Fuse_.

    .. _BRepAlgoAPI_Fuse: https://www.opencascade.com/doc/occt-7.1.0/refman/html/class_b_rep_algo_a_p_i___fuse.html

    Usage:

    >>> from afem.topology import *
    >>> e1 = EdgeByPoints((0., 0., 0.), (10., 0., 0.)).edge
    >>> e2 = EdgeByPoints((5., 1., 0.), (5., -1., 0.)).edge
    >>> bop = FuseShapes(e1, e2)
    >>> assert bop.is_done
    >>> shape = bop.shape
    >>> # Setting arguments and tools
    >>> bop = FuseShapes()
    >>> bop.set_args([e1])
    >>> bop.set_tools([e2])
    >>> bop.build()
    >>> assert bop.is_done
    """

    def __init__(self, shape1=None, shape2=None, parallel=True,
                 fuzzy_val=None):
        super(FuseShapes, self).__init__(shape1, shape2, parallel,
                                         fuzzy_val, BRepAlgoAPI_Fuse)


class CutShapes(BopAlgo):
    """
    Boolean cut operation.

    :param shape1: The first shape.
    :type shape1: OCC.TopoDS.TopoDS_Shape or None
    :param shape2: The second shape.
    :type shape2: OCC.TopoDS.TopoDS_Shape or None
    :param bool parallel: Option to run in parallel mode.
    :param float fuzzy_val: Fuzzy tolerance value.

    .. note::

        If *shape1* or *shape2* is *None* then the user is expected to manually
        set the arguments and tools and build the result.

    For more information see BRepAlgoAPI_Cut_.

    .. _BRepAlgoAPI_Cut: https://www.opencascade.com/doc/occt-7.1.0/refman/html/class_b_rep_algo_a_p_i___cut.html

    Usage:

    >>> from afem.topology import *
    >>> e1 = EdgeByPoints((0., 0., 0.), (10., 0., 0.)).edge
    >>> e2 = EdgeByPoints((5., 1., 0.), (5., -1., 0.)).edge
    >>> bop = CutShapes(e1, e2)
    >>> assert bop.is_done
    >>> shape = bop.shape
    >>> # Setting arguments and tools
    >>> bop = CutShapes()
    >>> bop.set_args([e1])
    >>> bop.set_tools([e2])
    >>> bop.build()
    >>> assert bop.is_done
    """

    def __init__(self, shape1=None, shape2=None, parallel=True,
                 fuzzy_val=None):
        super(CutShapes, self).__init__(shape1, shape2, parallel,
                                        fuzzy_val, BRepAlgoAPI_Cut)


class CommonShapes(BopAlgo):
    """
    Boolean common operation.

    :param shape1: The first shape.
    :type shape1: OCC.TopoDS.TopoDS_Shape or None
    :param shape2: The second shape.
    :type shape2: OCC.TopoDS.TopoDS_Shape or None
    :param bool parallel: Option to run in parallel mode.
    :param float fuzzy_val: Fuzzy tolerance value.

    .. note::

        If *shape1* or *shape2* is *None* then the user is expected to manually
        set the arguments and tools and build the result.

    For more information see BRepAlgoAPI_Common_.

    .. _BRepAlgoAPI_Common: https://www.opencascade.com/doc/occt-7.1.0/refman/html/class_b_rep_algo_a_p_i___common.html

    Usage:

    >>> from afem.topology import *
    >>> e1 = EdgeByPoints((0., 0., 0.), (10., 0., 0.)).edge
    >>> e2 = EdgeByPoints((5., 1., 0.), (5., -1., 0.)).edge
    >>> bop = CommonShapes(e1, e2)
    >>> assert bop.is_done
    >>> shape = bop.shape
    >>> # Setting arguments and tools
    >>> bop = CommonShapes()
    >>> bop.set_args([e1])
    >>> bop.set_tools([e2])
    >>> bop.build()
    >>> assert bop.is_done
    """

    def __init__(self, shape1=None, shape2=None, parallel=True,
                 fuzzy_val=None):
        super(CommonShapes, self).__init__(shape1, shape2, parallel,
                                           fuzzy_val, BRepAlgoAPI_Common)


class IntersectShapes(BopAlgo):
    """
    Boolean intersect operation.

    :param shape1: The first shape.
    :type shape1: OCC.TopoDS.TopoDS_Shape or afem.geometry.entities.Surface
    :param shape2: The second shape.
    :type shape2: OCC.TopoDS.TopoDS_Shape or afem.geometry.entities.Surface
    :param bool compute_pcurve1: Option to compute p-curves on shape 1.
    :param bool compute_pcurve2: Option to compute p-curves on shape 2.
    :param bool approximate: Option to approximate intersection curves.
    :param bool parallel: Option to run in parallel mode.
    :param float fuzzy_val: Fuzzy tolerance value.

    .. note::

        If *shape1* or *shape2* is *None* then the user is expected to manually
        set the arguments and tools and build the result.

    For more information see BRepAlgoAPI_Section_.

    .. _BRepAlgoAPI_Section: https://www.opencascade.com/doc/occt-7.1.0/refman/html/class_b_rep_algo_a_p_i___section.html

    Usage:

    >>> from afem.topology import *
    >>> e1 = EdgeByPoints((0., 0., 0.), (10., 0., 0.)).edge
    >>> e2 = EdgeByPoints((5., 1., 0.), (5., -1., 0.)).edge
    >>> bop = IntersectShapes(e1, e2)
    >>> assert bop.is_done
    >>> shape = bop.shape
    >>> # Setting arguments and tools
    >>> bop = IntersectShapes()
    >>> bop.set_args([e1])
    >>> bop.set_tools([e2])
    >>> bop.build()
    >>> assert bop.is_done
    """

    def __init__(self, shape1=None, shape2=None, compute_pcurve1=False,
                 compute_pcurve2=False, approximate=False, parallel=True,
                 fuzzy_val=None):
        super(IntersectShapes, self).__init__(None, None, parallel,
                                              fuzzy_val, BRepAlgoAPI_Section)

        self._bop.ComputePCurveOn1(compute_pcurve1)
        self._bop.ComputePCurveOn2(compute_pcurve2)
        self._bop.Approximation(approximate)

        build1, build2 = False, False
        if CheckShape.is_shape(shape1):
            self._bop.Init1(shape1)
            build1 = True
        elif CheckGeom.is_surface(shape1):
            self._bop.Init1(shape1.handle)
            build1 = True

        if CheckShape.is_shape(shape2):
            self._bop.Init2(shape2)
            build2 = True
        elif CheckGeom.is_surface(shape2):
            self._bop.Init2(shape2.handle)
            build2 = True

        if build1 and build2:
            self._bop.Build()

    def has_ancestor_face1(self, edge):
        """
        Get the ancestor face on the intersection edge on the first shape
        if available.

        :param OCC.TopoDS.TopoDS_Edge edge: The edge.

        :return: *True* and the face if available, *False* and *None* if not.
        :rtype: tuple(bool, OCC.TopoDS.TopoDS_Face or None)
        """
        f = TopoDS_Face()
        if self._bop.HasAncestorFaceOn1(edge, f):
            return True, f
        return False, None

    def has_ancestor_face2(self, edge):
        """
        Get the ancestor face on the intersection edge on the second shape
        if available.

        :param OCC.TopoDS.TopoDS_Edge edge: The edge.

        :return: *True* and the face if available, *False* and *None* if not.
        :rtype: tuple(bool, OCC.TopoDS.TopoDS_Face or None)
        """
        f = TopoDS_Face()
        if self._bop.HasAncestorFaceOn2(edge, f):
            return True, f
        return False, None


class SplitShapes(BopAlgo):
    """
    Split arbitrary shapes. This is a wrapper for the SALOME
    GEOMAlgo_Splitter tool.

    :param shape1: The first shape.
    :type shape1: OCC.TopoDS.TopoDS_Shape or None
    :param shape2: The second shape.
    :type shape2: OCC.TopoDS.TopoDS_Shape or None
    :param bool parallel: Option to run in parallel mode.
    :param float fuzzy_val: Fuzzy tolerance value.

    .. note::

        If *shape1* or *shape2* is *None* then the user is expected to manually
        set the arguments and tools and build the result.

    Usage:

    >>> from afem.topology import *
    >>> e1 = EdgeByPoints((0., 0., 0.), (10., 0., 0.)).edge
    >>> e2 = EdgeByPoints((5., 1., 0.), (5., -1., 0.)).edge
    >>> bop = SplitShapes(e1, e2)
    >>> assert bop.is_done
    >>> shape = bop.shape
    >>> # Setting arguments and tools
    >>> bop = SplitShapes()
    >>> bop.set_args([e1])
    >>> bop.set_tools([e2])
    >>> bop.build()
    >>> assert bop.is_done
    """

    def __init__(self, shape1=None, shape2=None, parallel=True,
                 fuzzy_val=None):
        super(SplitShapes, self).__init__(None, None, parallel,
                                          fuzzy_val, GEOMAlgo_Splitter)
        if CheckShape.is_shape(shape1) and CheckShape.is_shape(shape2):
            self._bop.AddArgument(shape1)
            self._bop.AddArgument(shape2)
            self._bop.Perform()

    def add_arg(self, shape):
        """
        Add an argument.

        :param OCC.TopoDS.TopoDS_Shape shape: The argument.

        :return: None.
        """
        self._bop.AddArgument(shape)

    def add_tool(self, shape):
        """
        Add a tool.

        :param OCC.TopoDS.TopoDS_Shape shape: The tool.

        :return: None.
        """
        self._bop.AddTool(shape)


class VolumesFromShapes(BopAlgo):
    """
    Build solids from a list of shapes.

    :param list[OCC.TopoDS.TopoDS_Shape] shapes: The shapes.
    :param bool intersect: Option to intersect the shapes before building
        solids.
    :param bool parallel: Option to run in parallel mode.
    :param float fuzzy_val: Fuzzy tolerance value.
    """

    def __init__(self, shapes, intersect=False, parallel=True,
                 fuzzy_val=None):
        super(VolumesFromShapes, self).__init__(None, None, parallel,
                                                fuzzy_val, BOPAlgo_MakerVolume)

        self.set_args(shapes)

        if intersect:
            self._bop.SetIntersect(True)
        else:
            self._bop.SetIntersect(False)

        self.build()

        self._solids = []
        for solid in ExploreShape.get_solids(self.shape):
            self._solids.append(topods_Solid(solid))

    @property
    def box(self):
        """
        :return: The bounding box of all provided shapes.
        :rtype: OCC.TopoDS.TopoDS_Solid
        """
        return self._bop.Box()

    @property
    def nsolids(self):
        """
        :return: The number of solids in the shape.
        :rtype: int
        """
        return len(self.solids)

    @property
    def solids(self):
        """
        :return: The list of solids.
        :rtype: list[OCC.TopoDS.TopoDS_Solid]
        """
        return self._solids


class CutCylindricalHole(BopAlgo):
    """
    Cut a cylindrical hole on a shape.

    :param shape: The shape.
    :type shape: OCC.TopoDS.TopoDS_Shape
    :param float radius: The radius of the hole.
    :param afem.geometry.entities.Axis1: The axis for the hole.
    :param bool parallel: Option to run in parallel mode.
    :param float fuzzy_val: Fuzzy tolerance value.

    Usage:

    >>> from afem.geometry import *
    >>> from afem.topology import *
    >>> pln = PlaneByAxes().plane
    >>> face = FaceByPlane(pln, -2., 2., -2., 2.).face
    >>> # The small offset is a workaround for planar cuts
    >>> p = Point(0., 0.1, 0.)
    >>> d = Direction(0., 1., 0.)
    >>> ax1 = Axis1(p, d)
    >>> bop = CutCylindricalHole(face, 1., ax1)
    >>> assert bop.is_done
    """

    def __init__(self, shape, radius, ax1, parallel=True, fuzzy_val=None):
        super(CutCylindricalHole, self).__init__(None, None, parallel,
                                                 fuzzy_val,
                                                 BRepFeat_MakeCylindricalHole)

        self._bop.Init(shape, ax1)
        self._bop.Perform(radius)


class LocalSplit(BopCore):
    """
    Perform a local split of a shape in the context of a basis shape. This tool
    only splits faces.

    :param OCC.TopoDS.TopoDS_Shape shape: The local shape.
    :param tool: The tool to split with.
    :type tool: OCC.TopoDS.TopoDS_Shape or afem.geometry.entities.Surface
    :param OCC.TopoDS.TopoDS_Shape basis_shape: The basis shape that the local
        shape is part of.
    :param bool approximate: Option to approximate intersection curves.
    :param bool parallel: Option to run in parallel mode.
    :param float fuzzy_val: Fuzzy tolerance value.
    """

    def __init__(self, shape, tool, basis_shape, approximate=False,
                 parallel=True, fuzzy_val=None):
        super(LocalSplit, self).__init__()

        # Intersect
        section = IntersectShapes(shape, tool, True, False, approximate,
                                  parallel, fuzzy_val)
        sec_edges = section.edges

        # Split
        self._bop = BRepFeat_SplitShape(basis_shape)
        for e in sec_edges:
            status, f = section.has_ancestor_face1(e)
            if status:
                self._bop.Add(e, f)
        self.build()


if __name__ == "__main__":
    import doctest

    doctest.testmod()
