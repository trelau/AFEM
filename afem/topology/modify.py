from OCC.BRepBuilderAPI import BRepBuilderAPI_Sewing
from OCC.ShapeBuild import ShapeBuild_ReShape
from OCC.ShapeFix import ShapeFix_Shape
from OCC.ShapeUpgrade import (ShapeUpgrade_ShapeDivideClosed,
                              ShapeUpgrade_ShapeDivideContinuity,
                              ShapeUpgrade_UnifySameDomain)
from OCC.TopoDS import topods_Edge

from afem.topology.create import CompoundByShapes
from afem.topology.explore import ExploreShape

__all__ = ["FixShape", "DivideClosedShape", "DivideC0Shape", "UnifyShape",
           "SewShape", "RebuildShapeWithShapes", "RebuildShapeByTool"]


class FixShape(object):
    """
    Attempt to fix the shape by applying a number of general fixes.

    :param OCC.TopoDS.TopoDS_Shape shape: The shape.
    :param float min_tol: Minimum allowed tolerance.
    :param float max_tol: Maximum allowed tolerance.
    """

    def __init__(self, shape, min_tol=None, max_tol=None):
        fix = ShapeFix_Shape(shape)
        if min_tol is not None:
            fix.SetMinTolerance(min_tol)
        if max_tol is not None:
            fix.SetMaxTolerance(max_tol)
        fix.Perform()
        self._shape = fix.Shape()

    @property
    def shape(self):
        """
        :return: The fixed shape.
        :rtype: OCC.TopoDS.TopoDS_Shape
        """
        return self._shape


class DivideClosedShape(object):
    """
    Divide all closed faces in a shape.

    :param OCC.TopoDS.TopoDS_Shape shape: The shape.
    """

    def __init__(self, shape):
        tool = ShapeUpgrade_ShapeDivideClosed(shape)
        tool.Perform()
        self._shape = tool.Result()

    @property
    def shape(self):
        """
        :return: The divided shape.
        :rtype: OCC.TopoDS.TopoDS_Shape
        """
        return self._shape


class DivideC0Shape(object):
    """
    Divide a shape at all C0 boundaries to form a C1 shape.

    :param OCC.TopoDS.TopoDS_Shape shape: The shape.
    """

    def __init__(self, shape):
        tool = ShapeUpgrade_ShapeDivideContinuity(shape)
        tool.Perform()
        self._shape = tool.Result()

    @property
    def shape(self):
        """
        :return: The divided shape.
        :rtype: OCC.TopoDS.TopoDS_Shape
        """
        return self._shape


class UnifyShape(object):
    """
    Unify edges and faces of a shape that lie on the same geometry.

    :param OCC.TopoDS.TopoDS_Shape shape: The shape.
    :param bool edges: Option to unify all possible edges.
    :param bool faces: Option to unify all possible faces.
    :param bool bsplines: Option to concatenate the curves of edges if they
        are C1 continuous.
    """

    def __init__(self, shape, edges=True, faces=True, bsplines=False):
        tool = ShapeUpgrade_UnifySameDomain(shape, edges, faces, bsplines)
        tool.Build()
        self._shape = tool.Shape()

    @property
    def shape(self):
        """
        :return: The unified shape.
        :rtype: OCC.TopoDS.TopoDS_Shape
        """
        return self._shape

    def generated(self, shape):
        """
        Get a shape generated from the old one.

        :param OCC.TopoDS.TopoDS_Shape shape: The old shape.

        :return: The new shape.
        :rtype: OCC.TopoDS.TopoDS_Shape
        """


class SewShape(object):
    """
    Sew the shape.

    :param OCC.TopoDS.TopoDS_Shape shape: The context shape to sew.
    :param float tol: Sewing tolerance. If *None* is provided then the
        average tolerance of the shape will be used. If no shape is
        provided, then a default value of 1.0e-7 is used.
    :param float min_tol: Minimum tolerance.
    :param float max_tol: Maximum tolerance.
    :param bool cut_free_edges: Option for cutting of free edges.
    :param bool non_manifold: Option for non-manifold processing.

    .. note::

        If *shape* is *None* then the user is expected to manually load the
        shape and perform the operation.

    For more information see BRepBuilderAPI_Sewing_.

    .. _BRepBuilderAPI_Sewing: https://www.opencascade.com/doc/occt-7.1.0/refman/html/class_b_rep_builder_a_p_i___sewing.html

    Usage:

    >>> from afem.topology import *
    >>> p1 = (0., 0., 0.)
    >>> p2 = (1., 0., 0.)
    >>> p3 = (1., 1., 0.)
    >>> p4 = (0., 1., 0.)
    >>> w1 = WireByPoints([p1, p2, p3, p4], True).wire
    >>> f1 = FaceByPlanarWire(w1).face
    >>> p1 = (0., 0., 0.)
    >>> p2 = (0., 0., 1.)
    >>> p3 = (0., 1., 1.)
    >>> p4 = (0., 1., 0.)
    >>> w2 = WireByPoints([p1, p2, p3, p4], True).wire
    >>> f2 = FaceByPlanarWire(w2).face
    >>> p1 = (0., 0.5, 0.)
    >>> p2 = (-1., 0.5, 0.)
    >>> p3 = (-1., 1.5, 0.)
    >>> p4 = (0., 1.5, 0.)
    >>> w3 = WireByPoints([p1, p2, p3, p4], True).wire
    >>> f3 = FaceByPlanarWire(w3).face
    >>> # Build an unconnected shell
    >>> shell = ShellByFaces([f1, f2, f3]).shell
    >>> tool = SewShape(shell, non_manifold=True)
    >>> shape = tool.sewed_shape
    >>> len(ExploreShape.get_faces(shape))
    3
    """

    def __init__(self, shape=None, tol=None, min_tol=None, max_tol=None,
                 cut_free_edges=False, non_manifold=False):
        if tol is None:
            if shape is None:
                tol = 1.0e-7
            else:
                tol = ExploreShape.get_tolerance(shape)

        self._tool = BRepBuilderAPI_Sewing(tol, True, True, cut_free_edges,
                                           non_manifold)

        if min_tol is not None:
            self._tool.SetMinTolerance(min_tol)
        if max_tol is not None:
            self._tool.SetMaxTolerance(max_tol)

        if shape is not None:
            self._tool.Load(shape)
            self._tool.Perform()

    def load(self, shape):
        """
        Load the context shape to sew.

        :param OCC.TopoDS.TopoDS_Shape shape: The shape.

        :return: None.
        """
        self._tool.Load(shape)

    def add(self, shape):
        """
        Add a shape to sew.

        :param OCC.TopoDS.TopoDS_Shape shape: The shape.

        :return: None.
        """
        self._tool.Add(shape)

    def perform(self):
        """
        Perform the sewing operation.

        :return: None.
        """
        self._tool.Perform()

    @property
    def sewed_shape(self):
        """
        :return: The sewed shape. May be a null shape if nothing is
            constructed.
        :rtype: OCC.TopoDS.TopoDS_Shape or OCC.TopoDS.TopoDS_Face or
            OCC.TopoDS.TopoDS_Shell or OCC.TopoDS.TopoDS_Solid or
            OCC.TopoDS.TopoDS_Compound
        """
        return self._tool.SewedShape()

    @property
    def n_free_edges(self):
        """
        :return: Number of free edges.
        :rtype: int
        """
        return self._tool.NbFreeEdges()

    @property
    def free_edges(self):
        """
        :return: Free edges.
        :rtype: list[OCC.TopoDS.TopoDS_Edge]
        """
        edges = []
        for i in range(1, self.n_free_edges + 1):
            e = topods_Edge(self._tool.FreeEdge(i))
            edges.append(e)
        return edges

    @property
    def n_multiple_edges(self):
        """
        :return: Number of edges connected to more than two faces.
        :rtype: int
        """
        return self._tool.NbMultipleEdges()

    @property
    def multiple_edges(self):
        """
        :return: Multiple edges.
        :rtype: list[OCC.TopoDS.TopoDS_Edge]
        """
        edges = []
        for i in range(1, self.n_free_edges + 1):
            e = topods_Edge(self._tool.MultipleEdge(i))
            edges.append(e)
        return edges

    @property
    def n_manifold_edges(self):
        """
        :return: Number of manifold edges.
        :rtype: int
        """
        return self._tool.NbContigousEdges()

    @property
    def manifold_edges(self):
        """
        :return: Manifold edges.
        :rtype: list[OCC.TopoDS.TopoDS_Edge]
        """
        edges = []
        for i in range(1, self.n_free_edges + 1):
            e = topods_Edge(self._tool.ContigousEdge(i))
            edges.append(e)
        return edges

    def is_modified(self, shape):
        """
        Check to see if input shape has been modified.

        :param OCC.TopoDS.TopoDS_Shape shape: The shape.

        :return: *True* if modified, *False* if not.
        :rtype: bool
        """
        self._tool.IsModified(shape)

    def modified(self, shape):
        """
        Get a modified shape.

        :param OCC.TopoDS.TopoDS_Shape shape: The shape.

        :return: The modified shape.
        :rtype: OCC.TopoDS.TopoDS_Shape
        """
        return self._tool.Modified(shape)

    def is_modified_subshape(self, subshape):
        """
        Check to see if input sub-shape has been modified.

        :param OCC.TopoDS.TopoDS_Shape subshape: The shape.

        :return: *True* if modified, *False* if not.
        :rtype: bool
        """
        self._tool.IsModifiedSubShape(subshape)

    def modified_subshape(self, subshape):
        """
        Get a modified sub-shape.

        :param OCC.TopoDS.TopoDS_Shape subshape: The shape.

        :return: The modified sub-shape.
        :rtype: OCC.TopoDS.TopoDS_Shape
        """
        return self._tool.ModifiedSubShape(subshape)


class RebuildShapeWithShapes(object):
    """
    Rebuild a shape by requesting substitutions on a shape.

    :param OCC.TopoDS.TopoDS_Shape old_shape: The old shape that will be
        rebuilt.
    """

    def __init__(self, old_shape):
        self._tool = ShapeBuild_ReShape()
        self._old_shape = old_shape

    def remove(self, old_shape):
        """
        Request to remove the old shape.

        :param OCC.TopoDS.TopoDS_Shape old_shape: The old shape. This is
            usually a sub-shape of the original old shape.

        :return: None.
        """
        self._tool.Remove(old_shape)

    def replace(self, old_shape, new_shapes):
        """
        Request to replace the old shape with a list of new shapes.

        :param OCC.TopoDS.TopoDS_Shape old_shape: The old shape. This is
            usually a sub-shape of the original old shape.
        :param list[OCC.TopoDS.TopoDS_Shape] new_shapes: The new shapes.

        :return: None.
        """
        cmp = CompoundByShapes(new_shapes).compound
        self._tool.Replace(old_shape, cmp)

    def apply(self, fix=False):
        """
        Apply the substitutions to the original old shape and return a new
        shape.

        :param bool fix: Option to use :class:`.FixShape` on the new shape in
            case the substitutions caused errors (e.g., like a solid is now a
            shell).

        :return: The new shape.
        :rtype: OCC.TopoDS.TopoDS_Shape
        """
        new_shape = self._tool.Apply(self._old_shape)
        if not fix:
            return new_shape
        return FixShape(new_shape).shape


class RebuildShapeByTool(object):
    """
    Rebuild a shape using a supported tool.

    :param OCC.TopoDS.TopoDS_Shape: The shape.
    :param tool: The tool.
    :type tool: afem.topology.bop.BopAlgo
    :param str replace_type: The level of shape to replace ('vertex',
        'edge', 'face').
    :param bool fix: Option to use :class:`.FixShape` on the new shape in
        case the substitutions caused errors (e.g., like a solid is now a
        shell).

    For more information see ShapeBuild_ReShape_.

    .. _ShapeBuild_ReShape: https://www.opencascade.com/doc/occt-7.1.0/refman/html/class_shape_build___re_shape.html

    Usage:

    >>> from afem.geometry import *
    >>> from afem.topology import *
    >>> pln1 = PlaneByAxes(axes='xy').plane
    >>> box1 = SolidByPlane(pln1, 10., 10., 10.).solid
    >>> pln2 = PlaneByAxes((1., 1., 1.), 'xy').plane
    >>> box2 = SolidByPlane(pln2, 5., 15., 5.).solid
    >>> cut = CutShapes(box1, box2)
    >>> assert cut.is_done
    >>> rebuild = RebuildShapeByTool(box1, cut, fix=True)
    >>> new_shape = rebuild.new_shape
    >>> CheckShape.is_solid(box1)
    True
    >>> CheckShape.is_shell(new_shape)
    True
    """

    def __init__(self, old_shape, tool, replace_type='face', fix=False):
        reshape = ShapeBuild_ReShape()

        # TODO Consider iterating through all shapes?

        # Old shapes
        replace_type = replace_type.lower()
        if replace_type in ['v', 'vertex', 'vertices']:
            old_shapes = ExploreShape.get_vertices(old_shape)
        elif replace_type in ['e', 'edge', 'edges']:
            old_shapes = ExploreShape.get_edges(old_shape)
        elif replace_type in ['f', 'face', 'faces']:
            old_shapes = ExploreShape.get_faces(old_shape)
        else:
            msg = 'Invalid replace type.'
            raise TypeError(msg)

        # Delete and replace
        for shape in old_shapes:
            # Deleted
            if tool.is_deleted(shape):
                reshape.Remove(shape)
                continue

            # Modified
            mod_shapes = tool.modified(shape)
            if mod_shapes:
                new_shape = CompoundByShapes(mod_shapes).compound
                reshape.Replace(shape, new_shape)

        new_shape = reshape.Apply(old_shape)
        if fix:
            self._new_shape = FixShape(new_shape).shape
        else:
            self._new_shape = new_shape

    @property
    def new_shape(self):
        """
        :return: The new shape after substitutions.
        :rtype: OCC.TopoDS.TopoDS_Shape
        """
        return self._new_shape


if __name__ == "__main__":
    import doctest

    doctest.testmod()
