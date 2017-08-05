from OCC.BRepBuilderAPI import BRepBuilderAPI_Sewing
from OCC.ShapeFix import ShapeFix_Shape
from OCC.ShapeUpgrade import (ShapeUpgrade_ShapeDivideClosed,
                              ShapeUpgrade_ShapeDivideContinuity,
                              ShapeUpgrade_UnifySameDomain)
from OCC.TopoDS import topods_Edge

from afem.topology.explore import ExploreShape

__all__ = ["FixShape", "DivideClosedShape", "DivideC0Shape", "UnifyShape",
           "SewShape"]


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
        self._shape = tool.Result

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
                 cut_free_edges=False, non_manifold=False, ):
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


if __name__ == "__main__":
    import doctest

    doctest.testmod()
