#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2017 Laughlin Research, L.L.C.
#
# This file is subject to the license agreement that was delivered
# with this source code.
#
# THE SOFTWARE AND INFORMATION ARE PROVIDED ON AN "AS IS" BASIS,
# WITHOUT ANY WARRANTIES OR REPRESENTATIONS EXPRESS, IMPLIED OR 
# STATUTORY; INCLUDING, WITHOUT LIMITATION, WARRANTIES OF QUALITY,
# PERFORMANCE, MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

from OCCT.BRepBuilderAPI import BRepBuilderAPI_Sewing
from OCCT.GeomAbs import GeomAbs_C1
from OCCT.ShapeBuild import ShapeBuild_ReShape
from OCCT.ShapeUpgrade import (ShapeUpgrade_ShapeDivideClosed,
                               ShapeUpgrade_ShapeDivideContinuity,
                               ShapeUpgrade_UnifySameDomain)
from OCCT.TopTools import (TopTools_DataMapOfShapeShape,
                           TopTools_IndexedMapOfShape)
from OCCT.TopoDS import TopoDS

from afem.occ.utils import to_lst_from_toptools_listofshape
from afem.topology.create import CompoundByShapes
from afem.topology.explore import ExploreShape

__all__ = ["DivideClosedShape", "DivideContinuityShape", "DivideC0Shape",
           "UnifyShape", "SewShape", "RebuildShapeWithShapes",
           "RebuildShapeByTool", "RebuildShapesByTool"]


class DivideClosedShape(object):
    """
    Divide all closed faces in a shape.

    :param OCCT.TopoDS.TopoDS_Shape shape: The shape.
    """

    def __init__(self, shape):
        tool = ShapeUpgrade_ShapeDivideClosed(shape)
        tool.Perform()
        self._shape = tool.Result()

    @property
    def shape(self):
        """
        :return: The divided shape.
        :rtype: OCCT.TopoDS.TopoDS_Shape
        """
        return self._shape


class DivideContinuityShape(object):
    """
    Divide a shape for a given continuity and tolerance.

    :param OCCT.TopoDS.TopoDS_Shape shape: The shape.
    :param float tol: The tolerance.
    :param OCCT.GeomAbs.GeomAbs_Shape continuity:
    """

    def __init__(self, shape, tol=0.001, continuity=GeomAbs_C1):
        tool = ShapeUpgrade_ShapeDivideContinuity(shape)
        tool.SetTolerance(tol)
        tool.SetBoundaryCriterion(continuity)
        tool.SetPCurveCriterion(continuity)
        tool.SetSurfaceCriterion(continuity)
        tool.Perform()
        self._shape = tool.Result()

    @property
    def shape(self):
        """
        :return: The divided shape.
        :rtype: OCCT.TopoDS.TopoDS_Shape
        """
        return self._shape


class DivideC0Shape(DivideContinuityShape):
    """
    Divide a shape at all C0 boundaries to form a C1 shape.

    :param OCCT.TopoDS.TopoDS_Shape shape: The shape.
    :param float tol: The tolerance.
    """

    def __init__(self, shape, tol=0.001):
        super(DivideC0Shape, self).__init__(shape, tol, GeomAbs_C1)


class UnifyShape(object):
    """
    Unify edges and faces of a shape that lie on the same geometry.

    :param OCCT.TopoDS.TopoDS_Shape shape: The shape.
    :param bool edges: Option to unify all possible edges.
    :param bool faces: Option to unify all possible faces.
    :param bool bsplines: Option to concatenate the curves of edges if they
        are C1 continuous.
    """

    def __init__(self, shape, edges=True, faces=True, bsplines=False):
        tool = ShapeUpgrade_UnifySameDomain(shape, edges, faces, bsplines)
        tool.Build()
        self._shape = tool.Shape()
        self._history = tool.History()

    @property
    def shape(self):
        """
        :return: The unified shape.
        :rtype: OCCT.TopoDS.TopoDS_Shape
        """
        return self._shape

    def modified(self, shape):
        """
        Return a list of shapes modified from the given shape.

        :param OCCT.TopoDS.TopoDS_Shape shape: The shape.

        :return: List of modified shapes.
        :rtype: list[OCCT.TopoDS.TopoDS_Shape]
        """
        return to_lst_from_toptools_listofshape(self._history.Modified(shape))

    def generated(self, shape):
        """
        Return a list of shapes generated from the given shape.

        :param OCCT.TopoDS.TopoDS_Shape shape: The shape.

        :return: List of generated shapes.
        :rtype: list[OCCT.TopoDS.TopoDS_Shape]
        """
        return to_lst_from_toptools_listofshape(self._history.Generated(shape))

    def is_deleted(self, shape):
        """
        Check to see if shape is deleted.

        :param OCCT.TopoDS.TopoDS_Shape shape: The shape.

        :return: *True* if deleted, *False* if not.
        :rtype: bool
        """
        return self._history.IsRemoved(shape)


class SewShape(object):
    """
    Sew the shape.

    :param OCCT.TopoDS.TopoDS_Shape shape: The context shape to sew.
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

    .. _BRepBuilderAPI_Sewing: https://www.opencascade.com/doc/occt-7.2.0/refman/html/class_b_rep_builder_a_p_i___sewing.html

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
                tol = ExploreShape.global_tolerance(shape)

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

        :param OCCT.TopoDS.TopoDS_Shape shape: The shape.

        :return: None.
        """
        self._tool.Load(shape)

    def add(self, shape):
        """
        Add a shape to be sewed or controlled.

        :param OCCT.TopoDS.TopoDS_Shape shape: The shape.

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
        :rtype: OCCT.TopoDS.TopoDS_Shape or OCCT.TopoDS.TopoDS_Face or
            OCCT.TopoDS.TopoDS_Shell or OCCT.TopoDS.TopoDS_Solid or
            OCCT.TopoDS.TopoDS_Compound
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
        :rtype: list[OCCT.TopoDS.TopoDS_Edge]
        """
        edges = []
        for i in range(1, self.n_free_edges + 1):
            e = TopoDS.Edge_(self._tool.FreeEdge(i))
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
        :rtype: list[OCCT.TopoDS.TopoDS_Edge]
        """
        edges = []
        for i in range(1, self.n_free_edges + 1):
            e = TopoDS.Edge_(self._tool.MultipleEdge(i))
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
        :rtype: list[OCCT.TopoDS.TopoDS_Edge]
        """
        edges = []
        for i in range(1, self.n_free_edges + 1):
            e = TopoDS.Edge_(self._tool.ContigousEdge(i))
            edges.append(e)
        return edges

    def is_modified(self, shape):
        """
        Check to see if input shape has been modified.

        :param OCCT.TopoDS.TopoDS_Shape shape: The shape.

        :return: *True* if modified, *False* if not.
        :rtype: bool
        """
        self._tool.IsModified(shape)

    def modified(self, shape):
        """
        Get a modified shape.

        :param OCCT.TopoDS.TopoDS_Shape shape: The shape.

        :return: The modified shape.
        :rtype: OCCT.TopoDS.TopoDS_Shape
        """
        return self._tool.Modified(shape)

    def is_modified_subshape(self, subshape):
        """
        Check to see if input sub-shape has been modified.

        :param OCCT.TopoDS.TopoDS_Shape subshape: The shape.

        :return: *True* if modified, *False* if not.
        :rtype: bool
        """
        self._tool.IsModifiedSubShape(subshape)

    def modified_subshape(self, subshape):
        """
        Get a modified sub-shape.

        :param OCCT.TopoDS.TopoDS_Shape subshape: The shape.

        :return: The modified sub-shape.
        :rtype: OCCT.TopoDS.TopoDS_Shape
        """
        return self._tool.ModifiedSubShape(subshape)


class RebuildShapeWithShapes(object):
    """
    Rebuild a shape by requesting substitutions on a shape.

    :param OCCT.TopoDS.TopoDS_Shape old_shape: The old shape that will be
        rebuilt.
    """

    def __init__(self, old_shape):
        self._tool = ShapeBuild_ReShape()
        self._old_shape = old_shape

    def remove(self, old_shape):
        """
        Request to remove the old shape.

        :param OCCT.TopoDS.TopoDS_Shape old_shape: The old shape. This is
            usually a sub-shape of the original old shape.

        :return: None.
        """
        self._tool.Remove(old_shape)

    def replace(self, old_shape, new_shapes):
        """
        Request to replace the old shape with a list of new shapes.

        :param OCCT.TopoDS.TopoDS_Shape old_shape: The old shape. This is
            usually a sub-shape of the original old shape.
        :param list[OCCT.TopoDS.TopoDS_Shape] new_shapes: The new shapes.

        :return: None.
        """
        cmp = CompoundByShapes(new_shapes).compound
        self._tool.Replace(old_shape, cmp)

    def apply(self):
        """
        Apply the substitutions to the original old shape and return a new
        shape.

        :return: The new shape.
        :rtype: OCCT.TopoDS.TopoDS_Shape
        """
        return self._tool.Apply(self._old_shape)


class RebuildShapeByTool(object):
    """
    Rebuild a shape using a supported tool. This method will first try to make
    substitutions on the faces of the shape. If not faces exist it will try
    the edges. If no edges exist it will try the vertices.

    :param OCCT.TopoDS.TopoDS_Shape old_shape: The old shape.
    :param tool: The tool.
    :type tool: afem.topology.bop.BopCore

    :raise ValueError: If there are no sub-shapes to substitute.

    For more information see ShapeBuild_ReShape_.

    .. _ShapeBuild_ReShape: https://www.opencascade.com/doc/occt-7.2.0/refman/html/class_shape_build___re_shape.html

    Usage:

    >>> from afem.geometry import *
    >>> from afem.topology import *
    >>> pln1 = PlaneByAxes(axes='xy').plane
    >>> box1 = SolidByPlane(pln1, 10., 10., 10.).solid
    >>> pln2 = PlaneByAxes((1., 1., 1.), 'xy').plane
    >>> box2 = SolidByPlane(pln2, 5., 15., 5.).solid
    >>> cut = CutShapes(box1, box2)
    >>> assert cut.is_done
    >>> rebuild = RebuildShapeByTool(box1, cut)
    >>> new_shape = rebuild.new_shape
    >>> CheckShape.is_solid(box1)
    True
    >>> shape = FixShape(new_shape).shape
    >>> CheckShape.is_shell(shape)
    True
    """

    def __init__(self, old_shape, tool):
        reshape = ShapeBuild_ReShape()

        # Old shapes
        old_shapes = ExploreShape.get_faces(old_shape)
        if not old_shapes:
            old_shapes = ExploreShape.get_edges(old_shape)
        if not old_shapes:
            old_shapes = ExploreShape.get_vertices(old_shape)
        if not old_shapes:
            raise ValueError('No sub-shapes to substitute.')

        # Delete and replace
        for shape in old_shapes:
            # Deleted
            if tool.is_deleted(shape):
                reshape.Remove(shape)
                continue

            # Modified
            mod_shapes = [s for s in tool.modified(shape)
                          if not s.IsSame(shape)]

            if mod_shapes:
                new_shape = CompoundByShapes(mod_shapes).compound
                reshape.Replace(shape, new_shape)

        self._new_shape = reshape.Apply(old_shape)

    @property
    def new_shape(self):
        """
        :return: The new shape after substitutions.
        :rtype: OCCT.TopoDS.TopoDS_Shape
        """
        return self._new_shape


class RebuildShapesByTool(object):
    """
    Rebuild multiple shapes using a supported tool. This method is intended
    to address the case where modified shapes from an old shape are not
    duplicated in adjacent shapes, like what would happen if rebuilding a
    single shape at a time without any context. If a modified shape has
    already been replaced in an old shape and is encountered again,
    it is not substituted in the later shape. This method will first try to
    make substitutions on the faces of the shape. If not faces exist it will
    try the edges. If no edges exist it will try the vertices.

    :param list[OCCT.TopoDS.TopoDS_Shape] old_shapes: The old shapes.
    :param tool: The tool.
    :type tool: afem.topology.bop.BopCore
    """

    def __init__(self, old_shapes, tool):
        reshape = ShapeBuild_ReShape()

        self._new_shapes = TopTools_DataMapOfShapeShape()
        index_map = TopTools_IndexedMapOfShape()

        for old_shape in old_shapes:
            # Old shapes
            shapes = ExploreShape.get_faces(old_shape)
            if not shapes:
                shapes = ExploreShape.get_edges(old_shape)
            if not shapes:
                shapes = ExploreShape.get_vertices(old_shape)
            if not shapes:
                continue

            # Delete and replace
            for shape in shapes:
                # Deleted
                if tool.is_deleted(shape):
                    reshape.Remove(shape)
                    continue

                # Modified considering shapes already used
                mod_shapes = tool.modified(shape)
                replace_shapes = []
                for mod_shape in mod_shapes:
                    if index_map.Contains(mod_shape):
                        continue
                    replace_shapes.append(mod_shape)
                    index_map.Add(mod_shape)

                if replace_shapes:
                    new_shape = CompoundByShapes(replace_shapes).compound
                    reshape.Replace(shape, new_shape)

            new_shape = reshape.Apply(old_shape)
            self._new_shapes.Bind(old_shape, new_shape)

    def new_shape(self, old_shape):
        """
        Get the new shape from the old shape.

        :param OCCT.TopoDS.TopoDS_Shape old_shape: The old shape provided in
            the initial inputs.

        :return: The new shape after substitutions.
        :rtype: OCCT.TopoDS.TopoDS_Shape

        :raises RuntimeError: If the old shape is not a key in the final
            results.
        """
        return self._new_shapes.Find(old_shape)


if __name__ == "__main__":
    import doctest

    doctest.testmod()
