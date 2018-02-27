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

from warnings import warn

from OCCT.BRepExtrema import (BRepExtrema_DistShapeShape, BRepExtrema_IsVertex,
                              BRepExtrema_IsOnEdge, BRepExtrema_IsInFace)

from afem.geometry.check import CheckGeom
from afem.geometry.entities import Point
from afem.topology.check import CheckShape

__all__ = ["DistanceShapeToShape", "DistanceShapeToShapes",
           "DistancePointToShapes"]


class DistanceShapeToShape(object):
    """
    Calculate minimum distance between two shapes.

    :param OCCT.TopoDS.TopoDS_Shape shape1: The first shape.
    :param OCCT.TopoDS.TopoDS_Shape shape2: The second shape.

    :raise RuntimeError: If OCC method fails.

    Usage:

    >>> from afem.topology import *
    >>> v1 = VertexByPoint((0., 0., 0.)).vertex
    >>> v2 = VertexByPoint((10., 0., 0.)).vertex
    >>> tool = DistanceShapeToShape(v1, v2)
    >>> tool.nsol
    1
    >>> tool.dmin
    10.0
    """

    def __init__(self, shape1, shape2, deflection=1.0e-7):
        self._tool = BRepExtrema_DistShapeShape(shape1, shape2, deflection)

    @property
    def is_done(self):
        """
        :return: *True* if algorithm is done, *False* if not.
        :rtype: bool
        """
        return self._tool.IsDone()

    @property
    def nsol(self):
        """
        :return: The number of solutions satisfying the minimum distance.
        :rtype:
        """
        return self._tool.NbSolution()

    @property
    def dmin(self):
        """
        :return: The minimum distance.
        :rtype: float
        """
        return self._tool.Value()

    @property
    def inner_solution(self):
        """
        :return: *True* if one of the shapes is a solid and the other is
            completely or partially inside the solid.
        :rtype: bool
        """
        return self._tool.InnerSolution()

    def point_on_shape1(self, n=1):
        """
        The point for the *n-th* solution on the first shape.

        :param int n: The index.

        :return: The point.
        :rtype: afem.geometry.entities.Point
        """
        gp_pnt = self._tool.PointOnShape1(n)
        return Point(gp_pnt.X(), gp_pnt.Z(), gp_pnt.Z())

    def point_on_shape2(self, n=1):
        """
        The point for the *n-th* solution on the second shape.

        :param int n: The index.

        :return: The point.
        :rtype: afem.geometry.entities.Point
        """
        gp_pnt = self._tool.PointOnShape2(n)
        return Point(gp_pnt.X(), gp_pnt.Z(), gp_pnt.Z())

    def support_type_shape1(self, n=1):
        """
        The type of support for the *n-th* solution on the first shape.

        :param int n: The index.

        :return: The support type.
        :rtype: OCCT.BRepExtrema.BRepExtrema_SupportType
        """
        return self._tool.SupportOnShape1(n)

    def is_vertex_shape1(self, n=1):
        """
        Check if support type is a vertex for the first shape.

        :param int n: The index.

        :return: *True* if a vertex, *False* otherwise.
        :rtype: bool
        """
        return self.support_type_shape1(n) == BRepExtrema_IsVertex

    def is_on_edge_shape1(self, n=1):
        """
        Check if support type is on an edge for the first shape.

        :param int n: The index.

        :return: *True* if on an edge, *False* otherwise.
        :rtype: bool
        """
        return self.support_type_shape1(n) == BRepExtrema_IsOnEdge

    def is_in_face_shape1(self, n=1):
        """
        Check if support type is in a face for the first shape.

        :param int n: The index.

        :return: *True* if in a face, *False* otherwise.
        :rtype: bool
        """
        return self.support_type_shape1(n) == BRepExtrema_IsInFace

    def support_type_shape2(self, n=1):
        """
        The type of support for the *n-th* solution on the second shape.

        :param int n: The index.

        :return: The support type.
        :rtype: OCCT.BRepExtrema.BRepExtrema_SupportType
        """
        return self._tool.SupportOnShape2(n)

    def is_vertex_shape2(self, n=1):
        """
        Check if support type is a vertex for the second shape.

        :param int n: The index.

        :return: *True* if a vertex, *False* otherwise.
        :rtype: bool
        """
        return self.support_type_shape2(n) == BRepExtrema_IsVertex

    def is_on_edge_shape2(self, n=1):
        """
        Check if support type is on an edge for the second shape.

        :param int n: The index.

        :return: *True* if on an edge, *False* otherwise.
        :rtype: bool
        """
        return self.support_type_shape2(n) == BRepExtrema_IsOnEdge

    def is_in_face_shape2(self, n=1):
        """
        Check if support type is in a face for the second shape.

        :param int n: The index.

        :return: *True* if in a face, *False* otherwise.
        :rtype: bool
        """
        return self.support_type_shape2(n) == BRepExtrema_IsInFace

    def support_on_shape1(self, n=1):
        """
        Get the shape where the *n-th* solution is on the first shape.

        :param int n: The index.

        :return: The support shape.
        :rtype: OCCT.TopoDS.TopoDS_Shape
        """
        return self._tool.SupportOnShape1(n)

    def support_on_shape2(self, n=1):
        """
        Get the shape where the *n-th* solution is on the second shape.

        :param int n: The index.

        :return: The support shape.
        :rtype: OCCT.TopoDS.TopoDS_Shape
        """
        return self._tool.SupportOnShape2(n)

    def par_on_edge_shape1(self, n=1):
        """
        Get the parameter of the *n-th* solution if it is on an edge of the
        first shape.

        :param int n: The index.

        :return: The parameter.
        :rtype: float
        """
        return self._tool.ParOnEdgeS1(n, 0.)

    def par_on_edge_shape2(self, n=1):
        """
        Get the parameter of the *n-th* solution if it is on an edge of the
        second shape.

        :param int n: The index.

        :return: The parameter.
        :rtype: float
        """
        return self._tool.ParOnEdgeS2(n, 0.)

    def par_on_face_shape1(self, n=1):
        """
        Get the parameters of the *n-th* solution if it is in a face of the
        first shape.

        :param int n: The index.

        :return: The parameters.
        :rtype: tuple(float, float)
        """
        return self._tool.ParOnFaceS1(n, 0., 0.)

    def par_on_face_shape2(self, n=1):
        """
        Get the parameters of the *n-th* solution if it is in a face of the
        second shape.

        :param int n: The index.

        :return: The parameters.
        :rtype: tuple(float, float)
        """
        return self._tool.ParOnFaceS2(n, 0., 0.)


class DistanceShapeToShapes(object):
    """
    Calculate the minimum distance between a shape and other shapes. Sort the
    results by distance.

    :param OCCT.TopoDS.TopoDS_Shape shape: The main shape.
    :param list[OCCT.TopoDS.TopoDS_Shape] other_shapes: The other shapes.

    :raises RuntimeWarning: If the distance between two shapes cannot be
        found. This shape will be ignored and process will continue.

    Usage:

    >>> from afem.topology import *
    >>> v1 = VertexByPoint((0., 0., 0.)).vertex
    >>> v2 = VertexByPoint((5., 0., 0.)).vertex
    >>> v3 = VertexByPoint((10., 0., 0.)).vertex
    >>> tool = DistanceShapeToShapes(v1, [v3, v2])
    >>> tool.dmin
    5.0
    >>> tool.dmax
    10.0
    >>> tool.sorted_distances
    [5.0, 10.0]
    """

    def __init__(self, shape, other_shapes):
        results = []
        for shape2 in other_shapes:
            dist = DistanceShapeToShape(shape, shape2)
            if dist.nsol == 0:
                warn("Could not calculate distance to a shape. Continuing...",
                     RuntimeWarning)
                continue
            results.append((dist.dmin, shape2))

        results.sort(key=lambda tup: tup[0])
        self._distances = [data[0] for data in results]
        self._shapes = [data[1] for data in results]

    @property
    def dmin(self):
        """
        :return: The minimum distance of all shapes.
        :rtype: float
        """
        return self._distances[0]

    @property
    def dmax(self):
        """
        :return: The maximum distance of all shapes.
        :rtype: float
        """
        return self._distances[-1]

    @property
    def sorted_distances(self):
        """
        :return: List of sorted distances.
        :rtype: list[float]
        """
        return self._distances

    @property
    def nearest_shape(self):
        """
        :return: The nearest shape.
        :rtype: OCCT.TopoDS.TopoDS_Shape
        """
        return self._shapes[0]

    @property
    def farthest_shape(self):
        """
        :return: The farthest shape.
        :rtype: OCCT.TopoDS.TopoDS_Shape
        """
        return self._shapes[-1]

    @property
    def sorted_shapes(self):
        """
        :return: List of shapes sorted by distance.
        :rtype: list[OCCT.TopoDS.TopoDS_Shape]
        """
        return self._shapes


class DistancePointToShapes(DistanceShapeToShapes):
    """
    Calculate the minimum distance between a point and other shapes. Sort the
    results by distance. This method converts the point to a vertex and then
    uses :class:`.DistanceShapeToShapes`.

    :param point_like pnt: The point.
    :param list[OCCT.TopoDS.TopoDS_Shape] other_shapes: The other shapes.

    :raise TypeError: If *pnt* cannot be converted to a point.
    """

    def __init__(self, pnt, other_shapes):
        pnt = CheckGeom.to_point(pnt)
        if not pnt:
            msg = 'Invalid point type provided.'
            raise TypeError(msg)

        v = CheckShape.to_vertex(pnt)
        super(DistancePointToShapes, self).__init__(v, other_shapes)
