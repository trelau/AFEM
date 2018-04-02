# This file is part of AFEM which provides an engineering toolkit for airframe
# finite element modeling during conceptual design.
#
# Copyright (C) 2016-2018  Laughlin Research, LLC (info@laughlinresearch.com)
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
from __future__ import division

from numpy import array, cross, linalg

from afem.geometry.entities import Point

__all__ = ["Node", "Element", "FaceSide"]


class Node(object):
    """
    Mesh node.

    :param OCCT.SMDS.SMDS_MeshNode the_node: The SMDS_MeshNode object.
    """

    def __init__(self, the_node):
        self._node = the_node

    def __str__(self):
        return 'Node {0}: ({1}, {2}, {3})'.format(self.id, *self.xyz)

    def __eq__(self, other):
        return self.id == other.nid

    def __hash__(self):
        return hash(self.id)

    @property
    def object(self):
        """
        :return: The underlying node.
        :rtype: OCCT.SMDS.SMDS_MeshNode
        """
        return self._node

    @property
    def id(self):
        """
        :return: The node ID.
        :rtype: int
        """
        return self._node.GetID()

    @property
    def x(self):
        """
        :return: The node x-location.
        :rtype: float
        """
        return self._node.X()

    @property
    def y(self):
        """
        :return: The node y-location.
        :rtype: float
        """
        return self._node.Y()

    @property
    def z(self):
        """
        :return: The node z-location.
        :rtype: float
        """
        return self._node.Z()

    @property
    def xyz(self):
        """
        :return: The node xyz-location.
        :rtype: numpy.ndarray
        """
        return array([self.x, self.y, self.z], dtype=float)


class Element(object):
    """
    Generic base class for elements.

    :param OCCT.SMDS.SMDS_MeshElement the_element: The element.
    """

    def __init__(self, the_element):
        self._elm = the_element

    def __str__(self):
        eid = 'Element {0}: '.format(str(self.id))
        nids = ' '.join([str(n) for n in self.nids])
        return ''.join([eid, nids])

    def __eq__(self, other):
        return self.id == other.eid

    def __hash__(self):
        return hash(self.id)

    @property
    def object(self):
        """
        :return: The underlying element.
        :rtype: OCCT.SMDS.SMDS_MeshElement
        """
        return self._elm

    @property
    def id(self):
        """
        :return: The element ID.
        :rtype: int
        """
        return self._elm.GetID()

    @property
    def is_quadratic(self):
        """
        :return: *True* if element is quadratic, *False* if not.
        :rtype: bool
        """
        return self._elm.IsQuadratic()

    @property
    def num_nodes(self):
        """
        :return: Number of nodes.
        :rtype: int
        """
        return self._elm.NbNodes()

    @property
    def num_edges(self):
        """
        :return: Number of edges.
        :rtype: int
        """
        return self._elm.NbEdges()

    @property
    def num_faces(self):
        """
        :return: Number of faces.
        :rtype: int
        """
        return self._elm.NbFaces()

    @property
    def num_corner_nodes(self):
        """
        :return: Number of corner nodes.
        :rtype: int
        """
        return self._elm.NbCornerNodes()

    @property
    def node_iter(self):
        """
        :return: Yield nodes of the element. The corner nodes will be first
            followed by medium nodes if quadratic.
        :rtype: collections.Iterable(afem.smesh.entities.Node)
        """
        iter_ = self._elm.nodeIterator()
        while iter_.more():
            yield Node(iter_.next())

    @property
    def point_iter(self):
        """
        :return: Yield nodes of the element as points. The corner points will be
            first followed by medium points if quadratic.
        :rtype: collections.Iterable(afem.geometry.entities.Point)
        """
        for n in self.node_iter:
            yield Point(n.x, n.y, n.z)

    @property
    def edge_iter(self):
        """
        :return: Yield edges of the element.
        :rtype: collections.Iterable(afem.smesh.entities.Element)
        """
        iter_ = self._elm.edgesIterator()
        while iter_.more():
            yield Element(iter_.next())

    @property
    def nids(self):
        """
        :return: The node ID's of the element.
        :rtype: list[int]
        """
        return [n.id for n in self.node_iter]

    @property
    def is_0d(self):
        """
        :return: *True* if the element has one node, *False* if not.
        :rtype: bool
        """
        return self.num_nodes == 1

    @property
    def is_1d(self):
        """
        :return: *True* if the element has two nodes, *False* if not.
        :rtype: bool
        """
        return self.num_nodes == 2

    @property
    def is_2d(self):
        """
        :return: *True* if the element has more than two nodes, *False* if
            not.
        :rtype: bool
        """
        return self.num_nodes > 2

    @property
    def is_tri(self):
        """
        :return: *True* if the element has three nodes, *False* if not.
        :rtype: bool
        """
        return self.num_nodes == 3

    @property
    def is_quad(self):
        """
        :return: *True* if the element has four nodes, *False* if not.
        :rtype: bool
        """
        return self.num_nodes == 4

    @property
    def length(self):
        """
        :return: Element length.
        :rtype: float
        """
        num_pnts = self.num_nodes
        if num_pnts not in [2, 3]:
            return 0.
        pnts = list(self.point_iter)
        d1 = pnts[0].distance(pnts[1])
        if num_pnts == 2:
            return d1
        return d1 + pnts[1].distance(pnts[2])

    @property
    def area(self):
        """
        :return: Element area.
        :rtype: float
        """
        area = 0.
        num_nodes = self.num_nodes
        if self.is_quadratic:
            num_nodes //= 2
        if num_nodes > 2:
            pnts = list(self.point_iter)
            v1 = pnts[1] - pnts[0]
            v2 = pnts[2] - pnts[0]
            vtot = cross(v1, v2)
            for i in range(3, num_nodes):
                v1 = pnts[i - 1] - pnts[0]
                v2 = pnts[i] - pnts[0]
                vtot += cross(v1, v2)
            area = 0.5 * linalg.norm(vtot)

        return area

    @property
    def min_angle(self):
        """
        :return: The minimum element angle in degrees.
        :rtype: float
        """
        raise NotImplementedError()

    @property
    def max_angle(self):
        """
        :return: The maximum element angle in degrees.
        :rtype: float
        """
        raise NotImplementedError()

    @property
    def aspect_ratio(self):
        """
        :return: The element aspect ratio.
        :rtype: float
        """
        raise NotImplementedError()

    @property
    def warp_angle(self):
        """
        :return: The element warping in degrees.
        :rtype: float
        """
        raise NotImplementedError()

    @property
    def taper_ratio(self):
        """
        :return: The element taper ratio.
        :rtype: float
        """
        raise NotImplementedError()

    @property
    def skew_angle(self):
        """
        :return: The element skew angle.
        :rtype: float
        """
        raise NotImplementedError()

    @property
    def jacobian(self):
        """
        :return: The element Jacobian.
        :rtype: float
        """
        raise NotImplementedError()

    def is_medium_node(self, node):
        """
        Check to see if the node is a medium node.

        :param afem.smesh.entities.Node node: A node.

        :return: *True* if medium node, *False* if not.
        :rtype: bool
        """
        return self._elm.IsMediumNode(node.object)

    def node_index(self, node):
        """
        Get the node index if it belongs to the element.

        :param afem.smesh.entities.Node node: A node.

        :return: Node index.
        :rtype: int
        """
        return self._elm.GetNodeIndex(node.object)

    def wrapped_index(self, indx):
        """
        Get a valid node index if given one needs fixed.

        :param int indx: A node index.

        :return: The wrapped node index.
        :rtype: int
        """
        return self._elm.WrappedIndex(indx)

    def get_node(self, indx):
        """
        Get the node at a given index.

        :param int indx: The node index.

        :return: The node.
        :rtype: afem.smesh.entities.Node
        """
        return Node(self._elm.GetNode(indx))

    def get_node_wrap(self, indx):
        """
        Get the node at a given index, fixing it if necessary.

        :param int indx: The node index.

        :return: The node.
        :rtype: afem.smesh.entities.Node
        """
        return Node(self._elm.GetNodeWrap(indx))


class FaceSide(object):
    """
    Entity the represents the side of a quasi-quadrilateral face. It can be
    composed of several edges and gives access to geometry and 1-D mesh.

    :param OCCT.StdMeshers.StdMeshers_FaceSide the_side: The StdMeshers_FaceSide
        instance.
    """

    def __init__(self, the_side):
        self._fside = the_side

    @property
    def object(self):
        """
        :return: The underlying object.
        :rtype: OCCT.StdMeshers.StdMeshers_FaceSide
        """
        return self._fside

    @property
    def num_edges(self):
        """
        :return: The number of edges.
        :rtype: int
        """
        return self._fside.NbEdges()

    @property
    def num_nodes(self):
        """
        :return: The number of nodes.
        :rtype: int
        """
        return self._fside.NbPoints()

    @property
    def num_segments(self):
        """
        :return: The number of mesh segments.
        :rtype: int
        """
        return self._fside.NbSegments()

    @property
    def missed_vertices(self):
        """
        :return: *True* if there are vertices without nodes.
        :rtype: bool
        """
        return self._fside.MissVertexNode()

    @property
    def ordered_nodes(self):
        """
        :return: List of ordered nodes along the side.
        :rtype: list[afem.smesh.entities.Node]
        """
        return [Node(n) for n in self._fside.GetOrderedNodes()]

    @property
    def is_closed(self):
        """
        :return: *True* if chaing of edges is closed.
        :rtype: bool
        """
        return self._fside.IsClosed()

    @property
    def length(self):
        """
        :return: Length of side.
        :rtype: float
        """
        return self._fside.Length()

    @property
    def edges(self):
        """
        :return: List of side edges.
        :rtype: list[OCCT.TopoDS.TopoDS_Edge]
        """
        return self._fside.Edges()

    @property
    def first_vertex(self):
        """
        :return: First vertex of side.
        :rtype: OCCT.TopoDS.TopoDS_Vertex
        """
        return self._fside.FirstVertex()

    @property
    def last_vertex(self):
        """
        :return: Last vertex of side.
        :rtype: OCCT.TopoDS.TopoDS_Vertex
        """
        return self._fside.LastVertex()

    def vertex_node(self, indx):
        """
        Get the node from a vertex.

        :param int indx: The vertex index (starts with 0).

        :return: The vertex node.
        :rtype: afem.smesh.entities.Node
        """
        return Node(self._fside.VertexNode(indx))
