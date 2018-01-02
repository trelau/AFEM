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

from OCCT.SMDS import SMDS_ListOfNodes, SMDS_ListOfElements
from OCCT.SMESH import SMESH_MeshEditor, SMESH_MesherHelper
from OCCT.TopAbs import TopAbs_EDGE

from afem.smesh.entities import Element, Node

__all__ = ["MeshEditor", "MeshHelper"]


class MeshEditor(object):
    """
    Mesh editor.

    :param afem.smesh.meshes.Mesh mesh: A mesh.
    """

    def __init__(self, mesh):
        self._editor = SMESH_MeshEditor(mesh.object)

    @property
    def object(self):
        """
        :return: The underlying mesh editor.
        :rtype: OCCT.SMESH.SMESH_MeshEditor
        """
        return self._editor

    def reorient(self, elm):
        """
        Reverse the element orientation.

        :param afem.smesh.entities.Element elm: The element.

        :return: *True* if reoriented, *False* if not.
        :rtype: bool
        """
        return self._editor.Reorient(elm.object)

    def smooth(self, elms=(), fixed_nodes=(), method='laplacian', iters=5,
               target_ar=1., in_2d=True):
        """
        Smooth the mesh.

        :param collections.Sequence(afem.smesh.entities.Element) elms: The
            elements to smooth. If empty then the whole mesh is smoothed.
        :param collections.Sequence(afem.smesh.entities.Node) fixed_nodes: The
            fixed nodes if any. Nodes on edges and boundaries are always fixed.
        :param str method: Either 'laplacian', 'l', 'centroidal', or 'c'.
        :param int iters: Number of iterations.
        :param float target_ar: Target aspect ratio.
        :param bool in_2d: Perform smoothing in 2-d parameters of nodes on face.

        :return: None.
        """
        elms = set([e.object for e in elms])
        fixed_nodes = set([n.object for n in fixed_nodes])
        if method.lower() in ['c', 'centroidal']:
            method = self._editor.CENTROIDAL
        else:
            method = self._editor.LAPLACIAN

        self._editor.Smooth(elms, fixed_nodes, method, iters, target_ar, in_2d)

    def find_coincident_nodes(self, nodes=(), tol=1.0e-7):
        """
        Find coincident nodes.

        :param collections.Sequence(afem.smesh.entities.Node) nodes: The nodes
            to search. If not provided then the whole mesh will be searched.
        :param float tol: Search tolerance.

        :return: Coincident nodes as a list of list of nodes. The length of the
            returned list will be the number of coincident nodes. Each row
            is a list of nodes that are coincident.
        :rtype: list[list[afem.smesh.entities.Node]]
        """
        nodes = set([n.object for n in nodes])
        smesh_list = self._editor.TListOfListOfNodes()
        self._editor.FindCoincidentNodes(nodes, tol, smesh_list, False)

        coincident_nodes = []
        for row in smesh_list:
            nodes = [Node(n) for n in row]
            coincident_nodes.append(nodes)

        return coincident_nodes

    def merge_nodes(self, nodes=None, avoid_making_holes=False, tol=1.0e-7):
        """
        Merge nodes.

        :param nodes: The nodes to merge. Each row is a list of nodes where the
            first node is kept and the others are replaced with the first. If
            *None* is provided then the entire mesh is first searched.
        :type nodes: list[list[afem.smesh.entities.Node]]
        :param bool avoid_making_holes: Avoid modifications that may spoil mesh
            topology.
        :param float tol: Search tolerance.

        :return: None.
        """
        if nodes is None:
            nodes = self.find_coincident_nodes(tol=tol)

        smesh_list = self._editor.TListOfListOfNodes()
        for row in nodes:
            smesh_list.push_back([n.object for n in row])

        self._editor.MergeNodes(smesh_list, avoid_making_holes)

    def find_equal_elements(self, elements=()):
        """
        Find equal elements.

        :param collections.Sequence(afem.smesh.entities.Element) elements: The
            elements to search. If not provided then the whole mesh will be
            searched.

        :return: Equal elements as a list of list of integers. The length of the
            returned list will be the number of equal elements. Each row
            is a list of element ID's that are equal.
        :rtype: list[list[int]]
        """
        elements = set([e.object for e in elements])
        smesh_list = self._editor.TListOfListOfElementsID()
        self._editor.FindEqualElements(elements, smesh_list)

        equal_ids = []
        for row in smesh_list:
            elements = [i for i in row]
            equal_ids.append(elements)

        return equal_ids

    def merge_elements(self, elements=None):
        """
        Merge elements.

        :param elements: The elements to merge. Each row is a list of element
            ID's where the first ID is kept and the others are replaced with the
            first. If *None* is provided then the entire mesh is first searched.
        :type elements: list[list[int]]

        :return: None.
        """
        if elements is None:
            elements = self.find_equal_elements()

        smesh_list = self._editor.TListOfListOfElementsID()
        for row in elements:
            smesh_list.push_back([i for i in row])

        self._editor.MergeElements(smesh_list)

    def merge_equal_elements(self):
        """
        Merge elements built on same nodes.

        :return: None.
        """
        self._editor.MergeEqualElements()

    def find_shape(self, elm):
        """
        Find the shape index that the element or node is on.

        :param elm: The node or element.
        :type elm: afem.smesh.entities.Node or afem.smesh.entities.Element

        :return: The shape index.
        :rtype: int
        """
        return self._editor.FindShape(elm.object)

    def is_medium(self, node):
        """
        Check if node is a medium.

        :param afem.smesh.entities.Node node: The node.

        :return: *True* if medium, *False* if not.
        :rtype: bool
        """
        return self._editor.IsMedium_(node.object)

    def double_elements(self, elms):
        """
        Double elements.

        :param collections.Sequence(afem.smesh.entities.Element) elms: The
            elements to double.

        :return: None.
        """
        elms = set([e.object for e in elms])
        self._editor.DoubleElements(elms)

    def double_nodes(self, nids):
        """
        Double nodes.

        :param collections.Sequence(int) nids: Sequence of node ID's.

        :return: *True* if doubled, *False* if not.
        :rtype: bool
        """
        nids = list(nids)
        return self._editor.DoubleNodes(nids)

    def transform(self, trsf, elements=(), copy=False, make_groups=False,
                  target_mesh=None):
        """
        Transform the elements.

        :param OCCT.gp.gp_Trsf trsf: The transformation.
        :param collection.Sequence(afem.smesh.entities.Element) elements: The
            elements to transform. If none are provided then the whole mesh is
            used.
        :param bool copy: Option to copy elements.
        :param bool make_groups: Option to make groups.
        :param afem.smesh.meshes.Mesh target_mesh: The target mesh to place
            elements.

        :return: List of element ID's.
        :rtype: list[int]
        """
        elements = set([e.object for e in elements])
        if target_mesh is None:
            target_mesh = self._editor.GetMesh()
        else:
            target_mesh = target_mesh.object
        # FIXME Why returning None and not actual type?
        if copy:
            return self._editor.Transform(elements, trsf, copy, make_groups,
                                          target_mesh)
        else:
            return self._editor.Transform(elements, trsf, copy, make_groups)

    def check_free_border(self, n1, n2, n3=None):
        """
        Check to see if the nodes are on a free border.

        :param afem.smesh.entities.Node n1: The first node.
        :param afem.smesh.entities.Node n2: The second node.
        :param afem.smesh.entities.Node n3: The third node.

        :return: *True* if on a free border, *False* if not.
        :rtype: bool
        """
        if n3 is None:
            return self._editor.CheckFreeBorderNodes_(n1.object, n2.object)
        else:
            return self._editor.CheckFreeBorderNodes_(n1.object, n2.object,
                                                      n3.object)

    def find_free_border(self, n1, n2, n3):
        """
        Return nodes and faces of a free border if found.

        :param afem.smesh.entities.Node n1: The first node.
        :param afem.smesh.entities.Node n2: The second node.
        :param afem.smesh.entities.Node n3: The third node.

        :return: List of nodes and elements of free borders.
        :rtype: tuple(list[afem.smesh.entities.Node],
            list[afem.smesh.entities.Element])
        """
        node_list = SMDS_ListOfNodes()
        elm_list = SMDS_ListOfElements()
        return self._editor.FindFreeBorder_(n1.object, n2.object, n3.object,
                                            node_list, elm_list)


class MeshHelper(object):
    """
    Mesh helper.

    :param afem.smesh.meshes.Mesh mesh: A mesh.
    """

    def __init__(self, mesh):
        self._helper = SMESH_MesherHelper(mesh.object)

    @property
    def object(self):
        """
        :return: The underlying object.
        :rtype: OCCT.SMESH.SMESH_MesherHelper
        """
        return self._helper

    @staticmethod
    def is_structured(submesh):
        """
        Check if a 2-D mesh on a face is structured.

        :param afem.smesh.meshes.SubMesh submesh: The submesh.

        :return: *True* if structured, *False* if not.
        :rtype: bool
        """
        return SMESH_MesherHelper.IsStructured_(submesh.object)

    @staticmethod
    def is_distorted2d(submesh, check_uv=False):
        """
        Check if a 2-D mesh on a face is distorted.

        :param afem.smesh.meshes.SubMesh submesh: The submesh.
        :param bool check_uv: Option to check using *uv* parameters.

        :return: *True* if distored, *False* if not.
        :rtype: bool
        """
        return SMESH_MesherHelper.IsDistorted2D_(submesh.object, check_uv)

    @staticmethod
    def subshape_by_node(node, mesh_ds):
        """
        Get the support shape of a node.

        :param afem.smesh.entities.Node node: The node.
        :param afem.smesh.meshes.MeshDS mesh_ds: The mesh data structure.

        :return: The support shape.
        :rtype: OCCT.TopoDS.TopoDS_Shape
        """
        return SMESH_MesherHelper.GetSubShapeByNode_(node.object,
                                                     mesh_ds.object)

    @staticmethod
    def common_ancestor(shape1, shape2, mesh, ancestor_type=TopAbs_EDGE):
        """
        Get a common ancestor between the two shapes.

        :param OCCT.TopoDS.TopoDS_Shape shape1: The first shape.
        :param OCCT.TopoDS.TopoDS_Shape shape2: The second shape.
        :param afem.smesh.meshes.Mesh mesh: The mesh.
        :param OCCT.TopAbs.TopAbs_ShapeEnum ancestor_type: The shape type.

        :return: The common ancestor.
        :rtype: OCCT.TopoDS.TopoDS_Edge
        """
        return SMESH_MesherHelper.GetCommonAncestor_(shape1, shape2,
                                                     mesh.object, ancestor_type)

    @staticmethod
    def is_subshape_by_shape(shape, main_shape):
        """
        Check to see if the shape is a sub-shape in a shape.

        :param OCCT.TopoDS.TopoDS_Shape shape: The shape.
        :param OCCT.TopoDS.TopoDS_Shape main_shape: The main shape.

        :return: *True* if a sub-shape, *False* if not.
        :rtype: bool
        """
        return SMESH_MesherHelper.IsSubShape_(shape, main_shape)

    @staticmethod
    def is_subshape_by_mesh(shape, mesh):
        """
        Check to see if the shape is a sub-shape in a mesh.

        :param OCCT.TopoDS.TopoDS_Shape shape: The shape.
        :param afem.smesh.meshes.Mesh mesh: The mesh.

        :return: *True* if a sub-shape, *False* if not.
        :rtype: bool
        """
        return SMESH_MesherHelper.IsSubShape_(shape, mesh.object)

    @staticmethod
    def get_angle(e1, e2, f, v):
        """
        Determine the angle between the two edges on a face.

        :param OCCT.TopoDS.TopoDS_Edge e1: The first edge.
        :param OCCT.TopoDS.TopoDS_Edge e2: The second edge.
        :param OCCT.TopoDS.TopoDS_Face f: The face the edges belong to.
        :param OCCT.TopoDS.TopoDS_Vertex v: The vertex connecting the two edges.

        :return: The angle between the egdges.
        :rtype: float
        """
        return SMESH_MesherHelper.GetAngle_(e1, e2, f, v)

    @staticmethod
    def is_closed_edge(edge):
        """
        Check if edge is closed.

        :param OCCT.TopoDS.TopoDS_Edge edge: The edge.

        :return: *True* if closed, *False* if not.
        :rtype: bool
        """
        return SMESH_MesherHelper.IsClosedEdge_(edge)

    @staticmethod
    def shape_by_hypothesis(hyp, shape, mesh):
        """
        Get a shape the hypothesis is applied to.

        :param afem.smesh.hypotheses.Hypothesis hyp: The hypothesis.
        :param OCCT.TopoDS.TopoDS_Shape shape: The shape.
        :param afem.smesh.meshes.Mesh mesh: The mesh.

        :return: The shape that the hypothesis is applied to.
        :rtype: OCCT.TopoDS.TopoDS_Shape
        """
        return SMESH_MesherHelper.GetShapeOfHypothesis_(hyp.object, shape,
                                                        mesh.object)

    def is_reversed_submesh(self, face):
        """
        Check to see if the elements have opposite orientation on the face.

        :param OCCT.TopoDS.TopoDS_Face face: The face.

        :return: *True* if reversed, *False* if not.
        :rtype: bool
        """
        return self._helper.IsReversedSubMesh(face)

    def set_subshape(self, shape):
        """
        Set the shape to make elements on.

        :param shape: The shape or the shape ID.
        :type shape: OCCT.TopoDS.TopoDS_Shape or int

        :return: None.
        """
        self._helper.SetSubShape(shape)

    def shape_to_index(self, shape):
        """
        Get the shape index.

        :param OCCT.TopoDS.TopoDS_Shape shape: The shape.

        :return: The shape index.
        :rtype: int
        """
        return self._helper.ShapeToIndex(shape)

    def add_node(self, x, y, z, id_=0, u=0., v=0.):
        """
        Create a node.

        :param float x: The x-location.
        :param float y: The y-location.
        :param float z: The z-location.
        :param int id_: The node ID. If zero then the ID will be assigned.
        :param float u: The node u-parameter.
        :param float v: The node v-parameter.

        :return: The created node.
        :rtype: afem.smesh.entities.Node
        """
        smesh_node = self._helper.AddNode(x, y, z, id_, u, v)
        return Node(smesh_node)

    def add_edge(self, n1, n2, id_=0, force3d=True):
        """
        Add an edge.

        :param afem.smesh.entities.Node n1: The first node.
        :param afem.smesh.entities.Node n2: The second node.
        :param int id_: The edge ID. If zero then the ID will be assigned.
        :param force3d: Unknown option.

        :return: The created edge.
        :rtype: afem.smesh.entities.Element
        """
        smesh_elm = self._helper.AddEdge(n1.object, n2.object, id_, force3d)
        return Element(smesh_elm)

    def add_face(self, n1, n2, n3, n4=None, id_=0, force3d=False):
        """
        Add a face.

        :param afem.smesh.entities.Node n1: The first node.
        :param afem.smesh.entities.Node n2: The second node.
        :param afem.smesh.entities.Node n3: The third node.
        :param afem.smesh.entities.Node n4: The fourth node. If provided then
            the face will be a quadrangle, if not provided then the face is
            a triangle.
        :param int id_: The face ID. If zero then the ID will be assigned.
        :param force3d: Unknown option.

        :return: The created face.
        :rtype: afem.smesh.entities.Element
        """
        if n4 is None:
            smesh_elm = self._helper.AddFace(n1.object, n2.object, n3.object,
                                             id_, force3d)
        else:
            smesh_elm = self._helper.AddFace(n1.object, n2.object, n3.object,
                                             n4.object, id_, force3d)
        return Element(smesh_elm)
