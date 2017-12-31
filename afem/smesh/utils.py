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

from OCCT.SMESH import SMESH_MeshEditor

from afem.smesh.entities import Node

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
            smesh_list.append([n.object for n in row])

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
            smesh_list.append([i for i in row])

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

    # TODO Find free border


class MeshHelper(object):
    # TODO MeshHelper
    pass
