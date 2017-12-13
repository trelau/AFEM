#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2017 Laughlin Research, L.L.C.
#
# This file is subject to the license agreement that was delivered
# with this source code.
#
# THE SOFTWARE AND INFORMATION ARE PROVIDED ON AN AS "AS IS" BASIS,
# WITHOUT ANY WARRANTIES OR REPRESENTATIONS EXPRESS, IMPLIED OR 
# STATUTORY; INCLUDING, WITHOUT LIMITATION, WARRANTIES OF QUALITY,
# PERFORMANCE, MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

from afem.fem.elements import Element
from afem.fem.hypotheses import HypothesisAPI, the_gen
from afem.fem.nodes import Node
from afem.topology.check import CheckShape

__all__ = ["Mesh", "SubMesh", "MeshAPI"]


class Mesh(object):
    """
    Mesh.

    :param str label: The name.
    """
    _all = {}
    _indx = 0

    def __init__(self, label):
        self._label = label
        Mesh._all[label] = self
        self._mesh = the_gen.CreateMesh(Mesh._indx, True)
        self._ds = self._mesh.GetMeshDS()
        self._id = Mesh._indx
        Mesh._indx += 1

    @property
    def handle(self):
        """
        :return: The underlying mesh object.
        :rtype: OCCT.SMESH.SMESH_Mesh
        """
        return self._mesh

    @property
    def id(self):
        """
        :return: The mesh ID.
        :rtype: int
        """
        return self._id

    @property
    def shape(self):
        """
        :return: The shape to mesh.
        :rtype: OCCT.TopoDS.TopoDS_Shape
        """
        return self.handle.GetShapeToMesh()

    @property
    def has_shape(self):
        """
        :return: ``True`` if the mesh has a shape, ``False`` if not.
        :rtype: bool
        """
        return self.handle.HasShapeToMesh()

    @property
    def nb_nodes(self):
        """
        :return: The number of nodes in the mesh.
        :rtype: int
        """
        return self._ds.NbNodes()

    @property
    def min_node_id(self):
        """
        :return: The minimum node ID in the mesh.
        :rtype: int
        """
        return self._ds.MinNodeID()

    @property
    def max_node_id(self):
        """
        :return: The maximum node ID in the mesh.
        :rtype: int
        """
        return self._ds.MaxNodeID()

    @property
    def min_elm_id(self):
        """
        :return: The minimum element ID in the mesh.
        :rtype: int
        """
        return self._ds.MinElementID()

    @property
    def max_elm_id(self):
        """
        :return: The maximum element ID in the mesh.
        :rtype: int
        """
        return self._ds.MaxElementID()

    @property
    def nodes(self):
        """
        :return: The mesh nodes.
        :rtype: list[afem.fem.nodes.Node]
        """
        return self.get_nodes()

    @property
    def elements(self):
        """
        :return: The mesh elements.
        :rtype: list[afem.fem.elements.Element]
        """
        elm_iter = self._ds.elementsIterator()
        elms = []
        while elm_iter.more():
            elm = Element(elm_iter.next())
            elms.append(elm)
        return elms

    @classmethod
    def get_mesh(cls, mesh):
        """
        Get a mesh.

        :param mesh: The mesh to get. If a mesh instance is given it is
            simply returned. If a string is given the mesh is retrieved by
            its label.
        :type mesh: afem.fem.meshes.Mesh or str

        :return: The mesh.
        :rtype: afem.fem.meshes.Mesh

        :raise KeyError: If a mesh cannot be found.
        """
        if isinstance(mesh, Mesh):
            return mesh

        return Mesh._all[mesh]

    def activate(self):
        """
        Activate this mesh.

        :return: None.
        """
        MeshAPI._active = self

    def shape_to_mesh(self, shape):
        """
        Set the shape to mesh.

        :param OCCT.TopoDS.TopoDS_Shape shape: The shape.

        :return: None
        """
        shape = CheckShape.to_shape(shape)
        self.handle.ShapeToMesh(shape)

    def add_hypothesis(self, hypothesis, shape=None):
        """
        Add a hypothesis to the shape.

        :param hypothesis: The hypothesis. Can be an instance or a label.
        :type hypothesis: afem.fem.hypotheses.Hypothesis or str
        :param shape: The shape the hypothesis applies to. This can be a
            sub-shape of the master shape. If not provided then the master
            shape is used.

        :return: None.

        :raise ValueError: If no shape is available to apply the hypothesis to.
        :raise ValueError: If no hypothesis is found.
        """
        shape = CheckShape.to_shape(shape)
        if not shape:
            if self.has_shape:
                shape = self.shape
            else:
                raise ValueError('No shape could be found.')

        hypothesis = HypothesisAPI.get_hypothesis(hypothesis)
        if not hypothesis:
            raise ValueError('No hypothesis could be found.')

        self.handle.AddHypothesis(shape, hypothesis.id)

    def compute(self):
        """
        Compute the mesh.

        :return: ``True`` if performed, ``False`` if not.
        :rtype: bool
        """
        return the_gen.Compute(self.handle, self.shape)

    def clear(self):
        """
        Clear all nodes and elements.

        :return: None.
        """
        self.handle.Clear()

    def get_submesh(self, sub_shape):
        """
        Get a SubMesh from a sub-shape.

        :param OCCT.TopoDS.TopoDS_Shape sub_shape: The sub-shape.

        :return: A sub-mesh.
        :rtype: afem.fem.meshes.SubMesh
        """
        the_mesh = self.handle.GetSubMesh(sub_shape)
        return SubMesh(the_mesh)

    def get_nodes(self, order=False):
        """
        Get the nodes of the mesh.

        :param bool order: Order nodes by their ID.

        :return: The nodes.
        :rtype: list[afem.fem.nodes.Node]
        """
        niter = self._ds.nodesIterator(order)
        nodes = []
        while niter.more():
            n = Node(niter.next())
            nodes.append(n)
        return nodes


class SubMesh(object):
    """
    SubMesh.

    :param OCCT.SMESH.SMESH_subMesh: The sub-mesh.
    """

    def __init__(self, the_submesh):
        self._mesh = the_submesh
        self._ds = the_submesh.GetSubMeshDS()

    @property
    def object(self):
        """
        :return: The underlying sub-mesh object.
        :rtype: OCCT.SMESH.SMESH_subMesh
        """
        return self._mesh

    @property
    def is_empty(self):
        """
        :return: ``True`` if the sub-mesh is empty, ``False`` if not.
        :rtype: bool
        """
        return self.handle.IsEmpty()

    @property
    def is_computed(self):
        """
        :return: ``True`` if the sub-mesh is computed, ``False`` if not.
        :rtype: bool
        """
        return self.handle.IsMeshComputed()

    @property
    def nb_nodes(self):
        """
        :return: The number of nodes in the sub-mesh.
        :rtype: int
        """
        return self._ds.NbNodes()

    @property
    def nodes(self):
        """
        :return: The sub-mesh nodes.
        :rtype: list[afem.fem.nodes.Node]
        """
        return self.get_nodes()

    def get_nodes(self, include_subshapes=True):
        """
        Get nodes from the sub-mesh.

        :param bool include_subshapes: Option to include sub-shapes when
            retrieving nodes.

        :return: The sub-mesh nodes.
        :rtype: list[afem.fem.nodes.Node]
        """
        # Return nodes on sub-shape only.
        if not include_subshapes:
            niter = self._ds.GetNodes()
            nodes = []
            while niter.more():
                n = Node(niter.next())
                nodes.append(n)
            return nodes

        # Can't figure out how to get nodes from sub-shapes, so get the
        # elements and then the nodes from them.
        # TODO Figure out how to access sub-nodes.
        elm_iter = self._ds.GetElements()
        node_set = set()
        while elm_iter.more():
            elm = elm_iter.next()
            niter = elm.nodeIterator()
            while niter.more():
                n = Node(niter.next())
                node_set.add(n)
        return list(node_set)


class MeshAPI(object):
    """
    Mesh API. This is used to manage meshes from one place.

    :var afem.fem.hypotheses.HypothesisAPI hypotheses: Access to
        HypothesisAPI interface.
    """
    _active = None
    hypotheses = HypothesisAPI()

    @classmethod
    def get_active(cls):
        """
        Get the active mesh.

        :return: The active mesh.
        :rtype: afem.fem.meshes.Mesh
        """
        return cls._active

    @classmethod
    def get_mesh(cls, mesh=None):
        """
        Get mesh.

        :param mesh: The mesh to get. If a mesh instance is given it is
            simply returned. If a string is given the mesh is retrieved by
            its label. If ``None`` is given then the active mesh is returned.
        :type mesh: afem.fem.meshes.Mesh or str or None

        :return: The mesh.
        :rtype: afem.fem.meshes.Mesh
        """
        if mesh is None:
            return cls._active

        return Mesh.get_mesh(mesh)

    @classmethod
    def make_active(cls, mesh):
        """
        Activate the mesh.

        :param mesh: The mesh to activate. If a mesh instance is given it is
            activated. If a string is given the mesh is retrieved by
            its label.
        :type mesh: afem.fem.meshes.Mesh or str

        :return: None.
        """
        mesh = cls.get_mesh(mesh)
        mesh.activate()

    @classmethod
    def create_mesh(cls, label, shape, active=True):
        """
        Create a mesh.

        :param str label: The label.
        :param OCCT.TopoDS.TopoDS_Shape shape: The shape.
        :param bool active: Option to make this mesh active.

        :return: The new mesh.
        :rtype: afem.fem.meshes.Mesh
        """
        mesh = Mesh(label)
        mesh.shape_to_mesh(shape)
        if active:
            mesh.activate()
        return mesh

    @classmethod
    def add_hypothesis(cls, hypothesis, shape=None, mesh=None):
        """
        Add a hypothesis to the shape in a mesh.

        :param hypothesis: The hypothesis. Can be an instance or a label.
        :type hypothesis: afem.fem.hypotheses.Hypothesis or str
        :param shape: The shape the hypothesis applies to. This can be a
            sub-shape of the master shape. If not provided then the master
            shape is used.
        :param mesh: The mesh to apply the hypothesis to. If ``None`` is
            provided then the active mesh is used.
        :type mesh: afem.fem.meshes.Mesh or str or None

        :return: None.
        """
        mesh = cls.get_mesh(mesh)
        mesh.add_hypothesis(hypothesis, shape)

    @classmethod
    def compute_mesh(cls, mesh=None):
        """
        Compute a mesh.

        :param mesh: The mesh to compute. If ``None`` is provided then the
            active mesh is used.
        :type mesh: afem.fem.meshes.Mesh or str or None

        :return: ``True`` if performed, ``False`` if not.
        :rtype: bool
        """
        mesh = cls.get_mesh(mesh)
        return mesh.compute()

    @classmethod
    def get_submesh(cls, sub_shape, mesh=None):
        """
        Get a SubMesh from the sub-shape.

        :param OCCT.TopoDS.TopoDS_Shape sub_shape: The sub-shape.
        :param mesh: The mesh to retrieve the sub-mesh from. If ``None`` is
            provided then the active mesh is used.
        :type mesh: afem.fem.meshes.Mesh or str or None

        :return: The sub-mesh.
        :rtype: afem.fem.meshes.SubMesh
        """
        mesh = cls.get_mesh(mesh)
        return mesh.get_submesh(sub_shape)

    @classmethod
    def get_nodes(cls, order=False, mesh=None):
        """
        Get nodes from the mesh.

        :param bool order: Order nodes order by ID.
        :param mesh: The mesh to retrieve the nodes from. If ``None`` is
            provided then the active mesh is used.
        :type mesh: afem.fem.meshes.Mesh or str or None

        :return: The nodes.
        :rtype: list[afem.fem.nodes.Node]
        """
        mesh = cls.get_mesh(mesh)
        return mesh.get_nodes(order)
