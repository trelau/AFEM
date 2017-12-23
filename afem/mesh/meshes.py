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

from OCCT.SMESH import SMESH_Gen

from afem.mesh.nodes import Node
from afem.topology.check import CheckShape

__all__ = ["MeshGen", "Mesh", "SubMesh"]


class MeshGen(object):
    """
    This class is the primary meshing database for a given instance.
    """

    def __init__(self):
        self._gen = SMESH_Gen()

    @property
    def object(self):
        """
        :return: The underlying mesh object.
        :rtype: OCCT.SMESH.SMESH_Gen
        """
        return self._gen

    @classmethod
    def new(cls, gen):
        """
        Create a new instance using an existing SMESH_Gen instance.

        :param OCCT.SMESH.SMESH_Gen gen: A SMESH_Gen instance.

        :return: The new instance.
        :rtype: afem.mesh.meshes.MeshGen
        """
        new_gen = cls.__new__(cls)
        new_gen._gen = gen
        return new_gen

    def new_id(self):
        """
        Generate a new unique ID within this generator.

        :return: A new unique ID.
        :rtype: int
        """
        return self.object.GetANewId()

    def create_mesh(self, is_embedded=False):
        """
        Create a mesh.

        :param bool is_embedded: Option for embedding mesh.

        :return: The mesh.
        :rtype: afem.mesh.meshes.Mesh
        """
        return Mesh(self, is_embedded)

    def check_algo_state(self, mesh, shape):
        """
        Check if computation would fail because of some bad algorithm state.

        :param afem.mesh.meshes.Mesh mesh: A mesh.
        :param OCCT.TopoDS.TopoDS_Shape shape: A shape.

        :return: *True* if ok, *False* if not.
        :rtype: bool
        """
        return self.object.CheckAlgoState(mesh.object, shape)

    def compute(self, mesh, shape):
        """
        Compute a mesh on a shape.

        :param afem.mesh.meshes.Mesh mesh: A mesh.
        :param OCCT.TopoDS.TopoDS_Shape shape: A shape.

        :return: *True* if computed, *False* if not.
        :rtype: bool
        """
        return self.object.Compute(mesh.object, shape)


class Mesh(object):
    """
    Mesh.

    :param afem.mesh.meshes.MeshGen gen: The MeshGen instance.
    :param bool is_embedded: Option for embedding mesh.
    """

    def __init__(self, gen, is_embedded=False):
        self._mesh = gen.object.CreateMesh(-1, is_embedded)

    @property
    def object(self):
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
        return self.object.GetId()

    @property
    def shape(self):
        """
        :return: The shape to mesh.
        :rtype: OCCT.TopoDS.TopoDS_Shape
        """
        return self.object.GetShapeToMesh()

    @property
    def has_shape(self):
        """
        :return: ``True`` if the mesh has a shape, ``False`` if not.
        :rtype: bool
        """
        return self.object.HasShapeToMesh()

    @classmethod
    def new(cls, mesh):
        """
        Create a new instance using an existing SMESH_Mesh instance.

        :param OCCT.SMESH.SMESH_Mesh mesh: A SMESH_Mesh instance.

        :return: The new instance.
        :rtype: afem.mesh.meshes.Mesh
        """
        new_mesh = cls.__new__(cls)
        new_mesh._mesh = mesh
        return new_mesh

    def shape_to_mesh(self, shape):
        """
        Set the shape to mesh.

        :param OCCT.TopoDS.TopoDS_Shape shape: The shape.

        :return: None
        """
        shape = CheckShape.to_shape(shape)
        self.object.ShapeToMesh(shape)

    def add_hypothesis(self, hypothesis, shape=None):
        """
        Add a hypothesis to the shape.

        :param afem.mesh.hypotheses.Hypothesis hypothesis: The hypothesis.
        :param shape: The shape the hypothesis applies to. This can be a
            sub-shape of the master shape. If not provided then the master
            shape is used.

        :return: None.

        :raise ValueError: If no shape is available to apply the hypothesis to.
        """
        if shape is None:
            if self.has_shape:
                shape = self.shape
            else:
                raise ValueError('No shape could be found.')

        self.object.AddHypothesis(shape, hypothesis.id)

    def clear(self):
        """
        Clear all nodes and elements.

        :return: None.
        """
        self.object.Clear()


class MeshDS(object):
    """
    Mesh data structure.
    """
    # TODO MeshDS


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
        return self.object.IsEmpty()

    @property
    def is_computed(self):
        """
        :return: ``True`` if the sub-mesh is computed, ``False`` if not.
        :rtype: bool
        """
        return self.object.IsMeshComputed()

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
        :rtype: list[afem.mesh.nodes.Node]
        """
        return self.get_nodes()

    def get_nodes(self, include_subshapes=True):
        """
        Get nodes from the sub-mesh.

        :param bool include_subshapes: Option to include sub-shapes when
            retrieving nodes.

        :return: The sub-mesh nodes.
        :rtype: list[afem.mesh.nodes.Node]
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
