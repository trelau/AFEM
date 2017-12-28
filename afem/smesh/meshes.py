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

from OCCT.SMESH import SMESH_Gen, SMESH_subMesh

from afem.smesh.entities import Node, Element
from afem.topology.check import CheckShape

__all__ = ["MeshGen", "Mesh", "MeshDS", "SubMesh", "SubMeshDS"]


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
        :rtype: afem.smesh.meshes.MeshGen
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
        return self._gen.GetANewId()

    def create_mesh(self, shape=None, is_embedded=False):
        """
        Create a mesh.

        :param OCCT.TopoDS.TopoDS_Shape shape: The shape to mesh.
        :param bool is_embedded: Option for embedding mesh.

        :return: The mesh.
        :rtype: afem.smesh.meshes.Mesh
        """
        the_mesh = Mesh(self, is_embedded)
        if shape is not None:
            the_mesh.shape_to_mesh(shape)
        return the_mesh

    def check_algo_state(self, mesh, shape):
        """
        Check if computation would fail because of some bad algorithm state.

        :param afem.smesh.meshes.Mesh mesh: A mesh.
        :param OCCT.TopoDS.TopoDS_Shape shape: A shape.

        :return: *True* if ok, *False* if not.
        :rtype: bool
        """
        return self._gen.CheckAlgoState(mesh.object, shape)

    def compute(self, mesh, shape):
        """
        Compute a mesh on a shape.

        :param afem.smesh.meshes.Mesh mesh: A mesh.
        :param OCCT.TopoDS.TopoDS_Shape shape: A shape.

        :return: *True* if computed, *False* if not.
        :rtype: bool
        """
        return self._gen.Compute(mesh.object, shape)


class Mesh(object):
    """
    Mesh.

    :param afem.smesh.meshes.MeshGen gen: The MeshGen instance.
    :param bool is_embedded: Option for embedding mesh.
    """

    def __init__(self, gen, is_embedded=False):
        self._mesh = gen.object.CreateMesh(-1, is_embedded)
        self._ds = MeshDS(self)

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
        return self._mesh.GetId()

    @property
    def shape(self):
        """
        :return: The shape to mesh.
        :rtype: OCCT.TopoDS.TopoDS_Shape
        """
        return self._mesh.GetShapeToMesh()

    @property
    def has_shape(self):
        """
        :return: ``True`` if the mesh has a shape, ``False`` if not.
        :rtype: bool
        """
        return self._mesh.HasShapeToMesh()

    @property
    def num_nodes(self):
        """
        :return: Number of nodes.
        :rtype: int
        """
        return self._mesh.NbNodes()

    @property
    def num_edges(self):
        """
        :return: Number of edges.
        :rtype: int
        """
        return self._mesh.NbEdges()

    @property
    def num_faces(self):
        """
        :return: Number of faces.
        :rtype: int
        """
        return self._mesh.NbFaces()

    @property
    def num_tris(self):
        """
        :return: Number of triangles.
        :rtype: int
        """
        return self._mesh.NbTriangles()

    @property
    def num_quads(self):
        """
        :return: Number of quadrangles.
        :rtype: int
        """
        return self._mesh.NbQuadrangles()

    @property
    def num_submesh(self):
        """
        :return: Number of sub-meshes.
        :rtype: int
        """
        return self._mesh.NbSubMesh()

    @property
    def num_group(self):
        """
        :return: Number of groups.
        :rtype: int
        """
        return self._mesh.NbGroup()

    @property
    def ds(self):
        """
        :return: The mesh data structure.
        :rtype: afem.smesh.meshes.MeshDS
        """
        return self._ds

    @classmethod
    def new(cls, mesh):
        """
        Create a new instance using an existing SMESH_Mesh instance.

        :param OCCT.SMESH.SMESH_Mesh mesh: A SMESH_Mesh instance.

        :return: The new instance.
        :rtype: afem.smesh.meshes.Mesh
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
        self._mesh.ShapeToMesh(shape)

    def add_hypothesis(self, hypothesis, shape=None):
        """
        Add a hypothesis to the shape.

        :param afem.smesh.hypotheses.Hypothesis hypothesis: The hypothesis.
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

        self._mesh.AddHypothesis(shape, hypothesis.id)

    def clear(self):
        """
        Clear all nodes and elements.

        :return: None.
        """
        self._mesh.Clear()

    def clear_submesh(self, shape):
        """
        Clear the nodes and elements from the shape.

        :param OCCT.TopoDS.TopoDS_Shape shape: The shape.

        :return: None.
        """
        shape_id = self.ds.shape_to_index(shape)
        self._mesh.ClearSubMesh(shape_id)

    def get_submesh(self, sub_shape):
        """
        Get a sub-mesh for the sub-shape.

        :param OCCT.TopoDS.TopoDS_Shape sub_shape: A sub-shape.

        :return: A sub-mesh.
        :rtype: afem.smesh.meshes.SubMesh
        """
        sub_mesh = self._mesh.GetSubMesh(sub_shape)
        return SubMesh.new(sub_mesh)

    def get_submesh_containing(self, sub_shape):
        """
        Get a sub-mesh for containing the sub-shape.

        :param OCCT.TopoDS.TopoDS_Shape sub_shape: A sub-shape.

        :return: A sub-mesh.
        :rtype: afem.smesh.meshes.SubMesh
        """
        sub_mesh = self._mesh.GetSubMeshContaining(sub_shape)
        return SubMesh.new(sub_mesh)


class MeshDS(object):
    """
    Mesh data structure.

    :param afem.smesh.meshes.Mesh mesh: A mesh.
    """

    def __init__(self, mesh):
        self._ds = mesh.object.GetMeshDS()

    @property
    def object(self):
        """
        :return: The underlying mesh object.
        :rtype: OCCT.SMESH.SMESHDS_Mesh
        """
        return self._ds

    @property
    def is_embedded_mode(self):
        """
        :return: *True* if embedded, *False* if not.
        :rtype: bool
        """
        return self._ds.IsEmbeddedMode()

    @property
    def id(self):
        """
        :return: The mesh ID.
        :rtype: int
        """
        return self._ds.GetPersistentId()

    @property
    def min_node_id(self):
        """
        :return: The minimum node id.
        :rtype: int
        """
        return self._ds.MinNodeID()

    @property
    def max_node_id(self):
        """
        :return: The maximum node id.
        :rtype: int
        """
        return self._ds.MaxNodeID()

    @property
    def min_elm_id(self):
        """
        :return: The minimum element id.
        :rtype: int
        """
        return self._ds.MinElementID()

    @property
    def max_elm_id(self):
        """
        :return: The maximum element id.
        :rtype: int
        """
        return self._ds.MaxElementID()

    @property
    def num_nodes(self):
        """
        :return: Number of nodes.
        :rtype: int
        """
        return self._ds.NbNodes()

    @property
    def num_elms(self):
        """
        :return: Number of elements.
        :rtype: int
        """
        return self._ds.NbElements()

    @property
    def num_edges(self):
        """
        :return: Number of edges.
        :rtype: int
        """
        return self._ds.NbEdges()

    @property
    def num_faces(self):
        """
        :return: Number of faces.
        :rtype: int
        """
        return self._ds.NbFaces()

    @property
    def node_iter(self):
        """
        :return: Yield nodes of the mesh.
        :rtype: collections.Iterable(afem.smesh.entities.Node)
        """
        iter_ = self._ds.nodesIterator(True)
        while iter_.more():
            yield Node(iter_.next())

    @property
    def edge_iter(self):
        """
        :return: Yield edges of the mesh.
        :rtype: collections.Iterable(afem.smesh.entities.Element)
        """
        iter_ = self._ds.edgesIterator(True)
        while iter_.more():
            yield Element(iter_.next())

    @property
    def faces_iter(self):
        """
        :return: Yield faces of the mesh.
        :rtype: collections.Iterable(afem.smesh.entities.Element)
        """
        iter_ = self._ds.facesIterator(True)
        while iter_.more():
            yield Element(iter_.next())

    def move_node(self, node, x, y, z):
        """
        Move node to given location.

        :param afem.smesh.nodes.Node node: The node.
        :param float x: The x-location.
        :param float y: The y-location.
        :param float z:  The z-location.

        :return: None.
        """
        self._ds.MoveNode(node.object, x, y, z)

    def renumber_nodes(self, start=1, step=1):
        """
        Renumber nodes.

        :param int start: Starting node id.
        :param int step: Step between id's.

        :return: None.
        """
        self._ds.Renumber(True, start, step)

    def renumber_elements(self, start=1, step=1):
        """
        Renumber elements.

        :param int start: Starting element id.
        :param int step: Step between id's.

        :return: None.
        """
        self._ds.Renumber(False, start, step)

    def has_elements(self, shape):
        """
        Check to see if the shape has mesh elements.

        :param OCCT.TopoDS.TopoDS_Shape shape: The shape.

        :return: *True* if the shape has elements, *False* if not.
        :rtype: bool
        """
        return self._ds.HasMeshElements(shape)

    def mesh_elements(self, shape):
        """
        Get elements of shape.

        :param OCCT.TopoDS.TopoDS_Shape shape: The shape.

        :return:
        :rtype:
        """
        sub_meshds = self._ds.MeshElements(shape)
        return SubMeshDS.new(sub_meshds)

    def shape_to_index(self, shape):
        """
        Get the shape index.

        :param OCCT.TopoDS.TopoDS_Shape shape: The shape.

        :return: The shape index.
        :rtype: int
        """
        return self._ds.ShapeToIndex(shape)

    def index_to_shape(self, indx):
        """
        Get the shape from an index.

        :param int indx: The index.

        :return: The shape.
        :rtype: OCCT.TopoDS.TopoDS_Shape
        """
        return self._ds.IndexToShape(indx)


class SubMesh(object):
    """
    Sub-mesh.

    :param afem.smesh.meshes.MeshGen gen: A generator.
    :param afem.smesh.meshes.Mesh mesh: A mesh.
    :param afem.smesh.meshes.MeshDS mesh_ds: A mesh data structure.
    :param OCCT.TopoDS.TopoDS_Shape sub_shape: The sub-shape.
    """

    def __init__(self, gen, mesh, mesh_ds, sub_shape):
        self._mesh = SMESH_subMesh(gen.new_id(), mesh.object, mesh_ds.object,
                                   sub_shape)
        self._ds = SubMeshDS(self)

    @property
    def object(self):
        """
        :return: The underlying sub-mesh object.
        :rtype: OCCT.SMESH.SMESH_subMesh
        """
        return self._mesh

    @property
    def id(self):
        """
        :return: The mesh ID.
        :rtype: int
        """
        return self._mesh.GetId()

    @property
    def shape(self):
        """
        :return: The sub-shape.
        :rtype: OCCT.TopoDS.TopoDS_Shape
        """
        return self._mesh.GetSubShape()

    @property
    def is_empty(self):
        """
        :return: *True* if the mesh is empty, *False* if not.
        :rtype: bool
        """
        return self._mesh.IsEmpty()

    @property
    def is_computed(self):
        """
        :return: *True* if the mesh is computed, *False* if not.
        :rtype: bool
        """
        return self._mesh.IsMeshComputed()

    @property
    def ds(self):
        """
        :return: The sub-mesh data structure.
        :rtype: afem.smesh.meshes.SubMeshDS
        """
        return self._ds

    @classmethod
    def new(cls, sub_mesh):
        """
        Create a new instance using an existing SMESH_subMesh instance.

        :param OCCT.SMESH.SMESH_subMesh sub_mesh: A SMESH_subMesh instance.

        :return: The new instance.
        :rtype: afem.smesh.meshes.SubMesh
        """
        new_mesh = cls.__new__(cls)
        new_mesh._mesh = sub_mesh
        new_mesh._ds = SubMeshDS(new_mesh)
        return new_mesh

    def can_add_hypothesis(self, hyp):
        """
        Check to see if the hypothesis can be attached.

        :param afem.smesh.hypotheses.Hypothesis hyp: The hypothesis.

        :return: *True* if can be attached, *False* if not.
        :rtype: bool
        """
        return self._mesh.CanAddHypothesis(hyp.object)

    def is_applicable_hypothesis(self, hyp):
        """
        Check to see if the hypothesis can be used to mesh.

        :param afem.smesh.hypotheses.Hypothesis hyp: The hypothesis.

        :return: *True* if applicable, *False* if not.
        :rtype: bool
        """
        return self._mesh.IsApplicableHypotesis(hyp.object)


class SubMeshDS(object):
    """
    Sub-mesh data structure.

    :param afem.smesh.meshes.SubMesh sub_mesh: A sub-mesh.
    """

    def __init__(self, sub_mesh):
        self._ds = sub_mesh.object.GetSubMeshDS()

    @property
    def object(self):
        """
        :return: The underlying mesh object.
        :rtype: OCCT.SMESH.SMESHDS_subMesh
        """
        return self._ds

    @property
    def id(self):
        """
        :return: The mesh ID.
        :rtype: int
        """
        return self._ds.GetID()

    @property
    def is_complex(self):
        """
        :return: *True* if sub-mesh is complex, *False* if not.
        :rtype: bool
        """
        return self._ds.IsComplexSubmesh()

    @property
    def is_quadratic(self):
        """
        :return: *True* if sub-mesh is quadratic, *False* if not.
        :rtype: bool
        """
        return self._ds.IsQuadratic()

    @property
    def num_nodes(self):
        """
        :return: Number of nodes.
        :rtype: int
        """
        return self._ds.NbNodes()

    @property
    def num_elms(self):
        """
        :return: Number of elements.
        :rtype: int
        """
        return self._ds.NbElements()

    @property
    def node_iter(self):
        """
        :return: Yield nodes of the sub-mesh.
        :rtype: collections.Iterable(afem.smesh.entities.Node)
        """
        iter_ = self._ds.GetNodes()
        while iter_.more():
            yield Node(iter_.next())

    @property
    def elm_iter(self):
        """
        :return: Yield elements of the sub-mesh.
        :rtype: collections.Iterable(afem.smesh.entities.Element)
        """
        iter_ = self._ds.GetElements()
        while iter_.more():
            yield Element(iter_.next())

    @classmethod
    def new(cls, sub_meshds):
        """
        Create a new instance using an existing SMESHDS_SubMesh instance.

        :param OCCT.SMESHDS.SMESHDS_SubMesh sub_meshds: A SMESHDS_SubMesh
            instance.

        :return: The new instance.
        :rtype: afem.smesh.meshes.SubMeshDS
        """
        new_ds = cls.__new__(cls)
        new_ds._ds = sub_meshds
        return new_ds

    def get_node(self, id_):
        """
        Get a node.

        :param int id_: The id of the node in the shape.

        :return: The node.
        :rtype: afem.smesh.entities.Node
        """
        return Node(self._ds.GetNode(id_))

    def get_element(self, id_):
        """
        Get an element.

        :param int id_: The id of the element in the shape.

        :return: The element.
        :rtype: afem.smesh.entities.Element
        """
        return Element(self._ds.GetElement(id_))

    def contains(self, elm):
        """
        Check to see if node or element is in the sub-mesh.

        :param elm: The node or element.
        :type elm: afem.smesh.entities.Node or afem.smesh.entities.Element

        :return: *True* if present, *False* if not.
        :rtype: bool
        """
        return self._ds.Contains(elm.object)

    def clear(self):
        """
        Clear all nodes and elements.

        :return: None.
        """
        self._ds.Clear()
