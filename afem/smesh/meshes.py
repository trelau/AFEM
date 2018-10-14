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
from OCCT.SMDSAbs import SMDSAbs_ElementType
from OCCT.SMESH import SMESH_Gen, SMESH_subMesh

from afem.smesh.entities import Node, Element
from afem.topology.entities import Shape

__all__ = ["MeshGen",
           "Mesh", "MeshDS",
           "SubMesh", "SubMeshDS",
           "MeshGroup"]


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
    def wrap(cls, gen):
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

        :param afem.topology.entities.Shape shape: The shape to mesh.
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
        :param afem.topology.entities.Shape shape: A shape.

        :return: *True* if ok, *False* if not.
        :rtype: bool
        """
        return self._gen.CheckAlgoState(mesh.object, shape.object)

    def compute(self, mesh, shape=None):
        """
        Compute a mesh on a shape.

        :param afem.smesh.meshes.Mesh mesh: A mesh.
        :param afem.topology.entities.Shape shape: The shape to compute mesh
            on. If not provided then the shape associated to the mesh is used.

        :return: *True* if computed, *False* if not.
        :rtype: bool

        :raise ValueError: If no shape is available to apply the hypothesis to.
        """
        if shape is None:
            if mesh.has_shape:
                shape = mesh.shape
            else:
                raise ValueError('No shape could be found.')
        return self._gen.Compute(mesh.object, shape.object)


class Mesh(object):
    """
    Mesh.

    :param afem.smesh.meshes.MeshGen gen: The MeshGen instance.
    :param bool is_embedded: Option for embedding mesh.

    :cvar OCCT.SMDSAbs.SMDSAbs_ElementType ALL: SMESH all elements type.
    :cvar OCCT.SMDSAbs.SMDSAbs_ElementType NODE: SMESH node type.
    :cvar OCCT.SMDSAbs.SMDSAbs_ElementType EDGE: SMESH edge type.
    :cvar OCCT.SMDSAbs.SMDSAbs_ElementType VOLUME: SMESH volume type.
    """
    ALL = SMDSAbs_ElementType.SMDSAbs_All
    NODE = SMDSAbs_ElementType.SMDSAbs_Node
    EDGE = SMDSAbs_ElementType.SMDSAbs_Edge
    FACE = SMDSAbs_ElementType.SMDSAbs_Face
    VOLUME = SMDSAbs_ElementType.SMDSAbs_Volume

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
        :rtype: afem.topology.entities.Shape
        """
        return Shape.wrap(self._mesh.GetShapeToMesh())

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
    def num_volumes(self):
        """
        :return: Number of volume elements.
        :rtype: int
        """
        return self._mesh.NbVolumes()

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

    def shape_to_mesh(self, shape):
        """
        Set the shape to mesh.

        :param afem.topology.entities.Shape shape: The shape.

        :return: None
        """
        self._mesh.ShapeToMesh(shape.object)

    def add_hypothesis(self, hypothesis, shape=None):
        """
        Add a hypothesis to the shape.

        :param afem.smesh.hypotheses.Hypothesis hypothesis: The hypothesis.
        :param shape: The shape the hypothesis applies to. This can be a
            sub-shape of the master shape. If not provided then the master
            shape is used.
        :type shape: afem.topology.entities.Shape

        :return: Status of adding hypothesis.
        :rtype: OCCT.SMESH.SMESH_Hypothesis.Hypothesis_Status

        :raise ValueError: If no shape is available to apply the hypothesis to.
        """
        if shape is None:
            if self.has_shape:
                shape = self.shape
            else:
                raise ValueError('No shape could be found.')

        return self._mesh.AddHypothesis(shape.object, hypothesis.id)

    def add_hypotheses(self, hypotheses, shape=None):
        """
        Add a hypotheses to the shape.

        :param hypotheses: The hypotheses to add.
        :type hypotheses:
            collections.Sequence(afem.smesh.hypotheses.Hypothesis)
        :param shape: The shape the hypothesis applies to. This can be a
            sub-shape of the master shape. If not provided then the master
            shape is used.
        :type shape: afem.topology.entities.Shape

        :return: Dictionary where the key is the hypothesis and the value is
            the hypothesis status.
        :rtype: dict

        :raise ValueError: If no shape is available to apply the hypothesis to.
        """
        status_dict = {}
        for hyp in hypotheses:
            status = self.add_hypothesis(hyp, shape)
            status_dict[hyp] = status
        return status_dict

    def clear(self):
        """
        Clear all nodes and elements.

        :return: None.
        """
        self._mesh.Clear()

    def clear_submesh(self, shape):
        """
        Clear the nodes and elements from the shape.

        :param afem.topology.entities.Shape shape: The shape.

        :return: None.
        """
        shape_id = self.ds.shape_to_index(shape)
        self._mesh.ClearSubMesh(shape_id)

    def get_submesh(self, sub_shape):
        """
        Get a sub-mesh for the sub-shape.

        :param afem.topology.entities.Shape sub_shape: A sub-shape.

        :return: A sub-mesh.
        :rtype: afem.smesh.meshes.SubMesh
        """
        sub_mesh = self._mesh.GetSubMesh(sub_shape.object)
        return SubMesh.wrap(sub_mesh)

    def get_submesh_containing(self, sub_shape):
        """
        Get a sub-mesh for containing the sub-shape.

        :param afem.topology.entities.Shape sub_shape: A sub-shape.

        :return: A sub-mesh.
        :rtype: afem.smesh.meshes.SubMesh
        """
        sub_mesh = self._mesh.GetSubMeshContaining(sub_shape.object)
        return SubMesh.wrap(sub_mesh)

    def create_group(self, name, type_, shape=None):
        """
        Create a group belonging to the mesh.

        :param str name: The name of the group.
        :param OCCT.SMDSAbs.SMDSAbs_ElementType type_: The element type for the
            group.
        :param afem.topology.entities.Shape shape: The shape to create the
            group on. If *None* provided, then the group is not on geometry.

        :return: The group.
        :rtype: afem.smesh.meshes.MeshGroup
        """
        return MeshGroup(self, name, type_, shape)

    def export_dat(self, fn):
        """
        Export the mesh to a DAT file.

        :param str fn: The output file.

        :return: None
        """
        self._mesh.ExportDAT(fn)

    def export_stl(self, fn, is_ascii=True):
        """
        Export the mesh to an STL file.

        :param str fn: The output file.
        :param bool is_ascii: ASCII text output or binary.

        :return: None.
        """
        self._mesh.ExportSTL(fn, is_ascii)

    def export_unv(self, fn):
        """
        Export the mesh to a UNV file.

        :param str fn: The output file.

        :return: None
        """
        self._mesh.ExportUNV(fn)

    def import_unv(self, fn):
        """
        Import a mesh from a UNV file.

        :param str fn: The input file.

        :return: None.
        """
        return self._mesh.UNVToMesh(fn)

    @classmethod
    def wrap(cls, mesh):
        """
        Create a new instance using an existing SMESH_Mesh instance.

        :param OCCT.SMESH.SMESH_Mesh mesh: A SMESH_Mesh instance.

        :return: The new instance.
        :rtype: afem.smesh.meshes.Mesh
        """
        new_mesh = cls.__new__(cls)
        new_mesh._mesh = mesh
        new_mesh._ds = MeshDS(new_mesh)
        return new_mesh


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

        :param afem.topology.entities.Shape shape: The shape.

        :return: *True* if the shape has elements, *False* if not.
        :rtype: bool
        """
        return self._ds.HasMeshElements(shape.object)

    def mesh_elements(self, shape):
        """
        Get elements of shape.

        :param afem.topology.entities.Shape shape: The shape.

        :return: Sub-mesh data structure of elements.
        :rtype: afem.smesh.meshes.SubMeshDS
        """
        sub_meshds = self._ds.MeshElements(shape.object)
        return SubMeshDS.wrap(sub_meshds)

    def shape_to_index(self, shape):
        """
        Get the shape index.

        :param afem.topology.entities.Shape shape: The shape.

        :return: The shape index.
        :rtype: int
        """
        return self._ds.ShapeToIndex(shape.object)

    def index_to_shape(self, indx):
        """
        Get the shape from an index.

        :param int indx: The index.

        :return: The shape.
        :rtype: afem.topology.entities.Shape
        """
        return Shape.wrap(self._ds.IndexToShape(indx))

    @classmethod
    def wrap(cls, mesh):
        """
        Create a new instance using an existing SMESHDS_Mesh instance.

        :param OCCT.SMESHDS.SMESHDS_Mesh mesh: A SMESHDS_Mesh instance.

        :return: The new instance.
        :rtype: afem.smesh.meshes.MeshDS
        """
        new_mesh = cls.__new__(cls)
        new_mesh._ds = mesh
        return new_mesh


class SubMesh(object):
    """
    Sub-mesh.

    :param afem.smesh.meshes.MeshGen gen: A generator.
    :param afem.smesh.meshes.Mesh mesh: A mesh.
    :param afem.smesh.meshes.MeshDS mesh_ds: A mesh data structure.
    :param afem.topology.entities.Shape sub_shape: The sub-shape.
    """

    def __init__(self, gen, mesh, mesh_ds, sub_shape):
        self._mesh = SMESH_subMesh(gen.new_id(), mesh.object, mesh_ds.object,
                                   sub_shape.object)
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
        :rtype: afem.topology.entities.Shape
        """
        return Shape.wrap(self._mesh.GetSubShape())

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
    def wrap(cls, sub_mesh):
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
    def wrap(cls, sub_meshds):
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


class MeshGroup(object):
    """
    Mesh group for nodes and elements.

    :param afem.smesh.meshes.Mesh mesh: The mesh the group belongs to.
    :param str name: The name of the group.
    :param OCCT.SMDSAbs.SMDSAbs_ElementType type_: The element type for the
        group.
    :param afem.topology.entities.Shape shape: The shape to create the
            group on. If *None* provided, then the group is not on geometry.
    """

    def __init__(self, mesh, name, type_, shape=None):
        if isinstance(shape, Shape):
            group = mesh.object.AddGroup(type_, name, -1, shape.object)
        else:
            group = mesh.object.AddGroup(type_, name, -1)
        self._group = group
        self._ds = group.GetGroupDS()
        self._mesh = mesh

    @property
    def object(self):
        """
        :return: The underlying SMESH_Group object.
        :rtype: OCCT.SMESH.SMESH_Group
        """
        return self._group

    @property
    def id(self):
        """
        :return: The group ID.
        :rtype: int
        """
        return self._ds.GetID()

    @property
    def type(self):
        """
        :return: The group type.
        :rtype: OCCT.SMDSAbs.SMDSAbs_ElementType
        """
        return self._ds.GetType()

    @property
    def mesh(self):
        """
        :return: The mesh the group belongs to.
        :rtype: afem.smesh.meshes.Mesh
        """
        return self._mesh

    @property
    def name(self):
        """
        :Getter: The group name.
        :Setter: Set the group name.
        :type: str
        """
        return self._group.GetName()

    @name.setter
    def name(self, name):
        self.set_name(name)

    @property
    def shape(self):
        """
        :Getter: The shape of the group.
        :Setter: Set the group shape.
        :type: afem.topology.entities.Shape

        :raise AttributeError: If the group is not associated with a shape.
        """
        return Shape.wrap(self._ds.GetShape())

    @shape.setter
    def shape(self, shape):
        self.set_shape(shape)

    @property
    def is_empty(self):
        """
        :return: Check if group is empty.
        :rtype: bool
        """
        return self._ds.IsEmpty()

    @property
    def size(self):
        """
        :return: The number of entities in the group.
        :rtype: int
        """
        return self._ds.Extent()

    @property
    def node_iter(self):
        """
        :return: Yield the nodes in the group.
        :rtype: collections.Iterable(afem.smesh.entities.Node)

        :raise TypeError: If group is not a node group.
        """
        if self.type != Mesh.NODE:
            raise TypeError('Group is not a node group.')

        it = self._ds.GetElements()
        while it.more():
            yield Node(it.next())

    @property
    def edge_iter(self):
        """
        :return: Yield the edges in the group.
        :rtype: collections.Iterable(afem.smesh.entities.Element)

        :raise TypeError: If group is not an edge group.
        """
        if self.type != Mesh.EDGE:
            raise TypeError('Group is not an edge group.')

        it = self._ds.GetElements()
        while it.more():
            yield Element(it.next())

    @property
    def face_iter(self):
        """
        :return: Yield the faces in the group.
        :rtype: collections.Iterable(afem.smesh.entities.Element)

        :raise TypeError: If group is not a face group.
        """
        if self.type != Mesh.FACE:
            raise TypeError('Group is not a face group.')

        it = self._ds.GetElements()
        while it.more():
            yield Element(it.next())

    def set_name(self, name):
        """
        Set the group name.

        :param str name: The group name.

        :return: None.
        """
        self._group.SetName(name)

    def set_shape(self, shape):
        """
        Set the shape of the group.

        :param afem.topology.entities.Shape shape: The shape.

        :return: None.
        """
        self._ds.SetShape(shape.object)

    def contains_id(self, eid):
        """
        Check if group contains the element by using its ID.

        :param int eid: The element ID.

        :return: *True* if in the group, *False* otherwise.
        :rtype: bool
        """
        return self._ds.Contains(eid)

    def contains_node(self, node):
        """
        Check if group contains the node.

        :param afem.smesh.entities.Node node: The node.

        :return: *True* if in the group, *False* otherwise.
        :rtype: bool

        :raise TypeError: If the group is not a node group.
        """
        if self.type != Mesh.NODE:
            raise TypeError('Group is not a node group.')

        return self._ds.Contains(node.object)

    def contains_elm(self, elm):
        """
        Check if group contains the element.

        :param afem.smesh.entities.Element elm: The element.

        :return: *True* if in the group, *False* otherwise.
        :rtype: bool

        :raise TypeError: If the group is not an element group.
        """
        if self.type == Mesh.NODE:
            raise TypeError('Group is not an element group.')

        return self._ds.Contains(elm.object)

    def union(self, other, name='union group'):
        """
        Union the entities of this group with another.

        :param afem.smesh.meshes.MeshGroup other: The other group.
        :param str name: The name of the new group.

        :return: New group.
        :rtype: afem.smesh.meshes.MeshGroup

        :raise TypeError: If the two groups are of different type or they do
            not share the same mesh.
        """
        if self.type != other.type:
            raise TypeError('Groups are not the same type.')
        if self.mesh.id != other.mesh.id:
            raise TypeError('Groups do not share the same mesh.')

        # Add elements of each group
        new_group = MeshGroup(self.mesh, name, self.type)
        it = self._ds.GetElements()
        while it.more():
            new_group._ds.Add(it.next())
        it = other._ds.GetElements()
        while it.more():
            new_group._ds.Add(it.next())

        return new_group

    def intersect(self, other, name='intersect group'):
        """
        Intersect the entities between this group and another.

        :param afem.smesh.meshes.MeshGroup other: The other group.
        :param str name: The name of the new group.

        :return: New group.
        :rtype: afem.smesh.meshes.MeshGroup

        :raise TypeError: If the two groups are of different type or they do
            not share the same mesh.
        """
        if self.type != other.type:
            raise TypeError('Groups are not the same type.')
        if self.mesh.id != other.mesh.id:
            raise TypeError('Groups do not share the same mesh.')

        # Add elements of first group if they are in other
        new_group = MeshGroup(self.mesh, name, self.type)
        it = self._ds.GetElements()
        while it.more():
            e = it.next()
            if other._ds.Contains(e):
                new_group._ds.Add(e)

        return new_group

    def subtract(self, other, name='Subtract group'):
        """
        Subtract the entities from this group by another.

        :param afem.smesh.meshes.MeshGroup other: The other group.
        :param str name: The name of the new group.

        :return: New group.
        :rtype: afem.smesh.meshes.MeshGroup

        :raise TypeError: If the two groups are of different type or they do
            not share the same mesh.
        """
        if self.type != other.type:
            raise TypeError('Groups are not the same type.')
        if self.mesh.id != other.mesh.id:
            raise TypeError('Groups do not share the same mesh.')

        # Add elements of first group if they are not in other
        new_group = MeshGroup(self.mesh, name, self.type)
        it = self._ds.GetElements()
        while it.more():
            e = it.next()
            if not other._ds.Contains(e):
                new_group._ds.Add(e)

        return new_group
