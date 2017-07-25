from OCC.SMESH import SMESH_Gen_get

from .hypotheses import HypothesisData
from .nodes import Node
from ..topology import ShapeTools

_mesh_gen = SMESH_Gen_get()

__all__ = ["Mesh", "SubMesh", "MeshData"]


class Mesh(object):
    """
    Mesh.
    """
    _all = {}
    _indx = 0

    def __init__(self, name):
        self._name = name
        Mesh._all[name] = self
        self._mesh = _mesh_gen.CreateMesh(Mesh._indx, True)
        self._ds = self._mesh.GetMeshDS()
        self._id = Mesh._indx
        Mesh._indx += 1

    @property
    def id(self):
        return self._id

    @property
    def shape(self):
        return self._mesh.GetShapeToMesh()

    @property
    def has_shape(self):
        return self._mesh.HasShapeToMesh()

    @property
    def smesh_obj(self):
        return self._mesh

    @property
    def nb_nodes(self):
        return self._ds.NbNodes()

    @property
    def min_node_id(self):
        return self._ds.MinNodeID()

    @property
    def max_node_id(self):
        return self._ds.MaxNodeID()

    @property
    def min_elm_id(self):
        return self._ds.MinElementID()

    @property
    def max_elm_id(self):
        return self._ds.MaxElementID()

    @property
    def nodes(self):
        return self.get_nodes()

    @classmethod
    def get_mesh(cls, mesh=None):
        """
        Get a mesh.

        :param mesh:

        :return:
        """
        if isinstance(mesh, Mesh):
            return mesh
        try:
            return Mesh._all[mesh]
        except KeyError:
            return None

    def activate(self):
        """
        Activate this mesh.

        :return:
        """
        MeshData._active = self
        return True

    def shape_to_mesh(self, shape):
        """
        Set the shape to mesh.

        :param shape:

        :return:
        """
        shape = ShapeTools.to_shape(shape)
        if not shape:
            return False
        self._mesh.ShapeToMesh(shape)
        return True

    def add_hypothesis(self, hypothesis, shape=None):
        """
        Add a hypothesis to the shape.

        :param hypothesis:
        :param shape:

        :return:
        """
        shape = ShapeTools.to_shape(shape)
        if not shape:
            if self.has_shape:
                shape = self.shape
            else:
                return False
        hypothesis = HypothesisData.get_hypothesis(hypothesis)
        if not hypothesis:
            return False
        self._mesh.AddHypothesis(shape, hypothesis.id)
        return True

    def compute(self):
        """
        Compute the mesh.

        :return:
        """
        return _mesh_gen.Compute(self._mesh, self.shape)

    def clear(self):
        """
        Clear all nodes and elements.

        :return:
        """
        self._mesh.Clear()
        return True

    def get_submesh(self, sub_shape):
        """
        Get a SubMesh from a sub-shape.

        :param sub_shape:

        :return:
        """
        mesh = self._mesh.GetSubMesh(sub_shape)
        return SubMesh(mesh)

    def get_nodes(self, order=False):
        """
        Get nodes from the mesh.

        :param bool order: Order nodes order by ID.

        :return:
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
    """

    def __init__(self, smesh_submesh):
        self._mesh = smesh_submesh
        self._ds = smesh_submesh.GetSubMeshDS()

    @property
    def is_empty(self):
        return self._mesh.IsEmpty()

    @property
    def is_computed(self):
        return self._mesh.IsMeshComputed()

    @property
    def nb_nodes(self):
        return self._ds.NbNodes()

    @property
    def nodes(self):
        return self.get_nodes()

    def get_nodes(self, include_subshapes=True):
        """
        Get nodes from the sub-mesh.

        :param bool include_subshapes: Option to include sub-shapes when
            retrieving nodes.

        :return:
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


class MeshData(object):
    """
    Mesh data manager.
    """
    _active = None
    hypotheses = HypothesisData()

    @classmethod
    def get_active(cls):
        """
        Get the active mesh.

        :return:
        """
        return cls._active

    @classmethod
    def get_mesh(cls, mesh=None):
        """
        Get mesh.

        :param mesh:

        :return:
        """
        mesh = Mesh.get_mesh(mesh)
        if mesh:
            return mesh
        return cls._active

    @classmethod
    def make_active(cls, mesh):
        """
        Activate the mesh.

        :param mesh:

        :return:
        """
        mesh = cls.get_mesh(mesh)
        if not mesh:
            return False
        mesh.activate()
        return True

    @classmethod
    def create_mesh(cls, name, shape, active=True):
        """
        Create a mesh

        :param name:
        :param shape:
        :param active:

        :return:
        """
        mesh = Mesh(name)
        mesh.shape_to_mesh(shape)
        if active:
            mesh.activate()
        return mesh

    @classmethod
    def add_hypothesis(cls, hypothesis, shape=None, mesh=None):
        """
        Add a hypothesis to the shape in a mesh.

        :param hypothesis:
        :param shape:
        :param mesh:

        :return:
        """
        mesh = cls.get_mesh(mesh)
        if not mesh:
            return False
        return mesh.add_hypothesis(hypothesis, shape)

    @classmethod
    def compute_mesh(cls, mesh=None):
        """
        Compute a mesh.

        :param mesh:

        :return:
        """
        mesh = cls.get_mesh(mesh)
        if not mesh:
            return False
        return mesh.compute()

    @classmethod
    def get_submesh(cls, sub_shape, mesh=None):
        """
        Get a SubMesh from the sub-shape.

        :param sub_shape:
        :param mesh:

        :return:
        """
        mesh = cls.get_mesh(mesh)
        if not mesh:
            return None
        return mesh.get_submesh(sub_shape)

    @classmethod
    def get_nodes(cls, order=False, mesh=None):
        """
        Get nodes from the mesh.

        :param bool order: Order nodes order by ID.
        :param mesh:

        :return:
        """
        mesh = cls.get_mesh(mesh)
        if not mesh:
            return Node
        return mesh.get_nodes(order)
