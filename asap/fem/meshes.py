from OCC.SMESH import SMESH_Gen_get

from .hypotheses import HypothesisData
from ..topology import ShapeTools

_mesh_gen = SMESH_Gen_get()


class Mesh(object):
    """
    Mesh class.
    """
    _all = {}
    _indx = 0

    def __init__(self, name):
        self._name = name
        Mesh._all[name] = self
        self._smesh_mesh = _mesh_gen.CreateMesh(Mesh._indx, True)
        self._smesh_ds = self._smesh_mesh.GetMeshDS()
        self._id = Mesh._indx
        Mesh._indx += 1

    @property
    def id(self):
        return self._id

    @property
    def shape(self):
        return self._smesh_mesh.GetShapeToMesh()

    @property
    def has_shape(self):
        return self._smesh_mesh.HasShapeToMesh()

    @property
    def smesh_obj(self):
        return self._smesh_mesh

    @property
    def nb_nodes(self):
        return self._smesh_ds.NbNodes()

    @property
    def min_node_id(self):
        return self._smesh_ds.MinNodeID()

    @property
    def max_node_id(self):
        return self._smesh_ds.MaxNodeID()

    @property
    def min_elm_id(self):
        return self._smesh_ds.MinElementID()

    @property
    def max_elm_id(self):
        return self._smesh_ds.MaxElementID()

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
        self._smesh_mesh.ShapeToMesh(shape)
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
        self._smesh_mesh.AddHypothesis(shape, hypothesis.id)
        return True

    def compute(self):
        """
        Compute the mesh.

        :return:
        """
        return _mesh_gen.Compute(self._smesh_mesh, self.shape)

    def clear(self):
        """
        Clear all nodes and elements.

        :return:
        """
        self._smesh_mesh.Clear()
        return True

    def get_submesh(self, subshape):
        """
        Get a SubMesh from a subshape.

        :param subshape:

        :return:
        """
        # TODO Get SubMeshes.


class SubMesh(object):
    pass


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
