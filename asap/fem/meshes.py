from OCC.SMESH import SMESH_Gen_get

from ..topology import ShapeTools

_mesh_gen = SMESH_Gen_get()


class Mesh(object):
    """
    Mesh class.
    """
    _meshes = {}
    _indx = 0

    def __init__(self, name):
        self._name = name
        self._meshes[name] = self
        self._smesh_mesh = _mesh_gen.CreateMesh(Mesh._indx, True)
        self._id = Mesh._indx
        Mesh._indx += 1

    @property
    def id(self):
        return self._id

    @property
    def shape(self):
        return self._smesh_mesh.GetShapeToMesh()

    @property
    def smesh_obj(self):
        return self._smesh_mesh

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

    def add_hypotheses(self, shape, *hypotheses):
        """
        Add hypotheses to the shape.

        :param hypotheses:
        :param shape:

        :return:
        """
        for hypo in hypotheses:
            self._smesh_mesh.AddHypothesis(shape, hypo.id)

    def compute(self):
        """
        Compute the mesh.

        :return:
        """
        return _mesh_gen.Compute(self._smesh_mesh, self.shape)


class SubMesh(object):
    pass


class MeshData(object):
    pass
