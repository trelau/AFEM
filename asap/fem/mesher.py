from OCC.TopAbs import TopAbs_EDGE, TopAbs_FACE

from .data import MeshData
from .methods.mesh_1d import mesh_edge
from .methods.mesh_2d import mesh_face
from ..topology import ShapeTools
from ..graphics.viewer import ViewableItem


class MeshShape(object):
    """
    Shape mesher.
    """

    @staticmethod
    def perform(shape, maxh=None, density=None, quad_dominated=True):
        """
        Mesh the shape.

        :param shape:
        :param maxh:
        :param density:
        :param quad_dominated:

        :return:
        """
        if not ShapeTools.is_shape(shape):
            return None

        if shape.ShapeType() == TopAbs_EDGE:
            return EdgeMesh(shape, maxh, density)

        if shape.ShapeType() == TopAbs_FACE:
            return FaceMesh(shape, maxh, density, quad_dominated)

        return None


class ShapeMesh(object):
    """
    Base class for shape meshes.
    """

    def __init__(self, shape):
        self._shape = shape
        self._nodes = []
        self._elms = []

    @property
    def shape(self):
        return self._shape

    @property
    def elements(self):
        return [e for e in self._elms]


class EdgeMesh(ShapeMesh):
    """
    Edge mesh.
    """

    def __init__(self, edge, maxh=None, density=None):
        super(EdgeMesh, self).__init__(edge)
        if ShapeTools.is_shape(edge) and edge.ShapeType() == TopAbs_EDGE:
            self._perform(edge, maxh, density)

    def _perform(self, edge, maxh, density):
        """
        Perform the edge mesh.
        """
        nodes = mesh_edge(edge, maxh, density)
        self._nodes = nodes
        MeshData.shape_to_mesh(edge, self)

    @property
    def nodes(self):
        return [n for n in self._nodes]


class FaceMesh(ShapeMesh, ViewableItem):
    """
    Face mesh.
    """

    def __init__(self, face, maxh=None, density=None, quad_dominated=True):
        super(FaceMesh, self).__init__(face)
        ViewableItem.__init__(self)
        if ShapeTools.is_shape(face) and face.ShapeType() == TopAbs_FACE:
            self._perform(face, maxh, density, quad_dominated)

    def _perform(self, face, maxh, density, quad_dominated):
        """
        Perform the face mesh.
        """
        # Mesh all the edges in the face.
        for e in ShapeTools.get_edges(face):
            if not MeshData.has_mesh(e):
                EdgeMesh(e, maxh, density)

        # Mesh the face.
        elms = mesh_face(face, quad_dominated)
        self._elms = elms
        MeshData.shape_to_mesh(face, self)

    @property
    def nodes(self):
        node_set = set()
        for e in self._elms:
            for n in e.nodes:
                node_set.add(n)
        return list(node_set)
