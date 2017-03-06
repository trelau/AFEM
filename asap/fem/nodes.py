from OCC.SMDS import SMDS_MeshNode
from numpy import array, float64


class Node(object):
    """
    Node.
    """

    def __init__(self, smesh_node):
        assert isinstance(smesh_node, SMDS_MeshNode)
        self._node = smesh_node

    def __str__(self):
        return 'Node {0}: ({1}, {2}, {3})'.format(self.nid, *self.xyz)

    def __eq__(self, other):
        return self.nid == other.nid

    def __hash__(self):
        return hash(self.nid)

    @property
    def nid(self):
        return self._node.GetID()

    @property
    def x(self):
        return self._node.X()

    @property
    def y(self):
        return self._node.Y()

    @property
    def z(self):
        return self._node.Z()

    @property
    def xyz(self):
        return array([self.x, self.y, self.z], dtype=float64)
