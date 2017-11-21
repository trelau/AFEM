from OCCT.SMDS import SMDS_MeshNode
from numpy import array

__all__ = ["Node"]


class Node(object):
    """
    Generic mesh node.

    :param OCCT.SMDS.SMDS_MeshNode the_node: The node.
    """

    def __init__(self, the_node):
        assert isinstance(the_node, SMDS_MeshNode)
        self._node = the_node

    def __str__(self):
        return 'Node {0}: ({1}, {2}, {3})'.format(self.nid, *self.xyz)

    def __eq__(self, other):
        return self.nid == other.nid

    def __hash__(self):
        return hash(self.nid)

    @property
    def handle(self):
        """
        :return: The underlying node.
        :rtype: OCCT.SMDS.SMDS_MeshNode
        """
        return self._node

    @property
    def nid(self):
        """
        :return: The node ID.
        :rtype: int
        """
        return self._node.GetID()

    @property
    def x(self):
        """
        :return: The node x-location.
        :rtype: float
        """
        return self._node.X()

    @property
    def y(self):
        """
        :return: The node y-location.
        :rtype: float
        """
        return self._node.Y()

    @property
    def z(self):
        """
        :return: The node z-location.
        :rtype: float
        """
        return self._node.Z()

    @property
    def xyz(self):
        """
        :return: The node xyz-location.
        :rtype: numpy.ndarray
        """
        return array([self.x, self.y, self.z], dtype=float)
