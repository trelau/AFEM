from numpy import array, float64
from numpy.linalg import norm

from ..config import Settings


class Node(object):
    """
    Node.
    """
    _ids = 1

    def __init__(self, t=None, u=None, v=None, pnt=None):
        self.t = t
        self.u = u
        self.v = v
        self.pnt = pnt
        self.h = None
        self.is_boundary = False
        self._nid = Node._ids
        Node._ids += 1

    def __str__(self):
        return 'Node {0}: ({1}, {2}, {3})'.format(self._nid, *self.pnt)

    @property
    def nid(self):
        return self._nid

    @property
    def uv(self):
        return array([self.u, self.v], dtype=float64)

    @property
    def x(self):
        return self.pnt[0]

    @property
    def y(self):
        return self.pnt[1]

    @property
    def z(self):
        return self.pnt[2]

    def distance(self, node):
        """
        Calculate the 3-D distance to another node.

        :param node: Other node.
        :type node: :class:`.Node`

        :return: Distance to other node.
        :rtype: float
        """
        return norm(self.pnt - node.pnt)

    def is_equal(self, node, tol=None):
        """
        Check to see if the nodes are coincident.

        :param node: Other node.
        :type node: :class:`.Node`
        :param float tol: Tolerance.

        :return: *True* if two nodes are coincident, *False* if not.
        :rtype: bool
        """
        if tol is None:
            tol = Settings.mtol
        return self.distance(node) <= tol
