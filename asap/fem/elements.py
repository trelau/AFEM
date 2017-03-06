from OCC.SMDS import SMDS_MeshElement

from .nodes import Node


class Element(object):
    """
    Base class for elements.
    """

    def __init__(self, smesh_element):
        assert isinstance(smesh_element, SMDS_MeshElement)
        self._elm = smesh_element

    def __str__(self):
        eid = 'Element {0}: '.format(str(self.eid))
        nids = ' '.join([str(n) for n in self.nids])
        return ''.join([eid, nids])

    def __eq__(self, other):
        return self.eid == other.eid

    def __hash__(self):
        return hash(self.eid)

    @property
    def eid(self):
        return self._elm.GetID()

    @property
    def nodes(self):
        nodes = []
        niter = self._elm.nodeIterator()
        while niter.more():
            n = Node(niter.next())
            nodes.append(n)
        return nodes

    @property
    def nids(self):
        return [n.nid for n in self.nodes]

    @property
    def ncount(self):
        return self._elm.NbNodes()

    @property
    def is_0d(self):
        return self.ncount == 1

    @property
    def is_1d(self):
        return self.ncount == 2

    @property
    def is_2d(self):
        return self.ncount > 2

    @property
    def is_tri(self):
        return self.ncount == 3

    @property
    def is_quad(self):
        return self.ncount == 4


class Elm0D(Element):
    """
    Generic 0-D element.
    """

    def __init__(self, smesh_element):
        super(Elm0D, self).__init__(smesh_element)


class Elm1D(Element):
    """
    Generic 1-D element.
    """

    def __init__(self, smesh_element):
        super(Elm1D, self).__init__(smesh_element)


class Elm2D(Element):
    """
    Generic 2-D element.
    """

    def __init__(self, smesh_element):
        super(Elm2D, self).__init__(smesh_element)
