from OCCT.SMDS import SMDS_MeshElement

from .nodes import Node

__all__ = ["Element", "Elm0D", "Elm1D", "Elm2D"]


class Element(object):
    """
    Base class for elements.

    :param OCCT.SMDS.SMDS_MeshElement the_element: The element.
    """

    def __init__(self, the_element):
        assert isinstance(the_element, SMDS_MeshElement)
        self._elm = the_element

    def __str__(self):
        eid = 'Element {0}: '.format(str(self.eid))
        nids = ' '.join([str(n) for n in self.nids])
        return ''.join([eid, nids])

    def __eq__(self, other):
        return self.eid == other.eid

    def __hash__(self):
        return hash(self.eid)

    @property
    def object(self):
        """
        :return: The underlying element.
        :rtype: OCCT.SMDS.SMDS_Element
        """
        return self._elm

    @property
    def eid(self):
        """
        :return: The element ID.
        :rtype: int
        """
        return self._elm.GetID()

    @property
    def nodes(self):
        """
        :return: The nodes of the element.
        :rtype: list[afem.fem.nodes.Node]
        """
        nodes = []
        niter = self._elm.nodeIterator()
        while niter.more():
            n = Node(niter.next())
            nodes.append(n)
        return nodes

    @property
    def nids(self):
        """
        :return: The node ID's of the element.
        :rtype: list[int]
        """
        return [n.nid for n in self.nodes]

    @property
    def ncount(self):
        """
        :return: The number of nodes in the element.
        :rtype: int
        """
        return self._elm.NbNodes()

    @property
    def is_0d(self):
        """
        :return: ``True`` if the element has one node, ``False`` if not.
        :rtype: bool
        """
        return self.ncount == 1

    @property
    def is_1d(self):
        """
        :return: ``True`` if the element has two nodes, ``False`` if not.
        :rtype: bool
        """
        return self.ncount == 2

    @property
    def is_2d(self):
        """
        :return: ``True`` if the element has more than two nodes, ``False`` if
            not.
        :rtype: bool
        """
        return self.ncount > 2

    @property
    def is_tri(self):
        """
        :return: ``True`` if the element has three nodes, ``False`` if not.
        :rtype: bool
        """
        return self.ncount == 3

    @property
    def is_quad(self):
        """
        :return: ``True`` if the element has four nodes, ``False`` if not.
        :rtype: bool
        """
        return self.ncount == 4


class Elm0D(Element):
    """
    Generic 0-D element.

    :param OCCT.SMDS.SMDS_MeshElement the_element: The SMESH element.
    """

    def __init__(self, the_element):
        super(Elm0D, self).__init__(the_element)


class Elm1D(Element):
    """
    Generic 1-D element.

    :param OCCT.SMDS.SMDS_MeshElement the_element: The SMESH element.
    """

    def __init__(self, the_element):
        super(Elm1D, self).__init__(the_element)


class Elm2D(Element):
    """
    Generic 2-D element.

    :param OCCT.SMDS.SMDS_MeshElement the_element: The SMESH element.
    """

    def __init__(self, the_element):
        super(Elm2D, self).__init__(the_element)
