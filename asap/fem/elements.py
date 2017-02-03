from itertools import count


class Element(object):
    """
    Base class for elements.
    """
    _ids = count(1)

    def __init__(self, nodes):
        self._nodes = nodes
        self._eid = self._ids.next()

    def __str__(self):
        eid = 'Element {0}: '.format(str(self.eid))
        nids = ' '.join([str(n) for n in self.nids])
        return ''.join([eid, nids])

    @property
    def eid(self):
        return self._eid

    @property
    def nodes(self):
        return [n for n in self._nodes]

    @property
    def nids(self):
        return [n.nid for n in self._nodes]

    @property
    def ncount(self):
        return len(self._nodes)

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

    def __init__(self, nodes):
        super(Elm0D, self).__init__(nodes)


class Elm1D(Element):
    """
    Generic 1-D element.
    """

    def __init__(self, nodes):
        super(Elm1D, self).__init__(nodes)


class Elm2D(Element):
    """
    Generic 2-D element.
    """

    def __init__(self, nodes):
        super(Elm2D, self).__init__(nodes)
