from OCC.gp import gp_Dir, gp_Vec
from numpy import add, array, float64, subtract
from numpy.linalg import norm

from .base import Geometry

__all__ = ["Direction", "Vector"]


class Direction(gp_Dir, Geometry):
    """
    3-D unit vector.
    """

    def __init__(self, *args):
        super(Direction, self).__init__(*args)
        Geometry.__init__(self)

    def __array__(self, dtype=float64, copy=True, order=None, subok=False,
                  ndmin=0):
        return array(self.ijk, dtype=dtype, copy=copy, order=order,
                     subok=subok, ndmin=ndmin)

    def __iter__(self):
        for elm in self.ijk:
            yield elm

    def __len__(self):
        return 3

    def __getitem__(self, item):
        return self.ijk[item]

    def __add__(self, other):
        return add(self, other)

    def __sub__(self, other):
        return subtract(self, other)

    @property
    def i(self):
        return self.X()

    @property
    def j(self):
        return self.Y()

    @property
    def k(self):
        return self.Z()

    @property
    def ijk(self):
        return array([self.i, self.j, self.k], dtype=float64)

    @property
    def xyz(self):
        return self.ijk

    @property
    def mag(self):
        return 1.


class Vector(gp_Vec, Geometry):
    """
    3-D vector.
    """

    def __init__(self, *args):
        super(Vector, self).__init__(*args)
        Geometry.__init__(self)

    def __array__(self, dtype=float64, copy=True, order=None, subok=False,
                  ndmin=0):
        return array(self.xyz, dtype=dtype, copy=copy, order=order,
                     subok=subok, ndmin=ndmin)

    def __iter__(self):
        for elm in self.xyz:
            yield elm

    def __len__(self):
        return 3

    def __getitem__(self, item):
        return self.xyz[item]

    def __add__(self, other):
        return add(self, other)

    def __sub__(self, other):
        return subtract(self, other)

    @property
    def x(self):
        return self.X()

    @property
    def y(self):
        return self.Y()

    @property
    def z(self):
        return self.Z()

    @property
    def xyz(self):
        return array([self.x, self.y, self.z], dtype=float64)

    @property
    def mag(self):
        return norm(self.xyz)

    @property
    def ijk(self):
        return self.xyz / self.mag

    def reverse(self):
        """
        Reverse the direction of the vector.

        :return: None.
        """
        self.Reverse()

    def normalize(self):
        """
        Normalize the vector.

        :return:
        """
        self.Normalize()

    def scale(self, scale):
        """
        Scale the vector.

        :param scale:

        :return:
        """
        self.Scale(scale)
