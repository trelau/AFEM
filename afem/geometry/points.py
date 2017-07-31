from OCC.gp import gp_Pnt, gp_Pnt2d, gp_XYZ
from numpy import add, array, float64, subtract

from .base import Geometry
from ..config import Settings
from ..utils.check import is_array_like

__all__ = ["Point", "Point2D"]


class Point(gp_Pnt, Geometry):
    """
    Point.
    """

    def __init__(self, *args):
        super(Point, self).__init__(*args)
        Geometry.__init__(self)

    def __str__(self):
        return 'Point = ({0}, {1}, {2})'.format(*self.xyz)

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
    def xyz(self):
        return array([self.X(), self.Y(), self.Z()], dtype=float64)

    @property
    def x(self):
        return self.X()

    @x.setter
    def x(self, x):
        self.SetX(x)

    @property
    def y(self):
        return self.Y()

    @y.setter
    def y(self, y):
        self.SetY(y)

    @property
    def z(self):
        return self.Z()

    @z.setter
    def z(self, z):
        self.SetZ(z)

    def copy(self):
        """
        Return a new copy of the point.

        :return:
        """
        return Point(*self.xyz)

    def set_xyz(self, xyz):
        """
        Set point coordinates.

        :param array_like xyz: Point coordinates.

        :return: *True* if set, *False* if not.
        :rtype: bool
        """
        if isinstance(xyz, gp_Pnt):
            self.SetXYZ(xyz.XYZ())
            return True
        if isinstance(xyz, gp_XYZ):
            self.SetXYZ(xyz)
            return True
        if is_array_like(xyz) and len(xyz) == 3:
            self.x, self.y, self.z = xyz
            return True
        return False

    def distance(self, other):
        """
        Compute the distance between two points.

        :param other: The other point.
        :type other: point_like

        :return: Distance to the other point.
        :rtype: float
        """
        if isinstance(other, gp_Pnt):
            return self.Distance(other)
        if is_array_like(other) and len(other) == 3:
            other = Point(*other)
            return self.Distance(other)
        return None

    def is_equal(self, other, tol=None):
        """
        Check for coincident points.

        :param other: The other point.
        :type other: point_like
        :param float tol: Tolerance for coincidence.

        :return: *True* if coincident, *False* if not.
        :rtype: bool
        """
        if tol is None:
            tol = Settings.gtol
        if isinstance(other, gp_Pnt):
            return self.IsEqual(other, tol)
        if is_array_like(other) and len(other) == 3:
            other = Point(*other)
            return self.IsEqual(other, tol)
        return False

    def translate(self, v):
        """
        Translate the point along the vector.

        :param v:

        :return:
        """
        self.Translate(v)


class Point2D(gp_Pnt2d, Geometry):
    """
    2-D point.
    """

    def __init__(self, *args):
        super(Point2D, self).__init__(*args)
        Geometry.__init__(self)

    def __str__(self):
        return 'Point 2-D = ({0}, {1})'.format(*self.xy)

    def __array__(self, dtype=float64, copy=True, order=None, subok=False,
                  ndmin=0):
        return array(self.xy, dtype=dtype, copy=copy, order=order,
                     subok=subok, ndmin=ndmin)

    def __iter__(self):
        for elm in self.xy:
            yield elm

    def __len__(self):
        return 2

    def __getitem__(self, item):
        return self.xy[item]

    def __add__(self, other):
        return add(self, other)

    def __sub__(self, other):
        return subtract(self, other)

    @property
    def xy(self):
        return array([self.X(), self.Y()], dtype=float64)

    @property
    def x(self):
        return self.X()

    @x.setter
    def x(self, x):
        self.SetX(x)

    @property
    def y(self):
        return self.Y()

    @y.setter
    def y(self, y):
        self.SetY(y)

    def copy(self):
        """
        Return a new copy of the point.

        :return:
        """
        return Point2D(*self.xy)

    def set_xy(self, xy):
        """
        Set point coordinates.

        :param array_like xy: Point coordinates.

        :return: *True* if set, *False* if not.
        :rtype: bool
        """
        if isinstance(xy, gp_Pnt):
            self.SetXY(xy.XYZ())
            return True
        if isinstance(xy, gp_XYZ):
            self.SetXY(xy)
            return True
        if is_array_like(xy) and len(xy) == 2:
            self.x, self.y = xy
            return True
        return False

    def distance(self, other):
        """
        Compute the distance between two points.

        :param other: The other point.
        :type other: point_like

        :return: Distance to the other point.
        :rtype: float
        """
        if isinstance(other, gp_Pnt2d):
            return self.Distance(other)
        if is_array_like(other) and len(other) == 2:
            other = Point2D(*other)
            return self.Distance(other)
        return None

    def is_equal(self, other, tol=None):
        """
        Check for coincident points.

        :param other: The other point.
        :type other: point_like
        :param float tol: Tolerance for coincidence.

        :return: *True* if coincident, *False* if not.
        :rtype: bool
        """
        if tol is None:
            tol = Settings.gtol
        if isinstance(other, gp_Pnt2d):
            return self.IsEqual(other, tol)
        if is_array_like(other) and len(other) == 2:
            other = Point2D(*other)
            return self.IsEqual(other, tol)
        return False
