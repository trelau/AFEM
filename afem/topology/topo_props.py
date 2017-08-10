from OCC.BRepGProp import (brepgprop_LinearProperties,
                           brepgprop_SurfaceProperties,
                           brepgprop_VolumeProperties)
from OCC.GProp import GProp_GProps
from numpy import array

from afem.geometry.geom_entities import Point
from afem.topology.topo_check import CheckShape

__all__ = ["ShapeProps", "LinearProps", "SurfaceProps", "VolumeProps",
           "LengthOfShapes", "AreaOfShapes"]


class ShapeProps(object):
    """
    Base class for shape properties.
    """

    def __init__(self):
        self._props = GProp_GProps()

    @property
    def mass(self):
        """
        :return: The mass of the shape. This corresponds to total length for
            linear properties, total area for surface properties, or total
            volume for volume properties.
        :rtype: float
        """
        return self._props.Mass()

    @property
    def cg(self):
        """
        :return: The center of gravity.
        :rtype: afem.geometry.geom_entities.Point
        """
        gp_pnt = self._props.CentreOfMass()
        return Point(gp_pnt.X(), gp_pnt.Y(), gp_pnt.Z())

    @property
    def static_moments(self):
        """
        :return: The static moments of inertia Ix, Iy, and Iz.
        :rtype: tuple(float)
        """
        return self._props.StaticMoments()

    @property
    def matrix_of_inertia(self):
        """
        :return: The 3 x 3 matrix of inertia.
        :rtype: numpy.ndarray
        """
        gp_mat = self._props.MatrixOfInertia()
        matrix = []
        for j in range(1, 4):
            row = []
            for i in range(1, 4):
                row.append(gp_mat.Value(i, j))
            matrix.append(row)
        return array(matrix, dtype=float)

    def moment_of_inertia(self, axis):
        """
        Compute the moment of inertia about the axis.

        :param afem.geometry.geom_entities.Axis1 axis: The axis.

        :return: The moment of inertia.
        :rtype: float
        """
        return self._props.MomentOfInertia(axis)


class LinearProps(ShapeProps):
    """
    Calculate linear properties of a shape.

    :param OCC.TopoDS.TopoDS_Shape shape: The shape.

    Usage:

    >>> from afem.topology import EdgeByPoints, LinearProps
    >>> e = EdgeByPoints((0., 0., 0.), (1., 0., 0.)).edge
    >>> props = LinearProps(e)
    >>> props.length
    1.0
    >>> props.cg
    Point(0.5, 0.0, 0.0)
    """

    def __init__(self, shape):
        super(LinearProps, self).__init__()
        brepgprop_LinearProperties(shape, self._props)

    @property
    def length(self):
        """
        :return: The total length of all edges of the shape.
        :rtype: float
        """
        return self.mass


class SurfaceProps(ShapeProps):
    """
    Calculate surface properties of a shape.

    :param OCC.TopoDS.TopoDS_Shape shape: The shape.

    Usage:

    >>> from afem.topology import EdgeByPoints, FaceByDrag, SurfaceProps
    >>> e = EdgeByPoints((0., 0., 0.), (1., 0., 0.)).edge
    >>> f = FaceByDrag(e, (0., 1., 0.)).face
    >>> props = SurfaceProps(f)
    >>> props.area
    0.9999999999999998
    >>> props.cg
    Point(0.5, 0.5, 0.0)
    """

    def __init__(self, shape):
        super(SurfaceProps, self).__init__()
        brepgprop_SurfaceProperties(shape, self._props)

    @property
    def area(self):
        """
        :return: The total area of all faces of the shape.
        :rtype: float
        """
        return self.mass


class VolumeProps(ShapeProps):
    """
    Calculate volume properties of a shape.

    :param OCC.TopoDS.TopoDS_Shape shape: The shape.

    Usage:

    >>> from afem.topology import *
    >>> e = EdgeByPoints((0., 0., 0.), (1., 0., 0.)).edge
    >>> f = FaceByDrag(e, (0., 1., 0.)).face
    >>> solid = SolidByDrag(f, (0., 0., 1.)).solid
    >>> props = VolumeProps(solid)
    >>> props.volume
    0.9999999999999998
    >>> props.cg
    Point(0.5, 0.5, 0.5)
    """

    def __init__(self, shape):
        super(VolumeProps, self).__init__()
        brepgprop_VolumeProperties(shape, self._props)

    @property
    def volume(self):
        """
        :return: The total volume of all solids of the shape.
        :rtype: float
        """
        return self.mass


class LengthOfShapes(object):
    """
    Calculate the total length of all edges of each shape and sort the results.

    :param list[OCC.TopoDS.TopoDS_Shape] shapes: The shapes.
    """

    def __init__(self, shapes):
        results = []
        for shape in shapes:
            shape = CheckShape.to_shape(shape)
            ls = LinearProps(shape).length
            results.append((ls, shape))

        results.sort(key=lambda tup: tup[0])
        self._lengths = [data[0] for data in results]
        self._shapes = [data[1] for data in results]

    @property
    def min_length(self):
        """
        :return: The minimum length.
        :rtype: float
        """
        return self._lengths[0]

    @property
    def max_length(self):
        """
        :return: The maximum length.
        :rtype: float
        """
        return self._lengths[-1]

    @property
    def sorted_lengths(self):
        """
        :return: List of sorted lengths.
        :rtype: list[float]
        """
        return self._lengths

    @property
    def shortest_shape(self):
        """
        :return: The shortest shape.
        :rtype: OCC.TopoDS.TopoDS_Shape
        """
        return self._shapes[0]

    @property
    def longest_shape(self):
        """
        :return: The longest shape.
        :rtype: OCC.TopoDS.TopoDS_Shape
        """
        return self._shapes[-1]

    @property
    def sorted_shapes(self):
        """
        :return: List of shapes sorted by length.
        :rtype: list[OCC.TopoDS.TopoDS_Shape]
        """
        return self._shapes


class AreaOfShapes(object):
    """
    Calculate the total area of each face for each shape and sort the results.

    :param list[OCC.TopoDS.TopoDS_Shape] shapes: The shapes.
    """

    def __init__(self, shapes):
        results = []
        for shape in shapes:
            shape = CheckShape.to_shape(shape)
            a = SurfaceProps(shape).area
            results.append((a, shape))

        results.sort(key=lambda tup: tup[0])
        self._areas = [data[0] for data in results]
        self._shapes = [data[1] for data in results]

    @property
    def min_area(self):
        """
        :return: The minimum area.
        :rtype: float
        """
        return self._areas[0]

    @property
    def max_area(self):
        """
        :return: The maximum area.
        :rtype: float
        """
        return self._areas[-1]

    @property
    def sorted_areas(self):
        """
        :return: List of sorted areas.
        :rtype: list[float]
        """
        return self._areas

    @property
    def smallest_shape(self):
        """
        :return: The smallest shape.
        :rtype: OCC.TopoDS.TopoDS_Shape
        """
        return self._shapes[0]

    @property
    def largest_shape(self):
        """
        :return: The largest shape.
        :rtype: OCC.TopoDS.TopoDS_Shape
        """
        return self._shapes[-1]

    @property
    def sorted_shape(self):
        """
        :return: List of shapes sorted by area.
        :rtype: list[OCC.TopoDS.TopoDS_Shape]
        """
        return self._shapes


if __name__ == "__main__":
    import doctest

    doctest.testmod()
