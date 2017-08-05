from OCC.BRepGProp import brepgprop_LinearProperties, \
    brepgprop_SurfaceProperties, brepgprop_VolumeProperties
from OCC.GProp import GProp_GProps

from afem.geometry.entities import Point

__all__ = ["ShapeProps", "LinearProps", "SurfaceProps", "VolumeProps"]


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
        :rtype: afem.geometry.entities.Point
        """
        gp_pnt = self._props.CentreOfMass()
        return Point(gp_pnt.X(), gp_pnt.Y(), gp_pnt.Z())


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


if __name__ == "__main__":
    import doctest

    doctest.testmod()
