from OCC.BRepGProp import brepgprop_LinearProperties, \
    brepgprop_SurfaceProperties, brepgprop_VolumeProperties
from OCC.GProp import GProp_GProps

from afem.geometry.entities import Point
from afem.topology.check import CheckShape

__all__ = ["ShapeProps", "LinearProps", "SurfaceProps", "VolumeProps",
           "LengthOfEdges", "AreaOfFaces"]


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


class LengthOfEdges(object):
    """
    Calculate the length of each edge and sort the results.

    :param list[OCC.TopoDS.TopoDS_Edge] edges: The edges.
    """

    def __init__(self, edges):
        results = []
        for e in edges:
            e = CheckShape.to_edge(e)
            le = LinearProps(e).length
            results.append((le, e))

        results.sort(key=lambda tup: tup[0])
        self._lengths = [data[0] for data in results]
        self._edges = [data[1] for data in results]

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
    def shortest_edge(self):
        """
        :return: The shortest edge.
        :rtype: OCC.TopoDS.TopoDS_Edge
        """
        return self._edges[0]

    @property
    def longest_edge(self):
        """
        :return: The longest edge.
        :rtype: OCC.TopoDS.TopoDS_Edge
        """
        return self._edges[-1]

    @property
    def sorted_edges(self):
        """
        :return: List of edges sorted by length.
        :rtype: list[OCC.TopoDS.TopoDS_Face]
        """
        return self._edges


class AreaOfFaces(object):
    """
    Calculate the area of each face and sort the results.

    :param list[OCC.TopoDS.TopoDS_Face] faces: The faces.
    """

    def __init__(self, faces):
        results = []
        for f in faces:
            f = CheckShape.to_face(f)
            a = SurfaceProps(f).area
            results.append((a, f))

        results.sort(key=lambda tup: tup[0])
        self._areas = [data[0] for data in results]
        self._faces = [data[1] for data in results]

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
    def smallest_face(self):
        """
        :return: The smallest face.
        :rtype: OCC.TopoDS.TopoDS_Face
        """
        return self._faces[0]

    @property
    def largest_face(self):
        """
        :return: The largest face.
        :rtype: OCC.TopoDS.TopoDS_Face
        """
        return self._faces[-1]

    @property
    def sorted_faces(self):
        """
        :return: List of faces sorted by area.
        :rtype: list[OCC.TopoDS.TopoDS_Face]
        """
        return self._faces


if __name__ == "__main__":
    import doctest

    doctest.testmod()
