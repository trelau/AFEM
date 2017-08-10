from math import sqrt

from OCC.BRepBndLib import brepbndlib_Add
from OCC.Bnd import Bnd_Box

from afem.geometry.geom_check import CheckGeom
from afem.geometry.geom_entities import Point

__all__ = ["BBox"]


class BBox(Bnd_Box):
    """
    Bounding box in 3-D space.

    For more information see Bnd_Box_.

    .. _Bnd_Box: https://www.opencascade.com/doc/occt-7.1.0/refman/html/class_bnd___box.html

    Usage:

    >>> from afem.topo_patch import *
    >>> e = EdgeByPoints((0., 0., 0.), (10., 0., 0.)).edge
    >>> bbox = BBox()
    >>> bbox.add_shape(e)
    >>> bbox.set_gap(0.)
    >>> bbox.gap
    0.0
    >>> bbox.pmin
    Point(0.0, 0.0, 0.0)
    >>> bbox.pmax
    Point(10.0, 0.0, 0.0)
    >>> bbox.diagonal
    10.0
    """

    def __init__(self):
        super(BBox, self).__init__()

    @property
    def is_void(self):
        """
        :return: *True* if bounding box is empty, *False* if not.
        :rtype: bool
        """
        return self.IsVoid()

    @property
    def pmin(self):
        """
        :return: Lower corner of bounding box. *None* if empty.
        :rtype: afem.geometry.geom_entities.Point
        """
        if self.is_void:
            return None
        return Point(self.CornerMin().XYZ())

    @property
    def pmax(self):
        """
        :return: Upper corner of bounding box. *None* if empty.
        :rtype: afem.geometry.geom_entities.Point
        """
        if self.is_void:
            return None
        return Point(self.CornerMax().XYZ())

    @property
    def xmin(self):
        """
        :return: Minimum x-component.
        :rtype: float
        """
        if self.is_void:
            return None
        return self.CornerMin().X()

    @property
    def xmax(self):
        """
        :return: Maximum x-component.
        :rtype: float
        """
        if self.is_void:
            return None
        return self.CornerMax().X()

    @property
    def ymin(self):
        """
        :return: Minimum y-component.
        :rtype: float
        """
        if self.is_void:
            return None
        return self.CornerMin().Y()

    @property
    def ymax(self):
        """
        :return: Maximum y-component.
        :rtype: float
        """
        if self.is_void:
            return None
        return self.CornerMax().Y()

    @property
    def zmin(self):
        """
        :return: Minimum z-component.
        :rtype: float
        """
        if self.is_void:
            return None
        return self.CornerMin().Z()

    @property
    def zmax(self):
        """
        :return: Maximum z-component.
        :rtype: float
        """
        if self.is_void:
            return None
        return self.CornerMax().Z()

    @property
    def gap(self):
        """
        :return: The gap of the bounding box.
        :rtype: float
        """
        return self.GetGap()

    @property
    def diagonal(self):
        """
        :return: The diagonal length of the box.
        :rtype: float
        """
        return sqrt(self.SquareExtent())

    def set_gap(self, gap):
        """
        Set the gap of the bounding box.

        :param float gap: The gap.

        :return: None.
        """
        self.SetGap(abs(gap))

    def enlarge(self, tol):
        """
        Enlarge the box with a tolerance value.

        :param float tol: The tolerance.

        :return: None.
        """
        self.Enlarge(tol)

    def add_box(self, bbox):
        """
        Add the other box to this one.

        :param afem.geometry.entities.BBox bbox: The other box.

        :return: None.

        :raise TypeError: If *bbox* cannot be converted to a bounding box.
        """
        if not isinstance(bbox, Bnd_Box):
            msg = 'Methods requires a BBox instance.'
            raise TypeError(msg)

        self.Add(bbox)

    def add_pnt(self, pnt):
        """
        Add a point to the bounding box.

        :param point_like pnt: The point.

        :return: None.

        :raise TypeError: If *pnt* cannot be converted to a point.
        """
        pnt = CheckGeom.to_point(pnt)
        if not pnt:
            msg = 'Invalid point type provided.'
            raise TypeError(msg)

        self.Add(pnt)

    def add_shape(self, shape):
        """
        Add shape to the bounding box.

        :param OCC.TopoDS.TopoDS_Shape shape: The shape.

        :return: None.
        """
        brepbndlib_Add(shape, self, True)

    def is_pnt_out(self, pnt):
        """
        Check to see if the point is outside the bounding box.

        :param point_like pnt: The point.

        :return: *True* if outside, *False* if not.
        :rtype: bool

        :raise TypeError: If *pnt* cannot be converted to a point.
        """
        pnt = CheckGeom.to_point(pnt)
        if not pnt:
            msg = 'Invalid point type provided.'
            raise TypeError(msg)

        return self.IsOut(pnt)

    def is_line_out(self, line):
        """
        Check to see if the line intersects the box.

        :param afem.geometry.geom_entities.Line line: The line.

        :return: *True* if outside, *False* if it intersects.
        :rtype: bool

        :raise TypeError: If *line* is not a line.
        """

        if not CheckGeom.is_line(line):
            msg = 'Methods requires a Line instance.'
            raise TypeError(msg)

        return self.IsOut(line)

    def is_pln_out(self, pln):
        """
        Check to see if the plane intersects the box.

        :param afem.geometry.geom_entities.Plane pln: The plane.

        :return: *True* if outside, *False* if it intersects.
        :rtype: bool

        :raise TypeError: If *pln* is not a plane.
        """

        if not CheckGeom.is_plane(pln):
            msg = 'Methods requires a Plane instance.'
            raise TypeError(msg)

        return self.IsOut(pln)

    def is_box_out(self, bbox):
        """
        Check to see if the bounding box intersects this one.

        :param afem.topology.topo_entities.BBox bbox: The other box.

        :return: *True* if outside, *False* if it intersects or is inside.
        :rtype: bool

        :raise TypeError: If *bbox* cannot be converted to a bounding box.
        """
        if not isinstance(bbox, Bnd_Box):
            msg = 'Methods requires a BBox instance.'
            raise TypeError(msg)

        return self.IsOut(bbox)

    def distance(self, bbox):
        """
        Calculate distance to other box.

        :param afem.geometry.entities.BBox bbox: The other box.

        :return: Distance to other box.
        :rtype: float

        :raise TypeError: If *bbox* cannot be converted to a bounding box.
        """
        if not isinstance(bbox, Bnd_Box):
            msg = 'Methods requires a BBox instance.'
            raise TypeError(msg)

        return self.Distance(bbox)


if __name__ == "__main__":
    import doctest

    doctest.testmod()
