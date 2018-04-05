# This file is part of AFEM which provides an engineering toolkit for airframe
# finite element modeling during conceptual design.
#
# Copyright (C) 2016-2018  Laughlin Research, LLC (info@laughlinresearch.com)
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
from math import sqrt

from OCCT.BRepBndLib import BRepBndLib
from OCCT.Bnd import Bnd_Box
from OCCT.TopAbs import TopAbs_ShapeEnum
from OCCT.TopExp import TopExp_Explorer
from OCCT.TopoDS import TopoDS

from afem.geometry.check import CheckGeom
from afem.geometry.entities import Point
from afem.graphics.display import ViewableItem

__all__ = ["Shape", "Vertex", "Edge", "Wire", "Face", "Shell", "Solid",
           "Compound", "CompSolid",
           "BBox"]


class Shape(ViewableItem):
    """
    Shape.

    :param OCCT.TopoDS.TopoDS_Shape shape: The underlying shape.
    """

    VERTEX = TopAbs_ShapeEnum.TopAbs_VERTEX
    EDGE = TopAbs_ShapeEnum.TopAbs_EDGE
    WIRE = TopAbs_ShapeEnum.TopAbs_WIRE
    FACE = TopAbs_ShapeEnum.TopAbs_FACE
    SHELL = TopAbs_ShapeEnum.TopAbs_SHELL
    SOLID = TopAbs_ShapeEnum.TopAbs_SOLID
    COMPOUND = TopAbs_ShapeEnum.TopAbs_COMPOUND
    COMPSOLID = TopAbs_ShapeEnum.TopAbs_COMPSOLID

    def __init__(self, shape):
        super(Shape, self).__init__()

        self._shape = shape

    def __hash__(self):
        """
        Use the hash code of the shape.
        """
        return hash(self.hash_code)

    @property
    def object(self):
        """
        :return: The underlying shape.
        :rtype: OCCT.TopoDS.TopoDS_Shape
        """
        return self._shape

    @property
    def hash_code(self):
        """
        :return: The hash code of the shape computed using the TShape and
            Location. Orientation is not used. The upper limit is 99,999.
        :rtype: int
        """
        return self.object.HashCode(99999)

    @property
    def is_null(self):
        """
        :return: *True* if shape is null, *False* if not.
        :rtype: bool
        """
        return self.object.IsNull()

    @property
    def shape_type(self):
        """
        :return: The shape type.
        :rtype: OCCT.TopAbs.TopAbs_ShapeEnum
        """
        return self.object.ShapeType()

    @property
    def is_vertex(self):
        """
        :return: *True* if a Vertex, *False* if not.
        :rtype: bool
        """
        return self.shape_type == TopAbs_ShapeEnum.TopAbs_VERTEX

    @property
    def is_edge(self):
        """
        :return: *True* if an Edge, *False* if not.
        :rtype: bool
        """
        return self.shape_type == TopAbs_ShapeEnum.TopAbs_EDGE

    @property
    def is_wire(self):
        """
        :return: *True* if a Wire, *False* if not.
        :rtype: bool
        """
        return self.shape_type == TopAbs_ShapeEnum.TopAbs_WIRE

    @property
    def is_face(self):
        """
        :return: *True* if a Face, *False* if not.
        :rtype: bool
        """
        return self.shape_type == TopAbs_ShapeEnum.TopAbs_FACE

    @property
    def is_shell(self):
        """
        :return: *True* if a Shell, *False* if not.
        :rtype: bool
        """
        return self.shape_type == TopAbs_ShapeEnum.TopAbs_SHELL

    @property
    def is_solid(self):
        """
        :return: *True* if a Solid, *False* if not.
        :rtype: bool
        """
        return self.shape_type == TopAbs_ShapeEnum.TopAbs_SOLID

    @property
    def is_compound(self):
        """
        :return: *True* if a Compound, *False* if not.
        :rtype: bool
        """
        return self.shape_type == TopAbs_ShapeEnum.TopAbs_COMPOUND

    @property
    def is_compsolid(self):
        """
        :return: *True* if a CompSolid, *False* if not.
        :rtype: bool
        """
        return self.shape_type == TopAbs_ShapeEnum.TopAbs_COMPSOLID

    @property
    def closed(self):
        """
        :return: The closed flag of the shape.
        :rtype: bool
        """
        return self.object.Closed()

    @property
    def infinite(self):
        """
        :return: The infinite flag of the shape.
        :rtype: bool
        """
        return self.object.Infinite()

    def reverse(self):
        """
        Reverse the orientation of the shape.

        :return: None.
        """
        self.object.Reverse()

    def is_partner(self, other):
        """
        Check if this shape shares the same TShape with the other. Locations
        and Orientations may differ.

        :param afem.topology.entities.Shape other: Other shape.

        :return: *True* if partner, *False* if not.
        :rtype: bool
        """
        return self.object.IsPartner(other.object)

    def is_same(self, other):
        """
        Check if this shape shares the same TShape and Location with the other.
        Orientations may differ.

        :param afem.topology.entities.Shape other: Other shape.

        :return: *True* if same, *False* if not.
        :rtype: bool
        """
        return self.object.IsSame(other.object)

    def is_equal(self, other):
        """
        Check if this shape is equal to the other. That is, they share the same
        TShape, Location, and Orientation.

        :param afem.topology.entities.Shape other: Other shape.

        :return: *True* if equal, *False* if not.
        :rtype: bool
        """
        return self.object.IsEqual(other.object)

    def downcasted(self):
        """
        Convert this shape to a more specific type based on its shape type.

        :return: The new shape. Returns itself if not converted.
        :rtype: afem.topology.entities.Shape
        """
        if self.is_vertex:
            v = TopoDS.Vertex_(self.object)
            return Vertex(v)
        if self.is_edge:
            e = TopoDS.Edge_(self.object)
            return Edge(e)
        if self.is_wire:
            w = TopoDS.Wire_(self.object)
            return Wire(w)
        if self.is_face:
            f = TopoDS.Face_(self.object)
            return Face(f)
        if self.is_shell:
            s = TopoDS.Shell_(self.object)
            return Shell(s)
        if self.is_solid:
            s = TopoDS.Solid_(self.object)
            return Solid(s)
        if self.is_compound:
            c = TopoDS.Compound_(self.object)
            return Compound(c)
        if self.is_compsolid:
            c = TopoDS.CompSolid_(self.object)
            return CompSolid(c)
        return self

    def get_shapes(self, type_):
        """
        Get sub-shapes of a specified type from the shape.

        :param OCCT.TopAbs.TopAbs_ShapeEnum type_: The shape type.

        :return: List of sub-shapes.
        :rtype: list(afem.topology.entities.Shape)
        """
        explorer = TopExp_Explorer(self.object, type_)
        shapes = []
        while explorer.More():
            si = Shape(explorer.Current()).downcasted()
            is_unique = True
            for s in shapes:
                if s.is_same(si):
                    is_unique = False
                    break
            if is_unique:
                shapes.append(si)
            explorer.Next()
        return shapes


class Vertex(Shape):
    """
    Vertex.
    """

    def __init__(self, vertex):
        super(Vertex, self).__init__(vertex)


class Edge(Shape):
    """
    Edge.
    """

    def __init__(self, edge):
        super(Edge, self).__init__(edge)


class Wire(Shape):
    """
    Wire.
    """

    def __init__(self, wire):
        super(Wire, self).__init__(wire)


class Face(Shape):
    """
    Face.
    """

    def __init__(self, face):
        super(Face, self).__init__(face)


class Shell(Shape):
    """
    Shell.
    """

    def __init__(self, shell):
        super(Shell, self).__init__(shell)


class Solid(Shape):
    """
    Solid.
    """

    def __init__(self, solid):
        super(Solid, self).__init__(solid)


class Compound(Shape):
    """
    Compound.
    """

    def __init__(self, compound):
        super(Compound, self).__init__(compound)


class CompSolid(Shape):
    """
    CompSolid.
    """

    def __init__(self, compsolid):
        super(CompSolid, self).__init__(compsolid)


class BBox(Bnd_Box):
    """
    Bounding box in 3-D space.

    For more information see Bnd_Box_.

    .. _Bnd_Box: https://www.opencascade.com/doc/occt-7.2.0/refman/html/class_bnd___box.html

    Usage:

    >>> from afem.topology import *
    >>> e = EdgeByPoints((0., 0., 0.), (10., 0., 0.)).edge
    >>> bbox = BBox()
    >>> bbox.add_shape(e)
    >>> bbox.set_gap(0.)
    >>> bbox.gap
    0.0
    >>> bbox.pmin
    Point(0.000, 0.000, 0.000)
    >>> bbox.pmax
    Point(10.000, 0.000, 0.000)
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
        :rtype: afem.geometry.entities.Point
        """
        if self.is_void:
            return None
        return Point(self.CornerMin().XYZ())

    @property
    def pmax(self):
        """
        :return: Upper corner of bounding box. *None* if empty.
        :rtype: afem.geometry.entities.Point
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

        :param OCCT.TopoDS.TopoDS_Shape shape: The shape.

        :return: None.
        """
        BRepBndLib.Add_(shape, self, True)

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

        :param afem.geometry.entities.Line line: The line.

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

        :param afem.geometry.entities.Plane pln: The plane.

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

        :param afem.topology.entities.BBox bbox: The other box.

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
