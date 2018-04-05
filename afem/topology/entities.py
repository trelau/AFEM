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
from itertools import product
from math import sqrt

from OCCT.BRep import BRep_Tool
from OCCT.BRepBndLib import BRepBndLib
from OCCT.BRepBuilderAPI import BRepBuilderAPI_Copy
from OCCT.BRepClass3d import BRepClass3d
from OCCT.BRepTools import BRepTools
from OCCT.Bnd import Bnd_Box
from OCCT.ShapeAnalysis import ShapeAnalysis_Edge, ShapeAnalysis_ShapeTolerance
from OCCT.TopAbs import TopAbs_ShapeEnum
from OCCT.TopExp import TopExp_Explorer
from OCCT.TopoDS import (TopoDS, TopoDS_Vertex, TopoDS_Edge, TopoDS_Wire,
                         TopoDS_Face, TopoDS_Shell, TopoDS_Solid,
                         TopoDS_Compound, TopoDS_CompSolid)

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

    SHAPE = TopAbs_ShapeEnum.TopAbs_SHAPE
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

        # The underlying OCCT shape
        self._shape = shape

    def __hash__(self):
        """
        Use the hash code of the shape.
        """
        return hash(self.hash_code)

    def __eq__(self, other):
        """
        Check equality using is_same.
        """
        return self.is_same(other)

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
        return self.shape_type == Shape.VERTEX

    @property
    def is_edge(self):
        """
        :return: *True* if an Edge, *False* if not.
        :rtype: bool
        """
        return self.shape_type == Shape.EDGE

    @property
    def is_wire(self):
        """
        :return: *True* if a Wire, *False* if not.
        :rtype: bool
        """
        return self.shape_type == Shape.WIRE

    @property
    def is_face(self):
        """
        :return: *True* if a Face, *False* if not.
        :rtype: bool
        """
        return self.shape_type == Shape.FACE

    @property
    def is_shell(self):
        """
        :return: *True* if a Shell, *False* if not.
        :rtype: bool
        """
        return self.shape_type == Shape.SHELL

    @property
    def is_solid(self):
        """
        :return: *True* if a Solid, *False* if not.
        :rtype: bool
        """
        return self.shape_type == Shape.SOLID

    @property
    def is_compound(self):
        """
        :return: *True* if a Compound, *False* if not.
        :rtype: bool
        """
        return self.shape_type == Shape.COMPOUND

    @property
    def is_compsolid(self):
        """
        :return: *True* if a CompSolid, *False* if not.
        :rtype: bool
        """
        return self.shape_type == Shape.COMPSOLID

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

    @property
    def vertices(self):
        """
        :return: The vertices of the shape.
        :rtype: list(afem.topology.entities.Vertex)
        """
        return self._get_shapes(self.VERTEX)

    @property
    def edges(self):
        """
        :return: The edges of the shape.
        :rtype: list(afem.topology.entities.Edge)
        """
        return self._get_shapes(self.EDGE)

    @property
    def wires(self):
        """
        :return: The wires of the shape.
        :rtype: list(afem.topology.entities.Wire)
        """
        return self._get_shapes(self.WIRE)

    @property
    def faces(self):
        """
        :return: The Face of the shape.
        :rtype: list(afem.topology.entities.Face)
        """
        return self._get_shapes(self.FACE)

    @property
    def shells(self):
        """
        :return: The shells of the shape.
        :rtype: list(afem.topology.entities.Shell)
        """
        return self._get_shapes(self.SHELL)

    @property
    def solids(self):
        """
        :return: The solids of the shape.
        :rtype: list(afem.topology.entities.Solid)
        """
        return self._get_shapes(self.SOLID)

    @property
    def compounds(self):
        """
        :return: The compounds of the shape.
        :rtype: list(afem.topology.entities.Compound)
        """
        return self._get_shapes(self.COMPOUND)

    @property
    def compsolids(self):
        """
        :return: The compsolids of the shape.
        :rtype: list(afem.topology.entities.CompSolid)
        """
        return self._get_shapes(self.COMPSOLID)

    @property
    def tol(self):
        """
        :return: The average global tolerance.
        :rtype: float
        """
        tol = ShapeAnalysis_ShapeTolerance()
        tol.AddTolerance(self.object)
        return tol.GlobalTolerance(0)

    @property
    def tol_min(self):
        """
        :return: The minimum global tolerance.
        :rtype: float
        """
        tol = ShapeAnalysis_ShapeTolerance()
        tol.AddTolerance(self.object)
        return tol.GlobalTolerance(-1)

    @property
    def tol_max(self):
        """
        :return: The minimum global tolerance.
        :rtype: float
        """
        tol = ShapeAnalysis_ShapeTolerance()
        tol.AddTolerance(self.object)
        return tol.GlobalTolerance(1)

    def _get_shapes(self, type_):
        """
        Get sub-shapes of a specified type from the shape.
        """
        explorer = TopExp_Explorer(self.object, type_)
        shapes = []
        while explorer.More():
            si = Shape.wrap(explorer.Current())
            is_unique = True
            for s in shapes:
                if s.is_same(si):
                    is_unique = False
                    break
            if is_unique:
                shapes.append(si)
            explorer.Next()
        return shapes

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

    def copy(self, geom=True):
        """
        Copy this shape.

        :param bool geom: Option to copy geometry.

        :return: The copied shape.
        :rtype: afem.topology.entities.Shape
        """
        return Shape.wrap(BRepBuilderAPI_Copy(self.object, geom).Shape())

    def shared_vertices(self, other):
        """
        Get shared vertices between this shape and the other.

        :param afem.topology.entities.Shape other: The other shape.

        :return: Shared vertices.
        :rtype: list(afem.topology.entities.Vertex)
        """
        verts1 = self.vertices
        verts2 = other.vertices
        if not verts1 or not verts2:
            return []

        shared_verts = []
        for v1, v2 in product(verts1, verts2):
            if v1.is_same(v2):
                unique = True
                for vi in shared_verts:
                    if vi.is_same(v1):
                        unique = False
                        break
                if unique:
                    shared_verts.append(v1)

        return shared_verts

    def shared_edges(self, other):
        """
        Get shared edges between this shape and the other.

        :param afem.topology.entities.Shape other: The other shape.

        :return: Shared edges.
        :rtype: list(afem.topology.entities.Edge)
        """
        edges1 = self.edges
        edges2 = other.edges
        if not edges1 or not edges2:
            return []

        shared_edges = []
        for e1, e2 in product(edges1, edges2):
            if e1.is_same(e2):
                unique = True
                for ei in shared_edges:
                    if ei.is_same(e1):
                        unique = False
                        break
                if unique:
                    shared_edges.append(e1)

        return shared_edges

    @staticmethod
    def wrap(shape):
        """
        Convert the shape to a more specific type based on its shape type.

        :param OCCT.TopoDS.TopoDS_Shape shape: The shape.

        :return: The new shape.
        :rtype: afem.topology.entities.Shape
        """
        if shape.ShapeType() == Shape.VERTEX:
            return Vertex(shape)
        if shape.ShapeType() == Shape.EDGE:
            return Edge(shape)
        if shape.ShapeType() == Shape.WIRE:
            return Wire(shape)
        if shape.ShapeType() == Shape.FACE:
            return Face(shape)
        if shape.ShapeType() == Shape.SHELL:
            return Shell(shape)
        if shape.ShapeType() == Shape.SOLID:
            return Solid(shape)
        if shape.ShapeType() == Shape.COMPOUND:
            return Compound(shape)
        if shape.ShapeType() == Shape.COMPSOLID:
            return CompSolid(shape)

        return Shape(shape)


class Vertex(Shape):
    """
    Vertex.

    :param OCCT.TopoDS.TopoDS_Vertex vertex: The vertex.
    """

    def __init__(self, vertex):
        if not vertex.ShapeType() == Shape.VERTEX:
            raise TypeError('Shape is not a TopoDS_Vertex.')
        if not isinstance(vertex, TopoDS_Vertex):
            vertex = TopoDS.Vertex_(vertex)
        super(Vertex, self).__init__(vertex)

    @property
    def point(self):
        """
        :return: The vertex point.
        :rtype: afem.geometry.entities.Point
        """
        gp_pnt = BRep_Tool.Pnt_(self.object)
        return Point(gp_pnt.X(), gp_pnt.Y(), gp_pnt.Z())

    def parameter(self, edge, face=None):
        pass


class Edge(Shape):
    """
    Edge.

    :param OCCT.TopoDS.TopoDS_Edge edge: The edge.
    """

    def __init__(self, edge):
        if not edge.ShapeType() == Shape.EDGE:
            raise TypeError('Shape is not a TopoDS_Edge.')
        if not isinstance(edge, TopoDS_Edge):
            edge = TopoDS.Edge_(edge)
        super(Edge, self).__init__(edge)

    @property
    def curve(self):
        pass

    @property
    def first_vertex(self):
        """
        :return: The first vertex of the edge.
        :rtype: afem.topology.entities.Vertex
        """
        return Vertex(ShapeAnalysis_Edge().FirstVertex(self.object))

    @property
    def last_vertex(self):
        """
        :return: The last vertex of the edge.
        :rtype: afem.topology.entities.Vertex
        """
        return Vertex(ShapeAnalysis_Edge().LastVertex(self.object))


class Wire(Shape):
    """
    Wire.

    :param OCCT.TopoDS.TopoDS_Wire wire: The wire.
    """

    def __init__(self, wire):
        if not wire.ShapeType() == Shape.WIRE:
            raise TypeError('Shape is not a TopoDS_Wire.')
        if not isinstance(wire, TopoDS_Wire):
            wire = TopoDS.Wire_(wire)
        super(Wire, self).__init__(wire)


class Face(Shape):
    """
    Face.

    :param OCCT.TopoDS.TopoDS_Face face: The face.
    """

    def __init__(self, face):
        if not face.ShapeType() == Shape.FACE:
            raise TypeError('Shape is not a TopoDS_Face.')
        if not isinstance(face, TopoDS_Face):
            face = TopoDS.Face_(face)
        super(Face, self).__init__(face)

    @property
    def surface(self):
        pass

    @property
    def outer_wire(self):
        """
        :return: The outer wire of the face.
        :rtype: afem.topology.entities.Wire
        """
        return Wire(BRepTools.OuterWire_(self.object))


class Shell(Shape):
    """
    Shell.

    :param OCCT.TopoDS.TopoDS_Shell shell: The shell.
    """

    def __init__(self, shell):
        if not shell.ShapeType() == Shape.SHELL:
            raise TypeError('Shape is not a TopoDS_Shell.')
        if not isinstance(shell, TopoDS_Shell):
            shell = TopoDS.Shell_(shell)
        super(Shell, self).__init__(shell)


class Solid(Shape):
    """
    Solid.

    :param OCCT.TopoDS.TopoDS_Solid solid: The solid.
    """

    def __init__(self, solid):
        if not solid.ShapeType() == Shape.SOLID:
            raise TypeError('Shape is not a TopoDS_Solid.')
        if not isinstance(solid, TopoDS_Solid):
            solid = TopoDS.Solid_(solid)
        super(Solid, self).__init__(solid)

    @property
    def outer_shell(self):
        """
        :return: The outer shell of the face.
        :rtype: afem.topology.entities.Shell
        """
        return Shell(BRepClass3d.OuterShell_(self.object))


class Compound(Shape):
    """
    Compound.

    :param OCCT.TopoDS.TopoDS_Compound compound: The compound.
    """

    def __init__(self, compound):
        if not compound.ShapeType() == Shape.COMPOUND:
            raise TypeError('Shape is not a TopoDS_Compound.')
        if not isinstance(compound, TopoDS_Compound):
            compound = TopoDS.Compound_(compound)
        super(Compound, self).__init__(compound)


class CompSolid(Shape):
    """
    CompSolid.

    :param OCCT.TopoDS.TopoDS_CompSolid compsolid: The compsolid.
    """

    def __init__(self, compsolid):
        if not compsolid.ShapeType() == Shape.COMPSOLID:
            raise TypeError('Shape is not a TopoDS_CompSolid.')
        if not isinstance(compsolid, TopoDS_CompSolid):
            compsolid = TopoDS.CompSolid_(compsolid)
        super(CompSolid, self).__init__(compsolid)


class BBox(Bnd_Box):
    """
    Bounding box in 3-D space.

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
