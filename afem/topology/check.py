from OCC.BRep import (BRep_Builder, BRep_Tool_SameParameter,
                      BRep_Tool_SameRange)
from OCC.BRepBuilderAPI import (BRepBuilderAPI_MakeEdge,
                                BRepBuilderAPI_MakeFace,
                                BRepBuilderAPI_MakeVertex,
                                BRepBuilderAPI_MakeWire)
from OCC.BRepCheck import BRepCheck_Analyzer
from OCC.ShapeAnalysis import ShapeAnalysis_Edge
from OCC.TopAbs import (TopAbs_COMPOUND, TopAbs_COMPSOLID, TopAbs_EDGE,
                        TopAbs_FACE, TopAbs_SHELL, TopAbs_SOLID, TopAbs_VERTEX,
                        TopAbs_WIRE)
from OCC.TopoDS import (TopoDS_CompSolid, TopoDS_Compound, TopoDS_Edge,
                        TopoDS_Face, TopoDS_Shape, TopoDS_Shell, TopoDS_Solid,
                        TopoDS_Vertex, TopoDS_Wire, topods_CompSolid,
                        topods_Compound, topods_Edge, topods_Face,
                        topods_Shell, topods_Solid, topods_Vertex, topods_Wire)
from OCC.gp import gp_Pnt

from afem.geometry.check import CheckGeom

__all__ = ["CheckShape"]


class CheckShape(object):
    """
    Check shape.
    """

    @staticmethod
    def is_shape(shape):
        """
        Check if entity is a shape.

        :param shape: The shape.

        :return: *True* if a shape, *False* if not.
        :rtype: bool
        """
        return isinstance(shape, TopoDS_Shape)

    @staticmethod
    def is_solid(shape):
        """
        Check if the shape is a solid.

        :param OCC.TopoDS.TopoDS_Shape shape: The shape.

        :return: *True* if a solid, *False* if not.
        :rtype: bool
        """
        try:
            return shape.ShapeType() == TopAbs_SOLID
        except AttributeError:
            return False

    @classmethod
    def is_valid(cls, shape):
        """
        Check the shape for errors.

        :param OCC.TopoDS.TopoDS_Shape shape: The shape.

        :return: *True* if valid, *False* if not.
        :rtype: bool
        """
        return BRepCheck_Analyzer(shape, True).IsValid()

    @staticmethod
    def is_seam(edge, face):
        """
        Check to see if the edge is a seam edge on the face.

        :param OCC.TopoDS.TopoDS_Edge edge: The edge.
        :param OCC.TopoDS.TopoDS_Face face: The face.

        :return: *True* if a seam, *False* if not.
        :rtype: bool
        """
        return ShapeAnalysis_Edge().IsSeam(edge, face)

    @classmethod
    def to_vertex(cls, entity):
        """
        Convert an entity to a vertex.

        :param entity: The entity.

        :return: A vertex.
        :rtype: OCC.TopoDS.TopoDS_Vertex

        :raise TypeError: If entity cannot be converted to a vertex.
        """
        if isinstance(entity, TopoDS_Vertex):
            return entity

        if isinstance(entity, gp_Pnt):
            return BRepBuilderAPI_MakeVertex(entity).Vertex()

        if CheckGeom.is_point_like(entity):
            p = gp_Pnt(entity[0], entity[1], entity[2])
            return BRepBuilderAPI_MakeVertex(p).Vertex()

        if cls.is_shape(entity) and entity.ShapeType() == TopAbs_VERTEX:
            return topods_Vertex(entity)

        raise TypeError('Failed to convert entity to a vertex.')

    @classmethod
    def to_edge(cls, entity):
        """
        Convert an entity to an edge.

        :param entity: The entity.

        :return: An edge.
        :rtype: OCC.TopoDS.TopoDS_Edge

        :raise TypeError: If entity cannot be converted to an edge.
        """
        if isinstance(entity, TopoDS_Edge):
            return entity

        if CheckGeom.is_curve_like(entity):
            return BRepBuilderAPI_MakeEdge(entity.handle).Edge()

        if cls.is_shape(entity) and entity.ShapeType() == TopAbs_EDGE:
            return topods_Edge(entity)

        raise TypeError('Failed to convert entity to an edge.')

    @classmethod
    def to_wire(cls, entity):
        """
        Convert an entity to a wire.

        :param entity: The entity.

        :return: A wire.
        :rtype: OCC.TopoDS.TopoDS_Wire

        :raise TypeError: If entity cannot be converted to a wire.
        """
        if isinstance(entity, TopoDS_Wire):
            return entity

        if CheckGeom.is_curve_like(entity):
            e = BRepBuilderAPI_MakeEdge(entity.handle).Edge()
            return BRepBuilderAPI_MakeWire(e).Wire()

        if cls.is_shape(entity) and entity.ShapeType() == TopAbs_EDGE:
            return BRepBuilderAPI_MakeWire(entity).Wire()

        if cls.is_shape(entity) and entity.ShapeType() == TopAbs_WIRE:
            return topods_Wire(entity)

        raise TypeError('Failed to convert entity to a wire.')

    @classmethod
    def to_face(cls, entity):
        """
        Convert an entity to a face.

        :param entity: The entity.

        :return: A face.
        :rtype: OCC.TopoDS.TopoDS_Face

        :raise TypeError: If entity cannot be converted to a face.
        """
        if isinstance(entity, TopoDS_Face):
            return entity

        if CheckGeom.is_surface_like(entity):
            return BRepBuilderAPI_MakeFace(entity.handle, 0.).Face()

        if cls.is_shape(entity) and entity.ShapeType() == TopAbs_FACE:
            return topods_Face(entity)

        raise TypeError('Failed to convert entity to a face.')

    @classmethod
    def to_shell(cls, entity):
        """
        Convert an entity to a shell.

        :param entity: The entity.

        :return: A shell.
        :rtype: OCC.TopoDS.TopoDS_Shell

        :raise TypeError: If entity cannot be converted to a shell.
        """
        if isinstance(entity, TopoDS_Shell):
            return entity

        if cls.is_shape(entity) and entity.ShapeType() == TopAbs_SHELL:
            return topods_Shell(entity)

        if cls.is_shape(entity) and entity.ShapeType() == TopAbs_FACE:
            shell = TopoDS_Shell()
            builder = BRep_Builder()
            builder.MakeShell(shell)
            builder.Add(shell, entity)
            return shell

        if CheckGeom.is_surface_like(entity):
            f = BRepBuilderAPI_MakeFace(entity.handle, 0.).Face()
            shell = TopoDS_Shell()
            builder = BRep_Builder()
            builder.MakeShell(shell)
            builder.Add(shell, f)
            return shell

        raise TypeError('Failed to convert entity to a shell.')

    @classmethod
    def to_solid(cls, entity):
        """
        Convert an entity to a solid.

        :param entity: The entity.

        :return: A solid.
        :rtype: OCC.TopoDS.TopoDS_Solid

        :raise TypeError: If entity cannot be converted to a solid.
        """
        if isinstance(entity, TopoDS_Solid):
            return entity

        if cls.is_shape(entity) and entity.ShapeType() == TopAbs_SOLID:
            return topods_Solid(entity)

        raise TypeError('Failed to convert entity to a solid.')

    @classmethod
    def to_compsolid(cls, entity):
        """
        Convert an entity to a compsolid.

        :param entity: The entity.

        :return: A compsolid.
        :rtype: OCC.TopoDS.TopoDS_CompSolid

        :raise TypeError: If entity cannot be converted to a compsolid.
        """
        if isinstance(entity, TopoDS_CompSolid):
            return entity

        if cls.is_shape(entity) and entity.ShapeType() == TopAbs_COMPSOLID:
            return topods_CompSolid(entity)

        raise TypeError('Failed to convert entity to a compsolid.')

    @classmethod
    def to_compound(cls, entity):
        """
        Convert an entity to a compound.

        :param entity: The entity.

        :return: A compound
        :rtype: OCC.TopoDS.TopoDS_Compound

        :raise TypeError: If entity cannot be converted to a compound.
        """
        if isinstance(entity, TopoDS_Compound):
            return entity

        if cls.is_shape(entity) and entity.ShapeType() == TopAbs_COMPOUND:
            return topods_Compound(entity)

        if cls.is_shape(entity):
            cp = TopoDS_Compound
            builder = BRep_Builder()
            builder.MakeCompound(cp)
            builder.Add(cp, entity)
            return cp

        raise TypeError('Failed to convert entity to a compound.')

    @classmethod
    def to_shape(cls, entity):
        """
        Convert the entity to a shape. This method tries to convert the
        entity to its most specific shape type.

        :param entity: The entity.

        :return: A shape.
        :rtype: OCC.TopoDS.TopoDS_Shape

        :raise TypeError: If entity cannot be converted to a shape.
        """
        if entity is None:
            return None

        # Shapes
        if isinstance(entity, TopoDS_Shape):
            if entity.ShapeType() == TopAbs_VERTEX:
                return topods_Vertex(entity)
            elif entity.ShapeType() == TopAbs_EDGE:
                return topods_Edge(entity)
            elif entity.ShapeType() == TopAbs_WIRE:
                return topods_Wire(entity)
            elif entity.ShapeType() == TopAbs_FACE:
                return topods_Face(entity)
            elif entity.ShapeType() == TopAbs_SHELL:
                return topods_Shell(entity)
            elif entity.ShapeType() == TopAbs_SOLID:
                return topods_Solid(entity)
            elif entity.ShapeType() == TopAbs_COMPSOLID:
                return topods_CompSolid(entity)
            elif entity.ShapeType() == TopAbs_COMPOUND:
                return topods_Compound(entity)
            else:
                raise TypeError('Failed to convert entity to a shape.')

        # Geometry
        if CheckGeom.is_point_like(entity):
            return cls.to_vertex(entity)
        if CheckGeom.is_curve_like(entity):
            return cls.to_edge(entity)
        if CheckGeom.is_surface_like(entity):
            return cls.to_face(entity)

        raise TypeError('Failed to convert entity to a shape.')

    @staticmethod
    def same_parameter(edge):
        """
        Returns the SameParameter flag for the edge.

        :param OCC.TopoDS.TopoDS_Edge edge: The edge.

        :return: The same parameter flag.
        :rtype: bool
        """
        return BRep_Tool_SameParameter(edge)

    @staticmethod
    def same_range(edge):
        """
        Returns the SameRange flag for the edge.

        :param OCC.TopoDS.TopoDS_Edge edge: The edge.

        :return: The same range flag.
        :rtype: bool
        """
        return BRep_Tool_SameRange(edge)
