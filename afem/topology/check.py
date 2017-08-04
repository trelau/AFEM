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

# TODO Fix import error
try:
    from afem.geometry.check import CheckGeom
except ImportError:
    pass

__all__ = ["CheckShape"]


class CheckShape(object):
    """
    Check shape.
    """

    @staticmethod
    def is_shape(shape):
        """
        Check if shape is a TopoDS_Shape.

        :param shape:

        :return:
        """
        return isinstance(shape, TopoDS_Shape)

    @staticmethod
    def is_solid(shape):
        """
        Check if the shape is a solid.

        :param shape:

        :return:
        """
        try:
            return shape.ShapeType() == TopAbs_SOLID
        except AttributeError:
            return False

    @classmethod
    def is_valid(cls, shape):
        """
        Check the shape for errors.

        :param shape:

        :return:
        """
        shape = cls.to_shape(shape)
        if not shape:
            return False
        return BRepCheck_Analyzer(shape, True).IsValid()

    @staticmethod
    def is_seam(edge, face):
        """
        Check to see if the edge is a seam edge on the face.

        :param edge:
        :param face:

        :return:
        """
        return ShapeAnalysis_Edge().IsSeam(edge, face)

    @classmethod
    def to_vertex(cls, entity):
        """
        Convert an entity to a vertex.

        :param entity:

        :return:
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

        return None

    @classmethod
    def to_edge(cls, entity):
        """
        Convert an entity to an edge.

        :param entity:

        :return:
        """
        if isinstance(entity, TopoDS_Edge):
            return entity

        if CheckGeom.is_curve_like(entity):
            return BRepBuilderAPI_MakeEdge(entity.handle).Edge()

        if cls.is_shape(entity) and entity.ShapeType() == TopAbs_EDGE:
            return topods_Edge(entity)

        return None

    @classmethod
    def to_wire(cls, entity):
        """
        Convert an entity to a wire.

        :param entity:

        :return:
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

        return None

    @classmethod
    def to_face(cls, entity):
        """
        Convert an entity to a face.

        :param entity:

        :return:
        """
        if isinstance(entity, TopoDS_Face):
            return entity

        if CheckGeom.is_surface_like(entity):
            return BRepBuilderAPI_MakeFace(entity.handle, 0.).Face()

        if cls.is_shape(entity) and entity.ShapeType() == TopAbs_FACE:
            return topods_Face(entity)

        return None

    @classmethod
    def to_shell(cls, entity):
        """
        Convert an entity to a shell.

        :param entity:

        :return:
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

        return None

    @classmethod
    def to_solid(cls, entity):
        """
        Convert an entity to a solid.

        :param entity:

        :return:
        """
        if isinstance(entity, TopoDS_Solid):
            return entity

        if cls.is_shape(entity) and entity.ShapeType() == TopAbs_SOLID:
            return topods_Solid(entity)

        return None

    @classmethod
    def to_compsolid(cls, entity):
        """
        Convert an entity to a compsolid.

        :param entity:

        :return:
        """
        if isinstance(entity, TopoDS_CompSolid):
            return entity

        if cls.is_shape(entity) and entity.ShapeType() == TopAbs_COMPSOLID:
            return topods_CompSolid(entity)

        return None

    @classmethod
    def to_compound(cls, entity):
        """
        Convert an entity to a compound.

        :param entity:

        :return:
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

        return None

    @classmethod
    def to_shape(cls, entity):
        """
        Convert the entity to a shape.

        :param entity:

        :return:
        """
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
                return None

        # Geometry
        if CheckGeom.is_point_like(entity):
            return cls.to_vertex(entity)
        if CheckGeom.is_curve_like(entity):
            return cls.to_edge(entity)
        if CheckGeom.is_surface_like(entity):
            return cls.to_face(entity)

        return None

    @staticmethod
    def same_parameter(edge):
        """
        Returns the SameParameter flag for the edge.

        :param edge:

        :return:
        """
        return BRep_Tool_SameParameter(edge)

    @staticmethod
    def same_range(edge):
        """
        Returns the SameRange flag for the edge.

        :param edge:

        :return:
        """
        return BRep_Tool_SameRange(edge)
