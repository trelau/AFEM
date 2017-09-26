from numpy import mean

from afem.structure.entities import SurfacePart
from afem.topology.bop import (FuseShapes, IntersectShapes,
                               SplitShapes)
from afem.topology.check import CheckShape
from afem.topology.create import CompoundByShapes, EdgeByCurve
from afem.topology.explore import ExploreShape
from afem.topology.modify import RebuildShapesByTool, SewShape

__all__ = ["FuseSurfaceParts", "FuseSurfacePartsByCref", "CutParts",
           "SewSurfaceParts", "SplitParts", "FuseAssemblies"]


class FuseSurfaceParts(object):
    """
    Fuse together multiple surface parts and rebuild their shapes.

    :param parts: The surface parts.
    :type parts: collections.Sequence[afem.structure.entities.SurfacePart]
    :param parts: The other surface parts.
    :type tools: collections.Sequence[afem.structure.entities.SurfacePart]
    :param float fuzzy_val: Fuzzy tolerance value.
    """

    def __init__(self, parts, tools, fuzzy_val=None):
        bop = FuseShapes(fuzzy_val=fuzzy_val)

        parts = list(parts)
        tools = list(tools)
        bop.set_args(parts)
        bop.set_tools(tools)
        bop.build()

        rebuild = RebuildShapesByTool(parts + tools, bop)
        for part in parts + tools:
            new_shape = rebuild.new_shape(part)
            part.set_shape(new_shape)

        self._is_done = bop.is_done
        self._fused_shape = bop.shape

    @property
    def is_done(self):
        """
        :return: *True* if operation is done, *False* if not.
        :rtype: bool
        """
        return self._is_done

    @property
    def fused_shape(self):
        """
        :return: The fused shape.
        :rtype: OCC.TopoDS.TopoDS_Shape
        """
        return self._fused_shape


class FuseSurfacePartsByCref(object):
    """
    Attempt to automatically fuse together adjacent parts based on the
    possible intersection of their reference curve. The part shapes are
    rebuilt in place.

    :param list[afem.structure.entities.SurfacePart] parts: The surface parts.
    :param float tol: The tolerance to use for checking possible
        intersections of the reference curves. Default is the maximum
        tolerance of the part shape.

    :raises TypeError: If a given part is not a surface part.
    """

    def __init__(self, parts, tol=None):
        self._is_done = False

        for part in parts:
            if not isinstance(part, SurfacePart):
                msg = 'Part is not a surface part.'
                raise TypeError(msg)

        # Test all combinations of parts for intersection of reference curve
        join_parts = []
        main_parts = []
        nparts = len(parts)
        for i in range(0, nparts - 1):
            main = parts[i]
            other_parts = []
            for j in range(i + 1, nparts):
                other = parts[j]
                if not main.has_cref or not other.has_cref:
                    continue
                if tol is None:
                    tol1 = ExploreShape.get_tolerance(main, 1)
                    tol2 = ExploreShape.get_tolerance(other, 1)
                    tol = max(tol1, tol2)
                e1 = EdgeByCurve(main.cref).edge
                e2 = EdgeByCurve(other.cref).edge
                bop = IntersectShapes(e1, e2, fuzzy_val=tol)
                if not bop.vertices:
                    continue
                # Store potential join
                other_parts.append(other)
            if other_parts:
                main_parts.append(main)
                join_parts.append(other_parts)

        # Join the parts
        for main, other_parts in zip(main_parts, join_parts):
            main.fuse(*other_parts)
            self._is_done = True

    @property
    def is_done(self):
        """
        :return: *True* if at least one fuse was performed, *False* if not.
        :rtype: bool
        """
        return self._is_done


class CutParts(object):
    """
    Cut each part with a shape and rebuild the part shape.

    :param parts: The parts to cut.
    :type parts: collections.Sequence[afem.structure.entities.Part]
    :param shape: The shape to cut with.
    :type shape: OCC.TopoDS.TopoDS_Shape or afem.geometry.entities.Surface
    """

    def __init__(self, parts, shape):
        parts = list(parts)

        shape = CheckShape.to_shape(shape)

        # TODO Really just loop through each one?
        for part in parts:
            part.cut(shape)


class SewSurfaceParts(object):
    """
    Sew edges of the surface parts and rebuild their shapes.

    :param parts: The parts to sew.
    :type parts: collections.Sequence[afem.structure.entities.SurfacePart]
    :param float tol: The tolerance. If not provided then the average
        tolerance from all part shapes will be used.
    :param float max_tol: Maximum tolerance. If not provided then the maximum
        tolerance from all part shapes will be used.
    """

    def __init__(self, parts, tol=None, max_tol=None):
        parts = list(parts)

        if tol is None:
            tol = mean([ExploreShape.get_tolerance(part, 0) for part in parts],
                       dtype=float)

        if max_tol is None:
            max_tol = max([ExploreShape.get_tolerance(part, 1) for part
                           in parts])

        sew = SewShape(tol=tol, max_tol=max_tol, cut_free_edges=True,
                       non_manifold=True)

        for part in parts:
            sew.add(part)
        sew.perform()

        for part in parts:
            if not sew.is_modified(part):
                continue
            mod_shape = sew.modified(part)
            part.set_shape(mod_shape)


class SplitParts(object):
    """
    Split part shapes and rebuild in place.

    :param list[afem.structure.entities.Part] parts: The parts that will be
        split and rebuilt.
    :param tools: The parts or shapes used to split the parts but are not
        modified.
    :type tools: list[OCC.TopoDS.TopoDS_Shape or afem.structure.entities.Part]
    :param float fuzzy_val: Fuzzy tolerance value.
    """

    def __init__(self, parts, tools=None, fuzzy_val=None):
        bop = SplitShapes(fuzzy_val=fuzzy_val)

        bop.set_args(parts)

        if tools is not None:
            bop.set_tools(tools)

        bop.build()

        rebuild = RebuildShapesByTool(parts, bop)
        for part in parts:
            new_shape = rebuild.new_shape(part)
            part.set_shape(new_shape)

        self._is_done = bop.is_done
        self._split_shape = bop.shape

    @property
    def is_done(self):
        """
        :return: *True* if operation is done, *False* if not.
        :rtype: bool
        """
        return self._is_done

    @property
    def split_shape(self):
        """
        :return: The split shape.
        :rtype: OCC.TopoDS.TopoDS_Shape
        """
        return self._split_shape


class FuseAssemblies(object):
    """
    Fuse assemblies and rebuild the part shapes. This tool puts all the part
    shapes into compounds before the Boolean operation.

    :param afem.structure.assembly.Assembly assy1: The first assembly.
    :param afem.structure.assembly.Assembly assy2: The second assembly.
    :param float fuzzy_val: Fuzzy tolerance value.
    :param bool include_subassy: Option to recursively include parts
            from any sub-assemblies.
    :param afem.structure.assembly.Assembly other_assys: Extra assemblies.
    """

    def __init__(self, assy1, assy2, fuzzy_val=None, include_subassy=True,
                 *other_assys):
        bop = FuseShapes(fuzzy_val=fuzzy_val)

        parts1 = assy1.get_parts(include_subassy)
        shape1 = CompoundByShapes(parts1).compound
        bop.set_args([shape1])

        tools = []
        other_parts = []
        for assy in [assy2] + list(other_assys):
            parts = assy.get_parts(include_subassy)
            other_parts += parts
            shape = CompoundByShapes(parts).compound
            tools.append(shape)
        bop.set_tools(tools)

        bop.build()

        all_parts = parts1 + other_parts
        rebuild = RebuildShapesByTool(all_parts, bop)
        for part in all_parts:
            new_shape = rebuild.new_shape(part)
            part.set_shape(new_shape)

        self._bop = bop

    @property
    def is_done(self):
        """
        :return: *True* if operation is done, *False* if not.
        :rtype: bool
        """
        return self._bop.is_done

    @property
    def fused_shape(self):
        """
        :return: The fused shape.
        :rtype: OCC.TopoDS.TopoDS_Shape
        """
        return self._bop.shape
