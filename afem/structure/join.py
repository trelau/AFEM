from afem.structure.entities import SurfacePart
from afem.topology.bop import FuseShapes, SplitShapes, IntersectShapes
from afem.topology.explore import ExploreShape
from afem.topology.modify import RebuildShapesByTool
from afem.topology.create import EdgeByCurve

__all__ = ["FuseSurfaceParts", "FuseSurfacePartsByCref", "SplitParts"]


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
    pass


class SewParts(object):
    pass


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
