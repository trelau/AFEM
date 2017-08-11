from afem.geometry.intersect import IntersectCurveCurve
from afem.structure.entities import SurfacePart
from afem.topology.bop import FuseShapes
from afem.topology.explore import ExploreShape
from afem.topology.modify import RebuildShapeByTool

__all__ = ["FuseSurfaceParts", "FuseSurfacePartsByCref"]


class FuseSurfaceParts(object):
    """
    Fuse together multiple surface parts and rebuild their shapes. The part
    shapes are rebuilt in place.

    :param list[afem.structure.entities.SurfacePart] parts: The surface parts.
    :param float fuzzy_val: Fuzzy tolerance value.
    """

    def __init__(self, parts, fuzzy_val=None):
        bop = FuseShapes(fuzzy_val=fuzzy_val)

        bop.set_args(parts)
        bop.build()

        for part in parts:
            rebuild = RebuildShapeByTool(part, bop, 'face')
            part.set_shape(rebuild.new_shape)

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
                cci = IntersectCurveCurve(main.cref, other.cref, tol)
                if not cci.success:
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
    pass
