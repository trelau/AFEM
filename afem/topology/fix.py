from OCC.ShapeFix import ShapeFix_Shape, ShapeFix_ShapeTolerance
from OCC.TopAbs import TopAbs_SHAPE

__all__ = ["FixTolerance", "FixShape"]


class FixTolerance(object):
    """
    Methods for adjusting shape tolerances.
    """

    @staticmethod
    def limit_tolerance(shape, tmin, tmax=0., styp=TopAbs_SHAPE):
        """
        Limit tolerances in a shape.

        :param OCC.TopoDS.TopoDS_Shape shape: The shape.
        :param float tmin: The minimum tolerance.
        :param float tmax: The maximum tolerance.
        :param OCC.TopAbs.TopAbs_ShapeEnum styp: The level of shape to set
            (i.e., only vertices, only edges, only faces, or all shapes).

        :return: *True* if at least one tolerance of a sub-shape was modified.
        :rtype: bool
        """
        return ShapeFix_ShapeTolerance().LimitTolerance(shape, tmin, tmax,
                                                        styp)

    @staticmethod
    def set_tolerance(shape, tol, styp=TopAbs_SHAPE):
        """
        Enforce tolerance on the given shape.

        :param OCC.TopoDS.TopoDS_Shape shape: The shape.
        :param float tol: The tolerance.
        :param OCC.TopAbs.TopAbs_ShapeEnum styp: The level of shape to set
            (i.e., only vertices, only edges, only faces, or all shapes).

        :return: None.
        """
        return ShapeFix_ShapeTolerance().SetTolerance(shape, tol, styp)


class FixShape(object):
    """
    Attempt to fix the shape by applying a number of general fixes.

    :param OCC.TopoDS.TopoDS_Shape shape: The shape.
    :param float min_tol: Minimum allowed tolerance.
    :param float max_tol: Maximum allowed tolerance.
    """

    def __init__(self, shape, min_tol=None, max_tol=None):
        # TODO Use context to fix sub-shapes
        # TODO Option to set specific tools
        fix = ShapeFix_Shape(shape)
        if min_tol is not None:
            fix.SetMinTolerance(min_tol)
        if max_tol is not None:
            fix.SetMaxTolerance(max_tol)
        fix.Perform()
        self._shape = fix.Shape()

    @property
    def shape(self):
        """
        :return: The fixed shape.
        :rtype: OCC.TopoDS.TopoDS_Shape
        """
        return self._shape
