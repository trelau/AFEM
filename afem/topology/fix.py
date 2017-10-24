from OCC.ShapeBuild import ShapeBuild_ReShape
from OCC.ShapeFix import ShapeFix_Shape, ShapeFix_ShapeTolerance
from OCC.TopAbs import TopAbs_SHAPE

__all__ = ["FixShape"]

_fix_tol = ShapeFix_ShapeTolerance()


class FixShape(object):
    """
    Attempt to fix the shape by applying a number of general fixes.

    :param OCC.TopoDS.TopoDS_Shape shape: The shape.
    :param float precision: Basic precision value.
    :param float min_tol: Minimum allowed tolerance.
    :param float max_tol: Maximum allowed tolerance.
    :param OCC.TopoDS.TopoDS_Shape context: The context shape.
    """

    def __init__(self, shape, precision=None, min_tol=None, max_tol=None,
                 context=None):
        self._tool = ShapeFix_Shape()

        if precision is not None:
            self._tool.SetPrecision(precision)
        if min_tol is not None:
            self._tool.SetMinTolerance(min_tol)
        if max_tol is not None:
            self._tool.SetMaxTolerance(max_tol)

        if context is not None:
            reshape = ShapeBuild_ReShape()
            reshape.Apply(context)
            self._tool.SetContext(reshape.GetHandle())

        self._tool.Init(shape)
        self._tool.Perform()

    @property
    def shape(self):
        """
        :return: The fixed shape.
        :rtype: OCC.TopoDS.TopoDS_Shape
        """
        return self._tool.Shape()

    @property
    def context(self):
        """
        :return: The context.
        :rtype: OCC.ShapeBuild.ShapeBuild_ReShape
        """
        return self._tool.Context().GetObject()

    def apply(self, shape):
        """
        Apply substitutions to the shape (or subshape) and get the result.

        :param OCC.TopoDS.TopoDS_Shape shape: The shape.

        :return: The new shape.
        :rtype: OCC.TopoDS.TopoDS_Shape
        """
        return self.context.Apply(shape)

    @staticmethod
    def limit_tolerance(shape, tol=1.0e-7, styp=TopAbs_SHAPE):
        """
        Limit tolerances in a shape.

        :param OCC.TopoDS.TopoDS_Shape shape: The shape.
        :param float tol: Target tolerance.
        :param OCC.TopAbs.TopAbs_ShapeEnum styp: The level of shape to set
            (i.e., only vertices, only edges, only faces, or all shapes).

        :return: *True* if at least one tolerance of a sub-shape was modified.
        :rtype: bool
        """
        return _fix_tol.LimitTolerance(shape, tol, 0., styp)

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
        return _fix_tol.SetTolerance(shape, tol, styp)
