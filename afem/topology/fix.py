from OCC.ShapeFix import ShapeFix_Shape


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
