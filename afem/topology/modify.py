from OCC.ShapeFix import ShapeFix_Shape
from OCC.ShapeUpgrade import (ShapeUpgrade_ShapeDivideClosed,
                              ShapeUpgrade_ShapeDivideContinuity,
                              ShapeUpgrade_UnifySameDomain)


class UnifyShape(object):
    """
    Unify edges and faces of a shape that lie on the same geometry.

    :param OCC.TopoDS.TopoDS_Shape shape: The shape.
    :param bool edges: Option to unify all possible edges.
    :param bool faces: Option to unify all possible faces.
    :param bool bsplines: Option to concatenate the curves of edges if they
        are C1 continuous.
    """

    def __init__(self, shape, edges=True, faces=True, bsplines=False):
        tool = ShapeUpgrade_UnifySameDomain(shape, edges, faces, bsplines)
        tool.Build()
        self._shape = tool.Shape()

    @property
    def shape(self):
        """
        :return: The unified shape.
        :rtype: OCC.TopoDS.TopoDS_Shape
        """
        return self._shape

    def generated(self, shape):
        """
        Get a shape generated from the old one.

        :param OCC.TopoDS.TopoDS_Shape shape: The old shape.

        :return: The new shape.
        :rtype: OCC.TopoDS.TopoDS_Shape
        """


class FixShape(object):
    """
    Attempt to fix the shape by applying a number of general fixes.

    :param OCC.TopoDS.TopoDS_Shape shape: The shape.
    :param float min_tol: Minimum allowed tolerance.
    :param float max_tol: Maximum allowed tolerance.
    """

    def __init__(self, shape, min_tol=None, max_tol=None):
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


class DivideClosedShape(object):
    """
    Divide all closed faces in a shape.

    :param OCC.TopoDS.TopoDS_Shape shape: The shape.
    """

    def __init__(self, shape):
        tool = ShapeUpgrade_ShapeDivideClosed(shape)
        tool.Perform()
        self._shape = tool.Result()

    @property
    def shape(self):
        """
        :return: The divided shape.
        :rtype: OCC.TopoDS.TopoDS_Shape
        """
        return self._shape


class DivideC0Shape(object):
    """
    Divide a shape at all C0 boundaries to form a C1 shape.

    :param OCC.TopoDS.TopoDS_Shape shape: The shape.
    """

    def __init__(self, shape):
        tool = ShapeUpgrade_ShapeDivideContinuity(shape)
        tool.Perform()
        self._shape = tool.Result

    @property
    def shape(self):
        """
        :return: The divided shape.
        :rtype: OCC.TopoDS.TopoDS_Shape
        """
        return self._shape
