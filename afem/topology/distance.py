from OCC.BRepExtrema import BRepExtrema_DistShapeShape

__all__ = ["DistanceShapeToShape"]


class DistanceShapeToShape(object):
    """
    Calculate minimum distance between two shapes.

    :param OCC.TopoDS.TopoDS_Shape shape1: The first shape.
    :param OCC.TopoDS.TopoDS_Shape shape2: The second shape.

    :raise RuntimeError: If OCC method fails.

    Usage:

    >>> from afem.topology import *
    >>> v1 = VertexByPoint((0., 0., 0.)).vertex
    >>> v2 = VertexByPoint((10., 0., 0.)).vertex
    >>> tool = DistanceShapeToShape(v1, v2)
    >>> tool.nsol
    1
    >>> tool.dmin
    10.0
    """

    def __init__(self, shape1, shape2):
        self._tool = BRepExtrema_DistShapeShape(shape1, shape2)
        if not self._tool.IsDone():
            msg = 'OCC BRepExtrema_DistShapeShape failed.'
            raise RuntimeError(msg)

    @property
    def nsol(self):
        """
        :return: The number of solutions satisfying the minimum distance.
        :rtype:
        """
        return self._tool.NbSolution()

    @property
    def dmin(self):
        """
        :return: The minimum distance.
        :rtype: float
        """
        return self._tool.Value()
