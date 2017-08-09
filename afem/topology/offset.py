from math import sqrt

from OCC.BRepOffsetAPI import BRepOffsetAPI_NormalProjection
from OCC.GeomAbs import GeomAbs_C2

from afem.geometry.create import occ_continuity
from afem.topology.explore import ExploreShape

__all__ = ["ProjectShape"]


class ProjectShape(object):
    """
    Project edges and wires onto a basis shape.

    :param OCC.TopoDS.TopoDS_Shape: The shape to project to.
    :param to_project: List of edges or wires to project.
    :type to_project: list[OCC.TopoDS.TopoDS_Edge or OCC.TopoDS.TopoDS_Wire]
    :param float tol3d: The 3-D tolerance.
    :param float tol2d: The 2-D tolerance. If not provided then
        *sqrt(tol3d)* is used.
    :param str continuity: Desired continuity ('C0', 'G1', 'C1', 'G2', 'C2',
        'C3').
    :param int max_degree: Max degree.
    :param int max_seg: Max segments.
    :param float max_dist: Max distance between target shape and shapes to
        project. If not satisfied then results for the corresponding shape
        are discarded.
    :param bool limit: Option to limit projected edges to the face boundaries.

    For more information see BRepOffsetAPI_NormalProjection_.

    .. _BRepOffsetAPI_NormalProjection: https://www.opencascade.com/doc/occt-7.1.0/refman/html/class_b_rep_offset_a_p_i___normal_projection.html

    Usage:

    >>> from afem.geometry import *
    >>> from afem.topology import *
    >>> pln = PlaneByAxes().plane
    >>> face = FaceByPlane(pln, -5., 5., -5., 5.).face
    >>> edge = EdgeByPoints((0., 1., 15.), (0., 1., -15.)).edge
    >>> proj = ProjectShape(face, [edge])
    >>> proj.is_done
    True
    >>> proj.nedges
    1
    """

    def __init__(self, shape, to_project, tol3d=1.0e-4, tol2d=None,
                 continuity='C2', max_degree=14, max_seg=16, max_dist=None,
                 limit=True):
        tool = BRepOffsetAPI_NormalProjection(shape)

        if tol2d is None:
            tol2d = sqrt(tol3d)
        try:
            cont = occ_continuity[continuity.upper()]
        except (KeyError, AttributeError):
            cont = GeomAbs_C2
        tool.SetParams(tol3d, tol2d, cont, max_degree, max_seg)

        if max_dist is not None:
            tool.SetMaxDistance(max_dist)

        if limit:
            tool.SetLimit(True)

        for item in to_project:
            tool.Add(item)

        tool.Build()
        self._tool = tool

    @property
    def is_done(self):
        """
        :return: *True* if projection was correctly built, *False* if not.
        :rtype: bool
        """
        return self._tool.IsDone()

    @property
    def projection(self):
        """
        :return: The projected shape. Tries to build the result as a
            compound of wires.
        :rtype: OCC.TopoDS.TopoDS_Shape
        """
        return self._tool.Projection()

    @property
    def nwires(self):
        """
        :return: The number of wires in the projection.
        :rtype: int
        """
        return len(self.wires)

    @property
    def wires(self):
        """
        :return: A list of wires in the projected shape.
        :rtype: list[OCC.TopoDS.TopoDS_Wire]
        """
        return ExploreShape.get_wires(self.projection)

    @property
    def nedges(self):
        """
        :return: The number of edges in the projection.
        :rtype: int
        """
        return len(self.edges)

    @property
    def edges(self):
        """
        :return: A list of edges in the projected shape.
        :rtype: list[OCC.TopoDS.TopoDS_Edge]
        """
        return ExploreShape.get_edges(self.projection)


class OffsetShape(object):
    # TODO OffsetShape
    pass


class LoftShape(object):
    # TODO LoftShape
    pass


class SweepShape(object):
    # TODO SweepShape
    pass


class SweepWithNormal(object):
    # TODO SweepWithNormal
    pass
