from math import sqrt

from OCC.BRepOffset import BRepOffset_Skin
from OCC.BRepOffsetAPI import (BRepOffsetAPI_MakeOffsetShape,
                               BRepOffsetAPI_MakePipe,
                               BRepOffsetAPI_MakePipeShell,
                               BRepOffsetAPI_NormalProjection,
                               BRepOffsetAPI_ThruSections)
from OCC.GeomAbs import GeomAbs_C2
from OCC.TopAbs import TopAbs_EDGE, TopAbs_VERTEX, TopAbs_WIRE

from afem.occ.occ_utils import (occ_continuity, occ_join_type, occ_parm_type,
                                occ_transition_mode)
from afem.topology.shape_check import CheckShape
from afem.topology.shape_explore import ExploreShape

__all__ = ["ProjectShape", "OffsetShape", "LoftShape", "SweepShape",
           "SweepShapeWithNormal"]


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
    """
    Offset a shape.

    :param OCC.TopoDS.TopoDS_Shape shape: The shape. It may be a face,
        shell, a solid, or a compound of these kinds.
    :param float offset: The offset value. The offset will be outside the
        shape if positive and inside if negative.
    :param float tol: Tolerance for coincidence for generated shapes. If not
        provided the average tolerance of the shape is used.
    :param str join_mode: Option for how to fill holes that may appear when
        offsetting two adjacent faces ('arc' or 'intersect').

    For more information see BRepOffsetAPI_MakeOffsetShape_.

    .. _BRepOffsetAPI_MakeOffsetShape: https://www.opencascade.com/doc/occt-7.1.0/refman/html/class_b_rep_offset_a_p_i___make_offset_shape.html

    Usage:

    >>> from afem.geometry import *
    >>> from afem.topology import *
    >>> pln = PlaneByAxes().plane
    >>> face = FaceByPlane(pln, -5., 5., -5., 5.).face
    >>> tool = OffsetShape(face, 5.)
    >>> shape = tool.shape
    """

    def __init__(self, shape, offset, tol=None, join_mode='arc'):
        if tol is None:
            tol = ExploreShape.get_tolerance(shape)

        join_mode = occ_join_type[join_mode.lower()]

        # TODO Remove internal edges option for OCC 7
        self._tool = BRepOffsetAPI_MakeOffsetShape(shape, offset, tol,
                                                   BRepOffset_Skin, False,
                                                   False, join_mode)
        self._tool.Build()
        self._shape = self._tool.Shape()

    @property
    def shape(self):
        """
        :return: The offset shape.
        :rtype: OCC.TopoDS.TopoDS_Shape
        """
        return self._shape


class LoftShape(object):
    """
    Loft a shape using a sequence of sections.

    :param sections: The sections of the loft. These
        are usually wires but the first and last section can be vertices.
        Edges are converted to wires before adding to the loft tool.
    :type sections: list[OCC.TopoDS.TopoDS_Vertex or OCC.TopoDS.TopoDS_Edge or
        OCC.TopoDS.TopoDS_Wire]
    :param bool is_solid: If *True* the tool will build a solid, otherwise
        it will build a shell.
    :param bool make_ruled: If *True* the faces between sections will be ruled
        surfaces, otherwise they are smoothed out by approximation.
    :param float pres3d: Defines the precision for the approximation algorithm.
    :param bool check_compatibility: Option to check the orientation of the
        sections to avoid twisted results and update to have the same number
        of edges.
    :param bool use_smoothing: Option to use approximation algorithm.
    :param str par_type: Parametrization type ('chord', 'uniform',
        'centripetal').
    :param str continuity: The desired continuity ('C0', 'G1', 'C1', 'G2',
        'C2', 'C3').
    :param int max_degree: The maximum degree for the approximation
        algorithm.

    :raise TypeError: If any of the sections cannot be added to the tool
        because they are of the wrong type.

    For more information see BRepOffsetAPI_ThruSections_.

    .. _BRepOffsetAPI_ThruSections: https://www.opencascade.com/doc/occt-7.1.0/refman/html/class_b_rep_offset_a_p_i___thru_sections.html

    Usage:

    >>> from afem.topology import *
    >>> pnts1 = [(0., 0., 0.), (5., 0., 5.), (10., 0., 0.)]
    >>> wire1 = WireByPoints(pnts1).wire
    >>> pnts2 = [(0., 10., 0.), (5., 10., -5.), (10., 10., 0.)]
    >>> wire2 = WireByPoints(pnts2).wire
    >>> loft = LoftShape([wire1, wire2])
    >>> loft.is_done
    True
    """

    def __init__(self, sections, is_solid=False, make_ruled=False,
                 pres3d=1.0e-6, check_compatibility=None,
                 use_smoothing=None, par_type=None, continuity=None,
                 max_degree=None):
        self._tool = BRepOffsetAPI_ThruSections(is_solid, make_ruled, pres3d)

        if check_compatibility is not None:
            self._tool.CheckCompatibility(check_compatibility)

        if use_smoothing is not None:
            self._tool.SetSmoothing(use_smoothing)

        if par_type is not None:
            par_type = occ_parm_type[par_type.lower()]
            self._tool.SetParType(par_type)

        if continuity is not None:
            continuity = occ_continuity[continuity.upper()]
            self._tool.SetContinuity(continuity)

        if max_degree is not None:
            self._tool.SetMaxDegree(max_degree)

        for section in sections:
            if section.ShapeType() == TopAbs_VERTEX:
                self._tool.AddVertex(section)
            elif section.ShapeType() == TopAbs_EDGE:
                wire = CheckShape.to_wire(section)
                self._tool.AddWire(wire)
            elif section.ShapeType() == TopAbs_WIRE:
                self._tool.AddWire(section)
            else:
                msg = 'Invalid shape type in loft.'
                raise TypeError(msg)

        self._tool.Build()

    @property
    def is_done(self):
        """
        :return: *True* if done, *False* if not.
        :rtype: bool
        """
        return self._tool.IsDone()

    @property
    def shape(self):
        """
        :return: The lofted shape.
        :rtype: OCC.TopoDS.TopoDS_Shape
        """
        return self._tool.Shape()

    @property
    def first_shape(self):
        """
        :return: The first/bottom shape of the loft if a solid was
            constructed.
        :rtype: OCC.TopoDS.TopoDS_Shape
        """
        return self._tool.FirstShape()

    @property
    def last_shape(self):
        """
        :return: The last/top shape of the loft if a solid was constructed.
        :rtype: OCC.TopoDS.TopoDS_Shape
        """
        return self._tool.LastShape()

    @property
    def max_degree(self):
        """
        :return: The max degree used in the approximation algorithm
        :rtype: int
        """
        return self._tool.MaxDegree()

    def generated_face(self, edge):
        """
        Get a face(s) generated by the edge. If the ruled option was used,
        then this returns each face generated by the edge. If the smoothing
        option was used, then this returns the face generated by the edge.

        :param OCC.TopoDS.TopoDS_Edge edge: The edge.

        :return: The face(s) generated by the edge.
        :rtype: OCC.TopoDS.TopoDS_Shape
        """
        return self._tool.GeneratedFace(edge)


class SweepShape(object):
    """
    Sweep a profile along a spine.

    :param OCC.TopoDS.TopoDS_Wire spine: The path for the sweep. This must
        be at least G1 continuous.
    :param OCC.TopoDS.TopoDS_Shape profile: The profile.

    For more information see BRepOffsetAPI_MakePipe_.

    .. _BRepOffsetAPI_MakePipe: https://www.opencascade.com/doc/occt-7.1.0/refman/html/class_b_rep_offset_a_p_i___make_pipe.html
    """

    def __init__(self, spine, profile):
        self._tool = BRepOffsetAPI_MakePipe(spine, profile)

        self._tool.Build()

    @property
    def is_done(self):
        """
        :return: *True* if done, *False* if not.
        :rtype: bool
        """
        return self._tool.IsDone()

    @property
    def shape(self):
        """
        :return: The swept shape.
        :rtype: OCC.TopoDS.TopoDS_Shape
        """
        return self._tool.Shape()

    @property
    def first_shape(self):
        """
        :return: The first/bottom shape of the sweep.
        :rtype: OCC.TopoDS.TopoDS_Shape
        """
        return self._tool.FirstShape()

    @property
    def last_shape(self):
        """
        :return: The last/top shape of the sweep.
        :rtype: OCC.TopoDS.TopoDS_Shape
        """
        return self._tool.LastShape()


class SweepShapeWithNormal(object):
    """
    Sweep sections along a spine using a support shape to define the
    local orientation.

    :param OCC.TopoDS.TopoDS_Wire spine: The spine.
    :param OCC.TopoDS.TopoDS_Shape spine_support: The shape that will define
        the normal during the sweeping algorithm. To be effective, each edge of
        the spine must have a representation on one face of the spine support.
    :param float tol3d: The 3-D tolerance.
    :param float tol_bound: The boundary tolerance.
    :param float tol_angular: The angular tolerance.
    :param int max_degree: The maximum degree allowed in the resulting surface.
    :param int max_segments: The maximum number of segments allowed in the
        resulting surface.
    :param bool force_c1: If *True*, the tool will attempt to approximate a
        C1 surface if a swept surface proved to be C0.
    :param str transition_mode: The transition mode to manage
        discontinuities on the swept shape ('transformed', 'right', 'round').

    For more information see BRepOffsetAPI_MakePipeShell_.

    .. _BRepOffsetAPI_MakePipeShell: https://www.opencascade.com/doc/occt-7.1.0/refman/html/class_b_rep_offset_a_p_i___make_pipe_shell.html
    """

    def __init__(self, spine, spine_support, tol3d=1.0e-4, tol_bound=1.0e-4,
                 tol_angular=1.0e-2, max_degree=None, max_segments=None,
                 force_c1=None, transition_mode='transformed', ):
        self._tool = BRepOffsetAPI_MakePipeShell(spine)
        self._tool.SetMode(spine_support)

        self._tool.SetTolerance(tol3d, tol_bound, tol_angular)

        if max_degree is not None:
            self._tool.SetMaxDegree(max_degree)

        if max_segments is not None:
            self._tool.SetMaxSegments(max_segments)

        if force_c1 is not None:
            self._tool.SetForceApproxC1(force_c1)

        if transition_mode is not None:
            transition_mode = occ_transition_mode[transition_mode.lower()]
            self._tool.SetTransitionMode(transition_mode)

    def add_profile(self, profile, with_contact=False, with_correction=False):
        """
        Add the profile to the tool.

        :param profile: The profile to add.
        :type profile: OCC.TopoDS.TopoDS_Vertex or OCC.TopoDS.TopoDS_Edge or
            OCC.TopoDS.TopoDS_Wire
        :param bool with_contact: If *True*, then the profile is translated
            to be in contact with the spine.
        :param bool with_correction: If *True*, then the profile is rotated
            to be orthogonal to the spine's tangent.

        :return: None.

        :raise TypeError: If the profile is not a vertex, edge, or wire.
        """
        if profile.ShapeType() == TopAbs_EDGE:
            profile = CheckShape.to_wire(profile)

        if profile.ShapeType() not in [TopAbs_VERTEX, TopAbs_WIRE]:
            msg = 'Invalid profile type.'
            raise TypeError(msg)

        self._tool.Add(profile, with_contact, with_correction)

    @property
    def is_ready(self):
        """
        :return: *True* if tool is ready to build the shape. *False* if not.
        :rtype: bool
        """
        return self._tool.IsReady()

    def build(self):
        """
        Build the resulting shape.

        :return: None.
        """
        self._tool.Build()

    def make_solid(self):
        """
        Attempts to transform the sweeping shell into a solid.

        :return: *True* if done, *False* if not.
        :rtype: bool
        """
        return self._tool.MakeSolid()

    @property
    def is_done(self):
        """
        :return: *True* if done, *False* if not.
        :rtype: bool
        """
        return self._tool.IsDone()

    @property
    def shape(self):
        """
        :return: The resulting shape.
        :rtype: OCC.TopoDS.TopoDS_Shape
        """
        return self._tool.Shape()

    @property
    def first_shape(self):
        """
        :return: The first/bottom shape of the sweep.
        :rtype: OCC.TopoDS.TopoDS_Shape
        """
        return self._tool.FirstShape()

    @property
    def last_shape(self):
        """
        :return: The last/top shape of the sweep.
        :rtype: OCC.TopoDS.TopoDS_Shape
        """
        return self._tool.LastShape()
