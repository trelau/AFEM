#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2017 Laughlin Research, L.L.C.
#
# This file is subject to the license agreement that was delivered
# with this source code.
#
# THE SOFTWARE AND INFORMATION ARE PROVIDED ON AN "AS IS" BASIS,
# WITHOUT ANY WARRANTIES OR REPRESENTATIONS EXPRESS, IMPLIED OR 
# STATUTORY; INCLUDING, WITHOUT LIMITATION, WARRANTIES OF QUALITY,
# PERFORMANCE, MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

try:
    from OCCT.BLSURFPlugin import BLSURFPlugin_BLSURF, BLSURFPlugin_Hypothesis

    has_mg = True
except ImportError:
    BLSURFPlugin_BLSURF, BLSURFPlugin_Hypothesis = None, None
    has_mg = False

from OCCT.NETGENPlugin import (NETGENPlugin_Hypothesis_2D,
                               NETGENPlugin_NETGEN_2D,
                               NETGENPlugin_NETGEN_2D_ONLY,
                               NETGENPlugin_SimpleHypothesis_2D)
from OCCT.SMESH import SMESH_Algo, SMESH_Hypothesis
from OCCT.StdMeshers import (StdMeshers_Adaptive1D,
                             StdMeshers_Deflection1D, StdMeshers_LocalLength,
                             StdMeshers_MaxLength,
                             StdMeshers_NumberOfSegments,
                             StdMeshers_QuadrangleParams,
                             StdMeshers_Quadrangle_2D, StdMeshers_Regular_1D,
                             StdMeshers_CompositeSegment_1D,
                             StdMeshers_QuadType)

from afem.geometry.check import CheckGeom
from afem.smesh.entities import Node

__all__ = ["Hypothesis", "Algorithm", "Regular1D", "CompositeSide1D",
           "MaxLength1D", "LocalLength1D", "NumberOfSegments1D", "Adaptive1D",
           "Deflection1D",
           "QuadrangleAlgo2D", "QuadrangleHypo2D",
           "NetgenAlgo2D", "NetgenAlgoOnly2D", "NetgenHypo2D", "NetgenSimple2D",
           "MeshGemsAlgo2D", "MeshGemsHypo2D"]


class Hypothesis(object):
    """
    Base class for all hypotheses.

    :param OCCT.SMESH.SMESH_Hypothesis hyp: The SMESH hypothesis.

    For more information see SMESH_Hypothesis_.

    .. _SMESH_Hypothesis: http://docs.salome-platform.org/latest/tui/SMESH/classSMESH__Hypothesis.html

    """

    def __init__(self, hyp):
        self._hyp = hyp

    @property
    def object(self):
        """
        :return: The underlying hypothesis.
        :rtype: OCCT.SMESH.SMESH_Hypothesis
        """
        return self._hyp

    @property
    def id(self):
        """
        :return: The hypothesis ID.
        :rtype: int
        """
        return self._hyp.GetID()

    @property
    def name(self):
        """
        :return: The hypothesis name.
        :rtype: str
        """
        return self._hyp.GetName()

    @property
    def dim(self):
        """
        :return: Hypothesis dimension.
        :rtype: int
        """
        return self._hyp.GetDim()


class Algorithm(Hypothesis):
    """
    Base class for algorithms.

    For more information see SMESH_Algo_.

    .. _SMESH_Algo: http://docs.salome-platform.org/latest/tui/SMESH/classSMESH__Algo.html

    """

    @staticmethod
    def edge_length(e):
        """
        Compute length of an edge.

        :param OCCT.TopoDS.TopoDS_Edge e: The edge.

        :return: Edge length.
        :rtype: float
        """
        return SMESH_Algo.EdgeLength_(e)

    @staticmethod
    def continuity(e1, e2):
        """
        Calculate the continuity of the two edges.

        :param OCCT.TopoDS.TopoDS_Edge e1: The first edge.
        :param OCCT.TopoDS.TopoDS_Edge e2: The second edge.

        :return: The continuity.
        :rtype: OCCT.GeomAbs.GeomAbs_Shape
        """
        return SMESH_Algo.Continuity_(e1, e2)

    @staticmethod
    def is_continuous(e1, e2):
        """
        Check if the two edges can be considered continuous.

        :param OCCT.TopoDS.TopoDS_Edge e1: The first edge.
        :param OCCT.TopoDS.TopoDS_Edge e2: The second edge.

        :return: *True* if continuous, *False* otherwise.
        :rtype: bool
        """
        return SMESH_Algo.IsContinuous_(e1, e2)

    @staticmethod
    def is_straight(e):
        """
        Check if the edge can be considered straight.

        :param OCCT.TopoDS.TopoDS_Edge e: The edge.

        :return: *True* if straight, *False* otherwise.
        :rtype: bool
        """
        return SMESH_Algo.IsStraight_(e)

    @staticmethod
    def is_degenerated(e):
        """
        Check if the edge has no 3-D curve.

        :param OCCT.TopoDS.TopoDS_Edge e: The edge.

        :return: *True* if no 3-D curve, *False* otherwise.
        :rtype: bool
        """
        return SMESH_Algo.isDegenerated_(e)

    @staticmethod
    def vertex_node(v, mesh):
        """
        Get a node built on a vertex.

        :param OCCT.TopoDS.TopoDS_Vertex v: The vertex.
        :param mesh: The mesh.
        :type mesh: afem.smesh.meshes.Mesh or afem.smesh.meshes.MeshDS

        :return: The node.
        :rtype: afem.smesh.entities.Node
        """
        return Node(SMESH_Algo.VertexNode_(v, mesh.object))

    @property
    def compatible_hypotheses(self):
        """
        :return: A list of OCCT hypotheses that are compatible with this
            algorithm.
        :rtype: list[str]
        """
        return self._hyp.GetCompatibleHypothesis()

    @property
    def compute_error(self):
        """
        :return: The compute error.
        :rtype: OCCT.SMESH.SMESH_ComputeError
        """
        return self._hyp.GetComputeError()

    @property
    def compute_error_name(self):
        """
        :return: The name of the compute error.
        :rtype: str
        """
        return self.compute_error.CommonName()

    def check_hypothesis(self, mesh, shape, status=SMESH_Hypothesis.HYP_OK):
        """
        Check the hypothesis in the given mesh and shape.

        :param afem.smesh.meshes.Mesh mesh: The mesh.
        :param OCCT.TopoDS.TopoDS_Shape shape: The shape.
        :param OCCT.SMESH.SMESH_Hypothesis.Hypothesis_Status status: The status
            to check.

        :return: *True* if check matches status, *False* otherwise.
        :rtype: bool
        """
        return self._hyp.CheckHypothesis(mesh.object, shape, status)

    def compute(self, mesh, shape):
        """
        Compute the mesh on a shape.

        :param afem.smesh.meshes.Mesh mesh: The mesh.
        :param OCCT.TopoDS.TopoDS_Shape shape: The shape.

        :return: *True* if completed, *False* if not.
        :rtype: bool
        """
        return self._hyp.Compute(mesh.object, shape)


class Regular1D(Algorithm):
    """
    Regular 1-D algorithm. Use this with a hypothesis to mesh edges.

    :param afem.smesh.meshes.MeshGen gen: A mesh generator.
    """

    def __init__(self, gen):
        hyp = StdMeshers_Regular_1D(gen.new_id(), -1, gen.object)
        super(Regular1D, self).__init__(hyp)


class CompositeSide1D(Algorithm):
    """
    Composite side 1-D discretization. This allows a 1-D hypothesis to be
    applied to a whole side of a geometrical face even if it is composed of
    several edges provided that they form C1 curve in all faces of the main
    shape.

    :param afem.smesh.meshes.MeshGen gen: A mesh generator.
    """

    def __init__(self, gen):
        hyp = StdMeshers_CompositeSegment_1D(gen.new_id(), -1, gen.object)
        super(CompositeSide1D, self).__init__(hyp)

    @staticmethod
    def get_face_side(mesh, e, f, ignore_meshed=False):
        """
        Return a face side the edge belongs to.

        :param afem.smesh.meshes.Mesh mesh: The mesh.
        :param OCCT.TopoDS.TopoDS_Edge e: The edge.
        :param OCCT.TopoDS.TopoDS_Face f: The face.
        :param bool ignore_meshed: Unclear what option does.

        :return: The face side.
        :rtype: OCCT.StdMeshers.StdMeshers_FaceSide
        """
        # TODO Wrap StdMeshers_FaceSide
        return StdMeshers_CompositeSegment_1D.GetFaceSide_(mesh.object, e, f,
                                                           ignore_meshed)


class MaxLength1D(Hypothesis):
    """
    Maximum length 1-D hypothesis.

    :param afem.smesh.meshes.MeshGen gen: A mesh generator.
    :param float max_length: Maximum edge length.
    """

    def __init__(self, gen, max_length):
        hyp = StdMeshers_MaxLength(gen.new_id(), -1, gen.object)
        super(MaxLength1D, self).__init__(hyp)

        self._hyp.SetLength(max_length)

    @property
    def max_length(self):
        """
        :return: The maximum edge length.
        :rtype: float
        """
        return self._hyp.GetLength()


class LocalLength1D(Hypothesis):
    """
    Local length 1-D hypothesis.

    :param afem.smesh.meshes.MeshGen gen: A mesh generator.
    :param float local_length: The desired edge length.
    :param float precision: The value used to control the rounding of the
        number of segments. Use 0.5 to provide rounding to nearest integer,
        1.0 for lower integer, 0.0 for higher integer.
    """

    def __init__(self, gen, local_length, precision=1.0e-7):
        hyp = StdMeshers_LocalLength(gen.new_id(), -1, gen.object)
        super(LocalLength1D, self).__init__(hyp)

        self._hyp.SetLength(local_length)
        self._hyp.SetPrecision(precision)

    @property
    def local_length(self):
        """
        :return: The local edge length.
        :rtype: float
        """
        return self._hyp.GetLength()

    @property
    def precision(self):
        """
        :return: The precision.
        :rtype: float
        """
        return self._hyp.GetPrecision()


class NumberOfSegments1D(Hypothesis):
    """
    Number of segments 1-D hypothesis.

    :param afem.smesh.meshes.MeshGen gen: A mesh generator.
    :param int nseg: The number of segments.
    """

    def __init__(self, gen, nseg):
        hyp = StdMeshers_NumberOfSegments(gen.new_id(), -1, gen.object)
        super(NumberOfSegments1D, self).__init__(hyp)

        self._hyp.SetNumberOfSegments(nseg)

    @property
    def nseg(self):
        """
        :return: The number of segments.
        :rtype: int
        """
        return self._hyp.GetNumberOfSegments()


class Adaptive1D(Hypothesis):
    """
    Adaptive length 1-D hypothesis.

    :param afem.smesh.meshes.MeshGen gen: A mesh generator.
    :param float min_size: Minimum edge size.
    :param float max_size: Maximum edge size.
    :param float deflection: Maximum distance from a segment to a curved edge.
    """

    def __init__(self, gen, min_size, max_size, deflection):
        hyp = StdMeshers_Adaptive1D(gen.new_id(), -1, gen.object)
        super(Adaptive1D, self).__init__(hyp)

        self._hyp.SetMinSize(min_size)
        self._hyp.SetMaxSize(max_size)
        self._hyp.SetDeflection(deflection)


class Deflection1D(Hypothesis):
    """
    Deflection length 1-D hypothesis.

    :param afem.smesh.meshes.MeshGen gen: A mesh generator.
    :param float deflection: Maximum distance from a segment to a curved edge.
    """

    def __init__(self, gen, deflection):
        hyp = StdMeshers_Deflection1D(gen.new_id(), -1, gen.object)
        super(Deflection1D, self).__init__(hyp)

        self._hyp.SetDeflection(deflection)


class QuadrangleAlgo2D(Algorithm):
    """
    Quadrangle 2-D algorithm.

    :param afem.smesh.meshes.MeshGen gen: A mesh generator.
    """

    def __init__(self, gen):
        hyp = StdMeshers_Quadrangle_2D(gen.new_id(), -1, gen.object)
        super(QuadrangleAlgo2D, self).__init__(hyp)

    @staticmethod
    def is_applicable(shape, check_all=True):
        """
        Check if this algorithm can mesh the shape.

        :param shape: The shape.
        :param bool check_all: If *True*, this check returns *True* if all
            shapes are applicable. If *False* this check returns *True* if at
            least one shape is ok.

        :return: Check whether algorithm is applicable.
        :rtype: bool
        """
        return StdMeshers_Quadrangle_2D.IsApplicable_(shape, check_all)


class QuadrangleHypo2D(Hypothesis):
    """
    Quadrangle 2-D parameters.

    :param afem.smesh.meshes.MeshGen gen: A mesh generator.
    :param OCCT.StdMeshers.StdMeshers_QuadType quad_type: The quadrangle
        preference for transitions.
    """

    def __init__(self, gen, quad_type=StdMeshers_QuadType.QUAD_STANDARD):
        hyp = StdMeshers_QuadrangleParams(gen.new_id(), -1, gen.object)
        super(QuadrangleHypo2D, self).__init__(hyp)

        self._hyp.SetQuadType(quad_type)

    def set_enforced_nodes(self, shapes, pnts):
        """
        Set enforced nodes on shapes.

        :param list[OCCT.TopoDS.TopoDS_Shape] shapes: List of shapes.
        :param list[point_like] pnts: List of points for enforced nodes.

        :return: None.
        """
        pnts = [CheckGeom.to_point(p) for p in pnts]
        self._hyp.SetEnforcedNodes(shapes, pnts)


class NetgenAlgo2D(Algorithm):
    """
    NETGEN 2-D algorithm.

    :param afem.smesh.meshes.MeshGen gen: A mesh generator.
    """

    def __init__(self, gen):
        hyp = NETGENPlugin_NETGEN_2D(gen.new_id(), -1, gen.object)
        super(NetgenAlgo2D, self).__init__(hyp)


class NetgenAlgoOnly2D(Algorithm):
    """
    NETGEN 2-D only algorithm. Takes into account pre-existing nodes on face
    boundaries.

    :param afem.smesh.meshes.MeshGen gen: A mesh generator.
    """

    def __init__(self, gen):
        hyp = NETGENPlugin_NETGEN_2D_ONLY(gen.new_id(), -1, gen.object)
        super(NetgenAlgoOnly2D, self).__init__(hyp)


class NetgenHypo2D(Hypothesis):
    """
    NETGEN 2-D hypothesis.

    :param afem.smesh.meshes.MeshGen gen: A mesh generator.
    :param float max_size: Maximum edge length.
    :param float min_size: Minimum edge length.
    :param bool allow_quads: Enable quad-dominated mesh.
    :param bool second_order: Enable second-order mesh.
    :param bool optimize: Enable mesh optimization.
    :param int fineness: Mesh fineness (0=very coarse, 1=coarse, 2=moderate,
        3=fine, 4=very fine, 5=user defined).
    :param float growth_rate: Growth rate between 0 to 1.
    :param int nseg_per_edge: Number of segments per edge.
    :param int nseg_per_radius: Number of segments per radius.
    :param bool surface_curvature: Enable surface curvature to define edge size.
    :param bool fuse_edges: Option to fuse edges.
    """

    def __init__(self, gen, max_size=1000., min_size=0.,
                 allow_quads=False, second_order=False, optimize=True,
                 fineness=2, growth_rate=0.3, nseg_per_edge=1,
                 nseg_per_radius=2, surface_curvature=False, fuse_edges=False):
        hyp = NETGENPlugin_Hypothesis_2D(gen.new_id(), -1, gen.object)
        super(NetgenHypo2D, self).__init__(hyp)

        self._hyp.SetMaxSize(max_size)
        self._hyp.SetMinSize(min_size)
        self._hyp.SetSecondOrder(second_order)
        self._hyp.SetOptimize(optimize)
        self._hyp.SetFineness(fineness)
        self._hyp.SetGrowthRate(growth_rate)
        self._hyp.SetNbSegPerEdge(nseg_per_edge)
        self._hyp.SetNbSegPerRadius(nseg_per_radius)
        self._hyp.SetQuadAllowed(allow_quads)
        self._hyp.SetSurfaceCurvature(surface_curvature)
        self._hyp.SetFuseEdges(fuse_edges)


class NetgenSimple2D(Hypothesis):
    """
    NETGEN 2-D simple hypothesis.

    :param afem.smesh.meshes.MeshGen gen: A mesh generator.
    :param float local_length: Local edge segment length.
    :param int nseg: Number of edge segments.
    :param bool allow_quads: Enable quad-dominated mesh.
    :param bool length_from_edges: Set maximum element area to be dependent on
        1-D discretization.
    :param float max_area: Set maximum element area.
    """

    def __init__(self, gen, local_length=None, nseg=None,
                 allow_quads=True, length_from_edges=False, max_area=0.):
        hyp = NETGENPlugin_SimpleHypothesis_2D(gen.new_id(), -1, gen.object)
        super(NetgenSimple2D, self).__init__(hyp)

        if local_length is not None:
            self._hyp.SetLocalLength(local_length)
        if nseg is not None:
            self._hyp.SetNumberOfSegments(nseg)
        self._hyp.SetAllowQuadrangles(allow_quads)
        if length_from_edges:
            self._hyp.LengthFromEdges()
        if max_area > 0.:
            self._hyp.SetMaxElementArea(max_area)


class MeshGemsAlgo2D(Algorithm):
    """
    MeshGems MGCAD-Surf algorithm.

    :param afem.smesh.meshes.MeshGen gen: A mesh generator.

    :raise NotImplementedError: If MeshGems is not available.
    """

    def __init__(self, gen):
        if not has_mg:
            raise NotImplementedError('MeshGems not available.')
        hyp = BLSURFPlugin_BLSURF(gen.new_id(), -1, gen.object, True)
        super(MeshGemsAlgo2D, self).__init__(hyp)


class MeshGemsHypo2D(Hypothesis):
    """
    MeshGems MGCAD-Surf hypothesis.

    :param afem.smesh.meshes.MeshGen gen: A mesh generator.
    :param float size: Desired global edge size.
    :param bool allow_quads: Enable quad-dominant mesh.

    :raise NotImplementedError: If MeshGems is not available.
    """

    def __init__(self, gen, size=None, allow_quads=True):
        if not has_mg:
            raise NotImplementedError('MeshGems not available.')
        hyp = BLSURFPlugin_Hypothesis(gen.new_id(), -1, gen.object, True)
        super(MeshGemsHypo2D, self).__init__(hyp)

        # Set a global physical size
        if size is not None:
            self._hyp.SetPhySize(size)
            self._hyp.SetMinSize(size)
            self._hyp.SetMaxSize(size)

        self._hyp.SetQuadAllowed(allow_quads)

    def set_physical_size(self, size, is_rel=False):
        """
        Set physical size.

        :param float size:
        :param bool is_rel:
        :return:
        """
        self._hyp.SetPhySize(size, is_rel)

    def set_min_size(self, size, is_rel=False):
        """
        Set minimum size.

        :param float size:
        :param bool is_rel:
        :return:
        """
        self._hyp.SetMinSize(size, is_rel)

    def set_max_size(self, size, is_rel=False):
        """
        Set maximum size.

        :param float size:
        :param bool is_rel:
        :return:
        """
        self._hyp.SetMaxSize(size, is_rel)

    def set_use_gradation(self, val=True):
        """

        :param bool val:
        :return:
        """
        self._hyp.SetUseGradation(val)

    def set_gradation(self, val):
        """

        :param float val:
        :return:
        """
        self._hyp.SetGradation(val)

    def set_quads_allowed(self, val):
        """

        :param bool val:
        :return:
        """
        self._hyp.SetQuadAllowed(val)

    def set_angle_mesh(self, val):
        """

        :param float val:
        :return:
        """
        self._hyp.SetAngleMesh(val)

    def set_chordal_error(self, val):
        """

        :param float val:
        :return:
        """
        self._hyp.SetChordalError(val)

    def set_anisotropic(self, val):
        """

        :param bool val:
        :return:
        """
        self._hyp.SetAnisotropic(val)

    def set_anisotropic_ratio(self, val):
        """

        :param float val:
        :return:
        """
        self._hyp.SetAnisotropicRatio(val)

    def set_remove_tiny_edges(self, val):
        """

        :param bool val:
        :return:
        """
        self._hyp.SetRemoveTinyEdges(val)

    def set_tiny_edge_length(self, val):
        """

        :param float val:
        :return:
        """
        self._hyp.SetTinyEdgeLength(val)

    def set_optimize_tiny_edges(self, val):
        """

        :param bool val:
        :return:
        """
        self._hyp.SetOptimiseTinyEdges(val)

    def set_tiny_edge_optimization_length(self, val):
        """

        :param float val:
        :return:
        """
        self._hyp.SetTinyEdgeOptimisationLength(val)

    def set_correct_surface_intersection(self, val):
        """

        :param bool val:
        :return:
        """
        self._hyp.SetCorrectSurfaceIntersection(val)

    def set_correct_surface_intersection_max_cost(self, val):
        """

        :param float val:
        :return:
        """
        self._hyp.SetCorrectSurfaceIntersectionMaxCost(val)

    def set_bad_element_removal(self, val):
        """

        :param bool val:
        :return:
        """
        self._hyp.SetBadElementRemoval(val)

    def set_bad_element_aspect_ratio(self, val):
        """

        :param float val:
        :return:
        """
        self._hyp.SetBadElementAspectRatio(val)

    def set_optimize_mesh(self, val):
        """

        :param bool val:
        :return:
        """
        self._hyp.SetOptimizeMesh(val)

    def set_respect_geometry(self, val):
        """

        :param bool val:
        :return:
        """
        self._hyp.SetRespectGeometry(val)

    def set_max_number_of_threads(self, val):
        """

        :param int val:
        :return:
        """
        self._hyp.SetMaxNumberOfThreads(val)

    def set_debug(self, val):
        """

        :param bool val:
        :return:
        """
        self._hyp.SetDebug(val)

    def set_required_entities(self, val):
        """

        :param str val: "respect", "ignore", "clear"
        :return:
        """
        self._hyp.SetRequiredEntities(val)

    def set_sewing_tolerance(self, val):
        """

        :param float val:
        :return:
        """
        self._hyp.SetSewingTolerance(val)

    def set_tags(self, val):
        """

        :param str val: "respect", "ignore", "clear"
        :return:
        """
        self._hyp.SetTags(val)

    def add_option(self, name, val):
        """

        :param str name:
        :param str val:
        :return:
        """
        self._hyp.AddOption(name, val)

    def add_hyperpatch(self, patch_ids):
        """

        :param list[set(int)] patch_ids:
        :return:
        """
        patch = [set(patch_ids)]
        self._hyp.SetHyperPatches(patch)
