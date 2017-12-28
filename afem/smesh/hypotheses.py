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
from OCCT.SMESH import SMESH_Hypothesis
from OCCT.StdMeshers import (QUAD_STANDARD, StdMeshers_Adaptive1D,
                             StdMeshers_Deflection1D, StdMeshers_LocalLength,
                             StdMeshers_MaxLength,
                             StdMeshers_NumberOfSegments,
                             StdMeshers_QuadrangleParams,
                             StdMeshers_Quadrangle_2D, StdMeshers_Regular_1D)

__all__ = ["Hypothesis", "Algorithm", "Regular1D", "MaxLength1D",
           "LocalLength1D", "NumberOfSegments1D", "Adaptive1D", "Deflection1D",
           "NetgenAlgo2D", "NetgenAlgoOnly2D", "NetgenHypo2D",
           "NetgenSimple2D", "QuadrangleAlgo2D", "QuadrangleHypo2D",
           "MeshGemsAlgo2D", "MeshGemsHypo2D"]


class Hypothesis(object):
    """
    Base class for all hypotheses.

    :param OCCT.SMESH.SMESH_Hypothesis hyp: The SMESH hypothesis.
    """

    def __init__(self, hyp):
        self._hyp = hyp

    @property
    def object(self):
        """
        :return: The underlying hypothesis.
        """
        return self._hyp

    @property
    def id(self):
        """
        :return: The hypothesis ID.
        :rtype: int
        """
        return self._hyp.GetID()


class Algorithm(Hypothesis):
    """
    Base class for algorithms.
    """

    def check_hypothesis(self, mesh, shape):
        """
        Check the hypothesis in the given mesh and shape.

        :param afem.smesh.meshes.Mesh mesh: The mesh.
        :param OCCT.TopoDS.TopoDS_Shape shape: The shape.

        :return: *True* if hypothesis is ok, *False* otherwise.
        :rtype: bool
        """
        return self.object.CheckHypothesis(mesh.object, shape,
                                           SMESH_Hypothesis.HYP_OK)

    def compute(self, mesh, shape):
        """
        Compute the mesh on a shape.

        :param afem.smesh.meshes.Mesh mesh: The mesh.
        :param OCCT.TopoDS.TopoDS_Shape shape: The shape.

        :return: *True* if completed, *False* if not.
        :rtype: bool
        """
        return self.object.Compute(mesh.object, shape)


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
    """
    # TODO Implement CompositeSide1D
    pass


class MaxLength1D(Hypothesis):
    """
    Maximum length 1-D hypothesis.
    """

    def __init__(self, gen, max_length):
        hyp = StdMeshers_MaxLength(gen.new_id(), -1, gen.object)
        super(MaxLength1D, self).__init__(hyp)

        self.object.SetLength(max_length)


class LocalLength1D(Hypothesis):
    """
    Local length 1-D hypothesis.
    """

    def __init__(self, gen, local_length):
        hyp = StdMeshers_LocalLength(gen.new_id(), -1, gen.object)
        super(LocalLength1D, self).__init__(hyp)

        self.object.SetLength(local_length)


class NumberOfSegments1D(Hypothesis):
    """
    Number of segments 1-D hypothesis.
    """

    def __init__(self, gen, nseg):
        hyp = StdMeshers_NumberOfSegments(gen.new_id(), -1, gen.object)
        super(NumberOfSegments1D, self).__init__(hyp)

        self.object.SetNumberOfSegments(nseg)


class Adaptive1D(Hypothesis):
    """
    Adaptive length 1-D hypothesis.
    """

    def __init__(self, gen, min_size, max_size, deflection):
        hyp = StdMeshers_Adaptive1D(gen.new_id(), -1, gen.object)
        super(Adaptive1D, self).__init__(hyp)

        self.object.SetMinSize(min_size)
        self.object.SetMaxSize(max_size)
        self.object.SetDeflection(deflection)


class Deflection1D(Hypothesis):
    """
    Deflection length 1-D hypothesis.
    """

    def __init__(self, gen, deflection):
        hyp = StdMeshers_Deflection1D(gen.new_id(), -1, gen.object)
        super(Deflection1D, self).__init__(hyp)

        self.object.SetDeflection(deflection)


class NetgenAlgo2D(Algorithm):
    """
    Netgen 2-D algorithm.
    """

    def __init__(self, gen):
        hyp = NETGENPlugin_NETGEN_2D(gen.new_id(), -1, gen.object)
        super(NetgenAlgo2D, self).__init__(hyp)


class NetgenAlgoOnly2D(Algorithm):
    """
    Netgen 2-D only algorithm.
    """

    def __init__(self, gen):
        hyp = NETGENPlugin_NETGEN_2D_ONLY(gen.new_id(), -1, gen.object)
        super(NetgenAlgoOnly2D, self).__init__(hyp)


class NetgenHypo2D(Hypothesis):
    """
    NETGEN 2-D hypothesis.
    """

    def __init__(self, gen, max_size=1000., min_size=0.,
                 allow_quads=False, second_order=False, optimize=True,
                 fineness=2, growth_rate=0.3, nseg_per_edge=1,
                 nseg_per_radius=2, surface_curvature=False, fuse_edges=False):
        hyp = NETGENPlugin_Hypothesis_2D(gen.new_id(), -1, gen.object)
        super(NetgenHypo2D, self).__init__(hyp)

        self.object.SetMaxSize(max_size)
        self.object.SetMinSize(min_size)
        self.object.SetSecondOrder(second_order)
        self.object.SetOptimize(optimize)
        self.object.SetFineness(fineness)
        self.object.SetGrowthRate(growth_rate)
        self.object.SetNbSegPerEdge(nseg_per_edge)
        self.object.SetNbSegPerRadius(nseg_per_radius)
        self.object.SetQuadAllowed(allow_quads)
        self.object.SetSurfaceCurvature(surface_curvature)
        self.object.SetFuseEdges(fuse_edges)


class NetgenSimple2D(Hypothesis):
    """
    NETGEN 2-D simple hypothesis.
    """

    def __init__(self, gen, local_length, allow_quads=True,
                 length_from_edges=False, max_area=0.):
        hyp = NETGENPlugin_SimpleHypothesis_2D(gen.new_id(), -1, gen.object)
        super(NetgenSimple2D, self).__init__(hyp)

        self.object.SetLocalLength(local_length)
        self.object.SetAllowQuadrangles(allow_quads)
        if length_from_edges:
            self.object.LengthFromEdges()
        if max_area > 0.:
            self.object.SetMaxElementArea(max_area)


class QuadrangleAlgo2D(Algorithm):
    """
    Quadrangle 2-D algorithm
    """

    def __init__(self, gen):
        hyp = StdMeshers_Quadrangle_2D(gen.new_id(), -1, gen.object)
        super(QuadrangleAlgo2D, self).__init__(hyp)

    def is_applicable(self, shape, check_all=True):
        """
        Check if this algorithm can mesh the shape.


        :param shape: The shape.
        :param bool check_all: If *True*, this check returns *True* if all
            shapes are applicable. If *False* this check returns *True* if at
            least one shape is ok.

        :return: Check whether algorithm is applicable.
        :rtype: bool
        """
        return self.object.IsApplicable_(shape, check_all)


class QuadrangleHypo2D(Hypothesis):
    """
    Quadrangle 2-D parameters.
    """

    def __init__(self, gen):
        hyp = StdMeshers_QuadrangleParams(gen.new_id(), -1, gen.object)
        super(QuadrangleHypo2D, self).__init__(hyp)

        self.object.SetQuadType(QUAD_STANDARD)


class MeshGemsAlgo2D(Algorithm):
    """
    MeshGems MGCAD-Surf algorithm.
    """

    def __init__(self, gen):
        if not has_mg:
            raise NotImplementedError('MeshGems not available.')
        hyp = BLSURFPlugin_BLSURF(gen.new_id(), -1, gen.object, True)
        super(MeshGemsAlgo2D, self).__init__(hyp)


class MeshGemsHypo2D(Hypothesis):
    """
    MeshGems MGCAD-Surf hypothesis.
    """

    def __init__(self, gen, size=None, allow_quads=True):
        if not has_mg:
            raise NotImplementedError('MeshGems not available.')
        hyp = BLSURFPlugin_Hypothesis(gen.new_id(), -1, gen.object, True)
        super(MeshGemsHypo2D, self).__init__(hyp)

        # Set a global physical size
        if size is not None:
            self.object.SetPhySize(size)
            self.object.SetMinSize(size)
            self.object.SetMaxSize(size)

        self.object.SetQuadAllowed(allow_quads)

    def set_physical_size(self, size, is_rel=False):
        """
        Set physical size.

        :param float size:
        :param bool is_rel:
        :return:
        """
        self.object.SetPhySize(size, is_rel)

    def set_min_size(self, size, is_rel=False):
        """
        Set minimum size.

        :param float size:
        :param bool is_rel:
        :return:
        """
        self.object.SetMinSize(size, is_rel)

    def set_max_size(self, size, is_rel=False):
        """
        Set maximum size.

        :param float size:
        :param bool is_rel:
        :return:
        """
        self.object.SetMaxSize(size, is_rel)

    def set_use_gradation(self, val=True):
        """

        :param bool val:
        :return:
        """
        self.object.SetUseGradation(val)

    def set_gradation(self, val):
        """

        :param float val:
        :return:
        """
        self.object.SetGradation(val)

    def set_quads_allowed(self, val):
        """

        :param bool val:
        :return:
        """
        self.object.SetQuadAllowed(val)

    def set_angle_mesh(self, val):
        """

        :param float val:
        :return:
        """
        self.object.SetAngleMesh(val)

    def set_chordal_error(self, val):
        """

        :param float val:
        :return:
        """
        self.object.SetChordalError(val)

    def set_anisotropic(self, val):
        """

        :param bool val:
        :return:
        """
        self.object.SetAnisotropic(val)

    def set_anisotropic_ratio(self, val):
        """

        :param float val:
        :return:
        """
        self.object.SetAnisotropicRatio(val)

    def set_remove_tiny_edges(self, val):
        """

        :param bool val:
        :return:
        """
        self.object.SetRemoveTinyEdges(val)

    def set_tiny_edge_length(self, val):
        """

        :param float val:
        :return:
        """
        self.object.SetTinyEdgeLength(val)

    def set_optimize_tiny_edges(self, val):
        """

        :param bool val:
        :return:
        """
        self.object.SetOptimiseTinyEdges(val)

    def set_tiny_edge_optimization_length(self, val):
        """

        :param float val:
        :return:
        """
        self.object.SetTinyEdgeOptimisationLength(val)

    def set_correct_surface_intersection(self, val):
        """

        :param bool val:
        :return:
        """
        self.object.SetCorrectSurfaceIntersection(val)

    def set_correct_surface_intersection_max_cost(self, val):
        """

        :param float val:
        :return:
        """
        self.object.SetCorrectSurfaceIntersectionMaxCost(val)

    def set_bad_element_removal(self, val):
        """

        :param bool val:
        :return:
        """
        self.object.SetBadElementRemoval(val)

    def set_bad_element_aspect_ratio(self, val):
        """

        :param float val:
        :return:
        """
        self.object.SetBadElementAspectRatio(val)

    def set_optimize_mesh(self, val):
        """

        :param bool val:
        :return:
        """
        self.object.SetOptimizeMesh(val)

    def set_respect_geometry(self, val):
        """

        :param bool val:
        :return:
        """
        self.object.SetRespectGeometry(val)

    def set_max_number_of_threads(self, val):
        """

        :param int val:
        :return:
        """
        self.object.SetMaxNumberOfThreads(val)

    def set_debug(self, val):
        """

        :param bool val:
        :return:
        """
        self.object.SetDebug(val)

    def set_required_entities(self, val):
        """

        :param str val: "respect", "ignore", "clear"
        :return:
        """
        self.object.SetRequiredEntities(val)

    def set_sewing_tolerance(self, val):
        """

        :param float val:
        :return:
        """
        self.object.SetSewingTolerance(val)

    def set_tags(self, val):
        """

        :param str val: "respect", "ignore", "clear"
        :return:
        """
        self.object.SetTags(val)

    def add_option(self, name, val):
        """

        :param str name:
        :param str val:
        :return:
        """
        self.object.AddOption(name, val)

    def add_hyperpatch(self, patch_ids):
        """

        :param list[set(int)] patch_ids:
        :return:
        """
        patch = [set(patch_ids)]
        self.object.SetHyperPatches(patch)
