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

# TODO Use SubMeshes instead of local hypotheses.

try:
    from OCCT.BLSURFPlugin import BLSURFPlugin_BLSURF, BLSURFPlugin_Hypothesis

    has_mg = True
except ImportError:
    BLSURFPlugin_BLSURF, BLSURFPlugin_Hypothesis = None, None
    has_mg = False

from OCCT.NETGENPlugin import (NETGENPlugin_Hypothesis, NETGENPlugin_NETGEN_2D,
                               NETGENPlugin_NETGEN_2D_ONLY,
                               NETGENPlugin_SimpleHypothesis_2D)
from OCCT.SMESH import SMESH_Gen, SMESH_Hypothesis
from OCCT.StdMeshers import (QUAD_STANDARD, StdMeshers_Adaptive1D,
                             StdMeshers_Deflection1D, StdMeshers_LocalLength,
                             StdMeshers_MaxLength,
                             StdMeshers_NumberOfSegments,
                             StdMeshers_QuadrangleParams,
                             StdMeshers_Quadrangle_2D, StdMeshers_Regular_1D)

__all__ = ["Hypothesis", "Algorithm", "Regular1D", "MaxLength1D",
           "LocalLength1D",
           "NumberOfSegments1D", "Adaptive1D", "Deflection1D",
           "NetgenHypothesis", "NetgenSimple2D", "NetgenAlgo2D",
           "NetgenAlgoOnly2D", "QuadrangleParams2D", "Quadrangle2D",
           "MeshGemsAlgo2D", "MeshGemsHypothesis2D", "HypothesisAPI"]

# Use a single instance of SMESH_Gen
the_gen = SMESH_Gen()


class Hypothesis(object):
    """
    Base class for all hypotheses.
    """
    _all = {}
    _indx = 0

    def __init__(self, label):
        self._label = label
        Hypothesis._all[label] = self
        self._id = Hypothesis._indx
        Hypothesis._indx += 1

    @property
    def label(self):
        """
        :return: The hypothesis label.
        :rtype: str
        """
        return self._label

    @property
    def id(self):
        """
        :return: The hypothesis ID.
        :rtype: int
        """
        return self._id

    @classmethod
    def get_hypothesis(cls, hypothesis):
        """
        Get a hypothesis.

        :param hypothesis: The hypothesis to get. If a hypothesis instance
            is provided it is simply returned. If a string is provided the
            label's of the hypotheses are used to find a match.
        :type hypothesis: afem.mesh.hypotheses.Hypothesis or str

        :return: The hypothesis.
        :rtype: afem.mesh.hypotheses.Hypothesis

        :raise KeyError: If hypothesis cannot be found.
        """
        if isinstance(hypothesis, Hypothesis):
            return hypothesis

        return Hypothesis._all[hypothesis]


class Algorithm(Hypothesis):
    """
    Base class for algorithms.
    """

    @property
    def object(self):
        """
        :return: The underlying hypothesis.
        :rtype: OCCT.SMESH.SMESH_Algo
        """
        return self._hypothesis

    def check_hypothesis(self, mesh, shape):
        """
        Check the hypothesis in the given mesh and shape.

        :param afem.mesh.meshes.Mesh mesh: The mesh.
        :param OCCT.TopoDS.TopoDS_Shape shape: The shape.

        :return: *True* if hypothesis is ok, *False* otherwise.
        :rtype: bool
        """
        return self.object.CheckHypothesis(mesh.object, shape,
                                           SMESH_Hypothesis.HYP_OK)

    def compute(self, mesh, shape):
        """
        Compute the mesh on a shape.

        :param afem.mesh.meshes.Mesh mesh: The mesh.
        :param OCCT.TopoDS.TopoDS_Shape shape: The shape.

        :return: *True* if completed, *False* if not.
        :rtype: bool
        """
        return self.object.Compute(mesh.object, shape)


class Regular1D(Algorithm):
    """
    Regular 1-D algorithm.
    """

    def __init__(self, label):
        self._hypothesis = StdMeshers_Regular_1D(Hypothesis._indx, 0,
                                                 the_gen)
        super(Regular1D, self).__init__(label)

    @property
    def object(self):
        """
        :return: The underlying hypothesis.
        :rtype: OCCT.StdMeshers.StdMeshers_Regular_1D
        """
        return self._hypothesis


class MaxLength1D(Hypothesis):
    """
    Maximum length 1-D hypothesis.
    """

    def __init__(self, label, max_length):
        self._hypothesis = StdMeshers_MaxLength(Hypothesis._indx, 0, the_gen)
        super(MaxLength1D, self).__init__(label)

        self.object.SetLength(max_length)

    @property
    def object(self):
        """
        :return: The underlying hypothesis.
        :rtype: OCCT.StdMeshers.StdMeshers_MaxLength
        """
        return self._hypothesis


class LocalLength1D(Hypothesis):
    """
    Local length 1-D hypothesis.
    """

    def __init__(self, label, local_length):
        self._hypothesis = StdMeshers_LocalLength(Hypothesis._indx, 0,
                                                  the_gen)
        super(LocalLength1D, self).__init__(label)

        self.object.SetLength(local_length)

    @property
    def object(self):
        """
        :return: The underlying hypothesis.
        :rtype: OCCT.StdMeshers.StdMeshers_LocalLength
        """
        return self._hypothesis


class NumberOfSegments1D(Hypothesis):
    """
    Number of segments 1-D hypothesis.
    """

    def __init__(self, label, nseg):
        self._hypothesis = StdMeshers_NumberOfSegments(Hypothesis._indx, 0,
                                                       the_gen)
        super(NumberOfSegments1D, self).__init__(label)

        self.object.SetNumberOfSegments(nseg)

    @property
    def object(self):
        """
        :return: The underlying hypothesis.
        :rtype: OCCT.StdMeshers.StdMeshers_NumberOfSegments
        """
        return self._hypothesis


class Adaptive1D(Hypothesis):
    """
    Adaptive length 1-D hypothesis.
    """

    def __init__(self, label, min_size, max_size, deflection):
        self._hypothesis = StdMeshers_Adaptive1D(Hypothesis._indx, 0,
                                                 the_gen)
        super(Adaptive1D, self).__init__(label)

        self.object.SetMinSize(min_size)
        self.object.SetMaxSize(max_size)
        self.object.SetDeflection(deflection)

    @property
    def object(self):
        """
        :return: The underlying hypothesis.
        :rtype: OCCT.StdMeshers.StdMeshers_Adaptive1D
        """
        return self._hypothesis


class Deflection1D(Hypothesis):
    """
    Deflection length 1-D hypothesis.
    """

    def __init__(self, label, deflection):
        self._hypothesis = StdMeshers_Deflection1D(Hypothesis._indx, 0,
                                                   the_gen)
        super(Deflection1D, self).__init__(label)

        self.object.SetDeflection(deflection)

    @property
    def object(self):
        """
        :return: The underlying hypothesis.
        :rtype: OCCT.StdMeshers.StdMeshers_Deflection1D
        """
        return self._hypothesis


class NetgenHypothesis(Hypothesis):
    """
    Netgen hypothesis.
    """

    def __init__(self, label, max_size=1000., min_size=0., allow_quads=False,
                 second_order=False, optimize=True, fineness=2,
                 growth_rate=0.3, nseg_per_edge=1, nseg_per_radius=2,
                 surface_curvature=False, fuse_edges=False):
        self._hypothesis = NETGENPlugin_Hypothesis(Hypothesis._indx, 0,
                                                   the_gen)
        super(NetgenHypothesis, self).__init__(label)

        self.object.SetMaxSize(max_size)
        self.object.SetMinSize(min_size)
        self.object.SetQuadAllowed(allow_quads)
        self.object.SetSecondOrder(second_order)
        self.object.SetOptimize(optimize)
        self.object.SetFineness(fineness)
        self.object.SetGrowthRate(growth_rate)
        self.object.SetNbSegPerEdge(nseg_per_edge)
        self.object.SetNbSegPerRadius(nseg_per_radius)
        self.object.SetSurfaceCurvature(surface_curvature)
        self.object.SetFuseEdges(fuse_edges)

    @property
    def object(self):
        """
        :return: The underlying hypothesis.
        :rtype: OCCT.NETGENPlugin.NETGENPlugin_Hypothesis
        """
        return self._hypothesis


class NetgenSimple2D(Hypothesis):
    """
    Netgen 2-D simple hypothesis.
    """

    def __init__(self, label, local_length, allow_quads=True,
                 length_from_edges=False, max_area=0.):
        self._hypothesis = NETGENPlugin_SimpleHypothesis_2D(Hypothesis._indx,
                                                            0, the_gen)
        super(NetgenSimple2D, self).__init__(label)

        self.object.SetLocalLength(local_length)
        self.object.SetAllowQuadrangles(allow_quads)
        if length_from_edges:
            self.object.LengthFromEdges()
        if max_area > 0.:
            self.object.SetMaxElementArea(max_area)

    @property
    def object(self):
        """
        :return: The underlying hypothesis.
        :rtype: OCCT.NETGENPlugin.NETGENPlugin_SimpleHypothesis_2D
        """
        return self._hypothesis


class NetgenAlgo2D(Algorithm):
    """
    Netgen 2-D algorithm.
    """

    def __init__(self, label):
        self._hypothesis = NETGENPlugin_NETGEN_2D(Hypothesis._indx, 0,
                                                  the_gen)
        super(NetgenAlgo2D, self).__init__(label)

    @property
    def object(self):
        """
        :return: The underlying hypothesis.
        :rtype: OCCT.NETGENPlugin.NETGENPlugin_NETGEN_2D
        """
        return self._hypothesis


class NetgenAlgoOnly2D(Algorithm):
    """
    Netgen 2-D only algorithm.
    """

    def __init__(self, label):
        self._hypothesis = NETGENPlugin_NETGEN_2D_ONLY(Hypothesis._indx, 0,
                                                       the_gen)
        super(NetgenAlgoOnly2D, self).__init__(label)

    @property
    def object(self):
        """
        :return: The underlying hypothesis.
        :rtype: OCCT.NETGENPlugin.NETGENPlugin_NETGEN_2D_ONLY
        """
        return self._hypothesis


class QuadrangleParams2D(Hypothesis):
    """
    Quadrangle 2-D parameters.
    """

    def __init__(self, label):
        self._hypothesis = StdMeshers_QuadrangleParams(Hypothesis._indx, 0,
                                                       the_gen)

        super(QuadrangleParams2D, self).__init__(label)

        self.object.SetQuadType(QUAD_STANDARD)

    @property
    def object(self):
        """
        :return: The underlying hypothesis.
        :rtype: OCCT.StdMeshers.StdMeshers_QuadrangleParams
        """
        return self._hypothesis


class Quadrangle2D(Algorithm):
    """
    Quadrangle 2-D algorithm
    """

    def __init__(self, label):
        self._hypothesis = StdMeshers_Quadrangle_2D(Hypothesis._indx, 0,
                                                    the_gen)
        super(Quadrangle2D, self).__init__(label)

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

    @property
    def object(self):
        """
        :return: The underlying hypothesis.
        :rtype: OCCT.StdMeshers.StdMeshers_Quadrangle_2D
        """
        return self._hypothesis


class MeshGemsAlgo2D(Algorithm):
    """
    MeshGems MGCAD-Surf algorithm.
    """

    def __init__(self, label):
        if not has_mg:
            raise NotImplementedError('MeshGems not available.')
        self._hypothesis = BLSURFPlugin_BLSURF(Hypothesis._indx, 0, the_gen,
                                               True)
        super(MeshGemsAlgo2D, self).__init__(label)

    @property
    def object(self):
        """
        :return: The underlying hypothesis.
        :rtype: OCCT.BLSURFPlugin.BLSURFPlugin_BLSURF
        """
        return self._hypothesis


class MeshGemsHypothesis2D(Hypothesis):
    """
    MeshGems MGCAD-Surf hypothesis.
    """

    def __init__(self, label, size=None, allow_quads=True):
        if not has_mg:
            raise NotImplementedError('MeshGems not available.')
        self._hypothesis = BLSURFPlugin_Hypothesis(Hypothesis._indx, 0, the_gen,
                                                   True)
        super(MeshGemsHypothesis2D, self).__init__(label)

        # Set a global physical size
        if size is not None:
            self.object.SetPhySize(size)
            self.object.SetMinSize(size)
            self.object.SetMaxSize(size)

        self.object.SetQuadAllowed(allow_quads)

    @property
    def object(self):
        """
        :return: The underlying hypothesis.
        :rtype: OCCT.BLSURFPlugin.BLSURFPlugin_Hypothesis
        """
        return self._hypothesis

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


class HypothesisAPI(object):
    """
    Hypothesis API. This is used to manage hypotheses from one place.
    """

    @staticmethod
    def get_hypothesis(hypothesis):
        """
        Get a hypothesis.

        :param hypothesis: The hypothesis to get. If a hypothesis instance
            is provided it is simply returned. If a string is provided the
            label's of the hypotheses are used to find a match.
        :type hypothesis: afem.mesh.hypotheses.Hypothesis or str

        :return: The hypothesis.
        :rtype: afem.mesh.hypotheses.Hypothesis

        :raise KeyError: If hypothesis cannot be found.
        """
        return Hypothesis.get_hypothesis(hypothesis)

    @staticmethod
    def create_regular_1d(label):
        """
        Create a Regular1D hypothesis.

        :param str label: The label.

        :return: The hypothesis.
        :rtype: afem.mesh.hypotheses.Regular1D
        """
        return Regular1D(label)

    @staticmethod
    def create_max_length_1d(label, max_length):
        """
        Create a MaxLength1D hypothesis.

        :param str label: The label.
        :param float max_length: The maximum length.

        :return: The hypothesis.
        :rtype: afem.mesh.hypotheses.MaxLength1D
        """
        return MaxLength1D(label, max_length)

    @staticmethod
    def create_local_length_1d(label, local_length):
        """
        Create a LocalLength1D hypothesis.

        :param str label: The label.
        :param float local_length: The local length.

        :return: The hypothesis.
        :rtype: afem.mesh.hypotheses.LocalLength1D
        """
        return LocalLength1D(label, local_length)

    @staticmethod
    def create_number_of_segments_1d(label, nseg):
        """
        Create a NumberOfSegments1D hypothesis.

        :param str label: The label.
        :param int nseg: The number of segments.

        :return: The hypothesis.
        :rtype: afem.mesh.hypotheses.NumberOfSegments1D
        """
        return NumberOfSegments1D(label, nseg)

    @staticmethod
    def create_adaptive_1d(label, min_size, max_size, deflection):
        """
        Create an Adaptive1D hypothesis.

        :param str label: The label.
        :param float min_size: The minimum size.
        :param float max_size: The maximum size.
        :param float deflection: The allowed deflection.

        :return: The hypothesis.
        :rtype: afem.mesh.hypotheses.Adaptive1D
        """
        return Adaptive1D(label, min_size, max_size, deflection)

    @staticmethod
    def create_deflection_1d(label, deflection):
        """
        Create Deflection1D hypothesis.

        :param str label: The label.
        :param float deflection: The deflection.

        :return: The hypothesis.
        :rtype: afem.mesh.hypotheses.Deflection1D
        """
        return Deflection1D(label, deflection)

    @staticmethod
    def create_netgen_hypothesis(label, max_size=1000., min_size=0.,
                                 allow_quads=False, second_order=False,
                                 optimize=True, fineness=2, growth_rate=0.3,
                                 nseg_per_edge=1, nseg_per_radius=2,
                                 surface_curvature=False, fuse_edges=False):
        """
        Create NetgenHypothesis.

        :param str label: The label.
        :param float min_size: The minimum size.
        :param float max_size: The maximum size.
        :param bool allow_quads: Option to allow quad-dominated mesh.
        :param bool second_order: Option to create second-order elements.
        :param bool optimize: Option to optimize mesh.
        :param int fineness: The fineness ratio.
        :param float growth_rate: The growth rate.
        :param int nseg_per_edge: The number of segments per edge.
        :param int nseg_per_radius: The number of segments per radius.
        :param bool surface_curvature: Option to mesh considering surface
            curvature.
        :param bool fuse_edges: Option to fuse edges with C1 continuity.

        :return: The hypothesis.
        :rtype: afem.mesh.hypotheses.NetgenHypothesis
        """
        return NetgenHypothesis(label, max_size, min_size, allow_quads,
                                second_order, optimize, fineness, growth_rate,
                                nseg_per_edge, nseg_per_radius,
                                surface_curvature, fuse_edges)

    @staticmethod
    def create_netgen_simple_2d(label, local_length, allow_quads=True,
                                length_from_edges=False, max_area=0.):
        """
        Create a NetgenSimple2D hypothesis.

        :param str label: The label.
        :param float local_length: The local length.
        :param bool allow_quads: Option to allow quad-dominated mesh.
        :param bool length_from_edges: Option to derive size of the elements
            from the element sizes on the edges.
        :param float max_area: The maximum area.

        :return: The hypothesis.
        :rtype: afem.mesh.hypotheses.NetgenSimple2D
        """
        return NetgenSimple2D(label, local_length, allow_quads,
                              length_from_edges, max_area)

    @staticmethod
    def create_netgen_algo_2d(label):
        """
        Create NetgenAlgo2D algorithm.

        :param str label: The label.

        :return: The algorithm.
        :rtype: afem.mesh.hypotheses.NetgenAlgo2D
        """
        return NetgenAlgo2D(label)

    @staticmethod
    def create_netgen_algo_only_2d(label):
        """
        Create NetgenAlgoOnly2D algorithm.

        :param str label: The label.

        :return: The algorithm.
        :rtype: afem.mesh.hypotheses.NetgenAlgoOnly2D
        """
        return NetgenAlgoOnly2D(label)

    @staticmethod
    def create_quadrangle_parameters(label):
        """
        Create QuadrangleParams2D hypothesis.

        :param str label: The label.

        :return: The hypothesis.
        :rtype: afem.mesh.hypotheses.QuadrangleParams2D
        """
        return QuadrangleParams2D(label)

    @staticmethod
    def create_quadrangle_aglo(label):
        """
        Create Quadrangle2D algorithm.

        :param str label: The label.

        :return: The hypothesis.
        :rtype: afem.mesh.hypotheses.Quadrangle2D
        """
        return Quadrangle2D(label)

    @staticmethod
    def create_meshgems_aglo_2d(label):
        """
        Create MeshGemsAlgo2D algorithm.

        :param str label: The label.

        :return: The hypothesis.
        :rtype: afem.mesh.hypotheses.MeshGemsAlgo2D
        """
        return MeshGemsAlgo2D(label)

    @staticmethod
    def create_meshgems_hypothesis_2d(label, physical_size, allow_quads=True):
        """
        Create MeshGemsHypothesis2D hypothesis.

        :param str label: The label.
        :param float physical_size:
        :param bool allow_quads:

        :return: The hypothesis.
        :rtype: afem.mesh.hypotheses.MeshGemsHypothesis2D
        """
        return MeshGemsHypothesis2D(label, physical_size, allow_quads)
