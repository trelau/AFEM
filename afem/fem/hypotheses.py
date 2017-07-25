from OCC.NETGENPlugin import NETGENPlugin_Hypothesis, NETGENPlugin_NETGEN_2D, \
    NETGENPlugin_NETGEN_2D_ONLY, NETGENPlugin_SimpleHypothesis_2D
from OCC.SMESH import SMESH_Gen_get
from OCC.StdMeshers import StdMeshers_Adaptive1D, StdMeshers_Deflection1D, \
    StdMeshers_LocalLength, StdMeshers_MaxLength, \
    StdMeshers_NumberOfSegments, StdMeshers_Regular_1D

_mesh_gen = SMESH_Gen_get()

__all__ = ["Regular1D", "MaxLength1D", "LocalLength1D", "NumberOfSegments1D",
           "Adaptive1D", "Deflection1D", "NetgenHypothesis", "NetgenSimple2D",
           "NetgenAlgo2D", "NetgenAlgoOnly2D", "HypothesisData"]


class Hypothesis(object):
    """
    Base class for all hypotheses.
    """
    _all = {}
    _indx = 0

    def __init__(self, name, smesh_hypo):
        self._name = name
        Hypothesis._all[name] = self
        self._smesh_hypo = smesh_hypo(Hypothesis._indx, 0, _mesh_gen)
        self._id = Hypothesis._indx
        Hypothesis._indx += 1

    @property
    def smesh_obj(self):
        return self._smesh_hypo

    @property
    def id(self):
        return self._id

    @classmethod
    def get_hypothesis(cls, hypothesis):
        """
        Get a hypothesis.

        :param hypothesis:

        :return:
        """
        if isinstance(hypothesis, Hypothesis):
            return hypothesis
        try:
            return Hypothesis._all[hypothesis]
        except KeyError:
            return None


class Regular1D(Hypothesis):
    """
    Regular 1-D algorithm
    """

    def __init__(self, name):
        smesh_hypo = StdMeshers_Regular_1D
        super(Regular1D, self).__init__(name, smesh_hypo)


class MaxLength1D(Hypothesis):
    """
    Maximum length 1-D hypothesis.
    """

    def __init__(self, name, max_length):
        smesh_hypo = StdMeshers_MaxLength
        super(MaxLength1D, self).__init__(name, smesh_hypo)
        self.smesh_obj.SetLength(max_length)


class LocalLength1D(Hypothesis):
    """
    Local length 1-D hypothesis.
    """

    def __init__(self, name, local_length):
        smesh_hypo = StdMeshers_LocalLength
        super(LocalLength1D, self).__init__(name, smesh_hypo)
        self.smesh_obj.SetLength(local_length)


class NumberOfSegments1D(Hypothesis):
    """
    Number of segments 1-D hypothesis.
    """

    def __init__(self, name, nseg):
        smesh_hypo = StdMeshers_NumberOfSegments
        super(NumberOfSegments1D, self).__init__(name, smesh_hypo)
        self.smesh_obj.SetNumberOfSegments(nseg)


class Adaptive1D(Hypothesis):
    """
    Adaptive length 1-D hypothesis.
    """

    def __init__(self, name, min_size, max_size, deflection):
        smesh_hypo = StdMeshers_Adaptive1D
        super(Adaptive1D, self).__init__(name, smesh_hypo)
        self.smesh_obj.SetMinSize(min_size)
        self.smesh_obj.SetMaxSize(max_size)
        self.smesh_obj.SetDeflection(deflection)


class Deflection1D(Hypothesis):
    """
    Deflection length 1-D hypothesis.
    """

    def __init__(self, name, deflection):
        smesh_hypo = StdMeshers_Deflection1D
        super(Deflection1D, self).__init__(name, smesh_hypo)
        self.smesh_obj.SetDeflection(deflection)


class NetgenHypothesis(Hypothesis):
    """
    Netgen hypothesis.
    """

    def __init__(self, name, max_size=1000., min_size=0., allow_quads=False,
                 second_order=False, optimize=True, fineness=2,
                 growth_rate=0.3, nseg_per_edge=1, nseg_per_radius=2,
                 surface_curvature=False, fuse_edges=False):
        smesh_hypo = NETGENPlugin_Hypothesis
        super(NetgenHypothesis, self).__init__(name, smesh_hypo)
        self.smesh_obj.SetMaxSize(max_size)
        self.smesh_obj.SetMinSize(min_size)
        self.smesh_obj.SetQuadAllowed(allow_quads)
        self.smesh_obj.SetSecondOrder(second_order)
        self.smesh_obj.SetOptimize(optimize)
        self.smesh_obj.SetFineness(fineness)
        self.smesh_obj.SetGrowthRate(growth_rate)
        self.smesh_obj.SetNbSegPerEdge(nseg_per_edge)
        self.smesh_obj.SetNbSegPerRadius(nseg_per_radius)
        self.smesh_obj.SetSurfaceCurvature(surface_curvature)
        self.smesh_obj.SetFuseEdges(fuse_edges)


class NetgenSimple2D(Hypothesis):
    """
    Netgen 2-D simple hypothesis.
    """

    def __init__(self, name, local_length, allow_quads=True,
                 length_from_edges=False, max_area=0.):
        smesh_hypo = NETGENPlugin_SimpleHypothesis_2D
        super(NetgenSimple2D, self).__init__(name, smesh_hypo)
        self.smesh_obj.SetLocalLength(local_length)
        self.smesh_obj.SetAllowQuadrangles(allow_quads)
        if length_from_edges:
            self.smesh_obj.LengthFromEdges()
        if max_area > 0.:
            self.smesh_obj.SetMaxElementArea(max_area)


class NetgenAlgo2D(Hypothesis):
    """
    Netgen 2-D algorithm.
    """

    def __init__(self, name):
        smesh_hypo = NETGENPlugin_NETGEN_2D
        super(NetgenAlgo2D, self).__init__(name, smesh_hypo)


class NetgenAlgoOnly2D(Hypothesis):
    """
    Netgen 2-D only algorithm.
    """

    def __init__(self, name):
        smesh_hypo = NETGENPlugin_NETGEN_2D_ONLY
        super(NetgenAlgoOnly2D, self).__init__(name, smesh_hypo)


class HypothesisData(object):
    """
    Hypothesis data manager.
    """

    @staticmethod
    def get_hypothesis(hypothesis):
        """
        Get a hypothesis.

        :param hypothesis:

        :return:
        """
        return Hypothesis.get_hypothesis(hypothesis)

    @staticmethod
    def create_regular_1d(name):
        """
        Create a Regular1D hypothesis.

        :param str name:

        :return:
        """
        return Regular1D(name)

    @staticmethod
    def create_max_length_1d(name, max_length):
        """
        Create a MaxLength1D hypothesis.

        :param name:
        :param max_length:

        :return:
        """
        return MaxLength1D(name, max_length)

    @staticmethod
    def create_local_length_1d(name, local_length):
        """
        Create a LocalLength1D hypothesis.

        :param name:
        :param local_length:

        :return:
        """
        return LocalLength1D(name, local_length)

    @staticmethod
    def create_number_of_segments_1d(name, nseg):
        """
        Create a NumberOfSegments1D hypothesis.

        :param name:
        :param nseg:

        :return:
        """
        return NumberOfSegments1D(name, nseg)

    @staticmethod
    def create_adaptive_1d(name, min_size, max_size, deflection):
        """
        Create an Adaptive1D hypothesis.

        :param name:
        :param min_size:
        :param max_size:
        :param deflection:

        :return:
        """
        return Adaptive1D(name, min_size, max_size, deflection)

    @staticmethod
    def create_deflection_1d(name, deflection):
        """
        Create Deflection1D hypothesis.

        :param name:
        :param deflection:

        :return:
        """
        return Deflection1D(name, deflection)

    @staticmethod
    def create_netgen_hypothesis(name, max_size=1000., min_size=0.,
                                 allow_quads=False, second_order=False,
                                 optimize=True, fineness=2, growth_rate=0.3,
                                 nseg_per_edge=1, nseg_per_radius=2,
                                 surface_curvature=False, fuse_edges=False):
        """
        Create NetgenHypothesis.

        :param name:
        :param max_size:
        :param min_size:
        :param allow_quads:
        :param second_order:
        :param optimize:
        :param fineness:
        :param growth_rate:
        :param nseg_per_edge:
        :param nseg_per_radius:
        :param surface_curvature:
        :param fuse_edges:

        :return:
        """
        return NetgenHypothesis(name, max_size, min_size, allow_quads,
                                second_order, optimize, fineness, growth_rate,
                                nseg_per_edge, nseg_per_radius,
                                surface_curvature, fuse_edges)

    @staticmethod
    def create_netgen_simple_2d(name, local_length, allow_quads=True,
                                length_from_edges=False, max_area=0.):
        """
        Create a NetgenSimple2D hypothesis.

        :param name:
        :param local_length:
        :param allow_quads:
        :param length_from_edges:
        :param max_area:

        :return:
        """
        return NetgenSimple2D(name, local_length, allow_quads,
                              length_from_edges, max_area)

    @staticmethod
    def create_netgen_algo_2d(name):
        """
        Create NetgenAlgo2D hypothesis.

        :param name:

        :return:
        """
        return NetgenAlgo2D(name)

    @staticmethod
    def create_netgen_algo_only_2d(name):
        """
        Create NetgenAlgoOnly2D hypothesis.

        :param name:

        :return:
        """
        return NetgenAlgoOnly2D(name)
