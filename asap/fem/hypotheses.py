from OCC.NETGENPlugin import NETGENPlugin_NETGEN_2D, \
    NETGENPlugin_SimpleHypothesis_2D
from OCC.SMESH import SMESH_Gen_get
from OCC.StdMeshers import StdMeshers_Adaptive1D, StdMeshers_Deflection1D, \
    StdMeshers_LocalLength, StdMeshers_MaxLength, \
    StdMeshers_NumberOfSegments, StdMeshers_Regular_1D

_mesh_gen = SMESH_Gen_get()


class Hypothesis(object):
    """
    Base class for all hypotheses.
    """
    _hypos = {}
    _indx = 0

    def __init__(self, name, smesh_hypo):
        self._name = name
        self._hypos[name] = self
        self._smesh_hypo = smesh_hypo(Hypothesis._indx, 0, _mesh_gen)
        self._id = Hypothesis._indx
        Hypothesis._indx += 1

    @property
    def smesh_obj(self):
        return self._smesh_hypo

    @property
    def id(self):
        return self._id


class Regular1d(Hypothesis):
    """
    Regular 1-D algorithm
    """

    def __init__(self, name):
        smesh_hypo = StdMeshers_Regular_1D
        super(Regular1d, self).__init__(name, smesh_hypo)


class MaxLength1d(Hypothesis):
    """
    Maximum length 1-D hypothesis.
    """

    def __init__(self, name, max_length):
        smesh_hypo = StdMeshers_MaxLength
        super(MaxLength1d, self).__init__(name, smesh_hypo)
        self.smesh_obj.SetLength(max_length)


class LocalLength1d(Hypothesis):
    """
    Local length 1-D hypothesis.
    """

    def __init__(self, name, local_length):
        smesh_hypo = StdMeshers_LocalLength
        super(LocalLength1d, self).__init__(name, smesh_hypo)
        self.smesh_obj.SetLength(local_length)


class NumberOfSegments1d(Hypothesis):
    """
    Number of segments 1-D hypothesis.
    """

    def __init__(self, name, nseg):
        smesh_hypo = StdMeshers_NumberOfSegments
        super(NumberOfSegments1d, self).__init__(name, smesh_hypo)
        self.smesh_obj.SetNumberOfSegments(nseg)


class Adaptive1d(Hypothesis):
    """
    Adaptive length 1-D hypothesis.
    """

    def __init__(self, name, min_size, max_size, deflection):
        smesh_hypo = StdMeshers_Adaptive1D
        super(Adaptive1d, self).__init__(name, smesh_hypo)
        self.smesh_obj.SetMinSize(min_size)
        self.smesh_obj.SetMaxSize(max_size)
        self.smesh_obj.SetDeflection(deflection)


class Deflection1d(Hypothesis):
    """
    Deflection length 1-D hypothesis.
    """

    def __init__(self, name, deflection):
        smesh_hypo = StdMeshers_Deflection1D
        super(Deflection1d, self).__init__(name, smesh_hypo)
        self.smesh_obj.SetDeflection(deflection)


class NetgenSimple2d(Hypothesis):
    """
    Netgen 2-D simple hypothesis.
    """

    def __init__(self, name, local_length, allow_quads=True):
        smesh_hypo = NETGENPlugin_SimpleHypothesis_2D
        super(NetgenSimple2d, self).__init__(name, smesh_hypo)
        self.smesh_obj.SetLocalLength(local_length)
        self.smesh_obj.SetAllowQuadrangles(allow_quads)


class NetgenAlgo2d(Hypothesis):
    """
    Netgen 2-D algorithm.
    """

    def __init__(self, name):
        smesh_hypo = NETGENPlugin_NETGEN_2D
        super(NetgenAlgo2d, self).__init__(name, smesh_hypo)
