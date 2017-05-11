from .curve_part import CurvePart


class Beam(CurvePart):
    """
    Beam.
    """

    def __init__(self, label, shape, cref=None):
        super(Beam, self).__init__(label, shape, cref)
