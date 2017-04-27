from .curve_part import CurvePart


class Beam(CurvePart):
    """
    Beam.
    """

    def __init__(self, label, curve_shape):
        super(Beam, self).__init__(label, curve_shape)
