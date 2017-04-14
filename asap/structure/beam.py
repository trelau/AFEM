from .curve_part import CurvePart


class Beam(CurvePart):
    """
    Beam.
    """

    def __init__(self, name, curve_shape):
        super(Beam, self).__init__(name, curve_shape)
