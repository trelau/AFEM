from .curve_part import CurvePart


class Stiffener(CurvePart):
    """
    Stiffener for surface parts.
    """

    def __init__(self, name, curve_shape):
        super(Stiffener, self).__init__(name, curve_shape, False)
