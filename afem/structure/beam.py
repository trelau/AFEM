from .curve_part import CurvePart

__all__ = ["Beam"]


class Beam(CurvePart):
    """
    Beam.
    """

    def __init__(self, label, shape, cref=None, assy=None):
        super(Beam, self).__init__(label, shape, cref, assy)
