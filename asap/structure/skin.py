from .surface_part import SurfacePart


class Skin(SurfacePart):
    """
    Skin part.
    """

    def __init__(self, name, shape, cref=None, sref=None, assy=None):
        super(Skin, self).__init__(name, shape, cref, sref, assy)
