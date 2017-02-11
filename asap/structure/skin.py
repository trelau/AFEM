from .surface_part import SurfacePart


class Skin(SurfacePart):
    """
    Skin part.
    """

    def __init__(self, name, rshape):
        super(Skin, self).__init__(name, rshape)
