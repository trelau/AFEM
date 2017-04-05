from .surface_part import SurfacePart


class Skin(SurfacePart):
    """
    Skin part.
    """

    def __init__(self, name, surface_shape):
        super(Skin, self).__init__(name, surface_shape)
