from OCC.BRepAdaptor import BRepAdaptor_Surface

from .tools import ShapeTools
from ..geometry.methods.create import create_nurbs_surface_from_occ


class ShapeAdaptorSurface(object):
    """
    Shape adaptor surface.
    """

    def __init__(self):
        self._shape = None
        self._surface = None

    @property
    def surface(self):
        return self._surface

    @property
    def shape(self):
        return self._shape

    def build(self, face, split_closed=True, split_c0=True):
        """
        Build the adaptor surface.

        :param face:
        :param split_closed:
        :param split_c0:

        :return:
        """
        face = ShapeTools.to_face(face)
        if not face:
            return False

        # Split face.
        if split_closed:
            shape = ShapeTools.divide_closed(face)
        else:
            shape = face
        if split_c0:
            shape = ShapeTools.divide_c0(shape)
        if not shape or shape.IsNull():
            return False

        # Build AFEM surface.
        adp_srf = BRepAdaptor_Surface(face)
        occ_srf = adp_srf.BSpline().GetObject()
        srf = create_nurbs_surface_from_occ(occ_srf)
        if not srf:
            return False

        self._shape = shape
        self._surface = srf
        return True
