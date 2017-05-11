from .methods.cut_parts import cut_wing_part_with_circle
from .methods.modify_parts import discard_wing_part_faces
from .surface_part import SurfacePart
from ..geometry import CreateGeom


class WingPart(SurfacePart):
    """
    Base class for wing parts.
    """

    def __init__(self, label, shape, cref=None, sref=None):
        super(WingPart, self).__init__(label, shape, cref, sref)

    def discard(self, shape=None, tol=None):
        """
        Discard faces of the part that are inside the solid or use automated
        if *solid* is *None*.

        :param shape:
        :param tol:

        :return:
        """
        if not shape:
            # Automatic method for wing parts.
            return discard_wing_part_faces(self)
        else:
            # Call original method.
            return super(WingPart, self).discard(shape, tol)

    def cut_hole(self, dx, r):
        """
        Cut a circular hole in the part (in development).

        :param dx:
        :param r:

        :return:
        """
        return cut_wing_part_with_circle(self, dx, r)

    def get_plane(self, u=None, dx=None, pnt=None, ref_pln=None):
        """
        Get a plane along the reference curve.

        :param u:
        :param dx:
        :param pnt:
        :param ref_pln:

        :return:
        """
        return CreateGeom.plane_on_curve(self.cref, u, dx, pnt, ref_pln)
