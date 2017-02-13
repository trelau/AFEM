from .methods.cut_parts import cut_wing_part_with_circle
from .methods.modify_parts import discard_wing_part_faces
from .surface_part import SurfacePart
from ..geometry import CheckGeom, CreateGeom, ProjectGeom
from ..topology import ShapeTools


class WingPart(SurfacePart):
    """
    Base class for wing parts.
    """

    def __init__(self, name, wing, rshape):
        super(WingPart, self).__init__(name, rshape)
        self._wing = wing
        self._cref = None

    @property
    def cref(self):
        return self._cref

    @property
    def p1(self):
        try:
            return self._cref.eval(self._cref.u1)
        except AttributeError:
            return None

    @property
    def p2(self):
        try:
            return self._cref.eval(self._cref.u2)
        except AttributeError:
            return None

    def set_cref(self, cref):
        """
        Set reference curve of frame.

        :param cref:

        :return:
        """
        if CheckGeom.is_curve_like(cref):
            self._cref = cref
            return True
        return False

    def local_to_global(self, u):
        """
        Convert local parameter from 0 <= u <= 1 to u1 <= u <= u2.

        :param u:

        :return:
        """
        try:
            return self._cref.local_to_global_param(u)
        except AttributeError:
            return None

    def eval(self, u):
        """
        Evaluate point on reference curve.

        :param u:
        :return:
        """
        try:
            return self._cref.eval(u)
        except AttributeError:
            return None

    def distribute_points(self, dx, s1=None, s2=None, u1=None, u2=None,
                          shape1=None, shape2=None):
        """
        Create evenly distributed points along the reference curve.

        :param dx:
        :param float s1:
        :param float s2:
        :param float u1:
        :param float u2:
        :param shape1:
        :param shape2:

        :return:
        """
        if isinstance(dx, int):
            npts = dx
            dx = None
        else:
            npts = None
        return ShapeTools.points_along_shape(self.cref, dx, npts, u1, u2,
                                             s1, s2, shape1, shape2)

    def project_points(self, pnts):
        """
        Project points to reference curve.
        """
        for p in pnts:
            ProjectGeom.point_to_geom(p, self.cref, True)

    def invert(self, pnt):
        """
        Invert the point on the reference curve.

        :param pnt:

        :return:
        """
        return ProjectGeom.invert(pnt, self.cref)

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
        Cut a circular hole in the part.

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
