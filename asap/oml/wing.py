from .body import Body
from .methods.extract import extract_wing_plane, extract_wing_ref_curve
from ..geometry import CheckGeom, CreateGeom, ProjectGeom
from ..topology.adaptors import ShapeAdaptorSurface


class Wing(Body):
    """
    Wing body.
    """

    def __init__(self, shape):
        super(Wing, self).__init__(shape)
        self._sref = None
        self._sref_adp = None

    @property
    def sref(self):
        return self._sref

    @property
    def sref_shape(self):
        return self._sref_adp.shape

    @property
    def u1(self):
        return self.sref.u1

    @property
    def v1(self):
        return self.sref.v1

    @property
    def u2(self):
        return self.sref.u2

    @property
    def v2(self):
        return self.sref.v2

    @property
    def uknots(self):
        return self.sref.uknots

    @property
    def vknots(self):
        return self.sref.vknots

    def set_sref(self, sref):
        """
        Set the wing reference surface.

        :param sref:
        :return:
        """
        if CheckGeom.is_surface_like(sref):
            self._sref = sref
            self._sref_adp = ShapeAdaptorSurface()
            self._sref_adp.build(sref)
            return True
        return False

    def eval(self, u, v):
        """
        Evaluate a point on the wing reference surface.

        :param float u: Reference surface parameter in u-direction.
        :param float v: Reference surface parameter in v-direction.

        :return: Point on reference surface.
        :rtype: :class:`.Point`
        """
        return self.sref.eval(u, v)

    def norm(self, u, v):
        """
        Evaluate the surface normal of the wing reference surface.

        :param float u: Reference surface parameter in u-direction.
        :param float v: Reference surface parameter in v-direction.

        :return: Reference surface normal.
        :rtype: :class:`.Vector`
        """
        return self.sref.norm(u, v)

    def extract_curve(self, uv1, uv2, rshape=None):
        """
        Extract curve along wing reference surface.

        :param uv1:
        :param uv2:
        :param rshape:

        :return:
        """
        return extract_wing_ref_curve(self, uv1, uv2, rshape)

    def extract_plane(self, uv1, uv2):
        """
        Extract a plane between the two points.  Points are provided as either
        wing parametric coordinates (u, v) or :class:`.Point` instances.

        :param uv1: Starting parameters (or point).
        :type uv1: :class:`.Point` or array_like
        :param uv2: Ending parameters (or point).
        :type uv2: :class:`.Point` or array_like

        :return: Plane oriented between start and end points.
        :rtype: :class:`.Plane`
        """
        return extract_wing_plane(self, uv1, uv2)

    def invert_point(self, p):
        """
        Find the wing reference surface parameters by inverting the point.

        :param p: The point.
        :type p: :class:`.Point` or array_like

        :return: Parameters on the wing reference surface (u, v). Will return
            (None, None) if the method fails.
        :rtype: tuple
        """
        p = CheckGeom.to_point(p)
        return ProjectGeom.invert(p, self.sref)

    def isocurve(self, u=None, v=None):
        """
        Extract isocurve in wing reference surface.

        :param u:
        :param v:

        :return:
        """
        return CreateGeom.isocurve(self.sref, u, v)
