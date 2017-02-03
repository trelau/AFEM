from OCC.BRepAdaptor import BRepAdaptor_Surface
from OCC.BRepBndLib import brepbndlib_Add
from OCC.BRepMesh import BRepMesh_IncrementalMesh
from OCC.BRepTools import breptools_UVBounds
from OCC.Bnd import Bnd_Box
from numpy import float64, linspace, sqrt, zeros
from numpy.linalg import norm
from scipy.interpolate import LinearNDInterpolator, RectBivariateSpline

from ...geometry.methods.create import create_nurbs_surface_from_occ


class FaceMap(object):
    """
    Face distance mapping class.
    """

    def __init__(self, face, build=True, n=10, h=None):
        self._f = face
        self._interp_udist = None
        self._interp_vdist = None
        self._interp_uparam = None
        self._interp_vparam = None
        if build:
            self.build_distance_map(n, h)

    def build_distance_map(self, n=10, h=None):
        """
        Build a distance map of the face.

        :param int n:
        :param float h:

        :return:
        """
        results = face_distance_map(self._f, n, h)
        self._interp_udist = results[0]
        self._interp_vdist = results[1]
        self._interp_uparam = results[2]
        self._interp_vparam = results[3]
        return True

    def eval_udist(self, u, v):
        """
        Evaluate u-distance given parameters.

        :param u:
        :param v:

        :return:
        """
        try:
            return float(self._interp_udist(u, v)[0, 0])
        except (AttributeError, IndexError):
            return None

    def eval_vdist(self, u, v):
        """
        Evaluate v-distance given parameters.

        :param u:
        :param v:

        :return:
        """
        try:
            return float(self._interp_vdist(u, v)[0, 0])
        except (AttributeError, IndexError):
            return None

    def eval_uparam(self, s, t):
        """
        Evaluate u-parameter given distances.

        :param s:
        :param t:

        :return:
        """
        return float(self._interp_uparam(s, t))

    def eval_vparam(self, s, t):
        """
        Evaluate t-parameter given distance.

        :param s:
        :param t:

        :return:
        """
        return float(self._interp_vparam(s, t))


def face_distance_map(face, n=None, h=None):
    """
    Build a distance map of the face.
    """
    if n is None and h is None:
        return None, None, None, None

    # Triangulate the face if not already done so a smaller bounding box is
    # calculated.
    BRepMesh_IncrementalMesh(face, 1., True)

    # Use length of bounding box diagonal to estimate grid steps.
    if h is not None:
        bbox = Bnd_Box()
        brepbndlib_Add(face, bbox)
        diag = sqrt(bbox.SquareExtent())
        ntemp = int(diag / h) + 1
        if ntemp > n:
            n = ntemp

    # Extract an adaptor surface from the face and convert.
    adp_srf = BRepAdaptor_Surface(face, True)
    surface = create_nurbs_surface_from_occ(adp_srf.BSpline().GetObject())

    # Parameter domain of the face on the surface.
    umin, umax, vmin, vmax = breptools_UVBounds(face)

    # Grid of uv values.
    ugrid = linspace(umin, umax, n)
    vgrid = linspace(vmin, vmax, n)

    # Number of grids.
    nu = len(ugrid)
    nv = len(vgrid)

    seval = surface.eval
    # Map u-direction.
    udist = zeros((nu, nv), dtype=float64)
    for j in range(0, nv):
        dv = vgrid[j]
        d = 0.
        for i in range(1, nu):
            du_0 = ugrid[i - 1]
            du_1 = ugrid[i]
            pi = seval(du_0, dv)
            pi1 = seval(du_1, dv)
            d += norm(pi1 - pi)
            udist[i, j] = d

    # Map v-direction.
    vdist = zeros((nu, nv), dtype=float64)
    for i in range(0, nu):
        du = ugrid[i]
        d = 0.
        for j in range(1, nv):
            dv_0 = vgrid[j - 1]
            dv_1 = vgrid[j]
            pi = seval(du, dv_0)
            pi1 = seval(du, dv_1)
            d += norm(pi1 - pi)
            vdist[i, j] = d

    # Build interpolation functions for parameters.
    xy = []
    uparam = []
    vparam = []
    for j in range(0, nv):
        for i in range(0, nu):
            xy.append([udist[i, j], vdist[i, j]])
            uparam.append(ugrid[i])
            vparam.append(vgrid[j])

    # Build interpolation functions for distance.
    interp_udist = RectBivariateSpline(ugrid, vgrid, udist, kx=1, ky=1)
    interp_vdist = RectBivariateSpline(ugrid, vgrid, vdist, kx=1, ky=1)
    interp_uparam = LinearNDInterpolator(xy, uparam, rescale=True)
    interp_vparam = LinearNDInterpolator(xy, vparam, rescale=True)

    return interp_udist, interp_vdist, interp_uparam, interp_vparam
