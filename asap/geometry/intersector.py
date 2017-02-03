from numpy import float64, inf, zeros
from scipy.spatial import KDTree

from .checker import CheckGeom
from .methods.distance import curve_nearest_point
from .methods.intersect import intersect_curve_curve, \
    intersect_curve_surface, intersect_surface_surface
from ..config import Settings


class IntersectGeom(object):
    """
    Intersect geometry.
    """

    @staticmethod
    def perform(geom1, geom2, itol=None):
        """
        Perform the intersection of two geometries.

        :param geom1: Geometry entity 1.
        :param geom2: Geometry entity 2.
        :param float itol: Intersection tolerance.

        :return: Intersector object depending on type of *geom1* and *geom2*.
            *None* if returned if entities are not supported.
        """
        # Curve intersection with...
        if CheckGeom.is_curve_like(geom1):
            # Curve
            if CheckGeom.is_curve_like(geom2):
                return IntersectCurveCurve(geom1, geom2, itol)
            # Surface/plane
            if CheckGeom.is_surface_like(geom2):
                return IntersectCurveSurface(geom1, geom2)

        # Surface intersection with...
        if CheckGeom.is_surface_like(geom1):
            # Curve
            if CheckGeom.is_curve_like(geom2):
                return IntersectCurveSurface(geom2, geom1)
            # Surface/plane
            if CheckGeom.is_surface_like(geom2):
                return IntersectSurfaceSurface(geom1, geom2, itol)

        # Return error if combination not supported.
        return IntersectError()


class IntersectError(object):
    """
    Class for handling intersection errors.
    """

    @property
    def success(self):
        return False

    @property
    def npts(self):
        return 0

    @property
    def ncrvs(self):
        return 0

    @property
    def points(self):
        return []

    @property
    def icrvs(self):
        return []


class CurveIntersector(object):
    """
    Base class for handling curve intersection methods and results.
    """

    def __init__(self, c1, c2):
        self._c1 = c1
        self._c2 = c2
        self._npts = 0
        self._results = []
        self._kdt = None

    def _set_results(self, npts, results):
        """
        Set curve intersection results.
        """
        if npts > 0:
            self._npts = npts
            # Replace ndarrays with point instances and build kd-tree.
            data = zeros((npts, 3), dtype=float64)
            for i in range(npts):
                data[i] = results[i][1]
            self._results = results
            self._kdt = KDTree(data)

    @property
    def npts(self):
        return self._npts

    @property
    def success(self):
        if self._npts > 0:
            return True
        return False

    @property
    def points(self):
        if self._npts <= 0:
            return []
        return [results[1] for results in self._results]

    @property
    def parameters(self):
        return [results[0] for results in self._results]

    def point(self, indx=1):
        """
        Return the point result by index.

        :param int indx: Index for point selection.

        :return: Intersection point at index.
        :rtype: :class:`.Point`
        """
        if indx > self._npts:
            return self._results[-1][1]
        return self._results[indx - 1][1]

    def query_point(self, p0, distance_upper_bound=inf):
        """
        Find the intersection result nearest to the provided point.

        :param p0: Point to search from.
        :type p0: :class:`.Point` or array_like
        :param distance_upper_bound: Return only results within this distance.
        :type distance_upper_bound: nonnegative float

        :return: Distance to nearest intersection result and its index (d, i).
            Returns (None, None) if no results are available.
        :rtype: tuple
        """
        if not self.success:
            return None, None
        if CheckGeom.is_point(p0):
            p0 = p0.xyz
        d, i = self._kdt.query(p0, 1, 0., 2, distance_upper_bound)
        return d, i

    def params_by_cref(self, cref):
        """
        Get the parameters of the intersection results by curve reference.

        :param cref: Reference curve (must have been using in the
            intersection method).

        :return: List of parameters for the reference curve. Returns empty
            list if  reference curve is not in the intersection.
        :rtype: list
        """
        if self._c1 is cref:
            return [results[0][0] for results in self._results]
        elif self._c2 is cref:
            if CheckGeom.is_surface_like(cref):
                return [results[0][1:] for results in self._results]
            return [results[0][1] for results in self._results]
        return []


class IntersectCurveCurve(CurveIntersector):
    """
    Curve-curve intersection.
    """

    def __init__(self, curve1, curve2, itol=None):
        super(IntersectCurveCurve, self).__init__(curve1, curve2)
        if CheckGeom.is_curve_like(curve1) and CheckGeom.is_curve_like(curve2):
            self._perform(curve1, curve2, itol)

    def _perform(self, curve1, curve2, itol):
        """
        Perform the curve-curve intersection.
        """
        npts, results = intersect_curve_curve(curve1, curve2, itol)
        self._set_results(npts, results)


class IntersectCurveSurface(CurveIntersector):
    """
    Curve-surface intersection.
    """

    def __init__(self, curve, surface):
        super(IntersectCurveSurface, self).__init__(curve, surface)
        if CheckGeom.is_curve_like(curve) and \
                CheckGeom.is_surface_like(surface):
            self._perform(curve, surface)

    def _perform(self, curve, surface):
        """
        Perform the curve-surface intersection.
        """
        npts, results = intersect_curve_surface(curve, surface)
        self._set_results(npts, results)

    @property
    def curve_parameters(self):
        return [results[0][0] for results in self._results]

    @property
    def surface_parameters(self):
        return [results[0][1:] for results in self._results]


class SurfaceIntersector(object):
    """
    Base class for handling surface intersection methods and results.
    """

    def __init__(self, itol=None):
        if itol is None:
            itol = Settings.gtol
        self._itol = itol
        self._icrvs = []

    @property
    def ncrvs(self):
        return len(self._icrvs)

    @property
    def success(self):
        if self.ncrvs > 0:
            return True
        return False

    @property
    def itol(self):
        return self._itol

    @property
    def icurves(self):
        return self._icrvs

    def get_icurve(self, indx=1):
        """
        Generate an :class:`.ICurve` for the specified intersection curve.

        :param int indx: Index of intersection curve.

        :return: Intersection curve.
        :rtype: :class:`.ICurve`
        """
        try:
            return self._icrvs[indx - 1]
        except IndexError:
            return None

    def get_icurves(self):
        """
        Get all intersection curves.

        :return: List of intersection curves.
        :rtype: list
        """
        return self._icrvs

    def curve_nearest_point(self, pref):
        """
        Get the index of the intersection curve that is nearest to the given
        reference point.

        :param array_like pref: Reference point.

        :return: Index of curve nearest point.
        :rtype: int
        """
        return curve_nearest_point(pref, self.icurves)


class IntersectSurfaceSurface(SurfaceIntersector):
    """
    Surface-surface intersection.
    """

    def __init__(self, surface1, surface2, itol=None):
        super(IntersectSurfaceSurface, self).__init__(itol)
        if CheckGeom.is_surface_like(surface1) and \
                CheckGeom.is_surface_like(surface2):
            self._perform(surface1, surface2)

    def _perform(self, surface1, surface2):
        """
        Perform the surface-surface intersection.
        """
        icrvs = intersect_surface_surface(surface1, surface2, self.itol)
        self._icrvs = icrvs
