from OCCT.TopoDS import TopoDS_Face

from afem.geometry import (PlaneByAxes, NurbsCurve2DByApprox,
                           NurbsCurve2DByInterp, Point2D)
from afem.topology import (CheckShape, FuseShapes, WiresByConnectedEdges,
                           FaceByPlanarWire)

__all__ = ["CrossSection"]


class CrossSection(object):
    """
    Planar cross section.

    :param afem.geometry.entities.Plane pln: The default construction plane. If
        *None* is provided, then the xy-plane is used.
    """

    def __init__(self, pln=None):
        if pln is None:
            pln = PlaneByAxes(axes='xy').plane

        self._pln = pln
        self._crvs = []
        self._shape = None
        self._wire_tool = None
        self._face = None

    @property
    def shape(self):
        """
        :return: The 3-D shape after building. If only one edge was available
            then an edge is returned, otherwise the edges are fused together
            and the resulting shape is returned.
        :rtype: OCCT.TopoDS.TopoDS_Edge or OCCT.TopoDS.TopoDS_Shape
        """
        return self._shape

    @property
    def nwires(self):
        """
        :return: Number of wires after building.
        :rtype: int
        """
        return self._wire_tool.nwires

    @property
    def wires(self):
        """
        :return: The wires after building.
        :rtype: list(OCCT.TopoDS.TopoDS_Wire)
        """
        return self._wire_tool.wires

    @property
    def has_face(self):
        """
        :return: Check if face is available.
        :rtype: bool
        """
        return isinstance(self._face, TopoDS_Face)

    @property
    def face(self):
        """
        :return: The face if available.
        :rtype: OCCT.TopoDS.TopoDS_Face
        """
        return self._face

    def _close(self, c):
        """
        Close a curve segment.
        """
        if c.is_closed:
            return None

        p1 = c.eval(c.u1)
        p2 = c.eval(c.u2)
        self.add_segment(p1, p2)

    def add_segment(self, p1, p2):
        """
        Add linear segment between the two points.

        :param point2d_like p1: The first point.
        :param point2d_like p2: The second point.

        :return: The 2-D curve.
        :rtype: afem.geometry.entities.NurbsCurve2D
        """
        return self.add_interp([p1, p2])

    def add_approx(self, pnts, close=False):
        """
        Add a curve by approximating the 2-D points.

        :param collections.Sequence(point2d_like) pnts: The points.
        :param bool close: Option to add a segment to close the 2-D curve.

        :return: The 2-D curve.
        :rtype: afem.geometry.entities.NurbsCurve2D
        """
        c = NurbsCurve2DByApprox(pnts).curve
        self._crvs.append(c)

        if close:
            self._close(c)

        return c

    def add_interp(self, pnts, close=False):
        """
        Add a curve by interpolating the 2-D points.

        :param collections.Sequence(point2d_like) pnts: The points.
        :param bool close: Option to add a segment to close the 2-D curve.

        :return: The 2-D curve.
        :rtype: afem.geometry.entities.NurbsCurve2D
        """
        c = NurbsCurve2DByInterp(pnts).curve
        self._crvs.append(c)

        if close:
            self._close(c)

        return c

    def rotate(self, angle, pnt=(0., 0.)):
        """
        Rotate each 2-D curve in the cross section.

        :param float angle: The angle in degrees.
        :param point2d_like pnt: The reference point for rotation.

        :return: None.
        """
        for c in self._crvs:
            c.rotate(pnt, angle)

    def scale(self, scale, pnt=(0., 0.)):
        """
        Scale each 2-D curve in the cross section.

        :param float scale: The scaling parameter.
        :param point2d_like pnt: The reference point for scaling.

        :return: None.
        """
        for c in self._crvs:
            c.scale(pnt, scale)

    def copy(self, pln=None):
        """
        Copy the cross section and underlying 2-D curves.

        :param afem.geometry.entities.Plane pln: The default plane for the
            new cross section. If *None*, then the default plane for this
            cross section is used.

        :return: A new cross section.
        :rtype: afem.sketch.entities.CrossSection
        """
        if pln is None:
            pln = self._pln

        new_cs = CrossSection(pln)
        for c in self._crvs:
            new_cs._crvs.append(c.copy())

        return new_cs

    def read_uiuc(self, fn, close=True):
        """
        Generate 2-D curves by reading an airfoil file from the UIUC database.

        :param str fn: The filename.
        :param bool close: Option to close the airfoil by adding a segment.

        :return: The 2-D curve.
        :rtype: afem.geometry.entities.NurbsCurve2D

        .. note::

            The UIUC airfoil is assumed to be in a particular format and this
            method may fail it's not. The points for the upper and lower
            surfaces from the file are combined to form a single curve that
            starts at the trailing edge, moves along the upper surface towards
            the leading edge, and then proceeds back to the trailing edge along
            the lower surface.
        """
        with open(fn, 'r') as fin:
            content = fin.read().splitlines()

        upr, lwr = [], []
        i = 3
        while i < len(content):
            line = content[i].strip().split()
            i += 1
            if not line:
                break
            x = float(line[0])
            z = float(line[1])
            upr.append(Point2D(x, z))

        while i < len(content):
            line = content[i].strip().split()
            i += 1
            if not line:
                break
            x = float(line[0])
            z = float(line[1])
            lwr.append(Point2D(x, z))

        upr.reverse()
        pnts = upr + lwr[1:]

        return self.add_approx(pnts, close)

    def build(self, pln=None, scale=None, rotate=None):
        """
        Build a shape in 3-D using the 2-D cross section curves. This method
        collects the 2-D curves, converts them into 3-D edges on the plane,
        fuses them together, and attempts to build wires and a face if
        applicable.

        :param pln: The plane to build on. If *None*, then the default plane
            is used.
        :param float scale: Scale the 3-D curve after construction on the
            plane. The reference point is the plane origin.
        :param float rotate: Rotate the 3-D curve after construction on the
            plane. The reference point is the plane origin.

        :return: *True* if successful, *False* otherwise.
        :rtype: bool
        """
        if pln is None:
            pln = self._pln
        if pln is None:
            raise RuntimeError('No plane is defined.')
        origin = pln.origin
        axis = pln.axis

        edges = []
        for c2d in self._crvs:
            c3d = c2d.to_3d(pln)
            if scale is not None:
                c3d.scale(origin, scale)
            if rotate is not None:
                c3d.rotate(axis, rotate)
            e = CheckShape.to_edge(c3d)
            edges.append(e)

        if len(edges) == 0:
            return False

        if len(edges) == 1:
            self._shape = edges[0]
        else:
            # Fuse the shape
            bop = FuseShapes()
            bop.set_args(edges[:-1])
            bop.set_tools(edges[-1:])
            bop.build()
            self._shape = bop.shape
            edges = bop.edges

        # Try and make wire(s)
        self._wire_tool = WiresByConnectedEdges(edges)

        # Try to make a face if one wire was created
        if self.nwires == 1:
            self._face = FaceByPlanarWire(self.wires[0]).face

        return True
