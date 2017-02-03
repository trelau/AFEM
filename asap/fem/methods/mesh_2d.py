import netgen.meshing as ng_mesh
from OCC.BRepAdaptor import BRepAdaptor_Curve2d, BRepAdaptor_Surface
from OCC.GeomAbs import GeomAbs_Plane
from OCC.TopAbs import TopAbs_REVERSED
from netgen.geom2d import SplineGeometry
from numpy import array, mean
from numpy.linalg import norm
from scipy.spatial import KDTree

from .map_face import FaceMap
from ..data import MeshData
from ..elements import Elm2D
from ..nodes import Node
from ...config import Settings
from ...geometry import CheckGeom
from ...topology import ShapeTools
from ...utils import pairwise


class Node2D(object):
    """
    2-D node for meshing faces.
    """

    def __init__(self, node3d, u=None, v=None, s=None, t=None, pnt=None):
        self.node3d = node3d
        self.u = u
        self.v = v
        self.s = s
        self.t = t
        self.pnt = pnt

    def __str__(self):
        return 'Node2D = ({0}, {1})'.format(self.u, self.v)

    def __eq__(self, other):
        return self.is_equal(other, Settings.mtol)

    @property
    def uv(self):
        return array([self.u, self.v])

    @property
    def st(self):
        return self.s, self.t

    def distance(self, other):
        """
        Distance to other node.
        """
        return norm(self.pnt - other.pnt)

    def is_equal(self, other, tol=None):
        """
        Check to see if nodes are equivalent.
        """
        if tol is None:
            tol = Settings.mtol
        if self.distance(other) <= tol:
            return True
        return False


def mesh_face(face, quad_dominated=False):
    """
    Mesh a face.
    """
    # Get the outer wire of the face and use a wire explorer to get the
    # boundary nodes.
    bnodes = []
    wire = ShapeTools.outer_wire(face)
    wexp = ShapeTools.wire_explorer(wire, face)
    while wexp.More():
        edge = ShapeTools.to_edge(wexp.Current())
        # Get 3-D nodes for the edge.
        edge_mesh = MeshData.mesh_from_shape(edge)
        nodes = edge_mesh.nodes

        # Reverse the node list if the edge is reversed.
        if wexp.Orientation() == TopAbs_REVERSED:
            nodes.reverse()

        # Parameters along edge.
        t = [n.t for n in nodes]

        # Get first and last parameters for the vertices of the edge and
        # update the list.
        u1, u2 = ShapeTools.parameters(edge, face)
        t[0], t[-1] = u1, u2

        # For each node except the last one, generate a 2-D node and add to
        # the list using an adaptor curve.
        adp_crv2d = BRepAdaptor_Curve2d(edge, face)
        for ni, ti in zip(nodes, t)[:-1]:
            gp_pnt2d = adp_crv2d.Value(ti)
            u, v = gp_pnt2d.X(), gp_pnt2d.Y()
            n2d = Node2D(ni, u, v, pnt=ni.pnt)
            bnodes.append(n2d)

        # Next edge.
        wexp.Next()

    if len(bnodes) < 3:
        return []

    # Check direction of boundary nodes are CCW.
    xy = [b.uv for b in bnodes]
    xy.append(xy[0])
    xy = array(xy)
    nm = xy.shape
    x0, y0 = xy[0]
    xn, yn = xy[-1]
    a = xn * y0 - x0 * yn
    for i in range(0, nm[0] - 1):
        a += xy[i, 0] * xy[i + 1, 1]
        a -= xy[i + 1, 0] * xy[i, 1]
    a *= 0.5
    # assert a > 0., "Boundary nodes need reversed"
    if a < 0:
        bnodes.reverse()

    # Min distance between boundary nodes for surface map.
    h = []
    for n1, n2 in pairwise(bnodes):
        h.append(n1.distance(n2))
    h.append(bnodes[-1].distance(bnodes[0]))
    h = mean(h)

    # Build face distance map.
    adp_srf = BRepAdaptor_Surface(face, True)
    smap = None
    if adp_srf.GetType() != GeomAbs_Plane:
        smap = FaceMap(face, True, 10, h)

    # Generate the geometry using Netgen points and line segments.
    geo = SplineGeometry()
    node_to_point = {}
    point_to_node = {}
    points = []
    for n in bnodes:
        if smap:
            s = smap.eval_udist(n.u, n.v)
            t = smap.eval_vdist(n.u, n.v)
        else:
            s, t = n.u, n.v
        p = geo.AppendPoint(s, t)
        node_to_point[n] = p
        point_to_node[p] = n
        points.append(p)
        n.s, n.t = s, t

    # Generate line segments from pairwise points.
    h = []
    for p1, p2 in pairwise(points):
        geo.Append(['line', p1, p2])
        n1 = point_to_node[p1]
        n2 = point_to_node[p2]
        h.append(norm(n1.pnt - n2.pnt))
    # Last segment.
    p1, p2 = points[-1], points[0]
    geo.Append(['line', p1, p2])
    n1 = point_to_node[p1]
    n2 = point_to_node[p2]
    h.append(norm(n1.pnt - n2.pnt))

    # Set warning level.
    if Settings.warnings:
        # noinspection PyUnresolvedReferences
        ng_mesh.SetMessageImportance(5)
    else:
        # noinspection PyUnresolvedReferences
        ng_mesh.SetMessageImportance(0)

    # Generate the mesh.
    h = max(h)
    # noinspection PyUnresolvedReferences
    mp = ng_mesh.MeshingParameters(maxh=h,
                                   quad_dominated=quad_dominated,
                                   optsteps2d=3,
                                   edge_subdivide=False)
    mesh = geo.GenerateMesh(mp=mp)

    # Get nodes.
    nodes = _get_nodes2d(bnodes, mesh, smap, adp_srf, Settings.mtol)

    # Get elements.
    elements = _get_elements(nodes, mesh)

    return elements


def _get_nodes2d(fixed_nodes2d, mesh, smap, adp_srf, tol):
    """
    Enumerate mesh nodes using a KDTree to correllate mesh nodes to original
    nodes. This is done since Netgen does not maintain numbering of input
    points.
    """
    # Build a KDTree of fixed nodes.
    data = [n.st for n in fixed_nodes2d]
    kdt = KDTree(data)

    # Convert 2-D mesh to 3-D points.
    p3d = []
    p_st = []
    for p in mesh.Points():
        s, t, _ = p.p
        p_st.append((s, t))
        if smap:
            u = smap.eval_uparam(s, t)
            v = smap.eval_vparam(s, t)
        else:
            u, v = s, t
        gp_pnt = adp_srf.Value(u, v)
        pnt = CheckGeom.to_point(gp_pnt)
        p3d.append(pnt)
    p_st = array(p_st)
    p3d = array(p3d)

    # Query KDTree using (s, t) points.
    _, indx = kdt.query(p_st, 1, 0., 2, tol)

    # For each mesh point, either use an existing node or create a new one.
    nodes2d = []
    n_nodes = len(fixed_nodes2d)
    n = 0
    for i in list(indx):
        if i < n_nodes:
            nodes2d.append(fixed_nodes2d[i])
        else:
            n3d = Node(pnt=p3d[n])
            n2d = Node2D(n3d)
            nodes2d.append(n2d)
        n += 1

    return nodes2d


def _get_elements(nodes, mesh):
    """
    Create elements from the mesh.
    """
    elms = []
    for e in mesh.Elements2D():
        # Element connectivity.
        conn = [pid.nr for pid in e.vertices]
        # Get nodes. Subract 1 from the index for Python indexing.
        ni = [nodes[i - 1].node3d for i in conn]
        # Create generic 2-D element.
        elms.append(Elm2D(ni))

    return elms
