# noinspection PyUnresolvedReferences
# import netgen.meshing as ng_mesh
from OCC.BRep import BRep_Tool
from OCC.BRepAdaptor import BRepAdaptor_Curve2d, BRepAdaptor_Surface
from OCC.GeomAbs import GeomAbs_Plane
from OCC.TopAbs import TopAbs_REVERSED
# noinspection PyUnresolvedReferences
# from netgen.geom2d import SplineGeometry
from numpy import array, mean
from numpy.linalg import norm
from scipy.spatial import KDTree

from .map_face import FaceMap
from ..elements import Elm2D
from ..mesh_mgr import MeshMgr
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

    @property
    def uv(self):
        return array([self.u, self.v])

    @property
    def st(self):
        return array([self.s, self.t])

    def distance(self, other):
        """
        Distance to other node.
        """
        return norm(self.uv - other.uv)

    def is_equal(self, other, tol=None):
        """
        Check to see if nodes are equivalent.
        """
        if tol is None:
            tol = Settings.ptol
        if self.distance(other) <= tol:
            return True
        return False


def mesh_face(face, quad_dominated=False):
    """
    Mesh a face.
    """
    # Average h.
    havg = []

    # Get the outer wire of the face and use a wire explorer to get the
    # boundary nodes.
    bwire = ShapeTools.outer_wire(face)
    bnodes2d = _process_closed_wire(bwire, face)

    if len(bnodes2d) < 3:
        return []

    # Check direction of boundary nodes.
    _check_order(bnodes2d, True)

    # Average distance between boundary nodes for surface map.
    hb = _havg(bnodes2d)
    havg.append(hb)

    # Explore for closed wires that define a hole.
    wires = ShapeTools.get_wires(face)
    internal_wires = []
    hnodes2d = []
    for wire in wires:
        # Make sure this is not the boundary wire.
        if bwire.IsSame(wire):
            continue

        # Check if wire is closed. If not, add it to internal wires for
        # later processing.
        if not BRep_Tool.IsClosed(wire):
            internal_wires.append(wire)
            continue

        # Process nodes on wire, check order, and calculate havg.
        nodes2d = _process_closed_wire(wire, face)
        _check_order(nodes2d, False)
        havg.append(_havg(nodes2d))

        # Add to hole list.
        hnodes2d.append(nodes2d)

    # Calculate average h.
    h = mean(havg)

    # Build face distance map.
    adp_srf = BRepAdaptor_Surface(face, True)
    smap = None
    if adp_srf.GetType() != GeomAbs_Plane:
        smap = FaceMap(face, True, 10, h)

    # Generate the geometry using Netgen points and line segments.
    geo = SplineGeometry()
    node2d_to_point = {}
    point_to_node2d = {}
    points = []

    # Boundary
    for n2d in bnodes2d:
        if smap:
            s = smap.eval_udist(n2d.u, n2d.v)
            t = smap.eval_vdist(n2d.u, n2d.v)
        else:
            s, t = n2d.u, n2d.v
        p = geo.AppendPoint(s, t)
        node2d_to_point[n2d] = p
        point_to_node2d[p] = n2d
        points.append(p)
        n2d.s, n2d.t = s, t

    # Generate line segments from pairwise points.
    h = []
    h2d = []
    for p1, p2 in pairwise(points):
        geo.Append(['line', p1, p2])
        n1 = point_to_node2d[p1]
        n2 = point_to_node2d[p2]
        h.append(norm(n1.pnt - n2.pnt))
        h2d.append(n1.distance(n2))
    # Last segment.
    p1, p2 = points[-1], points[0]
    geo.Append(['line', p1, p2])
    n1 = point_to_node2d[p1]
    n2 = point_to_node2d[p2]
    h.append(norm(n1.pnt - n2.pnt))
    h2d.append(n1.distance(n2))

    # Holes
    for hnode in hnodes2d:
        points = []
        for n2d in hnode:
            if smap:
                s = smap.eval_udist(n2d.u, n2d.v)
                t = smap.eval_vdist(n2d.u, n2d.v)
            else:
                s, t = n2d.u, n2d.v
            p = geo.AppendPoint(s, t)
            node2d_to_point[n2d] = p
            point_to_node2d[p] = n2d
            points.append(p)
            n2d.s, n2d.t = s, t

        for p1, p2 in pairwise(points):
            geo.Append(['line', p1, p2])
            n1 = point_to_node2d[p1]
            n2 = point_to_node2d[p2]
            h.append(norm(n1.pnt - n2.pnt))
            h2d.append(n1.distance(n2))
        # Last segment.
        p1, p2 = points[-1], points[0]
        geo.Append(['line', p1, p2])
        n1 = point_to_node2d[p1]
        n2 = point_to_node2d[p2]
        h.append(norm(n1.pnt - n2.pnt))
        h2d.append(n1.distance(n2))

    # Gather existing node 2-D instances.
    all_n2d = node2d_to_point.keys()

    # Calculate a tolerance in 2-D to check for equivalent nodes.
    tol2d = min(h2d) / 10.

    # Internal wires
    inodes2d = []
    for wire in internal_wires:
        # Process nodes on wire.
        nodes2d = _process_internal_wire(wire, face)

        # Check if any created 2-D node can be replaced with existing.
        for i, n2d1 in enumerate(nodes2d):
            for n2d2 in all_n2d:
                if n2d1.is_equal(n2d2, tol2d):
                    nodes2d[i] = n2d2
                    break

        # Add to list.
        inodes2d.append(nodes2d)

    # Generate internal line segments.
    for inode in inodes2d:
        points = []
        boundaries = []
        for n2d in inode:
            is_boundary = False
            if n2d in node2d_to_point:
                is_boundary = True
            elif smap:
                s = smap.eval_udist(n2d.u, n2d.v)
                t = smap.eval_vdist(n2d.u, n2d.v)
            else:
                s, t = n2d.u, n2d.v

            if is_boundary:
                p = node2d_to_point[n2d]
            else:
                p = geo.AppendPoint(s, t)
                node2d_to_point[n2d] = p
                point_to_node2d[p] = n2d
                n2d.s, n2d.t = s, t
            points.append(p)
            boundaries.append(is_boundary)

        for p1, p2 in pairwise(points):
            # The zero is for the boundary condition.
            geo.Append(['line', p1, p2], 0)

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
    all_n2d = []
    for n2d in bnodes2d:
        all_n2d.append(n2d)
    for hnode in hnodes2d:
        for n2d in hnode:
            all_n2d.append(n2d)
    for inode in inodes2d:
        for n2d in inode:
            if n2d not in all_n2d:
                all_n2d.append(n2d)

    nodes2d = _get_nodes2d(all_n2d, mesh, smap, adp_srf, Settings.mtol)

    # Get elements.
    elements = _get_elements(nodes2d, mesh)

    return elements


def _process_closed_wire(wire, face):
    """
    Process a closed wire.
    """
    wire_nodes2d = []
    wexp = ShapeTools.wire_explorer(wire, face)
    while wexp.More():
        edge = ShapeTools.to_edge(wexp.Current())
        # Get 3-D nodes for the edge.
        edge_mesh = MeshMgr.mesh_from_shape(edge)
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
        for ni, ti in list(zip(nodes, t))[:-1]:
            gp_pnt2d = adp_crv2d.Value(ti)
            u, v = gp_pnt2d.X(), gp_pnt2d.Y()
            n2d = Node2D(ni, u, v, pnt=ni.pnt)
            wire_nodes2d.append(n2d)

        # Next edge.
        wexp.Next()

    return wire_nodes2d


def _process_internal_wire(wire, face):
    """
    Process an internal wire.
    """
    wire_nodes2d = []
    wexp = ShapeTools.wire_explorer(wire, face)

    # Put all edges in ordered list. This is a hack since wire explorer
    # seems to crash if wire only has one edge...
    all_edges = ShapeTools.get_edges(wire, False)
    edges = []
    orienation = []
    if len(all_edges) < 1:
        return []
    elif len(all_edges) == 1:
        edge = ShapeTools.to_edge(all_edges[0])
        edges.append(edge)
        orienation.append(edge.Orientation())
    else:
        while wexp.More():
            edge = ShapeTools.to_edge(wexp.Current())
            edges.append(edge)
            orienation.append(wexp.Orientation())
            wexp.Next()

    nedges = len(edges)
    for i, edge in enumerate(edges):
        # Get 3-D nodes for the edge.
        edge_mesh = MeshMgr.mesh_from_shape(edge)
        nodes = edge_mesh.nodes

        # Reverse the node list if the edge is reversed.
        if orienation[i] == TopAbs_REVERSED:
            nodes.reverse()

        # Parameters along edge.
        t = [n.t for n in nodes]

        # Get first and last parameters for the vertices of the edge and
        # update the list.
        u1, u2 = ShapeTools.parameters(edge, face)
        t[0], t[-1] = u1, u2

        # For each node, generate a 2-D node and add to the list using an
        # adaptor curve.
        adp_crv2d = BRepAdaptor_Curve2d(edge, face)
        if i == nedges - 1:
            for ni, ti in list(zip(nodes, t)):
                gp_pnt2d = adp_crv2d.Value(ti)
                u, v = gp_pnt2d.X(), gp_pnt2d.Y()
                n2d = Node2D(ni, u, v, pnt=ni.pnt)
                wire_nodes2d.append(n2d)
        else:
            for ni, ti in list(zip(nodes, t))[:-1]:
                gp_pnt2d = adp_crv2d.Value(ti)
                u, v = gp_pnt2d.X(), gp_pnt2d.Y()
                n2d = Node2D(ni, u, v, pnt=ni.pnt)
                wire_nodes2d.append(n2d)

    return wire_nodes2d


def _check_order(nodes2d, is_boundary):
    """
    Check the order of the nodes2d along the wire.
    """
    xy = [n.uv for n in nodes2d]
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
    # Boundary nodes2d should be CCW for Netgen, CW for holes.
    if a < 0 and is_boundary:
        nodes2d.reverse()
    elif a > 0. and not is_boundary:
        nodes2d.reverse()


def _havg(nodes):
    """
    Calculate average distance between nodes in 3-D.
    """
    h = []
    for n1, n2 in pairwise(nodes):
        h.append(n1.pnt.distance(n2.pnt))
    h.append(nodes[-1].pnt.distance(nodes[0].pnt))
    return mean(h)


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
