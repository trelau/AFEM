from OCC.TopAbs import TopAbs_REVERSED
from OCC.BRepAdaptor import BRepAdaptor_Curve
from OCC.GCPnts import GCPnts_AbscissaPoint, GCPnts_UniformAbscissa
from numpy import floor

from ..mesh_mgr import MeshMgr
from ..nodes import Node
from ...config import Settings
from ...geometry import CheckGeom
from ...topology import ShapeTools


def mesh_edge(edge, maxh=None, density=None, face=None):
    """
    Mesh an edge.
    """
    if maxh is None and density is None:
        return []
    npts = None
    if density is not None:
        npts = density + 1

    # Adaptor curve for edge.
    if face is None:
        adp_crv = BRepAdaptor_Curve(edge)
    else:
        adp_crv = BRepAdaptor_Curve(edge, face)

    # Adjust maxh and/or density if needed.
    if maxh is not None:
        arc_len = GCPnts_AbscissaPoint.Length(adp_crv, Settings.gtol)
        nb_pts = int(arc_len / maxh) + 1
    else:
        nb_pts = int(npts)

    # Minimum number of points if maxh and npts are provided.
    if maxh is not None and npts is not None:
        if nb_pts < npts:
            nb_pts = int(npts)

    # Force at least one element.
    if nb_pts <= 1:
        nb_pts = 2

    pac = GCPnts_UniformAbscissa(adp_crv, nb_pts, Settings.gtol)
    if not pac.IsDone():
        return []

    # Gather points and parameters.
    pnts = []
    prms = []
    for i in range(1, pac.NbPoints() + 1):
        u = pac.Parameter(i)
        gp_pnt = adp_crv.Value(u)
        pnt = CheckGeom.to_point(gp_pnt)
        pnts.append(pnt)
        prms.append(u)
    if len(pnts) < 2:
        return []

    # Create nodes from the points without duplicates.
    npts = len(pnts)
    nodes = []
    v1 = ShapeTools.first_vertex(edge)
    v2 = ShapeTools.last_vertex(edge)
    # Flip the vertices if the edge is reversed because BRepAdapter_Curve's
    # always follow the orientation of the curve.
    if edge.Orientation() == TopAbs_REVERSED:
        v1, v2 = v2, v1
    n1 = MeshMgr.mesh_from_shape(v1)
    n2 = MeshMgr.mesh_from_shape(v2)
    for i in range(npts):
        if i == 0:
            if isinstance(n1, Node):
                node = n1
            else:
                node = Node(t=prms[i], pnt=pnts[i])
                MeshMgr.shape_to_mesh(v1, node)
        elif i == npts - 1:
            if isinstance(n2, Node):
                node = n2
            else:
                if pnts[i].is_equal(pnts[0], Settings.gtol):
                    node = nodes[0]
                else:
                    node = Node(t=prms[i], pnt=pnts[i])
                MeshMgr.shape_to_mesh(v2, node)
        else:
            node = Node(t=prms[i], pnt=pnts[i])
        nodes.append(node)

    return nodes
