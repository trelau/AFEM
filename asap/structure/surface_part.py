from OCC.TopAbs import TopAbs_COMPSOLID, TopAbs_SOLID
from OCC.TopoDS import TopoDS_Shape

from .methods.build_parts import build_surface_part
from .methods.explore_parts import get_shared_edges, get_shared_nodes
from .methods.form_parts import form_with_solid
from .methods.join_parts import join_surface_parts
from .methods.mesh_parts import mesh_surface_part
from .methods.modify_parts import cut_part, discard_faces_by_distance, \
    discard_faces_by_solid
from .part import Part
from ..fem.data import MeshData
from ..topology import ShapeTools


class SurfacePart(Part):
    """
    Base class for surface-based parts.
    """

    def __init__(self, name, rshape):
        super(SurfacePart, self).__init__(name)
        self._rshape = rshape
        self._fshapes = set()

    @property
    def rshape(self):
        return self._rshape

    @property
    def fshapes(self):
        return list(self._fshapes)

    @property
    def faces(self):
        return ShapeTools.get_faces(self)

    @property
    def edges(self):
        return ShapeTools.get_edges(self)

    @property
    def nfaces(self):
        return len(self.faces)

    @property
    def elements(self):
        elm_set = set()
        for f in self.faces:
            if not MeshData.has_mesh(f):
                continue
            mesh = MeshData.mesh_from_shape(f)
            for e in mesh.elements:
                elm_set.add(e)
        return elm_set

    @property
    def nodes(self):
        node_set = set()
        for e in self.elements:
            for n in e.nodes:
                node_set.add(n)
        return node_set

    def form(self, *bodies):
        """
        Form frame boundary with a body.

        :param bodies:

        :return:
        """
        status = {}
        for body in bodies:
            shape = form_with_solid(self._rshape, body)
            if not shape:
                status[body] = False
            else:
                self._fshapes.add(shape)
                status[body] = True
        return status

    def build(self, unify=False):
        """
        Build the frame shape.

        :param bool unify:

        :return:
        """
        return build_surface_part(self, unify)

    def join(self, *other_parts):
        """
        Joint with other part.

        :param other_parts:

        :return:
        """
        _other_parts = []
        for part in other_parts:
            if isinstance(part, TopoDS_Shape) and not part.IsNull():
                _other_parts.append(part)
        if not _other_parts:
            return False
        return join_surface_parts(self, *_other_parts)

    def cut(self, cutter):
        """
        Cut the part with a shape.

        :param cutter:

        :return:
        """
        return cut_part(self, cutter)

    def discard(self, shape, tol=None):
        """
        Discard faces of the part.

        :param shape:
        :param tol:

        :return:
        """
        shape = ShapeTools.to_shape(shape)
        if not shape:
            return False
        if shape.ShapeType() in [TopAbs_SOLID, TopAbs_COMPSOLID]:
            return discard_faces_by_solid(self, shape, tol)
        return discard_faces_by_distance(self, shape, tol)

    def shared_edges(self, other_part):
        """
        Get edges shared between the two parts.

        :param other_part:

        :return:
        """
        return get_shared_edges(self, other_part)

    def mesh(self, maxh=4., quad_dominated=True):
        """
        Mesh the part.

        :param maxh:
        :param quad_dominated:

        :return:
        """
        return mesh_surface_part(self, maxh, quad_dominated)

    def shared_nodes(self, other_part):
        """
        Get nodes shared between the two parts.

        :param other_part:

        :return:
        """
        return get_shared_nodes(self, other_part)
