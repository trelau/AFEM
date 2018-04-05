# This file is part of AFEM which provides an engineering toolkit for airframe
# finite element modeling during conceptual design.
#
# Copyright (C) 2016-2018  Laughlin Research, LLC (info@laughlinresearch.com)
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
from afem.exchange.xde import XdeDocument
from afem.geometry.create import PlaneByPoints, TrimmedCurveByParameters
from afem.geometry.entities import Surface
from afem.geometry.project import (ProjectPointToCurve,
                                   ProjectPointToSurface)
from afem.graphics.display import ViewableItem
from afem.topology import ExploreShape
from afem.topology.bop import IntersectShapes
from afem.topology.check import CheckShape
from afem.topology.create import FaceBySurface, WiresByConnectedEdges
from afem.topology.distance import DistancePointToShapes
from afem.topology.entities import BBox
from afem.topology.modify import DivideC0Shape, DivideClosedShape
from afem.topology.transform import mirror_shape

__all__ = ["Body"]


class Body(ViewableItem):
    """
    Generic class for solid bodies and encapsulating necessary
    information when creating structural components.

    :param solid: The solid.
    :type solid: OCCT.TopoDS.TopoDS_Solid or OCCT.TopoDS.TopoDS_Shape
    :param str label: The label.

    :raise TypeError: If ``solid`` is not a ``TopoDS_Solid`` or cannot be
        converted to one.
    """

    def __init__(self, solid, label=None):
        super(Body, self).__init__()
        solid = CheckShape.to_solid(solid)
        if solid is None:
            msg = 'Invalid shape provided to Body. Requires a TopoDS_Solid.'
            raise TypeError(msg)
        self._solid = solid
        self._metadata = {}
        self._label = label
        self._sref = None
        self._sref_shape = None

    @property
    def label(self):
        """
        :return: The label.
        :rtype: str
        """
        return self._label

    @property
    def solid(self):
        """
        :return: The solid.
        :rtype: OCCT.TopoDS.TopoDS_Solid
        """
        return self._solid

    @property
    def outer_shell(self):
        """
        :return: The outer shell.
        :rtype: OCCT.TopoDS.TopoDS_Shell
        """
        return ExploreShape.outer_shell(self._solid)

    @property
    def metadata(self):
        """
        :return: The metadata dictionary.
        :rtype: dict
        """
        return self._metadata

    @property
    def sref(self):
        """
        :return: The reference surface.
        :rtype: afem.geometry.entities.Surface
        """
        return self._sref

    @property
    def sref_shape(self):
        """
        :return: The reference shape.
        :rtype: OCCT.TopoDS.TopoDS_Shape
        """
        return self._sref_shape

    @property
    def u1(self):
        """
        :return: Lower u-parameter of reference surface.
        :rtype: float
        """
        return self._sref.u1

    @property
    def u2(self):
        """
        :return: Upper u-parameter of reference surface.
        :rtype: float
        """
        return self._sref.u2

    @property
    def v1(self):
        """
        :return: Lower v-parameter of reference surface.
        :rtype: float
        """
        return self._sref.v1

    @property
    def v2(self):
        """
        :return: Upper v-parameter of reference surface.
        :rtype: float
        """
        return self._sref.v2

    def set_label(self, label):
        """
        Set name of body.

        :param str label: The name.

        :return: None.
        """
        self._label = label

    def add_metadata(self, key, value):
        """
        Add metadata to the body.

        :param key: The key.
        :param value: The value.

        :return: None.
        """
        self._metadata[key] = value

    def get_metadata(self, key):
        """
        Get metadata.

        :param key: The key.

        :return: The key or *None* if not present.
        :rtype: object or None

        :raise KeyError: If the key is not in the dictionary.
        """
        return self._metadata[key]

    def set_solid(self, solid):
        """
        Set the solid for the body.

        :param OCCT.TopoDS.TopoDS_Solid solid: The solid.

        :return: None.

        :raise TypeError: If *solid* is not a solid.
        """
        if not CheckShape.is_solid(solid):
            msg = 'Invalid shape provided to Body. Requires a TopoDS_Solid.'
            raise TypeError(msg)
        self._solid = solid

    def set_sref(self, srf, divide_closed=True, divide_c0=True):
        """
        Set the reference surface.

        :param afem.geometry.entities.Surface srf: The surface.
        :param bool divide_closed: Option to divide the surface if closed
            when creating the reference shape.
        :param bool divide_c0: Option to divide the surface at C0 boundaries
            when creating the reference shape.

        :return: None.

        :raise TypeError: If *srf* is not a supported type.
        """
        if not isinstance(srf, Surface):
            msg = 'Unsupported surface type.'
            raise TypeError(msg)

        self._sref = srf

        shape = FaceBySurface(srf).face
        if divide_closed:
            shape = DivideClosedShape(shape).shape
        if divide_c0:
            shape = DivideC0Shape(shape).shape
        self._sref_shape = shape

    def eval(self, u, v):
        """
        Evaluate a point on the reference surface.

        :param float u: Parameter in u-direction.
        :param float v: Parameter in v-direction.

        :return: Point on reference surface.
        :rtype: afem.geometry.entities.Point
        """
        return self._sref.eval(u, v)

    def norm(self, u, v):
        """
        Evaluate the surface normal of the reference surface.

        :param float u: Parameter in u-direction.
        :param float v: Parameter in v-direction.

        :return: Reference surface normal.
        :rtype: afem.geometry.entities.Vector
        """
        return self._sref.norm(u, v)

    def invert(self, p):
        """
        Find the parameters on the reference surface by inverting the point.

        :param point_like p: The point.

        :return: Parameters on the wing reference surface (u, v).
        :rtype: tuple(float)

        :raise RuntimeError: If no points are found in the projection
            algorithm.
        """
        proj = ProjectPointToSurface(p, self._sref)
        if not proj.success:
            msg = 'Failed to invert point.'
            raise RuntimeError(msg)
        return proj.nearest_param

    def extract_plane(self, u1, v1, u2, v2):
        """
        Extract a plane between parameters on the reference surface. The
        plane will be defined by three points. The first point is at (u1, v1),
        the second point is at (u2, v2), and the third point will be offset
        from the first point in the direction of the reference surface
        normal. The points should not be collinear.

        :param float u1: First u-parameter.
        :param float v1: First v-parameter.
        :param float u2: Second u-parameter.
        :param float v2: Second v-parameter.

        :return: The plane.
        :rtype: afem.geometry.entities.Plane
        """
        p1 = self.eval(u1, v1)
        p2 = self.eval(u2, v2)
        vn = self.norm(u1, v1)
        p3 = p1.copy()
        p3.translate(vn)
        return PlaneByPoints(p1, p2, p3).plane

    def extract_curve(self, u1, v1, u2, v2, basis_shape=None):
        """
        Extract a trimmed curve within the reference surface between the
        parameters.

        :param float u1: First u-parameter.
        :param float v1: First v-parameter.
        :param float u2: Second u-parameter.
        :param float v2: Second v-parameter.
        :param basis_shape: The shape that will be used to intersect with
            the reference shape. If not provided a plane will be
            created using the *extract_plane()* method. The parameters
            should create points that are on or very near the intersection
            between these two shapes. If they are not they will be projected to
            the intersection which could yield unanticipated results.
        :type basis_shape: afem.geometry.entities.Surface or
            OCCT.TopoDS.TopoDS_Shape

        :return: The curve.
        :rtype: afem.geometry.entities.TrimmedCurve

        :raise RuntimeError: If method fails.
        """
        p1 = self.eval(u1, v1)
        p2 = self.eval(u2, v2)

        if basis_shape is None:
            basis_shape = self.extract_plane(u1, v1, u2, v2)
        basis_shape = CheckShape.to_shape(basis_shape)

        bop = IntersectShapes(basis_shape, self.sref_shape, approximate=True)
        shape = bop.shape

        edges = ExploreShape.get_edges(shape)
        builder = WiresByConnectedEdges(edges)
        if builder.nwires == 0:
            msg = 'Failed to extract any curves.'
            raise RuntimeError(msg)

        if builder.nwires == 1:
            wire = builder.wires[0]
        else:
            dist = DistancePointToShapes(p1, builder.wires)
            wire = dist.nearest_shape
        crv = ExploreShape.curve_of_shape(wire)

        proj = ProjectPointToCurve(p1, crv)
        if not proj.success:
            msg = 'Failed to project point to reference curve.'
            raise RuntimeError(msg)
        u1c = proj.nearest_param

        proj = ProjectPointToCurve(p2, crv)
        if not proj.success:
            msg = 'Failed to project point to reference curve.'
            raise RuntimeError(msg)
        u2c = proj.nearest_param

        if u1c > u2c:
            crv.reverse()
            u1c, u2c = crv.reversed_u(u1c), crv.reversed_u(u2c)

        return TrimmedCurveByParameters(crv, u1c, u2c).curve

    def bbox(self, tol=None):
        """
        Return a bounding box of the body.

        :param tol: Optional tolerance to enlarge the bounding box.
        :type tol: float or None

        :return: Bounding box of the body.
        :rtype: afem.topology.entities.BBox
        """
        bbox = BBox()
        bbox.add_shape(self._solid)
        if tol is not None:
            bbox.enlarge(tol)
        return bbox

    def mirrored(self, pln, label=None):
        """
        Mirror this Body using the plane.

        :param afem.geometry.entities.Plane pln: The plane.
        :param str label: The label of the new Body.

        :return: Mirrored Body.
        :rtype: afem.oml.entities.Body
        """
        solid = CheckShape.to_solid(mirror_shape(self.solid, pln))
        body = Body(solid, label)
        if self.sref is not None:
            sref = self.sref.copy()
            sref.mirror(pln)
            body.set_sref(sref)
        return body

    def save(self, fn, binary=True):
        """
        Save this body.

        :param str fn: The filename.
        :param bool binary: If *True*, the document will be saved in a binary
            format. If *False*, the document will be saved in an XML format.

        :return: *True* if saved, *False* otherwise.
        :rtype: bool

        .. note::

            This method only saves the label, shape, reference surface, and
            color of the Body. Other user-defined metadata is currently not
            saved.
        """
        return self.save_bodies(fn, [self], binary)

    def load(self, fn):
        """
        Load a document containing a single body.

        :param str fn: The filename. The extension should be either ".xbf" for
            a binary file or ".xml" for an XML file.

        :return: The body.
        :rtype: afem.oml.entities.Body
        """
        bodies = self.load_bodies(fn)
        return list(bodies.values())[0]

    @staticmethod
    def save_bodies(fn, bodies, binary=True):
        """
        Save Body instances.

        :param str fn: The filename.
        :param collections.Sequence(afem.oml.entities.Body) bodies: The Body
            instances
        :param bool binary: If *True*, the document will be saved in a binary
            format. If *False*, the document will be saved in an XML format.

        :return: *True* if saved, *False* otherwise.
        :rtype: bool

        .. note::

            This method only saves the label, shape, reference surface, and
            color of the Body. Other user-defined metadata is currently not
            saved.
        """
        # Create document
        doc = XdeDocument(binary)

        # Add the bodies
        for body in bodies:
            label = doc.add_shape(body.solid, body.label, False)
            label.set_string('Body')
            label.set_color(body.color)

            # Sref
            if body.sref is not None:
                face = FaceBySurface(body.sref).face
                label = doc.add_shape(face, body.label, False)
                label.set_string('SREF')

        return doc.save_as(fn)

    @staticmethod
    def load_bodies(fn):
        """
        Load saved Body instances.

        :param str fn: The filename. The extension should be either ".xbf" for
            a binary file or ".xml" for an XML file.

        :return: A dictionary where the key is the name and the value is the
            Body instance.
        :rtype: dict(str, afem.oml.entities.Body)

        :raise TypeError: If the file extension type is not supported.
        """
        if fn.endswith('.xbf'):
            binary = True
        elif fn.endswith('.xml'):
            binary = False
        else:
            raise TypeError('Document type not supported.')

        # Open document
        doc = XdeDocument(binary)
        doc.open(fn)

        # Get the main label and iterate on top-level children which
        # should be parts
        label = doc.shapes_label
        body_data = []
        sref_to_body = {}
        for current in label.children_iter:
            name = current.name
            shape = current.shape
            type_ = current.string
            color = current.color

            if None in [name, shape, type_]:
                continue

            # Check for reference geometry
            if type_ == 'SREF':
                face = CheckShape.to_face(shape)
                sref_to_body[name] = ExploreShape.surface_of_face(face)
                continue

            # Add part data
            body_data.append((name, shape, color))

        # Create bodies
        label_to_bodies = {}
        for label, shape, color in body_data:
            sref = None
            if label in sref_to_body:
                sref = sref_to_body[label]
            solid = CheckShape.to_solid(shape)
            body = Body(solid, label)
            if sref is not None:
                body.set_sref(sref)
            if color is not None:
                body.color = color
            label_to_bodies[label] = body

        return label_to_bodies
