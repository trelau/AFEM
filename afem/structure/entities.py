#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2017 Laughlin Research, L.L.C.
#
# This file is subject to the license agreement that was delivered
# with this source code.
#
# THE SOFTWARE AND INFORMATION ARE PROVIDED ON AN "AS IS" BASIS,
# WITHOUT ANY WARRANTIES OR REPRESENTATIONS EXPRESS, IMPLIED OR 
# STATUTORY; INCLUDING, WITHOUT LIMITATION, WARRANTIES OF QUALITY,
# PERFORMANCE, MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.
from OCCT.TopAbs import TopAbs_SHAPE, TopAbs_EDGE, TopAbs_FACE
from OCCT.TopoDS import TopoDS_Shape
from numpy import mean

from afem.config import logger
from afem.geometry.check import CheckGeom
from afem.geometry.create import (PlaneByNormal, PlaneFromParameter,
                                  PointFromParameter, PointsAlongCurveByNumber,
                                  TrimmedCurveByCurve)
from afem.geometry.entities import Axis1, Plane, TrimmedCurve
from afem.geometry.project import (ProjectPointToCurve,
                                   ProjectPointToSurface)
from afem.graphics.display import ViewableItem
from afem.oml.entities import Body
from afem.structure.group import GroupAPI
from afem.topology.bop import (CutCylindricalHole, CutShapes, FuseShapes,
                               IntersectShapes, LocalSplit, SplitShapes)
from afem.topology.check import CheckShape, ClassifyPointInSolid
from afem.topology.create import (CompoundByShapes, HalfspaceBySurface,
                                  PointAlongShape, PointsAlongShapeByDistance,
                                  PointsAlongShapeByNumber,
                                  ShellByFaces, WiresByShape, FaceByPlane,
                                  SolidByDrag)
from afem.topology.distance import DistanceShapeToShape
from afem.topology.explore import ExploreShape
from afem.topology.fix import FixShape
from afem.topology.modify import (RebuildShapeByTool,
                                  RebuildShapeWithShapes, RebuildShapesByTool,
                                  SewShape, UnifyShape)
from afem.topology.props import LengthOfShapes, LinearProps, SurfaceProps

__all__ = ["Part", "CurvePart", "Beam1D", "SurfacePart", "WingPart", "Spar",
           "Rib", "FuselagePart", "Bulkhead", "Floor", "Frame", "Skin",
           "Stiffener1D", "Stiffener2D", "Stringer", "Beam2D",
           "shape_of_entity"]


def shape_of_entity(entity):
    """
    Get the shape of the entity. This method is useful if method inputs can
    either be a part or a shape. If the entity is already a shape it will be
    returned. If the entity is part the shape of the part will be returned. If
    the entity is a curve or surface then it will be converted to a shape.

    :param entity: The entity.
    :type entity: afem.geometry.entities.Geometry or
        afem.structure.entities.Part or afem.oml.entities.Body or
        OCCT.TopoDS.TopoDS_Shape

    :return: The shape.
    :rtype: OCCT.TopoDS.TopoDS_Shape
    """
    if isinstance(entity, Part):
        return entity.shape
    elif isinstance(entity, Body):
        return entity.solid
    else:
        return CheckShape.to_shape(entity)


class Part(ViewableItem):
    """
    Base class for all parts.

    :param str label: The label.
    :param OCCT.TopoDS.TopoDS_Shape shape: The shape.
    :param cref: The reference curve. If it is not a :class:`.TrimmedCurve`,
        then it will be converted to one.
    :type cref: afem.geometry.entities.Curve or None
    :param sref: The reference surface.
    :type sref: afem.geometry.entities.Surface or None
    :param group: The group to add the part to. If not provided the part will
        be added to the active group.
    :type group: str or afem.structure.group.Group or None

    :raise RuntimeError: If *shape* is not a valid shape or cannot be
            converted to a shape.
    """
    _indx = 1
    _mesh = None
    _core_shape = TopAbs_SHAPE

    def __init__(self, label, shape, cref=None, sref=None, group=None):
        super(Part, self).__init__()

        self._cref, self._sref = None, None
        self._metadata = {}
        self._subparts = {}

        self._id = Part._indx
        Part._indx += 1

        self._label = label
        self._shape = shape

        if cref is not None:
            self.set_cref(cref)

        if sref is not None:
            self.set_sref(sref)

        # Add to group
        GroupAPI.add_parts(group, self)

        # Log
        msg = ' '.join(['Creating part:', label])
        logger.info(msg)

    @property
    def type(self):
        """
        :return: The type of the part.
        :rtype: str
        """
        return self.__class__.__name__

    @property
    def label(self):
        """
        :return: The part label.
        :rtype: str
        """
        return self._label

    @property
    def id(self):
        """
        :return: The unique part ID.
        :rtype: int
        """
        return self._id

    @property
    def shape(self):
        """
        :return: The part shape.
        :rtype: OCCT.TopoDS.TopoDS_Shape
        """
        return self._shape

    @property
    def is_null(self):
        """
        :return: *True* if part shape is null, *False* if not.
        :rtype: bool
        """
        return self._shape.IsNull()

    @property
    def tol(self):
        """
        :return: The average tolerance of the part shape.
        :rtype: float
        """
        return ExploreShape.global_tolerance(self._shape, 0)

    @property
    def max_tol(self):
        """
        :return: The maximum tolerance of the part shape.
        :rtype: float
        """
        return ExploreShape.global_tolerance(self._shape, 1)

    @property
    def min_tol(self):
        """
        :return: The minimum tolerance of the part shape.
        :rtype: float
        """
        return ExploreShape.global_tolerance(self._shape, -1)

    @property
    def metadata(self):
        """
        :return: The part metadata.
        :rtype: dict
        """
        return self._metadata

    @property
    def subparts(self):
        """
        :return: List of sub-parts associated to this part.
        :rtype: list[afem.structure.entities.Part]
        """
        return self._subparts.values()

    @property
    def cref(self):
        """
        :return: The part reference curve.
        :rtype: afem.geometry.entities.TrimmedCurve
        """
        return self._cref

    @property
    def sref(self):
        """
        :return: The part reference surface.
        :rtype: afem.geometry.entities.Surface
        """
        return self._sref

    @property
    def has_cref(self):
        """
        :return: *True* if part has a reference curve, *False* if not.
        :rtype: bool
        """
        return CheckGeom.is_curve(self._cref)

    @property
    def has_sref(self):
        """
        :return: *True* if part has a reference surface, *False* if not.
        :rtype: bool
        """
        return CheckGeom.is_surface(self._sref)

    @property
    def is_planar(self):
        """
        :return: *True* if the reference surface is a plane, *False* if not.
        :rtype: bool
        """
        return isinstance(self._sref, Plane)

    @property
    def plane(self):
        """
        :return: A plane if the reference surface is a plane.
        :rtype: afem.geometry.entities.Plane

        :raise TypeError: If the reference surface is not a plane.
        """
        if not self.is_planar:
            raise TypeError('Reference surface is not a plane.')
        return self._sref

    @property
    def p1(self):
        """
        :return: The first point of the part reference curve.
        :rtype: afem.geometry.entities.Point

        :raises ValueError: If the part has no reference curve.
        """
        if not self.has_cref:
            msg = 'Part has no reference curve.'
            raise ValueError(msg)
        return self._cref.p1

    @property
    def p2(self):
        """
        :return: The last point of the part reference curve.
        :rtype: afem.geometry.entities.Point

        :raises ValueError: If the part has no reference curve.
        """
        if not self.has_cref:
            msg = 'Part has no reference curve.'
            raise ValueError(msg)
        return self._cref.p2

    @property
    def nedges(self):
        """
        :return: The number of edges in the part shape.
        :rtype: int
        """
        return len(self.edges)

    @property
    def edges(self):
        """
        :return: All the edges of the part shape.
        :rtype: list[OCCT.TopoDS.TopoDS_Edge]
        """
        return ExploreShape.get_edges(self._shape)

    @property
    def edge_compound(self):
        """
        :return: A compound containing the part edges.
        :rtype: OCCT.TopoDS.TopoDS_Compound
        """
        return CompoundByShapes(self.edges).compound

    @property
    def nfaces(self):
        """
        :return: The number of faces in the part shape.
        :rtype: int
        """
        return len(self.faces)

    @property
    def faces(self):
        """
        :return: All the faces of the part shape.
        :rtype: list[OCCT.TopoDS.TopoDS_Face]
        """
        return ExploreShape.get_faces(self._shape)

    @property
    def face_compound(self):
        """
        :return: A compound containing the part faces.
        :rtype: OCCT.TopoDS.TopoDS_Compound
        """
        return CompoundByShapes(self.faces).compound

    @property
    def mesh(self):
        """
        :return: The active mesh.
        :rtype: afem.smesh.meshes.Mesh
        """
        return self._mesh

    @property
    def submesh(self):
        """
        :return: Sub-mesh for the part shape using the active mesh.
        :rtype: afem.smesh.meshes.SubMesh
        """
        # Put core shape types into a compound
        if isinstance(self, CurvePart):
            subshape = self.edge_compound
        else:
            subshape = self.face_compound
        return self._mesh.get_submesh(subshape)

    @property
    def elements(self):
        """
        :return: The elements of the part.
        :rtype: list[afem.mesh.entities.Element]
        """
        ds = self.submesh.ds
        return [e for e in ds.elm_iter]

    @property
    def nodes(self):
        """
        :return: The nodes of part.
        :rtype: list[afem.mesh.entities.Node]
        """
        ds = self.submesh.ds
        return [n for n in ds.node_iter]

    def set_cref(self, cref):
        """
        Set the part reference curve.

        :param afem.geometry.entities.Curve cref: The curve. If it is not a
            :class:`.TrimmedCurve`, then it will be converted to one. Access
            the original curve using the *basis_curve* property (i.e.,
            part.cref.basis_curve).
        :param cref: The reference curve.

        :return: None.

        :raise TypeError: If *cref* is an invalid curve.
        """
        if not CheckGeom.is_curve(cref):
            raise TypeError('Invalid curve type.')

        if isinstance(cref, TrimmedCurve):
            self._cref = cref
        else:
            self._cref = TrimmedCurveByCurve(cref).curve

    def set_u1(self, u1):
        """
        Set the first parameter of the reference curve.

        :param float u1: The parameter.

        :return: None.

        :raise ValueError: If the *u1* is greater than or equal to *u2*.
        """
        if u1 >= self._cref.u2:
            msg = ('First parameter greater than or equal to second '
                   'parameter of curve.')
            raise ValueError(msg)

        self._cref.set_trim(u1, self._cref.u2)

    def set_u2(self, u2):
        """
        Set the last parameter of the reference curve.

        :param float u2: The parameter.

        :return: None.

        :raise ValueError: If the *u2* is less than or equal to *u1*.
        """
        if u2 <= self._cref.u1:
            msg = ('Second parameter less than or equal to first '
                   'parameter of curve.')
            raise ValueError(msg)

        self._cref.set_trim(self._cref.u1, u2)

    def set_p1(self, p1):
        """
        Set the first parameter of the reference curve by inverting the point.

        :param point_like p1: The point.

        :return: None.

        :raise RuntimeError: If no projection can be found.
        """
        u1 = self.invert_cref(p1)
        if u1 is None:
            msg = 'No projection found on reference curve.'
            raise RuntimeError(msg)

        self.set_u1(u1)

    def set_p2(self, p2):
        """
        Set the last parameter of the reference curve by inverting the point.

        :param point_like p2: The point.

        :return: None.

        :raise RuntimeError: If no projection can be found.
        """
        u2 = self.invert_cref(p2)
        if u2 is None:
            msg = 'No projection found on reference curve.'
            raise RuntimeError(msg)

        self.set_u2(u2)

    def trim_u1(self, entity):
        """
        Trim the first parameter of the reference curve by interesting it with
        the entity.

        :param entity: The entity.
        :type entity: afem.geometry.entities.Geometry or
            OCCT.TopoDS.TopoDS_Shape

        :return: None.

        :raise RuntimeError: If an intersection with the reference curve cannot
            be found.
        """
        shape1 = CheckShape.to_shape(self.cref)
        shape2 = CheckShape.to_shape(entity)
        bop = IntersectShapes(shape1, shape2)
        if not bop.is_done:
            raise RuntimeError('Failed to intersect reference curve.')

        pnts = [ExploreShape.pnt_of_vertex(v) for v in bop.vertices]
        p = CheckGeom.nearest_point(self.p1, pnts)
        self.set_p1(p)

    def trim_u2(self, entity):
        """
        Trim the last parameter of the reference curve by interesting it with
        the entity.

        :param entity: The entity.
        :type entity: afem.geometry.entities.Geometry or
            OCCT.TopoDS.TopoDS_Shape

        :return: None.

        :raise RuntimeError: If an intersection with the reference curve cannot
            be found.
        """
        shape1 = CheckShape.to_shape(self.cref)
        shape2 = CheckShape.to_shape(entity)
        bop = IntersectShapes(shape1, shape2)
        if not bop.is_done:
            raise RuntimeError('Failed to intersect reference curve.')

        pnts = [ExploreShape.pnt_of_vertex(v) for v in bop.vertices]
        p = CheckGeom.nearest_point(self.p2, pnts)
        self.set_p2(p)

    def set_sref(self, sref):
        """
        Set the part reference surface.

        :param afem.geometry.entities.Surface sref: The surface.

        :return: None.

        :raise TypeError: If *sref* is an invalid surface.
        """
        if not CheckGeom.is_surface(sref):
            msg = 'Invalid surface type.'
            raise TypeError(msg)
        self._sref = sref

    def nullify(self):
        """
        Destroy reference of the part shape.

        :return: None
        """
        self._shape.Nullify()

    def add_metadata(self, key, value):
        """
        Add metadata to the part.

        :param key: The key.
        :param value: The value.

        :return: None.
        """
        self._metadata[key] = value

    def get_metadata(self, key):
        """
        Get metadata.

        :param key: The key.

        :return: The metadata.

        :raise KeyError: If key not present in the dictionary.
        """
        return self._metadata[key]

    def add_subpart(self, key, subpart):
        """
        Add a sub-part to the part.

        :param str key: The key.
        :param afem.structure.entities.Part subpart: The sub-part.

        :return: None.

        :raise TypeError: If *subpart* is not a part.
        """
        if not isinstance(subpart, Part):
            msg = 'Sub-part is not a part.'
            raise TypeError(msg)
        self._subparts[key] = subpart

    def get_subpart(self, key):
        """
        Get a sub-part.

        :param str key: The key.

        :return: The sub-part. Returns *None* if the key is not present.
        :rtype: afem.structure.entities.Part

        :raise KeyError: If key not present in the dictionary.
        """
        return self._subparts[key]

    def set_shape(self, shape):
        """
        Set the shape of the part.

        :param OCCT.TopoDS.TopoDS_Shape shape: The shape.

        :return: None.

        :raise TypeError: If *shape* is not a valid shape or cannot be
            converted to a shape.
        """
        shape = CheckShape.to_shape(shape)
        if not shape:
            msg = 'Invalid shape.'
            raise TypeError(msg)

        self._shape = shape

    def local_to_global_u(self, u):
        """
        Convert local parameter from 0 <= u <= 1 to u1 <= u <= u2 using the
        part reference curve.

        :param float u: Local u-parameter.

        :return: Global u-parameter.
        :rtype: float

        :raise AttributeError: If part does not have a reference curve.
        """
        if not self.has_cref:
            msg = 'Part does not have a reference curve.'
            raise AttributeError(msg)
        return self._cref.local_to_global_param(u)

    def point_on_cref(self, u):
        """
        Evaluate point on reference curve.

        :param float u: The parameter.

        :return: The point.
        :rtype: afem.geometry.entities.Point

        :raise AttributeError: If part does not have a reference curve.
        """
        if not self.has_cref:
            msg = 'Part does not have a reference curve.'
            raise AttributeError(msg)
        return self._cref.eval(u)

    def point_on_sref(self, u, v):
        """
        Evaluate point on reference surface.

        :param float u: The u-parameter.
        :param float v: The v-parameter.

        :return: The point.
        :rtype: afem.geometry.entities.Point

        :raise AttributeError: If part does not have a reference surface.
        """
        if not self.has_sref:
            msg = 'Part does not have a reference surface.'
            raise AttributeError(msg)
        return self._sref.eval(u, v)

    def point_from_parameter(self, ds, u0=None, is_rel=False):
        """
        Evaluate point on reference curve at a distance from a parameter.

        :param float ds: The distance.
        :param float u0: The parameter. If not provided the first parameter
            of the reference curve will be used.
        :param bool is_rel: Option specifying if the distance is absolute or
            a relative to the length of the reference curve. If relative, then
            *ds* is multiplied by the curve length to get the absolute value
            for the :class:`.PointFromParameter` method.

        :return: The point.
        :rtype: afem.geometry.entities.Point

        :raise AttributeError: If part does not have a reference curve.
        """
        if not self.has_cref:
            msg = 'Part does not have a reference curve.'
            raise AttributeError(msg)

        if u0 is None:
            u0 = self._cref.u1

        if is_rel:
            ds *= self._cref.length

        return PointFromParameter(self._cref, u0, ds).point

    def points_by_number(self, n, d1=None, d2=None, shape1=None,
                         shape2=None, tol=1.0e-7):
        """
        Create a specified number of points along the reference curve.

        :param int n: Number of points to create (*n* > 0).
        :param float d1: An offset distance for the first point. This is
            typically a positive number indicating a distance from *u1*
            towards *u2*.
        :param float d2: An offset distance for the last point. This is
            typically a negative number indicating a distance from *u2*
            towards *u1*.
        :param OCCT.TopoDS.TopoDS_Shape shape1: A shape to define the first
            point. This shape is intersected with the edge or wire.
        :param OCCT.TopoDS.TopoDS_Shape shape2: A shape to define the last
            point. This shape is intersected with the edge or wire.
        :param float tol: Tolerance.

        :return: The points.
        :rtype: list[afem.geometry.entities.Point]

        :raise AttributeError: If part does not have a reference curve.
        """
        if not self.has_cref:
            msg = 'Part does not have a reference curve.'
            raise AttributeError(msg)

        edge = CheckShape.to_edge(self._cref)
        builder = PointsAlongShapeByNumber(edge, n, d1, d2, shape1, shape2,
                                           tol)
        return builder.points

    def points_by_distance(self, maxd, nmin=0, d1=None, d2=None, shape1=None,
                           shape2=None, tol=1.0e-7):
        """
        Create a points along the reference curve by distance.

        :param float maxd: The maximum allowed spacing between points. The
            actual spacing will be adjusted to not to exceed this value.
        :param int nmin: Minimum number of points to create.
        :param float d1: An offset distance for the first point. This is
            typically a positive number indicating a distance from *u1*
            towards *u2*.
        :param float d2: An offset distance for the last point. This is
            typically a negative number indicating a distance from *u2*
            towards *u1*.
        :param OCCT.TopoDS.TopoDS_Shape shape1: A shape to define the first
            point. This shape is intersected with the edge or wire.
        :param OCCT.TopoDS.TopoDS_Shape shape2: A shape to define the last
            point. This shape is intersected with the edge or wire.
        :param float tol: Tolerance.

        :return: The points.
        :rtype: list[afem.geometry.entities.Point]

        :raise AttributeError: If part does not have a reference curve.
        """
        if not self.has_cref:
            msg = 'Part does not have a reference curve.'
            raise AttributeError(msg)

        edge = CheckShape.to_edge(self._cref)
        builder = PointsAlongShapeByDistance(edge, maxd, d1, d2, shape1,
                                             shape2, nmin, tol)
        return builder.points

    def point_to_cref(self, pnt, direction=None):
        """
        Project a point to reference curve.

        :param afem.geometry.entities.Point pnt: The point. Position will be
            updated.
        :param vector_like direction: Projection direction.

        :return: *True* if projected, *False* if not.
        :rtype: bool

        :raise AttributeError: If part does not have a reference curve.
        """
        if not self.has_cref:
            msg = 'Part does not have a reference curve.'
            raise AttributeError(msg)

        proj = ProjectPointToCurve(pnt, self._cref, direction, update=True)
        if not proj.success:
            return False
        return True

    def points_to_cref(self, pnts, direction=None):
        """
        Project points to the reference curve.

        :param list[afem.geometry.entities.Point] pnts: The points. Position
            will be updated.
        :param vector_like direction: Projection direction.

        :return: List of status for each point.
        :rtype: list[bool]

        :raise AttributeError: If part does not have a reference curve.
        """
        if not self.has_cref:
            msg = 'Part does not have a reference curve.'
            raise AttributeError(msg)

        success = []
        for p in pnts:
            status = self.point_to_cref(p, direction)
            success.append(status)
        return success

    def point_to_sref(self, pnt, direction=None):
        """
        Project a point to reference surface.

        :param afem.geometry.entities.Point pnt: The point. Position will be
            updated.
        :param vector_like direction: Projection direction.

        :return: *True* if projected, *False* if not.
        :rtype: bool

        :raise AttributeError: If part does not have a reference surface.
        """
        if not self.has_sref:
            msg = 'Part does not have a reference surface.'
            raise AttributeError(msg)

        proj = ProjectPointToSurface(pnt, self._sref, direction)
        if not proj.success:
            return False

        p = proj.nearest_point
        pnt.set_xyz(p.xyz)
        return True

    def points_to_sref(self, pnts, direction=None):
        """
        Project points to reference surface.

        :param list[afem.geometry.entities.Point] pnts: The points. Position
            will be updated.
        :param vector_like direction: Projection direction.

        :return: List of status for each point.
        :rtype: list[bool]

        :raise AttributeError: If part does not have a reference surface.
        """
        if not self.has_sref:
            msg = 'Part does not have a reference surface.'
            raise AttributeError(msg)

        success = []
        for p in pnts:
            status = self.point_to_sref(p, direction)
            success.append(status)
        return success

    def plane_from_parameter(self, ds, u0=None, is_rel=False, ref_pln=None,
                             tol=1.0e-7):
        """
        Get a plane along the reference curve.

        :param float ds: The distance.
        :param float u0: The parameter. If not provided the first parameter
            of the reference curve will be used.
        :param bool is_rel: Option specifying if the distance is absolute or
            a relative to the length of the reference curve. If relative, then
            *ds* is multiplied by the curve length to get the absolute value
            for the :class:`.PlaneFromParameter` method.
        :param afem.geometry.entities.Plane ref_pln: The normal of this plane
            will be used to define the normal of the new plane. If no plane is
            provided, then the first derivative of the curve will define the
            plane normal.
        :param float tol: Tolerance.

        :return: The plane.
        :rtype: afem.geometry.entities.Plane

        :raise AttributeError: If part does not have a reference curve.
        """
        if not self.has_cref:
            msg = 'Part does not have a reference curve.'
            raise AttributeError(msg)

        if u0 is None:
            u0 = self._cref.u1

        if is_rel:
            ds *= self.cref.length

        return PlaneFromParameter(self._cref, u0, ds, ref_pln, tol).plane

    def invert_cref(self, pnt):
        """
        Invert the point on the reference curve.

        :param point_like pnt: The point.

        :return: The parameter on the reference curve. Returns *None* if no
            projection is found.
        :rtype: float or None

        :raise AttributeError: If part does not have a reference curve.
        """
        if not self.has_cref:
            msg = 'Part does not have a reference curve.'
            raise AttributeError(msg)

        proj = ProjectPointToCurve(pnt, self._cref)
        if proj.success:
            return proj.nearest_param
        return None

    def invert_sref(self, pnt):
        """
        Invert the point on the reference surface.

        :param point_like pnt: The point.

        :return: The parameters on the reference surface (u, v). Returns
            *None* if no projection is found.
        :rtype: tuple(float) or tuple(None)

        :raise AttributeError: If part does not have a reference surface.
        """
        if not self.has_sref:
            msg = 'Part does not have a reference surface.'
            raise AttributeError(msg)

        proj = ProjectPointToSurface(pnt, self._sref)
        if proj.success:
            return proj.nearest_param
        return None, None

    def distance(self, other):
        """
        Find the minimum distance between the part and other shape.

        :param other: Other part or shape.
        :type other: OCCT.TopoDS.TopoDS_Shape or afem.structure.entities.Part

        :return: The minimum distance.
        :rtype: float
        """
        other = shape_of_entity(other)
        return DistanceShapeToShape(self._shape, other).dmin

    def check(self, raise_error=True):
        """
        Check the shape of the part.

        :param bool raise_error: Option to raise an error if the shape is
            not valid.

        :return: *True* if shape is valid, *False* if not.
        :rtype: bool

        :raise RuntimeError: If the check fails and *raise_error* is ``True``.
        """
        check = CheckShape(self._shape).is_valid

        if not raise_error:
            return check

        if check:
            return True
        msg = ' '.join(['The shape of the part is not valid. Label:',
                        self.label])
        raise RuntimeError(msg)

    def fix(self, precision=None, min_tol=None, max_tol=None, context=None,
            include_subgroup=True):
        """
        Attempt to fix the shape of the part using :class:`.FixShape`.

        :param float precision: Basic precision value.
        :param float min_tol: Minimum tolerance.
        :param float max_tol: Maximum tolerance.
        :param context: The context shape or group.
        :type context: OCCT.TopoDS.TopoDS_Shape or
            afem.structure.entities.Group or str
        :param bool include_subgroup: Option to recursively include parts
            from any subgroups.

        :return: None.
        """
        if context is not None:
            if not isinstance(context, TopoDS_Shape):
                context = GroupAPI.as_compound(context, include_subgroup)

        new_shape = FixShape(self._shape, precision, min_tol, max_tol,
                             context).shape
        self.set_shape(new_shape)

    def cut(self, cutter):
        """
        Cut the part shape and rebuild this part.

        :param cutter: The cutter. If geometry is provided it will be
            converted to a shape before the Boolean operation.
        :type cutter: OCCT.TopoDS.TopoDS_Shape or afem.structure.entities.Part
            or afem.geometry.entities.Geometry

        :return: *True* if shape was cut, *False* if not.
        :rtype: bool
        """
        cutter = shape_of_entity(cutter)
        cut = CutShapes(self._shape, cutter)
        if not cut.is_done:
            return False

        self.rebuild(cut)

        return True

    def split(self, splitter, rebuild_both=True):
        """
        Split the part shape and rebuild this part. Optionally rebuild the
        splitter if it is a part. This method should handle splitting with
        parts and shapes or different types (i.e., splitting a face with an
        edge).

        :param splitter: The splitter.
        :type splitter: OCCT.TopoDS.TopoDS_Shape or afem.structure.entities.Part
        :param bool rebuild_both: Option to rebuild both if *splitter* is a
            part.

        :return: *True* if split, *False* if not.
        :rtype: bool

        :raise TypeError: If this part is or the splitter not a curve or
            surface part.
        """
        other_part = None
        if isinstance(splitter, Part):
            other_part = splitter
        splitter = shape_of_entity(splitter)

        split = SplitShapes()
        if rebuild_both:
            split.set_args([self._shape, splitter])
        else:
            split.set_args([self._shape])
            split.set_tools([splitter])
        split.build()
        if not split.is_done:
            return False

        parts = [self]
        if rebuild_both and other_part is not None:
            parts.append(other_part)
        for part in parts:
            part.rebuild(split)
        return True

    def rebuild(self, tool):
        """
        Rebuild the part shape with a supported tool.

        :param afem.topology.bop.BopCore tool: The tool.

        :return: *True* if modified, *False* if not.
        :rtype: bool

        :raise TypeError: If this part is not a curve or surface part.
        """
        if isinstance(self, (CurvePart, SurfacePart)):
            rebuild = RebuildShapeByTool(self._shape, tool)
        else:
            msg = 'Invalid part type in rebuild operation.'
            raise TypeError(msg)

        new_shape = rebuild.new_shape
        self.set_shape(new_shape)
        return True

    def discard_by_solid(self, solid, tol=None):
        """
        Discard shapes of the part using a solid. Any shapes of the part that
        have centroids inside the solid will be removed. Edges are checked
        for curve parts and faces are checked for surface parts.

        :param OCCT.TopoDS.TopoDS_Solid solid: The solid.
        :param float tol: The tolerance. If not provided then the part
            tolerance will be used.

        :return: *True* if shapes were discarded, *False* if not.
        :rtype: bool

        :raise TypeError: If this part is not a curve or surface part.
        """
        if isinstance(self, CurvePart):
            shapes = self.edges
        elif isinstance(self, SurfacePart):
            shapes = self.faces
        else:
            msg = 'Invalid part type in discard operation.'
            raise TypeError(msg)

        if tol is None:
            tol = self.tol

        rebuild = RebuildShapeWithShapes(self._shape)
        classifer = ClassifyPointInSolid(solid, tol=tol)

        modified = False
        for shape in shapes:
            if isinstance(self, CurvePart):
                cg = LinearProps(shape).cg
            else:
                cg = SurfaceProps(shape).cg

            classifer.perform(cg, tol)
            if classifer.is_in:
                rebuild.remove(shape)
                modified = True

        if not modified:
            return False

        new_shape = rebuild.apply()
        self.set_shape(new_shape)
        return True

    def discard_by_dmax(self, entity, dmax):
        """
        Discard shapes of the part using a shape and a distance. If the
        distance between a shape of the part and the given shape is greater
        than *dmax*, then the shape is removed. Edges are checked
        for curve parts and faces are checked for surface parts.

        :param entity: The shape.
        :type entity: OCCT.TopoDS.TopoDS_Shape or
            afem.geometry.entities.Geometry
        :param float dmax: The maximum distance.

        :return: *True* if shapes were discarded, *False* if not.
        :rtype: bool

        :raise TypeError: If this part is not a curve or surface part.
        """
        entity = CheckShape.to_shape(entity)

        if isinstance(self, CurvePart):
            shapes = self.edges
        elif isinstance(self, SurfacePart):
            shapes = self.faces
        else:
            msg = 'Invalid part type in discard operation.'
            raise TypeError(msg)

        rebuild = RebuildShapeWithShapes(self._shape)

        modified = False
        for part_shape in shapes:
            dmin = DistanceShapeToShape(entity, part_shape).dmin
            if dmin > dmax:
                rebuild.remove(part_shape)
                modified = True

        if not modified:
            return False

        new_shape = rebuild.apply()
        self.set_shape(new_shape)
        return True

    def discard_by_dmin(self, entity, dmin):
        """
        Discard shapes of the part using a shape and a distance. If the
        distance between a shape of the part and the given shape is less
        than *dmin*, then the shape is removed. Edges are checked
        for curve parts and faces are checked for surface parts.

        :param entity: The shape.
        :type entity: OCCT.TopoDS.TopoDS_Shape or
            afem.geometry.entities.Geometry
        :param float dmin: The minimum distance.

        :return: *True* if shapes were discarded, *False* if not.
        :rtype: bool

        :raise TypeError: If this part is not a curve or surface part.
        """
        entity = CheckShape.to_shape(entity)

        if isinstance(self, CurvePart):
            shapes = self.edges
        elif isinstance(self, SurfacePart):
            shapes = self.faces
        else:
            msg = 'Invalid part type in discard operation.'
            raise TypeError(msg)

        rebuild = RebuildShapeWithShapes(self._shape)

        modified = False
        for part_shape in shapes:
            dmin_ = DistanceShapeToShape(entity, part_shape).dmin
            if dmin > dmin_:
                rebuild.remove(part_shape)
                modified = True

        if not modified:
            return False

        new_shape = rebuild.apply()
        self.set_shape(new_shape)
        return True

    def discard_by_cref(self, size=None):
        """
        Discard shapes of the part by using the reference curve. An infinite
        solid is created at each end of the reference curve using the curve
        tangent. Any shape that has a centroid in these solids is removed.
        For a curve part edges are discarded, for a SurfacePart faces are
        discarded.

        :param float size: Option to define a finite solid box which might be
            more robust than an infinite solid.

        :return: *True* if shapes were discard, *False* if not.
        :rtype: bool

        :raise AttributeError: If part does not have a reference curve.
        """
        if not self.has_cref:
            msg = 'Part does not have a reference curve.'
            raise AttributeError(msg)

        # Create vectors at each end of the reference curve pointing "out" of
        # the part
        u1, u2 = self._cref.u1, self._cref.u2
        v1 = self._cref.deriv(u1, 1)
        v2 = self._cref.deriv(u2, 1)
        # Reverse v1 so it's "out" of the part
        v1.reverse()

        # Create planes at each end
        p1 = self._cref.eval(u1)
        p2 = self._cref.eval(u2)
        pln1 = PlaneByNormal(p1, v1).plane
        pln2 = PlaneByNormal(p2, v2).plane

        # Translate points to define half space
        if size is None:
            pref1 = p1 + 100. * v1.ijk
            pref2 = p2 + 100. * v2.ijk
            hs1 = HalfspaceBySurface(pln1, pref1).solid
            hs2 = HalfspaceBySurface(pln2, pref2).solid
        else:
            w, h = 2 * [size / 2.]
            f1 = FaceByPlane(pln1, -w, w, -h, h).face
            f2 = FaceByPlane(pln2, -w, w, -h, h).face
            v1.normalize()
            v2.normalize()
            v1.scale(size)
            v2.scale(size)
            hs1 = SolidByDrag(f1, v1).solid
            hs2 = SolidByDrag(f2, v2).solid

        # Discard by solid
        status1 = self.discard_by_solid(hs1)
        status2 = self.discard_by_solid(hs2)
        return status1 or status2

    def shared_vertices(self, other, as_compound=False):
        """
        Get vertices shared between the two parts.

        :param other: The other part or shape.
        :type other: afem.structure.entities.Part or OCCT.TopoDS.TopoDS_Shape
        :param bool as_compound: Option to return the shared vertices in a
            compound.

        :return: Shared vertices.
        :rtype: list[OCCT.TopoDS.TopoDS_Vertex] or OCCT.TopoDS.TopoDS_Compound
        """
        other = shape_of_entity(other)
        verts = ExploreShape.get_shared_vertices(self._shape, other)
        if not as_compound:
            return verts
        return CompoundByShapes(verts).compound

    def shared_edges(self, other, as_compound=False):
        """
        Get edges shared between the two parts.

        :param other: The other part or shape.
        :type other: afem.structure.entities.Part or OCCT.TopoDS.TopoDS_Shape
        :param bool as_compound: Option to return the shared edges in a
            compound.

        :return: Shared edges.
        :rtype: list[OCCT.TopoDS.TopoDS_Edge] or OCCT.TopoDS.TopoDS_Compound
        """
        other = shape_of_entity(other)
        edges = ExploreShape.get_shared_edges(self._shape, other)
        if not as_compound:
            return edges
        return CompoundByShapes(edges).compound

    @classmethod
    def reset(cls):
        """
        Reset the part index counter back to 1 and the mesh to None.

        :return: None.
        """
        cls._indx = 1
        cls._mesh = None

    @classmethod
    def set_mesh(cls, mesh):
        """
        Set the active mesh for all parts.

        :param afem.smesh.meshes.Mesh mesh: The mesh.

        :return: None.
        """
        cls._mesh = mesh


class CurvePart(Part):
    """
    Base class for curve parts.
    """
    _core_shape = TopAbs_EDGE

    @property
    def length(self):
        """
        :return: The length of all the edges of the part.
        :rtype: float
        """
        return LinearProps(self._shape).length


class Beam1D(CurvePart):
    """
    Beam 1-D.
    """
    pass


class SurfacePart(Part):
    """
    Base class for surface parts.
    """
    _core_shape = TopAbs_FACE

    @property
    def length(self):
        """
        :return: The length of the reference curve if available. Otherwise
            it returns the length of all edges of the part.
        :rtype: float
        """
        if self.has_cref:
            return self.cref.length
        return LinearProps(self._shape).length

    @property
    def area(self):
        """
        :return: The area of all faces of the part.
        :rtype: float
        """
        return SurfaceProps(self._shape).area

    @property
    def stiffeners(self):
        """
        :return: List of stiffeners associated to the part.
        :rtype: list[afem.structure.entities.Stiffener1D]
        """
        return [part for part in self.subparts]

    def make_shell(self):
        """
        Attempt to make a shell from the faces of the part. This method
        only constructs the faces into a shell without checking for a valid
        shape.

        :return: A new shell from the shape of the part.
        :rtype: OCCT.TopoDS.TopoDS_Shell
        """
        return ShellByFaces(self.faces).shell

    def fuse(self, *other_parts):
        """
        Fuse with other surface parts and rebuild both.

        :param afem.structure.entities.SurfacePart other_parts: The other
            part(s).

        :return: *True* if fused, *False* if not.
        :rtype: bool
        """
        # Putting the other parts in a compound avoids fusing them to each
        # other
        other_shapes = [part.shape for part in other_parts]
        other_compound = CompoundByShapes(other_shapes).compound

        fuse = FuseShapes(self._shape, other_compound)
        if not fuse.is_done:
            return False

        # Rebuild the part shapes
        parts = [self] + list(other_parts)
        shapes = [part.shape for part in parts]
        rebuild = RebuildShapesByTool(shapes, fuse)
        for part in parts:
            new_shape = rebuild.new_shape(part.shape)
            part.set_shape(new_shape)

        return True

    def sew(self, *other_parts):
        """
        Sew with other parts and rebuild all parts.

        :param afem.structure.entities.SurfacePart other_parts: The other
            part(s).

        :return: *True* if sewed, *False* if not.
        :rtype: bool
        """
        parts = [self] + list(other_parts)
        shapes = [self._shape] + [part.shape for part in other_parts]

        tol = float(mean([ExploreShape.global_tolerance(shape, 0) for shape in
                          shapes], dtype=float))
        max_tol = max(
            [ExploreShape.global_tolerance(shape, 1) for shape in shapes])

        sew = SewShape(tol=tol, max_tol=max_tol, cut_free_edges=True,
                       non_manifold=True)
        for part in parts:
            sew.add(part.shape)
        sew.perform()

        for part in parts:
            if not sew.is_modified(part.shape):
                continue
            mod_shape = sew.modified(part.shape)
            part.set_shape(mod_shape)
        return True

    def merge(self, other, unify=False):
        """
        Merge other surface part or shape with this one.

        :param other: The other part or shape.
        :type other: afem.structure.entities.SurfacePart or
            OCCT.TopoDS.TopoDS_Shape
        :param bool unify: Option to attempt to unify same domains.

        :return: *True* if merged, *False* if not.
        :rtype: bool
        """
        # Fuse the parts
        fuse = FuseShapes(self._shape, other)
        if not fuse.is_done:
            return False

        # Reset the shape
        self.set_shape(fuse.shape)

        # Unify if possible
        if not unify:
            return True
        return self.unify()

    def unify(self, edges=True, faces=True, bsplines=False):
        """
        Attempt to unify the same domains of the part shape.

        :param bool edges: Option to unify all possible edges.
        :param bool faces: Option to unify all possible faces.
        :param bool bsplines: Option to concatenate the curves of edges if they
            are C1 continuous.

        :return: *True* if unified, *False* if not.
        :rtype: bool
        """
        unify = UnifyShape(self._shape, edges, faces, bsplines)
        new_shape = unify.shape
        self.set_shape(new_shape)

    def split_local(self, subshape, tool):
        """
        Locally split the faces of he sub-shape in the context of the part
        shape.

        :param OCCT.TopoDS.TopoDS_Shape subshape: The sub-shape.
        :param tool: The tool to split the sub-shape with.
        :type tool: OCCT.TopoDS.TopoDS_Shape or afem.geometry.entities.Surface

        :return: *True* if split, *False* if not.
        :rtype: bool
        """
        bop = LocalSplit(subshape, tool, self._shape)
        if not bop.is_done:
            return False

        self.rebuild(bop)
        return True

    def shared_nodes(self, other):
        """
        Get nodes shared between the two parts.

        :param afem.structure.entities.SurfacePart other: The other part.

        :return: Shared nodes.
        :rtype: list[afem.mesh.entities.Node]
        """
        nodes1 = set(self.nodes)
        nodes2 = set(other.nodes)
        return list(nodes1 & nodes2)

    def cut_hole(self, d, ds, u0=None, is_rel=False):
        """
        Cut a hole in the part and update its shape.

        :param float d: Diameter of hole.
        :param float ds: The distance.
        :param float u0: The parameter. If not provided the first parameter
            of the reference curve will be used.
        :param bool is_rel: Option specifying if the distance is absolute or
            a relative to the length of the reference curve. If relative, then
            *ds* is multiplied by the curve length to get the absolute value.

        :return: The status of the Boolean operation that will cut the hole.
            *True* if the operation was performed, *False* if not.
        :rtype: bool

        :raise AttributeError: If the part does not have a reference curve or
            surface.
        """
        if not self.has_cref:
            msg = 'Part does not have a reference curve.'
            raise AttributeError(msg)

        if not self.has_sref:
            msg = 'Part does not have a reference surface.'
            raise AttributeError(msg)

        # Intersect the shape with a plane to use for the hole height location
        pln = self.plane_from_parameter(ds, u0, is_rel)
        bop = IntersectShapes(self._shape, pln)
        wires = WiresByShape(bop.shape).wires
        los = LengthOfShapes(wires)
        max_length = los.max_length
        wire = los.longest_shape
        mid_length = max_length / 2.
        p = PointAlongShape(wire, mid_length).point
        u, v = self.invert_sref(p)
        v = self.sref.norm(u, v)
        v = CheckGeom.to_direction(v)
        ax1 = Axis1(p, v)

        # Cut the hole
        r = d / 2.
        bop = CutCylindricalHole(self._shape, r, ax1)
        self.set_shape(bop.shape)

        return bop.is_done

    def cut_holes(self, n, d):
        """
        Cut holes along the reference curve at evenly spaced intervals
        (experimental).

        :param int n: The number of holes.
        :param float d: The diameter.

        :return: None.

        :raise AttributeError: If the part does not have a reference curve or
            surface.
        """
        if not self.has_cref:
            msg = 'Part does not have a reference curve.'
            raise AttributeError(msg)

        pac = PointsAlongCurveByNumber(self.cref, n + 2)
        for u in pac.parameters[1:-1]:
            self.cut_hole(d, 0., u)


class WingPart(SurfacePart):
    """
    Base class for wing parts.
    """
    pass


class Spar(WingPart):
    """
    Wing spar.
    """
    pass


class Rib(WingPart):
    """
    Wing rib.
    """
    pass


class FuselagePart(SurfacePart):
    """
    Base class for fuselage parts.
    """
    pass


class Bulkhead(FuselagePart):
    """
    Bulkhead.
    """
    pass


class Floor(FuselagePart):
    """
    Floor.
    """
    pass


class Frame(FuselagePart):
    """
    Frame.
    """
    pass


class Skin(SurfacePart):
    """
    Skin.
    """
    pass


class Stiffener1D(CurvePart):
    """
    1-D stiffener for surface parts.
    """
    pass


class Stiffener2D(SurfacePart):
    """
    2-D stiffener for surface parts.
    """
    pass


class Stringer(SurfacePart):
    """
    Stringer.
    """
    pass


class Beam2D(SurfacePart):
    """
    Beam 2-D.
    """
    pass
