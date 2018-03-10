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

import os

from OCCT.BRepBuilderAPI import (BRepBuilderAPI_Transform,
                                 BRepBuilderAPI_MakeVertex,
                                 BRepBuilderAPI_MakeEdge,
                                 BRepBuilderAPI_MakeFace)
from OCCT.Quantity import (Quantity_TOC_RGB, Quantity_Color)
from OCCT.SMESH import SMESH_Mesh, SMESH_subMesh
from OCCT.TopoDS import TopoDS_Shape
from OCCT.Visualization import BasicViewer
from OCCT.gce import gce_MakeMirror
from numpy.random import rand

__all__ = ["ViewableItem", "Viewer"]

# Icon location
_icon = os.path.dirname(__file__) + '/resources/main.png'


class ViewableItem(object):
    """
    Base class for items that can be viewed.

    :var OCCT.Quantity.Quantity_Color color: The OCC color quantity. The
        color is set randomly during initialization.
    :var float transparency: The transparency level.
    :var afem.geometry.entities.Plane plane: The plane to mirror the object
        about. If provided then object will be mirrored about the plane for
        visualization purposes only.
    """

    def __init__(self):
        r, g, b = rand(1, 3)[0]
        self.color = Quantity_Color(r, g, b, Quantity_TOC_RGB)
        self.transparency = 0.
        self.mirror_plane = None

    def set_color(self, r, g, b):
        """
        Set color (0. <= r, g, b <= 1.).

        :param float r: Red.
        :param float g: Green.
        :param float b: Blue.

        :return: None.
        """
        if r > 1.:
            r /= 255.
        if g > 1.:
            g /= 255.
        if b > 1.:
            b /= 255.
        self.color = Quantity_Color(r, g, b, Quantity_TOC_RGB)

    def set_transparency(self, transparency):
        """
        Set the opacity for graphics.

        :param float transparency: Level of transparency (0 to 1).

        :return: None.
        """
        if transparency < 0.:
            transparency = 0.
        elif transparency > 1.:
            transparency = 1.
        self.transparency = transparency

    def set_mirror(self, pln):
        """
        Set a plane to mirror the item.

        :param afem.geometry.entities.Plane pln: The plane.

        :return: None.
        """
        self.mirror_plane = pln

    def get_mirrored(self):
        """
        Get the mirrored shape.

        :return: The mirrored shape.
        :rtype: OCCT.TopoDS.TopoDS_Shape
        """
        if not self.mirror_plane:
            return None

        trsf = gce_MakeMirror(self.mirror_plane.object.Pln()).Value()
        builder = BRepBuilderAPI_Transform(self, trsf, True)
        if not builder.IsDone():
            return None
        return builder.Shape()


class Viewer(BasicViewer):
    """
    Simple tool for viewing entities.
    """

    def __init__(self, width=800, height=600):
        super(Viewer, self).__init__(width, height)
        self.SetLabel('AFEM')

    def display_geom(self, geom, *args, **kwargs):
        """
        Display a geometric entity.

        :param afem.geometry.entities.Geometry geom: The geometry. Must be
            either a Point, Curve, or Surface.

        :return: The AIS_Shape created for the geometry. Returns *None* if the
            entity cannot be converted to a shape.
        :rtype: OCCT.AIS.AIS_Shape or None
        """
        from afem.geometry.entities import Point, Curve, Surface

        if isinstance(geom, Point):
            shape = BRepBuilderAPI_MakeVertex(geom).Vertex()
        elif isinstance(geom, Curve):
            shape = BRepBuilderAPI_MakeEdge(geom.object).Edge()
        elif isinstance(geom, Surface):
            shape = BRepBuilderAPI_MakeFace(geom.object, 1.0e-7).Face()
        else:
            return None

        return self.display_shape(shape, geom.color, geom.transparency)

    def display_body(self, body):
        """
        Display a body.

        :param afem.oml.entities.Body body: The body.

        :return: The AIS_Shape created for the body.
        :rtype: OCCT.AIS.AIS_Shape
        """
        return self.display_shape(body.solid, body.color, body.transparency)

    def display_part(self, part):
        """
        Display a part.

        :param afem.structure.entities.Part part: The part.

        :return: The AIS_Shape created for the part.
        :rtype: OCCT.AIS.AIS_Shape
        """
        return self.display_shape(part.shape, part.color, part.transparency)

    def display_parts(self, parts):
        """
        Display a sequence of parts.

        :param collections.Sequence[afem.structure.entities.Part] parts: The
            parts.

        :return: None.
        """
        for part in parts:
            self.display_part(part)

    def display_groups(self, group, include_subgroup=True):
        """
        Display all parts of a group.

        :param afem.structure.group.Group group: The group.
        :param bool include_subgroup: Option to recursively include parts
            from any subgroups.

        :return: None.
        """
        for part in group.get_parts(include_subgroup):
            self.display_part(part)

    def add(self, *items):
        """
        Add items to be displayed.

        :param items: The items.
        :type items: OCCT.TopoDS.TopoDS_Shape or
            afem.geometry.entities.Geometry or
            afem.oml.entities.Body or
            afem.structure.entities.Part or
            afem.structure.group.Group or
            OCCT.SMESH.SMESH_Mesh or
            OCCT.SMESH.SMESH_subMesh or
            afem.smesh.meshes.Mesh or
            afem.smesh.meshes.SubMesh

        :return: None.
        """
        from afem.geometry.entities import Geometry
        from afem.oml.entities import Body
        from afem.structure.group import Group
        from afem.structure.entities import Part
        from afem.smesh.meshes import Mesh, SubMesh

        for item in items:
            if isinstance(item, Part):
                self.display_part(item)
            elif isinstance(item, Group):
                self.display_groups(item)
            elif isinstance(item, Geometry):
                self.display_geom(item)
            elif isinstance(item, Body):
                self.display_body(item)
            elif isinstance(item, TopoDS_Shape):
                self.display_shape(item)
            elif isinstance(item, (Mesh, SubMesh)):
                self.display_mesh(item.object)
            elif isinstance(item, (SMESH_Mesh, SMESH_subMesh)):
                self.display_mesh(item)
