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
import os

from OCCT.BRepBuilderAPI import (BRepBuilderAPI_MakeVertex,
                                 BRepBuilderAPI_MakeEdge,
                                 BRepBuilderAPI_MakeFace)
from OCCT.Quantity import (Quantity_TOC_RGB, Quantity_Color)
from OCCT.SMESH import SMESH_Mesh, SMESH_subMesh
from OCCT.TopoDS import TopoDS_Shape
from OCCT.Visualization import BasicViewer
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
        return self.display_shape(body.solid.object, body.color, body.transparency)

    def display_part(self, part):
        """
        Display a part.

        :param afem.structure.entities.Part part: The part.

        :return: The AIS_Shape created for the part.
        :rtype: OCCT.AIS.AIS_Shape
        """
        return self.display_shape(part.shape.object, part.color, part.transparency)

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
        from afem.topology.entities import Shape

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
            elif isinstance(item, Shape):
                self.display_shape(item.object, item.color, item.transparency)
