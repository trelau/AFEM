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
from warnings import warn

from OCCT.Quantity import Quantity_TOC_RGB, Quantity_Color
from numpy.random import rand

__all__ = ["ViewableItem", "ShapeHolder"]


class ViewableItem(object):
    """
    Base class for types that can be viewed.

    :ivar OCCT.Quantity.Quantity_Color color: The OpenCASCADE color quantity.
    :ivar float transparency: The transparency level.
    """

    def __init__(self):
        r, g, b = rand(1, 3)[0]
        self.color = Quantity_Color(r, g, b, Quantity_TOC_RGB)
        self.transparency = 0.

    @property
    def displayed_shape(self):
        """
        :return: The shape to be displayed.
        :rtype: OCCT.TopoDS.TopoDS_Shape
        """
        msg = ('Viewing this type is not yet support. '
               'Implement displayed_shape property.')
        raise NotImplementedError(msg)

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


class ShapeHolder(ViewableItem):
    """
    Base class for types that store a shape.

    :param Type[afem.topology.entities.Shape] expected_type: The expected type.
    :param afem.topology.entities.Shape shape: The shape.
    """

    def __init__(self, expected_type, shape=None):
        super(ShapeHolder, self).__init__()
        self._type = expected_type
        self._type_name = expected_type.__name__
        self._shape = None
        if shape is not None:
            self.set_shape(shape)

    @property
    def shape(self):
        """
        :Getter: The shape.
        :Setter: Set the shape.
        :type: afem.topology.entities.Shape
        """
        return self._shape

    @shape.setter
    def shape(self, shape):
        self.set_shape(shape)

    @property
    def displayed_shape(self):
        """
        :return: The shape to be displayed.
        :rtype: OCCT.TopoDS.TopoDS_Shape
        """
        return self.shape.object

    def set_shape(self, shape):
        """
        Set the shape.

        :param afem.topology.entities.Shape shape: The shape.

        :return: None.
        """
        if not isinstance(shape, self._type):
            this = self.__class__.__name__
            other = shape.__class__.__name__
            msg = ('Invalid shape provided for a {} object. '
                   'Expected a {} but got a {}.'.format(this, self._type_name,
                                                        other))
            warn(msg, RuntimeWarning)

        self._shape = shape
