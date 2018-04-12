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
from collections import Sequence

from OCCT.Quantity import Quantity_TOC_RGB, Quantity_Color
from numpy.random import rand

from afem.config import logger

__all__ = ["Metadata", "NamedItem", "ViewableItem", "ShapeHolder"]


class Metadata(dict):
    """
    Simple class for storing data as (key, value) pairs.
    """

    def set(self, key, value):
        """
        Set data.

        :param key: The key.
        :param value: The value.

        :return: None.
        """
        self[key] = value


class NamedItem(object):
    """
    Base class for types that can be give a name.
    """

    def __init__(self, name='Item'):
        self._name = name
        self._metadata = Metadata()

    @property
    def name(self):
        """
        :Getter: The name.
        :Setter: Set name.
        :type: str
        """
        return self._name

    @name.setter
    def name(self, name):
        self.set_name(name)

    @property
    def metadata(self):
        """
        :return: The metadata dictionary.
        :rtype: afem.base.entities.Metadata
        """
        return self._metadata

    def set_name(self, name):
        """
        Set name.

        :param str name: The name.

        :return: None.
        """
        self._name = name


class ViewableItem(object):
    """
    Base class for types that can be viewed.
    """

    def __init__(self):
        self._color = None
        self._transparency = 0.

    @property
    def displayed_shape(self):
        """
        :return: The shape to be displayed.
        :rtype: OCCT.TopoDS.TopoDS_Shape
        """
        msg = ('Viewing this type is not yet support. '
               'Implement displayed_shape property.')
        raise NotImplementedError(msg)

    @property
    def color(self):
        """
        :return: The color or *None* if not set.
        :rtype: OCCT.Quantity.Quantity_Color or None
        """
        return self._color

    @property
    def transparency(self):
        """
        :return: The transparency.
        :rtype: float
        """
        return self._transparency

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
        self._color = Quantity_Color(r, g, b, Quantity_TOC_RGB)

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
        self._transparency = transparency

    def random_color(self):
        """
        Set the color randomly.

        :return: None.
        """
        r, g, b = rand(1, 3)[0]
        self._color = Quantity_Color(r, g, b, Quantity_TOC_RGB)


class ShapeHolder(ViewableItem):
    """
    Base class for types that store a shape.

    :param expected_type: The expected type(s).
    :type expected_type: Type[afem.topology.entities.Shape] or
        tuple(Type[afem.topology.entities.Shape])
    :param afem.topology.entities.Shape shape: The shape.
    """

    def __init__(self, expected_type, shape=None):
        super(ShapeHolder, self).__init__()
        if isinstance(expected_type, Sequence):
            self._types = expected_type
        else:
            self._types = (expected_type,)
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
        if not isinstance(shape, self._types):
            this = self.__class__.__name__
            other = shape.__class__.__name__

            expected = []
            for type_ in self._types:
                expected.append(type_.__name__)
            if len(expected) == 1:
                expected = 'a ' + expected[0]
            else:
                expected = 'one of (' + ', '.join(expected) + ')'
            name = 'Unknown'
            if isinstance(self, NamedItem):
                name = self.name
            msg = ('Invalid shape provided for a {} object with name {}. '
                   'Got a {} but expected {}.'.format(this, name, other,
                                                      expected))
            logger.warning(msg)

        self._shape = shape
