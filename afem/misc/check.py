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
try:
    from itertools import izip as _zip
except ImportError:
    _zip = zip

from collections import Sequence
from itertools import tee

from numpy import ndarray

__all__ = ["is_array_like", "is_array_type", "is_local_domain", "pairwise"]


def is_array_like(obj):
    """
    Test to see if an object is array_like (tuple, list, or ndarray).

    :param obj: Object to test.

    :return: *True* if array_like, *False* if not.
    :rtype: bool
    """
    return isinstance(obj, (Sequence, ndarray))


def is_array_type(rtype):
    """
    Test to see if return type parameter is a NumPy array.

    :param str rtype: Return type parameter.

    :return: *True* if return type parameter is a NumPy array, *False* if not.
    :rtype: bool
    """
    if rtype.lower() in ['ndarray', 'array', 'arr', 'np', 'a']:
        return True
    return False


def is_local_domain(domain):
    """
    Test to see if domain parameter is local or global.

    :param str domain: Domain parameter.

    :return: *True* if domain parameter is local, *False* if is global.
    :rtype: bool
    """
    if domain.lower() in ['l', 'local', 'loc']:
        return True
    return False


def pairwise(iterable):
    """
    s -> (s0,s1), (s1,s2), (s2, s3), ...
    """
    a, b = tee(iterable)
    next(b, None)
    return _zip(a, b)
