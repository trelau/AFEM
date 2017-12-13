#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2017 Laughlin Research, L.L.C.
#
# This file is subject to the license agreement that was delivered
# with this source code.
#
# THE SOFTWARE AND INFORMATION ARE PROVIDED ON AN AS "AS IS" BASIS,
# WITHOUT ANY WARRANTIES OR REPRESENTATIONS EXPRESS, IMPLIED OR
# STATUTORY; INCLUDING, WITHOUT LIMITATION, WARRANTIES OF QUALITY,
# PERFORMANCE, MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

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
