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

from OCCT.BRep import BRep_Builder
from OCCT.BRepTools import BRepTools
from OCCT.TopoDS import TopoDS_Shape


def write_brep(shape, fn):
    """
    Write a BREP file using the shape.

    :param OCCT.TopoDS.TopoDS_Shape shape: The shape.
    :param str fn: The filename.

    :return: None.
    """
    BRepTools.Write_(shape, fn)


def read_brep(fn):
    """
    Read a BREP file and return the shape.

    :param str fn: The filename.

    :return: The shape.
    :rtype: OCCT.TopoDS.TopoDS_Shape
    """
    shape = TopoDS_Shape()
    builder = BRep_Builder()
    BRepTools.Read_(shape, fn, builder)

    return shape
