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

from OCCT.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCCT.gce import gce_MakeMirror


def mirror_shape(shape, pln):
    """
    Mirror a shape about a plane.

    :param OCCT.TopoDS.TopoDS_Shape shape: The shape.
    :param afem.geometry.entities.Plane pln: The plane.

    :return: The mirrored shape.
    :rtype: OCCT.TopoDS.TopoDS_Shape

    :raise RuntimeError: If the transformation fails or is not done.
    """
    trsf = gce_MakeMirror(pln.object.Pln()).Value()
    builder = BRepBuilderAPI_Transform(shape, trsf, True)
    if not builder.IsDone():
        raise RuntimeError('Failed to mirror the shape. The builder is not done.')
    return builder.Shape()
