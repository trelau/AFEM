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
from OCCT.TColStd import (TColStd_Array1OfInteger, TColStd_Array1OfReal,
                          TColStd_Array2OfReal)
from OCCT.TColgp import (TColgp_Array1OfPnt, TColgp_Array1OfPnt2d,
                         TColgp_Array2OfPnt, TColgp_HArray1OfPnt,
                         TColgp_HArray1OfPnt2d)
from OCCT.TopAbs import (TopAbs_COMPOUND, TopAbs_COMPSOLID, TopAbs_EDGE,
                         TopAbs_FACE, TopAbs_SHELL, TopAbs_SOLID,
                         TopAbs_VERTEX,
                         TopAbs_WIRE)
from OCCT.TopTools import TopTools_ListOfShape
from OCCT.TopoDS import TopoDS
from OCCT.gp import gp_Pnt, gp_Pnt2d
from numpy import array as np_array, zeros

from afem.misc.check import is_array_like


def to_gp_pnt(p):
    """
    Convert the point_like entity to a gp_Pnt.
    """
    if isinstance(p, gp_Pnt):
        return p
    if is_array_like(p) and len(p) == 3:
        return gp_Pnt(*p)
    return None


def to_gp_pnt2d(p):
    """
    Convert the point_like entity to a gp_Pnt2d.
    """
    if isinstance(p, gp_Pnt2d):
        return p
    if is_array_like(p) and len(p) == 2:
        return gp_Pnt2d(*p)
    return None


def to_tcolgp_array1_pnt(pnts):
    """
    Convert the 1-D array of point_like entities to OCC data.

    :param array_like pnts: Array of points to convert.

    :return: OCC array of points.
    :rtype: TColgp_Array1OfPnt
    """
    gp_pnts = []
    for gp in pnts:
        gp = to_gp_pnt(gp)
        if not gp:
            continue
        gp_pnts.append(gp)

    n = len(gp_pnts)
    array = TColgp_Array1OfPnt(1, n)
    for i, gp in enumerate(gp_pnts, 1):
        array.SetValue(i, gp)

    return array


def to_tcolgp_array1_pnt2d(pnts):
    """
    Convert the 1-D array of point_like entities to OCC data.

    :param array_like pnts: Array of points to convert.

    :return: OCC array of points.
    :rtype: TColgp_Array1OfPnt2d
    """
    gp_pnts = []
    for gp in pnts:
        gp = to_gp_pnt2d(gp)
        if not gp:
            continue
        gp_pnts.append(gp)

    n = len(gp_pnts)
    array = TColgp_Array1OfPnt2d(1, n)
    for i, gp in enumerate(gp_pnts, 1):
        array.SetValue(i, gp)

    return array


def to_tcolgp_harray1_pnt(pnts):
    """
    Convert the 1-D array of point_like entities to OCC data.

    :param array_like pnts: Array of points to convert.

    :return: OCC array of points.
    :rtype: TColgp_HArray1OfPnt
    """
    gp_pnts = []
    for gp in pnts:
        gp = to_gp_pnt(gp)
        if not gp:
            continue
        gp_pnts.append(gp)

    n = len(gp_pnts)
    harray = TColgp_HArray1OfPnt(1, n)
    for i, gp in enumerate(gp_pnts, 1):
        harray.SetValue(i, gp)

    return harray


def to_tcolgp_harray1_pnt2d(pnts):
    """
    Convert the 1-D array of point2d_like entities to OCC data.

    :param point2d_like pnts: Array of 2-D points to convert.

    :return: OCC array of points.
    :rtype: TColgp_HArray1OfPnt2d
    """
    gp_pnts = []
    for gp in pnts:
        gp = to_gp_pnt2d(gp)
        if not gp:
            continue
        gp_pnts.append(gp)

    n = len(gp_pnts)
    harray = TColgp_HArray1OfPnt2d(1, n)
    for i, gp in enumerate(gp_pnts, 1):
        harray.SetValue(i, gp)

    return harray


def to_tcolstd_array1_real(array):
    """
    Convert the 1-D array of floats to OCC data.

    :param array_like array: Array of floats.

    :return: OCC array of floats.
    :rtype: TColStd_Array1OfReal
    """
    flts = [float(x) for x in array]
    n = len(flts)
    array = TColStd_Array1OfReal(1, n)
    for i, x in enumerate(flts, 1):
        array.SetValue(i, x)

    return array


def to_tcolstd_array1_integer(array):
    """
    Convert the 1-D array of integers to OCC data.

    :param array_like array: Array of integers.

    :return: OCC array of integers.
    :rtype: TColStd_Array1OfInteger
    """
    ints = [int(x) for x in array]
    n = len(ints)
    array = TColStd_Array1OfInteger(1, n)
    for i, x in enumerate(ints, 1):
        array.SetValue(i, x)

    return array


def to_tcolgp_array2_pnt(pnts):
    """
    Convert the 2-D array of point_like entities to OCC data.

    :param array_like pnts: Array of points to convert.

    :return: OCC array of points.
    :rtype: TColgp_Array2OfPnt
    """
    pnts = np_array(pnts, dtype=float)
    n, m = pnts.shape[0:2]
    gp_pnts = []
    for i in range(0, n):
        row = []
        for j in range(0, m):
            gp = to_gp_pnt(pnts[i, j])
            row.append(gp)
        gp_pnts.append(row)

    array = TColgp_Array2OfPnt(1, n, 1, m)
    for i, row in enumerate(gp_pnts, 1):
        for j, gp in enumerate(row, 1):
            array.SetValue(i, j, gp)

    return array


def to_tcolstd_array2_real(array):
    """
    Convert the 2-D array of floats to OCC data.

    :param array_like array: Array of floats.

    :return: OCC array of floats.
    :rtype: TColStd_Array2OfReal
    """
    flts = np_array(array, dtype=float)
    n, m = flts.shape
    array = TColStd_Array2OfReal(1, n, 1, m)
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            x = float(flts[i - 1, j - 1])
            array.SetValue(i, j, x)

    return array


def to_np_from_tcolstd_array1_real(tcol_array):
    """
    Convert OCC data to NumPy array.

    :param tcol_array: OCC array of floats.
    :type tcol_array: TColStd_Array1OfReal

    :return: NumPy array of floats.
    :rtype: ndarray
    """
    n = tcol_array.Length()
    array = zeros(n, dtype=float)
    for i in range(n):
        x = tcol_array.Value(i + 1)
        array[i] = x
    return array


def to_np_from_tcolstd_array1_integer(tcol_array):
    """
    Convert OCC data to NumPy array.

    :param tcol_array: OCC array of integers.
    :type tcol_array: TColStd_Array1OfInteger

    :return: NumPy array of integers.
    :rtype: ndarray
    """
    n = tcol_array.Length()
    array = zeros(n, dtype=int)
    for i in range(n):
        x = tcol_array.Value(i + 1)
        array[i] = x
    return array


def to_np_from_tcolgp_array1_pnt(tcol_array):
    """
    Convert OCC data to NumPy array.

    :param tcol_array: OCC array of points.
    :type tcol_array: TColgp_Array1OfPnt

    :return: NumPy array of points.
    :rtype: ndarray
    """
    n = tcol_array.Length()
    array = zeros((n, 3), dtype=float)
    for i in range(n):
        p = tcol_array.Value(i + 1)
        array[i, :] = p.X(), p.Y(), p.Z()
    return array


def to_np_from_tcolgp_array2_pnt(tcol_array):
    """
    Convert OCC data to NumPy array.

    :param tcol_array: OCC array of points.
    :type tcol_array: TColgp_Array2OfPnt

    :return: NumPy array of points.
    :rtype: ndarray
    """
    n, m = tcol_array.ColLength(), tcol_array.RowLength()
    array = zeros((n, m, 3), dtype=float)
    for i in range(n):
        for j in range(m):
            p = tcol_array.Value(i + 1, j + 1)
            array[i, j, :] = p.X(), p.Y(), p.Z()
    return array


def to_np_from_tcolstd_array2_real(tcol_array):
    """
    Convert OCC data to NumPy array.

    :param tcol_array: OCC array of floats.
    :type tcol_array: TColStd_Array2OfReal

    :return: NumPy array of floats.
    :rtype: ndarray
    """
    n, m = tcol_array.ColLength(), tcol_array.RowLength()
    array = zeros((n, m), dtype=float)
    for i in range(n):
        for j in range(m):
            x = tcol_array.Value(i + 1, j + 1)
            array[i, j] = x
    return array


def to_toptools_listofshape(shapes):
    """
    Create TopTools_ListOfShape from shapes.

    :param list[OCCT.TopoDS.TopoDS_Shape] shapes: List of shapes.

    :return: TopTools_ListOfShape
    :rtype: OCCT.TopTools.TopTools_ListOfShape
    """
    lst = TopTools_ListOfShape()
    for s in shapes:
        lst.Append(s)
    return lst


def to_lst_from_toptools_listofshape(toptools_list):
    """
    Create a list from a TopTools_ListOfShape.

    :param OCCT.TopTools.TopTools_ListOfShape toptools_list: The
        TopTools_ListOfShape.

    :return: A list of shapes.
    :rtype: list[OCCT.TopoDS.TopoDS_Shape]
    """
    if toptools_list.IsEmpty():
        return []
    lst = []
    while not toptools_list.IsEmpty():
        shp = toptools_list.First()
        lst.append(shp)
        toptools_list.RemoveFirst()
    return lst


def downcast_shape(shape):
    """
    Downcast a shape to its specific type.

    :param OCCT.TopoDS.TopoDS_Shape shape: The shape.

    :return: The downcasted shape.
    :rtype: OCCT.TopoDS.TopoDS_Vertex or OCCT.TopoDS.TopoDS_Edge or
        OCCT.TopoDS.TopoDS_Wire or OCCT.TopoDS.TopoDS_Face or
        OCCT.TopoDS.TopoDS_Shell or OCCT.TopoDS.TopoDS_CompSolid or
        OCCT.TopoDS.TopoDS_Compound
    """
    if shape.ShapeType() == TopAbs_VERTEX:
        return TopoDS.Vertex_(shape)
    if shape.ShapeType() == TopAbs_EDGE:
        return TopoDS.Edge_(shape)
    if shape.ShapeType() == TopAbs_WIRE:
        return TopoDS.Wire_(shape)
    if shape.ShapeType() == TopAbs_FACE:
        return TopoDS.Face_(shape)
    if shape.ShapeType() == TopAbs_SHELL:
        return TopoDS.Shell_(shape)
    if shape.ShapeType() == TopAbs_SOLID:
        return TopoDS.Solid_(shape)
    if shape.ShapeType() == TopAbs_COMPSOLID:
        return TopoDS.CompSolid_(shape)
    if shape.ShapeType() == TopAbs_COMPOUND:
        return TopoDS.Compound_(shape)
    raise ValueError('Could not downcast shape.')
