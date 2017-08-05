from OCC.TColStd import (TColStd_Array1OfInteger, TColStd_Array1OfReal,
                         TColStd_Array2OfReal)
from OCC.TColgp import (TColgp_Array1OfPnt, TColgp_Array1OfPnt2d,
                        TColgp_Array2OfPnt, TColgp_HArray1OfPnt)
from OCC.TopTools import TopTools_ListOfShape
from OCC.gp import gp_Pnt, gp_Pnt2d
from numpy import array as np_array, zeros

from afem.utils.check import is_array_like

__all__ = ["to_gp_pnt", "to_gp_pnt2d", "to_np_from_tcolgp_array1_pnt",
           "to_np_from_tcolgp_array2_pnt",
           "to_np_from_tcolstd_array1_integer",
           "to_np_from_tcolstd_array1_real", "to_np_from_tcolstd_array2_real",
           "to_tcolgp_array1_pnt", "to_tcolgp_array1_pnt2d",
           "to_tcolgp_array2_pnt", "to_tcolgp_harray1_pnt",
           "to_tcolstd_array1_integer", "to_tcolstd_array1_real",
           "to_tcolstd_array2_real", "to_toptools_listofshape",
           "to_lst_from_toptools_listofshape"]


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

    :param list[OCC.TopoDS.TopoDS_Shape] shapes: List of shapes.

    :return: TopTools_ListOfShape
    :rtype: OCC.TopTools.TopTools_ListOfShape
    """
    lst = TopTools_ListOfShape()
    for s in shapes:
        lst.Append(s)
    return lst


def to_lst_from_toptools_listofshape(toptools_list):
    """
    Create a list from a TopTools_ListOfShape.

    :param OCC.TopTools.TopTools_ListOfShape toptools_list: The
        TopTools_ListOfShape.

    :return: A list of shapes.
    :rtype: list[OCC.TopoDS.TopoDS_Shape]
    """
    if toptools_list.IsEmpty():
        return []
    lst = []
    while not toptools_list.IsEmpty():
        shp = toptools_list.First()
        lst.append(shp)
        toptools_list.RemoveFirst()
    return lst
