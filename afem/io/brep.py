from OCC.BRep import BRep_Builder
from OCC.BRepTools import breptools_Read, breptools_Write
from OCC.TopoDS import TopoDS_Shape


def write_brep(shape, fn):
    """
    Write a BREP file using the shape.

    :param OCC.TopoDS.TopoDS_Shape shape: The shape.
    :param str fn: The filename.

    :return: None.
    """
    breptools_Write(shape, fn)


def read_brep(fn):
    """
    Read a BREP file and return the shape.

    :param str fn: The filename.

    :return: The shape.
    :rtype: OCC.TopoDS.TopoDS_Shape
    """
    shape = TopoDS_Shape()
    builder = BRep_Builder()
    breptools_Read(shape, fn, builder)

    return shape
