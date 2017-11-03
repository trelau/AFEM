from OCC.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCC.gce import gce_MakeMirror


def mirror_shape(shape, pln):
    """
    Mirror a shape about a plane.

    :param OCC.TopoDS.TopoDS_Shape shape: The shape.
    :param afem.geometry.entities.Plane pln: The plane.

    :return: The mirrored shape.
    :rtype: OCC.TopoDS.TopoDS_Shape

    :raise RuntimeError: If the transformation fails or is not done.
    """
    trsf = gce_MakeMirror(pln.object.Pln()).Value()
    builder = BRepBuilderAPI_Transform(shape, trsf, True)
    if not builder.IsDone():
        raise RuntimeError('Failed to mirror the shape. The builder is not done.')
    return builder.Shape()
