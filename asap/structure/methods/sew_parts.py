from OCC.BRepBuilderAPI import BRepBuilderAPI_Sewing
from numpy import mean

from ...topology import ShapeTools


def sew_surface_parts(parts):
    """
    Sew surface parts.
    """
    if len(parts) == 0:
        return False

    # Put other parts into a list.
    parts = list(parts)

    # Calculate tolerances for sewing using the tolerances of the parts.
    tol = mean([ShapeTools.get_tolerance(part, 0) for part in parts])
    max_tol = max([ShapeTools.get_tolerance(part, 1) for part in parts])

    # Perform sewing
    sew = BRepBuilderAPI_Sewing(tol, True, True, True, True)
    sew.SetMaxTolerance(max_tol)
    for part in parts:
        sew.Add(part)
    sew.Perform()

    modified = False
    for part in parts:
        if not sew.IsModified(part):
            continue
        mod = sew.Modified(part)
        shape = ShapeTools.to_shape(mod)
        if shape:
            part.set_shape(shape)
            modified = True

    return modified
