from OCC.BRepAlgoAPI import BRepAlgoAPI_Fuse


def merge_surface_part(part, other):
    """
    Merge surface part with other part or shape.
    """

    if part.IsNull() or other.IsNull():
        return False

    # Fuse the part and the other part or shape.
    bop = BRepAlgoAPI_Fuse(part, other)
    if bop.ErrorStatus() != 0:
        return False

    # The resulting shape will be the original part now.
    part.set_shape(bop.Shape())
    return True
