from OCC.BOPAlgo import BOPAlgo_BOP, BOPAlgo_FUSE

from ...topology import ShapeTools


def build_surface_part(part, unify=False):
    """
    Build surface part shape.
    """
    # Get the formed shapes of the part. This should usually be a set of
    # compounds containing faces.
    fshapes = part.fshapes
    if not fshapes:
        return False

    if len(fshapes) > 1:
        # Use BOP Fuse algorithm to fuse together all formed shapes of the
        # frame.
        bop = BOPAlgo_BOP()
        bop.SetOperation(BOPAlgo_FUSE)
        bop.AddArgument(fshapes[0])
        for shape in fshapes[1:]:
            bop.AddTool(shape)
        bop.Perform()
        if bop.ErrorStatus() != 0:
            return False
        shape = bop.Shape()
    else:
        shape = fshapes[0]

    # Unify the frame if desired.
    if unify:
        shape = ShapeTools.unify_shape(shape, False, True, False)

    # Set the frame shape to reference the result.
    part.set_shape(shape)
    return True
