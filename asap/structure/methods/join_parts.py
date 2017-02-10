from OCC.BRepAlgoAPI import BRepAlgoAPI_Fuse

from .reshape_parts import reshape_surface_parts
from ...geometry import IntersectGeom
from ...topology import ShapeTools


def join_surface_parts(main_part, *other_parts):
    """
    Join parts using BOP Fuse.
    """

    if main_part.IsNull() or len(other_parts) == 0:
        return False

    # Put other parts into a compound.
    other_parts = list(other_parts)
    other_compound = ShapeTools.make_compound(other_parts)
    # Fuse the main part and the compound. Putting the other parts in a
    # compound avoids fusing them to each other.
    bop = BRepAlgoAPI_Fuse(main_part, other_compound)
    if bop.ErrorStatus() != 0:
        return False

    # Replace modified face(s) of result into original shapes.
    return reshape_surface_parts(bop, [main_part] + other_parts)


def join_wing_parts(parts, tol=None):
    """
    Attempt to automatically join wing parts.
    """
    # Test all combinations of the parts for potential intersection using
    # the part reference curve.
    join_parts = []
    main_parts = []
    nparts = len(parts)
    for i in range(0, nparts - 1):
        main = parts[i]
        other_parts = []
        for j in range(i + 1, nparts):
            other = parts[j]
            if None in [main.cref, other.cref]:
                continue
            if tol is None:
                tol1 = ShapeTools.get_tolerance(main, 1)
                tol2 = ShapeTools.get_tolerance(other, 1)
                tol = max(tol1, tol2)
            cci = IntersectGeom.perform(main.cref, other.cref, tol)
            if not cci.success:
                continue
            # Store potential join.
            other_parts.append(other)
        if other_parts:
            main_parts.append(main)
            join_parts.append(other_parts)

    # Join the parts using the key as the main part and the list as other
    # parts.
    status_out = False
    for main, other_parts in zip(main_parts, join_parts):
        if join_surface_parts(main, *other_parts):
            status_out = True

    return status_out
