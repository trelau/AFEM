from OCC.BRepAlgoAPI import BRepAlgoAPI_Fuse
from OCC.ShapeBuild import ShapeBuild_ReShape

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
    def _replace_faces(_part):
        reshape = ShapeBuild_ReShape()
        perform = False
        for f in _part.faces:
            mod = bop.Modified(f)
            if mod.IsEmpty():
                continue
            # Create a reshape operation to replace the face. Put all modified
            # shapes into a compound and replace.
            faces = []
            while not mod.IsEmpty():
                shape = mod.First()
                if not str(shape) in used_faces:
                    faces.append(shape)
                    used_faces.append(str(shape))
                mod.RemoveFirst()
            compound = ShapeTools.make_compound(faces)
            reshape.Replace(f, compound)
            perform = True
        if perform:
            new_shape = reshape.Apply(_part)
            _part.set_shape(new_shape)
        return perform

    # Only replace unique faces so they aren't overlapping between parts.
    # This should avoid duplicate/overlapping faces in splices.
    used_faces = []
    was_performed = []
    for part in [main_part] + other_parts:
        status = _replace_faces(part)
        was_performed.append(status)

    return True in was_performed


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
