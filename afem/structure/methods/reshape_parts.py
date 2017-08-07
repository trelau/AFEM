from OCC.BRep import BRep_Builder
from OCC.ShapeBuild import ShapeBuild_ReShape
from OCC.TopTools import TopTools_IndexedMapOfShape
from OCC.TopoDS import TopoDS_Compound

__all__ = []


def reshape_parts(tool, parts):
    """
    Update the part shape if modified by a tool.
    """
    status = []
    index_map = TopTools_IndexedMapOfShape()
    for part in parts:
        reshape = ShapeBuild_ReShape()
        performed = False
        for old_shape in part.reshapes:
            # Check deleted.
            if tool.IsDeleted(old_shape):
                # Remove the shape.
                reshape.Remove(old_shape)
                performed = True
                continue

            # TODO Handle generated?

            # Check modified.
            mod = tool.Modified(old_shape)
            if mod.IsEmpty():
                continue

            # Put modified shapes into a compound.
            new_shape = TopoDS_Compound()
            builder = BRep_Builder()
            builder.MakeCompound(new_shape)
            modified = False
            while not mod.IsEmpty():
                shape = mod.First()
                if not index_map.Contains(shape) and \
                        not old_shape.IsSame(shape):
                    builder.Add(new_shape, shape)
                    index_map.Add(shape)
                    modified = True
                mod.RemoveFirst()

            # Replace the old shape with modified shape(s).
            if modified:
                reshape.Replace(old_shape, new_shape)
                performed = True

        # Set the new shape.
        if performed:
            new_shape = reshape.Apply(part)
            part.set_shape(new_shape)
        status.append(performed)

        # Perform for sub-parts.
        for subpart in part.subparts:
            reshape_parts(tool, [subpart])

    return True in status
