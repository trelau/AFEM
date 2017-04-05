from OCC.ShapeBuild import ShapeBuild_ReShape

from ...topology import ShapeTools


def reshape_parts(tool, parts):
    """
    Update the part shape if modified by a tool.
    """
    used_shapes = []
    status = []
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

            # Check modified.
            mod = tool.Modified(old_shape)
            if mod.IsEmpty():
                continue
            # Replace the old shape with modified shape(s).
            new_shapes = []
            while not mod.IsEmpty():
                shape = mod.First()
                if not str(shape) in used_shapes:
                    new_shapes.append(shape)
                    used_shapes.append(str(shape))
                mod.RemoveFirst()
            if len(new_shapes) > 1:
                compound = ShapeTools.make_compound(new_shapes)
            else:
                compound = new_shapes[0]
            reshape.Replace(old_shape, compound)
            performed = True
        if performed:
            new_shape = reshape.Apply(part)
            part.set_shape(new_shape)
        status.append(performed)

    return True in status
