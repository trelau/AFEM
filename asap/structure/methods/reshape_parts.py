from OCC.ShapeBuild import ShapeBuild_ReShape

from ...topology import ShapeTools


def reshape_surface_parts(tool, parts):
    """
    Update the part shape.
    """
    used_faces = []
    status = []
    for part in parts:
        reshape = ShapeBuild_ReShape()
        performed = False
        for f in part.faces:
            # Check deleted.
            if tool.IsDeleted(f):
                # Remove the face.
                reshape.Remove(f)
                performed = True
                continue

            # Check modified.
            mod = tool.Modified(f)
            if mod.IsEmpty():
                continue
            # Replace the face with modified face(s).
            faces = []
            while not mod.IsEmpty():
                shape = mod.First()
                if not str(shape) in used_faces:
                    faces.append(shape)
                    used_faces.append(str(shape))
                mod.RemoveFirst()
            compound = ShapeTools.make_compound(faces)
            reshape.Replace(f, compound)
            performed = True
        if performed:
            new_shape = reshape.Apply(part)
            part.set_shape(new_shape)
        status.append(performed)

    return True in status
