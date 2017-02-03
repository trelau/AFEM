from ...topology import ShapeTools


def form_with_solid(shape, solid):
    """
    Form the frame with a solid using a Common operation.
    """
    # Perform BOP Common and get the resulting faces.
    faces = ShapeTools.bcommon(shape, solid, 'face')
    if not faces:
        return None

    # Put all faces in a compound
    compound = ShapeTools.make_compound(faces)
    return compound
