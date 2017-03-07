from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeFace

from ...topology import ShapeTools


def form_with_solid(shape, solid):
    """
    Form the frame with a solid using a Common operation.
    """
    # Perform BOP Common and get the resulting faces.
    faces = ShapeTools.bcommon(shape, solid, 'face')
    if not faces:
        faces = _hack_form_planar(shape, solid)
    if not faces:
        return None

    # Put all faces in a compound
    compound = ShapeTools.make_compound(faces)
    return compound


def _hack_form_planar(shape, solid):
    """
    Hopefully a temporary hack for a failing BOP with a planar shape.
    """
    edges = ShapeTools.bsection(shape, solid, 'edge')
    if not edges:
        return []

    wires = ShapeTools.connect_edges(edges)
    if not wires:
        return []

    faces = []
    for w in wires:
        if not w.Closed():
            continue
        builder = BRepBuilderAPI_MakeFace(w, True)
        if builder.IsDone():
            faces.append(builder.Face())
    return faces
