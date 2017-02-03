from .checker import CheckGeom
from .creator import CreateGeom
from .intersector import IntersectGeom
from .projector import ProjectGeom


class GeomTools(object):
    """
    Geometry tools.
    """
    check = CheckGeom()
    create = CreateGeom()
    project = ProjectGeom()
    intersect = IntersectGeom()
