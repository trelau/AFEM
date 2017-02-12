from .methods.create_parts import create_bulkhead_by_sref, \
    create_floor_by_sref, create_frame_by_sref, create_frames_at_shapes, \
    create_frames_between_planes, create_skin_from_body, create_surface_part, \
    create_wing_part_between_geom, create_wing_part_by_params, \
    create_wing_part_by_points, create_wing_part_by_sref, \
    create_wing_parts_between_planes


class CreateSpar(object):
    """
    Spar creator.
    """

    @staticmethod
    def by_parameters(name, wing, u1, v1, u2, v2, rshape=None, build=True):
        """
        Create a spar by wing parameters.

        :param name:
        :param wing:
        :param u1:
        :param v1:
        :param u2:
        :param v2:
        :param rshape:
        :param build:

        :return:
        """
        return create_wing_part_by_params('spar', name, wing, u1, v1, u2, v2,
                                          rshape, build)

    @staticmethod
    def by_points(name, wing, p1, p2, rshape=None, build=True):
        """
        Create a spar by points.

        :param name:
        :param wing:
        :param p1:
        :param p2:
        :param rshape:
        :param build:

        :return:
        """
        return create_wing_part_by_points('spar', name, wing, p1, p2,
                                          rshape, build)

    @staticmethod
    def by_sref(name, wing, rshape, build=True):
        """
        Create a spar using a reference shape.

        :param name:
        :param wing:
        :param rshape:
        :param build:

        :return:
        """
        return create_wing_part_by_sref('spar', name, wing, rshape, build)

    @staticmethod
    def between_geom(name, wing, geom1, geom2, rshape, build=True):
        """
        Create a spar between geometry.

        :param name:
        :param wing:
        :param geom1:
        :param geom2:
        :param rshape:
        :param build:

        :return:
        """
        return create_wing_part_between_geom('spar', name, wing, geom1, geom2,
                                             rshape, build)

    @staticmethod
    def between_planes(name, wing, planes, geom1, geom2, maxd=None,
                       nplns=None, indx=1):
        """
        Create spars evenly spaced between planes

        :param name:
        :param wing:
        :param planes:
        :param geom1:
        :param geom2:
        :param maxd:
        :param nplns:
        :param indx:

        :return:
        """
        return create_wing_parts_between_planes('spar', name, wing, planes,
                                                geom1, geom2, maxd, nplns,
                                                indx)


class CreateRib(object):
    """
    Rib creator.
    """

    @staticmethod
    def by_parameters(name, wing, u1, v1, u2, v2, rshape=None, build=True):
        """
        Create a rib by wing parameters.

        :param name:
        :param wing:
        :param u1:
        :param v1:
        :param u2:
        :param v2:
        :param rshape:
        :param build:

        :return:
        """
        return create_wing_part_by_params('rib', name, wing, u1, v1, u2, v2,
                                          rshape, build)

    @staticmethod
    def by_points(name, wing, p1, p2, rshape=None, build=True):
        """
        Create a rib by points.

        :param name:
        :param wing:
        :param p1:
        :param p2:
        :param rshape:
        :param build:

        :return:
        """
        return create_wing_part_by_points('rib', name, wing, p1, p2, rshape,
                                          build)

    @staticmethod
    def by_sref(name, wing, rshape, build=True):
        """
        Create a rib using a reference shape.

        :param name:
        :param wing:
        :param rshape:
        :param build:

        :return:
        """
        return create_wing_part_by_sref('rib', name, wing, rshape, build)

    @staticmethod
    def between_geom(name, wing, geom1, geom2, rshape, build=True):
        """
        Create a rib between geometry.

        :param name:
        :param wing:
        :param geom1:
        :param geom2:
        :param rshape:
        :param build:

        :return:
        """
        return create_wing_part_between_geom('rib', name, wing, geom1, geom2,
                                             rshape, build)

    @staticmethod
    def between_planes(name, wing, planes, geom1, geom2, maxd=None,
                       nplns=None, indx=1):
        """
        Create ribs evenly spaced between planes

        :param name:
        :param wing:
        :param planes:
        :param geom1:
        :param geom2:
        :param maxd:
        :param nplns:
        :param indx:

        :return:
        """
        return create_wing_parts_between_planes('rib', name, wing, planes,
                                                geom1, geom2, maxd, nplns,
                                                indx)


class CreateBulkhead(object):
    """
    Bulkhead creation.
    """

    @staticmethod
    def by_sref(name, fuselage, rshape, build=True):
        """
        Create a bulkhead by reference shape.

        :param name:
        :param fuselage:
        :param rshape:
        :param build:

        :return:
        """
        return create_bulkhead_by_sref(name, fuselage, rshape, build)


class CreateFloor(object):
    """
    Floor creation.
    """

    @staticmethod
    def by_sref(name, fuselage, rshape, build=True):
        """
        Create a floor by reference shape.

        :param name:
        :param fuselage:
        :param rshape:
        :param build:

        :return:
        """
        return create_floor_by_sref(name, fuselage, rshape, build)


class CreateFrame(object):
    """
    Frame creation.
    """

    @staticmethod
    def by_sref(name, fuselage, rshape, height):
        """
        Create a frame by reference shape.

        :param name:
        :param fuselage:
        :param rshape:
        :param height:

        :return:
        """
        return create_frame_by_sref(name, fuselage, rshape, height)

    @staticmethod
    def between_planes(name, fuselage, planes, height, maxd=None, nplns=None,
                       indx=1):
        """
        Create frames between planes.

        :param name:
        :param fuselage:
        :param planes:
        :param height:
        :param maxd:
        :param nplns:
        :param indx:

        :return:
        """
        return create_frames_between_planes(name, fuselage, planes, height,
                                            maxd, nplns, indx)

    @staticmethod
    def at_shapes(name, fuselage, shapes, height, indx=1):
        """
        Create frames at shapes.

        :param name:
        :param fuselage:
        :param shapes:
        :param height:
        :param indx:

        :return:
        """
        return create_frames_at_shapes(name, fuselage, shapes, height, indx)


class CreateSkin(object):
    """
    Skin creation.
    """

    @staticmethod
    def from_body(name, body):
        """
        Create skin from the outer shell of a body.

        :param name:
        :param body:

        :return:
        """
        return create_skin_from_body(name, body)


class CreatePart(object):
    """
    Part creator.
    """
    spar = CreateSpar()
    rib = CreateRib()
    bulkhead = CreateBulkhead()
    floor = CreateFloor()
    frame = CreateFrame()
    skin = CreateSkin()

    @staticmethod
    def surface_part(name, rshape, *bodies):
        """
        Create a surface frame.

        :param str name: Part name.
        :param surface_like rshape: Part reference shape.
        :param bodies:

        :return: Surface frame or *None* if method fails.
        :return: :class:`.SurfacePart`
        """
        return create_surface_part(name, rshape, *bodies)
