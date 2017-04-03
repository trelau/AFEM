from .methods.create_parts import create_bulkhead_by_sref, create_curve_part, \
    create_curve_part_by_section, create_floor_by_sref, create_frame_by_sref, \
    create_frames_at_shapes, create_frames_between_planes, \
    create_skin_from_body, create_skin_from_solid, create_surface_part, \
    create_wing_part_between_geom, create_wing_part_by_params, \
    create_wing_part_by_points, create_wing_part_by_sref, \
    create_wing_parts_along_curve, create_wing_parts_between_planes


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

    @staticmethod
    def along_curve(name, wing, curve, geom1, geom2, maxd=None, npts=None,
                    ref_pln=None, u1=None, u2=None, s1=None, s2=None, indx=1):
        """
        Create spars along a curve.

        :param name:
        :param wing:
        :param curve:
        :param geom1:
        :param geom2:
        :param maxd:
        :param npts:
        :param ref_pln:
        :param u1:
        :param u2:
        :param s1:
        :param s2:
        :param indx:

        :return:
        """
        return create_wing_parts_along_curve('spar', name, wing, curve, geom1,
                                             geom2, maxd, npts, ref_pln, u1,
                                             u2, s1, s2, indx)


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

    @staticmethod
    def along_curve(name, wing, curve, geom1, geom2, maxd=None, npts=None,
                    ref_pln=None, u1=None, u2=None, s1=None, s2=None, indx=1):
        """
        Create ribs along a curve.

        :param name:
        :param wing:
        :param curve:
        :param geom1:
        :param geom2:
        :param maxd:
        :param npts:
        :param ref_pln:
        :param u1:
        :param u2:
        :param s1:
        :param s2:
        :param indx:

        :return:
        """
        return create_wing_parts_along_curve('rib', name, wing, curve, geom1,
                                             geom2, maxd, npts, ref_pln, u1,
                                             u2, s1, s2, indx)


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
    def from_solid(name, solid):
        """
        Create skin from the outer shell of a solid.
        
        :param name: 
        :param solid:
         
        :return: 
        """
        return create_skin_from_solid(name, solid)

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
    def curve_part(name, curve_shape=None, shape1=None, shape2=None):
        """
        Create a curve part.
        
        :param str name: Part name. 
        :param curve_shape: Curve or shape for part.
        :param shape1:
        :param shape2:
        
        :return: Curve part of *None* if method fails.
        :rtype: :class:`.CurvePart`
        """
        if curve_shape:
            return create_curve_part(name, curve_shape)
        if shape1 and shape2:
            return create_curve_part_by_section(name, shape1, shape2)
        return None

    @staticmethod
    def surface_part(name, rshape, *bodies):
        """
        Create a surface part.

        :param str name: Part name.
        :param surface_like rshape: Part reference shape.
        :param bodies:

        :return: Surface part or *None* if method fails.
        :return: :class:`.SurfacePart`
        """
        return create_surface_part(name, rshape, *bodies)
