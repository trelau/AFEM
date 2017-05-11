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
    def by_parameters(label, wing, u1, v1, u2, v2, surface_shape=None,
                      bodies=()):
        """
        Create a spar by wing parameters.

        :param label:
        :param wing:
        :param u1:
        :param v1:
        :param u2:
        :param v2:
        :param surface_shape:
        :param bodies:

        :return:
        """
        return create_wing_part_by_params('spar', label, wing, u1, v1, u2, v2,
                                          surface_shape, bodies)

    @staticmethod
    def by_points(label, wing, p1, p2, surface_shape=None, bodies=()):
        """
        Create a spar by points.

        :param label:
        :param wing:
        :param p1:
        :param p2:
        :param surface_shape:
        :param bodies:

        :return:
        """
        return create_wing_part_by_points('spar', label, wing, p1, p2,
                                          surface_shape, bodies)

    @staticmethod
    def by_sref(label, wing, surface_shape, bodies=()):
        """
        Create a spar using a reference shape.

        :param label:
        :param wing:
        :param surface_shape:
        :param bodies:

        :return:
        """
        return create_wing_part_by_sref('spar', label, wing, surface_shape,
                                        bodies)

    @staticmethod
    def between_geom(label, wing, geom1, geom2, surface_shape, bodies=()):
        """
        Create a spar between geometry.

        :param label:
        :param wing:
        :param geom1:
        :param geom2:
        :param surface_shape:
        :param bodies:

        :return:
        """
        return create_wing_part_between_geom('spar', label, wing, geom1, geom2,
                                             surface_shape, bodies)

    @staticmethod
    def between_planes(label, wing, planes, geom1, geom2, maxd=None,
                       nplns=None, indx=1):
        """
        Create spars evenly spaced between planes

        :param label:
        :param wing:
        :param planes:
        :param geom1:
        :param geom2:
        :param maxd:
        :param nplns:
        :param indx:

        :return:
        """
        return create_wing_parts_between_planes('spar', label, wing, planes,
                                                geom1, geom2, maxd, nplns,
                                                indx)

    @staticmethod
    def along_curve(label, wing, curve, geom1, geom2, maxd=None, npts=None,
                    ref_pln=None, u1=None, u2=None, s1=None, s2=None, indx=1):
        """
        Create spars along a curve.

        :param label:
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
        return create_wing_parts_along_curve('spar', label, wing, curve, geom1,
                                             geom2, maxd, npts, ref_pln, u1,
                                             u2, s1, s2, indx)


class CreateRib(object):
    """
    Rib creator.
    """

    @staticmethod
    def by_parameters(label, wing, u1, v1, u2, v2, surface_shape=None,
                      bodies=()):
        """
        Create a rib by wing parameters.

        :param label:
        :param wing:
        :param u1:
        :param v1:
        :param u2:
        :param v2:
        :param surface_shape:
        :param bodies:

        :return:
        """
        return create_wing_part_by_params('rib', label, wing, u1, v1, u2, v2,
                                          surface_shape, bodies)

    @staticmethod
    def by_points(label, wing, p1, p2, surface_shape=None, bodies=()):
        """
        Create a rib by points.

        :param label:
        :param wing:
        :param p1:
        :param p2:
        :param surface_shape:
        :param bodies:

        :return:
        """
        return create_wing_part_by_points('rib', label, wing, p1, p2,
                                          surface_shape,
                                          bodies)

    @staticmethod
    def by_sref(label, wing, surface_shape, bodies=()):
        """
        Create a rib using a reference shape.

        :param label:
        :param wing:
        :param surface_shape:
        :param bodies:

        :return:
        """
        return create_wing_part_by_sref('rib', label, wing, surface_shape,
                                        bodies)

    @staticmethod
    def between_geom(label, wing, geom1, geom2, surface_shape, bodies=()):
        """
        Create a rib between geometry.

        :param label:
        :param wing:
        :param geom1:
        :param geom2:
        :param surface_shape:
        :param bodies:

        :return:
        """
        return create_wing_part_between_geom('rib', label, wing, geom1, geom2,
                                             surface_shape, bodies)

    @staticmethod
    def between_planes(label, wing, planes, geom1, geom2, maxd=None,
                       nplns=None, indx=1):
        """
        Create ribs evenly spaced between planes

        :param label:
        :param wing:
        :param planes:
        :param geom1:
        :param geom2:
        :param maxd:
        :param nplns:
        :param indx:

        :return:
        """
        return create_wing_parts_between_planes('rib', label, wing, planes,
                                                geom1, geom2, maxd, nplns,
                                                indx)

    @staticmethod
    def along_curve(label, wing, curve, geom1, geom2, maxd=None, npts=None,
                    ref_pln=None, u1=None, u2=None, s1=None, s2=None, indx=1):
        """
        Create ribs along a curve.

        :param label:
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
        return create_wing_parts_along_curve('rib', label, wing, curve, geom1,
                                             geom2, maxd, npts, ref_pln, u1,
                                             u2, s1, s2, indx)


class CreateBulkhead(object):
    """
    Bulkhead creation.
    """

    @staticmethod
    def by_sref(label, fuselage, surface_shape, bodies=()):
        """
        Create a bulkhead by reference shape.

        :param label:
        :param fuselage:
        :param surface_shape:
        :param bodies:

        :return:
        """
        return create_bulkhead_by_sref(label, fuselage, surface_shape, bodies)


class CreateFloor(object):
    """
    Floor creation.
    """

    @staticmethod
    def by_sref(label, fuselage, surface_shape, bodies=()):
        """
        Create a floor by reference shape.

        :param label:
        :param fuselage:
        :param surface_shape:
        :param bodies:

        :return:
        """
        return create_floor_by_sref(label, fuselage, surface_shape, bodies)


class CreateFrame(object):
    """
    Frame creation.
    """

    @staticmethod
    def by_sref(label, fuselage, surface_shape, height):
        """
        Create a frame by reference shape.

        :param label:
        :param fuselage:
        :param surface_shape:
        :param height:

        :return:
        """
        return create_frame_by_sref(label, fuselage, surface_shape, height)

    @staticmethod
    def between_planes(label, fuselage, planes, height, maxd=None, nplns=None,
                       indx=1):
        """
        Create frames between planes.

        :param label:
        :param fuselage:
        :param planes:
        :param height:
        :param maxd:
        :param nplns:
        :param indx:

        :return:
        """
        return create_frames_between_planes(label, fuselage, planes, height,
                                            maxd, nplns, indx)

    @staticmethod
    def at_shapes(label, fuselage, shapes, height, indx=1):
        """
        Create frames at shapes.

        :param label:
        :param fuselage:
        :param shapes:
        :param height:
        :param indx:

        :return:
        """
        return create_frames_at_shapes(label, fuselage, shapes, height, indx)


class CreateSkin(object):
    """
    Skin creation.
    """

    @staticmethod
    def from_solid(label, solid, copy=False):
        """
        Create skin from the outer shell of a solid.
        
        :param label: 
        :param solid:
        :param copy:
         
        :return: 
        """
        return create_skin_from_solid(label, solid, copy)

    @staticmethod
    def from_body(label, body, copy=False):
        """
        Create skin from the outer shell of a body.

        :param label:
        :param body:
        :param copy:

        :return:
        """
        return create_skin_from_body(label, body, copy)


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
    def curve_part(label, curve_shape=None, shape1=None, shape2=None):
        """
        Create a curve part.
        
        :param str label: Part label. 
        :param curve_shape: Curve or shape for part.
        :param shape1:
        :param shape2:
        
        :return: Curve part of *None* if method fails.
        :rtype: :class:`.CurvePart`
        """
        if curve_shape:
            return create_curve_part(label, curve_shape)
        if shape1 and shape2:
            return create_curve_part_by_section(label, shape1, shape2)
        return None

    @staticmethod
    def surface_part(label, surface_shape, bodies=()):
        """
        Create a surface part.

        :param str label: Part label.
        :param surface_like surface_shape: Part reference shape.
        :param bodies:

        :return: Surface part or *None* if method fails.
        :return: :class:`.SurfacePart`
        """
        return create_surface_part(label, surface_shape, bodies)
