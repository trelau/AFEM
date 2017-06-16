from .methods.create_parts import create_bulkhead_by_sref, create_curve_part, \
    create_curve_part_by_section, create_floor_by_sref, create_frame_by_sref, \
    create_frames_at_shapes, create_frames_between_planes, \
    create_skin_from_body, create_skin_from_solid, create_stiffener1d, \
    create_stiffener2d_by_section, create_stiffener2d_by_sections, \
    create_stiffener2d_by_wire, create_surface_part, \
    create_wing_part_between_geom, create_wing_part_by_params, \
    create_wing_part_by_points, create_wing_part_by_sref, \
    create_wing_parts_along_curve, create_wing_parts_between_planes


# TODO Check to see if parts are created at planes for between planes method.


class CreateSpar(object):
    """
    Spar creator.
    """

    @staticmethod
    def by_parameters(label, wing, u1, v1, u2, v2, surface_shape=None,
                      bodies=(), assy=None):
        """
        Create a spar by parameters.

        :param str label: Part label.
        :param wing: Wing to build part in.
        :type wing: :class:`.Wing`
        :param float u1: Parameter in u-direction at starting point
            (sref.u1 <= u1 <= sref.u2).
        :param float v1: Parameter in v-direction at starting point
            (sref.v1 <= v1 <= sref.v2).
        :param float u2: Parameter in u-direction at ending point
            (sref.u1 <= u2 <= sref.u2).
        :param float v2: Parameter in v-direction at ending point
            (sref.v1 <= v2 <= sref.v2).
        :param surface_like surface_shape: The basis shape to define the
            shape of the spar. If none is provided, then a plane will be
            defined between (u1, v1), (u2, v2), and a point translated from
            the reference surface normal at (u1, v1).
        :param iterable bodies: Extra bodies to be used to define the initial
            part shape.
        :param assy: The assembly to store the part in. If none is provided
            then the active assembly will be used by default. If *False* is
            provided then the part will not be added to any assembly.
        :type assy: str, :class:`.Assembly`, or bool

        :return: Spar part.
        :rtype: :class:`.Spar`

        **Notes**:

        - The wing component must have a reference surface (sref) available
          to evaluate the parameters. The u- and v-direction parameters are
          dependent on the parametrization of the reference surface.
          Typically, the u-direction is in the chord direction and the
          v-direction is along the span.

        - If any segment of the part is outside the wing reference surface,
          only the segment nearest the point at (u1, v1) is used.

        """
        return create_wing_part_by_params('spar', label, wing, u1, v1, u2, v2,
                                          surface_shape, bodies, assy)

    @staticmethod
    def by_points(label, wing, p1, p2, surface_shape=None, bodies=(),
                  assy=None):
        """
        Create a spar by points.

        :param str label: Part label.
        :param wing: Wing to build part in.
        :type wing: :class:`.Wing`
        :param p1: Starting point.
        :param p2: Ending point.
        :param surface_like surface_shape: The basis shape to define the
            shape of the spar. If none is provided, then a plane will be
            defined between p1, p2, and a point translated from
            the reference surface normal at p1.
        :param iterable bodies: Extra bodies to be used to define the initial
            part shape.
        :param assy: The assembly to store the part in. If none is provided
            then the active assembly will be used by default. If *False* is
            provided then the part will not be added to any assembly.
        :type assy: str, :class:`.Assembly`, or bool

        :return: Spar part.
        :rtype: :class:`.Spar`

        **Notes**:

        - The wing component must have a reference surface (sref) available
          to project the points to. This method first inverts the points and
          then uses the *by_parameters()* method.

        - If any segment of the part is outside the wing reference surface,
          only the segment nearest the point at p1 is used.

        """
        return create_wing_part_by_points('spar', label, wing, p1, p2,
                                          surface_shape, bodies, assy)

    @staticmethod
    def by_sref(label, wing, surface_shape, bodies=(), assy=None):
        """
        Create a spar by a basis shape.

        :param str label: Part label.
        :param wing: Wing to build part in.
        :type wing: :class:`.Wing`
        :param surface_like surface_shape: The basis shape to define the
            shape of the spar.
        :param iterable bodies: Extra bodies to be used to define the initial
            part shape.
        :param assy: The assembly to store the part in. If none is provided
            then the active assembly will be used by default. If *False* is
            provided then the part will not be added to any assembly.
        :type assy: str, :class:`.Assembly`, or bool

        :return: Spar part.
        :rtype: :class:`.Spar`

        **Notes**:

        - If any segment of the part is outside the wing reference surface
          only the longest segment will be used.

        - If available, the reference curve will be generated by the
          intersection of the basis shape and the wing reference surface.
          The orientation of this curve will be such that the starting point
          is nearest the (u1, v1) corner of the wing reference surface. For
          example, if the (u1, v1) corner is the root leading edge, then the
          part should be oriented such that its positive direction is
          outboard.

        """
        return create_wing_part_by_sref('spar', label, wing, surface_shape,
                                        bodies, assy)

    @staticmethod
    def between_geom(label, wing, geom1, geom2, surface_shape, bodies=(),
                     assy=None):
        """
        Create a spar between geometry.

        :param str label: Part label.
        :param wing: Wing to build part in.
        :type wing: :class:`.Wing`
        :param geom1: Starting geometry.
        :param geom2: Ending geometry.
        :param surface_like surface_shape: The basis shape to define the
            shape of the spar.
        :param iterable bodies: Extra bodies to be used to define the initial
            part shape.
        :param assy: The assembly to store the part in. If none is provided
            then the active assembly will be used by default. If *False* is
            provided then the part will not be added to any assembly.
        :type assy: str, :class:`.Assembly`, or bool

        :return: Spar part.
        :rtype: :class:`.Spar`

        **Notes**:

        - If any segment of the part is outside the wing reference surface
          only the segment nearest *geom1* will be used.

        - This method determine the part reference curve by intersecting the
          basis shape with the wing reference surface. It then intersects
          this curve with *geom1* and *geom2*, resulting in two points. These
          two points are used in the *by_points()* method.

        """
        return create_wing_part_between_geom('spar', label, wing, geom1, geom2,
                                             surface_shape, bodies, assy)

    @staticmethod
    def between_planes(label, wing, planes, geom1, geom2, maxd=None,
                       nplns=None, indx=1, assy=None):
        """
        Create planar spars between planes.

        :param str label: Part label.
        :param wing: Wing to build part in.
        :type wing: :class:`.Wing`
        :param iterable planes: Iterable of planes that parts will be
            created between. Parts are not created at the given planes.
        :param geom1: Starting geometry.
        :param geom2: Ending geometry.
        :param maxd: Maximum allowed spacing.
        :param nplns: Number of planes between the given planes.
        :param int indx: Index to append to label (space delimited). This
            creates unique part labels for multiple parts.
        :param assy: The assembly to store the part in. If none is provided
            then the active assembly will be used by default. If *False* is
            provided then the part will not be added to any assembly.
        :type assy: str, :class:`.Assembly`, or bool

        :return: List of :class:`.Spar` instances. Empty list is returned if
            no parts are created.
        :rtype: list

        **Notes**:

        - If both *maxd* and *nplns* are provided then *nplns* is used as a
          required minimum.

        """
        return create_wing_parts_between_planes('spar', label, wing, planes,
                                                geom1, geom2, maxd, nplns,
                                                indx, assy)

    @staticmethod
    def along_curve(label, wing, curve, geom1, geom2, maxd=None, npts=None,
                    ref_pln=None, u1=None, u2=None, s1=None, s2=None, indx=1,
                    assy=None):
        """
        Create planar spars along a curve.

        :param str label: Part label.
        :param wing: Wing to build part in.
        :type wing: :class:`.Wing`
        :param curve_like curve: Curve to generate parts along.
        :param geom1:
        :param geom2:
        :param maxd: Maximum allowed spacing.
        :param npts: Number of parts along curve.
        :param ref_pln: Reference plane for orienting planes along the
            curve. If none is provided then the local curve derivative is
            used to define the plane normal.
        :param u1: Starting parameter along curve.
        :param u2: Ending parameter along curve.
        :param s1: Offset distance for starting point.
        :param s2: Offset distane from end of curve (negative value).
        :param int indx: Index to append to label (space delimited). This
            creates unique part labels for multiple parts.
        :param assy: The assembly to store the part in. If none is provided
            then the active assembly will be used by default. If *False* is
            provided then the part will not be added to any assembly.
        :type assy: str, :class:`.Assembly`, or bool

        :return: List of :class:`.Spar` instances. Empty list is returned if
            no parts are created.
        :rtype: list

        **Notes**:

        - If both *maxd* and *npts* are provided then *npts* is used as a
          required minimum.

        """
        return create_wing_parts_along_curve('spar', label, wing, curve, geom1,
                                             geom2, maxd, npts, ref_pln, u1,
                                             u2, s1, s2, indx, assy)


class CreateRib(object):
    """
    Rib creator.
    """

    @staticmethod
    def by_parameters(label, wing, u1, v1, u2, v2, surface_shape=None,
                      bodies=(), assy=None):
        """
        Create a rib by parameters.

        :param str label: Part label.
        :param wing: Wing to build part in.
        :type wing: :class:`.Wing`
        :param float u1: Parameter in u-direction at starting point
            (sref.u1 <= u1 <= sref.u2).
        :param float v1: Parameter in v-direction at starting point
            (sref.v1 <= v1 <= sref.v2).
        :param float u2: Parameter in u-direction at ending point
            (sref.u1 <= u2 <= sref.u2).
        :param float v2: Parameter in v-direction at ending point
            (sref.v1 <= v2 <= sref.v2).
        :param surface_like surface_shape: The basis shape to define the
            shape of the rib. If none is provided, then a plane will be
            defined between (u1, v1), (u2, v2), and a point translated from
            the reference surface normal at (u1, v1).
        :param iterable bodies: Extra bodies to be used to define the initial
            part shape.
        :param assy: The assembly to store the part in. If none is provided
            then the active assembly will be used by default. If *False* is
            provided then the part will not be added to any assembly.
        :type assy: str, :class:`.Assembly`, or bool

        :return: Rib part.
        :rtype: :class:`.Rib`

        **Notes**:

        - The wing component must have a reference surface (sref) available
          to evaluate the parameters. The u- and v-direction parameters are
          dependent on the parametrization of the reference surface.
          Typically, the u-direction is in the chord direction and the
          v-direction is along the span.

        - If any segment of the part is outside the wing reference surface,
          only the segment nearest the point at (u1, v1) is used.

        """
        return create_wing_part_by_params('rib', label, wing, u1, v1, u2, v2,
                                          surface_shape, bodies, assy)

    @staticmethod
    def by_points(label, wing, p1, p2, surface_shape=None, bodies=(),
                  assy=None):
        """
        Create a rib by points.

        :param str label: Part label.
        :param wing: Wing to build part in.
        :type wing: :class:`.Wing`
        :param p1: Starting point.
        :param p2: Ending point.
        :param surface_like surface_shape: The basis shape to define the
            shape of the rib. If none is provided, then a plane will be
            defined between p1, p2, and a point translated from
            the reference surface normal at p1.
        :param iterable bodies: Extra bodies to be used to define the initial
            part shape.
        :param assy: The assembly to store the part in. If none is provided
            then the active assembly will be used by default. If *False* is
            provided then the part will not be added to any assembly.
        :type assy: str, :class:`.Assembly`, or bool

        :return: Rib part.
        :rtype: :class:`.Rib`

        **Notes**:

        - The wing component must have a reference surface (sref) available
          to project the points to. This method first inverts the points and
          then uses the *by_parameters()* method.

        - If any segment of the part is outside the wing reference surface,
          only the segment nearest the point at p1 is used.

        """
        return create_wing_part_by_points('rib', label, wing, p1, p2,
                                          surface_shape,
                                          bodies, assy)

    @staticmethod
    def by_sref(label, wing, surface_shape, bodies=(), assy=None):
        """
        Create a rib by a basis shape.

        :param str label: Part label.
        :param wing: Wing to build part in.
        :type wing: :class:`.Wing`
        :param surface_like surface_shape: The basis shape to define the
            shape of the rib.
        :param iterable bodies: Extra bodies to be used to define the initial
            part shape.
        :param assy: The assembly to store the part in. If none is provided
            then the active assembly will be used by default. If *False* is
            provided then the part will not be added to any assembly.
        :type assy: str, :class:`.Assembly`, or bool

        :return: Rib part.
        :rtype: :class:`.Rib`

        **Notes**:

        - If any segment of the part is outside the wing reference surface
          only the longest segment will be used.

        - If available, the reference curve will be generated by the
          intersection of the basis shape and the wing reference surface.
          The orientation of this curve will be such that the starting point
          is nearest the (u1, v1) corner of the wing reference surface. For
          example, if the (u1, v1) corner is the root leading edge, then the
          part should be oriented such that its positive direction is
          outboard.

        """
        return create_wing_part_by_sref('rib', label, wing, surface_shape,
                                        bodies, assy)

    @staticmethod
    def between_geom(label, wing, geom1, geom2, surface_shape, bodies=(),
                     assy=None):
        """
        Create a rib between geometry.

        :param str label: Part label.
        :param wing: Wing to build part in.
        :type wing: :class:`.Wing`
        :param geom1: Starting geometry.
        :param geom2: Ending geometry.
        :param surface_like surface_shape: The basis shape to define the
            shape of the rib.
        :param iterable bodies: Extra bodies to be used to define the initial
            part shape.
        :param assy: The assembly to store the part in. If none is provided
            then the active assembly will be used by default. If *False* is
            provided then the part will not be added to any assembly.
        :type assy: str, :class:`.Assembly`, or bool

        :return: Rib part.
        :rtype: :class:`.Rib`

        **Notes**:

        - If any segment of the part is outside the wing reference surface
          only the segment nearest *geom1* will be used.

        - This method determine the part reference curve by intersecting the
          basis shape with the wing reference surface. It then intersects
          this curve with *geom1* and *geom2*, resulting in two points. These
          two points are used in the *by_points()* method.

        """
        return create_wing_part_between_geom('rib', label, wing, geom1, geom2,
                                             surface_shape, bodies, assy)

    @staticmethod
    def between_planes(label, wing, planes, geom1, geom2, maxd=None,
                       nplns=None, indx=1, assy=None):
        """
        Create planar ribs between planes.

        :param str label: Part label.
        :param wing: Wing to build part in.
        :type wing: :class:`.Wing`
        :param iterable planes: Iterable of planes that parts will be
            created between. Parts are not created at the given planes.
        :param geom1: Starting geometry.
        :param geom2: Ending geometry.
        :param maxd: Maximum allowed spacing.
        :param nplns: Number of planes between the given planes.
        :param int indx: Index to append to label (space delimited). This
            creates unique part labels for multiple parts.
        :param assy: The assembly to store the part in. If none is provided
            then the active assembly will be used by default. If *False* is
            provided then the part will not be added to any assembly.
        :type assy: str, :class:`.Assembly`, or bool

        :return: List of :class:`.Rib` instances. Empty list is returned if
            no parts are created.
        :rtype: list

        **Notes**:

        - If both *maxd* and *nplns* are provided then *nplns* is used as a
          required minimum.

        """
        return create_wing_parts_between_planes('rib', label, wing, planes,
                                                geom1, geom2, maxd, nplns,
                                                indx, assy)

    @staticmethod
    def along_curve(label, wing, curve, geom1, geom2, maxd=None, npts=None,
                    ref_pln=None, u1=None, u2=None, s1=None, s2=None, indx=1,
                    assy=None):
        """
        Create planar ribs along a curve.

        :param str label: Part label.
        :param wing: Wing to build part in.
        :type wing: :class:`.Wing`
        :param curve_like curve: Curve to generate parts along.
        :param geom1:
        :param geom2:
        :param maxd: Maximum allowed spacing.
        :param npts: Number of parts along curve.
        :param ref_pln: Reference plane for orienting planes along the
            curve. If none is provided then the local curve derivative is
            used to define the plane normal.
        :param u1: Starting parameter along curve.
        :param u2: Ending parameter along curve.
        :param s1: Offset distance for starting point.
        :param s2: Offset distance from end of curve (negative value).
        :param int indx: Index to append to label (space delimited). This
            creates unique part labels for multiple parts.
        :param assy: The assembly to store the part in. If none is provided
            then the active assembly will be used by default. If *False* is
            provided then the part will not be added to any assembly.
        :type assy: str, :class:`.Assembly`, or bool

        :return: List of :class:`.Rib` instances. Empty list is returned if
            no parts are created.
        :rtype: list

        **Notes**:

        - If both *maxd* and *npts* are provided then *npts* is used as a
          required minimum.

        """
        return create_wing_parts_along_curve('rib', label, wing, curve, geom1,
                                             geom2, maxd, npts, ref_pln, u1,
                                             u2, s1, s2, indx, assy)


class CreateBulkhead(object):
    """
    Bulkhead creator.
    """

    @staticmethod
    def by_sref(label, fuselage, surface_shape, bodies=(), assy=None):
        """
        Create a bulkhead by a basis shape.

        :param str label: Part label.
        :param fuselage: Fuselage to build part in.
        :type fuselage: :class:`.Fuselage`
        :param surface_like surface_shape: The basis shape to define the
            shape of the bulkhead.
        :param iterable bodies: Extra bodies to be used to define the initial
            part shape.
        :param assy: The assembly to store the part in. If none is provided
            then the active assembly will be used by default. If *False* is
            provided then the part will not be added to any assembly.
        :type assy: str, :class:`.Assembly`, or bool

        :return: Bulkhead part.
        :rtype: :class:`.Bulkhead`
        """
        return create_bulkhead_by_sref(label, fuselage, surface_shape,
                                       bodies, assy)


class CreateFloor(object):
    """
    Floor creation.
    """

    @staticmethod
    def by_sref(label, fuselage, surface_shape, bodies=(), assy=None):
        """
        Create a floor by a basis shape.

        :param str label: Part label.
        :param fuselage: Fuselage to build part in.
        :type fuselage: :class:`.Fuselage`
        :param surface_like surface_shape: The basis shape to define the
            shape of the floor.
        :param iterable bodies: Extra bodies to be used to define the initial
            part shape.
        :param assy: The assembly to store the part in. If none is provided
            then the active assembly will be used by default. If *False* is
            provided then the part will not be added to any assembly.
        :type assy: str, :class:`.Assembly`, or bool

        :return: Floor part.
        :rtype: :class:`.Floor`
        """
        return create_floor_by_sref(label, fuselage, surface_shape, bodies,
                                    assy)


class CreateFrame(object):
    """
    Frame creator.
    """

    @staticmethod
    def by_sref(label, fuselage, surface_shape, height, assy=None):
        """
        Create a frame by a basis shape.

        :param str label: Part label.
        :param fuselage: Fuselage to build part in.
        :type fuselage: :class:`.Fuselage`
        :param surface_like surface_shape: The basis shape to define the
            shape of the frame.
        :param float height: Frame height.
        :param assy: The assembly to store the part in. If none is provided
            then the active assembly will be used by default. If *False* is
            provided then the part will not be added to any assembly.
        :type assy: str, :class:`.Assembly`, or bool

        :return: Frame part.
        :rtype: :class:`.Frame`
        """
        return create_frame_by_sref(label, fuselage, surface_shape, height,
                                    assy)

    @staticmethod
    def between_planes(label, fuselage, planes, height, maxd=None, nplns=None,
                       indx=1, assy=None):
        """
        Create planar frames between planes.

        :param str label: Part label.
        :param fuselage: Fuselage to build part in.
        :type fuselage: :class:`.Fuselage`
        :param iterable planes: Iterable of planes that parts will be
            created between. Parts are not created at the given planes.
        :param height: Frame height.
        :param maxd: Maximum allowed spacing.
        :param nplns: Number of planes between the given planes.
        :param int indx: Index to append to label (space delimited). This
            creates unique part labels for multiple parts.
        :param assy: The assembly to store the part in. If none is provided
            then the active assembly will be used by default. If *False* is
            provided then the part will not be added to any assembly.
        :type assy: str, :class:`.Assembly`, or bool

        :return: List of :class:`.Frame` instances. Empty list is returned if
            no parts are created.
        :rtype: list

        **Notes**:

        - If both *maxd* and *nplns* are provided then *nplns* is used as a
          required minimum.

        """
        return create_frames_between_planes(label, fuselage, planes, height,
                                            maxd, nplns, indx, assy)

    @staticmethod
    def at_shapes(label, fuselage, surface_shapes, height, indx=1, assy=None):
        """
        Create frames at shapes.

        :param str label: Part label.
        :param fuselage: Fuselage to build part in.
        :type fuselage: :class:`.Fuselage`
        :param iterable surface_shapes: Basis shapes where parts will be
            created.
        :param height: Frame height.
        :param int indx: Index to append to label (space delimited). This
            creates unique part labels for multiple parts.
        :param assy: The assembly to store the part in. If none is provided
            then the active assembly will be used by default. If *False* is
            provided then the part will not be added to any assembly.
        :type assy: str, :class:`.Assembly`, or bool

        :return: List of :class:`.Frame` instances. Empty list is returned if
            no parts are created.
        :rtype: list
        """
        return create_frames_at_shapes(label, fuselage, surface_shapes, height,
                                       indx, assy)


class CreateSkin(object):
    """
    Skin creator.
    """

    @staticmethod
    def from_solid(label, solid, copy=False, assy=None):
        """
        Create skin from the outer shell of a solid.
        
        :param str label: Part label.
        :param solid: Solid to extract outer shell from.
        :type solid: TopoDS_Solid
        :param bool copy: Option to first copy the shell before creating the
            part.
        :param assy: The assembly to store the part in. If none is provided
            then the active assembly will be used by default. If *False* is
            provided then the part will not be added to any assembly.
        :type assy: str, :class:`.Assembly`, or bool
         
        :return: Skin part.
        :rtype: :class:`.Skin`

        **Notes**:

        - The *copy* option is used to avoid modification of the part
          tolerance during boolean operations involving the body. This
          feature has been altered in new versions of OpenCASCADE and may
          be removed in the future.

        """
        return create_skin_from_solid(label, solid, copy, assy)

    @staticmethod
    def from_body(label, body, copy=False, assy=None):
        """
        Create skin from the outer shell of a body.

        :param str label: Part label.
        :param body: Body to extract outer shell from.
        :type body: :class:`.Body`
        :param bool copy: Option to first copy the shell before creating the
            part.
        :param assy: The assembly to store the part in. If none is provided
            then the active assembly will be used by default. If *False* is
            provided then the part will not be added to any assembly.
        :type assy: str, :class:`.Assembly`, or bool

        :return: Skin part.
        :rtype: :class:`.Skin`

        **Notes**:

        - The *copy* option is used to avoid modification of the part
          tolerance during boolean operations involving the body. This
          feature has been altered in new versions of OpenCASCADE and may
          be removed in the future.

        """
        return create_skin_from_body(label, body, copy, assy)


class CreateStringer(object):
    """
    Stringer creator.
    """

    @staticmethod
    def by_spine(label, surface_part, spine, h, runout_angle=30., assy=None):
        """
        Create a stringer on a surface part by a spine path.

        :param str label: Part label.
        :param surface_part: Surface part to create stringer on.
        :type surface_part: :class:`.SurfacePart`
        :param curve_like spine: Entity defining spine of stringer.
        :param float h: Stringer height.
        :param float runout_angle: Stringer run-out angles (degrees). This
            is the angle between the ti of the stringer and the surface part.
        :param assy: The assembly to store the part in. If none is provided
            then the active assembly will be used by default. If *False* is
            provided then the part will not be added to any assembly.
        :type assy: str, :class:`.Assembly`, or bool

        :return: Stringer part.
        :rtype: :class:`.Stringer`
        """
        return create_stiffener2d_by_wire('stringer', label, surface_part,
                                          spine, h, runout_angle,
                                          assy=assy)

    @staticmethod
    def by_section(label, surface_part, spine_shape, h, runout_angle=30.,
                   assy=None):
        """
        Create a stringer on a surface part by intersection.

        :param str label: Part label.
        :param surface_part: Surface part to create stringer on.
        :type surface_part: :class:`.SurfacePart`
        :param surface_like spine_shape: Shape used to intersect the surface
            and define the spine path.
        :param float h: Stringer height.
        :param float runout_angle: Stringer run-out angles (degrees). This
            is the angle between the ti of the stringer and the surface part.
        :param assy: The assembly to store the part in. If none is provided
            then the active assembly will be used by default. If *False* is
            provided then the part will not be added to any assembly.
        :type assy: str, :class:`.Assembly`, or bool

        :return: Stringer part.
        :rtype: :class:`.Stringer`

        **Notes**:

        - If more than one intersection segment is found then only the
          longest is used. Use *by_sections()* if all intersection segments
          should be used.

        """
        return create_stiffener2d_by_section('stringer', label, surface_part,
                                             spine_shape, h, runout_angle,
                                             assy=assy)

    @staticmethod
    def by_sections(label, surface_part, spine_shape, h, runout_angle=30.,
                    assy=None):
        """
        Create a stringer on a surface part by all intersections.

        :param label: 
        :param surface_part:
        :param spine_shape: 
        :param h: 
        :param runout_angle:
        :param assy:

        :param str label: Part label.
        :param surface_part: Surface part to create stringers on.
        :type surface_part: :class:`.SurfacePart`
        :param surface_like spine_shape: Shape used to intersect the surface
            and define the spine path.
        :param float h: Stringer height.
        :param float runout_angle: Stringer run-out angles (degrees). This
            is the angle between the ti of the stringer and the surface part.
        :param assy: The assembly to store the part in. If none is provided
            then the active assembly will be used by default. If *False* is
            provided then the part will not be added to any assembly.
        :type assy: str, :class:`.Assembly`, or bool

        :return: Stringer parts. Empty list is returned if no parts are
            created.
        :rtype: list
        """
        return create_stiffener2d_by_sections('stringer', label, surface_part,
                                              spine_shape, h, runout_angle,
                                              assy=assy)


class CreateStiffener2D(object):
    """
    Stiffener2D creation.
    """

    @staticmethod
    def by_spine(label, surface_part, spine, h, runout_angle=30.):
        """
        Create a 2-D stiffener on a surface part by a spine path.

        :param str label: Part label.
        :param surface_part: Surface part to create stiffener on.
        :type surface_part: :class:`.SurfacePart`
        :param curve_like spine: Entity defining spine of stiffener.
        :param float h: Stiffener height.
        :param float runout_angle: Stiffener run-out angles (degrees). This
            is the angle between the ti of the stiffener and the surface part.

        :return: Stiffener2D part.
        :rtype: :class:`.Stiffener2D`
        """
        return create_stiffener2d_by_wire('stiffener', label, surface_part,
                                          spine, h, runout_angle)

    @staticmethod
    def by_section(label, surface_part, spine_shape, h, runout_angle=30.):
        """
        Create a stiffener on a surface part by intersection.

        :param str label: Part label.
        :param surface_part: Surface part to create stiffener on.
        :type surface_part: :class:`.SurfacePart`
        :param surface_like spine_shape: Shape used to intersect the surface
            and define the spine path.
        :param float h: Stiffener height.
        :param float runout_angle: Stiffener run-out angles (degrees). This
            is the angle between the ti of the stiffener and the surface part.

        :return: Stiffener2D part.
        :rtype: :class:`.Stiffener2D`

        **Notes**:

        - If more than one intersection segment is found then only the
          longest is used. Use *by_sections()* if all intersection segments
          should be used.

        """
        return create_stiffener2d_by_section('stiffener', label, surface_part,
                                             spine_shape, h, runout_angle)

    @staticmethod
    def by_sections(label, surface_part, spine_shape, h, runout_angle=30.):
        """
        Create a stiffener on a surface part by all intersections.

        :param str label: Part label.
        :param surface_part: Surface part to create stiffeners on.
        :type surface_part: :class:`.SurfacePart`
        :param surface_like spine_shape: Shape used to intersect the surface
            and define the spine path.
        :param float h: Stiffener height.
        :param float runout_angle: Stiffener run-out angles (degrees). This
            is the angle between the ti of the stiffener and the surface part.

        :return: Stiffener2D parts. Empty list is returned if no parts are
            created.
        :rtype: list
        """
        return create_stiffener2d_by_sections('stiffener', label, surface_part,
                                              spine_shape, h, runout_angle)


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
    stringer = CreateStringer()
    stiffener2d = CreateStiffener2D()

    @staticmethod
    def curve_part(label, curve_shape=None, shape1=None, shape2=None,
                   assy=None):
        """
        Create a general curve part.
        
        :param str label: Part label. 
        :param curve_like curve_shape: Part basis shape.
        :param shape_like shape1: Starting shape.
        :param shape_like shape2: Ending shape.
        :param assy: The assembly to store the part in. If none is provided
            then the active assembly will be used by default. If *False* is
            provided then the part will not be added to any assembly.
        :type assy: str, :class:`.Assembly`, or bool
        
        :return: Curve part.
        :rtype: :class:`.CurvePart`
        """
        if curve_shape:
            return create_curve_part(label, curve_shape, assy)
        if shape1 and shape2:
            return create_curve_part_by_section(label, shape1, shape2, assy)
        return None

    @staticmethod
    def surface_part(label, surface_shape, bodies=(), assy=None):
        """
        Create a general surface part.

        :param str label: Part label.
        :param surface_like surface_shape: Part basis shape.
        :param iterable bodies: Extra bodies to be used to define the initial
            part shape.
        :param assy: The assembly to store the part in. If none is provided
            then the active assembly will be used by default. If *False* is
            provided then the part will not be added to any assembly.
        :type assy: str, :class:`.Assembly`, or bool

        :return: Surface part.
        :return: :class:`.SurfacePart`
        """
        return create_surface_part(label, surface_shape, bodies, assy)

    @staticmethod
    def create_stiffener1d(surface_part, label, stiffener):
        """
        Create a 1-D stiffener on the surface part.

        :param surface_part: Surface part to create stiffener on.
        :type surface_part: :class:`SurfacePart`
        :param str label: Part label.
        :param stiffener: Entity that defines the stiffener shape. This
            could be a curve, a curve like shape, or a Stiffener1D instance.
        :type stiffener: curve_like or :class:`.Stiffener1D`

        :return: Stiffener1D part.
        :rtype: :class:`.Stiffener1D`
        """
        return create_stiffener1d(surface_part, stiffener, label)
