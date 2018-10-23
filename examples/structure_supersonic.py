from afem.config import Settings
from afem.exchange import ImportVSP, StepWrite
from afem.geometry import *
from afem.graphics import Viewer
from afem.structure import *
from afem.topology import *


def build(wing, fuselage):
    """
    Build wing and fuselage structure. Parameters can be exposed as inputs to
    this method.

    :param afem.oml.entities.Body wing: The wing.
    :param afem.oml.entities.Body fuselage: The fuselage.
    """

    # Create to assemblies for storing parts in
    wa = GroupAPI.create_group('wing group')
    fa = GroupAPI.create_group('fuselage group')

    # Extract a bounding box of the fuselage enlarged by 6 units
    fuselage_bbox = fuselage.bbox(6.)

    # Get a list of v-direction knots of the wing reference surface. For
    # OpenVSP, each v-direction knot in the reference surface will coincide
    # with a cross section in the model. These parameters are used to define
    # structure at wing cross sections. The "vknots" of the wing reference
    # surface may not always be available, but since a NurbsSurface is used for
    # OpenVSP translation, they are in this use case.
    vknots = wing.sref.vknots

    # Center rib using an xz-plane
    xz_plane = PlaneByAxes().plane
    RibByShape('center rib', xz_plane, wing, wa)

    # Ribs 1, 2, and 3 are at specific key points in the wing and are defined
    # in the xz-plane. For this model, using the xz-plane makes the ribs align
    # vertically with the wing cross sections and avoids narrow faces.

    # Rib 1 is at the maximum y-location of the fuselage bounding box.
    # Remember that the bounding box has already been enlarged by 6 units.
    p = PointByXYZ(0., fuselage_bbox.ymax, 0.).point
    sref = PlaneByAxes(p, 'xz').plane
    rib1 = RibByShape('rib 1', sref, wing, wa).part

    # Rib 2 is at the third cross section of the wing. The v-direction
    # parameter is vknots[2] since Python indexing starts at 0.
    p = wing.sref.eval(0., vknots[2])
    sref = PlaneByAxes(p, 'xz').plane
    rib2 = RibByShape('rib 2', sref, wing, wa).part

    # Rib 3 is at the tip of the wing which is found using the parameters
    # (0., 1.) of the wing reference surface.
    p = wing.sref.eval(0., 1.)
    sref = PlaneByAxes(p, 'xz').plane
    rib3 = RibByShape('rib 3', sref, wing, wa).part

    # Spars are defined in the wing using a front and rear spar location and
    # then evenly spaced number of spars between them. At each location there
    # are three distinct spars:
    #
    #   1) The center spars which are in the yz-plane between the center rib
    #      and Rib 1.
    #   2) The inboard spars which are between Rib 1 and Rib 2.
    #   3) The outboard spars which are between Rib 2 and Rib 3.
    #
    # To maintain proper congruent edges, the orientation of the inboard and
    # outboard spars is defined using the intersection of the previous parts
    # via the PlaneByIntersectingShapes tool. If not used, spar interfaces may
    # be at different angles and result in invalid joints.
    #
    # Bulkheads that interface with the center spars are also created along the
    # way.

    # Front center spar at a relative distance along the center rib
    p1 = rib1.point_from_parameter(0.15, is_rel=True)
    sref = PlaneByAxes(p1, 'yz').plane
    fc_spar = SparBetweenShapes('center fspar', xz_plane, rib1.shape, wing,
                                sref, wa).part

    # Front inboard spar
    p2 = rib2.point_from_parameter(0.15, is_rel=True)
    sref = PlaneByIntersectingShapes(fc_spar.shape, rib1.shape, p2).plane
    inbd_fspar = SparByPoints('inbd fspar', fc_spar.cref.p2, p2, wing, sref,
                              wa).part

    # Front outboard spar
    p3 = rib3.point_from_parameter(0.15, is_rel=True)
    sref = PlaneByIntersectingShapes(inbd_fspar.shape, rib2.shape, p3).plane
    outbd_fspar = SparByPoints('outbd fspar', inbd_fspar.cref.p2, p3, wing,
                               sref, wa).part

    # Front spar bulkhead
    BulkheadByShape('fspar bh', fc_spar.plane, fuselage, fa)

    # Rear center spar at a relative distance along the center rib
    p1 = rib1.point_from_parameter(0.92, is_rel=True)
    sref = PlaneByAxes(p1, 'yz').plane
    rc_spar = SparBetweenShapes('center rspar', xz_plane, rib1.shape, wing,
                                sref, wa).part

    # The rear spar is a single spar between the rear center spar and Rib 3.
    # The orientation is defined by the rear center spar at a relative location
    # along Rib 3.
    p3 = rib3.point_from_parameter(0.80, is_rel=True)
    sref = PlaneByIntersectingShapes(rc_spar.shape, rib1.shape, p3).plane
    rspar = SparByPoints('rspar', rc_spar.cref.p2, p3, wing, sref, wa).part

    # Rear spar bulkhead
    BulkheadByShape('rspar bh', rc_spar.plane, fuselage, fa)

    # Generate evenly spaced points between the front and rear spars to define
    # intermediate spars. Start and end parameters are found along the key ribs
    # and used in the PointsAlongCurveByNumber tool. Only the interior points
    # are used so overlapping spars are not defined at the ends of the ribs.

    # Points for center spars along Rib 1
    u1 = rib1.cref.invert(fc_spar.cref.p2)
    u2 = rib1.cref.invert(rc_spar.cref.p2)
    pnts1 = PointsAlongCurveByNumber(rib1.cref, 6, u1, u2).interior_points

    # Points for inboard spars along Rib 2
    u1 = rib2.cref.invert(inbd_fspar.cref.p2)
    pi = IntersectCurveCurve(rib2.cref, rspar.cref).point(1)
    u2 = rib2.cref.invert(pi)
    pnts2 = PointsAlongCurveByNumber(rib2.cref, 6, u1, u2).interior_points

    # Points for outboard spars along Rib 3
    u1 = rib3.cref.invert(outbd_fspar.cref.p2)
    u2 = rib3.cref.invert(rspar.cref.p2)
    pnts3 = PointsAlongCurveByNumber(rib3.cref, 6, u1, u2).interior_points

    # Intermediate spars and bulkheads
    for p1, p2, p3 in zip(pnts1, pnts2, pnts3):
        # Bulkhead and center spar location
        sref = PlaneByAxes(p1, 'yz').plane
        BulkheadByShape('bh', sref, fuselage, fa)

        # Center spar
        spar1 = SparBetweenShapes('spar', xz_plane, rib1.shape, wing, sref,
                                  wa).part

        # Inboard spar
        sref = PlaneByIntersectingShapes(spar1.shape, rib1.shape, p2).plane
        spar2 = SparByPoints('spar', p1, p2, wing, sref, wa).part

        # Outboard spar
        sref = PlaneByIntersectingShapes(spar2.shape, rib2.shape, p3).plane
        SparByPoints('spar', p2, p3, wing, sref, wa)

    # Inboard ribs by number between Rib 1 and Rib 2. First define the planes
    # then define the ribs.
    plns = PlanesBetweenPlanesByNumber(rib1.plane, rib2.plane, 1).planes
    for pln in plns:
        RibByShape('inbd rib', pln, wing, wa)

    # Outboard ribs by distance between Rib 2 and Rib 3
    plns = PlanesBetweenPlanesByDistance(rib2.plane, rib3.plane, 60., 60.,
                                         -60.).planes
    for pln in plns:
        RibByShape('rib', pln, wing, wa)

    # Forward bulkhead at the second wing cross section in the yz-plane. The
    # origin of the plane is the leading edge (u=0.) at the second cross
    # section (vknots[1]).
    p = wing.sref.eval(0., vknots[1])
    sref = PlaneByAxes(p, 'yz').plane
    fbh = BulkheadByShape('forward bh', sref, fuselage, fa).part

    # Forward frames between an offset distance from the nose to the forward
    # bulkhead. Specify the spacing and the height.
    yz_plane = PlaneByAxes(axes='yz').plane
    frames = FramesBetweenPlanesByDistance('frame', yz_plane, fbh.plane, 60.,
                                           fuselage, 3., 84., group=fa).parts

    # Aft bulkheads between the rear spar bulkhead and a distance offset from
    # the aft end of the fuselage. Get the maximum x-distance (i.e., length of
    # fuselage) by again finding a fuselage bounding box but do not enlarge it.
    bbox = fuselage.bbox()
    aft_pln = PlaneByAxes((bbox.xmax, 0., 0.), 'yz').plane
    pln1 = PlaneByAxes(wing.sref.eval(1., 0.), 'yz').plane
    plns = PlanesBetweenPlanesByDistance(pln1, aft_pln, 60., 24., -24.).planes
    last_plane = plns[-1]
    for pln in plns:
        BulkheadByShape('bh', pln, fuselage, fa)

    # Gather lists of internal parts for the wing and fuselage structure
    internal_wing_parts = wa.get_parts()
    internal_fuse_parts = fa.get_parts()

    # Show the internal wing parts before fusing and discarding. The shapes of
    # parts may extend beyond their intended shape and must be adjusted. This
    # is because defining the initial shape of the part is done using a simple
    # Boolean common operation between the body (a solid) and the part basis
    # shape (using a plane or a face that extends beyond the body). Since this
    # Boolean operation results in shapes that do not yet represent the
    # intended part, we must devise a way to (at least semi-automatically)
    # refine the shape.
    gui.add(*internal_wing_parts)
    gui.start()
    gui.clear()

    # When a wing part is created and if a wing reference surface is available,
    # a part "reference curve" gets created and attached to the part. This
    # curve is the result of the intersection between the part basis shape and
    # the wing reference surface and is trimmed between the actual start and
    # end points of the part. These curves can be used to test for possible
    # intersection with other parts since they are trimmed. Also, faces of
    # shapes that may be beyond the intended part boundaries after trimming can
    # be identified and discarded with these curves. That is what the following
    # tools do.

    # Show the wing reference surface and the underlying part reference curves
    gui.add(wing.sref)
    for part in internal_wing_parts:
        gui.add(part.cref)
    gui.start()
    gui.clear()

    # Fuse the interfacing wing parts together first using their reference
    # curves. Then discard faces of the parts using the reference curves again.
    FuseSurfacePartsByCref(internal_wing_parts)
    DiscardByCref(internal_wing_parts)

    # Show the parts again to see the difference
    gui.add(*internal_wing_parts)
    gui.start()
    gui.clear()

    # Fuselage skin using the shell of the fuselage body
    fskin = SkinByBody('fuselage skin', fuselage, group=fa).part

    # Wing skin using the shell of the wing body
    wskin = SkinByBody('wing skin', wing, group=wa).part

    # Discard the "caps" of the wing skin since ribs have been defined there
    # and we don't want overlapping structure. For this, define an iso-curve
    # within the wing reference surface from root to tip. For this case, the
    # only faces that should be near this curve at the caps of the wing skin.
    # Use this curve and a minimum distance to discard these faces.
    gui.add(wskin)
    gui.start()

    cref = wing.sref.u_iso(0.5)
    wskin.discard_by_dmin(cref, 0.01)

    # Show after discarding
    gui.start()
    gui.clear()

    # Fuse wing skin and internal parts. This will trim the skin with the parts
    # and update all the shapes of each part.
    wskin.fuse(*internal_wing_parts)

    # Fuse fuselage skin and internal parts
    fskin.fuse(*internal_fuse_parts)

    # After some initial experimentation, it was found that the faces near the
    # forward and aft parts of the fuselage were causing issues in the meshing
    # process. These faces do not appear to be well-suited for CAD operations
    # because they look very small, narrow, and have sharp edges. For now,
    # discard these faces since they should have little effect on the global
    # behavior of the structural model.

    # Discard fuselage skin forward of first frame and aft of last bulkhead.
    # In this case solid boxes are defined that contain the faces that are to
    # be discarded. Any face that has a centroid inside the given solid will be
    # removed.
    gui.add(fskin)
    gui.start()

    # Forward faces
    pln = frames[0].plane
    face = FaceByPlane(pln, -1000, 1000, -1000, 1000).face
    box = SolidByDrag(face, (-1000, 0, 0)).solid
    fskin.discard_by_solid(box)

    # Aft faces
    face = FaceByPlane(last_plane, -1000, 1000, -1000, 1000).face
    box = SolidByDrag(face, (1000, 0, 0)).solid
    fskin.discard_by_solid(box)

    # Show after
    gui.start()
    gui.clear()

    # The fuselage structure to this point is defined entirely in the fuselage
    # body. If a "half-model" is desired, one option is to simply the cut the
    # fuselage structure with a box along the centerline as shown below.

    # Define a box along the center of the aircraft and cut all the parts with
    # it
    face = FaceByPlane(xz_plane, fuselage_bbox.xmin, fuselage_bbox.xmax,
                       fuselage_bbox.zmin, fuselage_bbox.zmax).face
    box = SolidByDrag(face, (0., -500., 0.)).solid
    CutParts(fa.get_parts(), box)

    # Show after
    gui.add(*fa.get_parts())
    gui.start()
    gui.clear()

    # During all the Boolean operations, it's possible that some shapes become
    # corrupt (poor tolerances, invalid topology, etc.) and need repaired
    # before further processing. The tool below is an automated approach for
    # fixing general issues in the shapes of each part in the context of an
    # assembly. Fixing part shapes in the context of their assembly is
    # important since the tool can look at the topology of the overall shape
    # rather than just the part shape by itself.
    FixGroup(wa)
    FixGroup(fa)

    # Fuse assemblies together. This tool will build two separate compounds
    # from the parts of each assembly and then do a Boolean fuse operation.
    # The shapes of each part are then updated.
    FuseGroups([wa, fa])

    # Fix master assembly, hopefully addressing global issues after a large
    # Boolean operation like fusing assemblies.
    master = GroupAPI.get_master()
    FixGroup(master)

    # Set some viewing options and show the model
    wskin.set_transparency(0.5)
    fskin.set_transparency(0.5)
    gui.add(*wa.get_parts())
    gui.add(*fa.get_parts())
    gui.start()

    # For now, just throw a global unstructured quad-dominant mesh with a
    # constance size at the whole thing. Local mesh control is possible but
    # not implemented for this example. Automated and streamlined tools and
    # methods are currently being developed to make the meshing process more
    # user friendly and intuitive.

    # Retrieve the global shape to mesh
    the_mesh = MeshVehicle(4.)

    # Compute the mesh
    print('Computing mesh...')
    the_mesh.compute()

    # View the mesh. Adjust these transparencies if you want to see the
    # internal mesh better.
    gui.clear()
    gui.add(the_mesh)
    gui.start()
    gui.clear()

    # Export the shape to a STEP file
    step = StepWrite('AP203', 'in')
    step.transfer(the_mesh.shape)
    step.write('supersonic.step')

    # Export the mesh (nodes and elements) to a bulk data file
    the_mesh.export_nastran('supersonic.bdf')


if __name__ == '__main__':
    # Set units to inch and log to console
    Settings.set_units('in')
    Settings.log_to_console()

    # Import OpenVSP model
    fname = r'..\models\supersonic.stp'
    vsp_import = ImportVSP(fname)

    # A few notes about importing an OpenVSP STEP model:
    #   1) If the model comes from the modified version that includes metadata
    #      then the OpenVSP components can be identified by type and retrieved
    #      by name. Otherwise the import methods just tries to make solid
    #      bodies from the components and gives them a generic name.
    #   2) The import process attempts to sew the faces of the OpenVSP
    #      components together to form solids.
    #   3) If a surface is found to be planar, it is replaced with a plane
    #      before sewing. This helps eliminate some degenerated edges.
    #   4) Faces sharing the same domain are unified. If flat wing caps are
    #      used, this usually results in a single face for the entire cap
    #      rather than two faces split between the upper and lower surface.

    # View the model
    gui = Viewer()
    gui.add(*vsp_import.all_bodies)
    gui.start()
    gui.clear()

    # Retrieve relative components by name and set transparency for viewing
    wing_ = vsp_import['wing']
    fuse_ = vsp_import['fuse']
    wing_.set_transparency(0.5)
    fuse_.set_transparency(0.5)

    # OpenVSP 3.5 was modified by Laughlin Research to construct and export
    # metadata and reference geometry in the STEP file. For a wing component,
    # the reference geometry includes a surface that is lofted through the
    # chord lines at each wing station. This surface is used to define
    # structure in terms of percent chord and/or semispan. Currently, the
    # u-direction of the surface is in the chordwise and the v-direction is
    # spanwise. Show this surface.
    gui.add(wing_)
    gui.add(wing_.sref)
    gui.start()
    gui.clear()

    # Build the structural model
    build(wing_, fuse_)
