from __future__ import print_function

import time

from asap.config import Settings
from asap.fem import MeshData
from asap.geometry import CreateGeom, ProjectGeom
from asap.graphics import Viewer
from asap.io import ImportVSP
from asap.structure import AssemblyData, CreatePart, PartTools
from asap.topology import ShapeTools


def build_wingbox(wing, params):
    """
    Simple 2 spar wing box with mid spar and aux spar.
    """
    # Set units to inch.
    Settings.set_units('in')

    # SETUP -------------------------------------------------------------------
    _default = {'rib spacing': 30.,
                'fspar chord': 0.15,
                'rspar chord': 0.65,
                'root span': 0.05,
                'tip span': 0.925,
                'assy name': 'RH main wingbox',
                'mid spar rib': 10,
                'root rib axes': 'xz',
                'part tol': 0.01,
                'build aux': True,
                'aux rib list': ['2', '5', '8'],
                'aux spar xloc': 0.85}

    for key in _default:
        if key not in params:
            params[key] = _default[key]

    # Unpack params.
    rib_spacing = params['rib spacing']
    fspar_chord = params['fspar chord']
    rspar_chord = params['rspar chord']
    root_span = params['root span']
    tip_span = params['tip span']
    assy_name = params['assy name']
    mid_spar_rib = params['mid spar rib']
    build_aux_spar = params['build aux']
    aux_rib_list = params['aux rib list']
    aux_spar_xloc = params['aux spar xloc']
    Settings.set_part_tol(params['part tol'])

    # BUILD -------------------------------------------------------------------
    AssemblyData.create_assy(assy_name)

    # Front spar
    fspar = CreatePart.spar.by_parameters('front spar', wing,
                                          fspar_chord, root_span,
                                          fspar_chord, tip_span)

    # Rear spar
    rspar = CreatePart.spar.by_parameters('rear spar', wing,
                                          rspar_chord, root_span,
                                          rspar_chord, tip_span)

    # Root rib
    p1 = fspar.p1
    if not build_aux_spar:
        p2 = rspar.p1
    else:
        p2 = wing.eval(aux_spar_xloc, root_span)
    root = CreatePart.rib.by_points('root rib', wing, p1, p2)

    # Tip rib
    p1 = fspar.p2
    p2 = rspar.p2
    tip = CreatePart.rib.by_points('tip rib', wing, p1, p2)

    # Generate points along rear spar and project to front spar to define ribs.
    prear = rspar.spaced_points(rib_spacing, rib_spacing, -rib_spacing)
    pfront = [p.copy() for p in prear]
    fspar.points_to_cref(pfront)
    i = 1
    ribs = []
    for pf, pr in zip(pfront, prear):
        if pf.is_equal(pr):
            continue
        name = ' '.join(['rib', str(i)])
        rib = CreatePart.rib.by_points(name, wing, pf, pr)
        if not rib:
            continue
        ribs.append(rib)
        i += 1

    # Build a mid spar.
    mspar = None
    if mid_spar_rib > 0:
        u = root.local_to_global_u(0.4)
        p1 = root.eval_cref(u)
        rib = ribs[mid_spar_rib - 1]
        u = rib.local_to_global_u(0.5)
        p2 = rib.eval_cref(u)
        mspar = CreatePart.spar.by_points('mid spar', wing, p1, p2)

    # Aux spar.
    aspar = None
    if build_aux_spar:
        p1 = root.p2
        u = rspar.local_to_global_u(0.25)
        p2 = rspar.eval_cref(u)
        # Find nearest rib point and set equal to aux spar.
        pnts = [rib.p2 for rib in ribs]
        p2 = CreateGeom.nearest_point(p2, pnts)
        indx = pnts.index(p2)
        rib = ribs[indx]
        # Use intersection of the rib and rear spar to define plane for aux
        # spar.
        sref = ShapeTools.plane_from_section(rspar, rib, p1)
        aspar = CreatePart.spar.by_points('aux spar', wing, p1, p2, sref)

    # Build ribs from root rib to front spar.
    root_ribs = []
    if mspar:
        # Fwd of mid spar
        u2 = root.invert_cref(mspar.p1)
        prib = root.spaced_points(3, u2=u2)[1:-1]
        pfront = [p.copy() for p in prib]
        fspar.points_to_cref(pfront)
        for pf, pr in zip(pfront, prib):
            if pf.is_equal(pr):
                continue
            name = ' '.join(['rib', str(i)])
            rib = CreatePart.rib.by_points(name, wing, pf, pr)
            if not rib:
                continue
            ribs.append(rib)
            root_ribs.append(rib)
            i += 1

        # Aft of mid spar
        u1 = root.invert_cref(mspar.p1)
        u2 = root.invert_cref(rspar.p1)
        prib = root.spaced_points(3, u1=u1, u2=u2)[1:-1]
        pfront = [p.copy() for p in prib]
        fspar.points_to_cref(pfront)
        for pf, pr in zip(pfront, prib):
            if pf.is_equal(pr):
                continue
            name = ' '.join(['rib', str(i)])
            rib = CreatePart.rib.by_points(name, wing, pf, pr)
            if not rib:
                continue
            ribs.append(rib)
            root_ribs.append(rib)
            i += 1

    # Construction geom for center structure.
    root_chord = wing.extract_curve((0, 0), (1, 0))

    # Front center spar.
    p2 = fspar.p1
    p1 = p2.copy()
    ProjectGeom.point_to_geom(p1, root_chord, True)
    CreatePart.spar.by_points('fc spar', wing, p1, p2)

    # Rear center spar.
    p2 = rspar.p1
    p1 = p2.copy()
    ProjectGeom.point_to_geom(p1, root_chord, True)
    CreatePart.spar.by_points('rc spar', wing, p1, p2)

    # Mid center spar
    if mid_spar_rib > 0:
        p2 = mspar.p1
        p1 = p2.copy()
        ProjectGeom.point_to_geom(p1, root_chord, True)
        CreatePart.spar.by_points('center mid spar', wing, p1, p2)

    # Center spar at each root rib intersection
    i = 1
    for rib in root_ribs:
        p2 = rib.p2
        p1 = p2.copy()
        ProjectGeom.point_to_geom(p1, root_chord, True)
        name = ' '.join(['center spar', str(i)])
        CreatePart.spar.by_points(name, wing, p1, p2)
        i += 1

    # Aux ribs
    if build_aux_spar:
        assy = AssemblyData.get_active()
        aux_rib_id = 1
        for rib_id in aux_rib_list:
            rib_name = ' '.join(['rib', rib_id])
            rib = assy.get_part(rib_name)
            if not rib:
                continue
            # Since the structure is not joined yet, intersect the rear spar
            # and rib shapes to find the edge(s). Use this edge to define a
            # plane so the aux ribs will line up with the main ribs.
            edges = ShapeTools.bsection(rspar, rib, 'e')
            if not edges:
                continue
            edge = edges[0]
            pnts = ShapeTools.points_along_edge(edge, 3)
            p1 = rib.p2
            p2 = p1.copy()
            aspar.points_to_cref([p2])
            sref = CreateGeom.fit_plane([p2] + pnts)
            if not sref:
                continue
            aux_rib_name = ' '.join(['aux rib', str(aux_rib_id)])
            CreatePart.rib.by_points(aux_rib_name, wing, p1, p2, sref)
            aux_rib_id += 1

    # JOIN --------------------------------------------------------------------
    # Fuse internal structure and discard faces
    internal_parts = AssemblyData.get_parts()
    PartTools.fuse_wing_parts(internal_parts)
    for part in internal_parts:
        part.discard()

    # SKIN --------------------------------------------------------------------
    skin = CreatePart.surface_part('wing skin', wing.shell)

    # Join the wing skin and internal structure
    skin.fuse(*internal_parts)

    # Discard faces touching reference surface.
    skin.discard(wing.sref)

    # Fix skin since it's not a single shell anymore, but a compound of two
    # shells (upper and lower skin).
    skin.fix()

    # Viewing
    skin.set_transparency(0.5)

    # MESH --------------------------------------------------------------------
    # Initialize
    shape_to_mesh = AssemblyData.prepare_shape_to_mesh()
    MeshData.create_mesh('wing-box mesh', shape_to_mesh)

    # Use a single global hypothesis based on local length.
    MeshData.hypotheses.create_netgen_simple_2d('netgen hypo', 4.)
    MeshData.hypotheses.create_netgen_algo_2d('netgen algo')
    MeshData.add_hypothesis('netgen hypo')
    MeshData.add_hypothesis('netgen algo')

    # Compute the mesh
    mesh_start = time.time()
    print('Computing mesh...')
    status = MeshData.compute_mesh()
    if not status:
        print('Failed to compute mesh')
    else:
        print('Meshing complete in ', time.time() - mesh_start, ' seconds.')

    # Uncomment this to export STEP file.
    # from asap.io import StepExport
    # step = StepExport()
    # step.transfer(shape_to_mesh)
    # step.write('wingbox.step')

    return AssemblyData.get_active()


if __name__ == '__main__':
    start = time.time()

    # Import model
    fname = r'..\models\777-200LR.stp'
    ImportVSP.step_file(fname)
    wing_in = ImportVSP.get_body('Wing')

    # Build wing box
    assy = build_wingbox(wing_in, {})

    print('Complete in ', time.time() - start, ' seconds.')

    Viewer.add_items(*assy.parts)
    Viewer.add_meshes(MeshData.get_active())
    Viewer.show()
