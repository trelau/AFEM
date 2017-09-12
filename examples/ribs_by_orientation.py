from __future__ import print_function

import time

from afem.config import Settings
from afem.geometry import *
from afem.graphics import Viewer
from afem.io import ImportVSP
from afem.structure import *

Settings.log_to_console(True)


def build_wingbox(wing, params):
    # Set units to inch.
    Settings.set_units('in')

    # SETUP -------------------------------------------------------------------
    _default = {'rib spacing': 72.,
                'fspar chord': 0.15,
                'rspar chord': 0.65,
                'root span': 0.05,
                'tip span': 0.925,
                'assy name': 'RH main wingbox'}

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

    # BUILD -------------------------------------------------------------------
    AssemblyAPI.create_assy(assy_name)

    # Front spar
    fspar = SparByParameters('front spar', fspar_chord, root_span,
                             fspar_chord, tip_span, wing).spar

    # Rear spar
    rspar = SparByParameters('rear spar', rspar_chord, root_span,
                             rspar_chord, tip_span, wing).spar

    # Root rib
    p1 = fspar.p1
    p2 = rspar.p1
    root = RibByPoints('root rib', p1, p2, wing).rib

    # Tip rib
    p1 = fspar.p2
    p2 = rspar.p2
    tip = RibByPoints('tip rib', p1, p2, wing).rib

    # Locally rotate a plane for a rib
    pln = PlaneByCurveAndSurface(rspar.cref, wing.sref, 1000.).plane
    pln.rotate_x(45.)
    RibBySurface('rotated rib', pln, wing)

    # # Ribs along spar using global orientation
    # ref_pln = PlaneByOrientation(gamma=-45.).plane
    # RibsAlongCurveByDistance('rib', rspar.cref, rib_spacing, fspar, rspar,
    #                          wing, d1=rib_spacing, d2=-rib_spacing,
    #                          ref_pln=ref_pln)

    # Ribs along spar using locally rotated planes.
    RibsAlongCurveAndSurfaceByDistance('rib', rspar.cref, wing.sref,
                                       rib_spacing, fspar, rspar, wing,
                                       d1=rib_spacing, d2=-rib_spacing,
                                       rot_x=30., rot_y=30.)

    # JOIN --------------------------------------------------------------------
    # Fuse internal structure and discard faces
    internal_parts = AssemblyAPI.get_parts(order=True)
    FuseSurfacePartsByCref(internal_parts)
    DiscardByCref(internal_parts)

    return AssemblyAPI.get_active()


if __name__ == '__main__':
    start = time.time()

    # Import model
    fname = r'..\models\777-200LR.stp'
    ImportVSP.step_file(fname)
    wing_in = ImportVSP.get_body('Wing')

    # Build wing box
    assy = build_wingbox(wing_in, {})

    print('Complete in ', time.time() - start, ' seconds.')

    Viewer.add(*assy.parts)

    Viewer.show()
