"""
Simple test of multi-panel
detector simulation
"""
from __future__ import absolute_import, division, print_function
import numpy as np
import copy

import simtbx.nanoBragg
from scitbx.matrix import sqr
from dxtbx.model.beam import BeamFactory
from dxtbx.model.crystal import CrystalFactory
from dxtbx.model.detector import DetectorFactory, Detector, Panel
from six.moves import range


# dxtbx beam model description
beam_descr = {'direction': (0.0, 0.0, 1.0),
             'divergence': 0.0,
             'flux': 1e12,
             'polarization_fraction': 1.,
             'polarization_normal': (0.0, 1.0, 0.0),
             'sigma_divergence': 0.0,
             'transmission': 1.0,
             'wavelength': 1.4}

# dxtbx crystal description
cryst_descr = {'__id__': 'crystal',
              'real_space_a': (79, 0, 0),
              'real_space_b': (0, 79, 0),
              'real_space_c': (0, 0, 38),
              'space_group_hall_symbol': '-P 4 2'}

# monolithic camera description
det_descr = {'panels':
               [{'fast_axis': (-1.0, 0.0, 0.0),
                 'gain': 1.0,
                 'identifier': '',
                 'image_size': (512, 256),
                 'mask': [],
                 'material': '',
                 'mu': 0.0,
                 'name': 'Panel',
                 'origin': (25.6, -12.8, -100),
                 'pedestal': 0.0,
                 'pixel_size': (0.1, 0.1),
                 'px_mm_strategy': {'type': 'SimplePxMmStrategy'},
                 'raw_image_offset': (0, 0),
                 'slow_axis': (0.0, 1.0, 0.0),
                 'thickness': 0.0,
                 'trusted_range': (0.0, 65536.0),
                 'type': ''}]}


beam = BeamFactory.from_dict(beam_descr)
whole_det = DetectorFactory.from_dict(det_descr)
cryst = CrystalFactory.from_dict(cryst_descr)

# ---------

SIM = simtbx.nanoBragg.nanoBragg( whole_det, beam)
SIM.Amatrix = sqr(cryst.get_A()).transpose().elems
SIM.default_F = 1
SIM.F000 = 10
SIM.oversample=2
SIM.Ncells_abc = (5,5,5)
SIM.interpolate=0
SIM.show_params()
SIM.add_nanoBragg_spots()
SIM.raw_pixels *= 2000
wholepix = SIM.raw_pixels.as_numpy_array()

# ---------

# now make 2 smaller detectors
# that equal the whole detector above
# then append them as a dials multi panel
origin0 = (25.6, -12.8, -100)
origin1 = (25.6, 0, -100)
panel_size =(512, 128)

pan0 = {'fast_axis': (-1.0, 0.0, 0.0),
        'gain': 1.0,
        'identifier': '',
        'image_size': (512, 128),
        'mask': [],
        'material': '',
        'mu': 0.0,
        'name': 'panel1',
        'origin': origin0,
        'pedestal': 0.0,
        'pixel_size': (0.1, 0.1),
        'px_mm_strategy': {'type': 'SimplePxMmStrategy'},
        'raw_image_offset': (0, 0),
        'slow_axis': (0.0, 1.0, 0.0),
        'thickness': 0.0,
        'trusted_range': (0.0, 65536.0),
        'type': ''}

pan1 = copy.deepcopy(pan0)
pan1["origin"] = origin1

# combine the panels into a multi-panel detector
parts_det = Detector()
parts_det.add_panel( Panel.from_dict(pan0))
parts_det.add_panel( Panel.from_dict(pan1))

SIMs = {}
for i_pan in range(2):
  SIM = simtbx.nanoBragg.nanoBragg( parts_det, beam, panel_id=i_pan)
  SIM.Amatrix = sqr(cryst.get_A()).transpose().elems
  SIM.default_F = 1
  SIM.F000 = 10
  SIM.oversample=2
  SIM.Ncells_abc = (5,5,5)
  SIM.interpolate=0
  SIM.add_nanoBragg_spots()
  SIM.raw_pixels *= 2000
  SIMs[i_pan] = SIM

# check comparison with whole camera
# ----------------------------------
# combine the two panels:
pix0 = SIMs[0].raw_pixels.as_numpy_array()
pix1 = SIMs[1].raw_pixels.as_numpy_array()
pix01 = np.vstack((pix0,pix1))

assert( np.allclose(wholepix, pix01))

if __name__=="__main__":
 print("OK")
