"""
We need a way to update geometric properties
via dxtbx objects, without re-instantiating the simulator
These unit tests check our functon update_dxtbx_geoms which allows
one to update those models for beam/detector post-instantiation
"""
from __future__ import division

import numpy as np
from simtbx.nanoBragg import sim_data
from dxtbx.model import Panel

# make a simple detector
det = sim_data.SimData.simple_detector(detector_distance_mm=150, pixelsize_mm=0.1, image_shape=(512, 512))

SIM = sim_data.SimData(use_default_crystal=True)
# snag the beam for later use
B = SIM.beam.nanoBragg_constructor_beam
# set the detector
SIM.detector = det
SIM.instantiate_diffBragg(auto_set_spotscale=True)
D = SIM.D
D.add_diffBragg_spots()
# get the simulated image
img = D.raw_pixels.as_numpy_array()
D.raw_pixels *= 0  # reset

# update the geom to the same thing as it already is (sanity check)
D.update_dxtbx_geoms(det, B, 0)
D.update_dxtbx_geoms(det, B, 0)  # do it a few times to ensure no beam center isnt walking
D.update_dxtbx_geoms(det, B, 0)
D.add_diffBragg_spots()
img2 = D.raw_pixels.as_numpy_array()
# verify the images are the same
assert np.allclose(img, img2)
D.raw_pixels *= 0  # reset the pixels

# make a new detector with shifted detector distance
det2 = sim_data.SimData.simple_detector(detector_distance_mm=151, pixelsize_mm=0.1, image_shape=(512, 512))
D.update_dxtbx_geoms(det2, B, 0)
# simulate the new pattern
D.add_diffBragg_spots()
img3 = D.raw_pixels.as_numpy_array()
D.free_all()  # free c++ objects

# make a new diffBragg instantiating directly using the shifted detector
# e.g. use the nanoBragg dxtbx constructor
SIM = sim_data.SimData(use_default_crystal=True)
SIM.detector = det2
SIM.instantiate_diffBragg(auto_set_spotscale=True)
# make the new pattern and verify its the same as img3
SIM.D.add_diffBragg_spots()
img4 = SIM.D.raw_pixels.as_numpy_array()
assert np.allclose(img3, img4)

# modify the detector dxtbx object directly to produce the same pattern:
# NOTE: this is how we can adjust the detector distance without re-instantiating
node = det[0]
node_d = node.to_dict()
O = node_d["origin"]
distance_shift = -1
O2 = O[0], O[1], O[2] + distance_shift
node_d["origin"] = O2
det[0] = Panel.from_dict(node_d)
assert det == det2
# NOTE: too bad dxtbx doesnt seem to have a set_origin or set_fdet

SIM.D.update_dxtbx_geoms(det, B, 0)
SIM.D.raw_pixels *= 0
SIM.D.add_diffBragg_spots()
img5 = SIM.D.raw_pixels.as_numpy_array()
assert np.allclose(img4, img5)

# sanity check:
assert np.allclose(det[0].get_beam_centre(B.get_s0()), SIM.D.beam_center_mm)

print("OK!")
