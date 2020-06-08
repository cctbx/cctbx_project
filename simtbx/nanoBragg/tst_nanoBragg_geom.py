"""
Tests nanoBragg geometry

usage:
libtbx.python tst_nanoBragg_geom.py --plot  # PASSES
libtbx.python tst_nanoBragg_geom.py --plot --rotations 0.1 0.1  # FAILS
libtbx.python tst_nanoBragg_geom.py --plot --rotations 0.1 0.1  --updatecenter # PASSES

Determines where a Bragg peak should be with sub-pixel accuracy
and then simulates an image and verifie that the centroid is within
0.25 pixels of its predicted spot
"""
from __future__ import absolute_import, division, print_function

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--plot", action="store_true", help="display the simulated image and spots")
parser.add_argument("--rotations", nargs=2, default=[0, 0],
                    help="two rotation angles (degrees) for tilting the detector panel (pitch/yaw)", type=float)
parser.add_argument("--updatecenter", action="store_true", help="update nanoBragg center after instantiation")
args = parser.parse_args()

import numpy as np
import simtbx.nanoBragg
from scitbx.matrix import sqr
from dxtbx.model.beam import BeamFactory
from dxtbx.model.crystal import CrystalFactory
from dxtbx.model.detector import DetectorFactory
from simtbx.nanoBragg import shapetype

# MAKE DXTBX BEAM
# dxtbx beam model description
beam_descr = {'direction': (0.0, 0.0, 1.0),
             'divergence': 0.0,
             'flux': 1e12,
             'polarization_fraction': 1.,
             'polarization_normal': (0.0, 1.0, 0.0),
             'sigma_divergence': 0.0,
             'transmission': 1.0,
             'wavelength': 1.3}
beam = BeamFactory.from_dict(beam_descr)


# MAKE DXTBX CRYSTAL
a,b,c = 80,60,40
# dxtbx crystal description
cryst_descr = {'__id__': 'crystal',
              'real_space_a': (a, 0, 0),
              'real_space_b': (0, b, 0),
              'real_space_c': (0, 0, c),
              'space_group_hall_symbol': ' P 2 2ab'}
cryst = CrystalFactory.from_dict(cryst_descr)


# MAKE DXTBX DETECTOR
# monolithic camera description
det_descr = {'panels':
               [{'fast_axis': (-1.0, 0.0, 0.0),
                 'gain': 1.0,
                 'identifier': '',
                 'image_size': (512, 512),
                 'mask': [],
                 'material': '',
                 'mu': 0.0,
                 'name': 'Panel',
                 'origin': (51.2, -51.2, -100),
                 'pedestal': 0.0,
                 'pixel_size': (0.1, 0.1),
                 'px_mm_strategy': {'type': 'SimplePxMmStrategy'},
                 'raw_image_offset': (0, 0),
                 'slow_axis': (0.0, 1.0, 0.0),
                 'thickness': 0.0,
                 'trusted_range': (0.0, 65536.0),
                 'type': ''}]}
whole_det = DetectorFactory.from_dict(det_descr)


# rotate the panel
panel = whole_det[0]
pitch, yaw = args.rotations
panel.rotate_around_origin(axis=(1,0,0), angle=pitch, deg=True)
panel.rotate_around_origin(axis=(0,1,0), angle=yaw, deg=True)

if args.plot:
    print("PLotting the detector model")
    from dxtbx.model import ExperimentList, Experiment
    import os
    Exper = Experiment()
    Exper.detector = whole_det
    Exper.beam = beam
    Elist = ExperimentList()
    Elist.append(Exper)
    Elist_file = "_tst_nanoBragg_geom_whole_detector_test.expt"
    Elist.as_json(Elist_file)
    os.system("dxtbx.plot_detector_models %s" % Elist_file)

# ---------

##
# find the Qx, Qy, Qz of the 10,-9,1 Bragg spot
a,b,c,_,_,_ = cryst.get_unit_cell().parameters()
Amat = [[a,0,0], [0,b,0], [0,0,c]]
Amat_inv = np.linalg.inv(Amat)
h,k,l = 10, -9, 1
qx, qy, qz = np.dot( Amat_inv, [h,k,l])

## given the detector and beam models find the q of every pixel
## on a oversampled detector grid (oversample by a factor of 10)
## this will give us the expected centroid of the spot in pixel units
panel = whole_det[0]
fx,fy,fz = panel.get_fast_axis()
sx,sy,sz = panel.get_slow_axis()
ox,oy,oz = panel.get_origin()
pixFast, pixSlow = panel.get_pixel_size()
assert pixFast == pixSlow
pix = pixFast

fast_dim, slow_dim = panel.get_image_size()

# pixel oversample factor
oversamp = 10

slow_coor, fast_coor = np.indices((slow_dim*oversamp, fast_dim*oversamp))
# NOTE: slow/fast_coor is the slow/fast coordinate of every pixel aranged on a 2D array
# e.g. slow_coor =
#  [0 0 0 0 ..
#   1 1 1 1 ..
#   2 2 2 2 ..

# fast_coor =
#  [ 0 1 2 3
#    0 1 2 3
#    0 1 2 3


# get the components of the scattering vector
kx = ox + fx*fast_coor*pix/oversamp + sx*slow_coor*pix/oversamp
ky = oy + fy*fast_coor*pix/oversamp + sy*slow_coor*pix/oversamp
kz = oz + fz*fast_coor*pix/oversamp + sz*slow_coor*pix/oversamp

knorm = np.sqrt(kx**2 + ky**2 + kz**2)

khat_x = kx/knorm
khat_y = ky/knorm
khat_z = kz/knorm

# forward scattering unit vector
ko_x, ko_y, ko_z = beam.get_unit_s0()

# Now we can determine which pixel in the image array
# corresponds to qx,qy,qz by solving  the system
#   det_qx = qx
#   det_qy = qy
#   det_qz = qz
per_wavelen = 1/beam.get_wavelength()
det_qx = per_wavelen * (khat_x - ko_x)
det_qy = per_wavelen * (khat_y - ko_y)
det_qz = per_wavelen * (khat_z - ko_z)


# The Bragg peak "peak" is where all pixel q coordinates are closest to the Bragg peak q
close_dist = np.sqrt( (det_qx-qx)**2 + (det_qy - qy)**2 + (det_qz-qz)**2 )
min_index = np.argmin(close_dist)
# get the centroid of the spot in pixel units
slow_peak, fast_peak = np.unravel_index(min_index, (slow_dim*oversamp, fast_dim*oversamp))
fast_peak /= oversamp
slow_peak /= oversamp


# simulate the spot
SIM = simtbx.nanoBragg.nanoBragg(whole_det, beam)
if args.updatecenter:
    SIM.beam_center_mm = whole_det[0].get_beam_centre(beam.get_s0())
SIM.Amatrix = sqr(cryst.get_A()).transpose().elems
SIM.xtal_shape = shapetype.Gauss
SIM.default_F = 1e3
SIM.F000 = 0
SIM.oversample=2
SIM.Ncells_abc = (5,5,5)
SIM.interpolate=0
SIM.show_params()
SIM.add_nanoBragg_spots()
SIM.raw_pixels *= 2000
wholepix = SIM.raw_pixels.as_numpy_array()

i1,i2,j1,j2 = map(int, [fast_peak-10, fast_peak+10, slow_peak-10, slow_peak+10])
subimg = wholepix[j1:j2, i1:i2]

Ivals =subimg.ravel()
Isum = Ivals.sum()
Y,X = np.indices(subimg.shape)
X = X + i1 + 0.5
Y = Y + j1 + 0.5
xcom = (X.ravel()*Ivals).sum() / Isum
ycom = (Y.ravel()*Ivals).sum() / Isum

x_dev = xcom-fast_peak
y_dev = ycom-slow_peak

pix_offset = np.sqrt(x_dev**2 + y_dev**2)

print("True position on the camera is at fast=%.4f slow=%.4f" % (fast_peak, slow_peak))
print("Centroid spot position on the camera is at fast=%.4f slow=%.4f" % (xcom, ycom))
print("offset = %.4f pixels" % pix_offset)

if args.plot:
    import pylab as plt
    plt.figure(1)
    plt.imshow(wholepix, extent=(0,fast_dim, slow_dim, 0))
    plt.plot([fast_peak], [slow_peak], 'ro', ms=2, label="true center")
    plt.plot([xcom], [ycom], 'ks', ms=2, label='center of mass')
    plt.legend()
    plt.show()

assert pix_offset < 0.3
print("OK!")
