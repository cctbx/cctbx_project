from __future__ import division
from simtbx.nanoBragg import sim_data

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--plotimages", action='store_true')
parser.add_argument("--panel", choices=["x", "y", "z"], default="x", help="which origin coordinate to check")
parser.add_argument("--plotlines", action='store_true')
args = parser.parse_args()
SIM = sim_data.SimData(use_default_crystal=True)

det = sim_data.SimData.simple_detector(detector_distance_mm=150, pixelsize_mm=0.1, image_shape=(512, 512))

if args.panel == "x":
  refine_idx = 15  # id of origin coordinate in diffBragg
elif args.panel == "y":
  refine_idx = 16  # id of origin coordinate in diffBragg
else:  # args.panel==z
  refine_idx = 10

B = SIM.beam.nanoBragg_constructor_beam

# set the detector
SIM.detector = det
SIM.instantiate_diffBragg(auto_set_spotscale=True)
D = SIM.D
D.oversample_omega = True  #False
D.nopolar = True

D.refine(refine_idx)

D.initialize_managers()
D.add_diffBragg_spots()

# get the simulated image and derivative
img = D.raw_pixels_roi.as_numpy_array()
deriv = D.get_derivative_pixels(refine_idx).as_numpy_array()
D.raw_pixels_roi *= 0
D.raw_pixels *= 0

node = det[0]
node_d = node.to_dict()
O = node_d["origin"][0], node_d["origin"][1], node_d["origin"][2]
OX,OY,OZ = O

# update the detector distance
import numpy as np
all_shifts = []
all_errors = []
shifts_mm = [2*i*(1e-7) for i in range(1,30,2)]
print(shifts_mm)
import pylab as plt

from dxtbx.model import Panel

for i_shift, delta_shift in enumerate(shifts_mm):
  # update the detector model

  if args.panel=="x":
    shifted = OX+delta_shift, OY, OZ
  elif args.panel=="y":
    shifted = OX, OY+delta_shift, OZ
  else:
    delta_shift = delta_shift*10  # larger shift for Z
    shifted = OX, OY, OZ + delta_shift

  node_d["origin"] = shifted
  det[0] = Panel.from_dict(node_d)

  D.update_dxtbx_geoms(det, B, 0)
  D.add_diffBragg_spots()
  img_forward = D.raw_pixels_roi.as_numpy_array()
  D.raw_pixels_roi *= 0
  D.raw_pixels *= 0
  delta_shift_meters = delta_shift*1e-3
  fdiff = (img_forward-img)/delta_shift_meters

  bragg = img > 1e-2

  error = (np.abs(fdiff[bragg] - deriv[bragg])).mean()
  all_errors.append(error)

  all_shifts.append(delta_shift_meters)

  print("Error=%2.7g shift=%2.7g um" % (error, delta_shift*1000))
  if args.plotimages:
    plt.clf()
    y = slice(40, 65, 1)
    x = slice(415, 437, 1)
    plt.subplot(121)
    plt.imshow(fdiff[y, x])
    plt.title("finite diff")
    plt.subplot(122)
    plt.imshow(deriv[y, x])
    plt.title("analytical")
    plt.draw()
    plt.suptitle("delta=%f mm, Shift %d / %d"
                 % ( delta_shift, i_shift + 1, len(shifts_mm)))
    plt.pause(0.3)

if args.plotlines:
  plt.close()

  plt.plot(all_shifts, all_errors, 'o')
  plt.suptitle("finite 1st difference error", fontsize=16)
  plt.xlabel("h", fontsize=16)
  plt.show()

from scipy.stats import linregress
l = linregress(all_shifts, all_errors)
print("finite diff l.rvalue=%10.7g" % l.rvalue)
assert l.rvalue > .99
assert l.slope > 0
assert l.pvalue < 1e-6

print("OK!")
