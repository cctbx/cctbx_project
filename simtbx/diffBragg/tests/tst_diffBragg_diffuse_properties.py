"""
This test checks the setter and getter for Ncells parameter
"""
from __future__ import division

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--cuda", action="store_true")
parser.add_argument("--plot", action='store_true')
parser.add_argument("--sigma", default=3.16, type=float)
parser.add_argument("--gamma", default=50, type=float)

args = parser.parse_args()
if args.cuda:
    import os
    os.environ["DIFFBRAGG_USE_CUDA"] = "1"

from simtbx.nanoBragg import sim_data
from pylab import *


S = sim_data.SimData(use_default_crystal=True)
S.detector = sim_data.SimData.simple_detector(180, 0.1, (1024, 1024))
S.instantiate_diffBragg(verbose=0, oversample=0, auto_set_spotscale=True)
S.D.spot_scale = 100000
S.D.use_diffuse = True
gamma = args.gamma
sigma = args.sigma
S.D.diffuse_gamma = gamma, gamma, gamma
S.D.diffuse_sigma = sigma, sigma, sigma

assert np.allclose(S.D.diffuse_gamma, (gamma, gamma, gamma))
assert np.allclose(S.D.diffuse_sigma, (sigma, sigma, sigma))

S.D.initialize_managers()
S.D.add_diffBragg_spots()
det_sh = 1024,1024
img = S.D.raw_pixels_roi.as_numpy_array().reshape(det_sh)

print("OK")
if args.plot:
    S.D.use_diffuse = False
    S.D.raw_pixels_roi*=0
    S.D.add_diffBragg_spots()
    img0 = S.D.raw_pixels_roi.as_numpy_array().reshape(det_sh)

    S.D.use_diffuse = True
    S.D.diffuse_gamma = (2*gamma, 2*gamma, 2*gamma)
    S.D.raw_pixels_roi*=0
    S.D.add_diffBragg_spots()
    img2 = S.D.raw_pixels_roi.as_numpy_array().reshape(det_sh)

    S.D.use_diffuse = True
    S.D.diffuse_gamma = (4*gamma, 4*gamma, 4*gamma)
    S.D.raw_pixels_roi*=0
    S.D.add_diffBragg_spots()
    img4 = S.D.raw_pixels_roi.as_numpy_array().reshape(det_sh)

    S.D.diffuse_gamma = (gamma, gamma, gamma)
    S.D.diffuse_sigma = (2*sigma, 2*sigma,2*sigma)
    S.D.raw_pixels_roi*=0
    S.D.add_diffBragg_spots()
    img2_sigma = S.D.raw_pixels_roi.as_numpy_array().reshape(det_sh)

    S.D.diffuse_sigma = (4*sigma, 4*sigma,4*sigma)
    S.D.raw_pixels_roi*=0
    S.D.add_diffBragg_spots()
    img4_sigma = S.D.raw_pixels_roi.as_numpy_array().reshape(det_sh)




    FS=8
    subplot(231)
    title("no diffuse", fontsize=FS)
    imshow(img0)
    gca().set_xticks([])
    gca().set_yticks([])

    subplot(232)
    title("sigma=%.1f, gamma=%.1f" % (sigma, gamma), fontsize=FS)
    imshow(img)
    gca().set_xticks([])
    gca().set_yticks([])

    subplot(233)
    title("sigma=%.1f, gamma=%.1f" % (sigma, 2*gamma), fontsize=FS)
    imshow(img2)
    gca().set_xticks([])
    gca().set_yticks([])

    subplot(234)
    title("sigma=%.1f, gamma=%.1f" % (sigma, 4*gamma), fontsize=FS)
    imshow(img4)
    gca().set_xticks([])
    gca().set_yticks([])

    subplot(235)
    title("sigma=%.1f, gamma=%.1f" % (2*sigma, gamma), fontsize=FS)
    imshow(img2_sigma)
    gca().set_xticks([])
    gca().set_yticks([])

    subplot(236)
    title("sigma=%.1f, gamma=%.1f" % (4*sigma, gamma), fontsize=FS)
    imshow(img4_sigma)
    gca().set_xticks([])
    gca().set_yticks([])

    show()
