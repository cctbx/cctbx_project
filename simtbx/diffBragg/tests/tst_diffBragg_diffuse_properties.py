"""
This test checks the setter and getter for Ncells parameter
"""
from __future__ import division

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--cuda", action="store_true")
parser.add_argument("--plot", action='store_true')
parser.add_argument("--onlyDiffuse", action='store_true')
parser.add_argument("--sigma", default=0.3, type=float)
parser.add_argument("--gamma", default=1, type=float)

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

det_sh = 1024,1024
print("OK")
ONLY_DIFFUSE = args.onlyDiffuse
if args.plot:
    all_imgs = []
    for iso in [True, False]:
        F, axs = plt.subplots(nrows=2, ncols=3)
        F.suptitle("isotropic gamma=%s" % iso)
        axs = [a for ax in axs for a in ax]
        plots = [(False, None, None, False),
                 (True, sigma, gamma, ONLY_DIFFUSE),
                 (True, sigma, 10*gamma, ONLY_DIFFUSE),
                 (True, sigma, 100*gamma, ONLY_DIFFUSE),
                 (True, 2*sigma, gamma, ONLY_DIFFUSE),
                 (True, 4*sigma, gamma, ONLY_DIFFUSE)]
        for i_ax, (use_diffuse, sig, gam, only_diffuse) in enumerate(plots):
            S.D.use_diffuse = use_diffuse
            S.D.raw_pixels_roi*=0
            if use_diffuse:
                if iso:
                    ga1=ga2=ga3=gam
                else:
                    ga1,ga2,ga3 = 6*gam, 3*gam, gam

                S.D.diffuse_gamma = ga1, ga2, ga3
                S.D.diffuse_sigma = sig, sig, sig
            S.D.only_diffuse = False #only_diffuse
            S.D.verbose = 2
            S.D.add_diffBragg_spots()
            img = S.D.raw_pixels_roi.as_numpy_array().reshape(det_sh)
            all_imgs.append(img)
            m = img.mean()
            s = img.std()
            vmin = m-s
            vmax = m+s
            axs[i_ax].imshow(img, vmin=vmin, vmax=vmax)

            if use_diffuse:
                axs[i_ax].set_title("sig=%.2f, gam=%.1f %.1f %.1f" % (sig,ga1, ga2, ga3), fontsize=8)
            else:
                axs[i_ax].set_title("no diffuse", fontsize=8)

    show()
