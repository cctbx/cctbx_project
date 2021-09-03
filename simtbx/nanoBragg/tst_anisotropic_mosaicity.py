from __future__ import division
import random
import numpy as np
from simtbx.nanoBragg.tst_gaussian_mosaicity2 import check_distributions
from scitbx.matrix import sqr, col
import math
import pylab as plt

from simtbx.nanoBragg.anisotropic_mosaicity import (
  search_directions, _compute,
  d_eta_tensor_a, d_eta_tensor_b, d_eta_tensor_c, d_eta_tensor_isotropic )

class AnisoUmats:
    def __init__(self, num_random_samples=500, seed=8675309):
        """
        num_random_samples, number of mosaic domain umatrices to generate (should be even)
        seed, random seed used for permuting the 1-to-1 rotation axis / rotation angle mapping
        """
        if num_random_samples %2 == 1:
            raise ValueError("Num random samples should be an even number")
        Nrand = int(num_random_samples/2)
        self.hemisph_samples = search_directions(Nrand)
        R = random.Random()
        R.seed(seed)
        self.angle_indices = [(float(i)+0.5) / Nrand for i in range(Nrand)]
        R.shuffle(self.angle_indices) # evenly spaced numbers between 0 and 1, in random order

    def generate_isotropic_Umats(self, eta, compute_derivs=True):
        assert compute_derivs==False, "True not implemented"

        if eta < 0:
                raise ValueError("Mosaicities need to be >= 0")

        derivs = [d_eta_tensor_isotropic]

        all_U = []
        all_Uprime = []
        all_Udblprime = []

        default_C = col((1,0,0))
        for i, pt in enumerate(self.hemisph_samples):
            rot_ax = col(pt)

            ang_idx = self.angle_indices[i]
            U, Up, Udp = _compute(rot_ax, ang_idx, eta, default_C,
                                  derivs=derivs if compute_derivs else None)
            all_U += U
            if compute_derivs:
                all_Uprime += Up
                all_Udblprime += Udp

        if compute_derivs:
            assert len(all_U) == len(all_Uprime) == len(all_Udblprime)

        return all_U, all_Uprime, all_Udblprime

    def generate_Umats(self, eta_tensor, crystal=None, plot=None,
                       how=1, compute_derivs=True, verbose=False):
        """
        :param crystal: dxtbx crystal model
        :param eta_tensor: sequence of 9 numbers specifying the mosacitiy tensor
        :param num_axes: how many points to sample the unit hemisphere with
        :param num_angles_per_axis:  produces 2x this number of angles per axes to sample the angle distribution
        :param plot: string, if not None then produce a plot with this value as the title
        :param how: if 0, then do the full treatment (6 Umat derivatives)
                    if 1, then do the diagonal only (3 Umat derivatives)
                    if 2, then there is no anisotropy (1 Umat derivative)
        :param num_random_samples, if an integer, ignore num_axes and num_angles_per_axis and
            use this number to generate random samples
        :param compute_derivs, boolean, if False, only compute the Umats
        :return:
        """
        assert compute_derivs==False, "True not implemented"
        assert how==1, "other methods not implemented"
        if plot is not None:
            from mpl_toolkits.mplot3d import Axes3D # noqa
        if how == 1:
            for i in [1, 2, 3, 5, 6, 7]:
                assert eta_tensor[i] == 0
            for idiag in [0,4,8]:
              eta_tensor[idiag] = 1./(eta_tensor[idiag] * eta_tensor[idiag])
              #units are inverse degrees-squared
        elif how == 2:
            assert eta_tensor[0] == eta_tensor[4] == eta_tensor[8]
        else:
            assert eta_tensor[1] == eta_tensor[3]
            assert eta_tensor[2] == eta_tensor[6]
            assert eta_tensor[5] == eta_tensor[7]

        if how in [0, 1]:
            assert crystal is not None

        for val in eta_tensor:
            if val < 0:
                raise ValueError("Mosaicities need to be >= 0")

        if crystal is not None:
            a, b, c = map(col, crystal.get_real_space_vectors())
            unit_a = a.normalize()
            unit_b = b.normalize()
            unit_c = c.normalize()

        if plot is not None:
            f = plt.figure()
            ax = f.add_subplot(111, projection='3d')
            f2  = plt.figure()
            ax2 = f2.add_subplot(111, projection='3d')
            ax.set_title(plot)
            ax2.set_title(plot)

        if plot is not None:
            x,y,z = self.hemisph_samples.T
            ax2.scatter(x,y,z,s=5,marker='s', color='r', alpha=0.5)

        eta_tensor = sqr(eta_tensor)

        if how == 0:  # full treatment
            derivs = [d_eta_tensor_a, d_eta_tensor_b, d_eta_tensor_c,
                      d_eta_tensor_d, d_eta_tensor_e, d_eta_tensor_f]
        elif how == 1:  # diag only treatment
            derivs = [d_eta_tensor_a, d_eta_tensor_b, d_eta_tensor_c]
        else:  # isotropic treatment
            derivs = [d_eta_tensor_isotropic]

        all_U = []
        all_Uprime = []
        all_Udblprime = []
        for i, pt in enumerate(self.hemisph_samples):
            rot_ax = col(pt)

            if crystal is not None:
                Ca = rot_ax.dot(unit_a)
                Cb = rot_ax.dot(unit_b)
                Cc = rot_ax.dot(unit_c)
                C = col((Ca, Cb, Cc)).normalize() # XXX Fix Me this is the rot ax in lab space
                                                  # but it needs to go into skew space to apply to eta
            else:
                C = col((1, 0, 0))  # arbitrary

            # effective mosaic rotation dependent on eta tensor
            eta_eff = 1./math.sqrt(C.dot(eta_tensor*C))

            ang_idx = self.angle_indices[i]
            U, Up, Udp = _compute(rot_ax, ang_idx, eta_eff, C,
                                  derivs=derivs if compute_derivs else None)
            all_U += U
            if compute_derivs:
                all_Uprime += Up
                all_Udblprime += Udp

        if plot is not None:
            A = col((1,0,0))
            Anew = []
            for umat in all_U:
                Anew.append( umat * A)
            x,y,z = np.array(Anew).T
            ax.scatter(x,y,z,s=2, alpha=1)
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            ax.set_zlabel("z")
            plt.show()

        if compute_derivs and how == 0:
            assert 6*len(all_U) == len(all_Uprime) == len(all_Udblprime)
        elif compute_derivs and how == 1:
            assert 3 * len(all_U) == len(all_Uprime) == len(all_Udblprime)
        elif compute_derivs and how == 2:
            assert len(all_U) == len(all_Uprime) == len(all_Udblprime)

        if verbose:
            nm_angles = check_distributions.get_angular_rotation(all_U)
            nm_rms_angle = np.sqrt(np.mean(nm_angles * nm_angles))
            print("Normal rms angle is ", nm_rms_angle)

        return all_U, all_Uprime, all_Udblprime

