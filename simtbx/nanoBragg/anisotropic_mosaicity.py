from __future__ import division
import random
from scipy import special

try:
    from collections.abc import Iterable
except ModuleNotFoundError:
    from collections import Iterable

import numpy as np
from simtbx.nanoBragg.tst_gaussian_mosaicity2 import check_distributions
from scitbx.matrix import sqr, col

def search_directions(N=1000):
    """
    See Journal of Magnetic Resonance 138, 288–297 (1999)
    equation A6
    :param N: number of points on hemisphere
    :return: Nx3 numpy array of unit vectors
    """
    ti = (np.arange(1, N+1) - 0.5) / N
    THETA = np.arccos(ti)
    PHI = np.sqrt(np.pi*N) * np.arcsin(ti)

    u_vecs = np.zeros((N, 3))
    x = np.sin(THETA) * np.cos(PHI)
    y = np.sin(THETA) * np.sin(PHI)
    z = np.cos(THETA)
    u_vecs[:N, 0] = x
    u_vecs[:N, 1] = y
    u_vecs[:N, 2] = z

    return u_vecs


def _compute(rot_ax, ang_idx, eta_eff, Cvec, derivs=None, second_derivs=None):
    """

    :param rot_ax: scitbx.matrix.col rotation axis
    :param ang_idx: list of indices corresponding to uniform samplings of the cummulative distribution function
    :param eta_eff: float, effective mosaicity in degrees
    :param Cvec: scitbx.matrix.col the vector of rotation axis along the projects
    :param derivs: list of d_etaEffective_d_eta derivative tensors , should be len 1 for isotropic models, else len 3
    :param second_derivs: same as derivs, yet dsquared_detaEffecitve_d_eta_squared
    :return: Umatrices, and there first and second derivatives w.r.t. eta
    """
    # store positive and negative rotation matrices in a list
    Us, Uprimes, Udblprimes = [],[],[]

    # rotation amount
    factor = np.sqrt(2) * special.erfinv(ang_idx) * np.pi / 180.
    rot_ang = eta_eff * factor

    # first derivatives
    if derivs is not None:
        d_theta_d_eta = []
        dsquared_theta_d_eta_squared = []
        common_term = -(0.5 *eta_eff**3) * factor
        for d, d2 in zip(derivs, second_derivs):
            G = np.dot(Cvec, np.dot(d, Cvec))
            d_theta_d_eta.append(common_term*G)

            G2 = np.dot(Cvec, np.dot(d2, Cvec))
            dsquared_theta_d_eta_squared.append(common_term*(-1.5* eta_eff**2 * G**2 + G2))

    # do for both postivie and negative rotations for even distribution of Umats
    for rot_sign in [1, -1]:
        U = rot_ax.axis_and_angle_as_r3_rotation_matrix(rot_sign*rot_ang, deg=False)
        Us.append(U)

        if derivs is not None:
            dU_d_theta = rot_ax.axis_and_angle_as_r3_derivative_wrt_angle(rot_sign*rot_ang, deg=False) # 1st deriv
            d2U_d_theta2 = rot_ax.axis_and_angle_as_r3_derivative_wrt_angle(rot_sign*rot_ang, deg=False, second_order=True)  # second deriv
            for d, d2 in zip(d_theta_d_eta, dsquared_theta_d_eta_squared):
                dU_d_eta = rot_sign*dU_d_theta*d
                d2U_d_eta2 = d2U_d_theta2*(d**2) + dU_d_theta*d2
                Uprimes.append(dU_d_eta)
                Udblprimes.append(d2U_d_eta2)

    return Us, Uprimes, Udblprimes


class AnisoUmats:

    def __init__(self, num_random_samples=500, seed=8675309):
        """
        num_random_samples, number of mosaic domain umatrices to generate with generate_Umats method (should be even)
        seed, random seed used for permuting the 1-to-1 rotation axis / rotation angle mapping
        """
        if num_random_samples %2 == 1:
            raise ValueError("Num random samples should be an even number")
        Nrand = int(num_random_samples/2)
        self.hemisph_samples = search_directions(Nrand)
        R = random.Random()
        R.seed(seed)
        self.angle_indices = [float(i) / Nrand for i in range(Nrand)]
        R.shuffle(self.angle_indices)

    def generate_Umats(self, eta, crystal=None, transform_eta=False,
                       how=1, compute_derivs=True, verbose=False ):
        """
        :param eta: float or 3-tuple specfying the mosaicity in degrees
        :param crystal: dxtbx crystal model
        :param transform_eta: bool, if True, then form an eta tensor that uses the
        crystal axes as its basis (experimental, not sure if its correct)
        references:
        https://en.wikipedia.org/wiki/Ellipsoid#As_a_quadric
        https://math.stackexchange.com/a/1119690/721977

        :param num_axes: how many points to sample the unit hemisphere with
        :param num_angles_per_axis:  produces 2x this number of angles per axes to sample the angle distribution
        :param how: if 0, then do the full treatment (6 Umat derivatives)
                    if 1, then do the diagonal only (3 Umat derivatives)
                    if 2, then there is no anisotropy (1 Umat derivative)
        :param num_random_samples, if an integer, ignore num_axes and num_angles_per_axis and
            use this number to generate random samples
        :param compute_derivs, boolean, if False, only compute the Umats
        :return:
        """
        if how==0:
            raise NotImplementedError("Still working out details and use cases for 6-parameter mosaicity model.")

        if isinstance(eta, Iterable):
            assert len(eta) == 3
            eta_a, eta_b, eta_c = eta
        else:
            if not how==2:
                raise ValueError("passing in a float for eta assumes how=2 (isotropic model)")
            eta_a = eta_b = eta_c = eta

        eta_tensor = np.diag(1./np.array([eta_a, eta_b, eta_c])**2)

        if how==1 and transform_eta:
            a,b,c = map(lambda x: col(x).normalize(), crystal.get_real_space_vectors())
            # S is a matrix whose columns are a,b,c
            S = np.reshape( sqr(a.elems+b.elems+c.elems).transpose(), (3,3))
            Sinv = np.linalg.inv(S)
            # re-write eta_tensor in crystal basis ?
            eta_tensor = np.dot(S, np.dot(eta_tensor, Sinv))

        if how==1:
            d_eta_tensor_a = np.diag([-2*eta_a**-3,0,0])
            d_eta_tensor_b = np.diag([0,-2*eta_b**-3,0])
            d_eta_tensor_c = np.diag([0,0, -2*eta_c**-3])
            if transform_eta:
                d_eta_tensor_a = np.dot(S, np.dot(d_eta_tensor_a, Sinv))
                d_eta_tensor_b = np.dot(S, np.dot(d_eta_tensor_b, Sinv))
                d_eta_tensor_c = np.dot(S, np.dot(d_eta_tensor_c, Sinv))
            derivs = d_eta_tensor_a, d_eta_tensor_b, d_eta_tensor_c

            d2_eta_tensor_a = np.diag([6*eta_a**-4,0,0])
            d2_eta_tensor_b = np.diag([0,6*eta_b**-4,0])
            d2_eta_tensor_c = np.diag([0,0,6*eta_c**-4])
            if transform_eta:
                d2_eta_tensor_a = np.dot(S, np.dot(d2_eta_tensor_a, Sinv))
                d2_eta_tensor_b = np.dot(S, np.dot(d2_eta_tensor_b, Sinv))
                d2_eta_tensor_c = np.dot(S, np.dot(d2_eta_tensor_c, Sinv))
            second_derivs = d2_eta_tensor_a, d2_eta_tensor_b, d2_eta_tensor_c
        elif how==2:
            # eta_a = eta_b = eta_c
            derivs = [np.diag([-2*eta_a**-3]*3)]
            second_derivs = [np.diag([6*eta_a**-4]*3)]

        all_U = []
        all_Uprime = []
        all_Udblprime = []
        for i, pt in enumerate(self.hemisph_samples):
            rot_ax = col(pt)
            C = rot_ax

            # effective mosaic rotation dependent on eta tensor
            C_eta_C = np.dot(C, np.dot(eta_tensor, C))
            # NOTE: added in the abs to protect sqrt, prob doesnt matter since we sample +- rotation for every axis...
            eta_eff = 1/np.sqrt(np.abs(C_eta_C))

            ang_idx = self.angle_indices[i]
            U, Up, Udp = _compute(rot_ax, ang_idx, eta_eff, C,
                                  derivs=derivs if compute_derivs else None,
                                  second_derivs=second_derivs if compute_derivs else None)
            all_U += U
            if compute_derivs:
                all_Uprime += Up
                all_Udblprime += Udp

        if compute_derivs and how == 1:
            assert 3 * len(all_U) == len(all_Uprime) == len(all_Udblprime)
        elif compute_derivs and how == 2:
            assert len(all_U) == len(all_Uprime) == len(all_Udblprime)

        if verbose:
            nm_angles = check_distributions.get_angular_rotation(all_U)
            nm_rms_angle = np.sqrt(np.mean(nm_angles * nm_angles))
            print("Normal rms angle is ", nm_rms_angle)

        return all_U, all_Uprime, all_Udblprime
