from __future__ import division
from scipy import special
import numpy as np
from simtbx.nanoBragg.tst_gaussian_mosaicity2 import check_distributions
from scitbx.matrix import sqr, col
import pylab as plt


d_eta_tensor_a = sqr((1, 0, 0,
                      0, 0, 0,
                      0, 0, 0))
d_eta_tensor_b = sqr((0, 0, 0,
                      0, 1, 0,
                      0, 0, 0))
d_eta_tensor_c = sqr((0, 0, 0,
                      0, 0, 0,
                      0, 0, 1))
d_eta_tensor_d = sqr((0, 1, 0,
                      1, 0, 0,
                      0, 0, 0))
d_eta_tensor_e = sqr((0, 0, 0,
                      0, 0, 1,
                      0, 1, 0))
d_eta_tensor_f = sqr((0, 0, 1,
                      0, 0, 0,
                      1, 0, 0))

d_eta_tensor_isotropic = sqr((1,0,0,
                              0,1,0,
                              0,0,1))


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


def _compute(rot_ax, ang_idx, eta_eff, Cvec, derivs=None):

    # store positive and negative rotation matrices in a list
    Us, Uprimes, Udblprimes = [],[],[]

    # rotation amount
    factor = np.sqrt(2) * special.erfinv(ang_idx) * np.pi / 180.
    rot_ang = eta_eff * factor

    # first derivatives
    if derivs is not None:
        d_theta_d_etas = []
        for d_eta_tensor in derivs:
            d_theta_d_etas.append(Cvec.dot(d_eta_tensor*Cvec) * factor)
        # second deriv of theta w.r.t eta is 0!

    # do for both postivie and negative rotations for even distribution of Umats
    for rot_sign in [1, -1]:
        U = rot_ax.axis_and_angle_as_r3_rotation_matrix(rot_sign*rot_ang, deg=False)
        Us.append(U)

        if derivs is not None:
            dU_d_theta = rot_ax.axis_and_angle_as_r3_derivative_wrt_angle(rot_sign*rot_ang, deg=False) # 1st deriv
            d2U_d_theta2 = rot_ax.axis_and_angle_as_r3_derivative_wrt_angle(rot_sign*rot_ang, deg=False, second_order=True)  # second deriv
            for d_theta_d_eta in d_theta_d_etas:
                dU_d_eta = rot_sign*dU_d_theta*d_theta_d_eta
                d2U_d_eta2 = d2U_d_theta2*(d_theta_d_eta**2)
                Uprimes.append(dU_d_eta)
                Udblprimes.append(d2U_d_eta2)

    return Us, Uprimes, Udblprimes


def generate_Umats(eta_tensor, crystal=None, plot=None,
                   how=1, num_random_samples=500, compute_derivs=True, verbose=False):
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
    if plot is not None:
        from mpl_toolkits.mplot3d import Axes3D # noqa
    if how == 1:
        for i in [1, 2, 3, 5, 6, 7]:
            assert eta_tensor[i] == 0
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

    if num_random_samples %2 == 1:
        raise ValueError("Num random samples should be an even number")
    Nrand = int(num_random_samples/2)
    hemisph_samples = search_directions(Nrand)
    np.random.seed(8675309)
    angle_indices = np.random.permutation(np.arange(Nrand)/ Nrand)

    if plot is not None:
        x,y,z = hemisph_samples.T
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
    for i, pt in enumerate(hemisph_samples):
        rot_ax = col(pt)

        if crystal is not None:
            Ca = rot_ax.dot(unit_a)
            Cb = rot_ax.dot(unit_b)
            Cc = rot_ax.dot(unit_c)
            C = col((Ca, Cb, Cc)).normalize()
        else:
            C = col((1, 0, 0))  # arbitrary

        # effective mosaic rotation dependent on eta tensor
        eta_eff = C.dot(eta_tensor*C)

        ang_idx = angle_indices[i]
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


if __name__ == "__main__":
    cryst_dict = dict([('__id__', 'crystal'), ('real_space_a', (-48.93914505851325, -61.4985726090971, 0.23980318971727585)), ('real_space_b', (-27.63556200961052, 72.26768337463876, 13.81410546001183)), ('real_space_c', (-42.92524538136074, 33.14788397044063, -259.2845460893375)), ('space_group_hall_symbol', '-P 6 2'), ('ML_half_mosaicity_deg', 0.02676231907923616), ('ML_domain_size_ang', 4646.073492432425)])
    from dxtbx.model import Crystal
    cryst = Crystal.from_dict(cryst_dict)

    # mosaicities in degrees
    a, b, c = 0.025, 0.025, 0.075
    d, e, f = 0.01, 0.05, 0.09

    # spherical cap model
    etas = (a, 0, 0,
            0, a, 0,
            0, 0, a)

    Nmos = 1000
    # Generate mosaic models with randomized sampling, no derivatives
    U, Uprime, Udblprime =generate_Umats(etas, crystal=cryst, num_random_samples=Nmos, plot="random angle sampling along spiral", how=2, verbose=True, compute_derivs=True)
    from simtbx.nanoBragg.tst_gaussian_mosaicity import plotter2
    etas_2 = (a/2.,0,0, 0,a/2.,0,0,0,a/2.)
    U_by2, _,_ =generate_Umats(etas_2, crystal=cryst, num_random_samples=Nmos, plot="random angle sampling along spiral", how=2, verbose=True, compute_derivs=False)
    plotter2(U, U_by2, False)

    # isotropic case should not depend on the crystal model
    Unoxtal, _, _ =generate_Umats(etas, crystal=None, num_random_samples=Nmos, plot=None, how=2, compute_derivs=False)
    for i in range(len(U)):
        assert np.allclose(U[i].elems, Unoxtal[i].elems)

    eps = 0.00001  # degrees
    etas_shifted = (a+eps, 0, 0,
                  0, a+eps, 0,
                  0, 0, a+eps)
    etas_shifted_minus = (a-eps, 0, 0,
                          0, a-eps, 0,
                          0, 0, a-eps)
    Ushift, _,_ = generate_Umats(etas_shifted, num_random_samples=Nmos, plot=None, how=2)
    Ushift_minus, _,_ =generate_Umats(etas_shifted_minus, num_random_samples=Nmos, plot=None, how=2)
    #from simtbx.nanoBragg.tst_gaussian_mosaicity2 import check_finite_second_order, check_finite
    from simtbx.nanoBragg.tst_gaussian_mosaicity2 import check_finite
    check_finite(U, Ushift, Uprime, eps)
    print("Finite differences check!")
    #check_finite_second_order(U, Ushift, Ushift_minus, Udblprime, eps)
    #print("Finite second differences check!")

    # anisotropic cap model
    etas = (a, 0, 0,
            0, b, 0,
            0, 0, c)
    generate_Umats(etas, cryst, num_random_samples=Nmos,
                   plot="3-parameter anisotropic random sampling along spiral", how=1, compute_derivs=False)

    # fully anisotropic cap model
    etas = (a, d, f,
            d, b, e,
            f, e, c)
    generate_Umats(etas, cryst, num_random_samples=Nmos, plot="6-parameter anisotropic sampling along spiral", how=0)
