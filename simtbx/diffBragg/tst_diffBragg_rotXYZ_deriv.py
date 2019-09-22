
from simtbx.diffBragg.tst_diffBragg_rotXYZ import args

from simtbx.diffBragg.tst_diffBragg_rotXYZ import  get_diffBragg_instance

def main(rot_idx):
    import numpy as np
    from scitbx.matrix import col, sqr

    vals = []
    theta_vals = [0.00001 * (2**i) for i in range(1, 25, 1)]
    for theta in theta_vals:
        theta = theta * np.pi / 180

        D = get_diffBragg_instance()
        D.vectorize_umats()

        # STEP 1: simulate the un-perturbed image:
        D.add_diffBragg_spots()
        img0 = D.raw_pixels_roi.as_numpy_array()

        # STEP 2: simulate the same crystal , directly perturbed by the rotation matrix in python

        if rot_idx == 0:
            rot_ax = col((-1, 0, 0))
        elif rot_idx == 1:
            rot_ax = col((0, -1, 0))
        elif rot_idx == 2:
            rot_ax = col((0, 0, -1))
        else:
            assert False, "Rot idx should be 0,1 or 2"

        R = rot_ax.axis_and_angle_as_r3_rotation_matrix(theta, deg=False)

        # apply the rotation matrix to the crystal Amatrix:
        Umat_orig = D.Umatrix

        #Arecip_orig = sqr(D.Amatrix)
        #Areal = Arecip_orig.inverse()
        #Areal = R * Areal
        #Arecip = Areal.inverse()
        ## put back in diffBragg:
        #D.Amatrix = Arecip

        # simulate the scattering in the rotated crystal:
        D.raw_pixels_roi *= 0
        D.Umatrix = R
        D.add_diffBragg_spots()
        img = D.raw_pixels_roi.as_numpy_array()

        # STEP3 : compute finite differenceL
        finite_diff = (img-img0)/theta

        # STEP4 : compute the analytical derivative evaluated at theta=0
        D.refine(rot_idx)
        D.initialize_managers()

        D.Umatrix = Umat_orig
        #D.zero_raw_pixel_rois()
        D.raw_pixels_roi *= 0
        D.set_value(rot_idx, 0)
        D.add_diffBragg_spots()

        ana_deriv = D.get_derivative_pixels(rot_idx).as_numpy_array()

        if args.plot:
            import pylab as plt
            plt.subplot(121)
            plt.imshow(finite_diff)
            plt.title("finite diff.")
            plt.subplot(122)
            plt.imshow(ana_deriv)
            plt.title("analytical")
            #plt.gca().images[0].set_clim(plt.gcf().axes[0].images[0].get_clim())
            plt.suptitle("Theta = %f deg. " % (theta * 180 / np.pi))
            plt.draw()
            plt.pause(0.2)

        val = abs(ana_deriv-finite_diff)
        print ("Theta = %.4f deg, Average discrepancy = %.3g" %
                (theta*180 / np.pi, val.mean()))

        #print  abs(ana_deriv.max()-finite_diff.max()) / (.5*ana_deriv.max() + .5*finite_diff.max())
        vals.append(val)

    # STEP6: for the smallest perturbation, assert finite difference
    # is equivalent to analytical derivative within 1e-6 units ( per pixel)
    assert np.all(vals[0] < 1e-3)

    from scipy import stats
    x = theta_vals[:16]
    y = [v.mean() for v in vals][:16]
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    assert slope > 0
    assert r_value > 0.9
    assert p_value < 1e-5

    # TODO: use second derivative and finite difference error model to check error scaling

    if args.plot:
        plt.close()
        plt.plot(theta_vals, [v.mean() for v in vals], '.')
        ax = plt.gca()
        ax.set_xlabel("theta (degrees)")
        ax.set_xscale("log")
        ax.set_ylabel(r"$\langle |\,$finite_diff - analytical$\,| \rangle$", fontsize=14)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        plt.subplots_adjust(left=0.12)
        plt.show()


if __name__ == "__main__":
    main(rot_idx=0)
    main(rot_idx=1)
    main(rot_idx=2)
    print("OK!")
