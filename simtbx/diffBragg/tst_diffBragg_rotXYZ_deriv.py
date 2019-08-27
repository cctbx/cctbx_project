
from simtbx.diffBragg.tst_diffBragg_rotXYZ import args

from simtbx.diffBragg.tst_diffBragg_rotXYZ import  get_diffBragg_instance

def main():
    import numpy as np

    for thetaX in np.linspace(0.00001,20, 50):
        thetaX = thetaX * np.pi / 180

        D = get_diffBragg_instance()
        D.vectorize_umats()

        # STEP 1: simulate the un-perturbed image:
        D.add_diffBragg_spots()
        imgX0 = D.raw_pixels.as_numpy_array()

        # STEP 2: simulate the same crystal , directly perturbed by the rotation matrix in python
        D.raw_pixels*=0
        from scitbx.matrix import col, sqr
        x = col((-1, 0, 0))  # rotate about X by amount thetaX
        RX = x.axis_and_angle_as_r3_rotation_matrix(thetaX, deg=False)
        # apply the rotation matrix to the crystal Amatrix:
        Arecip_orig = sqr(D.Amatrix)
        Areal = Arecip_orig.inverse()
        Areal = RX * Areal
        Arecip = Areal.inverse()
        # put back in diffBragg:
        D.Amatrix = Arecip
        # simulate the scattering in the rotated crystal:
        D.add_diffBragg_spots()
        imgX = D.raw_pixels.as_numpy_array()

        # STEP3 : compute finite differenceL
        finite_diff = (imgX-imgX0)/thetaX

        # STEP4 : compute the analytical derivative evaluated at thetaX=0
        rotX = 0
        D.refine(rotX)
        D.initialize_managers()

        D.raw_pixels *= 0
        D.set_value(rotX, 0)
        D.Amatrix = Arecip_orig
        D.add_diffBragg_spots()

        ana_deriv = D.get_derivative_pixels(rotX).as_numpy_array()

        if args.plot:
            import pylab as plt
            plt.subplot(121)
            plt.imshow(finite_diff)
            plt.subplot(122)
            plt.imshow(ana_deriv)
            #plt.gca().images[0].set_clim(plt.gcf().axes[0].images[0].get_clim())
            plt.draw()
            plt.pause(.7)

        #print ( "Theta = %.4f deg, SUM discrepancy = %.3g" %
        #        (thetaX*180 / np.pi,  abs(ana_deriv.max()-finite_diff.max())))  #np.abs(finite_diff - ana_deriv).sum())  )

        print  abs(ana_deriv.max()-finite_diff.max()) / (.5*ana_deriv.max() + .5*finite_diff.max())
        #from IPython import embed
        #embed()


if __name__ == "__main__":
    main()
    print("OK")
