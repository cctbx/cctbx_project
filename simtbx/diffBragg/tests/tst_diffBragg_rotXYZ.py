
from argparse import ArgumentParser
parser = ArgumentParser("diffBragg tests")
parser.add_argument("--plot", action='store_true')
args = parser.parse_args()


def main():
    import numpy as np
    from simtbx.diffBragg.utils import get_diffBragg_instance
    from scitbx.matrix import col, sqr

    n_trials = 10
    np.random.seed(n_trials)
    angles_XYZ = np.random.random((n_trials, 3)) * 3 * np.pi / 180.
    print (angles_XYZ*180 / np.pi)

    D = get_diffBragg_instance()

    rotX, rotY, rotZ = 0, 1, 2
    D.refine(rotX)  # rotX
    D.refine(rotY)  # rotY
    D.refine(rotZ)  # rotZ

    D.initialize_managers()
    D.vectorize_umats()

    x = col((-1, 0, 0))
    y = col((0, -1, 0))
    z = col((0, 0, -1))
    Uorig = sqr(D.Umatrix)

    if args.plot:
        import pylab as plt
        plt.figure()
        axA = plt.subplot(121)
        axB = plt.subplot(122)
    for i_ang, (thetaX, thetaY, thetaZ) in enumerate(angles_XYZ):

        RX = x.axis_and_angle_as_r3_rotation_matrix(thetaX, deg=False)
        RY = y.axis_and_angle_as_r3_rotation_matrix(thetaY, deg=False)
        RZ = z.axis_and_angle_as_r3_rotation_matrix(thetaZ, deg=False)

        Misset = RX*RY*RZ
        U = Misset*Uorig
        D.raw_pixels *= 0
        D.set_value(rotX,0)
        D.set_value(rotY,0)
        D.set_value(rotZ,0)
        D.Umatrix = U.elems
        D.add_diffBragg_spots()
        imgA = D.raw_pixels.as_numpy_array()

        D.raw_pixels *= 0
        D.set_value(rotX, thetaX)
        D.set_value(rotY, thetaY)
        D.set_value(rotZ, thetaZ)
        D.Umatrix = Uorig
        D.add_diffBragg_spots()
        imgB = D.raw_pixels.as_numpy_array()
        if args.plot:
            m = imgA[ imgA > 0].mean()
            s = imgA[ imgA > 0].std()

            axA.clear()
            axB.clear()
            axA.imshow( imgA, vmin=0, vmax=m+2*s, cmap='gnuplot')
            axB.imshow( imgB, vmin=0, vmax=m+2*s, cmap='gnuplot')
            plt.draw()
            plt.pause(1.)

        assert(np.allclose(imgA, imgB, atol=1e-4))
        print("OK (%d / %d)" % (i_ang+1, len(angles_XYZ)))


if __name__ == "__main__":
    main()
    print("OK")
