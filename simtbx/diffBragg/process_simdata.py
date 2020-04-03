
# NOTE: Used for testsing only

def process_simdata(plot=False, angles=None, perturb="rotXYZ"):
    """
    returns data needed for refinement script as a 5-tuple:
    --> spot_roi, spot_hkl, perturbed Amatrix, ground truth Amatrix
        spot_abc (initial guess), simulated image
    """
    import numpy as np
    from scitbx.matrix import sqr, col
    from simtbx.diffBragg import utils
    from simtbx.diffBragg.sim_data import SimData
    if plot:
        import pylab as plt

    # STEP 1: simulate an image
    #<><><><><><><><><><><><><>
    S = SimData()
    S.instantiate_diffBragg()
    #S._add_diffBragg_spots()
    img = S.generate_simulated_image()

    # STEP 2: find the strong spots
    #<><><><><><><><><><><><><><><>
    spot_data = utils.get_spot_data(img, thresh=200)
    ss_spot, fs_spot = map(np.array, zip(*spot_data["maxIpos"]))  # slow/fast  scan coords of strong spots
    num_spots = len(ss_spot)
    if plot:
        plt.imshow(img, vmax=200)
        plt.plot( fs_spot,ss_spot,'o', mfc='none', mec='r')
        plt.title("Simulated image with strong spots marked")
        plt.show()

    # STEP 3 : perturb the known ground truth Amatrix intentionally
    #<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    Areal_GrnTru = sqr(S.D.Amatrix).inverse()
    if perturb=="rotXYZ":
        # in order to index the spots we need the Amatrix from the simulation

        # perturb the Amatrix to simulate uncertainty
        x = col((-1, 0, 0))
        y = col((0, -1, 0))
        z = col((0, 0, -1))

        if angles is None:
            np.random.seed(1)
            angles = np.random.uniform(0, 0.04, 3)
            angles = 0.3, 0.05, 0.3
        RX = x.axis_and_angle_as_r3_rotation_matrix(angles[0], deg=True)
        RY = y.axis_and_angle_as_r3_rotation_matrix(angles[1], deg=True)
        RZ = z.axis_and_angle_as_r3_rotation_matrix(angles[2], deg=True)
        misset = RX*RY*RZ
        Areal_approx = misset*Areal_GrnTru
    #elif perturb=="ucell":
    #    ucell_truth = S.Fhkl.unit_cell()
    #    from cctbx import crystal
    #    a, b, c, al, be, ga = ucell_truth.parameters()
    #    symbol = S.Fhkl.space_group_info().type().lookup_symbol()
    #    symm_perturb = crystal.symmetry("%f,%f,%f,%f,%f,%f" % (a * 0.99, b * 1.01, c, al, be, ga), symbol)
    #    ortho_matrix_approx = symm_perturb.unit_cell().orthogonalization_matrix()

    #    from IPython import embed
    #    embed()

    if plot:
        from itertools import cycle
        S.D.Amatrix = Areal_approx.inverse()
        S.D.raw_pixels*=0
        S._add_nanoBragg_spots()
        S._add_background()
        S._add_noise()
        img_approx = S.D.raw_pixels.as_numpy_array()
        imgs = cycle( (img, img_approx))
        titles = cycle(("ground truth", "perturbed"))
        count = 0
        while count < 5:
            plt.cla()
            plt.imshow(imgs.next(), vmax=200)
            title = titles.next()
            title += "\nContinuing in %d ..." % (5-count)
            plt.title(title)
            plt.draw()
            plt.pause(1)
            count += 1

    # STEP 4 : assign miller indices to strong spots using perturbed A matrix
    #<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    # get the momentum transfer vectors of each strong spot
    q_vecs = utils.x_y_to_q(fs_spot, ss_spot, detector=S.detector, beam=S.beam.xray_beams[0])

    # multiply with the real-space Amatrix in order to get fractional miller indices
    hkl = np.dot(Areal_approx.as_numpy_array(), q_vecs.T)
    hkli = list(map(lambda h: np.ceil(h - 0.5).astype(int), hkl))
    hkli = np.vstack(hkli).T

    if plot:
        plt.close()
        plt.imshow(img, vmin=0, vmax=200)
        ax = plt.gca()
        for i_hkl, (i,j) in enumerate(zip(fs_spot, ss_spot)):
            hkl_str = "%d %d %d" % tuple(hkli[i_hkl])
            ax.text(i-2,j-2, s=hkl_str)
        plt.show()

    # STEP 5:  get initial estimates of the a,b,c background plane parameters
    #<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    # make a mask of the strong spot (simple, just the smallest box fitting around the spot)
    is_bg_pixel = np.ones( img.shape, bool)
    for bb_ss, bb_fs in spot_data["bboxes"]:
        is_bg_pixel[bb_ss, bb_fs] = False

    # now fit tilting planes
    shoebox_sz = 16
    tilt_abc = np.zeros((num_spots, 3))
    spot_roi = np.zeros((num_spots, 4), int)
    if plot:
        patches = []
    img_shape = S.detector[0].get_image_size()
    for i_spot, (x_com, y_com) in enumerate(zip(fs_spot, ss_spot)):
        i1 = int(max(x_com - shoebox_sz / 2., 0))
        i2 = int(min(x_com + shoebox_sz / 2., img_shape[0]))
        j1 = int(max(y_com - shoebox_sz / 2., 0))
        j2 = int(min(y_com + shoebox_sz / 2., img_shape[1]))

        shoebox_img = img[j1:j2, i1:i2]
        shoebox_mask = is_bg_pixel[j1:j2, i1:i2]

        tilt, bgmask, coeff, _ = utils.tilting_plane(
            shoebox_img,
            mask=shoebox_mask,  # mask specifies which spots are bg pixels...
            zscore=2)

        tilt_abc[i_spot] = coeff[1], coeff[2], coeff[0]  # store as fast-scan coeff, slow-scan coeff, offset coeff

        spot_roi[i_spot] = i1, i2, j1, j2
        if plot:
            R = plt.Rectangle(xy=(x_com-shoebox_sz/2, y_com-shoebox_sz/2.),
                          width=shoebox_sz,
                          height=shoebox_sz,
                          fc='none', ec='r')
            patches.append(R)

    if plot:
        patch_coll = plt.mpl.collections.PatchCollection(patches, match_original=True)
        plt.imshow( img, vmin=0, vmax=200)
        plt.gca().add_collection(patch_coll)
        plt.show()

    return spot_roi, hkli, Areal_approx, Areal_GrnTru, tilt_abc, img, misset

if __name__=="__main__":
    out = process_simdata(plot=True)
