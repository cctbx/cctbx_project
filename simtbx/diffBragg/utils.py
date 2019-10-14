from scipy import ndimage
from scipy.optimize import minimize
import numpy as np
import pylab as plt


def get_spot_data(img, thresh=0, filter=None, **kwargs):
    """
    Kwargs are passed to the filter used to smear the spots
    :param img: numpy image, assumed to be simulated
    :param thresh: minimum value, this should be  >= the minimum intensity separating spots..
    :param filter: a filter to apply to the data, typically one of scipy.ndimage
        the kwargs will be passed along to this filter
    :return: useful spot dictionary, numpy version of a reflection table..
    """

    if filter is not None:
        labimg, nlab = ndimage.label( filter(img, **kwargs) > thresh)
    else:
        labimg, nlab = ndimage.label( img > thresh)

    if nlab == 0:
        return None

    bboxes = ndimage.find_objects(labimg)

    comIpos = ndimage.center_of_mass(img, labimg, range(1, nlab+1))
    maxIpos = ndimage.maximum_position(img, labimg, range(1, nlab+1))
    maxI = ndimage.maximum(img, labimg, range(1, nlab+1))
    meanI = ndimage.mean(img, labimg, range(1,nlab+1))
    varI = ndimage.variance(img, labimg, range(1,nlab+1))

    return {'comIpos': comIpos,
            'bboxes': bboxes,
            'maxIpos': maxIpos,
            'maxI': maxI,
            'meanI': meanI,
            'varI': varI}


def x_y_to_q(x,y, detector, beam):
    """
    :param x:  fast scan coordinate of spot
    :param y:  slow scan coordinate of spots
    :param detector:  dxtbx detector model
    :param beam:  dxtbx beam model
    :return: numpy array of q-vectors corresponding to the x,y coords
    """
    orig = np.array(detector[0].get_origin())
    fs = np.array(detector[0].get_fast_axis())
    ss = np.array(detector[0].get_slow_axis())
    fs_pixsize, ss_pixsize = detector[0].get_pixel_size()

    q_vecs = []
    for i_fs, i_ss in zip(x,y):
        s1 = orig + i_fs*fs*fs_pixsize + i_ss*ss*ss_pixsize  # scattering vector
        s1 = s1 / np.linalg.norm(s1) / beam.get_wavelength()
        q_vecs.append(s1-beam.get_s0())  # momentum transfer

    return np.vstack(q_vecs)


def is_outlier(points, thresh=3.5):
    """http://stackoverflow.com/a/22357811/2077270"""
    if len(points.shape) == 1:
        points = points[:, None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median) ** 2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh


def tilting_plane(img, mask=None, zscore=2 ):
    """
    :param img:  numpy image
    :param mask:  boolean mask, same shape as img, True is good pixels
        mask should include strong spot pixels and bad pixels, e.g. zingers
    :param zscore: modified z-score for outlier detection, lower increases number of outliers
    :return: tilting plane, same shape as img
    """
    Y, X = np.indices(img.shape)
    YY, XX = Y.ravel(), X.ravel()

    img1d = img.ravel()

    if mask is None:
        mask = np.ones(img.shape, bool)
    mask1d = mask.ravel()

    out1d = np.zeros(mask1d.shape, bool)
    out1d[mask1d] = is_outlier(img1d[mask1d].ravel(), zscore)
    out2d = out1d.reshape(img.shape)

    fit_sel = np.logical_and(~out2d, mask)  # fit plane to these points, no outliers, no masked
    x, y, z = X[fit_sel], Y[fit_sel], img[fit_sel]
    guess = np.array([np.ones_like(x), x, y]).T
    coeff, r, rank, s = np.linalg.lstsq(guess, z, rcond=-1)
    ev = (coeff[0] + coeff[1]*XX + coeff[2]*YY)
    return ev.reshape(img.shape), out2d, coeff, True


def _positive_plane(x, xcoord, ycoord, data):
    return np.sum((np.exp(x[0]) + np.exp(x[1])*xcoord + np.exp(x[2])*ycoord - data)**2)


def positive_tilting_plane(img, mask=None, zscore=2):
    """
    :param img:  numpy image
    :param mask:  boolean mask, same shape as img, True is good pixels
        mask should include strong spot pixels and bad pixels, e.g. zingers
    :param zscore: modified z-score for outlier detection, lower increases number of outliers
    :return: tilting plane, same shape as img
    """
    Y, X = np.indices(img.shape)
    YY, XX = Y.ravel(), X.ravel()

    img1d = img.ravel()

    if mask is None:
        mask = np.ones(img.shape, bool)
    mask1d = mask.ravel()

    out1d = np.zeros(mask1d.shape, bool)
    out1d[mask1d] = is_outlier(img1d[mask1d].ravel(), zscore)
    out2d = out1d.reshape(img.shape)

    fit_sel = np.logical_and(~out2d, mask)  # fit plane to these points, no outliers, no masked
    x, y, z = X[fit_sel], Y[fit_sel], img[fit_sel]
    out = minimize(
        fun=_positive_plane,
        x0=np.array([1e-2, 1e-2, 1e-2]),
        args=(x, y, z),
        method='Nelder-Mead')
        # bounds=((0, 1e10), (0, 1e10), (0, 1e10)),
        # method='L-BFGS-B')
    coeff = out.x
    ev = (np.exp(coeff[0]) + np.exp(coeff[1])*XX + np.exp(coeff[2])*YY)
    return ev.reshape(img.shape), out2d, coeff, out.success


def refine_model_from_angles(dxcryst, angles=(0, 0, 0)):
    from simtbx.diffBragg.nanoBragg_crystal import nanoBragg_crystal

    C = nanoBragg_crystal()

    C.dxtbx_crystal = dxcryst
    angles = np.array(angles) * 180 / np.pi
    from scitbx.matrix import col
    x = col((-1, 0, 0))
    y = col((0, -1, 0))
    z = col((0, 0, -1))

    RX = x.axis_and_angle_as_r3_rotation_matrix(angles[0], deg=True)
    RY = y.axis_and_angle_as_r3_rotation_matrix(angles[1], deg=True)
    RZ = z.axis_and_angle_as_r3_rotation_matrix(angles[2], deg=True)
    C.missetting_matrix = RX*RY*RZ

    dxcryst_refined = C.dxtbx_crystal_with_missetting()

    return dxcryst_refined


def lower_triangle(matrix_sqr):
    return [matrix_sqr[i] for i in [0, 3, 6, 4, 7, 8]]


def upper_triangle(matrix_sqr):
    return [matrix_sqr[i] for i in [0, 1, 2, 4, 5, 8]]


def get_diffBragg_instance():
    """
    Simple method to get a diffBRagg instance
    USED IN TESTING
    Returns an instance of diffBragg
    """
    from dxtbx.model.crystal import CrystalFactory
    from dxtbx.model.detector import DetectorFactory
    from dxtbx.model.beam import BeamFactory
    from simtbx.nanoBragg.tst_nanoBragg_basic import fcalc_from_pdb
    from simtbx.nanoBragg import shapetype
    from simtbx.diffBragg import diffBragg

    wavelen = 1.24
    flux = 1e15
    SHAPE = shapetype.Gauss

    NCELLS_ABC = 15  # means (15, 15, 15)

    beam_descr = {'direction': (0.0, 0.0, 1.0),
                  'divergence': 0.0,
                  'flux': 5e11,
                  'polarization_fraction': 1.,
                  'polarization_normal': (0.0, 1.0, 0.0),
                  'sigma_divergence': 0.0,
                  'transmission': 1.0,
                  'wavelength': wavelen}

    cryst_descr = {'__id__': 'crystal',
                   'real_space_a': (50, 0, 0),
                   'real_space_b': (0, 60, 0),
                   'real_space_c': (0, 0, 70),
                   'space_group_hall_symbol': '-P 4 2'}

    det_descr = {'panels':
                     [{'fast_axis': (-1.0, 0.0, 0.0),
                       'gain': 1.0,
                       'identifier': '',
                       'image_size': (196, 196),
                       'mask': [],
                       'material': '',
                       'mu': 0.0,
                       'name': 'Panel',
                       'origin': (19.6, -19.6, -550),
                       'pedestal': 0.0,
                       'pixel_size': (0.1, 0.1),
                       'px_mm_strategy': {'type': 'SimplePxMmStrategy'},
                       'raw_image_offset': (0, 0),
                       'slow_axis': (0.0, 1.0, 0.0),
                       'thickness': 0.0,
                       'trusted_range': (0.0, 65536.0),
                       'type': ''}]}

    DET = DetectorFactory.from_dict(det_descr)
    BEAM = BeamFactory.from_dict(beam_descr)
    crystal = CrystalFactory.from_dict(cryst_descr)

    Fhkl = fcalc_from_pdb(resolution=4, algorithm="fft", wavelength=wavelen)

    D = diffBragg(DET, BEAM, verbose=0, panel_id=0)
    D.xtal_shape = SHAPE
    D.Ncells_abc = NCELLS_ABC
    D.wavelength_A = wavelen
    D.flux = flux
    D.mosaic_spread_deg = 0.01
    D.mosaic_domains = 10
    D.Fhkl = Fhkl
    D.Bmatrix = crystal.get_B()
    D.Umatrix = crystal.get_U()
    return D


def process_simdata(spots, img, thresh=20, plot=False, shoebox_sz=20, edge_reflections=True):
    """
    This is a helper function for some of the tests
    :param spots:  simulated image without background/noise
    :param img: simulated image with background/noise
    :param thresh: spot threshold
    :param plot: whether to make plots
    :return: spot_rois as an Nx4 array and abc background planes as an Nx3 array
    """
    spot_data = get_spot_data(spots, thresh=thresh)
    ss_spot, fs_spot = map(np.array, zip(*spot_data["maxIpos"]))  # slow/fast  scan coords of strong spots
    num_spots = len(ss_spot)
    if plot:
        plt.imshow(img, vmax=200)
        plt.plot(fs_spot, ss_spot, 'o', mfc='none', mec='r')
        plt.title("Simulated image with strong spots marked")
        plt.xlim(-.5, img.shape[1]-.5)
        plt.ylim(img.shape[0]-.5, -.5)
        plt.show()

    is_bg_pixel = np.ones(img.shape, bool)
    for bb_ss, bb_fs in spot_data["bboxes"]:
        is_bg_pixel[bb_ss, bb_fs] = False

    # now fit tilting planes
    #tilt_abc = np.zeros((num_spots, 3))
    #spot_roi = np.zeros((num_spots, 4), int)
    tilt_abc = []
    spot_roi = []
    if plot:
        patches = []
    img_shape = img.shape  # TODO: verify fast slow scan
    successes = []
    for i_spot, (x_com, y_com) in enumerate(zip(fs_spot, ss_spot)):
        i1 = int(max(x_com - shoebox_sz / 2., 0))
        i2 = int(min(x_com + shoebox_sz / 2., img_shape[0]-1))
        j1 = int(max(y_com - shoebox_sz / 2., 0))
        j2 = int(min(y_com + shoebox_sz / 2., img_shape[1]-1))

        if not edge_reflections:
            if i2-i1 < shoebox_sz-1 or j2-j1 < shoebox_sz-1:
                print("Skipping reflection on the detector edge!")
                continue

        shoebox_img = img[j1:j2, i1:i2]
        shoebox_mask = is_bg_pixel[j1:j2, i1:i2]

        tilt, bgmask, coeff, success = tilting_plane(
            shoebox_img,
            mask=shoebox_mask,  # mask specifies which spots are bg pixels...
            zscore=2)
        success = True
        #tilt, bgmask, coeff, success = positive_tilting_plane(
        #    shoebox_img,
        #    mask=shoebox_mask,  # mask specifies which spots are bg pixels...
        #    zscore=2)
        successes.append(success)

        tilt_abc.append((coeff[1], coeff[2], coeff[0]))  # store as fast-scan coeff, slow-scan coeff, offset coeff
        #tilt_abc[i_spot] = coeff[1], coeff[2], coeff[0]  # store as fast-scan coeff, slow-scan coeff, offset coeff

        spot_roi.append((i1, i2, j1, j2))
        #spot_roi[i_spot] = i1, i2, j1, j2
        if plot:
            R = plt.Rectangle(xy=(x_com-shoebox_sz/2, y_com-shoebox_sz/2.),
                          width=shoebox_sz,
                          height=shoebox_sz,
                          fc='none', ec='r')
            patches.append(R)

    if plot:
        patch_coll = plt.mpl.collections.PatchCollection(patches, match_original=True)
        plt.imshow(img, vmin=0, vmax=200)
        plt.gca().add_collection(patch_coll)
        plt.show()

    spot_roi = np.array(spot_roi).astype(int)
    tilt_abc = np.array(tilt_abc)
    spot_roi = spot_roi[successes]
    tilt_abc = tilt_abc[successes]
    return spot_roi, tilt_abc

