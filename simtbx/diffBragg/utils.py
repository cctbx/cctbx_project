from scipy import ndimage
import numpy as np


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
    from cxid9114 import utils
    Y, X = np.indices(img.shape)
    YY, XX = Y.ravel(), X.ravel()

    img1d = img.ravel()

    if mask is None:
        mask = np.ones(img.shape, bool)
    mask1d = mask.ravel()

    out1d = np.zeros(mask1d.shape, bool)
    out1d[mask1d] = utils.is_outlier(img1d[mask1d].ravel(), zscore)
    out2d = out1d.reshape(img.shape)

    fit_sel = np.logical_and(~out2d, mask)  # fit plane to these points, no outliers, no masked
    x, y, z = X[fit_sel], Y[fit_sel], img[fit_sel]
    guess = np.array([np.ones_like(x), x, y ] ).T
    coeff, r, rank, s = np.linalg.lstsq(guess, z)
    ev = (coeff[0] + coeff[1]*XX + coeff[2]*YY )
    return ev.reshape(img.shape), out2d, coeff


def refine_model_from_angles(dxcryst, angles=(0, 0, 0)):
    from simtbx.diffBragg.nanoBragg_crystal import nanoBragg_crystal
    from copy import deepcopy

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

    dxcryst_refined = deepcopy(dxcryst)
    dxcryst_refined.set_A(C.Amatrix_realspace.inverse())

    return dxcryst_refined

