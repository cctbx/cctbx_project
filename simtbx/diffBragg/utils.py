from __future__ import absolute_import, division, print_function
import re
import os

from scipy import fft
from scipy import optimize
from scipy.ndimage import binary_dilation
from itertools import zip_longest
import math
import pickle
from simtbx.diffBragg.refiners.crystal_systems import OrthorhombicManager, TetragonalManager, MonoclinicManager, HexagonalManager
from scipy.spatial import cKDTree
from scipy.optimize import minimize
from dials.algorithms.image.filter import convolve
from scipy.interpolate import SmoothBivariateSpline
from cctbx.array_family import flex
from cctbx import miller, sgtbx
from cctbx.crystal import symmetry

import numpy as np
from simtbx.nanoBragg.utils import ENERGY_CONV
from dxtbx.model import Detector, Panel
from scipy import ndimage
from simtbx.nanoBragg.sim_data import SimData
from simtbx.nanoBragg.nanoBragg_beam import NBbeam
from simtbx.nanoBragg.nanoBragg_crystal import NBcrystal
from dxtbx.imageset import MemReader
from dxtbx.imageset import ImageSet, ImageSetData
from dxtbx.model.experiment_list import ExperimentListFactory
from dxtbx.model import ExperimentList

import logging
MAIN_LOGGER = logging.getLogger("main")


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


def label_background_pixels(roi_img, thresh=3.5, iterations=1, only_high=True):
    """
    roi_img, a two dimensional pixel image numpy array
    thresh, median deviation Z-score; pixels with Z-score below thresh are flagged as background
    iterations, lots of bright pixels can skew the stats, so iteratively label background pixels
              and recompute Z-scores
    """
    img_shape = roi_img.shape
    img1 = roi_img.copy().ravel()   # 1-D version
    background_pixels = None
    while iterations > 0:
        if background_pixels is None:
            outliers = is_outlier(img1, thresh)
            m = np.median(img1[~outliers])
            if only_high:
                outliers = np.logical_and(outliers, img1 > m)
            background_pixels = ~outliers
        else:
            where_bg = np.where(background_pixels)[0]
            outliers = is_outlier(img1[background_pixels], thresh)
            m = np.median(img1[background_pixels][~outliers])
            if only_high:
                outliers = np.logical_and(outliers, img1[background_pixels] > m)
            background_pixels[where_bg[outliers]] = False
        iterations = iterations - 1

    return background_pixels.reshape(img_shape)


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


def tilting_plane(img, mask=None, zscore=np.inf, spline=False,
                old_style=True, return_resid=False, return_error=False,
                integration=False, sigma_readout=0):
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
    if old_style:
        coeff, r, rank, s = np.linalg.lstsq(guess, z, rcond=-1)
        ev = (coeff[0] + coeff[1]*XX + coeff[2]*YY)
        coeff_errors = None,None,None
        if spline:
            sp = SmoothBivariateSpline(x, y, z, kx=1, ky=1)
            tilt = sp.ev(XX, YY).reshape(img.shape)
        else:
            tilt = ev.reshape(img.shape)
    else:
        rho = z
        W = np.diag(1/( sigma_readout**2 + rho))
        A = guess
        AWA = np.dot(A.T, np.dot(W, A))
        AWA_inv = np.linalg.inv(AWA)
        AtW = np.dot( A.T,W)
        a,b,c = np.dot(np.dot(AWA_inv, AtW), rho)
        coeff = (a,b,c)
        tilt = (XX*b + YY*c + a).reshape(img.shape)

        # vector of residuals
        r = y-np.dot(A,(a,b,c))
        r_fact = np.dot(r.T,np.dot( W,r)) / (len(rho)-3)
        error = AWA_inv * r_fact
        coeff_errors = np.sqrt(error[0][0]), np.sqrt(error[1][1]), np.sqrt(error[2][2])

    return_packet = [tilt, out2d, coeff, True]
    if return_resid:
        return_packet.append(r)
    if return_error:
        return_packet.append(coeff_errors)
    return return_packet


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
    from simtbx.nanoBragg.nanoBragg_crystal import NBcrystal

    C = NBcrystal()

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

    D = diffBragg(DET, BEAM, verbose=0)
    D.xtal_shape = SHAPE
    D.Ncells_abc = NCELLS_ABC
    D.wavelength_A = wavelen
    D.flux = flux
    D.mosaic_spread_deg = 0.01
    D.mosaic_domains = 10
    from simtbx.nanoBragg.nanoBragg_crystal import NBcrystal
    braggC = NBcrystal()
    braggC.miller_array = Fhkl
    idx = braggC.miller_array.indices()
    amps = braggC.miller_array.data()
    D.Fhkl_tuple = idx, amps, None
    D.Bmatrix = crystal.get_B()
    D.Umatrix = crystal.get_U()
    return D


def perturb_miller_array(F, factor, perturb_log_vals=True):
    if not F.is_xray_amplitude_array():
        F = F.amplitudes()
    Fdat = np.array( F.data())

    if perturb_log_vals:
        Fdat = np.log(Fdat)

    try:
        Fdat = np.random.uniform(Fdat-factor, Fdat+factor)
    except OverflowError:
        return None
    if perturb_log_vals:
        Fdat = np.exp(Fdat)

    Fdat = flex.double(np.ascontiguousarray(Fdat))

    mset = miller.set(F.crystal_symmetry(), indices=F.indices(), anomalous_flag=True)
    return miller.array(miller_set=mset, data=Fdat).set_observation_type_xray_amplitude()


def map_hkl_list(Hi_lst, anomalous_flag=True, symbol="P43212"):
    from cctbx import sgtbx
    from dials.array_family import flex as dials_flex
    sg_type = sgtbx.space_group_info(symbol=symbol).type()
    # necessary for py3 to type cast the ints
    type_casted_Hi_lst = tuple([(int(x), int(y), int(z)) for x, y, z in Hi_lst])
    Hi_flex = dials_flex.miller_index(type_casted_Hi_lst)
    miller.map_to_asu(sg_type, anomalous_flag, Hi_flex)
    return list(Hi_flex)


def compare_with_ground_truth(a, b, c, dxcryst_models, symbol="C121", verbose=False):
    """
    ported from LS49
    :param a: ground truth a (real space)
    :param b: ground truth b (real space)
    :param c: gournd truth c (real space)
    :param dxcryst_models: P1 models or whatever
    :param symbol: target space group
    :return: list of angles, one for each model
    """
    from cctbx import sgtbx
    from dxtbx.model import MosaicCrystalSauter2014
    from cctbx import crystal_orientation
    from libtbx.test_utils import approx_equal
    from scitbx.array_family import flex
    from scitbx.matrix import sqr

    sgi = sgtbx.space_group_info(symbol=symbol)
    CB_OP_C_P = sgi.change_of_basis_op_to_primitive_setting()

    icount = 0
    angles = []
    icount += 1
    for crystal_model in dxcryst_models:
        if crystal_model.get_space_group().info().type().lookup_symbol() == "P 1":
            crystal_model = crystal_model.change_basis(CB_OP_C_P.inverse())
        mosaic_model = MosaicCrystalSauter2014(crystal_model)

        direct_A = mosaic_model.get_A_inverse_as_sqr()
        sim_compatible = direct_A
        integrated_Ori = crystal_orientation.crystal_orientation(
            sim_compatible,
            crystal_orientation.basis_type.direct)

        header_Ori = crystal_orientation.crystal_orientation(
            tuple(a)+tuple(b)+tuple(c),
            crystal_orientation.basis_type.direct)

        C2_ground_truth = header_Ori.change_basis(CB_OP_C_P.inverse())
        if verbose:
            C2_ground_truth.show(legend="C2_ground_truth")

        # align integrated model with ground truth
        cb_op_align = integrated_Ori.best_similarity_transformation(C2_ground_truth, 50, 1)
        aligned_Ori = integrated_Ori.change_basis(sqr(cb_op_align))
        if verbose:
            aligned_Ori.show(legend="integrated, aligned")
            MAIN_LOGGER.debug("alignment matrix", cb_op_align)

        U_integrated = aligned_Ori.get_U_as_sqr()
        U_ground_truth = C2_ground_truth.get_U_as_sqr()
        missetting_rot = U_integrated * U_ground_truth.inverse()
        assert approx_equal(missetting_rot.determinant(), 1.0)

        # now calculate the angle as mean a_to_a,b_to_b,c_to_c
        aoff = aligned_Ori.a.angle(C2_ground_truth.a, deg=True)
        boff = aligned_Ori.b.angle(C2_ground_truth.b, deg=True)
        coff = aligned_Ori.c.angle(C2_ground_truth.c, deg=True)

        hyp = flex.mean(flex.double((aoff, boff, coff)))
        angles.append(hyp)

    return angles


def fcalc_from_pdb(resolution, algorithm=None, wavelength=0.9, anom=True, ucell=None, symbol=None, as_amplitudes=True):
    # borrowed from tst_nanoBragg_basic
    pdb_lines = """HEADER TEST
CRYST1   50.000   60.000   70.000  90.00  90.00  90.00 P 1
ATOM      1  O   HOH A   1      56.829   2.920  55.702  1.00 20.00           O
ATOM      2  O   HOH A   2      49.515  35.149  37.665  1.00 20.00           O
ATOM      3  O   HOH A   3      52.667  17.794  69.925  1.00 20.00           O
ATOM      4  O   HOH A   4      40.986  20.409  18.309  1.00 20.00           O
ATOM      5  O   HOH A   5      46.896  37.790  41.629  1.00 20.00           O
ATOM      6 SED  MSE A   6       1.000   2.000   3.000  1.00 20.00          SE
END
"""
    from iotbx import pdb
    pdb_inp = pdb.input(source_info=None, lines=pdb_lines)
    xray_structure = pdb_inp.xray_structure_simple()
    if ucell is not None:
        assert symbol is not None
        from cctbx.xray import structure
        from cctbx import crystal
        crystal_sym = crystal.symmetry(unit_cell=ucell, space_group_symbol=symbol)
        xray_structure = structure(scatterers=xray_structure.scatterers(), crystal_symmetry=crystal_sym)
    # take a detour to insist on calculating anomalous contribution of every atom
    scatterers = xray_structure.scatterers()
    if anom:
        from cctbx.eltbx import henke
        for sc in scatterers:
            expected_henke = henke.table(sc.element_symbol()).at_angstrom(wavelength)
            sc.fp = expected_henke.fp()
            sc.fdp = expected_henke.fdp()
    # how do we do bulk solvent?
    primitive_xray_structure = xray_structure.primitive_setting()
    P1_primitive_xray_structure = primitive_xray_structure.expand_to_p1()
    fcalc = P1_primitive_xray_structure.structure_factors(
        d_min=resolution, anomalous_flag=anom, algorithm=algorithm).f_calc()
    if as_amplitudes:
        fcalc = fcalc.amplitudes().set_observation_type_xray_amplitude()
    return fcalc


def nearest_non_zero(lst, idx):
    # https: // codereview.stackexchange.com / a / 172121 / 78230
    if lst[idx] > 0:
        return lst[idx]
    before, after = lst[:idx], lst[idx+1:]
    for b_val, a_val in zip_longest(reversed(before), after, fillvalue=0):
        # N.B. I applied `reversed` inside `zip_longest` here. This
        # ensures that `before` and `after` are the same type, and that
        # `before + [lst[idx]] + after == lst`.
        if b_val > 0:
            return b_val
        if a_val > 0:
            return a_val
    else:
        return 0  # all zeroes in this list


def update_miller_array_at_indices(miller_array, indices, new_values):
    if not miller_array.is_xray_amplitude_array():
        raise ValueError("Miller array is assumed to be an amplitude array")

    if len(indices) != len(new_values):
        raise ValueError("Indices (length %d) should be the same length as new_values (length %d)"
                         % (len(indices), len(new_values)))
    mset = miller_array.set()
    fdata_map = {tuple(h): val for h,val in zip(miller_array.indices(), miller_array.data())}
    for i_h, h in enumerate(indices):
        equivs = [idx.h() for idx in miller.sym_equiv_indices(mset.space_group(), h).indices()]
        n_missing = 0
        for equiv_h in equivs:
            if equiv_h in fdata_map:
                fdata_map[h] = new_values[i_h]
            else:
                n_missing += 1
        if n_missing == len(equivs):
            raise KeyError("Trying to update index %s which is not in the miller array" % " ".join(h))

    new_data = flex.double(list(fdata_map.values()))
    new_indices = flex.miller_index(tuple(fdata_map.keys()))
    new_mset = miller.set(mset.crystal_symmetry(), indices=new_indices, anomalous_flag=mset.anomalous_flag())
    new_miller_array = miller.array(miller_set=new_mset, data=new_data).set_observation_type_xray_amplitude()
    return new_miller_array


def fiber2D_integ(x,y,g):
    return math.atan((x*y)/(g*math.sqrt(g*g + x*x + y*y)))/(2.0*math.pi)


def makeMoffat_integPSF(fwhm_pixel, sizeX, sizeY):
    ''' Integral form of contribution of moffat PSF to signal recorded in a square pixel'''
    g = fwhm_pixel*0.65238
    psf = np.zeros((sizeX, sizeY))
    sx = int(sizeX/2)
    sy = int(sizeY/2)
    for y in range(-sy, sy+1):
      for x in range(-sx, sx+1):
        # Holton 2012 paper says x,y should be pixel center; this does not seem right ?
        psf[x+sx,y+sy] = fiber2D_integ(x+1./2,y+1./2,g)-fiber2D_integ(x+1./2,y-1./2,g)-fiber2D_integ(x-1./2,y+1./2,g)+fiber2D_integ(x-1./2,y-1./2,g)
        #psf[x+sx, -y+sy] = fiber2D_integ(x+1./2,y+1./2,g)-fiber2D_integ(x+1./2,y-1./2,g)-fiber2D_integ(x-1./2,y+1./2,g)+fiber2D_integ(x-1./2,y-1./2,g)
        # Trying to get pixel center instead
        #psf[x+sx, -y+sy] = fiber2D_integ(x+1,y+1,g)-fiber2D_integ(x+1,y,g)-fiber2D_integ(x,y+1,g)+fiber2D_integ(x,y,g)
    psf = psf/psf.sum()
    psf = psf.tolist()
    psf = flex.double(psf)
    return psf


def convolve_padded_img(img, psf, sz=5):
    img = np.array(img)
    iY, iX = img.shape
    pY, pX = psf.focus()

    new_iY = iY
    if pY >= iY - sz:
        new_iY = pY + sz
    new_iX = iX
    if pX >= iX - sz:
        new_iX = pX + sz

    assert new_iX >= iX
    assert new_iY >= iY
    padX = new_iX - iX
    padY = new_iY - iY

    x = int(padX/2)
    y = int(padY/2)
    img = np.pad(img, ((y, y+1), (x, x+1)), mode='median')
    assert img.shape[0] >= pY + sz
    assert img.shape[1] >= pX + sz

    conv_img = convolve(flex.double(img), psf)
    conv_img = conv_img.as_numpy_array()[y:y+iY, x:x+iX]
    return conv_img


def convolve_with_psf(image_data, fwhm=27.0, pixel_size=177.8, psf_radius=7, sz=5, psf=None):
    ''' Given a 2D numpy array of image data, convolve with a PSF. '''
    # Currently only supporting fiber PSF i.e power law form as proposed in Holton et. al 2012, Journal of Synchotron Radiation
    if psf is None:
        xpsf=2*psf_radius+1
        ypsf=2*psf_radius+1
        fwhm_pixel=fwhm/pixel_size
        psf = makeMoffat_integPSF(fwhm_pixel, xpsf, ypsf)
    img_shape = image_data.shape
    psf_shape = psf.focus()
    if psf_shape[0] > img_shape[0] - sz or psf_shape[1] > img_shape[1] - sz:
        convolved_image = convolve_padded_img(image_data, psf, sz)
    else:
        convolved_image = convolve(flex.double(image_data), psf)
        convolved_image = convolved_image.as_numpy_array()
    return convolved_image


def get_roi_from_spot(refls, fdim, sdim, shoebox_sz=10):
    fs_spot, ss_spot, _ = zip(*refls['xyzobs.px.value'])
    rois = []
    is_on_edge = []
    for i_spot, (x_com, y_com) in enumerate(zip(fs_spot, ss_spot)):
        x_com = x_com - 0.5
        y_com = y_com - 0.5
        i1 = int(max(x_com - shoebox_sz / 2., 0))
        i2 = int(min(x_com + shoebox_sz / 2., fdim-1))
        j1 = int(max(y_com - shoebox_sz / 2., 0))
        j2 = int(min(y_com + shoebox_sz / 2., sdim-1))

        if i2-i1 < shoebox_sz-1 or j2-j1 < shoebox_sz-1:
            is_on_edge.append(True)
        else:
            is_on_edge.append(False)
        roi = i1, i2, j1, j2
        rois.append(roi)
    return rois, is_on_edge


def add_rlp_column(refls, experiment):
    keys = list(refls[0].keys())
    if "s1" in keys:
        s1 = refls['s1']
        s1_norm = np.array(s1) / np.linalg.norm(s1,axis=1)[:,None]
        wavelen = experiment.beam.get_wavelength()
        s0 = np.array([list(experiment.beam.get_unit_s0())]*len(refls))
        q_vecs = 1/wavelen*(s1_norm-s0)
        refls['rlp'] = flex.vec3_double(tuple(map(tuple, q_vecs)))
    else:
        raise KeyError("Need rlp or s1 column in refl table!")


def get_roi_deltaQ(refls, delta_Q, experiment):
    """
    :param refls: reflection table (needs rlp column)
    :param delta_Q:  width of the ROI in inverse Angstromg (e.g. 0.05)
    :param experiment:
    :return:
    """
    nref = len(refls)
    assert nref >0

    keys = list(refls[0].keys())
    if "rlp" not in  keys:
        add_rlp_column(refls, experiment)

    beam = experiment.beam
    detector = experiment.detector
    assert beam is not None and experiment is not None

    rois = []
    is_on_edge = []
    for i_refl in range(nref):
        roi, on_edge = determine_shoebox_ROI(detector, delta_Q, beam.get_wavelength(), refls[i_refl])
        rois.append(roi)
        is_on_edge.append( on_edge)
    return rois, is_on_edge


def get_roi_background_and_selection_flags(refls, imgs, shoebox_sz=10, reject_edge_reflections=False,
                                   reject_roi_with_hotpix=True, background_mask=None, hotpix_mask=None,
                                   bg_thresh=3.5, set_negative_bg_to_zero=False,
                                   pad_for_background_estimation=None, use_robust_estimation=True, sigma_rdout=3.,
                                   min_trusted_pix_per_roi=4, deltaQ=None, experiment=None, weighted_fit=True,
                                   tilt_relative_to_corner=False, ret_cov=False, allow_overlaps=False):

    # TODO handle divide by 0 warning that happens in is_outlier, when labeling background pix?
    npan, sdim, fdim = imgs.shape

    if hotpix_mask is not None:
        assert hotpix_mask.shape == imgs.shape

    if background_mask is not None:
        assert background_mask.shape == imgs.shape

    if deltaQ is None:  # then use fixed size ROIS determined by shoebox_sz
        rois, is_on_edge = get_roi_from_spot(refls, fdim, sdim, shoebox_sz=shoebox_sz)
    else:
        assert experiment is not None
        if len(refls) == 0:
            return
        rois, is_on_edge = get_roi_deltaQ(refls, deltaQ, experiment)

    #for i_ref, (x1,x2,y1,y2) in enumerate(rois):
    #    x,y,_ = refls[i_ref]["xyzobs.px.value"]
    #    assert x1 <= x <= x2, "refl %d, %f %f %f" % (i_ref, x1, x, x2)
    #    assert y1 <= y <= y2, "refl %d, %f %f %f" % (i_ref, y1, y, y2)


    tilt_abc = []
    kept_rois = []
    panel_ids = []
    all_cov =[]
    selection_flags = []
    num_roi_negative_bg = 0
    num_roi_nan_bg = 0
    background = np.ones(imgs.shape)*-1
    i_roi = 0
    #import pylab as ply
    #patches = []
    while i_roi < len(rois):
    #for i_roi, roi in enumerate(rois):
        roi = rois[i_roi]
        i1, i2, j1, j2 = roi
        #rect = plt.Rectangle(xy=(i1,j1), width=i2-i1, height=j2-j1, fc='none', ec='w')
        is_selected = True
        if is_on_edge[i_roi] and reject_edge_reflections:
            MAIN_LOGGER.debug("Reflection %d bounded by x1=%d,x2=%d,y1=%d,y2=%d is on edge" % (i_roi, i1,i1,j2,j2))
            is_selected = False
        pid = refls[i_roi]['panel']

        if hotpix_mask is not None:
            is_hotpix = hotpix_mask[pid, j1:j2, i1:i2]
            num_hotpix = is_hotpix.sum()
            if num_hotpix > 0 and reject_roi_with_hotpix:
                MAIN_LOGGER.debug("reflection %d has hot pixel" % i_roi)
                is_selected = False
            if num_hotpix > min_trusted_pix_per_roi:
                is_selected = False

        dimY, dimX = imgs[pid].shape
        j1_nopad = j1
        i1_nopad = i1
        j2_nopad = j2
        i2_nopad = i2
        if pad_for_background_estimation is not None:
            j1 = max(0, j1-pad_for_background_estimation)
            i1 = max(0, i1-pad_for_background_estimation)
            j2 = min(dimY, j2+pad_for_background_estimation)
            i2 = min(dimX, i2+pad_for_background_estimation)

        shoebox = imgs[pid, j1:j2, i1:i2]

        if not isinstance(sigma_rdout, float) and not isinstance(sigma_rdout, int):
            shoebox_sigma_readout = sigma_rdout[pid, j1:j2, i1:i2]
        else:
            shoebox_sigma_readout = sigma_rdout

        if background_mask is not None:
            is_background = background_mask[pid, j1:j2, i1:i2]
        else:
            is_background = label_background_pixels(shoebox,thresh=bg_thresh, iterations=2)

        Ycoords, Xcoords = np.indices((j2-j1, i2-i1))

        if use_robust_estimation:
            bg_pixels = shoebox[is_background]
            bg_signal = np.median(bg_pixels)
            if bg_signal < 0:
                num_roi_negative_bg += 1
                if set_negative_bg_to_zero:
                    bg_signal = 0
                else:
                    is_selected = False
            #if i_roi==0:
            #    bg_signal = 3.1936748
            #if i_roi==1:
            #    bg_signal = 5.1090536
            tilt_a, tilt_b, tilt_c = 0, 0, bg_signal
            covariance = None

            tilt_plane = tilt_a * Xcoords + tilt_b * Ycoords + tilt_c
        else:
            if tilt_relative_to_corner:
                Xcoords += i1
                Ycoords += j1
            fit_results = fit_plane_equation_to_background_pixels(
                shoebox_img=shoebox, fit_sel=is_background, sigma_rdout=shoebox_sigma_readout,
                weighted=weighted_fit)
            if fit_results is None:
                tilt_a = tilt_b = tilt_c = covariance = 0
                is_selected = False
                MAIN_LOGGER.debug("tilt fit failed for reflection %d, probably too few pixels" % i_roi)
                tilt_plane = np.zeros_like(Xcoords)
            else:
                (tilt_a, tilt_b, tilt_c), covariance = fit_results
                #if XY is not None:
                #    tilt_plane = tilt_a * (Xcoords+XY[0]) + tilt_b * (Ycoords+XY[1]) + tilt_c
                #else:
                tilt_plane = tilt_a * Xcoords + tilt_b * Ycoords + tilt_c
                if np.any(np.isnan(tilt_plane)) and is_selected:
                    is_selected = False
                    num_roi_nan_bg += 1
                if np.min(tilt_plane) < 0:  # dips below
                    num_roi_negative_bg += 1
                    is_selected = False

        is_overlapping = not np.all(background[pid, j1_nopad:j2_nopad, i1_nopad:i2_nopad] == -1)
        if not allow_overlaps and is_overlapping:
            # NOTE : move away from this option, it potentially moves the pixel centroid outside of the ROI (in very rare instances)
            MAIN_LOGGER.debug("region of interest already accounted for roi size= %d %d" % (i2_nopad-i1_nopad, j2_nopad-j1_nopad))
            rois[i_roi] = i1_nopad+1,i2_nopad,j1_nopad+1,j2_nopad
            continue

        # unpadded ROI dimension
        roi_dimY = j2_nopad-j1_nopad
        roi_dimX = i2_nopad-i1_nopad

        if roi_dimY < 2 or roi_dimX < 2:
            is_selected = False
        #assert np.all(background[pid, j1_nopad:j2_nopad, i1_nopad:i2_nopad] == -1), "region of interest already accounted for"

        j1 = j1_nopad-j1
        i1 = i1_nopad-i1
        plane_in_unpadded_roi = tilt_plane[j1: j1 + roi_dimY, i1: i1 + roi_dimX]
        background[pid, j1_nopad:j2_nopad, i1_nopad:i2_nopad] = plane_in_unpadded_roi
        tilt_abc.append((tilt_a, tilt_b, tilt_c))
        all_cov.append(covariance)
        kept_rois.append(roi)
        panel_ids.append(pid)
        #patches.append(rect)
        selection_flags.append(is_selected)
        #if pid == 141 and 150 < ii < 160 and 112 < jj < 122:
        #    from IPython import embed
        #    embed()
        i_roi += 1

    #assert len(kept_rois) == len(refls)
    #for i_ref, (x1,x2,y1,y2) in enumerate(kept_rois):
    #    x,y,_ = refls[i_ref]["xyzobs.px.value"]
    #    assert x1 <= x <= x2, "refl %d, %f %f %f" % (i_ref, x1, x, x2)
    #    assert y1 <= y <= y2, "refl %d, %f %f %f" % (i_ref, y1, y, y2)

    #plt.imshow(imgs[0], vmax=100)
    #for p in patches:
    #    plt.gca().add_patch(p)
    #plt.show()
    MAIN_LOGGER.debug("Number of ROI with negative BGs: %d / %d" % (num_roi_negative_bg, len(rois)))
    MAIN_LOGGER.debug("Number of ROI with NAN in BGs: %d / %d" % (num_roi_nan_bg, len(rois)))
    if ret_cov:
        return kept_rois, panel_ids, tilt_abc, selection_flags, background, all_cov
    else:
        return kept_rois, panel_ids, tilt_abc, selection_flags, background


def determine_shoebox_ROI(detector, delta_Q, wavelength_A, refl, centroid="obs"):
    panel = detector[refl["panel"]]
    width = get_width_of_integration_shoebox(panel, delta_Q, wavelength_A, refl['rlp'])
    wby2 = int(round(width/2.))
    fdim,sdim = panel.get_image_size()
    if centroid=='obs':  #TODO: give this more options
        i_com, j_com,_ = refl['xyzobs.px.value']
    else:
        i1,i2,j1,j2,_,_ = refl['bbox']  # super weird funky spots can skew bbox such that its a bad measure of centroid
        i_com = (i1+i2) * .5
        j_com = (j1+j2) * .5
    i1, i2, j1, j2 = i_com - wby2, i_com + wby2, j_com - wby2, j_com + wby2
    on_edge = False
    if i1 < 0 or j1 < 0:
        on_edge = True
    elif i2 >= fdim  or j2 >= sdim:
        on_edge = True
    i1 = max(int(i1), 0)
    i2 = min(int(i2), fdim)
    j1 = max(int(j1), 0)
    j2 = min(int(j2), sdim)
    return (i1, i2, j1, j2), on_edge


def get_width_of_integration_shoebox(detector_panel, delta_Q, ave_wavelength_A,  rlp):
    Qmag = 2 * np.pi * np.linalg.norm(rlp)  # magnitude momentum transfer of the RLP in physicist convention
    detdist_mm = detector_panel.get_distance()
    pixsize_mm = detector_panel.get_pixel_size()[0]
    rad1 = (detdist_mm / pixsize_mm) * np.tan(2 * np.arcsin((Qmag - delta_Q * .5) * ave_wavelength_A / 4 / np.pi))
    rad2 = (detdist_mm / pixsize_mm) * np.tan(2 * np.arcsin((Qmag + delta_Q * .5) * ave_wavelength_A / 4 / np.pi))
    bbox_extent = (rad2 - rad1) / np.sqrt(2)  # rad2 - rad1 is the diagonal across the bbox
    return bbox_extent


def check_if_tilt_dips_below_zero(tilt_coefs, dim_slowfast):
    slow_dim, fast_dim = dim_slowfast
    tX, tY, tZ = tilt_coefs
    tilt_plane = fast_dim*tX + slow_dim*tY + tZ
    dips_below_zero = False
    if np.min(tilt_plane) < 0:
        dips_below_zero = True
    return dips_below_zero


def fit_plane_equation_to_background_pixels(shoebox_img, fit_sel, sigma_rdout=3, weighted=True, relative_XY=None):
    """
    :param shoebox_img: 2D pixel image (usually less than 100x100 pixels)
    :param fit_sel: pixels to be fit to plane (usually background pixels)
    :param sigma_rdout: float or 2D image , readout noise term ( in shoebox_img pixel units)
    :return: coefficients of tilt plane (fast-scan, slow-scan, offset), as well as a covariance estimate matrix
    """
    # fast scan pixels, slow scan pixels, pixel values (corrected for gain)
    Y, X = np.indices(shoebox_img.shape)
    if relative_XY is not None:
        X += relative_XY[0]
        Y += relative_XY[1]
    fast, slow, rho_bg = X[fit_sel], Y[fit_sel], shoebox_img[fit_sel]


    try:
        sigma_rdout = sigma_rdout[fit_sel]
    except TypeError:
        pass

    # do the fit of the background plane
    A = np.array([fast, slow, np.ones_like(fast)]).T
    # weights matrix:
    if weighted:
        W = np.diag(1 / (sigma_rdout ** 2 + rho_bg))
        W[W <0] = 0
        W = np.sqrt(W)
    else:
        W = np.eye(len(rho_bg))
    AWA = np.dot(A.T, np.dot(W, A))
    try:
        AWA_inv = np.linalg.inv(AWA)
    except np.linalg.LinAlgError:
        MAIN_LOGGER.warning("WARNING: Fit did not work.. investigate reflection, number of pixels used in fit=%d" % len(fast))
        return None
    AtW = np.dot(A.T, W)
    t1, t2, t3 = np.dot(np.dot(AWA_inv, AtW), rho_bg)

    # get the covariance of the tilt coefficients:
    # vector of residuals
    r = rho_bg - np.dot(A, (t1, t2, t3))
    Nbg = len(rho_bg)
    r_fact = np.dot(r.T, np.dot(W, r)) / (Nbg - 3)  # 3 parameters fit
    var_covar = AWA_inv * r_fact  # TODO: check for correlations in the off diagonal elems

    return (t1, t2, t3), var_covar


def strip_thickness_from_detector(detector):
    """warning: this overrides detector hierarchy"""
    thin_detector = Detector()
    for i_pan in range(len(detector)):
        panel = detector[i_pan]
        panel_dict = panel.to_dict()
        # NOTE: necessary to set origin and fast/slow because panel.to_dict doesnt care about hierarchy
        panel_dict["origin"] = panel.get_origin()
        panel_dict["fast_axis"] = panel.get_fast_axis()
        panel_dict["slow_axis"] = panel.get_slow_axis()
        panel_dict["mu"] = 0
        panel_dict["thickness"] = 0
        panel_dict["px_mm_strategy"] = {'type': 'SimplePxMmStrategy'}
        thin_panel = Panel.from_dict(panel_dict)
        thin_detector.add_panel(thin_panel)
    return thin_detector


def get_14bit_from_jungfrau(expt):
    iset = expt.imageset
    F = iset.get_format_class()
    if len(iset.paths()) != 1:
        raise ValueError("imageset should have exactly 1 path")
    fclass = F.get_instance(iset.paths()[0])
    return fclass.get_14bit_component(iset.indices()[0])


def get_pedestalRMS_from_jungfrau(expt, gain_modes_too=False):
    iset = expt.imageset
    F = iset.get_format_class()
    if len(iset.paths()) != 1:
        raise ValueError("imageset should have exactly 1 path")
    fclass = F.get_instance(iset.paths()[0])
    return fclass.get_pedestal_rms(iset.indices()[0], return_gain_modes=gain_modes_too)

def gain_data_from_expt(expt):
    gainexp_path = expt.gainmap_path
    MAIN_LOGGER.info("Loading gainmap data from experiment %s" % gainexp_path)
    gainexp = ExperimentListFactory.from_json_file(gainexp_path, check_format=True)[0]
    iset = gainexp.imageset
    gain_data = [rawdat.as_numpy_array() >> 14 for rawdat in iset.get_raw_data(0)]
    return np.array(gain_data)


def image_data_from_expt(expt, as_double=True):
    iset = expt.imageset
    if len(iset) == 0:
        raise ValueError("imageset should have 1 shot")
    if len(iset) > 1:
        raise ValueError("imageset should have only 1 shot. This expt has imageset with %d shots" % len(iset))
    try:
        flex_data = iset.get_raw_data(0)
    except Exception as err:
        assert str(type(err)) == "<class 'Boost.Python.ArgumentError'>", "something weird going on with imageset data"
        flex_data = iset.get_raw_data()
    if not isinstance(flex_data, tuple):
        flex_data = (flex_data, )
    img_data = np.array([data.as_numpy_array() for data in flex_data])
    if as_double:
        img_data = img_data.astype(np.float64)
    return img_data


def simulator_from_expt_and_params(expt, params=None, oversample=0, device_id=0, init_scale=1, total_flux=1e12,
                                   ncells_abc=(10,10,10), has_isotropic_ncells=True, mosaicity=0, num_mosaicity_samples=1, mtz_name=None,
                                   mtz_column=None, default_F=0, dmin=1.5, dmax=30, spectra_file=None, spectra_stride=1,
                                   complex_F=None):

    if params is not None:
        oversample = params.simulator.oversample
        device_id = params.simulator.device_id
        init_scale = params.simulator.init_scale
        total_flux = params.simulator.total_flux

        ncells_abc = params.simulator.crystal.ncells_abc
        ncells_def = params.simulator.crystal.ncells_def
        has_isotropic_ncells = params.simulator.crystal.has_isotropic_ncells
        mosaicity = params.simulator.crystal.mosaicity
        num_mosaicity_samples = params.simulator.crystal.num_mosaicity_samples

        mtz_name = params.simulator.structure_factors.mtz_name
        mtz_column = params.simulator.structure_factors.mtz_column
        default_F = params.simulator.structure_factors.default_F
        dmin = params.simulator.structure_factors.dmin
        dmax = params.simulator.structure_factors.dmax

        spectra_file = params.simulator.spectrum.filename
        spectra_stride = params.simulator.spectrum.stride
        aniso_mos_spread = params.simulator.crystal.anisotropic_mosaicity
    else:
        ncells_def = None
        aniso_mos_spread = None

    if has_isotropic_ncells:
        if len(set(ncells_abc)) != 1 :
            raise ValueError("`isotropic_ncells=True`, so `ncells_abc` should all be the same, not %d %d %d" % tuple(ncells_abc))

    # make a simulator instance
    SIM = SimData()
    if params.simulator.detector.force_zero_thickness:
        SIM.detector = strip_thickness_from_detector(expt.detector)
    else:
        SIM.detector = expt.detector

    # create nanoBragg crystal
    crystal = NBcrystal()
    crystal.isotropic_ncells = has_isotropic_ncells
    if params.simulator.crystal.rotXYZ_ucell is not None:
        rotXYZ = params.simulator.crystal.rotXYZ_ucell[:3]
        ucell_p = params.simulator.crystal.rotXYZ_ucell[3:]
        crystal.dxtbx_crystal = modify_crystal(rotXYZ, ucell_p, expt.crystal)
    else:
        crystal.dxtbx_crystal = expt.crystal
    crystal.thick_mm = 0.1  # hard code a thickness, will be over-written by the scale
    crystal.Ncells_abc = tuple(ncells_abc)  #params.simulator.init_ncells_abc
    if ncells_def is not None:
        crystal.Ncells_def = tuple(ncells_def)  #params.simulator.init_ncells_abc
    crystal.anisotropic_mos_spread_deg = aniso_mos_spread
    crystal.n_mos_domains = num_mosaicity_samples
    crystal.mos_spread_deg = mosaicity
    if params is not None:
        crystal.mos_angles_per_axis = params.simulator.crystal.mos_angles_per_axis
        crystal.num_mos_axes = params.simulator.crystal.num_mos_axes
    if mtz_name is None:
        miller_data = make_miller_array(
            symbol=expt.crystal.get_space_group().info().type().lookup_symbol(),
            unit_cell=expt.crystal.get_unit_cell(), d_min=dmin, d_max=dmax)
    else:
        miller_data = open_mtz(mtz_name, mtz_column)
    if complex_F is not None:
        miller_data = complex_F

    crystal.miller_array = miller_data
    SIM.crystal = crystal
    #if params is not None:
    #    SIM.Umats_method = params.simulator.crystal.mosaicity_method

    # create a nanoBragg beam
    beam = NBbeam()
    beam.size_mm = params.simulator.beam.size_mm
    beam.unit_s0 = expt.beam.get_unit_s0()
    if spectra_file is not None:
        init_spectrum = load_spectra_file(spectra_file, total_flux, spectra_stride, as_spectrum=True)
    else:
        assert total_flux is not None
        init_spectrum = [(expt.beam.get_wavelength(), total_flux)]
    beam.spectrum = init_spectrum
    SIM.beam = beam

    SIM.panel_id = 0
    SIM.instantiate_diffBragg(oversample=oversample, device_Id=device_id, default_F=default_F)
    if init_scale is not None:
        #TODO phase this parameter out since its redundant?
        SIM.update_nanoBragg_instance("spot_scale", init_scale)
    return SIM


def get_flux_and_energy(beam=None, spec_file=None, total_flux=1e12, pinkstride=None):
    if spec_file is not None:
        FLUX, energies = load_spectra_file(spec_file, total_flux=total_flux, pinkstride=pinkstride)
    else:
        assert beam is not None
        FLUX = [total_flux]
        energies = [ENERGY_CONV / beam.get_wavelength()]

    return FLUX, energies


def open_mtz(mtzfname, mtzlabel=None, verbose=False):
    if mtzlabel is None:
        mtzlabel = "fobs(+)fobs(-)"
    if verbose:
        MAIN_LOGGER.info("Opening mtz file %s , label %s" % (mtzfname, mtzlabel))
    from iotbx.reflection_file_reader import any_reflection_file
    miller_arrays = any_reflection_file(mtzfname).as_miller_arrays()

    possible_labels = []
    foundlabel = False
    for ma in miller_arrays:
        label = ma.info().label_string()
        possible_labels.append(label)
        if label == mtzlabel:
            foundlabel = True
            break

    assert foundlabel, "MTZ Label not found... \npossible choices: %s" % (" ".join(possible_labels))
    ma = ma.as_amplitude_array()
    return ma


def make_miller_array(symbol, unit_cell, defaultF=1000, d_min=1.5, d_max=999):
    sgi = sgtbx.space_group_info(symbol)
    # TODO: allow override of ucell
    symm = symmetry(unit_cell=unit_cell, space_group_info=sgi)
    miller_set = symm.build_miller_set(anomalous_flag=True, d_min=d_min, d_max=d_max)
    # NOTE does build_miller_set automatically expand to p1 ? Does it obey systematic absences ?
    # Note how to handle sys absences here ?
    Famp = flex.double(np.ones(len(miller_set.indices())) * defaultF)
    mil_ar = miller.array(miller_set=miller_set, data=Famp).set_observation_type_xray_amplitude()
    return mil_ar

def save_spectra_file(spec_file, wavelengths, weights):
    data = np.array([wavelengths, weights])
    np.savetxt(spec_file, data.T, delimiter=',', header="wavelengths, weights")


def load_spectra_file(spec_file, total_flux=None, pinkstride=1, as_spectrum=False):
    wavelengths, weights = np.loadtxt(spec_file, float, delimiter=',', skiprows=1).T
    if isinstance(wavelengths, float) and isinstance(weights, float):
        # the file had one entry:
        wavelengths = np.array([wavelengths])
        weights = np.array([weights])
    if pinkstride > len(wavelengths) or pinkstride == 0:
        raise ValueError("Incorrect value for pinkstride")
    wavelengths = wavelengths[::pinkstride]
    weights = weights[::pinkstride]
    energies = ENERGY_CONV/wavelengths
    if total_flux is not None:
        weights = weights / weights.sum() * total_flux
    if as_spectrum:
        return list(zip(list(wavelengths), list(weights)))
    else:
        return weights, energies


def load_mask(maskfile):
    if maskfile is None:
        return None
    with open(maskfile, 'rb') as o:
        mask = pickle.load(o)
    if isinstance(mask, tuple):
        mask = np.array([m.as_numpy_array() for m in mask])
    else:
        mask = mask.as_numpy_array()
    return mask


def unitcell_sigmas(unitcell_manager, unitcell_sigmas):
    name_mapping = {'a_Ang': 0, 'b_Ang': 1, 'c_Ang': 2, 'alpha_rad': 3, 'beta_rad': 4, 'gamma_rad': 5}
    variable_sigmas = []
    for name in unitcell_manager.variable_names:
        sig = unitcell_sigmas[name_mapping[name]]
        variable_sigmas.append(sig)
    return variable_sigmas


def manager_from_crystal(crystal):
    """
    :param crystal:  dxtbx crystal model
    :return:
    """
    params = crystal.get_unit_cell().parameters()
    return manager_from_params(params)


def manager_from_params(ucell_p):

    a, b, c, al, be, ga = ucell_p
    if np.isclose(a,b) and not np.isclose(a,c) and np.allclose([al, be, ga], [90]*3):
        manager = TetragonalManager(a=a, c=c)

    elif not np.isclose(a,b) and not np.isclose(b,c) and not np.isclose(a,c) and np.allclose([al, be, ga], [90]*3):
        manager = OrthorhombicManager(a=a, b=b, c=c)

    elif np.isclose(a,b) and not np.isclose(a,c) and np.allclose([al, be], [90]*2) and np.allclose([ga], [120]):
        manager = HexagonalManager(a=a, c=c)

    elif not np.isclose(a,b) and not np.isclose(b,c) and not np.isclose(a,c) and np.allclose([al, ga], [90]*2) and not np.allclose([be], [120]):
        manager = MonoclinicManager(a=a, b=b, c=c, beta=be*np.pi/180.)

    else:
        raise NotImplementedError("Not yet implemented for crystal model")

    return manager


def refls_from_sims(panel_imgs, detector, beam, thresh=0, filter=None, panel_ids=None,
                    max_spot_size=1000, **kwargs):
    """
    This is for converting the centroids in the noiseless simtbx images
    to a multi panel reflection table

    :param panel_imgs: list or 3D array of detector panel simulations
    :param detector: dxtbx  detector model of a caspad
    :param beam:  dxtxb beam model
    :param thresh: threshol intensity for labeling centroids
    :param filter: optional filter to apply to images before
        labeling threshold, typically one of scipy.ndimage's filters
    :param pids: panel IDS , else assumes panel_imgs is same length as detector
    :param kwargs: kwargs to pass along to the optional filter
    :return: a reflection table of spot centroids
    """
    from dials.algorithms.spot_finding.factory import FilterRunner
    from dials.model.data import PixelListLabeller, PixelList
    from dials.algorithms.spot_finding.finder import pixel_list_to_reflection_table

    if panel_ids is None:
        panel_ids = np.arange(len(detector))
    pxlst_labs = []
    for i, pid in enumerate(panel_ids):
        plab = PixelListLabeller()
        img = panel_imgs[i]
        if filter is not None:
            mask = filter(img, **kwargs) > thresh
        else:
            mask = img > thresh
        img_sz = detector[int(pid)].get_image_size()  # for some reason the int cast is necessary in Py3
        flex_img = flex.double(img)
        flex_img.reshape(flex.grid(img_sz))

        flex_mask = flex.bool(mask)
        flex_mask.resize(flex.grid(img_sz))
        pl = PixelList(0, flex.double(img), flex.bool(mask))
        plab.add(pl)

        pxlst_labs.append(plab)

    El = explist_from_numpyarrays(panel_imgs, detector, beam)
    iset = El.imagesets()[0]
    refls = pixel_list_to_reflection_table(
        iset, pxlst_labs,
        min_spot_size=1,
        max_spot_size=max_spot_size,  # TODO: change this ?
        filter_spots=FilterRunner(),  # must use a dummie filter runner!
        write_hot_pixel_mask=False)[0]

    return refls


class FormatInMemory:
    """
    this class is a special image type
    necessary to create dxtbx imagesets and
    datablocks from numpy array images
    and masks.
    """
    def __init__(self, image, mask=None):
        self.image = image
        if image.dtype != np.float64:
            self.image = self.image.astype(np.float64)
        if mask is None:
            self.mask = np.ones_like(self.image).astype(np.bool)
        else:
            assert (mask.shape == image.shape)
            assert(mask.dtype == bool)
            self.mask = mask

    def get_raw_data(self):
        if len(self.image.shape)==2:
            return flex.double(self.image)
        else:
            return tuple([flex.double(panel) for panel in self.image])

    def get_mask(self, goniometer=None):
        if len(self.image.shape)==2:
            return flex.bool(self.mask)
        else:
            return tuple([flex.bool(panelmask) for panelmask in self.mask])


def explist_from_numpyarrays(image, detector, beam, mask=None):
    """
    So that one can do e.g.
    >> dblock = datablock_from_numpyarrays( image, detector, beam)
    >> refl = flex.reflection_table.from_observations(dblock, spot_finder_params)
    without having to utilize the harddisk

    :param image:  numpy array image, or list of numpy arrays
    :param mask:  numpy mask, should be same shape format as numpy array
    :param detector: dxtbx detector model
    :param beam: dxtbx beam model
    :return: datablock for the image
    """
    if isinstance( image, list):
        image = np.array( image)
    if mask is not None:
        if isinstance( mask, list):
            mask = np.array(mask).astype(bool)
    I = FormatInMemory(image=image, mask=mask)
    reader = MemReader([I])
    iset_Data = ImageSetData(reader, None) # , masker)
    iset = ImageSet(iset_Data)
    iset.set_beam(beam)
    iset.set_detector(detector)
    explist = ExperimentListFactory.from_imageset_and_crystal(iset, None)
    return explist


def parse_panel_input_string(panel_str):
    panel_range_str = panel_str.split(",")
    all_pids =[]
    for s in panel_range_str:
        if s == "":
            continue
        if "-" in s:
            if s.count("-") > 1:
                raise ValueError("panel str %s is not formatted correctly, too many '-' between commas" % panel_str )
            pid1, pid2 = map(int, s.split("-"))
            if not pid1 < pid2:
                raise ValueError("panel str %s is not formatted correctly, in 'x-y' x should be less than y" % panel_str )
            pids = [pid for pid in range(pid1, pid2+1)]
        else:
            pids = [int(s)]
        all_pids += pids
    if len(all_pids) != len(set(all_pids)):
        MAIN_LOGGER.warning("WARNING: duplicate pids found in panel input string %s, please double check!" % panel_str)
    return list(set(all_pids))


def get_ncells_mask_from_string(mask_string):
    assert mask_string in ["000", "101", "110", "011", "111"]
    s0, s1, s2 = mask_string
    mask = int(s0), int(s1), int(s2)
    return mask


def load_panel_group_file(panel_group_file):
    lines = open(panel_group_file, 'r').readlines()
    groups = {}
    for l in lines:
        try:
            panel, group = map(int, l.strip().split())
        except ValueError:
            continue
        groups[panel] = group
    return groups


def spots_from_pandas(pandas_frame, mtz_file=None, mtz_col=None,
                      oversample_override=None,
                      Ncells_abc_override=None,
                      pink_stride_override=None,
                      cuda=False, device_Id=0, time_panels=False,
                      d_max=999, d_min=1.5, defaultF=1000.,
                      njobs=1,
                      output_img=None, omp=False, norm_by_spectrum=False,
                      symbol_override=None, quiet=False, reset_Bmatrix=False):
    if time_panels and quiet:
        MAIN_LOGGER.debug("NOTE: quiet=True will suppress panel simulation timing print output")
    from joblib import Parallel, delayed
    from simtbx.nanoBragg.utils import flexBeam_sim_colors

    df = pandas_frame

    if not quiet:MAIN_LOGGER.info("Loading experiment models")
    expt_name = df.opt_exp_name.values[0]
    El = ExperimentListFactory.from_json_file(expt_name, check_format=False)
    expt = El[0]

    if "detz_shift_mm" in list(df):  # NOTE, this could also be inside expt_name directly
        expt.detector = shift_panelZ(expt.detector, df.detz_shift_mm.values[0])

    if reset_Bmatrix:
        ucell_params = df[["a", "b", "c", "al", "be", "ga"]].values[0]
        ucell_man = manager_from_params(ucell_params)
        expt.crystal.set_B(ucell_man.B_recipspace)
    if not quiet:MAIN_LOGGER.info("Done loading models!")
    if not quiet:MAIN_LOGGER.info("Crystal model:")
    if not quiet:El[0].crystal.show()
    assert len(df) == 1
    Ncells_abc = tuple(map(lambda x: int(round(x)), df.ncells.values[0]))
    if Ncells_abc_override is not None:
        Ncells_abc = Ncells_abc_override
    spot_scale = df.spot_scales.values[0]
    beamsize_mm = df.beamsize_mm.values[0]
    total_flux = df.total_flux.values[0]
    oversample = df.oversample.values[0]
    if oversample_override is not None:
        oversample = oversample_override

    # get the optimized spectra
    if "spectrum_filename" in list(df) and df.spectrum_filename.values[0] is not None:
        spectrum_file = df.spectrum_filename.values[0]
        pink_stride = df.spectrum_stride.values[0]
        if norm_by_spectrum:
            nspec = load_spectra_file(spectrum_file)[0].shape[0]
            spot_scale = spot_scale / nspec
        if pink_stride_override is not None:
            pink_stride = pink_stride_override
        fluxes, energies = load_spectra_file(spectrum_file, total_flux=total_flux,
                                             pinkstride=pink_stride)
    else:
        assert total_flux is not None
        fluxes = np.array([total_flux])
        energies = np.array([ENERGY_CONV/expt.beam.get_wavelength()])
        if not quiet: MAIN_LOGGER.info("Running MONO sim")
        nspec = 1
    lam0 = df.lam0.values[0]
    lam1 = df.lam1.values[0]
    if lam0 == -1:
        lam0 = 0
    if lam1 == -1:
        lam1 = 1
    wavelens = ENERGY_CONV / energies
    wavelens = lam0 + lam1*wavelens
    energies = ENERGY_CONV / wavelens

    if mtz_file is not None:
        assert mtz_col is not None
        Famp = open_mtz(mtz_file, mtz_col)
    else:
        Famp = make_miller_array_from_crystal(expt.crystal, dmin=d_min, dmax=d_max, defaultF=defaultF, symbol=symbol_override)

    crystal = expt.crystal
    # TODO opt exp should already include the optimized Amatrix, so no need to set it twice
    crystal.set_A(df.Amats.values[0])

    panel_list = list(range(len(expt.detector)))
    pids_per_job = np.array_split(panel_list, njobs)

    def main(pids):
        results = flexBeam_sim_colors(CRYSTAL=expt.crystal, DETECTOR=expt.detector, BEAM=expt.beam, Famp=Famp,
                                      fluxes=fluxes, energies=energies, beamsize_mm=beamsize_mm,
                                      Ncells_abc=Ncells_abc, spot_scale_override=spot_scale,
                                      cuda=cuda, device_Id=device_Id, oversample=oversample, time_panels=time_panels and not quiet,
                                      pids=pids, omp=omp, show_params=not quiet)
        return results

    results = Parallel(n_jobs=njobs)(delayed(main)(pids_per_job[jid]) for jid in range(njobs))
    results = [result for job_results in results for result in job_results]

    if output_img is not None:
        save_model_to_image(expt, results, output_img, save_experiment_data=False)

    pids, imgs = zip(*results)
    order = np.argsort(pids)
    results = np.array([imgs[i] for i in order])

    return results


def roi_spots_from_pandas(pandas_frame,  rois_per_panel, mtz_file=None, mtz_col=None,
                      oversample_override=None,
                      Ncells_abc_override=None,
                      pink_stride_override=None,
                      spectrum_override=None,
                      cuda=False, device_Id=0, time_panels=False,
                      d_max=999, d_min=1.5, defaultF=1e3,
                      omp=False,
                      norm_by_spectrum=False,
                      symbol_override=None, quiet=False, reset_Bmatrix=False, nopolar=False,
                      force_no_detector_thickness=False, printout_pix=None, norm_by_nsource=False):
    if time_panels and quiet:
        MAIN_LOGGER.info("NOTE: quiet=True will suppress panel simulation timing print output")
    from simtbx.nanoBragg.utils import flexBeam_sim_colors

    df = pandas_frame

    if not quiet: MAIN_LOGGER.info("Loading experiment models")
    expt_name = df.opt_exp_name.values[0]
    El = ExperimentListFactory.from_json_file(expt_name, check_format=False)
    expt = El[0]
    if "detz_shift_mm" in list(df):  # NOTE, this could also be inside expt_name directly
        expt.detector = shift_panelZ(expt.detector, df.detz_shift_mm.values[0])

    if force_no_detector_thickness:
        expt.detector = strip_thickness_from_detector(expt.detector)
    if reset_Bmatrix:
        ucell_params = df[["a", "b", "c", "al", "be", "ga"]].values[0]
        ucell_man = manager_from_params(ucell_params)
        expt.crystal.set_B(ucell_man.B_recipspace)
    if not quiet:MAIN_LOGGER.info("Done loading models!")
    if not quiet:MAIN_LOGGER.info("Crystal model:")
    if not quiet:expt.crystal.show()
    assert len(df) == 1
    Ncells_abc = df.ncells.values[0]  #tuple(map(lambda x: int(round(x)), df.ncells.values[0]))
    if Ncells_abc_override is not None:
        Ncells_abc = Ncells_abc_override
    spot_scale = df.spot_scales.values[0]
    beamsize_mm = df.beamsize_mm.values[0]
    total_flux = df.total_flux.values[0]
    oversample = df.oversample.values[0]
    if oversample_override is not None:
        oversample = oversample_override

    # get the optimized spectra
    if spectrum_override is None:
        if "spectrum_filename" in list(df) and df.spectrum_filename.values[0] is not None:
            spectrum_file = df.spectrum_filename.values[0]
            pink_stride = df.spectrum_stride.values[0]
            if norm_by_spectrum:
                nspec = load_spectra_file(spectrum_file)[0].shape[0]
                spot_scale = spot_scale / nspec
            if pink_stride_override is not None:
                pink_stride = pink_stride_override
            fluxes, energies = load_spectra_file(spectrum_file, total_flux=total_flux,
                                                 pinkstride=pink_stride)
        else:
            fluxes = np.array([total_flux])
            energies = np.array([ENERGY_CONV/expt.beam.get_wavelength()])
            if not quiet: MAIN_LOGGER.info("Running MONO sim")

    else:
        wavelens, fluxes = map(np.array, zip(*spectrum_override))
        energies = ENERGY_CONV / wavelens

    lam0 = df.lam0.values[0]
    lam1 = df.lam1.values[0]
    if lam0 == -1:
        lam0 = 0
    if lam1 == -1:
        lam1 = 1
    wavelens = ENERGY_CONV / energies
    wavelens = lam0 + lam1*wavelens
    energies = ENERGY_CONV / wavelens

    if mtz_file is not None:
        assert mtz_col is not None
        Famp = open_mtz(mtz_file, mtz_col)
    else:
        Famp = make_miller_array_from_crystal(expt.crystal, dmin=d_min, dmax=d_max, defaultF=defaultF, symbol=symbol_override)

    results = flexBeam_sim_colors(CRYSTAL=expt.crystal, DETECTOR=expt.detector, BEAM=expt.beam, Famp=Famp,
                                  fluxes=fluxes, energies=energies, beamsize_mm=beamsize_mm,
                                  Ncells_abc=Ncells_abc, spot_scale_override=spot_scale,
                                  cuda=cuda, device_Id=device_Id, oversample=oversample, time_panels=time_panels and not quiet,
                                  pids=list(rois_per_panel.keys()),
                                  rois_perpanel=rois_per_panel,
                                  omp=omp, show_params=not quiet, nopolar=nopolar,
                                  printout_pix=printout_pix)
    if norm_by_nsource:
        return {pid: image/len(energies) for pid, image in results}
    else:
        return {pid: image for pid, image in results}


def spots_from_pandas_and_experiment(expt, pandas_pickle, mtz_file=None, mtz_col=None,
                                     spectrum_file=None, total_flux=1e12, pink_stride=1,
                                     beamsize_mm=0.001, oversample=0, d_max=999, d_min=1.5, defaultF=1e3,
                                     cuda=False, device_Id=0, time_panels=True,
                                     output_img=None, njobs=1, save_expt_data=False,
                                     as_numpy_array=False, omp=False):
    import pandas
    from joblib import Parallel, delayed
    from simtbx.nanoBragg.utils import flexBeam_sim_colors

    MAIN_LOGGER.info("Loading experiment")
    if isinstance(expt, str):
        El = ExperimentListFactory.from_json_file(expt, check_format=save_expt_data)
        expt = El[0]
    else:
        assert "Experiment" in str(type(expt))
    MAIN_LOGGER.info("Done loading!")
    df = pandas.read_pickle(pandas_pickle)
    assert len(df) == 1
    Ncells_abc = tuple(map(lambda x : int(round(x)), df.ncells.values[0]))
    spot_scale = df.spot_scales.values[0]
    if mtz_file is not None:
        assert mtz_col is not None
        Famp = open_mtz(mtz_file, mtz_col)
    else:
        Famp = make_miller_array_from_crystal(expt.crystal, dmin=d_min, dmax=d_max, defaultF=defaultF)

    # get the optimized spectra
    if spectrum_file is not None:
        fluxes, energies = load_spectra_file(spectrum_file, total_flux=total_flux, pinkstride=pink_stride)
    else:
        fluxes = np.array([total_flux])
        energies = np.array([ENERGY_CONV/expt.beam.get_wavelength()])
    lam0 = df.lam0.values[0]
    lam1 = df.lam1.values[0]
    if lam0 == -1:
        lam0 = 0
    if lam1 == -1:
        lam1 = 1
    wavelens = ENERGY_CONV / energies
    wavelens = lam0 + lam1*wavelens
    energies = ENERGY_CONV / wavelens

    panel_list = list(range(len(expt.detector)))
    pids_per_job = np.array_split(panel_list, njobs)

    def main(pids):
        results = flexBeam_sim_colors(CRYSTAL=expt.crystal, DETECTOR=expt.detector, BEAM=expt.beam, Famp=Famp,
                                      fluxes=fluxes, energies=energies, beamsize_mm=beamsize_mm,
                                      Ncells_abc=Ncells_abc, spot_scale_override=spot_scale,
                                      cuda=cuda, device_Id=device_Id, oversample=oversample, time_panels=time_panels,
                                      pids=pids, omp=omp)
        return results

    results = Parallel(n_jobs=njobs)(delayed(main)(pids_per_job[jid]) for jid in range(njobs))
    results = [result for job_results in results for result in job_results]

    if output_img is not None:
        save_model_to_image(expt, results, output_img, save_experiment_data=save_expt_data)
    if as_numpy_array:
        pids, imgs = zip(*results)
        order = np.argsort(pids)
        assert np.all(order == np.arange(len(pids)))
        results = np.array([imgs[i] for i in order])

    return results


def make_miller_array_from_crystal(Crystal, dmin, dmax, defaultF=1000, symbol=None):
    if symbol is None:
        symbol = Crystal.get_space_group().info().type().lookup_symbol()
    Famp = make_miller_array(
        symbol=symbol,
        unit_cell=Crystal.get_unit_cell(), d_min=dmin, d_max=dmax, defaultF=defaultF)
    return Famp


def save_model_to_image(expt, model_results, output_img_file, save_experiment_data=False):
    ordered_img = {}
    npanels = len(expt.detector)
    n_img = 1
    if save_experiment_data:
        n_img = 2
    for pid, img in model_results:
        ordered_img[pid] = img
    ordered_img = np.array([ordered_img[pid] for pid in range(npanels)])
    panelX, panelY = expt.detector[0].get_image_size()
    from simtbx.nanoBragg.utils import H5AttributeGeomWriter

    with H5AttributeGeomWriter(filename=output_img_file, image_shape=(npanels, panelY, panelX), num_images=n_img,
                               beam=expt.beam, detector=expt.detector) as writer:
        writer.add_image(ordered_img)
        if save_experiment_data:
            exp_data = image_data_from_expt(expt)
            writer.add_image(exp_data)
    MAIN_LOGGER.info("Wrote model to image %s" % output_img_file)


def index_refls(refls, exper, tolerance=0.333):
    from dials.algorithms.indexing import assign_indices
    refls['id'] = flex.int(len(refls), -1)
    refls['imageset_id'] = flex.int(len(refls), 0)
    El = ExperimentList()
    El.append(exper)
    refls.centroid_px_to_mm(El)
    refls.map_centroids_to_reciprocal_space(El)
    idx_assign = assign_indices.AssignIndicesGlobal(tolerance=tolerance)
    idx_assign(refls, El)
    return refls


def indexed_from_model(strong_refls, model_images, expt, thresh=1, tolerance=0.333, Qdist_cutoff=0.003):
    model_refls = refls_from_sims(model_images, expt.detector, expt.beam, thresh=thresh)
    model_refls['id'] = flex.int(len(model_refls), 0)
    model_refls['imageset_id'] = flex.int(len(model_refls), 0)

    El = ExperimentList()
    El.append(expt)

    model_refls = index_refls(model_refls, expt, tolerance=tolerance)

    #_=refls_to_hkl(model_refls, expt.detector, expt.beam, expt.crystal,
    #                                update_table=True)
    MAIN_LOGGER.info("Indexed %d / %d spots from the model using the nominal wavelength"
          % (sum(model_refls['id'] == 0), len(model_refls)))
    if strong_refls is None:
        model_refls['xyzcal.mm'] = model_refls['xyzobs.mm.value']
        model_refls['xyzcal.px'] = model_refls['xyzobs.px.value']
        return model_refls

    strong_refls.centroid_px_to_mm(El)
    strong_refls.map_centroids_to_reciprocal_space(El)

    pids = set(strong_refls['panel'])
    Rindexed = None
    for pid in pids:
        Rmodel = model_refls.select(model_refls['panel'] == pid)
        Rstrong = strong_refls.select(strong_refls['panel'] == pid)
        if len(Rmodel) == 0:
            continue
        Qxyz_model = np.array(Rmodel['rlp'])
        Qxyz_data = np.array(Rstrong['rlp'])
        tree = cKDTree(Qxyz_model)
        dists, nearest = tree.query(Qxyz_data, k=1)   # find the nearest model spot to each data spot

        miller_index = [Rmodel['miller_index'][n] for n in nearest]
        xyzcal_px = [Rmodel['xyzobs.px.value'][n] for n in nearest]
        xyzcal_mm = [Rmodel['xyzobs.mm.value'][n] for n in nearest]

        Rstrong['miller_index'] = flex.miller_index(miller_index)
        Rstrong['xyzcal.px'] = flex.vec3_double(xyzcal_px)
        Rstrong['xyzcal.mm'] = flex.vec3_double(xyzcal_mm)
        within_cutoff = flex.bool(dists <= Qdist_cutoff)
        Rstrong = Rstrong.select(within_cutoff)

        if Rindexed is None:
            Rindexed = Rstrong
        else:
            Rindexed.extend(Rstrong)

    return Rindexed


def remove_multiple_indexed(R):
    """R: dials reflection table instance"""
    import pandas
    from dials.array_family import flex
    obs = R['xyzobs.px.value']
    cal = R['xyzcal.px']
    dists = (obs-cal).norms()
    h, k, l = zip(*R['miller_index'])
    df = pandas.DataFrame({"dist": dists, "h": h, "k": k, "l": l})
    gb = df.groupby(["h", "k", "l"])
    hkls = set(R["miller_index"])
    selected_rows = []
    for hkl in hkls:
        g = gb.get_group(hkl)
        MAIN_LOGGER.debug( hkl, len(g))
        if len(g) > 1:
            row = g.iloc[g.dist.argmin()].name
        else:
            row = g.iloc[0].name
        selected_rows.append(row)

    sel = flex.bool([i in selected_rows for i in range(len(R))])
    R = R.select(sel)
    return R


def load_spectra_from_dataframe(df):
    total_flux = df.total_flux.values[0]
    spectrum_file = df.spectrum_filename.values[0]
    pink_stride = df.spectrum_stride.values[0]
    spec = load_spectra_file(spectrum_file, total_flux=total_flux,
                            pinkstride=pink_stride, as_spectrum=True)
    return spec


def save_numpy_mask_as_flex(numpymask, outfile):
    import pickle
    from dials.array_family import flex
    flexmask = tuple((flex.bool(m) for m in numpymask))
    with open(outfile, "wb") as f:
        pickle.dump(flexmask, f)


def load_dials_flex_mask(path_to_mask, as_numpy=True):
    import pickle
    mask = pickle.load(open(path_to_mask, "rb"))
    if as_numpy:
        mask = np.array([m.as_numpy_array() for m in mask])
    return mask


def refls_to_q(refls, detector, beam, update_table=False):

    orig_vecs = {}
    fs_vecs = {}
    ss_vecs = {}
    u_pids = set(refls['panel'])
    for pid in u_pids:
        orig_vecs[pid] = np.array(detector[pid].get_origin())
        fs_vecs[pid] = np.array(detector[pid].get_fast_axis())
        ss_vecs[pid] = np.array(detector[pid].get_slow_axis())

    s1_vecs = []
    q_vecs = []
    panels = refls["panel"]
    n_refls = len(refls)
    for i_r in range(n_refls):
        r = refls[i_r]
        pid = r['panel']
        i_fs, i_ss, _ = r['xyzobs.px.value']
        panel = detector[pid]
        orig = orig_vecs[pid] #panel.get_origin()
        fs = fs_vecs[pid] #panel.get_fast_axis()
        ss = ss_vecs[pid] #panel.get_slow_axis()

        fs_pixsize, ss_pixsize = panel.get_pixel_size()
        s1 = orig + i_fs*fs*fs_pixsize + i_ss*ss*ss_pixsize  # scattering vector
        s1 = s1 / np.linalg.norm(s1) / beam.get_wavelength()
        s1_vecs.append(s1)
        q_vecs.append(s1-beam.get_s0())

    if update_table:
        refls['s1'] = flex.vec3_double(tuple(map(tuple, s1_vecs)))
        refls['rlp'] = flex.vec3_double(tuple(map(tuple, q_vecs)))

    return np.vstack(q_vecs)


def strong_spot_mask(refl_tbl, detector, dilate=None):
    """

    :param refl_tbl:
    :param detector:
    :return:  strong spot mask same shape as detector (multi panel)
    """
    # dxtbx detector size is (fast scan x slow scan)
    nfast, nslow = detector[0].get_image_size()
    n_panels = len(detector)
    is_strong_pixel = np.zeros((n_panels, nslow, nfast), bool)
    panel_ids = refl_tbl["panel"]
    # mask the strong spots so they dont go into tilt plane fits
    for i_refl in range(len(refl_tbl)): # (x1, x2, y1, y2, _, _) in enumerate(refl_tbl['bbox']):
        i_panel = panel_ids[i_refl]
        sb = refl_tbl[i_refl]["shoebox"]
        x1, x2, y1, y2, _, _ = sb.bbox
        is_strong = sb.mask.as_numpy_array()[0] == 5
        bb_ss = slice(y1, y2, 1)
        bb_fs = slice(x1, x2, 1)
        is_strong_pixel[i_panel, bb_ss, bb_fs] = is_strong

    if dilate is not None:
        is_strong_pixel = np.array([binary_dilation(m, iterations=dilate) for m in is_strong_pixel])

    return is_strong_pixel


def modify_crystal(anglesXYZ, ucell_params, starting_crystal):
    from scitbx.matrix import col
    from copy import deepcopy
    x = col((-1, 0, 0))
    y = col((0, -1, 0))
    z = col((0, 0, -1))
    RX = x.axis_and_angle_as_r3_rotation_matrix(anglesXYZ[0], deg=False)
    RY = y.axis_and_angle_as_r3_rotation_matrix(anglesXYZ[1], deg=False)
    RZ = z.axis_and_angle_as_r3_rotation_matrix(anglesXYZ[2], deg=False)
    M = RX * RY * RZ
    q = M.r3_rotation_matrix_as_unit_quaternion()
    rot_ang, rot_ax = q.unit_quaternion_as_axis_and_angle(deg=True)

    UCELL_MAN = manager_from_params(ucell_params)
    B = UCELL_MAN.B_recipspace

    C = deepcopy(starting_crystal)
    C.set_B(B)
    if rot_ang > 0:
        C.rotate_around_origin(rot_ax, rot_ang)
    return C


def extract_parmaeters(experiment, refls, use_back_computed_wavelen=False):
    from xfel.mono_simulation import max_like
    from dials.array_family import flex
    from scitbx.matrix import sqr, col
    from dials.algorithms.refinement.prediction.managed_predictors import ExperimentsPredictorFactory
    from dxtbx.model import Experiment
    from copy import deepcopy
    from dxtbx_model_ext import SimplePxMmStrategy

    crystal = experiment.crystal
    Amat = crystal.get_A()
    beam = experiment.beam
    s0_nominal = col(beam.get_s0())
    det = experiment.detector
    det = strip_thickness_from_detector(det)

    refls["local_id"] = flex.int_range(len(refls))

    all_new_delpsi = []
    d_i = []

    has_wave = False
    if "wavelength_Ang" in list(refls.keys()):
        has_wave = True

    waves2 =[]
    for i_ref in range(len(refls)):
        if has_wave:
            wave = refls[i_ref]["wavelength_Ang"]
        else:
            wave = beam.get_wavelength()
        pid = refls[i_ref]['panel']
        panel = det[pid]
        hkl = col(refls[i_ref]['miller_index'])

        if set(hkl)=={0}:
            continue

        # is this necessary ?
        panel.set_px_mm_strategy(SimplePxMmStrategy())
        panel.set_mu(0)
        panel.set_thickness(0)

        res = 1 / np.linalg.norm(refls[i_ref]['rlp'])
        d_i.append(res)
        xcal, ycal, _ = refls[i_ref]['xyzcal.mm']

        if use_back_computed_wavelen:
            cent = panel.get_beam_centre_lab(beam.get_unit_s0())
            x, y, _ = refls[i_ref]['xyzobs.mm.value']
            xyz = panel.get_lab_coord((x, y))
            rad_dist = np.sqrt(np.sum((np.array(xyz) - np.array(cent)) ** 2))
            det_dist = panel.get_distance()
            theta = 0.5 * np.arctan(rad_dist / det_dist)
            back_computed_wavelen = 2 * res * np.sin(theta)
            waves2.append(back_computed_wavelen)
            wave = back_computed_wavelen
            # re-predict the spot with new wavelength
            El = ExperimentList()
            E = deepcopy(experiment)
            E.beam.set_wavelength(back_computed_wavelen)
            E.detector = det
            El.append(E)
            ref_predictor = ExperimentsPredictorFactory.from_experiments(El, force_stills=True)
            this_refl = refls.select(refls["local_id"] == i_ref)
            assert len(this_refl) == 1
            this_refl['id'] = flex.int(1, 0)  # reset id so ref_predictor can match to El
            this_refl = ref_predictor(this_refl)  # this seems to modify inplace, so = statement not necessary?
            xcal, ycal, _ = this_refl[0]['xyzcal.mm']  # this is the new position of the prediction


        # get the s1 vector with nominal shot wavelength
        xyzcal = col(panel.get_lab_coord((xcal, ycal)))
        s1_spot = xyzcal.normalize() / wave
        s0_recalc = s0_nominal.normalize() / wave

        # get the rlp
        qhkl = sqr(Amat) * hkl
        q0 = qhkl.normalize()

        # r-vector with shot wavelength
        qcal = s1_spot - s0_recalc

        # compute using q1,q0 projection
        e1 = (q0.cross(s0_recalc)).normalize()
        q1 = q0.cross(e1)
        delpsi = -np.arctan2(qcal.dot(q1), qcal.dot(q0))
        all_new_delpsi.append(delpsi)

    Deff = 2000
    eta_rad = 0.001
    out = max_like.minimizer(flex.double(d_i), flex.double(all_new_delpsi), eta_rad, Deff)
    D = 2 / out.x[0]
    eta = out.x[1] * 180 / np.pi
    return D, eta, waves2


def parse_reso_string(s):
    """
    :param s:  a strong formated as %f-%f,%f-%f,%f-%f etc.
    :return: two floats as a tuple
    """
    vals =[]
    try:
        for subs in s.strip().split(","):
            a, b = map(float, subs.strip().split("-"))
            assert a < b
            vals.append((a,b))
    except Exception as error:
        MAIN_LOGGER.error("Failed to parse string!", error)
        raise ValueError("Wrong string format, see error above")
    return vals


def integrate_subimg(subimg, is_bragg_peak, tilt_abc, sigma_rdout_in_photon, covar, is_hotpix=None):
    Y,X = np.indices(subimg.shape)
    sel = is_bragg_peak

    if is_hotpix is not None:
        sel = np.logical_and(is_bragg_peak, ~is_hotpix)
    p = X[sel]  # fast scan coords
    q = Y[sel]  # slow scan coords
    rho_peak = subimg[sel]  # pixel values

    t1, t2, t3 = tilt_abc
    Isum = np.sum(rho_peak - t1*p - t2*q - t3)  # summed spot intensity

    var_rho_peak = sigma_rdout_in_photon**2 + rho_peak  # include readout noise in the variance
    Ns = len(rho_peak)  # number of integrated peak pixels

    # variance propagated from tilt plane constants
    var_a_term = covar[0, 0] * ((np.sum(p)) ** 2)
    var_b_term = covar[1, 1] * ((np.sum(q)) ** 2)
    var_c_term = covar[2, 2] * (Ns ** 2)

    # total variance of the spot
    var_Isum = np.sum(var_rho_peak) + var_a_term + var_b_term + var_c_term

    return Isum, var_Isum


def strong_spot_mask2(refl_tbl, detector):
    """
    :param refl_tbl:
    :param detector:
    :return:  strong spot mask same shape as detector (multi panel)
    """
    # dxtbx detector size is (fast scan x slow scan)
    nfast, nslow = detector[0].get_image_size()
    n_panels = len(detector)
    is_strong_pixel = np.zeros((n_panels, nslow, nfast), bool)
    panel_ids = refl_tbl["panel"]
    # mask the strong spots so they dont go into tilt plane fits
    for i_refl, (x1, x2, y1, y2, _, _) in enumerate(refl_tbl['bbox']):
        i_panel = panel_ids[i_refl]
        bb_ss = slice(y1, y2, 1)
        bb_fs = slice(x1, x2, 1)
        is_model = refl_tbl[i_refl]['shoebox'].mask.as_numpy_array()[0] == 5
        is_strong_pixel[i_panel,bb_ss, bb_fs] = is_model
    return is_strong_pixel

def q_to_hkl(q_vecs, crystal):
    from scitbx.matrix import sqr
    Ai = sqr(crystal.get_A()).inverse()
    Ai = Ai.as_numpy_array()
    HKL = np.dot( Ai, q_vecs.T)
    HKLi = map( lambda h: np.ceil(h-0.5).astype(int), HKL)
    return np.vstack(HKL).T, np.vstack(HKLi).T


def refls_to_pixelmm(refls, detector):
    """
    returns the mm position of the spot in pixels
    referenced to the panel origin, which I think is the
    center of the first pixel in memory for the panel
    :param reflection table
    :param dxtbx detecotr
    :return: np.array Nreflsx2 for fast,slow mm coord referenced to panel origin
        in the plane of the panel
    """
    ij_mm = np.zeros( (len(refls),2))
    for i_r,r in enumerate(refls):
        panel = detector[r['panel']]
        i,j,_ = r['xyzobs.px.value']
        ij_mm[i_r] = panel.pixel_to_millimeter( (i,j) )
    return ij_mm


def refls_to_hkl(refls, detector, beam, crystal,
                 update_table=False, returnQ=False):
    """
    convert pixel panel reflections to miller index data
    :param refls:  reflecton table for a panel or a tuple of (x,y)
    :param detector:  dxtbx detector model
    :param beam:  dxtbx beam model
    :param crystal: dxtbx crystal model
    :param update_table: whether to update the refltable
    :param returnQ: whether to return intermediately computed q vectors
    :return: if as_numpy two Nx3 numpy arrays are returned
        (one for fractional and one for whole HKL)
        else dictionary of hkl_i (nearest) and hkl (fractional)
    """
    from scitbx.matrix import sqr
    if 'rlp' not in list(refls.keys()):
        q_vecs = refls_to_q(refls, detector, beam, update_table=update_table)
    else:
        q_vecs = np.vstack([refls[i_r]['rlp'] for i_r in range(len(refls))])
    Ai = sqr(crystal.get_A()).inverse()
    Ai = Ai.as_numpy_array()
    HKL = np.dot( Ai, q_vecs.T)
    HKLi = np.ceil(HKL-0.5) #map( lambda h: np.ceil(h-0.5).astype(int), HKL)
    if update_table:
        #refls['miller_index'] = flex.miller_index(len(refls),(0,0,0))
        refls['miller_index'] = flex.miller_index(list(map(tuple, HKLi.T.astype(np.int32))))
        #from IPython import embed
        #embed();exit()
        #mil_idx = flex.vec3_int(tuple(map(tuple, np.vstack(HKLi).T)))
        #for i in range(len(refls)):
        #    refls['miller_index'][i] = mil_idx[i]
    if returnQ:
        return np.vstack(HKL).T, np.vstack(HKLi).T, q_vecs
    else:
        return np.vstack(HKL).T, np.vstack(HKLi).T


def refls_to_q(refls, detector, beam, update_table=False):

    orig_vecs = {}
    fs_vecs = {}
    ss_vecs = {}
    u_pids = set(refls['panel'])
    for pid in u_pids:
        orig_vecs[pid] = np.array(detector[pid].get_origin())
        fs_vecs[pid] = np.array(detector[pid].get_fast_axis())
        ss_vecs[pid] = np.array(detector[pid].get_slow_axis())

    s1_vecs = []
    q_vecs = []
    for i_r in range(len(refls)):
        r = refls[i_r]
        pid = r['panel']
        i_fs, i_ss, _ = r['xyzobs.px.value']
        panel = detector[pid]
        orig = orig_vecs[pid] #panel.get_origin()
        fs = fs_vecs[pid] #panel.get_fast_axis()
        ss = ss_vecs[pid] #panel.get_slow_axis()

        fs_pixsize, ss_pixsize = panel.get_pixel_size()
        s1 = orig + i_fs*fs*fs_pixsize + i_ss*ss*ss_pixsize  # scattering vector
        s1 = s1 / np.linalg.norm(s1) / beam.get_wavelength()
        s1_vecs.append(s1)
        q_vecs.append(s1-beam.get_s0())

    if update_table:
        refls['s1'] = flex.vec3_double(tuple(map(tuple,s1_vecs)))
        refls['rlp'] = flex.vec3_double(tuple(map(tuple,q_vecs)))

    return np.vstack(q_vecs)


def get_panels_fasts_slows(expt, pids, rois):
    npan = len(expt.detector)
    nfast, nslow = expt.detector[0].get_image_size()
    MASK = np.zeros((npan, nslow, nfast), bool)
    ROI_ID = np.zeros((npan, nslow, nfast), 'uint16')
    #ROI_ID = NP_ONES((npan, nslow, nfast), 'uint16') * mx
    nspots = len(rois)
    for i_spot in range(nspots):
        x1, x2, y1, y2 = rois[i_spot]
        if x2-x1 == 0 or y2-y1 == 0:
            continue
        pid = pids[i_spot]
        MASK[pid, y1:y2, x1:x2] = True
        ROI_ID[pid, y1:y2, x1:x2] = i_spot
    p,s,f = np.where(MASK)
    roi_id = ROI_ID[p,s,f]
    pan_fast_slow = np.ascontiguousarray((np.vstack([p,f,s]).T).ravel())
    pan_fast_slow = flex.size_t(pan_fast_slow)
    return pan_fast_slow, roi_id

def refls_by_panelname(refls):
    Nrefl = len(refls)
    panel_names = np.array([refls["panel"][i] for i in range( Nrefl)])
    uniq_names = np.unique(panel_names)
    refls_by_panel = {name: refls.select(flex.bool(panel_names == name))
                      for name in uniq_names }
    return refls_by_panel


def f_double_prime(energies, a,b,c,d, deriv=None):
    """
    generate a 4 parameter f double prime curve based on the sigmoid function
    :param energies:
    :param a: offset
    :param b: amplitude
    :param c: center
    :param d: slope
    :param deriv: string flag, c,d or None,
        eg 'c' returns deriv of f_dbl_prime w.r.t. c, None returns function f_dbl_prime
    :return: f double prime as a function of energy
    """
    exp_arg = -d * (energies - c)
    e_term = np.exp(exp_arg)
    if deriv=="c":
        return -b*(1+e_term)**-2 * e_term * d
    elif deriv=="d":
        return -b*(1+e_term)**-2 * e_term * (c-energies)
    else:
        return a + b / (1+e_term)


def f_prime(f_double_prime, S=None, padn=5000):
    """
    generate an f_prime from an f_double_prime curve
    using the kramers kronig relationship
    :param f_double_prime: function or its derivative
    :param S:
    :param padn:
    :return:
    """
    if S is None:
        S = np.sin(np.linspace(-np.pi / 2, np.pi / 2, padn)) * 0.5 + 0.5
    else:
        padn = S.shape[0]
    F = f_double_prime
    Fin = np.hstack((F[0] * S,
                     F,
                     F[-1] * (1 - S)))  # sin padding as used in Sherrell thesis, TODO window function ?
    Ft = fft.fft(Fin)
    iFt = -1 * fft.ifft(1j * np.sign(fft.fftfreq(Ft.shape[0])) * Ft).real
    return iFt[padn:-padn]


def make_gauss(sig, mu):
    return lambda x: 1 / (sig * (2*np.pi)**.5) * np.e ** (-(x-mu)**2/(2 * sig**2))

def make_lrntz(gam, mu):
    return lambda x: gam / (np.pi*((x-mu)**2 + gam**2))

def make_voigt(sig, gam, mu, x):
    #
    #
    from scipy.interpolate import interp1d
    gauss = make_gauss(sig, mu)(x)
    lrntz = make_lrntz(gam, mu)(x)
    voigt = np.convolve(gauss, lrntz)
    voigt /= abs(voigt.max())
    r = np.linspace(x[0], x[-1], len(voigt))
    voigt_intrp = interp1d(r, voigt, kind='linear')(x)
    return voigt, voigt_intrp


def shift_panelZ(D, shift):
    """
    :param D:  dxtbx detector
    :param shift:  shift in mm for origin Z component
    :return: new detector with shift applied to each panel origin
    """
    newD = Detector()
    for pid in range(len(D)):
        panel = D[pid]
        pan_dict = panel.to_dict()
        x,y,z = panel.get_origin()
        pan_dict["origin"] = x,y,z+shift
        pan_dict["fast_axis"] = panel.get_fast_axis()
        pan_dict["slow_axis"] = panel.get_slow_axis()
        newP = Panel.from_dict(pan_dict)
        newD.add_panel(newP)
    return newD

9
def safe_makedirs(name):
    if not os.path.exists(name):
        os.makedirs(name)


def get_rs(name):
    s = re.search("run[0-9]+_shot[0-9]+", name)
    return name[s.start():s.end()]


def read_exp_ref_spec(fname):
    return [l.strip().split() for l in open(fname, 'r').readlines()]


def _rfactor_minimizer_target(k, F1, F2):
    return F1.r1_factor(F2, scale_factor=k[0])


def compute_scale_to_minmize_r_factor(F1, F2, anom=True):
    """
    F1, F2 : two miller arrays (should be of type amplitude)
    scale factor should be applied to F2 such that the R factor is minmized between the two arrays!"""
    if min(F1.data()) < 0:
        F1 = F1.as_amplitude_array()
    if min(F2.data()) < 0:
        F2 = F2.as_amplitude_array()

    indices_common = set(F1.indices()).intersection(F2.indices())
    indices_common = flex.miller_index(list(indices_common))
    F1 = F1.select_indices(indices_common)
    F2_mset = F1.miller_set(indices_common, anom)
    F2_map = {h:d for h,d in zip(F2.indices(), F2.data())}
    F2_data = flex.double([F2_map[h] for h in indices_common])
    F2 = miller.array(F2_mset, F2_data)
    F1=F1.sort(by_value='packed_indices')
    F2=F2.sort(by_value='packed_indices')

    res = minimize(_rfactor_minimizer_target,
                       x0=[1], args=(F1, F2),
                   method='Nelder-Mead')

    assert res.success
    r1_scale = res.x[0]
    MAIN_LOGGER.debug("Optimization successful!, using scale factor=%f" % r1_scale)
    return r1_scale


def get_ave_FF(F, w=500, use_fft=True):
    d = np.array(F.d_spacings().data())
    order = np.argsort(d)
    amps = np.array(F.data())**2
    MA = moving_average
    if use_fft:
        MA = moving_average_fft
    amps_ave = MA(amps[order], w)
    d_ave = MA(d[order], w)
    return d_ave, amps_ave


def moving_average(x, w):  # stack overflow recipe
    return np.convolve(x, np.ones(w), 'valid') / w


def moving_average_fft(x, w):
    fftx = np.fft.fft(x)
    fftw = np.fft.fft(np.ones(w))
    move_ave = np.fft.ifft( fftx*fftw ).real
    return move_ave / w


#def moving_average_fft(x, w):
#    fftx = np.fft.fft(x)
#
#    nx = len(x)
#    fw = np.zeros(nx)
#    start = int(nx/2) - int(w/2)
#    fw[start: start+w] = 1
#    fftw = np.fft.fft(fw)
#    move_ave = np.fft.ifft( fftx[:len(x)]*fftw[:len(x)] ).real
#    return move_ave / w

def smooth(x, beta=10.0, window_size=11):
    # make sure the window size is odd
    if window_size % 2 == 0:
        window_size += 1

    # apply the smoothing function
    s = np.r_[x[window_size-1:0:-1], x, x[-1:-window_size:-1]]
    w = np.kaiser(window_size, beta)
    y = np.convolve( w/w.sum(), s, mode='valid' )

    # remove the extra array length convolve adds
    b = int((window_size-1) / 2.)
    smoothed = y[b:len(y)-b]

    return smoothed


class WilsonUpdater:
    def __init__(self, dspacing_at_i_fcell, nbins=100):
        self.order = np.argsort(dspacing_at_i_fcell)
        self.dsort = dspacing_at_i_fcell[self.order]
        self.bins = [b[0]-1e-6 for b in np.array_split(self.dsort, nbins)] + [self.dsort[-1]+1e-6]
        #assert min(self.bins) < min(dspacing_at_i_fcell)
        #assert max(self.bins) > max(dspacing_at_i_fcell)
        self.unsorted_digs = np.digitize(dspacing_at_i_fcell, self.bins)
        self.sorted_digs = np.digitize(self.dsort, self.bins)
        self.nbins = nbins
        self.dspacing_at_i_fcell = dspacing_at_i_fcell

    def num_in_i_fcells_res_bin(self):
        """returns the number of Fhkl in the res bin belonging to i_fcell"""
        mapper = {}
        for i in range(1, self.nbins+1):
            num_in_bin = np.sum(self.sorted_digs == i)
            mapper[i] = num_in_bin
        mapped_vals = np.array([mapper[i] for i in self.unsorted_digs])
        return mapped_vals

    def get_reference(self, fhkl_current_at_i_fcell):
        datsort = fhkl_current_at_i_fcell[self.order]**2
        #valmeans = []
        mapper = {}
        for i in range(1, self.nbins+1):
            val_mean = np.mean(datsort[self.sorted_digs == i])
            #val_mean = np.median(datsort[self.sorted_digs == i])
            #valmeans.append(val_mean)
            #d_mean = self.bins[i-1].mean()
            #dmeans.append(d_mean)
            mapper[i] = val_mean
        mapped_vals = np.array([mapper[i] for i in self.unsorted_digs])
        return mapped_vals


def show_diffBragg_state(D, debug_pixel_panelfastslow):
    """
    D, diffBragg instance
    debug_pixel_panelfastslow, 3-tuple of ints, panelId, fast coord, slow coord
    """
    # TODO be careful with zero-ing the pixels - is this really what we want to do ?
    # TODO, rather than print the state, pickle the state
    D.show_params()
    MAIN_LOGGER.info("internal spot scale=%f" % D.spot_scale)
    D.raw_pixels*=0
    p, f, s = debug_pixel_panelfastslow
    D.printout_pixel_fastslow = f, s
    D.add_diffBragg_spots((p, f, s))
    D.raw_pixels*=0


def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))


def fit_gauss( peak_pro, xdata, mu=0, sig=1, amp=100):

    guess = amp, mu,  sig

    try:
        fit,success = optimize.curve_fit(gauss,
                                         xdata=xdata, ydata=peak_pro, p0 =guess )
    except RuntimeError:
        return None
    return fit, success, gauss(xdata, *fit)
