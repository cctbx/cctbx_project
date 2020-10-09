from itertools import zip_longest
import math
import pickle
from simtbx.diffBragg.refiners.crystal_systems import OrthorhombicManager, TetragonalManager, MonoclinicManager, HexagonalManager
from scipy.optimize import minimize
from dials.algorithms.image.filter import convolve
from scipy.interpolate import SmoothBivariateSpline
from cctbx.array_family import flex
from cctbx import miller, sgtbx
from cctbx.crystal import symmetry
import numpy as np
from scipy.ndimage.morphology import binary_dilation
from simtbx.nanoBragg.utils import ENERGY_CONV
from dxtbx.model import Detector, Panel
from scipy import ndimage
import pylab as plt
from simtbx.diffBragg.sim_data import SimData
from simtbx.diffBragg.nanoBragg_beam import nanoBragg_beam
from simtbx.diffBragg.nanoBragg_crystal import nanoBragg_crystal
from dxtbx.imageset import MemReader #, MemMasker
from dxtbx.imageset import ImageSet, ImageSetData
from dxtbx.model.experiment_list import ExperimentListFactory


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


def label_background_pixels(roi_img, thresh=3.5, iterations=1):
    """
    iteratively determine background pixels in a subimg
    """
    img_shape = roi_img.shape
    img1 = roi_img.copy().ravel()   # 1-D version
    background_pixels = None
    while iterations > 0:
        if background_pixels is None:
            outliers = is_outlier(img1, thresh)
            m = np.median(img1[~outliers])
            out_and_high = np.logical_and(outliers, img1 > m)
            background_pixels = ~out_and_high
        else:
            where_bg = np.where(background_pixels)[0]
            outliers = is_outlier(img1[background_pixels], thresh)
            m = np.median(img1[background_pixels][~outliers])
            out_and_high = np.logical_and(outliers, img1[background_pixels] > m)
            background_pixels[where_bg[out_and_high]] = False
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
        m = np.median(img)
        s = np.std(img)
        vmax=m+2*s
        vmin=m-2*s
        plt.imshow(img, vmax=vmax, vmin=vmin)
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
        plt.imshow(img, vmin=vmin, vmax=vmax)
        plt.gca().add_collection(patch_coll)
        plt.show()

    spot_roi = np.array(spot_roi).astype(int)
    tilt_abc = np.array(tilt_abc)
    spot_roi = spot_roi[successes]
    tilt_abc = tilt_abc[successes]
    return spot_roi, tilt_abc


def perturb_miller_array(F, factor, perturb_log_vals=True):
    if not F.is_xray_amplitude_array():
        F = F.amplitudes()
    Fdat = np.array( F.data())

    if perturb_log_vals:
        Fdat = np.log(Fdat)

    try:
        Fdat = np.random.uniform(Fdat-factor, Fdat+factor)
    except OverflowError:
        print ('You suck')
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
            print("alignment matrix", cb_op_align)

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
    #for y in range(sy, -(sy+1), -1):
      for x in range(-sx, sx+1):
        #for x in range(-sx, sx+1):
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


def get_roi_background_and_selection_flags(refls, imgs, shoebox_sz=10, reject_edge_reflections=False,
                                   reject_roi_with_hotpix=True, background_mask=None, hotpix_mask=None,
                                   bg_thresh=3.5, set_negative_bg_to_zero=False,
                                   pad_for_background_estimation=None):

    npan, sdim, fdim = imgs.shape

    if hotpix_mask is not None:
        assert hotpix_mask.shape == imgs.shape

    if background_mask is not None:
        assert background_mask.shape == imgs.shape

    rois, is_on_edge = get_roi_from_spot(refls, fdim, sdim, shoebox_sz=shoebox_sz)
    tilt_abc = []
    kept_rois = []
    panel_ids = []
    selection_flags = []
    num_roi_negative_bg = 0
    for i_roi, roi in enumerate(rois):
        i1, i2, j1, j2 = roi
        is_selected = True
        if is_on_edge[i_roi] and reject_edge_reflections:
            print("is on edge")
            is_selected = False
        pid = refls[i_roi]['panel']

        if hotpix_mask is not None:
            is_hotpix = hotpix_mask[pid, j1:j2, i1:i2]
            num_hotpix = is_hotpix.sum()
            if num_hotpix > 0 and reject_roi_with_hotpix:
                print("has hot pixel")
                is_selected = False

        dimY, dimX = imgs[pid].shape
        if pad_for_background_estimation is not None:
            j1 = max(0, j1-pad_for_background_estimation)
            i1 = max(0, i1-pad_for_background_estimation)
            j2 = min(dimY, j2+pad_for_background_estimation)
            i2 = min(dimX, i2+pad_for_background_estimation)

        shoebox = imgs[pid, j1:j2, i1:i2]

        if background_mask is not None:
            is_background = background_mask[pid, j1:j2, i1:i2]
        else:
            is_background = label_background_pixels(shoebox,thresh=bg_thresh, iterations=2)

        bg_pixels = shoebox[is_background]
        bg_signal = np.median(bg_pixels)

        if bg_signal < 0:
            print("neg")
            num_roi_negative_bg += 1
            if set_negative_bg_to_zero:
                bg_signal = 0
            else:
                print("background is negative")
                is_selected = False
        tilt_abc.append((0, 0, bg_signal))
        kept_rois.append(roi)
        panel_ids.append(pid)
        selection_flags.append(is_selected)

    print("Number of ROI with negative BGs: %d / %d" % (num_roi_negative_bg, len(rois)))

    return kept_rois, panel_ids, tilt_abc, selection_flags


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
                                   mtz_column=None, default_F=0, dmin=1.5, dmax=30, spectra_file=None, spectra_stride=1):

    if params is not None:
        oversample = params.simulator.oversample
        device_id = params.simulator.device_id
        init_scale = params.simulator.init_scale
        total_flux = params.simulator.total_flux

        ncells_abc = params.simulator.crystal.ncells_abc
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
    crystal = nanoBragg_crystal()
    crystal.isotropic_ncells = has_isotropic_ncells
    crystal.dxtbx_crystal = expt.crystal
    crystal.thick_mm = 0.1  # hard code a thickness, will be over-written by the scale
    crystal.Ncells_abc = tuple(ncells_abc)  #params.simulator.init_ncells_abc
    crystal.n_mos_domains = num_mosaicity_samples
    crystal.mos_spread_deg = mosaicity
    if mtz_name is None:
        miller_data = make_miller_array(
            symbol=expt.crystal.get_space_group().info().type().lookup_symbol(),
            unit_cell=expt.crystal.get_unit_cell(), d_min=dmin, d_max=dmax)
    else:
        miller_data = open_mtz(mtz_name, mtz_column)
    crystal.miller_array = miller_data
    SIM.crystal = crystal

    # create a nanoBragg beam
    beam = nanoBragg_beam()
    beam.size_mm = params.simulator.beam.size_mm
    beam.unit_s0 = expt.beam.get_unit_s0()
    if spectra_file is not None:
        init_spectrum = load_spectra_file(spectra_file, total_flux, spectra_stride, as_spectrum=True)
    else:
        init_spectrum = [(expt.beam.get_wavelength(), total_flux)]
    beam.spectrum = init_spectrum
    SIM.beam = beam

    SIM.panel_id = 0
    SIM.instantiate_diffBragg(oversample=oversample, device_Id=device_id, default_F=default_F)
    if init_scale is not None:
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
        print("Opening mtz file %s , label %s" % (mtzfname, mtzlabel))
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


def load_spectra_file(spec_file, total_flux=1., pinkstride=1, as_spectrum=False):
    wavelengths, weights = np.loadtxt(spec_file, float, delimiter=',', skiprows=1).T
    wavelengths = wavelengths[::pinkstride]
    weights = weights[::pinkstride]
    energies = ENERGY_CONV/wavelengths
    FLUXES = weights / weights.sum() * total_flux
    if as_spectrum:
        return list(zip(list(wavelengths), list(FLUXES)))
    else:
        return FLUXES, energies


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

    a, b, c, al, be, ga = crystal.get_unit_cell().parameters()
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


def refls_from_sims(panel_imgs, detector, beam, thresh=0, filter=None, panel_ids=None, **kwargs):
    """
    This class is for converting the centroids in the noiseless simtbx images
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
    from dials.algorithms.spot_finding.finder import PixelListToReflectionTable

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

    pixlst_to_reftbl = PixelListToReflectionTable(
        min_spot_size=1,
        max_spot_size=194 * 184,  # TODO: change this ?
        filter_spots=FilterRunner(),  # must use a dummie filter runner!
        write_hot_pixel_mask=False)

    El = explist_from_numpyarrays(panel_imgs, detector, beam)
    iset = El.imagesets()[0]
    refls = pixlst_to_reftbl(iset, pxlst_labs)[0]

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
    #masker = MemMasker([I])
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
        print("WARNING: duplicate pids found in panel input string %s, please double check!" % panel_str)
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


#def tally_local_statistics(spot_rois):
#    # tally up all miller indices in this refinement
#    total_pix = 0
#    for x1, x2, y1, y2 in spot_rois:
#        total_pix += (x2 - x1) * (y2 - y1)
#
#    nspots = len(spot_rois)
#
#    # total_pix = self.all_pix
#    # Per image we have 3 rotation angles to refine
#    n_rot_param = 3
#
#    # by default we assume each shot refines its own ncells param (mosaic domain size Ncells_abc in nanoBragg)
#    n_per_image_ncells_param = 1 if params.simulator.crystal.has_isotroic_ncells else 3
#
#    # by default each shot refines its own unit cell parameters (e.g. a,b,c,alpha, beta, gamma)
#    n_per_image_ucell_param = n_ucell_param
#
#    # 1 crystal scale factor refined per shot (overall scale)
#    n_per_image_scale_param = 1
#
#    self.n_param_per_image = [n_rot_param + n_per_image_ncells_param + n_per_image_ucell_param +
#                          n_per_image_scale_param + 3 * n_spot_per_image[i]
#                              for i in range(n_images)]
#
#    total_per_image_unknowns = sum(self.n_param_per_image)
#
#    # NOTE: local refers to per-image
#    self.n_local_unknowns = total_per_image_unknowns
#
#    mem = self._usage()  # get memory usage
#    # note: roi para
#
#    # totals across ranks
#    if has_mpi:
#        n_images = comm.reduce(n_images, MPI.SUM, root=0)
#        n_spot_tot = comm.reduce(n_spot_tot, MPI.SUM, root=0)
#        total_pix = comm.reduce(total_pix, MPI.SUM, root=0)
#        mem = comm.reduce(mem, MPI.SUM, root=0)
#
#    # Gather so that each rank knows exactly how many local unknowns are on the other ranks
#    if has_mpi:
#        local_unknowns_per_rank = comm.gather(self.n_local_unknowns)
#    else:
#        local_unknowns_per_rank = [self.n_local_unknowns]
#
#    if rank == 0:
#        total_local_unknowns = sum(local_unknowns_per_rank)  # across all ranks
#    else:
#        total_local_unknowns = None
#
#    self.local_unknowns_across_all_ranks = total_local_unknowns
#    if has_mpi:
#        self.local_unknowns_across_all_ranks = comm.bcast(self.local_unknowns_across_all_ranks, root=0)
#
#    # TODO: what is the 2 for (its gain and detector distance which are not currently refined...
#    self.n_global_params = 2 + n_global_ucell_param + n_global_ncells_param + self.num_hkl_global  # detdist and gain + ucell params
#
#    self.n_total_unknowns = self.local_unknowns_across_all_ranks + self.n_global_params  # gain and detdist (originZ)
#
#    # also report total memory usage
#    # mem_tot = mem
#    # if has_mpi:
#    #    mem_tot = comm.reduce(mem_tot, MPI.SUM, root=0)
#
#    if has_mpi:
#        comm.Barrier()
#    if rank == 0:
#        print("\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>")
#        print("MPIWORLD TOTALZ: images=%d, spots=%d, pixels=%2.2g, Nlocal/Nglboal=%d/%d, usage=%2.2g GigaBytes"
#              % (n_images, n_spot_tot, total_pix, total_local_unknowns, self.n_global_params, mem))
#        print("Total time elapsed= %.4f seconds" % (time.time() - self.time_load_start))
#        print("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n")
#
#        # determine where in the global parameter array does this rank
#        # parameters begin
#        self.starts_per_rank = {}
#        xpos = 0
#        for _rank, n_unknown in enumerate(local_unknowns_per_rank):
#            self.starts_per_rank[_rank] = xpos
#            xpos += n_unknown
#    else:
#        self.starts_per_rank = None
#
#    if has_mpi:
#        self.starts_per_rank = comm.bcast(self.starts_per_rank, root=0)
#
#
