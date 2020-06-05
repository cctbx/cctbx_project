from scipy import ndimage
from itertools import zip_longest
from scipy.optimize import minimize
import numpy as np
import pylab as plt
from scipy.interpolate import SmoothBivariateSpline
from cctbx import miller
from cctbx.array_family import flex


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
    sg_type = sgtbx.space_group_info(symbol=symbol).type()
    # necessary for py3 to type cast the ints
    type_casted_Hi_lst = tuple([(int(x), int(y), int(z)) for x, y, z in Hi_lst])
    Hi_flex = flex.miller_index(type_casted_Hi_lst)
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
