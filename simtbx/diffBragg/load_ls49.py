

def strong_spot_mask(refl_tbl, img_size):
    """note only works for strong spot reflection tables
    img_size is slow-scan, fast-scan"""
    import numpy as np
    from dials.algorithms.shoebox import MaskCode
    Nrefl = len(refl_tbl)
    masks = [refl_tbl[i]['shoebox'].mask.as_numpy_array()
             for i in range(Nrefl)]
    code = MaskCode.Foreground.real

    x1, x2, y1, y2, z1, z2 = zip(*[refl_tbl[i]['shoebox'].bbox
                                   for i in range(Nrefl)])
    spot_mask = np.zeros(img_size, bool)
    for i1, i2, j1, j2, M in zip(x1, x2, y1, y2, masks):
        slcX = slice(i1, i2, 1)
        slcY = slice(j1, j2, 1)
        spot_mask[slcY, slcX] = M & code == code
    return spot_mask


def process_ls49_image(ls49_data_dir="/Users/dermen/crystal/modules/cctbx_project/simtbx/diffBragg/LS49_sim2"):
    import os
    import numpy as np
    import dxtbx
    from dxtbx.model.experiment_list import ExperimentListFactory
    from simtbx.diffBragg import utils
    from dials.array_family import flex

    os.chdir(ls49_data_dir)

    loader = dxtbx.load("ls49_0.npz")
    img = loader.get_raw_data().as_numpy_array()

    exp_list = ExperimentListFactory.from_json_file("idx-ls49_0_refined.expt", check_format=False)
    exp = exp_list[0]

    C = exp.crystal
    B = exp.beam
    D = exp.detector

    refls = flex.reflection_table.from_file("idx-ls49_0_integrated.refl")

    # NOTE: filter reflections
    #sel = np.zeros(len(refls), bool)
    #n_keep_refls = 30
    #sel[np.random.permutation(len(refls))[:n_keep_refls]] = True
    #refls = refls.select(flex.bool(sel))

    snr = refls["intensity.sum.value"]/flex.sqrt(refls["intensity.sum.variance"])
    order = np.argsort(snr)[::-1]
    refls = refls.select(snr > snr[order[20]])

    bboxes = [list(refls["shoebox"][i].bbox)[:4] for i in range(len(refls))]
    bboxes = np.array(bboxes)
    bboxes[bboxes > 3000] = 2999
    mill_idx = [list(refls["miller_index"][i]) for i in range(len(refls))]

    R2 = flex.reflection_table.from_file("idx-ls49_0_indexed.refl")
    strong_mask = strong_spot_mask(refl_tbl=R2, img_size=img.shape)

    is_bg_pixel = np.logical_not(strong_mask)

    # now fit tilting planes
    num_spots = len(refls)
    tilt_abc = np.zeros((num_spots, 3))

    for i_spot, (i1, i2, j1, j2) in enumerate(bboxes):
        shoebox_img = img[j1:j2, i1:i2]
        shoebox_mask = is_bg_pixel[j1:j2, i1:i2]

        tilt, bgmask, coeff = utils.tilting_plane(
            shoebox_img,
            mask=shoebox_mask,  # mask specifies which spots are bg pixels...
            zscore=2)

        tilt_abc[i_spot] = coeff[1], coeff[2], coeff[0]  # store as fast-scan coeff, slow-scan coeff, offset coeff

    data = np.load("LS49_data0.npz")
    spectrum = zip(data["wavelens"][33:66], data["fluxes"][33:66])
    sfall = data["sfall"][()]

    return {"dxcrystal": C, "dxdetector": D, "dxbeam": B, "mill_idx": mill_idx, "data_img": img,
            "bboxes_x1x2y1y2": bboxes, "tilt_abc": tilt_abc, "spectrum": spectrum, "sfall": sfall}


def process_ls49_image_real(tstamp="20180501143559313",
                            ls49_data_dir="/Users/dermen/crystal/modules/cctbx_project/simtbx/diffBragg/LS49_real_data2"):
    import os
    import pickle
    import numpy as np
    from scipy.interpolate import interp1d
    import dxtbx
    from dxtbx.model.experiment_list import ExperimentListFactory
    from simtbx.diffBragg import utils
    from dials.array_family import flex
    from iotbx import mtz

    os.chdir(ls49_data_dir)

    GAIN = 0.46

    loader = dxtbx.load("idx-%s.cbf"%tstamp)
    img = loader.get_raw_data().as_numpy_array() / GAIN

    exp_list = ExperimentListFactory.from_json_file("idx-%s_refined.expt" % tstamp, check_format=False)
    exp = exp_list[0]

    C = exp.crystal
    B = exp.beam
    D = exp.detector

    refls = flex.reflection_table.from_file("idx-%s_integrated.refl" %tstamp)

    snr = refls["intensity.sum.value"]/flex.sqrt(refls["intensity.sum.variance"])
    order = np.argsort(snr)[::-1]
    Nstrongest = 10
    refls = refls.select(snr > snr[order[Nstrongest]])

    bboxes = [list(refls["shoebox"][i].bbox)[:4] for i in range(len(refls))]
    bboxes = np.array(bboxes)
    bboxes[bboxes > 960] = 959  # TODO: fixme in the code
    mill_idx = [list(refls["miller_index"][i]) for i in range(len(refls))]

    R2 = flex.reflection_table.from_file("idx-%s_indexed.refl" % tstamp)
    strong_mask = strong_spot_mask(refl_tbl=R2, img_size=img.shape)

    is_bg_pixel = np.logical_not(strong_mask)

    is_BAD_pixel = np.logical_not(pickle.load(open("mask_r4.pickle", "r"))[0].as_numpy_array())
    is_bg_pixel[is_BAD_pixel] = False

    # now fit tilting planes
    num_spots = len(refls)
    tilt_abc = np.zeros((num_spots, 3))
    tilts = []
    for i_spot, (i1, i2, j1, j2) in enumerate(bboxes):
        shoebox_img = img[j1:j2, i1:i2]
        shoebox_mask = is_bg_pixel[j1:j2, i1:i2]

        tilt, bgmask, coeff = utils.tilting_plane(
            shoebox_img,
            mask=shoebox_mask,  # mask specifies which spots are bg pixels...
            zscore=2)
        tilts.append( tilt)
        tilt_abc[i_spot] = coeff[1], coeff[2], coeff[0]  # store as fast-scan coeff, slow-scan coeff, offset coeff

    chann_lambda, channI = np.array(pickle.load(open("fee_data_r0222.pickle", "r"))[tstamp]).T
    I = interp1d(chann_lambda, channI)
    interp_energies = np.arange(7120, 7150, 0.5)  # half eV resolution
    interp_fluxes = I(interp_energies)

    interp_fluxes /= interp_fluxes.sum()
    interp_fluxes *= 1e12
    spectrum = zip(12398.419739640716/interp_energies, interp_fluxes)

    M = mtz.object("ls49_oxy_2.5_s0_mark0.mtz")
    sfall = M.as_miller_arrays_dict()[('crystal', 'dataset', 'Iobs')]
    sfall = sfall.as_amplitude_array()

    return {"dxcrystal": C, "dxdetector": D, "dxbeam": B, "mill_idx": mill_idx, "data_img": img,
            "bboxes_x1x2y1y2": bboxes, "tilt_abc": tilt_abc, "spectrum": spectrum, "sfall": sfall}
