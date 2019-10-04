# uncompyle6 version 3.4.0
# Python bytecode 2.7 (62211)
# Decompiled from: Python 2.7.15 (default, Nov 27 2018, 21:24:58) 
# [GCC 4.2.1 Compatible Apple LLVM 10.0.0 (clang-1000.11.45.5)]
# Embedded file name: /Users/dermen/crystal/modules/cctbx_project/simtbx/diffBragg/load_ls49.py
# Compiled at: 2019-09-18 08:58:23


def strong_spot_mask(refl_tbl, img_size):
    """note only works for strong spot reflection tables
    img_size is slow-scan, fast-scan"""
    import numpy as np
    from dials.algorithms.shoebox import MaskCode
    Nrefl = len(refl_tbl)
    masks = [ refl_tbl[i]['shoebox'].mask.as_numpy_array() for i in range(Nrefl)
            ]
    code = MaskCode.Foreground.real
    x1, x2, y1, y2, z1, z2 = zip(*[ refl_tbl[i]['shoebox'].bbox for i in range(Nrefl)
                                  ])
    spot_mask = np.zeros(img_size, bool)
    for i1, i2, j1, j2, M in zip(x1, x2, y1, y2, masks):
        slcX = slice(i1, i2, 1)
        slcY = slice(j1, j2, 1)
        spot_mask[(slcY, slcX)] = M & code == code

    return spot_mask


def process_ls49_image(ls49_data_dir='/Users/dermen/crystal/modules/cctbx_project/simtbx/diffBragg/LS49_sim2'):
    import os, numpy as np, dxtbx
    from dxtbx.model.experiment_list import ExperimentListFactory
    from simtbx.diffBragg import utils
    from dials.array_family import flex
    os.chdir(ls49_data_dir)
    loader = dxtbx.load('ls49_0.npz')
    img = loader.get_raw_data().as_numpy_array()
    exp_list = ExperimentListFactory.from_json_file('idx-ls49_0_refined.expt', check_format=False)
    exp = exp_list[0]
    C = exp.crystal
    B = exp.beam
    D = exp.detector
    refls = flex.reflection_table.from_file('idx-ls49_0_integrated.refl')
    snr = refls['intensity.sum.value'] / flex.sqrt(refls['intensity.sum.variance'])
    order = np.argsort(snr)[::-1]
    refls = refls.select(snr > snr[order[20]])
    bboxes = [ list(refls['shoebox'][i].bbox)[:4] for i in range(len(refls)) ]
    bboxes = np.array(bboxes)
    bboxes[bboxes > 3000] = 2999
    mill_idx = [ list(refls['miller_index'][i]) for i in range(len(refls)) ]
    R2 = flex.reflection_table.from_file('idx-ls49_0_indexed.refl')
    strong_mask = strong_spot_mask(refl_tbl=R2, img_size=img.shape)
    is_bg_pixel = np.logical_not(strong_mask)
    num_spots = len(refls)
    tilt_abc = np.zeros((num_spots, 3))
    for i_spot, (i1, i2, j1, j2) in enumerate(bboxes):
        shoebox_img = img[j1:j2, i1:i2]
        shoebox_mask = is_bg_pixel[j1:j2, i1:i2]
        tilt, bgmask, coeff = utils.tilting_plane(shoebox_img, mask=shoebox_mask, zscore=2)
        tilt_abc[i_spot] = (
         coeff[1], coeff[2], coeff[0])

    data = np.load('LS49_data0.npz')
    spectrum = zip(data['wavelens'][33:66], data['fluxes'][33:66])
    sfall = data['sfall'][()]
    return {'dxcrystal': C, 'dxdetector': D, 'dxbeam': B, 'mill_idx': mill_idx, 'data_img': img, 'bboxes_x1x2y1y2': bboxes, 
       'tilt_abc': tilt_abc, 'spectrum': spectrum, 'sfall': sfall}


def process_ls49_image_real(tstamp='20180501143555114', #tstamp='20180501143559313',
                            ls49_data_dir='/Users/dermen/crystal/modules/cctbx_project/simtbx/diffBragg/LS49_real_data2',
                            Nstrongest = 10,
                            resmax=3.5, resmin=2.5):
    import os, pickle, numpy as np
    from scipy.interpolate import interp1d
    import dxtbx
    from dxtbx.model.experiment_list import ExperimentListFactory
    from simtbx.diffBragg import utils
    from dials.array_family import flex
    from iotbx import mtz
    os.chdir(ls49_data_dir)
    GAIN = 0.75
    loader = dxtbx.load('idx-%s.cbf' % tstamp)
    img = loader.get_raw_data().as_numpy_array() / GAIN
    exp_list = ExperimentListFactory.from_json_file('idx-%s_refined.expt' % tstamp, check_format=False)
    exp = exp_list[0]
    C = exp.crystal
    B = exp.beam
    D = exp.detector
    #refls = flex.reflection_table.from_file('idx-%s_indexed.refl' % tstamp)
    refls = flex.reflection_table.from_file('idx-%s_integrated.refl' % tstamp)
    Nbefore = len(refls)
    refls = refls.select(flex.bool([resmin < d < resmax for d in refls['d']]))
    print("Kept %d out of %d refls in the res range %2.2f to %2.2f"
          % (len(refls), Nbefore, resmin, resmax))
    snr = refls['intensity.sum.value'] / flex.sqrt(refls['intensity.sum.variance'])
    order = np.argsort(snr)[::-1]
    refls = refls.select(snr > snr[order[Nstrongest]])
    snr2 = refls['intensity.sum.value'] / flex.sqrt(refls['intensity.sum.variance'])

    bboxes = [list(refls['shoebox'][i].bbox)[:4] for i in range(len(refls)) ]
    bboxes = np.array(bboxes)
    bboxes[bboxes > 960] = 959
    mill_idx = [ list(refls['miller_index'][i]) for i in range(len(refls)) ]
    R2 = flex.reflection_table.from_file('idx-%s_indexed.refl' % tstamp)
    strong_mask = strong_spot_mask(refl_tbl=R2, img_size=img.shape)
    is_bg_pixel = np.logical_not(strong_mask)
    is_BAD_pixel = np.logical_not(pickle.load(open('mask_r4.pickle', 'r'))[0].as_numpy_array())
    is_bg_pixel[is_BAD_pixel] = False
    num_spots = len(refls)
    tilt_abc = np.zeros((num_spots, 3))
    tilts = []
    for i_spot, (i1, i2, j1, j2) in enumerate(bboxes):
        shoebox_img = img[j1:j2, i1:i2]
        shoebox_mask = is_bg_pixel[j1:j2, i1:i2]
        tilt, bgmask, coeff, _ = utils.tilting_plane(shoebox_img, mask=shoebox_mask, zscore=2)
        tilts.append(tilt)
        tilt_abc[i_spot] = (coeff[1], coeff[2], coeff[0])

    chann_lambda, channI = np.array(pickle.load(open('fee_data_r0222.pickle', 'r'))[tstamp]).T
    I = interp1d(chann_lambda, channI)
    interp_energies = np.arange(7120, 7150, 0.5)
    interp_fluxes = I(interp_energies)
    interp_fluxes /= interp_fluxes.sum()
    interp_fluxes *= 1000000000000.0
    spectrum = zip(12398.419739640716 / interp_energies, interp_fluxes)


    #M = mtz.object('ls49_oxy_2.5_s0_mark0.mtz')
    #sfall = M.as_miller_arrays_dict()[('crystal', 'dataset', 'Iobs')]

    M = mtz.object('anom_ls49_oxy_2.5_s0_mark0.mtz')
    sfall = M.as_miller_arrays_dict()[('crystal', 'dataset', 'IMEAN')]
    sfall = sfall.as_amplitude_array()
    return {'dxcrystal': C, 'dxdetector': D, 'dxbeam': B, 'mill_idx': mill_idx, 'data_img': img, 'bboxes_x1x2y1y2': bboxes, 
       'tilt_abc': tilt_abc, 'spectrum': spectrum, 'sfall': sfall,
            'mask': is_BAD_pixel}
# okay decompiling load_ls49.pyc


if __name__ == "__main__":
    from argparse import ArgumentParser
    from simtbx.diffBragg.refiners import RefineMissetAndUcell
    from simtbx.diffBragg.sim_data import SimData
    import numpy as np
    from simtbx.diffBragg.refiners.crystal_systems import MonoclinicManager
    from simtbx.diffBragg import nanoBragg_crystal, nanoBragg_beam
    from copy import deepcopy
    from IPython import embed

    parser = ArgumentParser()
    parser.add_argument("--plot", action='store_true')
    parser.add_argument("--scaleonly", action='store_true')
    args = parser.parse_args()

    if args.scaleonly:
        data = process_ls49_image_real(Nstrongest=14, resmin=2.5, resmax=12.2)
    else:
        data = process_ls49_image_real(Nstrongest=13, resmin=2.5, resmax=4.2)

    C = data["dxcrystal"]
    D = data["dxdetector"]
    B = data["dxbeam"]

    a,b,c,_,_,_ = C.get_unit_cell().parameters()
    Deff = 304 #1000
    Ncells_abc = np.power(4/3*np.pi*Deff**3/a/b/c, 1/3.)
    nbcryst = nanoBragg_crystal.nanoBragg_crystal()
    nbcryst.Ncells_abc = Ncells_abc, Ncells_abc, Ncells_abc
    nbcryst.mos_spread_deg = 0.01
    nbcryst.n_mos_domains = 1
    nbcryst.thick_mm = 0.005
    nbcryst.miller_array = data["sfall"]
    nbcryst.dxtbx_crystal = C

    nbbeam = nanoBragg_beam.nanoBragg_beam()
    nbbeam.unit_s0 = B.get_unit_s0()
    nbbeam.spectrum = data["spectrum"]

    SIM = SimData()
    SIM.crystal = nbcryst
    SIM.detector = D
    SIM.beam = nbbeam

    SIM.instantiate_diffBragg(adc_offset=0,
                              oversample=0,
                              interpolate=0,
                              verbose=0)
    SIM.D.show_params()

    ucell = C.get_unit_cell().parameters()
    UcellMan = MonoclinicManager(a=ucell[0], b=ucell[1], c=ucell[2], beta=ucell[4]*np.pi / 180.)
    import pylab as plt
    if args.plot:
        I = data["data_img"].copy()
        M = ~data["mask"]
        I*=M
        m = I[I > 0].mean()
        s = I[I > 0].std()
        vmin = m-s
        vmax=m+2.5*s

        plt.imshow(I, vmin=vmin, vmax=vmax)
        for x1, x2, y1, y2 in data["bboxes_x1x2y1y2"]:
            patch = plt.Rectangle(
                width=x2-x1,
                height=y2-y1,
                xy=(x1, y1),
                ec='r', fc='none')
            plt.gca().add_patch(patch)
        plt.show()

    RUC = RefineMissetAndUcell(
        spot_rois=data["bboxes_x1x2y1y2"],
        abc_init=data["tilt_abc"],
        img=data["data_img"],
        SimData_instance=SIM,
        plot_images=args.plot,
        ucell_manager=UcellMan,
        init_gain=1,
        init_scale=-0.3164287)
    RUC.trad_conv = True
    RUC.trad_conv_eps = 1e-5
    RUC.max_calls = 250
    RUC.refine_background_planes = False
    RUC.refine_Amatrix = True
    if args.scaleonly:
        RUC.refine_background_planes = False
        RUC.refine_Amatrix = False
    RUC.refine_crystal_scale = True
    RUC.refine_gain_fac = True
    RUC.run()
    print("Done.")
    print("Refined scale =%f", RUC.x[-1])

    C2 = deepcopy(C)
    ang, ax = RUC.get_correction_misset(as_axis_angle_deg=True)
    C2.rotate_around_origin(ax, ang)
    C2.set_B(RUC.get_refined_Bmatrix())

    # refined unit cell parameters
    ucell_ref = C2.get_unit_cell().parameters()

    print("ground truth unit cell: %2.7g,%2.7g,%2.7g,%2.7g,%2.7g,%2.7g" % ucell)
    print("refined unit cell: %2.7g,%2.7g,%2.7g,%2.7g,%2.7g,%2.7g" % ucell_ref)
    print("")

    C2.show()
    #embed()
    exit()












































    SIM.D.Umatrix = C2.get_U()
    SIM.D.Bmatrix = C2.get_B()
    SIM.D.Ncells_abc = 30, 30, 30
    s = RUC.x[-1]

    RUC = RefineMissetAndUcell(
        spot_rois=data["bboxes_x1x2y1y2"],
        abc_init=data["tilt_abc"],
        img=data["data_img"],
        SimData_instance=SIM,
        plot_images=args.plot,
        ucell_manager=UcellMan,
        init_scale=s)
    RUC.trad_conv = True
    RUC.trad_conv_eps = 1e-5
    RUC.max_calls = 250
    RUC.refine_background_planes = False
    RUC.refine_Amatrix = False #True
    RUC.refine_crystal_scale = True
    RUC.run()
    print("Done.")

    C2 = deepcopy(C)
    ang, ax = RUC.get_correction_misset(as_axis_angle_deg=True)
    C2.rotate_around_origin(ax, ang)
    C2.set_B(RUC.get_refined_Bmatrix())

    # refined unit cell parameters
    ucell_ref = C2.get_unit_cell().parameters()

    print("ground truth unit cell: %2.7g,%2.7g,%2.7g,%2.7g,%2.7g,%2.7g" % ucell)
    print("refined unit cell: %2.7g,%2.7g,%2.7g,%2.7g,%2.7g,%2.7g" % ucell_ref)
    print("")

    C2.show()
    embed()

