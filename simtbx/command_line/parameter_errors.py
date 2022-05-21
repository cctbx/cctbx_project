from __future__ import division

import os
import pandas
import numpy as np
from dials.array_family import flex

from shutil import copyfile
from simtbx.diffBragg import hopper_utils, utils
from dxtbx.model import ExperimentList


def load_expt_from_df(df, opt=False):
    """

    :param df: a hopper-formatted pandas dataframe with a single row
    :return: experiment
    """
    if opt:
        data_expt_name = df.opt_exp_name.values[0]
    else:
        data_expt_name = df.exp_name.values[0]
    assert os.path.exists(data_expt_name)
    data_expt = ExperimentList.from_file(data_expt_name)[0]
    return data_expt


def get_errors(phil_file,expt_name, refl_name, pkl_name, outfile_prefix=None, verbose=False, devid=0):
    """

    :param phil_file:
    :param expt_name:
    :param refl_name:
    :param pkl_name:
    :param outfile_prefix:
    :param verbose:
    :return:
    """
    params = utils.get_extracted_params_from_phil_sources(phil_file)
    Mod = hopper_utils.DataModeler(params)
    if not Mod.GatherFromExperiment(expt_name, refl_name):
        return
    df = pandas.read_pickle(pkl_name)
    Mod.SimulatorFromExperiment(df)
    if params.spectrum_from_imageset:
        data_expt = load_expt_from_df(df)
        spec = hopper_utils.downsamp_spec_from_params(params, data_expt)
    elif df.spectrum_filename.values[0] is not None:
        spec = utils.load_spectra_from_dataframe(df)
    else:
        data_expt = load_expt_from_df(df)
        spec = [(data_expt.beam.get_wavelength(), df.total_flux.values[0])]
    Mod.SIM.beam.spectrum = spec
    Mod.SIM.D.xray_beams = Mod.SIM.beam.xray_beams
    Mod.SIM.D.device_Id = devid
    target = hopper_utils.TargetFunc(Mod.SIM)
    # set up the refinement flags
    num_param = len(Mod.SIM.P)
    x0 = np.ones(num_param)

    vary = np.ones(num_param, bool)
    for p in Mod.SIM.P.values():
        if not p.refine:
            vary[p.xpos] = False

    target.vary = vary  # fixed flags
    target.x0 = np.array(x0, np.float64)

    if Mod.SIM.P["RotXYZ0_xtal0"].refine:
        Mod.SIM.D.refine(hopper_utils.ROTX_ID)
        Mod.SIM.D.refine(hopper_utils.ROTY_ID)
        Mod.SIM.D.refine(hopper_utils.ROTZ_ID)
    if Mod.SIM.P["Nabc0"].refine:
        Mod.SIM.D.refine(hopper_utils.NCELLS_ID)
    if Mod.SIM.P["ucell0"].refine:
        for i_ucell in range(len(Mod.SIM.ucell_man.variables)):
            Mod.SIM.D.refine(hopper_utils.UCELL_ID_OFFSET + i_ucell)
    if Mod.SIM.P["eta_abc0"].refine:
        Mod.SIM.D.refine(hopper_utils.ETA_ID)
    if Mod.SIM.P["detz_shift"].refine:
        Mod.SIM.D.refine(hopper_utils.DETZ_ID)
    if Mod.SIM.D.use_diffuse:
        Mod.SIM.D.refine(hopper_utils.DIFFUSE_ID)

    model_bragg, Jac = hopper_utils.model(x0, Mod.SIM, Mod.pan_fast_slow,compute_grad=True, dont_rescale_gradient=True)
    model_pix = model_bragg + Mod.all_background

    u = Mod.all_data - model_pix  # residuals, named "u" in notes

    sigma_rdout = params.refiner.sigma_r / params.refiner.adu_per_photon
    v = model_pix + sigma_rdout**2
    one_by_v = 1/v
    G = 1-2*u - u*u*one_by_v
    coef = one_by_v*(one_by_v*G - 2  - 2*u*one_by_v -u*u*one_by_v*one_by_v)

    coef_t = coef[Mod.all_trusted]
    Jac_t = Jac[:,Mod.all_trusted]
    # if we are only optimizing Fhkl, then the Hess is diagonal matrix
    diag_Hess = -.5*np.sum(coef_t*(Jac_t)**2, axis=1)
    with np.errstate(divide='ignore', invalid='ignore'):
        variance_s = 1/diag_Hess

    ## if we optimized per-shot scale along with Fhkl scales, then the Hess is an arrow matrix (diagonal with elem in first row/col)
    #name_to_i_Hess = {}
    #name_to_i_Hess["G_xtal0"] = 0
    #i_Hess = 1
    #for name in Mod.SIM.P:
    #    if name.startswith("scale_roi"):
    #        name_to_i_Hess[name] = i_Hess
    #        i_Hess += 1
    #Hess = np.zeros((len(name_to_i_Hess), len(name_to_i_Hess)))
    #scale_p = Mod.SIM.P["G_xtal0"]
    #overall_scale = scale_p.get_val(x0[scale_p.xpos])
    #name_from_i_Hess = {i:name for name,i in name_to_i_Hess.items()}

    #for name in name_to_i_Hess:
    #    p = Mod.SIM.P[name]
    #    xpos = p.xpos
    #    i_Hess = name_to_i_Hess[name]
    #    val = diag_Hess[xpos]
    #    Hess[i_Hess, i_Hess] = val

    ## offdiagonal terms
    #jac_coef_t = (.5*one_by_v*G)[Mod.all_trusted]
    #for i_Hess in range(1, len(name_to_i_Hess)):
    #    name = name_from_i_Hess[i_Hess]
    #    p = Mod.SIM.P[name]

    #    val_off_diag = jac_coef_t*Jac_t[p.xpos]
    #    val_off_diag = val_off_diag.sum() / overall_scale

    #    Hess[0, i_Hess] = val_off_diag
    #    Hess[i_Hess, 0] = val_off_diag

    F = Mod.SIM.crystal.miller_array
    Fmap = {h:amp for h,amp in zip(F.indices(), F.data())}
    all_I = []
    all_s = []
    all_varI = []
    #assert len(Mod.roi_id_unique) == len(Mod.refls)
    flex_varI = flex.double(len(Mod.refls),0)
    flex_I = flex.double(len(Mod.refls),0)
    sel = flex.bool(len(Mod.refls), False)

    Mod.set_slices("all_refls_idx")
    #for roi_id in Mod.roi_id_unique:
    for refl_idx in Mod.all_refls_idx_unique:
        refl_idx = int(refl_idx)
        data_slc = Mod.all_refls_idx_slices[refl_idx]
        assert len(data_slc)==1
        data_slc = data_slc[0]
        roi_id = int(Mod.roi_id[data_slc][0])
        p = Mod.SIM.P["scale_roi%d" % roi_id]
        # TODO : double check scale evaluated from x=1
        scale = p.get_val(1)
        var_s = variance_s[p.xpos]
        hkl = Mod.hi_asu_perpix[data_slc][0]
        if hkl not in Fmap:
            continue
        amp = Fmap[hkl]
        I_hkl = amp**2
        var_I = I_hkl **2 * var_s
        if var_I <= 1e-6 or var_I > 1e10:
            continue
        I = scale*I_hkl
        h,k,l=hkl
        if verbose:
            print("hkl=%d,%d,%d . I=%f +- %f" %(h,k,l, I, var_I))
        all_I.append(I)
        all_varI.append(var_I)
        all_s.append(scale)

        #refl_idx = int(Mod.all_refls_idx[data_slc][0])
        sel[refl_idx] = True
        flex_I[refl_idx] = I
        flex_varI[refl_idx] = var_I

        refl = Mod.refls[refl_idx]
        assert refl["scale_factor"] == scale

    Mod.refls["intensity.sum.value"] = flex_I
    Mod.refls["intensity.sum.variance"] = flex_varI
    Mod.refls["xyzobs.px.value"] = Mod.refls["xyzcal.px"]
    integ_refls = Mod.refls.select(sel)
    #all_s = np.array(all_s)
    #all_I = np.array(all_I)
    #all_varI = np.array(all_varI)
    #from IPython import embed;embed();exit()

    hopper_utils.free_SIM_mem(Mod.SIM)
    if outfile_prefix is not None:
        integ_refls.as_file(outfile_prefix+".refl")
        copyfile(expt_name, outfile_prefix+".expt")
    if verbose:
        print("Done.")


if __name__=="__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("phil", type=str, help="path to a simtbx.diffBragg phil file")
    parser.add_argument("expt", type=str, help="path to an experiment list file")
    parser.add_argument("refl", type=str, help="path to a reflection table")
    parser.add_argument("pkl", type=str,
                        help="path to a pandas pkl created by simtbx.diffBragg.hopper or simtbx.diffBragg.hopper_process")
    parser.add_argument("--loud", action="store_true", help="show more output")
    args = parser.parse_args()
    import logging

    logger = logging.getLogger("diffBragg.main")
    logger.setLevel(logging.DEBUG)
    get_errors(args.phil, args.expt, args.refl, args.pkl)
