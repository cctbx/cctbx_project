from __future__ import print_function, division

# LIBTBX_SET_DISPATCHER_NAME fcalc_refiner

import numpy as np
from copy import deepcopy
from dials.array_family import flex
from dxtbx.model import ExperimentList
from simtbx.command_line.hopper import single_expt_pandas
import sys
from simtbx.command_line import geometry_refiner
import os
from cctbx import miller, crystal
from libtbx.mpi4py import  MPI
COMM = MPI.COMM_WORLD

import logging
MAIN_LOGGER = logging.getLogger("diffBragg.main")

from simtbx.diffBragg import hopper_utils, ensemble_refine_launcher
from simtbx.diffBragg.refiners.parameters import RangedParameter, Parameters
import pandas
import glob
from scipy.optimize import basinhopping
if COMM.rank > 0:
    sys.tracebacklimit = 0

import time
import pylab as plt


class Target(geometry_refiner.Target):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __call__(self, x, *args, **kwargs):
        self.iternum += 1
        t = time.time()
        self.x0[self.vary] = x

        f, self.g, self.sigmaZ = target_and_grad(self.x0, self.ref_params, *args, **kwargs)
        t = time.time()-t
        if COMM.rank==0:
            self.all_times.append(t)
            time_per_iter = np.mean(self.all_times)
            pred_offset_str = ", ".join(map(lambda x: "%.4f" %x, self.med_offsets))
            print("Iteration %d:\n\tResid=%f, sigmaZ %f, t-per-iter=%.4f sec, pred_offsets=%s"
                  % (self.iternum, f, self.sigmaZ, time_per_iter, pred_offset_str), flush=True)
        if self.iternum % self.save_state_freq==0 and self.iternum >0:
            if not self.overwrite_state:
                params = args[-1]  # phil params
                temp_pandas_dir = params.geometry.pandas_dir
                params.geometry.pandas_dir=params.geometry.pandas_dir + "-iter%d" % self.iternum
            med_offset = write_output_files(self.x0, self.ref_params, *args, **kwargs)
            self.med_offsets.append(med_offset)
            self.med_iternums.append(self.iternum)
            if self.plot:
                self.ax.clear()
                self.ax.plot(self.med_iternums, self.med_offsets)
                self.ax.set_ylabel("median |xobs-xcal| (pixels)")
                self.ax.set_xlabel("iteration #")
                plt.draw()
                plt.pause(0.01)
            if not self.overwrite_state:
                params.geometry.pandas_dir=temp_pandas_dir
        return f



class FcalcParameters:

    def __init__(self, phil_params, launcher):
        self.parameters = []
        for fcell_idx in launcher.asu_from_idx:
            name = "scale_fcell%d" % fcell_idx
            p = RangedParameter(name=name,
                                minval=0, maxval=1e12, fix=False,
                                center=1, beta=1e12, init=1)
            self.parameters.append(p)


class CrystalParameters:

    def __init__(self, phil_params, data_modelers):
        self.phil = phil_params
        self.parameters = []
        for i_shot in data_modelers:
            Mod = data_modelers[i_shot]

            for i_N in range(3):
                p = Mod.PAR.Nabc[i_N]
                ref_p = RangedParameter(name="rank%d_shot%d_Nabc%d" % (COMM.rank, i_shot, i_N),
                                          minval=p.minval, maxval=p.maxval, fix=self.phil.fix.Nabc, init=p.init,
                                          center=p.center, beta=p.beta)
                self.parameters.append(ref_p)

            for i_rot in range(3):
                p = Mod.PAR.RotXYZ_params[i_rot]
                ref_p = RangedParameter(name="rank%d_shot%d_RotXYZ%d" % (COMM.rank, i_shot, i_rot),
                                          minval=p.minval, maxval=p.maxval, fix=self.phil.fix.RotXYZ, init=p.init,
                                          center=p.center, beta=p.beta)
                self.parameters.append(ref_p)

            p = Mod.PAR.Scale
            ref_p = RangedParameter(name="rank%d_shot%d_Scale" % (COMM.rank, i_shot),
                                      minval=p.minval, maxval=p.maxval, fix=self.phil.fix.G, init=p.init,
                                      center=p.center, beta=p.beta)
            self.parameters.append(ref_p)

            for i_uc in range(len(Mod.PAR.ucell)):
                p = Mod.PAR.ucell[i_uc]
                ref_p = RangedParameter(name="rank%d_shot%d_Ucell%d" % (COMM.rank, i_shot, i_uc),
                                          minval=p.minval, maxval=p.maxval, fix=self.phil.fix.ucell, init=p.init,
                                          center=p.center, beta=p.beta)
                self.parameters.append(ref_p)


def model(x, ref_params, i_shot, Modeler, SIM, return_model=False):
    """

    :param x: rescaled parameter array (global)
    :param ref_params: simtbx.diffBragg.refiners.parameters.Parameters() instance
    :param i_shot: shot index for this data model,
        the simtbx.diffBragg.refiners.parameters.RangerParameter objs stored in ref_params
        have names which include i_shot
    :param Modeler: DataModeler for i_shot
    :param SIM: instance of sim_data.SimData
    :param return_model: if true, bypass the latter half of the method and return the Bragg scattering model
    :return: either the Bragg scattering model (if return_model), or else a 3-tuple of
        (float, dict of float, float)
        (negative log likelihood, gradient of negative log likelihood, average sigmaZ for the shot)
    """

    rotX = ref_params["rank%d_shot%d_RotXYZ%d" % (COMM.rank, i_shot, 0)]
    rotY = ref_params["rank%d_shot%d_RotXYZ%d" % (COMM.rank, i_shot, 1)]
    rotZ = ref_params["rank%d_shot%d_RotXYZ%d" % (COMM.rank, i_shot, 2)]
    Na = ref_params["rank%d_shot%d_Nabc%d" % (COMM.rank, i_shot, 0)]
    Nb = ref_params["rank%d_shot%d_Nabc%d" % (COMM.rank, i_shot, 1)]
    Nc = ref_params["rank%d_shot%d_Nabc%d" % (COMM.rank, i_shot, 2)]
    G = ref_params["rank%d_shot%d_Scale" % (COMM.rank, i_shot)]
    num_uc_p = len(Modeler.ucell_man.variables)
    ucell_pars = [ref_params["rank%d_shot%d_Ucell%d" % (COMM.rank, i_shot, i_uc)] for i_uc in range(num_uc_p)]

    # update the photon energy spectrum for this shot
    SIM.beam.spectrum = Modeler.spectra
    SIM.D.xray_beams = SIM.beam.xray_beams

    # update the Bmatrix
    Modeler.ucell_man.variables = [p.get_val(x[p.xpos]) for p in ucell_pars]
    Bmatrix = Modeler.ucell_man.B_recipspace
    SIM.D.Bmatrix = Bmatrix
    for i_ucell in range(len(ucell_pars)):
        SIM.D.set_ucell_derivative_matrix(
            i_ucell + hopper_utils.UCELL_ID_OFFSET,
            Modeler.ucell_man.derivative_matrices[i_ucell])

    # update the Umat rotation matrix and the RotXYZ perturbation
    SIM.D.Umatrix = Modeler.PAR.Umatrix
    SIM.D.set_value(hopper_utils.ROTX_ID, rotX.get_val(x[rotX.xpos]))
    SIM.D.set_value(hopper_utils.ROTY_ID, rotY.get_val(x[rotY.xpos]))
    SIM.D.set_value(hopper_utils.ROTZ_ID, rotZ.get_val(x[rotZ.xpos]))

    # update the mosaic block size
    SIM.D.set_ncells_values((Na.get_val(x[Na.xpos]),
                             Nb.get_val(x[Nb.xpos]),
                             Nc.get_val(x[Nc.xpos])))

    npix = int(len(Modeler.pan_fast_slow)/3.)

    # calculate the forward Bragg scattering and gradients
    SIM.D.add_diffBragg_spots(Modeler.pan_fast_slow)

    # set the scale factors per ROI, store them in dict for use when computing gradient
    perRoiScaleFactors = {}
    for fcell_idx in Modeler.fcell_idx_unique:
        p = ref_params["scale_fcell%d" % fcell_idx]
        scale_fac = p.get_val(x[p.xpos])

        # there are multiple slices because the same ASU can be present multiple times on a single image
        slices = Modeler.fcell_idx_slices[fcell_idx]
        for slc in slices:
            Modeler.per_roi_scales_per_pix[slc] = scale_fac
        perRoiScaleFactors[fcell_idx] = (scale_fac, p)

    bragg_no_scale = (SIM.D.raw_pixels_roi[:npix]).as_numpy_array()

    # get the per-shot scale factor
    scale = G.get_val(x[G.xpos])

    #combine the per-shot scale factor with the per-roi scale factors
    all_bragg_scales = scale*Modeler.per_roi_scales_per_pix

    # scale the bragg scattering
    bragg = all_bragg_scales*bragg_no_scale

    # this is the total forward model:
    model_pix = bragg + Modeler.all_background

    if return_model:
        return model_pix

    # compute the negative log Likelihood
    resid = (Modeler.all_data - model_pix)
    resid_square = resid ** 2
    V = model_pix + Modeler.nominal_sigma_rdout ** 2
    neg_LL = (.5*(np.log(2*np.pi*V) + resid_square / V))[Modeler.all_trusted].sum()

    # compute the z-score sigma as a diagnostic
    zscore_sigma = np.std((resid / np.sqrt(V))[Modeler.all_trusted])

    # store the gradients
    J = {}
    # this term is a common factor in all of the gradients
    common_grad_term = (0.5 / V * (1 - 2 * resid - resid_square / V))

    if perRoiScaleFactors:
        # the gradient in this case is the bragg scattering, scaled by only the total shot scale (G in the literature)
        bragg_no_roi = bragg_no_scale*scale

        for fcell_idx in perRoiScaleFactors:
            scale_fac, p = perRoiScaleFactors[fcell_idx]
            slices = Modeler.fcell_idx_slices[fcell_idx]
            for slc in slices:
                d = p.get_deriv(x[p.xpos], bragg_no_roi[slc])
                d_trusted = Modeler.all_trusted[slc]
                common_term_slc = common_grad_term[slc]
                if p.name in J:
                    J[p.name] += (common_term_slc*d)[d_trusted].sum()
                else:
                    J[p.name] = (common_term_slc*d)[d_trusted].sum()



    # scale factor gradients
    if not G.fix:
        bragg_no_roi_scale = bragg_no_scale*Modeler.per_roi_scales_per_pix
        scale_grad = G.get_deriv(x[G.xpos], bragg_no_roi_scale)
        J[G.name] = (common_grad_term*scale_grad)[Modeler.all_trusted].sum()

    # Umat gradients
    for i_rot, rot in enumerate([rotX, rotY, rotZ]):
        if not rot.fix:
            rot_db_id = geometry_refiner.ROTXYZ_ID[i_rot]
            rot_grad = scale*SIM.D.get_derivative_pixels(rot_db_id).as_numpy_array()[:npix]
            rot_grad = rot.get_deriv(x[rot.xpos], rot_grad)
            J[rot.name] = (common_grad_term*rot_grad)[Modeler.all_trusted].sum()

    # mosaic block size gradients
    if not Na.fix:
        Nabc_grad = SIM.D.get_ncells_derivative_pixels()
        for i_N, N in enumerate([Na, Nb, Nc]):
            N_grad = scale*(Nabc_grad[i_N][:npix].as_numpy_array())
            N_grad = N.get_deriv(x[N.xpos], N_grad)
            J[N.name] = (common_grad_term*N_grad)[Modeler.all_trusted].sum()

    # unit cell gradients
    if not ucell_pars[0].fix:
        for i_ucell, uc_p in enumerate(ucell_pars):
            d = scale*SIM.D.get_derivative_pixels(hopper_utils.UCELL_ID_OFFSET+i_ucell).as_numpy_array()[:npix]
            d = uc_p.get_deriv(x[uc_p.xpos], d)
            J[ucell_pars[i_ucell].name] = (common_grad_term*d)[Modeler.all_trusted].sum()

    # detector model gradients
    detector_derivs = []
    for diffbragg_parameter_id in geometry_refiner.PAN_OFS_IDS+geometry_refiner.PAN_XYZ_IDS:
        try:
            d = common_grad_term*scale*(SIM.D.get_derivative_pixels(diffbragg_parameter_id).as_numpy_array()[:npix])
        except ValueError:
            d = None
        detector_derivs.append(d)
    names = "RotOrth", "RotFast", "RotSlow", "ShiftX", "ShiftY", "ShiftZ"
    for group_id in Modeler.unique_panel_group_ids:
        for name in names:
            J["group%d_%s" % (group_id, name)] = 0
        for pixel_rng in Modeler.group_id_slices[group_id]:
            trusted_pixels = Modeler.all_trusted[pixel_rng]
            for i_name, name in enumerate(names):
                par_name = "group%d_%s" % (group_id, name)
                det_param = ref_params[par_name]
                if det_param.fix:
                    continue
                pixderivs = detector_derivs[i_name][pixel_rng][trusted_pixels]
                pixderivs = det_param.get_deriv(x[det_param.xpos], pixderivs)
                J[par_name] += pixderivs.sum()

    return neg_LL, J, model_pix, zscore_sigma



def target_and_grad(x, ref_params, data_modelers, SIM, params):
    """
    Returns the target functional and the gradients
    :param x: float array of parameter values as seen by scipt.optimize (rescaled)
    :param ref_params: refinement parameter objects (diffBragg.refiners.parameters.Parameters() )
    :param data_modelers: dict of data modelers (one per experiment)
    :param SIM: sim_data instance
    :param params: phil parameters
    :return: 2-tuple, target and gradients
    """
    target_functional = 0
    grad = np.zeros(len(x))

    save_name = params.geometry.optimized_detector_name
    if not all(params.geometry.fix.panel_rotations) and not all(params.geometry.fix.panel_rotations):
        geometry_refiner.update_detector(x, ref_params, SIM, save_name)

    all_shot_sigZ = []
    for i_shot in data_modelers:
        Modeler = data_modelers[i_shot]

        neg_LL, neg_LL_grad, model_pix, per_shot_sigZ = model(x, ref_params, i_shot, Modeler, SIM)
        all_shot_sigZ.append(per_shot_sigZ)

        # accumulate the target functional for this rank/shot
        target_functional += neg_LL

        if params.use_restraints:
            for name in ref_params:
                par = ref_params[name]
                if not par.is_global and not par.fix:
                    val = par.get_restraint_val(x[par.xpos])
                    target_functional += val

        # accumulate the gradients for this rank/shot
        for name in ref_params:
            if name in neg_LL_grad:
                par = ref_params[name]
                grad[par.xpos] += neg_LL_grad[name]
                # for restraints only update the per-shot restraint gradients here
                if params.use_restraints and not par.is_global and not par.fix:
                    grad[par.xpos] += par.get_restraint_deriv(x[par.xpos])

    # sum the target functional and the gradients across all ranks
    target_functional = COMM.bcast(COMM.reduce(target_functional))
    grad = COMM.bcast(COMM.reduce(grad))

    if params.use_restraints and params.geometry.betas.close_distances is not None:
        target_functional += np.std(SIM.D.close_distances) / params.geometry.betas.close_distances

    ## add in the detector parameter restraints
    if params.use_restraints:
        for name in ref_params:
            par = ref_params[name]
            if par.is_global and not par.fix:
                target_functional += par.get_restraint_val(x[par.xpos])
                grad[par.xpos] += par.get_restraint_deriv(x[par.xpos])

    all_shot_sigZ = COMM.reduce(all_shot_sigZ)
    if COMM.rank == 0:
        all_shot_sigZ = np.median(all_shot_sigZ)

    return target_functional, grad, all_shot_sigZ


def load_pkl(params):
    if params.geometry.input_pkl is not None:
        df = pandas.read_pickle(params.geometry.input_pkl)
    else:
        assert params.geometry.input_pkl_glob is not None
        fnames = glob.glob(params.geometry.input_pkl_glob)
        dfs = []
        for i_f, f in enumerate(fnames):
            if i_f % COMM.size != COMM.rank:
                continue
            if COMM.rank==0:
                print("Loaing hopper pkl %d / %d" %(i_f+1, len(fnames)), flush=True)
            df_i = pandas.read_pickle(f)
            dfs.append(df_i)
        dfs = COMM.reduce(dfs)
        if COMM.rank==0:
            df = pandas.concat(dfs)
        else:
            df = None
        df = COMM.bcast(df)

    if params.skip is not None:
        df = df.iloc[params.skip:]
    if params.first_n is not None:
        df = df.iloc[:params.first_n]
    return df


def toggle_refined(params, launcher):
    # configure diffBragg instance for gradient computation
    if not params.fix.RotXYZ:
        for i_rot in range(3):
            launcher.SIM.D.refine(geometry_refiner.ROTXYZ_ID[i_rot])
    if not params.fix.Nabc:
        launcher.SIM.D.refine(hopper_utils.NCELLS_ID)
    if not params.fix.ucell:
        for i_ucell in range(launcher.SIM.num_ucell_param):
            launcher.SIM.D.refine(hopper_utils.UCELL_ID_OFFSET + i_ucell)
    for i, diffbragg_id in enumerate(geometry_refiner.PAN_OFS_IDS):
        if not params.geometry.fix.panel_rotations[i]:
            launcher.SIM.D.refine(diffbragg_id)

    for i, diffbragg_id in enumerate(geometry_refiner.PAN_XYZ_IDS):
        if not params.geometry.fix.panel_translations[i]:
            launcher.SIM.D.refine(diffbragg_id)


def fcalc_min(params):
    """
    :param params: phil parameters (simtbx/diffBragg/phil.py)
    """

    launcher = ensemble_refine_launcher.RefineLauncher(params)
    df = load_pkl(params)

    pdir = params.geometry.pandas_dir
    assert pdir is not None, "provide a pandas_dir where output files will be generated"
    params.geometry.optimized_detector_name = os.path.join(pdir, os.path.basename(params.geometry.optimized_detector_name))
    if COMM.rank==0:
        if not os.path.exists(pdir):
            os.makedirs(pdir)
    if COMM.rank == 0:
        print("Will optimize using %d experiments" %len(df))
    launcher.load_inputs(df, refls_key=params.geometry.refls_key)
    launcher.SIM.asu_from_idx = launcher.asu_from_idx

    for i_shot in launcher.Modelers:
        Modeler = launcher.Modelers[i_shot]

        # each pixel corresponds to 1 hkl (roughly), fcell_idx is a unique global (all shot) mapping  for hkls
        Modeler.fcell_idx = np.array([launcher.idx_from_asu[h] for h in Modeler.hi_asu_perpix])
        Modeler.set_slices("fcell_idx")
        Modeler.per_roi_scales_per_pix = np.ones_like(Modeler.all_data)

        geometry_refiner.set_group_id_slices(Modeler, launcher.panel_group_from_id)

    # same on every rank:
    det_params = geometry_refiner.DetectorParameters(params, launcher.panel_groups_refined, launcher.n_panel_groups)

    fcalc_params = FcalcParameters(params, launcher)

    # different on each rank
    crystal_params = CrystalParameters(params,launcher.Modelers)
    crystal_params.parameters = COMM.bcast(COMM.reduce(crystal_params.parameters))

    LMP = Parameters()
    for p in crystal_params.parameters + det_params.parameters + fcalc_params.parameters:
        LMP.add(p)

    # attached some objects to SIM for convenience
    launcher.SIM.panel_reference_from_id = launcher.panel_reference_from_id
    launcher.SIM.panel_group_from_id = launcher.panel_group_from_id
    launcher.SIM.panel_groups_refined = launcher.panel_groups_refined

    # set the GPU device
    launcher.SIM.D.device_Id = COMM.rank % params.refiner.num_devices
    npx_str = "(rnk%d, dev%d): %d pix" %(COMM.rank, launcher.SIM.D.device_Id, launcher.NPIX_TO_ALLOC)
    npx_str = COMM.gather(npx_str)
    if COMM.rank==0:
        print("How many pixels each rank will allocate for on its device:")
        print("; ".join(npx_str))
    launcher.SIM.D.Npix_to_allocate = launcher.NPIX_TO_ALLOC
    toggle_refined(params, launcher)

    # do a barrel roll!
    target = Target(LMP, save_state_freq=params.geometry.save_state_freq, overwrite_state=params.geometry.save_state_overwrite)
    fcn_args = (launcher.Modelers, launcher.SIM, params)
    lbfgs_kws = {"jac": target.jac,
                 "method": "L-BFGS-B",
                 "args": fcn_args,
                 "options":  {"ftol": params.ftol, "gtol": 1e-10, "maxfun":1e5, "maxiter":params.lbfgs_maxiter}}

    result = basinhopping(target, target.x0[target.vary],
                 niter=params.niter,
                 minimizer_kwargs=lbfgs_kws,
                 T=params.temp,
                 callback=target.at_min_callback,
                 disp=False,
                 stepsize=params.stepsize)

    target.x0[target.vary] = result.x
    Xopt = target.x0  # optimized, rescaled parameters

    if params.geometry.optimized_results_tag is not None:
        write_output_files(Xopt, LMP, launcher.Modelers, launcher.SIM, params)

    if COMM.rank == 0:
        geometry_refiner.save_opt_det(params, target.x0, target.ref_params, launcher.SIM)


def write_output_files(Xopt, LMP, Modelers, SIM, params):
    """
    Writes refl and exper files for each experiment modeled during
    the ensemble refiner
    :param Xopt: float array of optimized rescaled parameter values
    :param LMP: simtbx.diffBragg.refiners.parameters.Parameters() object
    :param Modelers: data modelers (launcher.Modleers
    :param SIM: instance of sim_data (launcher.SIM)
    :param params: phil params, simtbx.diffBragg.phil.py
    """
    opt_det = geometry_refiner.get_optimized_detector(Xopt, LMP, SIM)

    # Store the hessian of negative log likelihood for error estimation
    # must determine total number of refined Fhkls and then create a vector of 0's of that length
    num_fhkl_param = 0
    for name in LMP:
        if "fcell" in name:
            num_fhkl_param += 1
    diag_hess = np.zeros(num_fhkl_param)

    if params.geometry.pandas_dir is not None and COMM.rank == 0:
        if not os.path.exists(params.geometry.pandas_dir):
            os.makedirs(params.geometry.pandas_dir)
        refdir = os.path.join(params.geometry.pandas_dir, "refls")
        expdir = os.path.join(params.geometry.pandas_dir, "expts")
        for dname in [refdir, expdir]:
            if not os.path.exists(dname):
                os.makedirs(dname)

    all_shot_pred_offsets = []
    for i_shot in Modelers:
        Modeler = Modelers[i_shot]
        # these are in simtbx.diffBragg.refiners.parameters.RangedParameter objects
        rotX = LMP["rank%d_shot%d_RotXYZ%d" % (COMM.rank, i_shot, 0)]
        rotY = LMP["rank%d_shot%d_RotXYZ%d" % (COMM.rank, i_shot, 1)]
        rotZ = LMP["rank%d_shot%d_RotXYZ%d" % (COMM.rank, i_shot, 2)]
        num_uc_p = len(Modeler.ucell_man.variables)
        ucell_pars = [LMP["rank%d_shot%d_Ucell%d" % (COMM.rank, i_shot, i_uc)] for i_uc in range(num_uc_p)]

        # convert rotation angles back to radians (thats what the parameters.RangedParamter.get_val method does)
        rotXYZ = rotX.get_val(Xopt[rotX.xpos]), \
                 rotY.get_val(Xopt[rotY.xpos]), \
                 rotZ.get_val(Xopt[rotZ.xpos])

        # ucell_man is an instance of
        # simtbx.diffBragg.refiners.crystal_systems.manager.Manager()
        # (for the correct xtal system)
        Modeler.ucell_man.variables = [p.get_val(Xopt[p.xpos]) for p in ucell_pars]
        ucpar = Modeler.ucell_man.unit_cell_parameters

        new_crystal = hopper_utils.new_cryst_from_rotXYZ_and_ucell(rotXYZ, ucpar, Modeler.E.crystal)
        new_exp = deepcopy(Modeler.E)
        new_exp.crystal = new_crystal
        wave, wt = map(np.array, zip(*Modeler.spectra))
        ave_wave = (wave*wt).sum()/wt.sum()
        new_exp.beam.set_wavelength(ave_wave)
        new_exp.detector = opt_det

        Modeler.best_model = model(Xopt, LMP, i_shot, Modeler, SIM, return_model=True)
        Modeler.best_model_includes_background = True

        # Get the bragg-only component of model in order to compute hessian terms
        bragg = Modeler.best_model - Modeler.all_background

        # store the updated per-roi scale factors in the new refl table
        roi_scale_factor = flex.double(len(Modeler.refls), 1)
        for ii, fcell_idx in enumerate(Modeler.fcell_idx_unique):
            p = LMP["scale_fcell%d" % fcell_idx]
            scale_fac = p.get_val(Xopt[p.xpos])
            slices = Modeler.fcell_idx_slices[fcell_idx]
            for slc in slices:
                # update the refl table column
                roi_refl_ids = Modeler.all_refls_idx[slc]
                unique_refl_ids = np.unique(roi_refl_ids)
                for refl_idx in unique_refl_ids:
                    roi_scale_factor[refl_idx] = scale_fac

                # update the hessian of the log likelihood
                # first derivative is the Bragg component of the model divided by the scale factor
                # TODO what if scale_fac is close to 0 ?
                first_deriv = bragg[slc] / scale_fac
                u = Modeler.all_data[slc] - Modeler.best_model[slc]
                v = Modeler.best_model[slc] + Modeler.nominal_sigma_rdout**2
                one_by_v = 1 / v
                G = 1 - 2 * u - u * u * one_by_v
                hessian_coef = one_by_v * (one_by_v * G - 2 - 2 * u * one_by_v - u * u * one_by_v * one_by_v)
                trusted_slc = Modeler.all_trusted[slc]
                diag_hess[fcell_idx] += -0.5*(hessian_coef * (first_deriv**2))[trusted_slc].sum()

        Modeler.refls["global_scale_factor"] = roi_scale_factor

        # get the new refls
        new_refl = hopper_utils.get_new_xycalcs(Modeler, new_exp, old_refl_tag="before_geom_ref")

        new_refl_fname, refl_ext = os.path.splitext(Modeler.refl_name)
        new_refl_fname = "rank%d_%s_%s%s" % (COMM.rank, os.path.basename(new_refl_fname), params.geometry.optimized_results_tag, refl_ext)
        if not new_refl_fname.endswith(".refl"):
            new_refl_fname += ".refl"
        new_refl_fname = os.path.join(params.geometry.pandas_dir,"refls",  new_refl_fname)
        new_refl.as_file(new_refl_fname)
        shot_pred_offsets = geometry_refiner.get_dist_from_R(new_refl)
        all_shot_pred_offsets += list(shot_pred_offsets)

        new_expt_fname, expt_ext = os.path.splitext(Modeler.exper_name)
        new_expt_fname = "rank%d_%s_%s%s" % (COMM.rank, os.path.basename(new_expt_fname), params.geometry.optimized_results_tag, expt_ext)

        if not new_expt_fname.endswith(".expt"):
            new_expt_fname += ".expt"

        new_expt_fname = os.path.join(params.geometry.pandas_dir,"expts", new_expt_fname)
        new_exp_lst = ExperimentList()
        new_exp_lst.append(new_exp)
        new_exp_lst.as_file(new_expt_fname)

        if params.geometry.pandas_dir is not None:
            a,b,c,al,be,ga = ucpar
            ncells_p = [LMP["rank%d_shot%d_Nabc%d" % (COMM.rank, i_shot, i)] for i in range(3)]
            Na,Nb,Nc = [p.get_val(Xopt[p.xpos]) for p in ncells_p]
            scale_p = LMP["rank%d_shot%d_Scale" %(COMM.rank, i_shot)]
            scale = scale_p.get_val(Xopt[scale_p.xpos])

            _,fluxes = zip(*SIM.beam.spectrum)
            eta_a = eta_b = eta_c = np.nan
            df= single_expt_pandas(xtal_scale=scale, Amat=new_crystal.get_A(),
                                   ncells_abc=(Na, Nb, Nc), ncells_def=(0,0,0),
                                   eta_abc=(eta_a, eta_b, eta_c),
                                   diff_gamma=(np.nan, np.nan, np.nan),
                                   diff_sigma=(np.nan, np.nan, np.nan),
                                   detz_shift=0,
                                   use_diffuse=params.use_diffuse_models,
                                   gamma_miller_units=params.gamma_miller_units,
                                   eta=np.nan,
                                   rotXYZ=tuple(rotXYZ),
                                   ucell_p = (a,b,c,al,be,ga),
                                   ucell_p_init=(np.nan, np.nan, np.nan, np.nan, np.nan, np.nan),
                                   lam0_lam1 = (np.nan, np.nan),
                                   spec_file=Modeler.spec_name,
                                   spec_stride=params.simulator.spectrum.stride,
                                   flux=sum(fluxes), beamsize_mm=SIM.beam.size_mm,
                                   orig_exp_name=Modeler.exper_name,
                                   opt_exp_name=os.path.abspath(new_expt_fname),
                                   spec_from_imageset=params.spectrum_from_imageset,
                                   oversample=SIM.D.oversample,
                                   opt_det=params.opt_det, stg1_refls=Modeler.refl_name, stg1_img_path=None)
            pandas_name = os.path.splitext(os.path.basename(new_expt_fname))[0] + ".pkl"
            pandas_name = os.path.join(params.geometry.pandas_dir, pandas_name)
            df.to_pickle(pandas_name)
            modeler_name = pandas_name.replace(".pkl", ".npy")
            np.save(modeler_name, Modeler)

    all_shot_pred_offsets = COMM.reduce(all_shot_pred_offsets)
    if COMM.rank==0:
        median_pred_offset = np.median(all_shot_pred_offsets)
    else:
        median_pred_offset = None
    median_pred_offset = COMM.bcast(median_pred_offset)

    # reduce the hessian over all shots then compute the errors of the structure factors
    diag_hess = COMM.reduce(diag_hess)

    uc_p = np.zeros(6)
    nshot = 0
    for i_shot in Modelers:
        Mod = Modelers[i_shot]
        num_uc_p = len(Mod.ucell_man.variables)
        ucell_pars = [LMP["rank%d_shot%d_Ucell%d" % (COMM.rank, i_shot, i_uc)] for i_uc in range(num_uc_p)]
        Mod.ucell_man.variables = [p.get_val(Xopt[p.xpos]) for p in ucell_pars]
        uc_p += np.array(Mod.ucell_man.unit_cell_parameters)
        nshot += 1
    nshot = COMM.reduce(nshot)
    uc_p = COMM.reduce(uc_p)


    if COMM.rank==0:

        ave_uc_p = uc_p / nshot

        fhkl_file = os.path.join(params.geometry.pandas_dir, "final_merge.mtz")

        F = SIM.crystal.miller_array
        Fmap = {h: amp for h, amp in zip(F.indices(), F.data())}

        with np.errstate(divide='ignore', invalid='ignore'):
            scale_variance = 1 / diag_hess

        indices = flex.miller_index()
        data = flex.double()
        sigmas = flex.double()
        for fcell_idx in range(num_fhkl_param):

            pname = "scale_fcell%d" % fcell_idx
            p = LMP[pname]
            scale = p.get_val(Xopt[p.xpos])

            hkl = SIM.asu_from_idx[fcell_idx]
            F_no_scale = Fmap[hkl]
            Ihkl = scale* F_no_scale**2
            Fhkl = np.sqrt(Ihkl)
            var_scale = scale_variance[fcell_idx]
            if var_scale <= 0:
                continue
            sig_F = 0.5*F_no_scale / np.sqrt(scale) * np.sqrt(var_scale)
            if np.isinf(sig_F):
                continue
            indices.append(hkl)
            data.append(Fhkl)
            sigmas.append(sig_F)
        # store an optimized mtz, and a numpy array with the same information

        sym = crystal.symmetry(tuple(ave_uc_p), SIM.crystal.symbol)
        mset = miller.set(sym, indices, True)
        ma = miller.array(mset, data, sigmas)
        ma = ma.set_observation_type_xray_amplitude().as_anomalous_array()
        ma.as_mtz_dataset(column_root_label="F").mtz_object().write(fhkl_file)

    return median_pred_offset


if __name__ == "__main__":
    from argparse import ArgumentParser
    from libtbx.phil import parse
    from simtbx.diffBragg.phil import philz, hopper_phil

    parser = ArgumentParser()
    parser.add_argument("--phil", type=str, required=True, help="path to a phil string")
    parser.add_argument("--cmdlinePhil", nargs="+", default=None, type=str, help="command line phil params")
    progargs = parser.parse_args()

    phil_scope = parse(philz+hopper_phil)
    arg_interp = phil_scope.command_line_argument_interpreter(home_scope="")

    phil_file = open(progargs.phil, "r").read()
    user_phil = parse(phil_file)
    phil_sources = [user_phil]

    if progargs.cmdlinePhil is not None:
        command_line_phils = [arg_interp.process(phil) for phil in progargs.cmdlinePhil]
        phil_sources += command_line_phils

    working_phil, unused = phil_scope.fetch(sources=phil_sources, track_unused_definitions=True)
    for loc in unused:
        print("WARNING: unused phil:", loc)

    params = working_phil.extract()
    fcalc_min(params)
