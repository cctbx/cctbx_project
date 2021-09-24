
import numpy as np
import lmfit
from libtbx.mpi4py import  MPI
COMM = MPI.COMM_WORLD

import logging
MAIN_LOGGER = logging.getLogger("main")

from simtbx.diffBragg import hopper_utils, ensemble_refine_launcher

# diffBragg internal indices for derivative manager
ROTXYZ_ID = hopper_utils.ROTXYZ_IDS
PAN_O_ID = 14
PAN_F_ID = 17
PAN_S_ID = 18
PAN_X_ID = 15
PAN_Y_ID = 16
PAN_Z_ID = 10
PAN_OFS_IDS = PAN_O_ID, PAN_F_ID, PAN_S_ID
PAN_XYZ_IDS = PAN_X_ID, PAN_Y_ID, PAN_Z_ID

DEG_TO_PI = np.pi/180.


class DetectorParameters:

    def __init__(self, phil_params, panel_groups_refined, num_panel_groups ):

        self.parameters = []
        GEO = phil_params.geometry
        for i_group in range(num_panel_groups):
            group_has_data = i_group in panel_groups_refined
            if not group_has_data:
                continue
            vary_rots = [not fixed_flag and group_has_data for fixed_flag in GEO.fix.panel_rotations]
            o = lmfit.Parameter(name="group%d_RotOrth" % i_group, value=0, #GEO.init.panel_rotations[2],
                                min=GEO.min.panel_rotations[0]*DEG_TO_PI, max=GEO.max.panel_rotations[0]*DEG_TO_PI,
                                vary=vary_rots[0])
            f = lmfit.Parameter(name="group%d_RotFast" % i_group, value=0, #GEO.init.panel_rotations[0],
                                min=GEO.min.panel_rotations[1]*DEG_TO_PI, max=GEO.max.panel_rotations[1]*DEG_TO_PI,
                                vary=vary_rots[1])
            s = lmfit.Parameter(name="group%d_RotSlow" % i_group, value=0, #GEO.init.panel_rotations[1],
                                min=GEO.min.panel_rotations[2]*DEG_TO_PI, max=GEO.max.panel_rotations[2]*DEG_TO_PI,
                                vary=vary_rots[2])

            vary_shifts = [not fixed_flag and group_has_data for fixed_flag in GEO.fix.panel_translations]
            x = lmfit.Parameter(name="group%d_ShiftX" % i_group, value=0, #GEO.init.panel_translations[0],
                                min=GEO.min.panel_translations[0]*1e-3, max=GEO.max.panel_translations[0]*1e-3,
                                vary=vary_shifts[0])
            y = lmfit.Parameter(name="group%d_ShiftY" % i_group, value=0, #GEO.init.panel_translations[1],
                                min=GEO.min.panel_translations[1]*1e-3, max=GEO.max.panel_translations[1]*1e-3,
                                vary=vary_shifts[1])
            z = lmfit.Parameter(name="group%d_ShiftZ" % i_group, value=0, #GEO.init.panel_translations[2],
                                min=GEO.min.panel_translations[2]*1e-3, max=GEO.max.panel_translations[2]*1e-3,
                                vary=vary_shifts[2])

            self.parameters += [f, s, o, x, y, z]


class CrystalParameters:

    def __init__(self, data_modelers):

        self.parameters = []
        for i_shot in data_modelers:
            Mod = data_modelers[i_shot]

            for i_N in range(3):
                p = Mod.PAR.Nabc[i_N]
                lmfit_p = lmfit.Parameter("rank%d_shot%d_Nabc%d" % (COMM.rank, i_shot, i_N),
                                          min=p.minval, max=p.maxval, vary=not p.fix, value=p.init)
                self.parameters.append(lmfit_p)

            for i_rot in range(3):
                p = Mod.PAR.RotXYZ_params[i_rot]
                lmfit_p = lmfit.Parameter("rank%d_shot%d_RotXYZ%d" % (COMM.rank, i_shot, i_rot),
                                          min=p.minval, max=p.maxval, vary=not p.fix, value=p.init)
                self.parameters.append(lmfit_p)

            p = Mod.PAR.Scale
            lmfit_p = lmfit.Parameter("rank%d_shot%d_Scale" % (COMM.rank, i_shot),
                                      min=p.minval, max=p.maxval, vary=not p.fix, value=p.init)
            self.parameters.append(lmfit_p)

            for i_uc in range(len(Mod.PAR.ucell)):
                p = Mod.PAR.ucell[i_uc]
                lmfit_p = lmfit.Parameter("rank%d_shot%d_Ucell%d" % (COMM.rank, i_shot, i_uc),
                                          min=p.minval, max=p.maxval, vary=not p.fix, value=p.init)
                self.parameters.append(lmfit_p)


class Target:
    def __init__(self):
        pass

    def callbk(self, lmfit_params, iter, resid, *fcn_args, **fcn_kws):
        if COMM.rank==0:
            print("Iteration %d:\n\tResid=%f, sigmaZ %f" % (iter, resid, self.sigmaZ))
            print("\tAverage det shift=%f, average det rot=%f" % (-1,-1)) # TODO

    def __call__(self, lmfit_params, *fcn_args, **fcn_kws):
        f, self.g, self.sigmaZ = target_and_grad(lmfit_params, *fcn_args, **fcn_kws)
        return f

    def jac(self, lmfit_params, *args, **kwargs):
        return self.g


def model(x, i_shot, Modeler, SIM, compute_grad=True):

    rotX = x["rank%d_shot%d_RotXYZ%d" % (COMM.rank, i_shot, 0)]
    rotY = x["rank%d_shot%d_RotXYZ%d" % (COMM.rank, i_shot, 1)]
    rotZ = x["rank%d_shot%d_RotXYZ%d" % (COMM.rank, i_shot, 2)]
    Na = x["rank%d_shot%d_Nabc%d" % (COMM.rank, i_shot, 0)]
    Nb = x["rank%d_shot%d_Nabc%d" % (COMM.rank, i_shot, 1)]
    Nc = x["rank%d_shot%d_Nabc%d" % (COMM.rank, i_shot, 2)]
    G = x["rank%d_shot%d_Scale" % (COMM.rank, i_shot)]
    num_uc_p = len(Modeler.ucell_man.variables)
    ucell_pars = [x["rank%d_shot%d_Ucell%d" % (COMM.rank, i_shot, i_uc)] for i_uc in range(num_uc_p)]

    # update the Bmatrix
    Modeler.ucell_man.variables = [p.value for p in ucell_pars]
    Bmatrix = Modeler.ucell_man.B_recipspace
    SIM.D.Bmatrix = Bmatrix
    if compute_grad:
        for i_ucell in range(len(ucell_pars)):
            SIM.D.set_ucell_derivative_matrix(
                i_ucell + hopper_utils.UCELL_ID_OFFSET,
                Modeler.ucell_man.derivative_matrices[i_ucell])

    # update the Umat rotation matrix
    SIM.D.set_value(hopper_utils.ROTX_ID, rotX.value)
    SIM.D.set_value(hopper_utils.ROTY_ID, rotY.value)
    SIM.D.set_value(hopper_utils.ROTZ_ID, rotZ.value)

    # update the mosaic block size
    SIM.D.set_ncells_values(tuple([Na.value, Nb.value, Nc.value]))
    npix = int(len(Modeler.pan_fast_slow)/3.)

    # get the forward Bragg scatterint
    SIM.D.add_diffBragg_spots(Modeler.pan_fast_slow)
    bragg_no_scale = SIM.D.raw_pixels_roi[:npix]
    bragg_no_scale = bragg_no_scale.as_numpy_array()

    # apply the per-shot scale factor
    scale = G.value
    bragg = scale*bragg_no_scale

    # add the background
    model_pix = bragg + Modeler.all_background

    # compute the negative log Likelihood
    resid = (Modeler.all_data - model_pix)
    resid_square = resid ** 2
    V = model_pix + Modeler.sigma_rdout ** 2
    neg_LL = (.5*(np.log(2*np.pi*V) + resid_square / V))[Modeler.all_trusted].sum()

    # compute the z-score sigma as a diagnostic
    zscore_sigma = np.std(resid / np.sqrt(V))

    # store the gradients
    J = {}
    if compute_grad:
        # this term is a common factor in all of the gradients
        common_grad_term = (0.5 / V * (1 - 2 * resid - resid_square / V))

        # scale factor gradients
        if G.vary:
            scale_grad = bragg_no_scale
            J[G.name] = (common_grad_term*scale_grad)[Modeler.all_trusted].sum()

        # Umat gradients
        for i_rot, rot in enumerate([rotX, rotY, rotZ]):
            if rot.vary:
                rot_db_id = ROTXYZ_ID[i_rot]
                rot_grad = scale*SIM.D.get_derivative_pixels(rot_db_id).as_numpy_array()[:npix]
                J[rot.name] = (common_grad_term*rot_grad)[Modeler.all_trusted].sum()

        # mosaic block size gradients
        if Na.vary:
            Nabc_grad = SIM.D.get_ncells_derivative_pixels()
            for i_N, N in enumerate([Na, Nb, Nc]):
                N_grad = scale*(Nabc_grad[i_N][:npix].as_numpy_array())
                J[N.name] = (common_grad_term*N_grad)[Modeler.all_trusted].sum()

        # unit cell gradients
        if ucell_pars[0].vary:
            for i_ucell in range(len(ucell_pars)):
                d = scale*SIM.D.get_derivative_pixels(hopper_utils.UCELL_ID_OFFSET+i_ucell).as_numpy_array()[:npix]
                J[ucell_pars[i_ucell].name] = (common_grad_term*d)[Modeler.all_trusted].sum()

        # detector model gradients
        detector_derivs = []
        for diffbragg_parameter_id in PAN_OFS_IDS+PAN_XYZ_IDS:
            d = common_grad_term*scale*(SIM.D.get_derivative_pixels(diffbragg_parameter_id).as_numpy_array()[:npix])
            detector_derivs.append(d)
        names = "RotOrth", "RotFast", "RotSlow", "ShiftX", "ShiftY", "ShiftZ"
        for group_id in Modeler.unique_panel_group_ids:
            for name in names:
                J["group%d_%s" % (group_id, name)] = 0
            for pixel_rng in Modeler.group_id_slices[group_id]:
                trusted_pixels = Modeler.all_trusted[pixel_rng]
                #grad_component = common_grad_term[pixel_rng]
                for i_name, name in enumerate(names):
                    d = detector_derivs[i_name][pixel_rng]
                    #J["group%d_%s" % (group_id, name)] += (grad_component*d)[trusted_pixels].sum()
                    J["group%d_%s" % (group_id, name)] += d[trusted_pixels].sum()

    return neg_LL, J, model_pix, zscore_sigma


def set_group_id_slices(Modeler, group_id_from_panel_id):
    """finds the boundaries for each panel group ID in the 1-D array of per-shot data
    Modeler: DataModeler instance with loaded data
    group_id_from_panel_id : dict where key is panel id and value is group id
    """
    Modeler.all_group_id = [group_id_from_panel_id[pid] for pid in Modeler.all_pid]
    splitter = np.where(np.diff(Modeler.all_group_id) != 0)[0]+1
    npix = len(Modeler.all_data)
    slices = [slice(V[0], V[-1]+1, 1) for V in np.split(np.arange(npix), splitter)]
    group_ids = [V[0] for V in np.split(np.array(Modeler.all_group_id), splitter)]
    group_id_slices = {}
    for i_group, slc in zip(group_ids, slices):
        if i_group not in group_id_slices:
            group_id_slices[i_group] = [slc]
        else:
            group_id_slices[i_group].append(slc)
    Modeler.unique_panel_group_ids = set(Modeler.all_group_id)
    logging.debug("Modeler has data on %d unique panel groups" % (len(Modeler.unique_panel_group_ids)))
    Modeler.group_id_slices = group_id_slices


def update_detector(x, SIM):
    det = SIM.detector
    for pid in range(len(det)):
        group_id = SIM.panel_group_from_id[pid]
        if not group_id in SIM.panel_groups_refined:
            continue
        Oang = x["group%d_RotOrth" % group_id].value
        Fang = x["group%d_RotFast" % group_id].value
        Sang = x["group%d_RotSlow" % group_id].value
        Xdist = x["group%d_ShiftX" % group_id].value
        Ydist = x["group%d_ShiftY" % group_id].value
        Zdist = x["group%d_ShiftZ" % group_id].value

        origin_of_rotation = SIM.panel_reference_from_id[pid]
        SIM.D.reference_origin = origin_of_rotation
        SIM.D.update_dxtbx_geoms(det, SIM.beam.nanoBragg_constructor_beam, pid,
                                  Oang, Fang, Sang, Xdist, Ydist, Zdist,
                                  force=False)


def target_and_grad(x, x_mapping, data_modelers, SIM, compute_grad=True):
    target_functional = 0
    grad = np.zeros(len(x)) if compute_grad else None

    update_detector(x, SIM)

    all_shot_sigZ = []
    for i_shot in data_modelers:
        Modeler = data_modelers[i_shot]

        neg_LL, neg_LL_grad, model_pix, per_shot_sigZ = model(x, i_shot, Modeler, SIM, compute_grad)
        all_shot_sigZ.append(per_shot_sigZ)

        # accumulate the target functional for this rank/shot
        target_functional += neg_LL

        # accumulate the gradients for this rank/shot
        if compute_grad:
            for par_name in neg_LL_grad:
                grad_idx = x_mapping[par_name]
                grad[grad_idx] += neg_LL_grad[par_name]

        # TODO add in the restraints

    # sum the target functional and the gradients across all ranks
    target_functional = COMM.bcast(COMM.reduce(target_functional))
    if compute_grad:
        grad = COMM.bcast(COMM.reduce(grad))

    all_shot_sigZ = COMM.reduce(all_shot_sigZ)
    if COMM.rank==0:
        all_shot_sigZ = np.median(all_shot_sigZ)

    return target_functional, grad, all_shot_sigZ


def geom_min(params):
    import pandas
    launcher = ensemble_refine_launcher.RefineLauncher(params)
    df = pandas.read_pickle(params.geometry.input_pkl)
    launcher.load_inputs(df, refls_key="stage1_refls")

    # same on every rank:
    det_params = DetectorParameters(params, launcher.panel_groups_refined, launcher.n_panel_groups)

    # different on each rank
    crystal_params = CrystalParameters(launcher.Modelers)
    crystal_params.parameters = COMM.bcast(COMM.reduce(crystal_params.parameters))

    LMP = lmfit.Parameters()
    LMP.add_many(*(crystal_params.parameters + det_params.parameters))
    LMP_index_mapping = {name: i for i, name in enumerate(LMP.keys())}

    for i_shot in launcher.Modelers:
        Modeler = launcher.Modelers[i_shot]
        set_group_id_slices(Modeler, launcher.panel_group_from_id)

    # attached some objects to SIM for convenience
    launcher.SIM.panel_reference_from_id = launcher.panel_reference_from_id
    launcher.SIM.panel_group_from_id = launcher.panel_group_from_id
    launcher.SIM.panel_groups_refined = launcher.panel_groups_refined

    # compute gradients, depending on the refinement method
    do_grads = params.geometry.optimize_method == "lbfgsb"
    if not do_grads:
        assert params.geometry.optimize_method == "nelder"

    # configure diffBragg instance for gradient computation
    # TODO: fix flags currently unsupported in lmfit with gradients? One can always "fix" a parameter by
    #       setting the range in DetectorParameters/CrystalParameters to be infinitesimal, e.g. +-1e-10
    if do_grads:
        if not params.fix.RotXYZ:
            for i_rot in range(3):
                launcher.SIM.D.refine(ROTXYZ_ID[i_rot])
        if not params.fix.Nabc:
            launcher.SIM.D.refine(hopper_utils.NCELLS_ID)
        if not params.fix.ucell:
            for i_ucell in range(launcher.SIM.num_ucell_param):
                launcher.SIM.D.refine(hopper_utils.UCELL_ID_OFFSET + i_ucell)
        for i, diffbragg_id in enumerate(PAN_OFS_IDS):
            if not params.geometry.fix.panel_rotations[i]:
                launcher.SIM.D.refine(diffbragg_id)

        for i, diffbragg_id in enumerate(PAN_XYZ_IDS):
            if not params.geometry.fix.panel_translations[i]:
                launcher.SIM.D.refine(diffbragg_id)

    # do a barrel roll!
    target = Target()
    fcn_args = [LMP_index_mapping, launcher.Modelers, launcher.SIM, do_grads]
    fcn_kws = {}
    lbfgs_kws = {}
    if do_grads:
        lbfgs_kws = {"jac": target.jac,
                    "options":  {"ftol": params.ftol, "gtol": 1e-10, "maxfun":1e5, "maxiter":1e5}}

    minzer = lmfit.Minimizer(userfcn=target, params=LMP, fcn_args=fcn_args, fcn_kws=fcn_kws, iter_cb=target.callbk,
                             scale_covar=False, calc_covar=False)
    result = minzer.minimize(method=params.geometry.optimize_method, params=LMP, **lbfgs_kws)

    opt_det = get_optimized_detector(result.params, launcher.SIM)
    from dxtbx.model import Experiment, ExperimentList
    El = ExperimentList()
    E = Experiment()
    E.detector = opt_det
    El.append(E)
    El.as_file(params.geometry.optimized_detector_name)

    print_parmams(result.params)


def get_optimized_detector(x, SIM):
    from dxtbx.model import Detector, Panel
    new_det = Detector()
    for pid in range(len(SIM.detector)):
        panel = SIM.detector[pid]
        panel_dict = panel.to_dict()
        group_id = SIM.panel_group_from_id[pid]
        if group_id in SIM.panel_groups_refined:
            Oang = x["group%d_RotOrth" % group_id].value
            Fang = x["group%d_RotFast" % group_id].value
            Sang = x["group%d_RotSlow" % group_id].value
            Xdist = x["group%d_ShiftX" % group_id].value
            Ydist = x["group%d_ShiftY" % group_id].value
            Zdist = x["group%d_ShiftZ" % group_id].value

            origin_of_rotation = SIM.panel_reference_from_id[pid]
            SIM.D.reference_origin = origin_of_rotation
            SIM.D.update_dxtbx_geoms(SIM.detector, SIM.beam.nanoBragg_constructor_beam, pid,
                                     Oang, Fang, Sang, Xdist, Ydist, Zdist,
                                     force=False)
            fdet = SIM.D.fdet_vector
            sdet = SIM.D.sdet_vector
            origin = SIM.D.get_origin()
        else:
            fdet = panel.get_fast_axis()
            sdet = panel.get_slow_axis()
            origin = panel.get_origin()
        panel_dict["fast_axis"] = fdet
        panel_dict["slow_axis"] = sdet
        panel_dict["origin"] = origin

        new_det.add_panel(Panel.from_dict(panel_dict))

    return new_det


def print_parmams(params):
    for name in params:
        print("%s: %f" % (name, params[name].value))


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("--phil", type=str, required=True, help="path to a phil string")
    parser.add_argument("--cmdlinePhil", nargs="+", default=None, type=str, help="command line phil params")
    progargs = parser.parse_args()

    from libtbx.phil import parse
    from simtbx.diffBragg.phil import philz, hopper_phil

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
    geom_min(params)
