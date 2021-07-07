from __future__ import absolute_import, division, print_function
import socket
import glob
from copy import deepcopy
from simtbx.diffBragg.hopper_utils import look_at_x, model, get_param_from_x, DataModeler, get_data_model_pairs
import h5py
from dxtbx.model.experiment_list import ExperimentList
try:
    import pandas
except ImportError:
    print("Please intsall pandas, libtbx.python -m pip install pandas")
    exit()
from scitbx.matrix import sqr, col

ROTX_ID = 0
ROTY_ID = 1
ROTZ_ID = 2
NCELLS_ID = 9
UCELL_ID_OFFSET = 3
DETZ_ID = 10


# LIBTBX_SET_DISPATCHER_NAME simtbx.diffBragg.hopper


import numpy as np
np.seterr(invalid='ignore')
import os
from dials.array_family import flex
from libtbx.mpi4py import MPI

COMM = MPI.COMM_WORLD
from libtbx.phil import parse

from simtbx.diffBragg import utils
from simtbx.diffBragg.phil import philz

hopper_phil = """
use_float32 = False
  .type = bool
  .help = store pixel data and background models in 32bit arrays
test_gathered_file = False
  .type = bool
  .help = run a quick test to ensure the gathered data file preserves information
load_data_from_refls = False
  .type = bool
  .help = load image data, background etc from reflection tables
gathered_output_file = None
  .type = str
  .help = optional file for storing a new hopper input file which points to the gathered data dumps
only_dump_gathers = False
  .type = bool
  .help = only reads in image data, fits background planes, and dumps
  .help = results to disk, writes a new exper refl file at the end
gathers_dir = None
  .type = str
  .help = folder where gathered data reflection tables
  .help = will be writen (if dump_gathers=True)
dump_gathers = False
  .type = bool
  .help = optionally dump the loaded experimental data to reflection tables
  .help = for portability
spectrum_from_imageset = False
  .type = bool
  .help = if True, load the spectrum from the imageset in the experiment, then probably downsample it
downsamp_spec {
  filt_freq = 0.07
    .type = float
    .help = low pass filter frequency in units of inverse spectrometer pixels (??)
  filt_order = 3
    .type = int
    .help = order for bandpass butter filter
  tail = 50
    .type = int
    .help = endpoints of the spectrum that are used in background estimation
  delta_en = 0.5
    .type = float
    .help = final resolution of downsampled spectrum in eV
}
apply_best_crystal_model = False
  .type = bool
  .help = depending on what experiments in the exper refl file, one may want
  .help = to apply the optimal crystal transformations (this parameter only matters
  .help = if params.best_pickle is not None)
filter_unpredicted_refls_in_output = True
  .type = bool
  .help = filter reflections in the output refl table for which there was no model bragg peak
  .help = after stage 1 termination
tag = simplex
  .type = str
  .help = output name tag
ignore_existing = False
  .type = bool
  .help = experimental, ignore expts that already have optimized models in the output dir
global_method = *basinhopping annealing
  .type = choice
  .help = the method of global optimization to use
nelder_mead_maxfev = 60
  .type = int
  .help = multiplied by total number of modeled pixels to get max number of iterations
niter_per_J = 3
  .type = int
  .help = if using gradient descent, compute gradients
  .help = every niter_per_J iterations .
rescale_params = True
  .type = bool
  .help = use rescaled range parameters
use_likelihood_target = False
  .type = bool
  .help = if True, then use negative log Likelihood derived from a gaussian noise model
  .help = as opposed to the least squares target equation
best_pickle = None
  .type = str
  .help = path to a pandas pickle containing the best models for the experiments
betas {
  detz_shift = 10
    .type = float
    .help = restraint variance for detector shift target
  ucell = [0,0,0,0,0,0]
    .type = floats
    .help = beta values for unit cell constants
  RotXYZ = 0
    .type = float
    .help = restraint factor for the rotXYZ restraint
  Nabc = [0,0,0]
    .type = floats(size=3)
    .help = restraint factor for the ncells abc
  G = 0
    .type = float
    .help = restraint factor for the scale G
  B = 0
    .type = float
    .help = restraint factor for Bfactor
}
dual {
  initial_temp = 5230
    .type = float
    .help = init temp for dual annealing
  no_local_search = False
    .type = bool
    .help = whether to try local search procedure with dual annealing
    .help = if False, then falls back on classical simulated annealing
  visit = 2.62
    .type = float
    .help = dual_annealing visit param, see scipy optimize docs
  accept = -5
    .type = float
    .help = dual_annealing accept param, see scipy optimize docs
}
centers {
  detz_shift = 0
    .type = float
    .help = restraint target for detector shift along z-direction
  ucell = [63.66, 28.87, 35.86, 1.8425]
    .type = floats
    .help = centers for unit cell constants
  RotXYZ = [0,0,0]
    .type = floats(size=3)
    .help = restraint target for Umat rotations
  Nabc = [100,100,100]
    .type = floats(size=3)
    .help = restraint target for Nabc
  G = 100
    .type = float
    .help = restraint target for scale G
  B = 0
    .type = float
    .help = restraint target for Bfactor
}
levmar {
  damper = 1e-5
    .type = float
    .help = damping coefficient
  maxiter = 100
    .type = int
    .help = max iterations
  up = 10
    .type = float
    .help = scale-up factor
  down = 10
    .type = float
    .help = scale-down factor
  eps4 = 1e-3
    .type = float
    .help = metric improvement threshold for accepting parameter shift
}
skip = None
  .type = int
  .help = skip this many exp
hess = None
  .type = str
  .help = scipy minimize hessian argument, 2-point, 3-point, cs, or None
stepsize = 0.5
  .type = float
  .help = basinhopping stepsize
temp = 1
  .type = float
  .help = temperature for basin hopping algo
niter = 100
  .type = int
  .help = number of basin hopping iters
lsq = True
  .type = bool
  .help = minimizes least squares, if False, minimizes likelihood
weights = True
  .type = bool
  .help = use weights in the target function
exp_ref_spec_file = None
  .type = str
  .help = path to 3 col txt file containing file names for exper, refl, spectrum (.lam)
method = None
  .type = str
  .help = minimizer method
opt_det = None
  .type = str
  .help = path to experiment with optimized detector model
number_of_xtals = 1
  .type = int
  .help = number of crystal domains to model per shot
sanity_test_input = True
  .type = bool
  .help = sanity test input
outdir = True
  .type = str
  .help = output folder
quiet = False
  .type = bool
  .help = silence most output
max_process = -1
  .type = int
  .help = max exp to process
atols = [0.0001, 0.0001]
  .type = floats(size=2)
  .help = atol and fatol to be passed to nelder mead minimizer (termination params)
plot_at_end = False
  .type = bool
  .help = plot subimgs at end of minimize
embed_at_end = False
  .type = bool
  .help = embed to ipython at end of minimize
sigmas {
  detz_shift = 1
    .type = float
    .help = sensitivity shift for the overall detector shift along z-direction
  Nabc = [1,1,1]
    .type = floats(size=3)
    .help = sensitivity for Nabc
  Ndef = [1,1,1]
    .type = floats(size=3)
    .help = sensitivity for Ndef
  RotXYZ = [1,1,1]
    .type = floats(size=3)
    .help = sensitivity for RotXYZ
  G = 1
    .type = float
    .help = sensitivity for scale factor
  B = 1
    .type = float
    .help = sensitivity for Bfactor
  ucell = [1,1,1,1,1,1]
    .type = floats
    .help = sensitivity for unit cell params
  Fhkl = 1
    .type = float
    .help = sensitivity for structure factors
}
init {
  detz_shift = 0
    .type = float
    .help = initial value for the detector position overall shift along z-direction in millimeters
  Nabc = [100,100,100]
    .type = floats(size=3)
    .help = init for Nabc
  Ndef = [0,0,0]
    .type = floats(size=3)
    .help = init for Ndef
  RotXYZ = [0,0,0]
    .type = floats(size=3)
    .help = init for RotXYZ
  G = 1
    .type = float
    .help = init for scale factor
  B = 0
    .type = float
    .help = init for B factor
}
mins {
  detz_shift = -10
    .type = float
    .help = min value for detector z-shift in millimeters
  Nabc = [3,3,3]
    .type = floats(size=3)
    .help = min for Nabc
  Ndef = [-200,-200,-200]
    .type = floats(size=3)
    .help = min for Ndef
  RotXYZ = [-1,-1,-1]
    .type = floats(size=3)
    .help = min for rotXYZ in degrees
  G = 0
    .type = float
    .help = min for scale G
  B = 0
    .type = float
    .help = min for Bfactor
  Fhkl = 0
    .type = float
    .help = min for structure factors
}
maxs {
  detz_shift = 10
    .type = float
    .help = max value for detector z-shift in millimeters
  eta = 0.1
    .type = float
    .help = maximum mosaic spread in degrees
  Nabc = [300,300,300]
    .type = floats(size=3)
    .help = max for Nabc
  Ndef = [200,200,200]
    .type = floats(size=3)
    .help = max for Ndef
  RotXYZ = [1,1,1]
    .type = floats(size=3)
    .help = max for rotXYZ in degrees
  G = 1e12
    .type = float
    .help = max for scale G
  B = 1e3
    .type = float
    .help = max for Bfactor
  Fhkl = 1e6
    .type = float
    .help = max for structure factors
}
fix {
  G = False
    .type = bool
    .help = fix the Bragg spot scale during refinement
  RotXYZ = False
    .type = bool
    .help = fix the misorientation matrix during refinement
  Nabc = False
    .type = bool
    .help = fix the diagonal mosaic domain size parameters during refinement
  ucell = False
    .type = bool
    .help = fix the unit cell during refinement
  detz_shift = False
    .type = bool
    .help = fix the detector distance shift during refinement
}
relative_tilt = True
  .type = bool
  .help = fit tilt coef relative to roi corner
num_mosaic_blocks = 1
  .type = int
  .help = number of mosaic blocks making up mosaic spread dist (not implemented)
ucell_edge_perc = 10 
  .type = float
  .help = precentage for allowing ucell to fluctuate during refinement
ucell_ang_abs = 5
  .type = float
  .help = absolute angle deviation in degrees for unit cell angles to vary during refinement
no_Nabc_scale = False
  .type = bool
  .help = toggle Nabc scaling of the intensity
"""

philz = hopper_phil + philz
phil_scope = parse(philz)


class Script:
    def __init__(self):
        from dials.util.options import OptionParser

        self.params = self.parser = None
        if COMM.rank == 0:
            self.parser = OptionParser(
                usage="",  # stage 1 (per-shot) diffBragg refinement",
                sort_options=True,
                phil=phil_scope,
                read_experiments=True,
                read_reflections=True,
                check_format=False,
                epilog="PyCuties")
        self.parser = COMM.bcast(self.parser)
        if COMM.rank == 0:
            self.params, _ = self.parser.parse_args(show_diff_phil=True)
        self.params = COMM.bcast(self.params)

    def run(self):
        assert os.path.exists(self.params.exp_ref_spec_file)
        input_lines = None
        best_models = None
        if COMM.rank == 0:
            input_lines = open(self.params.exp_ref_spec_file, "r").readlines()
            if self.params.sanity_test_input:
                for line in input_lines:
                    for fname in line.strip().split():
                        if not os.path.exists(fname):
                            raise FileNotFoundError("File %s not there " % fname)
            if self.params.best_pickle is not None:
                if not self.params.quiet: print("reading pickle %s" % self.params.best_pickle)
                best_models = pandas.read_pickle(self.params.best_pickle)

            if self.params.dump_gathers:
                if self.params.gathers_dir is None:
                    raise ValueError("Need to provide a file dir path in order to dump_gathers")
                utils.safe_makedirs(self.params.gathers_dir)
        input_lines = COMM.bcast(input_lines)
        best_models = COMM.bcast(best_models)

        input_lines = input_lines[self.params.skip:]
        if self.params.ignore_existing:
            exp_names_already =None
            if COMM.rank==0:
                exp_names_already = {os.path.basename(f) for f in glob.glob("%s/expers/rank*/*.expt" % self.params.outdir)}
            exp_names_already = COMM.bcast(exp_names_already)

        exp_gatheredRef_spec = []  # optional list of expt, refls, spectra
        for i_exp, line in enumerate(input_lines):
            if i_exp == self.params.max_process:
                break
            if i_exp % COMM.size != COMM.rank:
                continue

            print("COMM.rank %d on shot  %d / %d" % (COMM.rank, i_exp + 1, len(input_lines)))
            exp, ref, spec = line.strip().split()

            if self.params.ignore_existing:
                basename = os.path.splitext(os.path.basename(exp))[0]
                opt_exp = "%s_%s_%d.expt" % (self.params.tag, basename, i_exp)
                if opt_exp in exp_names_already:
                    continue

            best = None
            if best_models is not None:
                best = best_models.query("exp_name=='%s'" % exp)
                if len(best) == 0:
                    best = best_models.query("opt_exp_name=='%s'" % exp)
                if len(best) != 1:
                    raise ValueError("Should be 1 entry for exp %s in best pickle %s" % (exp, self.params.best_pickle))
            self.params.simulator.spectrum.filename = spec
            Modeler = DataModeler(self.params)
            if self.params.load_data_from_refls:
                gathered = Modeler.GatherFromReflectionTable(exp, ref)
            else:
                gathered = Modeler.GatherFromExperiment(exp, ref)
            if not gathered:
                print("No refls in %s; CONTINUE; COMM.rank=%d" % (ref, COMM.rank))
                continue
            if self.params.dump_gathers:
                output_name = os.path.splitext(os.path.basename(exp))[0]
                output_name += "_withData.refl"
                output_name = os.path.join(self.params.gathers_dir, output_name)
                Modeler.dump_gathered_to_refl(output_name, do_xyobs_sanity_check=True)  # NOTE do this is modelin strong spots only
                if self.params.test_gathered_file:
                    all_data = Modeler.all_data.copy()
                    all_roi_id = Modeler.roi_id.copy()
                    all_bg = Modeler.all_background.copy()
                    all_trusted = Modeler.all_trusted.copy()
                    all_pids = np.array(Modeler.pids)
                    all_rois = np.array(Modeler.rois)
                    new_Modeler = DataModeler(self.params)
                    assert new_Modeler.GatherFromReflectionTable(exp, output_name)
                    assert np.allclose(new_Modeler.all_data, all_data)
                    assert np.allclose(new_Modeler.all_background, all_bg)
                    assert np.allclose(new_Modeler.rois, all_rois)
                    assert np.allclose(new_Modeler.pids, all_pids)
                    assert np.allclose(new_Modeler.all_trusted, all_trusted)
                    assert np.allclose(new_Modeler.roi_id, all_roi_id)

                exp_gatheredRef_spec.append((exp, os.path.abspath(output_name), spec))
                if self.params.only_dump_gathers:
                    continue

            Modeler.SimulatorFromExperiment(best)
            if self.params.use_float32:
                Modeler.all_data = Modeler.all_data.astype(np.float32)
                Modeler.all_background = Modeler.all_background.astype(np.float32)

            if self.params.refiner.randomize_devices:
                dev = np.random.choice(self.params.refiner.num_devices)
                print("Rank %d will use random device %d on host %s" % (COMM.rank, dev, socket.gethostname()), flush=True)
            else:
                dev = COMM.rank % self.params.refiner.num_devices
                print("Rank %d will use fixed device %d on host %s" % (COMM.rank, dev, socket.gethostname()), flush=True)

            Modeler.SIM.D.device_Id = dev

            # initial parameters (all set to 1, 7 parameters (scale, rotXYZ, Ncells_abc) per crystal (sausage) and then the unit cell parameters
            nparam = 7 * Modeler.SIM.num_xtals + len(Modeler.SIM.ucell_man.variables) + 1
            if self.params.rescale_params:
                x0 = [1] * nparam
            else:
                x0 = [np.nan]*nparam
                for i_xtal in range(Modeler.SIM.num_xtals):
                    x0[7*i_xtal] = Modeler.SIM.Scale_params[i_xtal].init
                    x0[7*i_xtal+1] = Modeler.SIM.RotXYZ_params[3*i_xtal].init
                    x0[7*i_xtal+2] = Modeler.SIM.RotXYZ_params[3*i_xtal+1].init
                    x0[7*i_xtal+3] = Modeler.SIM.RotXYZ_params[3*i_xtal+2].init
                    x0[7*i_xtal+4] = Modeler.SIM.Nabc_params[3*i_xtal].init
                    x0[7*i_xtal+5] = Modeler.SIM.Nabc_params[3*i_xtal+1].init
                    x0[7*i_xtal+6] = Modeler.SIM.Nabc_params[3*i_xtal+2].init

                nucell = len(Modeler.SIM.ucell_man.variables)
                for i_ucell in range(nucell):
                    x0[7*Modeler.SIM.num_xtals+i_ucell] = Modeler.SIM.ucell_params[i_ucell].init
                    #x0[7 * Modeler.SIM.num_xtals + i_ucell] = Modeler.SIM.ucell_params[i_ucell].init
                x0[7*Modeler.SIM.num_xtals+nucell] = Modeler.SIM.DetZ_param.init

            x = Modeler.Minimize(x0)
            save_up(Modeler, x, exp, i_exp, ref)

        if self.params.dump_gathers and self.params.gathered_output_file is not None:
            exp_gatheredRef_spec = COMM.reduce(exp_gatheredRef_spec)
            if COMM.rank==0:
                o = open(self.params.gathered_output_file, "w")
                for e,r,s in exp_gatheredRef_spec:
                    o.write("%s %s %s\n" % (e,r,s))
                o.close()


def save_up(Modeler, x, exp, i_exp, input_refls):
    best_model,_ = model(x, Modeler.SIM, Modeler.pan_fast_slow, compute_grad=False)
    print("Optimized:")
    look_at_x(x,Modeler.SIM)

    rank_imgs_outdir = os.path.join(Modeler.params.outdir, "imgs", "rank%d" % COMM.rank)
    if not os.path.exists(rank_imgs_outdir):
        os.makedirs(rank_imgs_outdir)

    rank_refls_outdir = os.path.join(Modeler.params.outdir, "refls", "rank%d" % COMM.rank)
    if not os.path.exists(rank_refls_outdir):
        os.makedirs(rank_refls_outdir)

    basename = os.path.splitext(os.path.basename(exp))[0]

    img_path = os.path.join(rank_imgs_outdir, "%s_%s_%d.h5" % (Modeler.params.tag, basename, i_exp))

    if Modeler.SIM.num_xtals == 1:
        save_to_pandas(x, Modeler.SIM, exp, Modeler.params, Modeler.E, i_exp, input_refls, img_path)

    new_refls_file = os.path.join(rank_refls_outdir, "%s_%s_%d.refl" % (Modeler.params.tag, basename, i_exp))
    # save_model_Z(img_path, all_data, best_model, pan_fast_slow, sigma_rdout)

    data_subimg, model_subimg, strong_subimg, bragg_subimg = get_data_model_pairs(Modeler.rois, Modeler.pids, Modeler.roi_id, best_model, Modeler.all_data, background=Modeler.all_background)

    comp = {"compression": "lzf"}
    new_refls = deepcopy(Modeler.refls)
    new_refls['dials.xyzcal.px'] = deepcopy(new_refls['xyzcal.px'])
    new_xycalcs = flex.vec3_double(len(Modeler.refls), (0,0,0))
    h5_roi_id = flex.int(len(Modeler.refls), -1)
    with h5py.File(img_path, "w") as h5:
        for i_roi in range(len(data_subimg)):
            h5.create_dataset("data/roi%d" % i_roi, data=data_subimg[i_roi], **comp)
            h5.create_dataset("model/roi%d" % i_roi, data=model_subimg[i_roi], **comp)
            if bragg_subimg[0] is not None:
                h5.create_dataset("bragg/roi%d" % i_roi, data=bragg_subimg[i_roi], **comp)
                com = np.nan, np.nan, np.nan
                if np.any(bragg_subimg[i_roi]>0):
                    I = bragg_subimg[i_roi]
                    Y,X = np.indices(bragg_subimg[i_roi].shape)
                    x1,_,y1,_ = Modeler.rois[i_roi]
                    X += x1
                    Y += y1
                    Isum = I.sum()
                    xcom = (X*I).sum() / Isum
                    ycom = (Y*I).sum() / Isum
                    com = xcom+.5, ycom+.5, 0

                ref_idx = Modeler.refls_idx[i_roi]
                h5_roi_id[ref_idx] = i_roi
                new_xycalcs[ref_idx] = com


        h5.create_dataset("rois", data=Modeler.rois)
        h5.create_dataset("pids", data=Modeler.pids)
        h5.create_dataset("sigma_rdout", data=Modeler.sigma_rdout)

    new_refls["xyzcal.px"] = new_xycalcs
    new_refls["h5_roi_idx"] = h5_roi_id
    if Modeler.params.filter_unpredicted_refls_in_output:
        sel = [not np.isnan(x) for x,y,z in new_xycalcs]
        new_refls = new_refls.select(flex.bool(sel))
    new_refls.as_file(new_refls_file)

    if Modeler.params.embed_at_end:
        from IPython import embed
        embed()

    Modeler.SIM.D.free_all()
    Modeler.SIM.D.free_Fhkl2()
    Modeler.SIM.D.gpu_free()


def save_to_pandas(x, SIM, orig_exp_name, params, expt, rank_exp_idx, stg1_refls, stg1_img_path):
    rank_exper_outdir = os.path.join(params.outdir, "expers", "rank%d" % COMM.rank)
    rank_pandas_outdir = os.path.join(params.outdir, "pandas", "rank%d" % COMM.rank)
    for d in [rank_exper_outdir, rank_pandas_outdir]:
        if not os.path.exists(d):
            os.makedirs(d)

    if SIM.num_xtals > 1:
        raise NotImplemented("cant save pandas for multiple crystals yet")
    scale, rotX, rotY, rotZ, Na, Nb, Nc,a,b,c,al,be,ga,detz_shift = get_param_from_x(x, SIM)
    shift = np.nan
    #if SIM.shift_param is not None:
    #    shift = SIM.shift_param.get_val(x[-1])
    xtal_scales = [scale]
    eta_a = eta_b = eta_c = 0
    a_init, b_init, c_init, al_init, be_init, ga_init = SIM.crystal.dxtbx_crystal.get_unit_cell().parameters()

    xax = col((-1, 0, 0))
    yax = col((0, -1, 0))
    zax = col((0, 0, -1))
    ## update parameters:
    RX = xax.axis_and_angle_as_r3_rotation_matrix(rotX, deg=False)
    RY = yax.axis_and_angle_as_r3_rotation_matrix(rotY, deg=False)
    RZ = zax.axis_and_angle_as_r3_rotation_matrix(rotZ, deg=False)
    M = RX * RY * RZ
    U = M * sqr(SIM.crystal.dxtbx_crystal.get_U())
    SIM.crystal.dxtbx_crystal.set_U(U)

    ucparam = a, b, c, al, be, ga
    ucman = utils.manager_from_params(ucparam)
    SIM.crystal.dxtbx_crystal.set_B(ucman.B_recipspace)

    Amats = [SIM.crystal.dxtbx_crystal.get_A()]
    ncells_def_vals = [(0, 0, 0)]
    ncells_vals = [(Na, Nb, Nc)]
    eta = [0]
    lam0 = [-1]
    lam1 = [-1]
    df = pandas.DataFrame({
        # "panX": list(panX), "panY": list(panY), "panZ": list(panZ),
        # "panO": list(panO), "panF": list(panF), "panS": list(panS),
        "spot_scales": xtal_scales, "Amats": Amats, "ncells": ncells_vals,
        "eta_abc": [(eta_a, eta_b, eta_c)],
        "detz_shift_mm": [detz_shift * 1e3],
        "ncells_def": ncells_def_vals,
        "fp_fdp_shift": [shift],
        # "bgplanes": bgplanes, "image_corr": image_corr,
        # "init_image_corr": init_img_corr,
        # "fcell_xstart": fcell_xstart,
        # "ucell_xstart": ucell_xstart,
        # "init_misorient": init_misori, "final_misorient": final_misori,
        # "bg_coef": bg_coef,
        "eta": eta,
        "rotX": rotX,
        "rotY": rotY,
        "rotZ": rotZ,
        "a": a, "b": b, "c": c, "al": al, "be": be, "ga": ga,
        "a_init": a_init, "b_init": b_init, "c_init": c_init, "al_init": al_init,
        "lam0": lam0, "lam1": lam1,
        "be_init": be_init, "ga_init": ga_init})
    # "scale_xpos": scale_xpos,
    # "ncells_xpos": ncells_xstart,
    # "bgplanes_xpos": bgplane_xpos})

    basename = os.path.splitext(os.path.basename(orig_exp_name))[0]
    opt_exp_path = os.path.join(rank_exper_outdir, "%s_%s_%d.expt" % (params.tag, basename, rank_exp_idx))
    pandas_path = os.path.join(rank_pandas_outdir, "%s_%s_%d.pkl" % (params.tag, basename, rank_exp_idx))
    expt.crystal = SIM.crystal.dxtbx_crystal
    # expt.detector = refiner.get_optimized_detector()
    new_exp_list = ExperimentList()
    new_exp_list.append(expt)
    new_exp_list.as_file(opt_exp_path)

    df["spectrum_filename"] = os.path.abspath(params.simulator.spectrum.filename)
    df["spectrum_stride"] = params.simulator.spectrum.stride
    df["total_flux"] = params.simulator.total_flux
    df["beamsize_mm"] = SIM.beam.size_mm
    df["exp_name"] = os.path.abspath(orig_exp_name)
    df["opt_exp_name"] = os.path.abspath(opt_exp_path)
    df["oversample"] = params.simulator.oversample
    if params.opt_det is not None:
        df["opt_det"] = params.opt_det
    df["stage1_refls"] = stg1_refls
    df["stage1_output_img"] = stg1_img_path

    df.to_pickle(pandas_path)


if __name__ == '__main__':
    from dials.util import show_mail_on_error

    with show_mail_on_error():
        script = Script()
        script.run()
