from __future__ import absolute_import, division, print_function
import socket
import glob
from copy import deepcopy
from simtbx.diffBragg import hopper_utils
import h5py
from dxtbx.model.experiment_list import ExperimentList
try:
    import pandas
except ImportError:
    print("Please intsall pandas, libtbx.python -m pip install pandas")
    exit()
from scitbx.matrix import sqr, col

try:
    from line_profiler import LineProfiler
except ImportError:
    LineProfiler = None

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
# TODO, figure out why next 3 lines are sometimes necessary?!
if not hasattr(COMM, "rank"):
    COMM.rank = 0
    COMM.size = 1
from libtbx.phil import parse

from simtbx.diffBragg import utils
from simtbx.diffBragg.phil import philz

import logging
from simtbx.diffBragg import mpi_logger
from simtbx.diffBragg.phil import hopper_phil


philz = hopper_phil + philz
phil_scope = parse(philz)


def make_rank_outdir(root, subfolder):
    rank_imgs_outdir = os.path.join(root, subfolder, "rank%d" % COMM.rank)
    if not os.path.exists(rank_imgs_outdir):
        os.makedirs(rank_imgs_outdir)
    return rank_imgs_outdir


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
            assert self.params.outdir is not None
        self.params = COMM.bcast(self.params)
        if COMM.rank == 0:
            if not os.path.exists(self.params.outdir):
                utils.safe_makedirs(self.params.outdir)
        COMM.barrier()

        if self.params.logging.logname is None:
            self.params.logging.logname = "main_stage1.log"
        if self.params.profile_name is None:
            self.params.profile_name = "prof_stage1.log"
        mpi_logger.setup_logging_from_params(self.params)

    def run(self):
        assert os.path.exists(self.params.exp_ref_spec_file)
        input_lines = None
        best_models = None
        if COMM.rank == 0:
            input_lines = open(self.params.exp_ref_spec_file, "r").readlines()
            if self.params.skip is not None:
                input_lines = input_lines[self.params.skip:]
            if self.params.first_n is not None:
                input_lines = input_lines[:self.params.first_n]
            if self.params.sanity_test_input:
                hopper_utils.sanity_test_input_lines(input_lines)

            if self.params.best_pickle is not None:
                logging.info("reading pickle %s" % self.params.best_pickle)
                best_models = pandas.read_pickle(self.params.best_pickle)

            if self.params.dump_gathers:
                if self.params.gathers_dir is None:
                    raise ValueError("Need to provide a file dir path in order to dump_gathers")
                utils.safe_makedirs(self.params.gathers_dir)
        input_lines = COMM.bcast(input_lines)
        best_models = COMM.bcast(best_models)

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

            logging.info("COMM.rank %d on shot  %d / %d" % (COMM.rank, i_exp + 1, len(input_lines)))
            line_fields = line.strip().split()
            assert len(line_fields) in [2, 3]
            if len(line_fields) == 2:
                exp, ref = line_fields
                spec = None
            else:
                exp, ref, spec = line_fields

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
            Modeler = hopper_utils.DataModeler(self.params)
            if self.params.load_data_from_refls:
                gathered = Modeler.GatherFromReflectionTable(exp, ref)
            else:
                gathered = Modeler.GatherFromExperiment(exp, ref)
            if not gathered:
                logging.warning("No refls in %s; CONTINUE; COMM.rank=%d" % (ref, COMM.rank))
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
                    new_Modeler = hopper_utils.DataModeler(self.params)
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
            if self.params.refiner.verbose is not None and COMM.rank==0:
                Modeler.SIM.D.verbose = self.params.refiner.verbose
            if self.params.profile:
                Modeler.SIM.record_timings = True
            if self.params.use_float32:
                Modeler.all_data = Modeler.all_data.astype(np.float32)
                Modeler.all_background = Modeler.all_background.astype(np.float32)

            if self.params.refiner.randomize_devices:
                dev = np.random.choice(self.params.refiner.num_devices)
                logging.info("Rank %d will use randomly chosen device %d on host %s" % (COMM.rank, dev, socket.gethostname()))
            else:
                dev = COMM.rank % self.params.refiner.num_devices
                logging.info("Rank %d will use device %d on host %s" % (COMM.rank, dev, socket.gethostname()))

            Modeler.SIM.D.device_Id = dev

            nparam = len(Modeler.SIM.P)
            x0 = [1] * nparam
            x = Modeler.Minimize(x0)
            if self.params.profile:
                Modeler.SIM.D.show_timings(COMM.rank) #, out)
            save_up(Modeler, x, exp, i_exp, ref)

        if self.params.dump_gathers and self.params.gathered_output_file is not None:
            exp_gatheredRef_spec = COMM.reduce(exp_gatheredRef_spec)
            if COMM.rank == 0:
                o = open(self.params.gathered_output_file, "w")
                for e, r, s in exp_gatheredRef_spec:
                    if s is not None:
                        o.write("%s %s %s\n" % (e,r,s))
                    else:
                        o.write("%s %s\n" % (e,r))
                o.close()


def save_up(Modeler, x, exp, i_exp, input_refls):
    LOGGER = logging.getLogger("refine")
    Modeler.best_model, _ = hopper_utils.model(x, Modeler.SIM, Modeler.pan_fast_slow, compute_grad=False)
    LOGGER.info("Optimized values for i_exp %d:" % i_exp)
    hopper_utils.look_at_x(x, Modeler.SIM)

    rank_imgs_outdir = make_rank_outdir(Modeler.params.outdir, "imgs")
    rank_SIMlog_outdir = make_rank_outdir(Modeler.params.outdir, "simulator_state")
    rank_refls_outdir = make_rank_outdir(Modeler.params.outdir, "refls")
    rank_spectra_outdir = make_rank_outdir(Modeler.params.outdir, "spectra")

    basename = os.path.splitext(os.path.basename(exp))[0]
    img_path = os.path.join(rank_imgs_outdir, "%s_%s_%d.h5" % (Modeler.params.tag, basename, i_exp))
    SIMlog_path = os.path.join(rank_SIMlog_outdir, "%s_%s_%d.txt" % (Modeler.params.tag, basename, i_exp))
    new_refls_file = os.path.join(rank_refls_outdir, "%s_%s_%d.refl" % (Modeler.params.tag, basename, i_exp))
    spectra_path = os.path.join(rank_spectra_outdir, "%s_%s_%d.lam" % (Modeler.params.tag, basename, i_exp))

    if Modeler.SIM.num_xtals == 1:
        save_to_pandas(x, Modeler.SIM, exp, Modeler.params, Modeler.E, i_exp, input_refls, img_path)

    data_subimg, model_subimg, trusted_subimg, bragg_subimg = Modeler.get_data_model_pairs()

    comp = {"compression": "lzf"}
    new_refls = deepcopy(Modeler.refls)
    has_xyzcal = 'xyzcal.px' in list(new_refls.keys())
    if has_xyzcal:
        new_refls['dials.xyzcal.px'] = deepcopy(new_refls['xyzcal.px'])
    new_xycalcs = flex.vec3_double(len(Modeler.refls), (np.nan,np.nan,np.nan))
    h5_roi_id = flex.int(len(Modeler.refls), -1)
    with h5py.File(img_path, "w") as h5:
        sigmaZs = []
        for i_roi in range(len(data_subimg)):
            dat = data_subimg[i_roi]
            fit = model_subimg[i_roi]
            trust = trusted_subimg[i_roi]
            sig = np.sqrt(fit + Modeler.sigma_rdout**2)
            Z = (dat - fit) / sig
            sigmaZ = Z[trust].std()
            sigmaZs.append(sigmaZ)
            h5.create_dataset("data/roi%d" % i_roi, data=data_subimg[i_roi], **comp)
            h5.create_dataset("model/roi%d" % i_roi, data=model_subimg[i_roi], **comp)
            if bragg_subimg[0] is not None:
                h5.create_dataset("bragg/roi%d" % i_roi, data=bragg_subimg[i_roi], **comp)
                if np.any(bragg_subimg[i_roi] > 0):
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
        h5.create_dataset("sigmaZ_vals", data=sigmaZs)
        if Modeler.Hi_asu is not None:
            h5.create_dataset("Hi_asu", data=Modeler.Hi_asu)

    new_refls["xyzcal.px"] = new_xycalcs
    new_refls["h5_roi_idx"] = h5_roi_id
    if Modeler.params.filter_unpredicted_refls_in_output:
        sel = [not np.isnan(x) for x,y,z in new_xycalcs]
        new_refls = new_refls.select(flex.bool(sel))
    new_refls.as_file(new_refls_file)

    if Modeler.params.refiner.debug_pixel_panelfastslow is not None:
        # TODO separate diffBragg logger
        utils.show_diffBragg_state(Modeler.SIM.D, Modeler.params.refiner.debug_pixel_panelfastslow)
        print("Refiner scale=%f" % Modeler.SIM.Scale_params[0].get_val(x[0]))

    hopper_utils.finalize_SIM(Modeler.SIM, log=SIMlog_path, lam=spectra_path)


def save_to_pandas(x, SIM, orig_exp_name, params, expt, rank_exp_idx, stg1_refls, stg1_img_path):
    LOGGER = logging.getLogger("refine")
    rank_exper_outdir = make_rank_outdir(params.outdir, "expers")
    rank_pandas_outdir =make_rank_outdir(params.outdir, "pandas")

    #if SIM.num_xtals > 1:
    #    raise NotImplemented("cant save pandas for multiple crystals yet")
    scale, rotX, rotY, rotZ, Na, Nb, Nc,diff_gam_a, diff_gam_b, diff_gam_c, diff_sig_a, diff_sig_b, diff_sig_c, a,b,c,al,be,ga,detz_shift = hopper_utils.get_param_from_x(x, SIM)
    if params.isotropic.diffuse_gamma:
        diff_gam_b = diff_gam_c = diff_gam_a
    if params.isotropic.diffuse_sigma:
        diff_sig_b = diff_sig_c = diff_sig_a
    shift = np.nan
    #if SIM.shift_param is not None:
    #    shift = SIM.shift_param.get_val(x[-1])
    xtal_scales = [scale]
    eta_a, eta_b, eta_c = hopper_utils.get_mosaicity_from_x(x, SIM)
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
        "diffuse_gamma": [(diff_gam_a, diff_gam_b, diff_gam_c)],
        "diffuse_sigma": [(diff_sig_a, diff_sig_b, diff_sig_c)],
        "fp_fdp_shift": [shift],
        "use_diffuse_models": [params.use_diffuse_models],
        "gamma_miller_units": [params.gamma_miller_units],
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
    LOGGER.debug("saved opt_exp %s with wavelength %f" % (opt_exp_path, expt.beam.get_wavelength()))

    spec_file = None
    if params.simulator.spectrum.filename is not None:
        spec_file = os.path.abspath(params.simulator.spectrum.filename)
    df["spectrum_filename"] = spec_file
    df["spectrum_stride"] = params.simulator.spectrum.stride

    df["total_flux"] = SIM.D.flux  # params.simulator.total_flux
    df["beamsize_mm"] = SIM.beam.size_mm
    df["exp_name"] = os.path.abspath(orig_exp_name)
    df["opt_exp_name"] = os.path.abspath(opt_exp_path)
    df["spectrum_from_imageset"] = params.spectrum_from_imageset
    df["oversample"] = params.simulator.oversample
    if params.opt_det is not None:
        df["opt_det"] = params.opt_det
    df["stage1_refls"] = stg1_refls
    df["stage1_output_img"] = stg1_img_path

    df.to_pickle(pandas_path)
    return df


if __name__ == '__main__':
    from dials.util import show_mail_on_error

    with show_mail_on_error():
        script = Script()
        RUN = script.run
        lp = None
        if LineProfiler is not None and script.params.profile:
            lp = LineProfiler()
            lp.add_function(hopper_utils.model)
            lp.add_function(hopper_utils.target_func)
            RUN = lp(script.run)
        elif script.params.profile:
            print("Install line_profiler in order to use logging: libtbx.python -m pip install line_profiler")

        RUN()

        if lp is not None:
            stats = lp.get_stats()
            hopper_utils.print_profile(stats, ["model", "target_func"])
