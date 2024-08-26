from __future__ import absolute_import, division, print_function
import socket
from copy import deepcopy
import glob
from simtbx.diffBragg import utils, hopper_utils
from dxtbx.model.experiment_list import ExperimentListFactory
import time
from xfel.merging.application.input.file_loader import create_experiment_identifier
from xfel.poly.recompute_mosaic_params import extract_mosaic_parameters_using_lambda_spread
from dials.algorithms.shoebox import MaskCode
import sys
try:
    import pandas
except ImportError:
    print("Please install pandas, libtbx.python -m pip install pandas")
    exit()

try:
    from line_profiler import LineProfiler
except ImportError:
    LineProfiler = None

from simtbx.diffBragg.device import DeviceWrapper

ROTX_ID = 0
ROTY_ID = 1
ROTZ_ID = 2
NCELLS_ID = 9
UCELL_ID_OFFSET = 3
DETZ_ID = 10

# LIBTBX_SET_DISPATCHER_NAME simtbx.diffBragg.hopper
# LIBTBX_SET_DISPATCHER_NAME hopper

import numpy as np
np.seterr(invalid='ignore')
from scipy.stats import pearsonr
import os
from libtbx.mpi4py import MPI

COMM = MPI.COMM_WORLD
# TODO, figure out why next 3 lines are sometimes necessary?!
if not hasattr(COMM, "rank"):
    COMM.rank = 0
    COMM.size = 1
from libtbx.phil import parse

from simtbx.diffBragg import utils
from simtbx.diffBragg.phil import philz
from simtbx.modeling import predictions
from dials.model.data import Shoebox
from dials.array_family import flex
from simtbx.command_line import integrate
from dxtbx.model import ExperimentList

import logging
from simtbx.diffBragg.phil import hopper_phil

philz = hopper_phil + philz
phil_scope = parse(philz)



class Script:
    def __init__(self):
        from dials.util.options import ArgumentParser

        self.params = None
        if COMM.rank == 0:
            self.parser = ArgumentParser(
                usage="",  # stage 1 (per-shot) diffBragg refinement",
                sort_options=True,
                phil=phil_scope,
                read_experiments=True,
                read_reflections=True,
                check_format=False,
                epilog="PyCuties")
            self.params, _ = self.parser.parse_args(show_diff_phil=True)
            assert self.params.outdir is not None
            utils.safe_makedirs(self.params.outdir)
            ts = time.strftime("%Y%m%d-%H%M%S")
            diff_phil_outname = os.path.join(self.params.outdir, "diff_phil_run_at_%s.txt" % ts)
            with open(diff_phil_outname, "w") as o:
                o.write("command line:\n%s\n" % (" ".join(sys.argv)))
                o.write("workding directory: \n%s\n" %os.getcwd())
                o.write("diff phil:\n")
                o.write(self.parser.diff_phil.as_str())
            just_diff_phil_outname = os.path.join(self.params.outdir, "diff.phil")
            with open(just_diff_phil_outname, "w") as o:
                o.write(self.parser.diff_phil.as_str())
        self.params = COMM.bcast(self.params)

        self.dev = COMM.rank % self.params.refiner.num_devices
        logging.info("Rank %d will use device %d on host %s" % (COMM.rank, self.dev, socket.gethostname()))

        if self.params.logging.logname is None:
            self.params.logging.logname = "main_stage1.log"
        if self.params.profile_name is None:
            self.params.profile_name = "prof_stage1.log"
        from simtbx.diffBragg import mpi_logger
        mpi_logger.setup_logging_from_params(self.params)

    def run(self):
        MAIN_LOGGER = logging.getLogger("diffBragg.main")
        assert os.path.exists(self.params.exp_ref_spec_file)
        input_lines = None
        best_models = None
        pd_dir = os.path.join(self.params.outdir, "pandas")
        expt_ref_dir = os.path.join(self.params.outdir, "expers_refls")
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

            utils.safe_makedirs(pd_dir)
            utils.safe_makedirs(expt_ref_dir)

        COMM.barrier()
        input_lines = COMM.bcast(input_lines)
        best_models = COMM.bcast(best_models)

        if self.params.ignore_existing:
            exp_names_already =None
            refl_names_already = None
            if COMM.rank==0:
                exp_names_already = {os.path.basename(f) for f in glob.glob("%s/expers/rank*/*.expt" % self.params.outdir)}
                refl_names_already = {os.path.basename(f) for f in glob.glob("%s/refls/rank*/*.refl" % self.params.outdir)}
            exp_names_already = COMM.bcast(exp_names_already)
            refl_names_already = COMM.bcast(refl_names_already)

        exp_gatheredRef_spec = []  # optional list of expt, refls, spectra
        trefs = []
        this_rank_dfs = []  # dataframes storing the modeling results for each shot
        this_rank_Elints = ExperimentList()
        this_rank_Rints = None
        this_rank_Ridxs = None
        chunk_id = 0
        shots_per_chunk=self.params.shots_per_chunk
        CHECKER = None
        try:
            from simtbx.tests import roi_check
            CHECKER = roi_check.roiCheck(self.params.roi.filter_scores.state_file)
        except:
            CHECKER = None
        for i_shot, line in enumerate(input_lines):
            time_to_refine = time.time()
            if i_shot == self.params.max_process:
                break
            if i_shot % COMM.size != COMM.rank:
                continue

            logging.info("COMM.rank %d on shot  %d / %d" % (COMM.rank, i_shot + 1, len(input_lines)))
            exp, ref, exp_idx, spec = hopper_utils.split_line(line)

            if self.params.ignore_existing:
                basename = os.path.splitext(os.path.basename(exp))[0]
                exists = False
                for ii in [i_shot, 0]:
                    opt_exp = "%s_%s_%d_%d.expt" % (self.params.tag, basename, exp_idx, ii)
                    opt_refl = opt_exp.replace(".expt", ".refl")
                    if opt_exp in exp_names_already and opt_refl in refl_names_already:
                        exists = True
                        break
                if exists:
                    print("Found existing!! %d" % i_shot)
                    continue

            best = None
            if best_models is not None:
                if "exp_idx" not in list(best_models):
                    best_models["exp_idx"]= 0
                best = best_models.query("exp_name=='%s'" % os.path.abspath(exp)).query("exp_idx==%d" % exp_idx)

                if len(best) != 1:
                    raise ValueError("Should be 1 entry for exp %s in best pickle %s" % (exp, self.params.best_pickle))
            self.params.simulator.spectrum.filename = spec
            Modeler = hopper_utils.DataModeler(self.params)
            Modeler.exper_name = exp
            Modeler.exper_idx = exp_idx
            Modeler.refl_name = ref
            Modeler.rank = COMM.rank
            Modeler.i_shot = i_shot
            # Optional prediction step?
            if self.params.load_data_from_refls:
                gathered = Modeler.GatherFromReflectionTable(exp, ref, sg_symbol=self.params.space_group)
            else:
                gathered = Modeler.GatherFromExperiment(exp, ref,
                                                        remove_duplicate_hkl=self.params.remove_duplicate_hkl,
                                                        sg_symbol=self.params.space_group,
                                                        exp_idx=exp_idx)
            if not gathered:
                logging.warning("No refls in %s; CONTINUE; COMM.rank=%d" % (ref, COMM.rank))
                continue
            MAIN_LOGGER.info("Modeling %s (%d refls)" % (exp, len(Modeler.refls)))
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

            if self.params.refiner.reference_geom is not None:
                detector = ExperimentListFactory.from_json_file(self.params.refiner.reference_geom, check_format=False)[0].detector
                Modeler.E.detector = detector

            # here we support inputting an experiment list with multiple crystals
            # the first crystal in the exp list is used to instantiate a diffBragg instance,
            # the remaining crystals are added to the sim_data instance for use during hopper_utils modeling
            # best pickle is not supported yet for multiple crystals
            # also, if number of crystals is >1 , then the params.number_of_xtals flag will be overridden
            exp_list = ExperimentListFactory.from_json_file(exp, False)
            xtals = exp_list.crystals()  # TODO: fix as this is broken now that we allow multi image experiments
            if self.params.consider_multicrystal_shots and len(xtals) > 1:
                assert best is None, "cannot pass best pickle if expt list has more than one crystal"
                assert self.params.number_of_xtals==1, "if expt list has more than one xtal, leave number_of_xtals as the default"
                self.params.number_of_xtals = len(xtals)
                MAIN_LOGGER.debug("Found %d xtals with unit cells:" %len(xtals))
                for xtal in xtals:
                    MAIN_LOGGER.debug("%.4f %.4f %.4f %.4f %.4f %.4f" % xtal.get_unit_cell().parameters())
            if self.params.record_device_timings and COMM.rank >0:
                self.params.record_device_timings = False  # only record for rank 0 otherwise there's too much output
            SIM = hopper_utils.get_simulator_for_data_modelers(Modeler)
            Modeler.set_parameters_for_experiment(best)
            MAIN_LOGGER.debug("Set parameters for experiment")
            Modeler.Umatrices = [Modeler.E.crystal.get_U()]

            # TODO: move this to SimulatorFromExperiment
            # TODO: fix multi crystal shot mode
            if best is not None and "other_spotscales" in list(best) and "other_Umats" in list(best):
                Modeler.Umatrices[0] = Modeler.E.get_U()
                assert len(xtals) == len(best.other_spotscales.values[0])+1
                for i_xtal in range(1, len(xtals),1):
                    scale_xt = best.other_spotscales.values[0][i_xtal]
                    Umat_xt = best.other_Umats.values[0][i_xtal]
                    Modeler.Umatrices[i_xtal] = Umat_xt
                    Modeler.P["G_xtal%d" %i_xtal] = scale_xt

            SIM.D.store_ave_wavelength_image = self.params.store_wavelength_images
            if self.params.refiner.verbose is not None and COMM.rank==0:
                SIM.D.verbose = self.params.refiner.verbose
            if self.params.profile:
                SIM.record_timings = True
            if self.params.use_float32:
                Modeler.all_data = Modeler.all_data.astype(np.float32)
                Modeler.all_background = Modeler.all_background.astype(np.float32)

            SIM.D.device_Id = self.dev

            nparam = len(Modeler.P)
            if SIM.refining_Fhkl:
                nparam += SIM.Num_ASU*SIM.num_Fhkl_channels
            x0 = [1] * nparam
            tref = time.time()
            MAIN_LOGGER.info("Beginning refinement of shot %d / %d" % (i_shot+1, len(input_lines)))
            try:
                x = Modeler.Minimize(x0, SIM, i_shot=i_shot)
                for i_rep in range(self.params.filter_after_refinement.max_attempts):
                    if not self.params.filter_after_refinement.enable:
                        continue
                    final_sigz = Modeler.target.all_sigZ[-1]
                    niter = len(Modeler.target.all_sigZ)
                    too_few_iter = niter < self.params.filter_after_refinement.min_prev_niter
                    too_high_sigz = final_sigz > self.params.filter_after_refinement.max_prev_sigz
                    if too_few_iter or too_high_sigz:
                        Modeler.filter_pixels(self.params.filter_after_refinement.threshold)
                        x = Modeler.Minimize(x0, SIM, i_shot=i_shot)

                if self.params.perRoi_finish:
                    old_params = deepcopy(self.params)
                    old_P = deepcopy(Modeler.P)
                    # fix all of the refinement variables except for perRoiScale
                    for fix_name in dir(self.params.fix):
                        if fix_name.startswith("_"):
                            continue
                        setattr(self.params.fix, fix_name, True)
                    self.params.fix.perRoiScale = False
                    Modeler.params = self.params
                    Modeler.set_parameters_for_experiment(best)
                    new_x = np.array([1.]*len(Modeler.P))
                    for name in old_P:
                        new_p = Modeler.P[name]
                        old_p = old_P[name]
                        new_x[new_p.xpos] = x[old_p.xpos]
                        assert not new_p.refine

                    x = Modeler.Minimize(new_x, SIM, i_shot=i_shot)

                    # reset the params
                    self.params = old_params

                    # Filter poor fits
                    if self.params.roi.filter_scores.enable and CHECKER is not None:
                        Modeler.best_model, _ = hopper_utils.model(x, Modeler, SIM, compute_grad=False)
                        Modeler.best_model_includes_background = False
                        num_good = Modeler.filter_bad_scores(CHECKER)
                        if num_good == 0:
                            Modeler.clean_up(SIM)
                            continue

                    # Then repeat minimization, fixing roi fits
                    Modeler.params = self.params
                    old_P = deepcopy(Modeler.P)
                    Modeler.set_parameters_for_experiment(best)
                    new_x = np.array([1.] * len(Modeler.P))
                    for name in Modeler.P:
                        new_p = Modeler.P[name]
                        old_p = old_P[name]
                        new_x[new_p.xpos] = x[old_p.xpos]
                    x = Modeler.Minimize(new_x, SIM, i_shot=i_shot)

            except StopIteration:
                x = Modeler.target.x0
            tref = time.time()-tref
            sigz = niter = None
            try:
                sigz, niter, _ = Modeler.get_best_hop()
            except Exception:
                pass

            trefs.append(tref)
            print_s = "Finished refinement of shot %d / %d in %.4f sec. (rank mean t/im=%.4f sec.)" \
                        % (i_shot+1, len(input_lines), tref, np.mean(trefs))
            if sigz is not None and niter is not None:
                print_s += " Ran %d iterations. Final sigmaZ = %.1f," % (niter, sigz)
            if COMM.rank==0:
                MAIN_LOGGER.info(print_s)
            else:
                MAIN_LOGGER.debug(print_s)
            if self.params.profile:
                SIM.D.show_timings(COMM.rank)

            dbg = self.params.debug_mode
            if dbg and COMM.rank > 0 and self.params.debug_mode_rank0_only:
                dbg = False
            shot_df = Modeler.save_up(x, SIM, rank=COMM.rank, i_shot=i_shot,
                            save_fhkl_data=dbg, save_refl=dbg, save_modeler_file=dbg,
                            save_sim_info=dbg, save_pandas=dbg, save_traces=dbg, save_expt=dbg)

            #if self.params.predictions.integrate_phil is not None:
            do_integrate = self.params.predictions.integrate_phil is not None
            Rstrong = None
            #Modeler.params.predictions.threshold = 10
            #Modeler.params.predictions.use_peak_detection = True
            #Modeler.params.predictions.use_diffBragg_mtz = True
            #pid = Modeler.pids[0]
            #x1, x2, y1, y2 = Modeler.rois[0]
            #fast = int((x2 + x1) / 2)
            #slow = int((y2 + y1) / 2)
            #Modeler.params.predictions.printout_pix = (pid, fast, slow)
            if do_integrate:
                pred_out = predictions.get_predict(Modeler.E,
                                Rstrong, Modeler.params, dev=self.dev, df=shot_df, filter_dupes=True, return_pix=True)
                pred = imgs = None
                if pred_out is None:
                    Modeler.clean_up(SIM)
                    continue
                pred, imgs = pred_out

            pfs = Modeler.pan_fast_slow.as_numpy_array()
            p, f, s = pfs[0::3], pfs[1::3], pfs[2::3]
            if do_integrate:
                predict_model = imgs[p, s, f]
            predict_subimgs = []
            hopper_subimgs = []

            ccs = []
            shoeboxes = []
            scores = []
            Ridx = flex.reflection_table()
            for roi, slc in Modeler.roi_id_slices.items():
                hopper_pix = Modeler.best_model[slc[0]]

                x1, x2, y1, y2 = Modeler.rois[int(roi)]
                ydim = y2 - y1
                xdim = x2 - x1

                hopper_subimg = hopper_pix.reshape((ydim, xdim))

                if do_integrate:
                    predict_pix = predict_model[slc[0]]
                    if len(np.unique(predict_pix)) == 1 or len(np.unique(hopper_pix)) == 1:
                        continue
                    cc = pearsonr(predict_pix, hopper_pix)[0]

                    predict_subimg = predict_pix.reshape((ydim, xdim))

                hopper_trust = Modeler.all_trusted[slc[0]].reshape((ydim,xdim))
                if not np.any(hopper_trust):
                    continue

                hopper_subimgs.append(hopper_subimg)
                hopper_bg = Modeler.all_background[slc[0]].reshape((ydim, xdim))
                if do_integrate:
                    ccs.append(cc)
                    predict_subimgs.append(predict_subimg)

                    sb = Shoebox((x1, x2, y1, y2, 0, 1))
                    sb.allocate()
                    sb.data = flex.float(np.ascontiguousarray(hopper_subimg[None]))
                    sb.background = flex.float(np.ascontiguousarray(hopper_bg[None]))

                    dials_mask = np.zeros((ydim, xdim)).astype(np.int32)
                    dials_mask[hopper_trust] = dials_mask[hopper_trust] + MaskCode.Valid
                    braggMask = (hopper_subimg > 10) & hopper_trust
                    dials_mask[braggMask] = dials_mask[braggMask] + MaskCode.Strong
                    bgMask = (hopper_subimg < 1) & hopper_trust
                    dials_mask[bgMask] = dials_mask[bgMask] + MaskCode.Background
                    sb.mask = flex.int(np.ascontiguousarray(dials_mask[None]))
                    shoeboxes.append(sb)

                refl = Modeler.refls[roi:roi+1]
                Ridx.extend(refl)
                if CHECKER is not None:
                    data_subimg = Modeler.all_data[slc[0]].reshape((ydim, xdim))
                    score = CHECKER.score(data_subimg, hopper_subimg + hopper_bg)
                    scores.append(score)
            if len(Ridx) == 0:
                Modeler.clean_up(SIM)
                continue

            if scores:
                Ridx["model_score"] = flex.double(scores)
            if shoeboxes and 'shoebox' != Ridx:
                Ridx['shoebox'] = flex.shoebox(shoeboxes)
            Ridx.set_flags(flex.bool(len(Ridx), True), Ridx.flags.indexed)
            from IPython import embed;embed()


            Elint = ExperimentList()
            Eref = deepcopy(Modeler.E)
            Eref.crystal.set_A(shot_df.Amats.values[0])
            Eref.detector = SIM.detector
            Elint.append(Eref)
            if do_integrate:
                Elint, Rint = integrate.integrate(self.params.predictions.integrate_phil, Elint, Ridx, pred)
                spot_ev, delpsi_rad, Deff, eta_est = extract_mosaic_parameters_using_lambda_spread(Elint[0], Rint, verbose=False)
                Rint["spot_ev"] = spot_ev
                Rint["delpsical.rad"] = delpsi_rad
                # here is where to do the diffBragg fit of intensity

            if Modeler.E.identifier is not None:
                ident = Modeler.E.identifier
            else:
                ident = create_experiment_identifier(Modeler.E, Modeler.exper_name, Modeler.exper_idx)

            Rint_id = Ridx_id = len(this_rank_Elints)
            if do_integrate:
                eid_int = Rint.experiment_identifiers()
                for k in eid_int.keys():
                    del eid_int[k]
                eid_int[Rint_id] = ident
                Rint['id'] = flex.int(len(Rint), Rint_id)

            eid_idx = Ridx.experiment_identifiers()
            for k in eid_idx.keys():
                del eid_idx[k]
            eid_idx[Ridx_id] = ident
            Ridx['id'] = flex.int(len(Ridx), Ridx_id)

            Elint[0].identifier = ident

            #ccs = [cc for cc in ccs if not np.isnan(cc)]
            #if not np.allclose(ccs, 1):
            #    from IPython import embed;embed()
            # verify the prediction is identical to the best model
            shot_df['identifier'] = ident
            shot_df["hopper_time"] = time.time()-time_to_refine
            shot_df["hopper_line"] = line  # exp_ref_spec file line for re-running
            if scores:
                shot_df["scores"] = np.mean(scores)

            # here is where to re-launch if doing predictions

            if Modeler.params.refiner.debug_pixel_panelfastslow is not None:
                # TODO separate diffBragg logger
                utils.show_diffBragg_state(SIM.D, Modeler.params.refiner.debug_pixel_panelfastslow)

            # TODO verify this works:
            if SIM.D.record_timings:
                SIM.D.show_timings(COMM.rank)
            Modeler.clean_up(SIM)
            del SIM.D  # TODO: is this necessary ?

            this_rank_dfs.append(shot_df)
            this_rank_Elints.extend(Elint)
            if do_integrate:
                if 'shoebox' in Rint:
                    del Rint['shoebox']
                if this_rank_Rints is None:
                    this_rank_Rints = Rint
                else:
                    this_rank_Rints.extend(Rint)
            #if "shoebox" in Ridx:
            #    del Ridx['shoebox']
            #TODO fill in new xyzcal.px column ?
            if this_rank_Ridxs is None:
                this_rank_Ridxs = Ridx
            else:
                this_rank_Ridxs.extend(Ridx)



            if len(this_rank_dfs) == shots_per_chunk:
                save_composite_files(this_rank_dfs, this_rank_Elints, this_rank_Ridxs, this_rank_Rints,
                                     pd_dir, expt_ref_dir, chunk_id)
                chunk_id += 1
                this_rank_dfs = []  # dataframes storing the modeling results for each shot
                this_rank_Elints = ExperimentList()
                this_rank_Rints = None
                this_rank_Ridxs = None

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

        save_composite_files(this_rank_dfs, this_rank_Elints, this_rank_Ridxs, this_rank_Rints,
                             pd_dir, expt_ref_dir, chunk_id)


def save_composite_files(dfs, expts, refls, refls_int, pd_dir, exp_ref_dir, chunk=0):
    df_name = os.path.join(pd_dir, "hopper_results_rank%d_chunk%d.pkl" % (COMM.rank, chunk))
    expt_name = os.path.join(exp_ref_dir, "hopper_rank%d_chunk%d.expt" % (COMM.rank,chunk))
    int_name = os.path.join(exp_ref_dir, "hopper_rank%d_chunk%d_integrated.refl" % (COMM.rank,chunk))
    idx_name = os.path.join(exp_ref_dir, "hopper_rank%d_chunk%d.refl" % (COMM.rank,chunk))
    if dfs:
        this_rank_dfs = pandas.concat(dfs).reset_index(drop=True)
        this_rank_dfs.to_pickle(df_name)
    if expts:
        expts.as_file(expt_name)
    if refls:
        refls.as_file(idx_name)
    if refls_int:
        refls_int.as_file(int_name)

        #MAIN_LOGGER.info("MPI-Gathering data frames across ranks")
        #all_rank_dfs = COMM.gather(this_rank_dfs)
        #if COMM.rank==0:
        #    all_rank_dfs = pandas.concat(all_rank_dfs)
        #    all_rank_dfs.reset_index(inplace=True, drop=True)
        #    all_df_name = os.path.join(self.params.outdir, "hopper_results.pkl")
        #    all_rank_dfs.to_pickle(all_df_name)


if __name__ == '__main__':
    from dials.util import show_mail_on_error

    with show_mail_on_error():
        print("Starting script")
        script = Script()
        RUN = script.run
        lp = None
        if LineProfiler is not None and script.params.profile:
            lp = LineProfiler()
            lp.add_function(hopper_utils.model)
            lp.add_function(hopper_utils.target_func)
            lp.add_function(RUN)
            lp.add_function(utils.get_roi_background_and_selection_flags)
            lp.add_function(hopper_utils.DataModeler.GatherFromExperiment)
            lp.add_function(utils.simulator_for_refinement)
            lp.add_function(utils.simulator_from_expt_and_params)
            lp.add_function(hopper_utils.get_simulator_for_data_modelers)
            RUN = lp(script.run)
        elif script.params.profile:
            print("Install line_profiler in order to use logging: libtbx.python -m pip install line_profiler")

        with DeviceWrapper(script.dev) as _:
            #with np.errstate(all='raise'):
            try:
                RUN()
            except Exception as err:
                err_file = os.path.join(script.params.outdir, "rank%d_hopper_fail.err" % COMM.rank)
                with open(err_file, "w") as o:
                    from traceback import format_tb
                    _, _, tb = sys.exc_info()
                    tb_s = "".join(format_tb(tb))
                    #tb_s = tb_s.replace("\n", "\nRANK%04d" % COMM.rank)
                    err_s = str(err) + "\n" + tb_s
                    o.write(err_s)
                raise err
        COMM.barrier()

        if lp is not None:
            stats = lp.get_stats()
            hopper_utils.print_profile(stats, ["model", "target_func", "run", "get_roi_background_and_selection_flags", "GatherFromExperiment",
                                               "simulator_for_refinement", "simulator_from_expt_and_params", "get_simulator_for_data_modelers"])
