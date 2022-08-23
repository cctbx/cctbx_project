from __future__ import absolute_import, division, print_function
import socket
import glob
from simtbx.diffBragg import hopper_utils
from dxtbx.model.experiment_list import ExperimentListFactory
import time
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
            refl_names_already = None
            if COMM.rank==0:
                exp_names_already = {os.path.basename(f) for f in glob.glob("%s/expers/rank*/*.expt" % self.params.outdir)}
                refl_names_already = {os.path.basename(f) for f in glob.glob("%s/refls/rank*/*.refl" % self.params.outdir)}
            exp_names_already = COMM.bcast(exp_names_already)
            refl_names_already = COMM.bcast(refl_names_already)

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
                exists = False
                for ii in [i_exp, 0]:
                    opt_exp = "%s_%s_%d.expt" % (self.params.tag, basename, ii)
                    opt_refl = opt_exp.replace(".expt", ".refl")
                    if opt_exp in exp_names_already and opt_refl in refl_names_already:
                        exists = True
                        break
                if exists:
                    print("Found existing!! %d" % i_exp)
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
            Modeler.exper_name = exp
            Modeler.refl_name = ref
            Modeler.rank = COMM.rank
            Modeler.i_exp = i_exp
            if self.params.load_data_from_refls:
                gathered = Modeler.GatherFromReflectionTable(exp, ref, sg_symbol=self.params.space_group)
            else:
                gathered = Modeler.GatherFromExperiment(exp, ref,
                                                        remove_duplicate_hkl=self.params.remove_duplicate_hkl, sg_symbol=self.params.space_group)
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
            xtals = exp_list.crystals()
            if len(xtals) > 1:
                assert best is None, "cannot pass best pickle if expt list has more than one crystal"
                assert self.params.number_of_xtals==1, "if expt list has more than one xtal, leave number_of_xtals as the default"
                self.params.number_of_xtals = len(xtals)
                MAIN_LOGGER.debug("Found %d xtals with unit cells:" %len(xtals))
                for xtal in xtals:
                    MAIN_LOGGER.debug("%.4f %.4f %.4f %.4f %.4f %.4f" % xtal.get_unit_cell().parameters())
            Modeler.SimulatorFromExperiment(best)
            Modeler.SIM.Umatrices = [xtal.get_U() for xtal in xtals]
            # TODO, move this to SimulatorFromExperiment
            if best is not None and "other_spotscales" in list(best) and "other_Umats" in list(best):
                Modeler.SIM.Umatrices[0] = self.E.get_U()
                assert len(xtals) == len(best.other_spotscales.values[0])+1
                for i_xtal in range(1, len(xtals),1):
                    scale_xt = best.other_spotscales.values[0][i_xtal]
                    Umat_xt = best.other_Umats.values[0][i_xtal]
                    Modeler.SIM.Umatrices[i_xtal] = Umat_xt
                    Modeler.SIM.P["G_xtal%d" %i_xtal] = scale_xt

            Modeler.SIM.D.store_ave_wavelength_image = self.params.store_wavelength_images
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
            #print("RANK%d, DEV=%d" % (COMM.rank, dev))

            Modeler.SIM.D.device_Id = dev

            nparam = len(Modeler.SIM.P)
            x0 = [1] * nparam
            try:
                x = Modeler.Minimize(x0)
            except StopIteration:
                x = Modeler.target.x0
            if self.params.profile:
                Modeler.SIM.D.show_timings(COMM.rank)

            Modeler.save_up(x, rank=COMM.rank, i_exp=i_exp)
            if Modeler.params.refiner.debug_pixel_panelfastslow is not None:
                # TODO separate diffBragg logger
                utils.show_diffBragg_state(Modeler.SIM.D, Modeler.params.refiner.debug_pixel_panelfastslow)


            Modeler.clean_up()

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
