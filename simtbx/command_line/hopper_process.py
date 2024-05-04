from __future__ import absolute_import, division, print_function
import time
import sys

# LIBTBX_SET_DISPATCHER_NAME diffBragg.stills_process
# LIBTBX_SET_DISPATCHER_NAME diffBragg.hopper_process

from dials.command_line.stills_process import Processor
from xfel.small_cell.small_cell import small_cell_index_detail
from simtbx.diffBragg import hopper_utils
from dials.array_family import flex
from simtbx.diffBragg.hopper_io import save_to_pandas
from dxtbx.model import ExperimentList
import numpy as np
import socket
from libtbx.mpi4py import MPI
COMM = MPI.COMM_WORLD

import logging
from libtbx.utils import Sorry
from simtbx.modeling import predictions

logger = logging.getLogger("dials.command_line.stills_process")


#include scope dials.command_line.stills_process.phil_scope
phil_str = """
include scope xfel.small_cell.command_line.small_cell_process.phil_scope
diffBragg {
  include scope simtbx.command_line.hopper.phil_scope
}
skip_hopper=False
  .type = bool
  .help = if True, then skip the hopper refinement, i.e. just run stills
  .help = process without refinement like usual
silence_dials_loggers = False
  .type = bool
  .help = if True, reduce the output clutter from the dials.refine and dials.index methods
save_pandas = True
  .type = bool
  .help = save the model parameters in a pandas frame
combine_pandas = True
  .type = bool
  .help = if True, combine all pandas frames after hopper_process completes
partial_correct = False
  .type = bool
  .help = if True, compute partialities from the prediction models  and use them to normalize measurements
  .help = for examining the modeling results, can be loaded using modeler = np.load(fname, allow_pickle=True)[()]
refspec = None
  .type = str
  .help = path to a reference .lam file to use as the spectra for each shot
reidx_obs = False
  .type = bool
  .help = optionally reindex the strong spot observations after running sills indexer refinement
db_loglevel = 0 1 *2
  .type = choice
  .help = Log level for diffBragg main logger
  .help = 0=critical (less verbose), 1=info (verbose), 2=debug (most verbose)
refine_predictions = False
  .type = bool
  .help = optionally refine the list of predicted reflections before integrating
use_small_cell_indexing = False
  .type = bool
  .help = use this to index sparse patterns from small unit cell xtals with at least 3 reflections
"""
import os
from libtbx.phil import parse
phil_scope = parse(phil_str, process_includes=True)


class Hopper_Processor(Processor):

    def __init__(self, *args, **kwargs):
        super(Hopper_Processor, self).__init__(*args, **kwargs)
        command_file = os.path.join(self.params.output.output_dir, "hopper_process_cmdline.txt")
        if COMM.rank==0:
            with open(command_file, "w") as o:
                o.write(" ".join(sys.argv))
                o.close()
            hop_proc_params_file = os.path.join(self.params.output.output_dir, "hopperProcessParams.npy")
            np.save(hop_proc_params_file, self.params)
        if self.params.output.composite_output:
            print("WARNING!! COMPOSITE OUTPUT NOT YET SUPPORTED, disabling for now")
        self.params.output.composite_output = False

        self.stage1_df = None  # the pandas dataframe containing model parameters after running stage 1 (see method refine)
        self.stage1_modeler = None  # the data modeler used during stage 1 refinement
        self.modeler_dir = None  # path for writing the data modelers
        self.known_crystal_models = None  # default flag, should be defined in stills_process init
        self.SIM = None # place holder for simtbx/nanoBragg/sim_data.SimData instance

        if self.params.silence_dials_loggers:
            #dials.algorithms.indexing.indexer: model
            logging.getLogger("dials.algorithms.indexing.nave_parameters").setLevel(logging.ERROR)
            logging.getLogger("dials.algorithms.indexing.stills_indexer").setLevel(logging.ERROR)
            logging.getLogger("dials.algorithms.refinement.refiner").setLevel(logging.ERROR)
            logging.getLogger("dials.algorithms.refinement.reflection_manager").setLevel(logging.ERROR)
            logging.getLogger("dials.algorithms.refinement.reflection_manager").setLevel(logging.ERROR)
        # configure the diffBragg logger
        dblog = logging.getLogger("diffBragg.main")
        H = logging.StreamHandler()
        if self.params.db_loglevel=='2':
            dblog.setLevel(logging.DEBUG)
            H.setLevel(logging.DEBUG)
        elif self.params.db_loglevel=='1':
            dblog.setLevel(logging.INFO)
            H.setLevel(logging.INFO)
        else:
            dblog.setLevel(logging.CRITICAL)
            H.setLevel(logging.CRITICAL)
        dblog.addHandler(H)

        self._create_modeler_dir()

    def index(self, experiments, reflections):
        """optionally do small cell indexing , else default to stills process"""
        if self.params.use_small_cell_indexing:
            max_clique_len, experiments, indexed = small_cell_index_detail(experiments, reflections, self.params, write_output=False)
        else:
            experiments, indexed = super(Hopper_Processor, self).index(experiments, reflections)

        return experiments, indexed

    def _create_modeler_dir(self):
        """makes the directory where data modelers, logs, and spectra files will be written"""
        if self.params.diffBragg.debug_mode:
            self.modeler_dir = os.path.join(self.params.output.output_dir, "modelers")
            if COMM.rank == 0:
                if not os.path.exists(self.modeler_dir):
                    os.makedirs(self.modeler_dir)
            COMM.barrier()

    @property
    def device_id(self):
        dev = COMM.rank % self.params.diffBragg.refiner.num_devices
        print("Rank %d will use fixed device %d on host %s" % (COMM.rank, dev, socket.gethostname()), flush=True)
        return dev

    def find_spots(self, experiments):
        st = time.time()

        logger.info("*" * 80)
        logger.info("Finding Strong Spots")
        logger.info("*" * 80)

        # Find the strong spots
        observed = flex.reflection_table.from_observations(
            experiments, self.params, is_stills=True
        )

        # Reset z coordinates for dials.image_viewer; see Issues #226 for details
        xyzobs = observed["xyzobs.px.value"]
        for i in range(len(xyzobs)):
            xyzobs[i] = (xyzobs[i][0], xyzobs[i][1], 0)
        bbox = observed["bbox"]
        for i in range(len(bbox)):
            bbox[i] = (bbox[i][0], bbox[i][1], bbox[i][2], bbox[i][3], 0, 1)

        if self.params.output.composite_output:
            n = len(self.all_strong_reflections.experiment_identifiers())
            for i, experiment in enumerate(experiments):
                refls = observed.select(observed["id"] == i)
                refls["id"] = flex.int(len(refls), n)
                del refls.experiment_identifiers()[i]
                refls.experiment_identifiers()[n] = experiment.identifier
                self.all_strong_reflections.extend(refls)
                n += 1
        else:
            # Save the reflections to file
            logger.info("\n" + "-" * 80)
            if self.params.output.strong_filename:
                self.save_reflections(observed, self.params.output.strong_filename)

        logger.info("")
        logger.info("Time Taken = %f seconds", time.time() - st)
        self.observed = observed  # note this is the only change needed to dials.stills_process.find_spots
        return observed

    def refine(self, exps, ref, refining_predictions=False, best=None):
        exps_out = exps
        if not self.params.skip_hopper:
            if self.params.dispatch.refine:
                exps, ref = super(Hopper_Processor, self).refine(exps, ref)
                print("WARNING: hopper_process will always run its own refinement, ignoring dials.refine phil scope")
            self.params.dispatch.refine = False
            assert len(exps) == 1
            if self.params.reidx_obs:
                exps, ref = self._reindex_obs(exps, self.observed)

            exp, ref, self.stage1_modeler, self.SIM, x = hopper_utils.refine(exps[0], ref,
                                               self.params.diffBragg,
                                               spec=self.params.refspec,
                                               gpu_device=self.device_id, return_modeler=True, best=best)
            orig_exp_name = os.path.abspath(self.params.output.refined_experiments_filename)
            refls_name = os.path.abspath(self.params.output.indexed_filename)
            self.params.diffBragg.outdir = self.params.output.output_dir
            # TODO: what about composite mode ?

            #def save_to_pandas(x, Mod, SIM, orig_exp_name, params, expt, rank_exp_idx, stg1_refls, stg1_img_path=None,
            #                   rank=0, write_expt=True, write_pandas=True, exp_idx=0):
            self.stage1_df = save_to_pandas(x, self.stage1_modeler, self.SIM, orig_exp_name, self.params.diffBragg,
                                            self.stage1_modeler.E, rank_exp_idx=0, stg1_refls=refls_name, stg1_img_path=None, rank=COMM.rank,
                                            write_expt=False, write_pandas=True, exp_idx=0)
            exps_out = ExperimentList()
            exps_out.append(exp)

            basename = os.path.splitext(os.path.basename(refls_name))[0]
            self._save_modeler_info(basename)

        out = super(Hopper_Processor, self).refine(exps_out, ref)
        return out

    def _reindex_obs(self, exps, ref):
        """use known_crystal_models indexing method to add more indexed spots which can
        then be refined using diffBragg
        This is useful when dials.stills_indexer happens to prefer a lower resolution cutoff
        to obtain a descent xtal model (controlled by refinement_protocol.d_min_start)
        but one still wants to try refining the higher resolution spots obtained with spot finder
        """
        # NOTE: this method assumes known_crystal_models is not being used in any other way ...
        self.known_crystal_models = exps.crystals()
        # cache these parameters in case they are set for stills_indexer
        tmp_tol = self.params.indexing.index_assignment.simple.hkl_tolerance
        tmp_prot = self.params.indexing.refinement_protocol.mode
        self.params.indexing.index_assignment.simple.hkl_tolerance = 0.5  # go for broke!
        self.params.indexing.refinement_protocol.mode = "repredict_only"  # no more refinement from dials
        exps, ref = self.index(exps, ref)

        self.params.indexing.index_assignment.simple.hkl_tolerance = tmp_tol
        self.params.indexing.refinement_protocol.mode = tmp_prot
        self.known_crystal_models = None
        return exps, ref

    def _save_modeler_info(self, basename):
        if self.params.diffBragg.debug_mode:
            modeler_fname = "%s_%s" % (basename, "modeler.npy")
            modeler_fname = os.path.abspath(os.path.join(self.modeler_dir, modeler_fname))
            np.save(modeler_fname, self.stage1_modeler)  # pickle the modeler, set __setstate__

            spectra_fname = "%s_%s" % (basename, "spectra.npy")
            spectra_fname = os.path.abspath(os.path.join(self.modeler_dir, spectra_fname))
            SIM_state_fname = "%s_%s" % (basename, "SIM_state.npy")
            SIM_state_fname = os.path.abspath(os.path.join(self.modeler_dir, SIM_state_fname))
            hopper_utils.write_SIM_logs(self.SIM, log=SIM_state_fname, lam=spectra_fname)

    def integrate(self, experiments, indexed):
        if self.params.skip_hopper:
            return super(Hopper_Processor, self).integrate(experiments, indexed)
        st = time.time()

        logger.info("*" * 80)
        logger.info("Integrating Reflections")
        logger.info("*" * 80)

        indexed, _ = self.process_reference(indexed)

        if self.params.integration.integration_only_overrides.trusted_range:
            for detector in experiments.detectors():
                for panel in detector:
                    panel.set_trusted_range(
                        self.params.integration.integration_only_overrides.trusted_range
                    )

        if self.params.dispatch.coset:
            from dials.algorithms.integration.sublattice_helper import integrate_coset

            integrate_coset(self, experiments, indexed)

        # Get the integrator from the input parameters
        logger.info("Configuring integrator from input parameters")
        from dials.algorithms.integration.integrator import create_integrator
        from dials.algorithms.profile_model.factory import ProfileModelFactory

        # Compute the profile model
        # Predict the reflections
        # Match the predictions with the reference
        # Create the integrator
        experiments = ProfileModelFactory.create(self.params, experiments, indexed)
        new_experiments = ExperimentList()
        new_reflections = flex.reflection_table()
        for expt_id, expt in enumerate(experiments):
            if (
                self.params.profile.gaussian_rs.parameters.sigma_b_cutoff is None
                or expt.profile.sigma_b()
                < self.params.profile.gaussian_rs.parameters.sigma_b_cutoff
            ):
                refls = indexed.select(indexed["id"] == expt_id)
                refls["id"] = flex.int(len(refls), len(new_experiments))
                # refls.reset_ids()
                del refls.experiment_identifiers()[expt_id]
                refls.experiment_identifiers()[len(new_experiments)] = expt.identifier
                new_reflections.extend(refls)
                new_experiments.append(expt)
            else:
                logger.info(
                    "Rejected expt %d with sigma_b %f"
                    % (expt_id, expt.profile.sigma_b())
                )
        experiments = new_experiments
        indexed = new_reflections
        if len(experiments) == 0:
            raise RuntimeError("No experiments after filtering by sigma_b")
        logger.info("")
        logger.info("=" * 80)
        logger.info("")
        logger.info("Predicting reflections")
        logger.info("")
        # NOTE: this is the only changed needed to dials.stills_process
        # TODO: multi xtal
        # TODO: add in normal dials predictions as an option
        predicted, model = predictions.get_predicted_from_pandas(
            self.stage1_df, self.params.diffBragg, self.observed,
            experiments[0].identifier, self.device_id,
            spectrum_override=self.SIM.beam.spectrum)
        if self.params.refine_predictions:
            experiments, rnd2_refls = self.refine(experiments, predicted, refining_predictions=True, best=self.stage1_df)
            # TODO: match rnd2_refls with indexed.refl and re-save indexed.refl
            predicted, model = predictions.get_predicted_from_pandas(
                self.stage1_df, self.params.diffBragg, self.observed,
                experiments[0].identifier, self.device_id,
                spectrum_override=self.SIM.beam.spectrum)

        predicted.match_with_reference(indexed)
        integrator = create_integrator(self.params, experiments, predicted)

        # Integrate the reflections
        integrated = integrator.integrate()

        if self.params.partial_correct:
            integrated = predictions.normalize_by_partiality(
                integrated, model, default_F=self.params.diffBragg.predictions.default_Famplitude,
                gain=self.params.diffBragg.refiner.adu_per_photon)

        # correct integrated intensities for absorption correction, if necessary
        for abs_params in self.params.integration.absorption_correction:
            if abs_params.apply:
                if abs_params.algorithm == "fuller_kapton":
                    from dials.algorithms.integration.kapton_correction import (
                        multi_kapton_correction,
                    )
                elif abs_params.algorithm == "kapton_2019":
                    from dials.algorithms.integration.kapton_2019_correction import (
                        multi_kapton_correction,
                    )

                experiments, integrated = multi_kapton_correction(
                    experiments, integrated, abs_params.fuller_kapton, logger=logger
                )()

        if self.params.significance_filter.enable:
            from dials.algorithms.integration.stills_significance_filter import (
                SignificanceFilter,
            )

            sig_filter = SignificanceFilter(self.params)
            filtered_refls = sig_filter(experiments, integrated)
            accepted_expts = ExperimentList()
            accepted_refls = flex.reflection_table()
            logger.info(
                "Removed %d reflections out of %d when applying significance filter",
                len(integrated) - len(filtered_refls),
                len(integrated),
                )
            for expt_id, expt in enumerate(experiments):
                refls = filtered_refls.select(filtered_refls["id"] == expt_id)
                if len(refls) > 0:
                    accepted_expts.append(expt)
                    refls["id"] = flex.int(len(refls), len(accepted_expts) - 1)
                    accepted_refls.extend(refls)
                else:
                    logger.info(
                        "Removed experiment %d which has no reflections left after applying significance filter",
                        expt_id,
                    )

            if len(accepted_refls) == 0:
                raise Sorry("No reflections left after applying significance filter")
            experiments = accepted_expts
            integrated = accepted_refls

        # Delete the shoeboxes used for intermediate calculations, if requested
        if self.params.integration.debug.delete_shoeboxes and "shoebox" in integrated:
            del integrated["shoebox"]

        if self.params.output.composite_output:
            if (
                self.params.output.integrated_experiments_filename
                or self.params.output.integrated_filename
            ):
                assert (
                    self.params.output.integrated_experiments_filename is not None
                    and self.params.output.integrated_filename is not None
                )

                n = len(self.all_integrated_experiments)
                self.all_integrated_experiments.extend(experiments)
                for i, experiment in enumerate(experiments):
                    refls = integrated.select(integrated["id"] == i)
                    refls["id"] = flex.int(len(refls), n)
                    del refls.experiment_identifiers()[i]
                    refls.experiment_identifiers()[n] = experiment.identifier
                    self.all_integrated_reflections.extend(refls)
                    n += 1
        else:
            # Dump experiments to disk
            if self.params.output.integrated_experiments_filename:

                experiments.as_json(self.params.output.integrated_experiments_filename)

            if self.params.output.integrated_filename:
                # Save the reflections
                self.save_reflections(
                    integrated, self.params.output.integrated_filename
                )

        self.write_integration_pickles(integrated, experiments)
        from dials.algorithms.indexing.stills_indexer import (
            calc_2D_rmsd_and_displacements,
        )

        rmsd_indexed, _ = calc_2D_rmsd_and_displacements(indexed)
        log_str = "RMSD indexed (px): %f\n" % rmsd_indexed
        for i in range(6):
            bright_integrated = integrated.select(
                (
                    integrated["intensity.sum.value"]
                    / flex.sqrt(integrated["intensity.sum.variance"])
                )
                >= i
            )
            if len(bright_integrated) > 0:
                rmsd_integrated, _ = calc_2D_rmsd_and_displacements(bright_integrated)
            else:
                rmsd_integrated = 0
            log_str += (
                "N reflections integrated at I/sigI >= %d: % 4d, RMSD (px): %f\n"
                % (i, len(bright_integrated), rmsd_integrated)
            )

        for crystal_model in experiments.crystals():
            if hasattr(crystal_model, "get_domain_size_ang"):
                log_str += ". Final ML model: domain size angstroms: {:f}, half mosaicity degrees: {:f}".format(
                    crystal_model.get_domain_size_ang(),
                    crystal_model.get_half_mosaicity_deg(),
                )

        logger.info(log_str)

        logger.info("")
        logger.info("Time Taken = %f seconds", time.time() - st)
        return integrated


def run(args=None):
    from dials.command_line import stills_process
    stills_process.Processor = Hopper_Processor
    stills_process.phil_scope = phil_scope
    script = stills_process.Script()
    from simtbx.diffBragg.device import DeviceWrapper
    params, options, all_paths = script.parser.parse_args(
        args, show_diff_phil=True, return_unhandled=True, quick_parse=True
    )
    if params.output.composite_output:
        raise NotImplementedError("diffBragg.stills_process currently does not support composite_output mode")
    dev = COMM.rank % params.diffBragg.refiner.num_devices
    with DeviceWrapper(dev) as _:
        script.run(args)
    return script


if __name__ == "__main__":
    script_that_was_run = run()
    if COMM.rank==0:
        try:
            params = script_that_was_run.params
        except AttributeError as err:
            print(err)
            print("Looks like the program never launched, check input paths, image files, phil files, current working dir etc.. ")
            sys.exit()
        if params.combine_pandas:
            if not params.save_pandas:
                print("No pandas tables saved, so will not combine")
                exit()
            import pandas
            import glob
            fnames = glob.glob("%s/pandas/rank*/*pkl" % params.output.output_dir)
            logging.info("There are %d pandas output files to combine" % len(fnames))
            if fnames:
                df = pandas.concat([pandas.read_pickle(f) for f in fnames])
                combined_table = os.path.join(params.output.output_dir, "hopper_process_summary.pkl")
                df.to_pickle(combined_table)
                logging.info("Saved summary pandas table: %s" % combined_table)
