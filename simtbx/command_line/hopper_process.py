from __future__ import absolute_import, division, print_function
import time

# LIBTBX_SET_DISPATCHER_NAME simtbx.diffBragg.hopper_process

from dials.command_line.stills_process import Processor
from simtbx.diffBragg.hopper_utils import refine
from dials.array_family import flex
from simtbx.command_line.hopper import save_to_pandas
from dxtbx.model import ExperimentList
import numpy as np
import socket
from libtbx.mpi4py import MPI
COMM = MPI.COMM_WORLD

import logging
from libtbx.utils import Sorry
from simtbx.modeling import predictions

logger = logging.getLogger("dials.command_line.stills_process")


phil_str = """
include scope dials.command_line.stills_process.phil_scope
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
"""
import os
from libtbx.phil import parse
phil_scope = parse(phil_str, process_includes=True)


class Hopper_Processor(Processor):

    def __init__(self, *args, **kwargs):
        super(Hopper_Processor, self).__init__(*args, **kwargs)

        self.stage1_df = None  # the pandas dataframe containing model parameters after running stage 1 (see method refine)

        if self.params.silence_dials_loggers:
            logging.getLogger("dials.algorithms.indexing.nave_parameters").setLevel(logging.ERROR)
            logging.getLogger("dials.algorithms.indexing.stills_indexer").setLevel(logging.ERROR)
            logging.getLogger("dials.algorithms.refinement.refiner").setLevel(logging.ERROR)
            logging.getLogger("dials.algorithms.refinement.reflection_manager").setLevel(logging.ERROR)
            logging.getLogger("dials.algorithms.refinement.reflection_manager").setLevel(logging.ERROR)
        logging.basicConfig(level=logging.DEBUG)

    @property
    def device_id(self):
        if self.params.diffBragg.refiner.randomize_devices:
            dev = np.random.choice(self.params.diffBragg.refiner.num_devices)
            print("Rank %d will use random device %d on host %s" % (COMM.rank, dev, socket.gethostname()), flush=True)
        else:
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

    def refine(self, exps, ref):
        exps_out = exps
        if not self.params.skip_hopper:
            if self.params.dispatch.refine:
                print("WARNING: hopper_process will always run its own refinement, ignoring dials.refine phil scope")
            self.params.dispatch.refine = False
            assert len(exps)==1
            # TODO MPI select GPU device

            exp, ref, data_modeler, x = refine(exps[0], ref, self.params.diffBragg, gpu_device=self.device_id, return_modeler=True)
            orig_exp_name = os.path.abspath(self.params.output.refined_experiments_filename)
            refls_name = os.path.abspath(self.params.output.indexed_filename)
            self.params.diffBragg.outdir = self.params.output.output_dir
            # TODO: what about composite mode ?
            self.stage1_df = save_to_pandas(x, data_modeler.SIM, orig_exp_name, self.params.diffBragg, data_modeler.E, 0, refls_name, None)
            exps_out = ExperimentList()
            exps_out.append(exp)
        return super(Hopper_Processor, self).refine(exps_out, ref)

    def integrate(self, experiments, indexed):
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
            from xfel.util.sublattice_helper import integrate_coset

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
            experiments[0].identifier, self.device_id)

        predicted.match_with_reference(indexed)
        integrator = create_integrator(self.params, experiments, predicted)

        # Integrate the reflections
        integrated = integrator.integrate()

        if self.params.partial_correct:
            integrated = predictions.normalize_by_partiality(
                integrated, model, default_F=self.params.diffBragg.predictions.default_Famplitude)

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
        log_str = f"RMSD indexed (px): {rmsd_indexed:f}\n"
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
    script.run(args)
    return script


if __name__ == "__main__":
    script_that_was_run = run()
    if COMM.rank==0:
        params = script_that_was_run.params
        if params.combine_pandas:
            if not params.save_pandas:
                print("No pandas tables saved, so will not combine")
                exit()
            import pandas
            import glob
            fnames = glob.glob("%s/pandas/rank*/*pkl" % params.output.output_dir)
            logging.info("There are %d pandas output files to combine")
            df = pandas.concat([pandas.read_pickle(f) for f in fnames])
            combined_table = os.path.join(params.output.output_dir, "hopper_process_summary.pkl")
            df.to_pickle(combined_table)
            logging.info("Saved summary pandas table: %s" % combined_table)
