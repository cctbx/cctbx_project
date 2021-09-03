from __future__ import absolute_import, division, print_function

# LIBTBX_SET_DISPATCHER_NAME simtbx.diffBragg.stage_two

from libtbx.mpi4py import MPI
from simtbx.command_line.hopper import hopper_phil
import time
import logging
from simtbx.diffBragg import mpi_logger

COMM = MPI.COMM_WORLD

if COMM.rank > 0:
    import warnings
    warnings.filterwarnings("ignore")


from simtbx.diffBragg.phil import philz
from simtbx.diffBragg import ensemble_refine_launcher
from libtbx.phil import parse
from dials.util import show_mail_on_error

help_message = "stage 2 (global) diffBragg refinement"

script_phil = """
pandas_table = None
  .type = str 
  .help = path to an input pandas table (usually output by simtbx.diffBragg.predictions)
prep_time = 60
  .type = float
  .help = Time spent optimizing order of input dataframe to better divide shots across ranks
  .help = Unit is seconds, 1-2 minutes of prep might save a lot of time during refinement!
"""

philz = script_phil + philz + hopper_phil
phil_scope = parse(philz)


class Script:

    def __init__(self):
        from dials.util.options import OptionParser

        self.parser = None
        if COMM.rank == 0:
            self.parser = OptionParser(
                usage="",  # stage 1 (per-shot) diffBragg refinement",
                sort_options=True,
                phil=phil_scope,
                read_experiments=False,
                read_reflections=False,
                check_format=False,
                epilog=help_message)
        self.parser = COMM.bcast(self.parser)
        self.params, _ = self.parser.parse_args(show_diff_phil=COMM.rank==0)

    def run(self):
        #self.params = None
        #if COMM.rank == 0:
        #self.params = COMM.bcast(self.params)
        if self.params.pandas_table is None:
            raise ValueError("Pandas table input required")

        refine_starttime = time.time()
        if not self.params.refiner.randomize_devices:
            self.params.simulator.device_id = COMM.rank % self.params.refiner.num_devices
        refiner = ensemble_refine_launcher.global_refiner_from_parameters(self.params)
        print("Time to refine experiment: %f" % (time.time()- refine_starttime))

        #TODO save MTZ


if __name__ == '__main__':
    try:
        from line_profiler import LineProfiler
    except ImportError:
        LineProfiler = None
    from simtbx.diffBragg.refiners import local_refiner
    from simtbx.diffBragg import hopper_utils
    with show_mail_on_error():
        script = Script()
        RUN = script.run
        lp = None
        if LineProfiler is not None and script.params.profile:
            lp = LineProfiler()
            lp.add_function(ensemble_refine_launcher.RefineLauncher.launch_refiner)
            lp.add_function(local_refiner.LocalRefiner._compute_functional_and_gradients)
            lp.add_function(local_refiner.LocalRefiner._run_diffBragg_current)
            lp.add_function(local_refiner.LocalRefiner._update_Fcell)
            lp.add_function(local_refiner.LocalRefiner._extract_pixel_data)
            lp.add_function(local_refiner.LocalRefiner._Fcell_derivatives)
            lp.add_function(local_refiner.LocalRefiner._mpi_aggregation)
            lp.add_function(local_refiner.LocalRefiner._setup)
            lp.add_function(hopper_utils.DataModeler.GatherFromExperiment)
            RUN = lp(script.run)

        if script.params.outdir is None:
            od = script.params.refiner.io.output_dir
            script.params.outdir = od if od is not None else '.'

        if script.params.logging.logname is None:
            script.params.logging.logname = "main_stage2.log"
        if script.params.profile_name is None:
            script.params.profile_name = "prof_stage2.log"
        if script.params.logging.disable:
            logging.disable(level=logging.CRITICAL)  # disables CRITICAL and below
        else:
            mpi_logger.setup_logging_from_params(script.params)

        RUN()

        if lp is not None:
            stats = lp.get_stats()
            from simtbx.diffBragg import hopper_utils
            hopper_utils.print_profile(stats,
                    ["launch_refiner", "_compute_functional_and_gradients", "_run_diffBragg_current",
                     "_update_Fcell", "_extract_pixel_data", "_Fcell_derivatives", "_mpi_aggregation",
                     "GatherFromExperiment", "_setup"])
