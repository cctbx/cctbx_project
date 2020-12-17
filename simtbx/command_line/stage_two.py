from __future__ import absolute_import, division, print_function

# LIBTBX_SET_DISPATCHER_NAME simtbx.diffBragg.stage_two

from libtbx.mpi4py import MPI
from dxtbx.model import ExperimentList
import os
import time

COMM = MPI.COMM_WORLD

if COMM.rank > 0:
    import warnings
    warnings.filterwarnings("ignore")


from simtbx.diffBragg.phil import philz
from simtbx.diffBragg import ensemble_refine_launcher
from dxtbx.model.experiment_list import ExperimentListFactory
from dials.array_family import flex
from libtbx.phil import parse
from dials.util import show_mail_on_error

help_message = "stage 2 (global) diffBragg refinement"

script_phil = """
show_timing = True
  .type = bool
  .help = print a refinement duration for each iteration experiment
d_min = 2
 .type = float
 .help = high res lim for binner
d_max = 999
 .type = float
 .help = low res lim for binner
n_bins = 10
  .type = int
  .help = number of binner bins
"""

philz = script_phil + philz
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

    def run(self):
        self.params = None
        if COMM.rank == 0:
            self.params, _ = self.parser.parse_args(show_diff_phil=True)
        self.params = COMM.bcast(self.params)

        refine_starttime = time.time()
        if not self.params.refiner.randomize_devices:
            self.params.simulator.device_id = COMM.rank % self.params.refiner.num_devices
        refiner = ensemble_refine_launcher.global_refiner_from_parameters(refls, pandas_list, self.params)
        if self.params.show_timing:
            print("Time to refine experiment: %f" % (time.time()- refine_starttime))

        #TODO save MTZ


if __name__ == '__main__':
    with show_mail_on_error():
        script = Script()
        script.run()