from __future__ import absolute_import, division, print_function


# LIBTBX_SET_DISPATCHER_NAME simtbx.diffBragg.hopper_process

from dials.command_line.stills_process import Processor
from simtbx.diffBragg.hopper_utils import refine
from dxtbx.model import ExperimentList
import numpy as np
import socket
from libtbx.mpi4py import MPI
COMM = MPI.COMM_WORLD

import logging

logger = logging.getLogger("dials.command_line.stills_process")


phil_str = """
include scope dials.command_line.stills_process.phil_scope
hopper {
  include scope simtbx.command_line.hopper.phil_scope
}
skip_hopper=False
  .type = bool
  .help = if True, then skip the hopper refinement, i.e. just run stills
  .help = process without refinement like usual
"""
from libtbx.phil import parse
phil_scope = parse(phil_str, process_includes=True)


class Hopper_Processor(Processor):

    def __init__(self, *args, **kwargs):
        super(Hopper_Processor, self).__init__(*args, **kwargs)
        if self.params.hopper.quiet:
            logging.getLogger("dials.algorithms.indexing.nave_parameters").setLevel(logging.ERROR)
            logging.getLogger("dials.algorithms.indexing.stills_indexer").setLevel(logging.ERROR)
            logging.getLogger("dials.algorithms.refinement.refiner").setLevel(logging.ERROR)
            logging.getLogger("dials.algorithms.refinement.reflection_manager").setLevel(logging.ERROR)
            logging.getLogger("dials.algorithms.refinement.reflection_manager").setLevel(logging.ERROR)

    def refine(self, exps, ref):
        exps_out = exps
        if not self.params.skip_hopper:
            if self.params.dispatch.refine:
                print("WARNING: hopper_process will always run its own refinement, ignoring dials.refine phil scope")
            self.params.dispatch.refine = False
            assert len(exps)==1
            # TODO MPI select GPU device

            if self.params.hopper.refiner.randomize_devices:
                dev = np.random.choice(self.params.hopper.refiner.num_devices)
                print("Rank %d will use random device %d on host %s" % (COMM.rank, dev, socket.gethostname()), flush=True)
            else:
                dev = COMM.rank % self.params.hopper.refiner.num_devices
                print("Rank %d will use fixed device %d on host %s" % (COMM.rank, dev, socket.gethostname()), flush=True)

            exp, ref = refine(exps[0], ref, self.params.hopper, gpu_device=dev)
            exps_out = ExperimentList()
            exps_out.append(exp)
        return super(Hopper_Processor, self).refine(exps_out, ref)


if __name__=="__main__":
    from dials.command_line import stills_process
    stills_process.Processor = Hopper_Processor
    stills_process.phil_scope = phil_scope
    stills_process.run()
