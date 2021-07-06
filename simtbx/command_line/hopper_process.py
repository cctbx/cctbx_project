from __future__ import absolute_import, division, print_function

import os

# LIBTBX_SET_DISPATCHER_NAME simtbx.diffBragg.hopper_process

from dials.command_line.stills_process import Processor
from libtbx.utils import Abort, Sorry
from simtbx.diffBragg.hopper_utils import refine
from dxtbx.model import Experiment, ExperimentList
from dxtbx.model.experiment_list import ExperimentListFactory

import logging
logger = logging.getLogger("dials.command_line.stills_process")


phil_str = """
include scope dials.command_line.stills_process.phil_scope
hopper {
  include scope simtbx.command_line.hopper.phil_scope
}
"""
from libtbx.phil import parse
phil_scope = parse(phil_str, process_includes=True)



class Hopper_Processor(Processor):

    def refine(self, exps, ref):
        if self.params.dispatch.refine:
            print("WARNING: hopper_process will always run its own refinement, ignoring dials.refine phil scope")
        self.params.dispatch.refine = False
        assert len(exps)==1
        # TODO MPI select GPU device
        spec='/global/cfs/cdirs/m3562/der/braggnanimous/7534_new_lams_p05/run802_shot343_newlam.lam'
        exp, ref = refine(exps[0], ref, self.params.hopper, spec=spec, gpu_device=0)
        El =ExperimentList()
        El.append(exp)
        return super(Hopper_Processor, self).refine(El, ref)


if __name__=="__main__":
    from dials.command_line import stills_process
    stills_process.Processor = Hopper_Processor
    stills_process.phil_scope = phil_scope
    stills_process.run()
