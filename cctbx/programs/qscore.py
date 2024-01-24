from __future__ import absolute_import, division, print_function
from libtbx.program_template import ProgramTemplate
from libtbx import group_args
from cctbx.maptbx import qscore

import numpy as np

# =============================================================================

class Program(ProgramTemplate):

    description = """
  Perform a Qscore analysis for map-model fit
  """

    datatypes = ['phil', 'model', 'real_map']

    master_phil_str = """
    include scope cctbx.maptbx.qscore.master_phil_str
    """

    def validate(self):
        pass

    def run(self):
        print("Running")
        mmm = self.data_manager.get_map_model_manager()
        if self.params.qscore.selection != None:
            selection = np.where(mmm.model().selection(self.params.qscore.selection).as_numpy_array())[0]
            if len(selection)==0:
                print("Finished... nothing selected")
                self.result = group_args()
                return
            self.params.qscore.selection = selection
        qscore_per_atom= qscore.run_qscore(
            mmm,
            self.params.qscore,
            log=self.logger)

        print("Finished")
        self.result = group_args(qscore_per_atom=qscore_per_atom)


    def get_results(self):
        return self.result
