from __future__ import absolute_import, division, print_function
from phenix.program_template import ProgramTemplate
from libtbx import group_args
from cctbx.maptbx.qscore import qscore_np

import numpy as np


# =============================================================================

class Program(ProgramTemplate):

    description = """
  Perform a Qscore analysis for map-model fit
  """

    datatypes = ['phil', 'model', 'real_map']

    master_phil_str = """
  nproc = 8
      .type = int
      .help = Number of processors to use
      .short_caption = Number of processors to use
      .expert_level = 1
  n_probes = 32
      .type = int
      .help = Number of radial probes to use
      .short_caption = Number of radial probes to use
      .expert_level = 1
  selection = None
    .type = str
    .help = Only test atoms within this selection
    .short_caption = Only test atoms within this selection
    .expert_level = 1

  shell_radius_start = 0.1
    .type = float
    .help = Start testing density at this radius from atom
    .short_caption = Start testing density at this radius from atom
    .expert_level = 1

  shell_radius_stop = 2
    .type = float
    .help = Stop testing density at this radius from atom
    .short_caption = Stop testing density at this radius from atom
    .expert_level = 1

  shell_radius_num = 20
    .type = int
    .help = The number of radial shells
    .short_caption = The number of radial shells (includes start/stop, so minimum 2)
    .expert_level = 1

  probe_allocation_method = precalculate
    .type = str
    .help = The method used to allocate radial probes
    .short_caption = Either 'progressive' or 'precalculate'. Progressive is the original method \
                     where probes are proposed and rejected iteratively. \
                     Precalculate is a method where probes are pre-allocated and \
                     rejected once. Parallelization is done by radial shell. \
                     Precalculate is much faster but will yield slightly different results.


  """

    def validate(self):
        pass

    def run(self):
        print("Running")
        shells = np.linspace(self.params.shell_radius_start,
                             self.params.shell_radius_stop,
                             num=self.params.shell_radius_num,
                             endpoint=True)
        mmm = self.data_manager.get_map_model_manager()
        version = 2 if self.params.probe_allocation_method == "precalculate" else 1  # 'progressive'

        # TODO: move param unpacking to the library function
        result = qscore_np(
            mmm,
            selection=self.params.selection,
            n_probes=self.params.n_probes,
            shells=shells,
            nproc=self.params.nproc,
            version=version,
            log=self.logger)
        print("Result:")
        for val in result:
            print(str(val)+",")
        self.result = group_args(qscore=result)

    def get_results(self):
        return self.result
