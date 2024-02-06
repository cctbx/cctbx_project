from __future__ import absolute_import, division, print_function
<<<<<<< Updated upstream
from phenix.program_template import ProgramTemplate
from libtbx import group_args
from cctbx.maptbx.qscore import qscore_np

=======

import json
from pathlib import Path

from libtbx.program_template import ProgramTemplate
from libtbx import group_args
from cctbx.maptbx.qscore import (
    calc_qscore,
    calc_qscore_flex,
    cctbx_atoms_to_df,
    write_bild_spheres
)
>>>>>>> Stashed changes
import numpy as np


# =============================================================================

class Program(ProgramTemplate):

  description = """
  Perform a Qscore analysis for map-model fit
  """

  datatypes = ['phil', 'model', 'real_map']

<<<<<<< Updated upstream
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


=======
  master_phil_str = """
  include scope cctbx.maptbx.qscore.master_phil_str
>>>>>>> Stashed changes
  """

  def validate(self):
    assert self.params.qscore.backend in [
      "numpy","flex"
      ], "Provide one of 'numpy', 'flex'"

<<<<<<< Updated upstream
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
=======
    assert self.params.qscore.probe_allocation_method in [
      "progressive", "precalculate"
    ], "Provide one of 'progressive' or 'precalculate'"

  def run(self):
    print("Running")

    # get initial data
    mmm = self.data_manager.get_map_model_manager()

    # calculate shells
    if len(self.params.qscore.shells) ==0 :
        start = self.params.qscore.shell_radius_start
        stop = self.params.qscore.shell_radius_stop
        num = self.params.qscore.shell_radius_num
        shells = list(np.linspace(
        start,
        stop,
        num,
        endpoint=True))
        for shell in shells:
          self.params.qscore.shells.append(shell)


    # ignore hydrogens
    model = mmm.model()
    model = model.select(model.selection("not element H"))

    # make mmm
    mmm.set_model(model,overwrite=True)


    # run qscore
    backend = self.params.qscore.backend
    calc_func = calc_qscore if backend == "numpy" else calc_qscore_flex
    qscore_result= calc_func(
        mmm,
        self.params.qscore,
        log=self.logger)


    self.result = group_args(**qscore_result)

    # # save as dataframe
    # df = cctbx_atoms_to_df(model.get_atoms())
    # df["Q-score"] = self.result.qscore_per_atom
    # self.result.qscore_dataframe =



    self.write_results()

  def get_results(self):
    return self.result

  def get_results_as_JSON(self):
    results_dict = {
      "flat_results" : self.result.qscore_dataframe.to_dict(orient="records")
    }
    return json.dumps(results_dict,indent=2)


  def write_results(self):
    with open("qscore_results.json","w") as fh:
      fh.write(self.get_results_as_JSON())

    # write bild files
    if self.params.qscore.debug:
      print("Writing probe debug files...Using a small selection is recommended",
            file=self.logger)
      debug_path = Path("qscore_debug")
      debug_path.mkdir(exist_ok=True)
      for i,shell in enumerate(self.params.qscore.shells):
        shell = str(round(shell,2))
        probe_xyz = self.result.probe_xyz[i]
        n_shells, n_atoms,n_probes,_ = self.result.probe_xyz.shape
        probe_xyz_flat = probe_xyz.reshape((n_atoms*n_probes,3))
        out_file = Path(debug_path,f"probes_shell_{shell}.bild")
        write_bild_spheres(probe_xyz_flat,str(out_file),r=0.2)
>>>>>>> Stashed changes
