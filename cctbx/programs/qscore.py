from __future__ import absolute_import, division, print_function
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
    assert self.params.qscore.backend in [
      "numpy","flex"
      ], "Provide one of 'numpy', 'flex'"

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
        out_file = Path(debug_path,"probes_shell_"+shell+".bild")
        write_bild_spheres(probe_xyz_flat,str(out_file),r=0.2)
