from __future__ import absolute_import, division, print_function
import json
from pathlib import Path

from cctbx.array_family import flex
from libtbx.program_template import ProgramTemplate
from libtbx import group_args
from cctbx.maptbx.qscore import (
    calc_qscore,
    cctbx_atoms_to_df,
    write_bild_spheres
)
from libtbx.utils import Sorry
import numpy as np
import pandas as pd

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
    # test for sane parameters
    if not (4<=self.params.qscore.n_probes<=512):
      raise Sorry("Provide n_probe values in the range 4-512")

    if  not (4<=self.params.qscore.shell_radius_num<=128):
      raise Sorry("Provide shell_radius_num values in range 4-128")



  def run(self):
    # set up logging
    log = open("qscore.log", "w")

    self.logger.register("log", log, atexit_send_to=None)

    # get initial data
    mmm = self.data_manager.get_map_model_manager()

    # calculate shells
    self.shells = []
    # add a range of shells
    start = self.params.qscore.shell_radius_start
    stop = self.params.qscore.shell_radius_stop
    num = self.params.qscore.shell_radius_num
    shells = list(np.linspace(start,stop,num,endpoint=True))

    for shell in reversed(shells):
      self.shells.insert(0,shell)


    # ignore hydrogens
    model = mmm.model()
    model.add_crystal_symmetry_if_necessary()
    model = model.remove_hydrogens()

    # make mmm
    mmm.set_model(model,overwrite=True)

    # print output
    self._print("Running Q-score:")
    self._print("\nRadial shells (Radii used to place probes):")
    self._print(",".join([str(round(shell,2)) for shell in shells]))
    self._print("\n Number of probes per radial shell: "+str(self.params.qscore.n_probes))

    # run qscore
    qscore_result= calc_qscore(
        mmm,
        selection=self.params.qscore.selection,
        n_probes=self.params.qscore.n_probes,
        rtol=self.params.qscore.rtol,
        shells=self.shells,
        nproc=self.params.qscore.nproc,
        log=self.logger)


    self.result = group_args(**qscore_result)
    # calculate some metrics
    df = self.result.qscore_dataframe
    # round
    df_numeric = df.select_dtypes(include=['number']).round(2)
    df[df_numeric.columns] = df_numeric

    if self.params.qscore.selection is not None:
      model = model.select(model.selection(self.params.qscore.selection))
    assert model.get_number_of_atoms()==len(df)


    self._print("\nFinished running.\n\n Q-score results:")
    sel_mc = "protein and (name C or name N or name CA or name O or name CB)"
    sel_mc = model.selection(sel_mc)
    sel_sc = ~sel_mc

    # Side chains
    df["Side Chain"] = pd.NA
    if sel_sc.count(True)>0:
      q_sc = flex.mean(self.result.qscore_per_atom.select(sel_sc))
      q_sc = round(q_sc,2)
      df.loc[sel_sc.as_numpy_array(),"Side Chain"] = np.array(self.result.qscore_per_atom.select(sel_sc))
    else:
      q_sc = None
    df["Side Chain"] = df["Side Chain"].astype("Float64")

    # Main chain
    df["Main Chain"] = pd.NA
    if sel_mc.count(True)>0:
      q_mc = flex.mean(self.result.qscore_per_atom.select(sel_mc))
      q_mc = round(q_mc,2)
      df.loc[sel_mc.as_numpy_array(),"Main Chain"] = np.array(self.result.qscore_per_atom.select(sel_mc))
    else:
      q_mc = None
    df["Main Chain"] = df["Main Chain"].astype("Float64")
    # All
    q_all = flex.mean(self.result.qscore_per_atom)
    q_all = round(q_all,2)

    # Per chain
    q_chains = df.groupby("chain_id").agg('mean',numeric_only=True)
    q_chains.drop(columns=["id","resseq","x","y","z","Q-Residue"],inplace=True)
    q_chains = q_chains.round(2)
    self._print("\n\nBy residue:")
    self._print("----------------------------------------")
    pd.set_option('display.max_rows', 20)
    df.drop(columns=["x","y","z","id","Q-score"],inplace=True)
    self._print(df.groupby(["model_id","chain_id","resseq","altloc"]).agg("mean"))
    self._print("\n\nBy chain:")
    self._print("----------------------------------------")
    self._print(q_chains)
    self._print("\n\nAll:")
    self._print("----------------------------------------")
    df = pd.DataFrame({"Main Chain":[q_mc],"Side Chain":[q_sc],"Overall":[q_all]})
    self._print(df)


    # store in results
    self.result.q_score_chain_df = q_chains
    self.result.q_score_side_chain = q_sc
    self.result.q_score_main_chain = q_mc
    self.result.q_score_overall = q_all

    # write out
    if self.params.qscore.write_probes:
      self.write_bild_spheres()

    if self.params.qscore.write_to_bfactor_pdb:
      self.write_to_bfactor_pdb(model,self.result.qscore_per_atom)


  def get_results(self):
    return self.result

  def get_results_as_JSON(self):
    results_dict = {
      "flat_results" : self.result.qscore_dataframe.to_dict(orient="records")
    }
    return json.dumps(results_dict,indent=2)

  def write_to_bfactor_pdb(self,model,qscore_per_atom):
    model.set_b_iso(qscore_per_atom)

    with open("qscore_bfactor_field.pdb","w") as fh:
      fh.write(model.model_as_pdb())

  def write_bild_spheres(self):
    # write bild files
    if self.params.qscore.write_probes:
      self._print("Writing probe debug files...Using a small selection is recommended",
            file=self.logger)
      debug_path = Path("qscore_debug")
      debug_path.mkdir(exist_ok=True)
      for i,shell in enumerate(self.shells):
        shell = str(round(shell,2))
        probe_xyz = self.result.probe_xyz[i]
        n_shells, n_atoms,n_probes,_ = self.result.probe_xyz.shape
        probe_xyz_flat = probe_xyz.reshape((n_atoms*n_probes,3))
        out_file = Path(debug_path,"probes_shell_"+shell+".bild")
        write_bild_spheres(probe_xyz_flat,str(out_file),r=0.2)
