"""Calculates ideal CB atom position based on mainchain geometry,
then displays deviation of modeled CB from that ideal positon.
"""
from __future__ import absolute_import, division, print_function

import os
from mmtbx.validation.cbetadev import cbetadev
from libtbx.program_template import ProgramTemplate
#from libtbx.utils import Sorry
from datetime import datetime

class Program(ProgramTemplate):
  prog = os.getenv('LIBTBX_DISPATCHER_NAME')
  description="""
  %(prog)s file.pdb [params.eff] [options ...]

Calculates ideal CB atom position based on mainchain geometry,
then displays deviation of modeled CB from that ideal positon.
Deviations of >= 0.25A are considered outliers

Options:

  model=input_file      input PDB file
  output=text, kin, or bullseye    select type of output
  json=False            Outputs results as JSON compatible dictionary
  outliers_only=False   suppress non-outlier results

Example:

  %(prog)s model=1ubq.pdb
""" % locals()

  master_phil_str = """
  include scope mmtbx.validation.molprobity_cmdline_phil_str
    json = False
      .type = bool
      .help = "Prints results as JSON format dictionary"
    cbetadev {
      output = *text kin bullseye
        .type = choice
        .help = '''choose output type'''
      apply_phi_psi_correction = False
        .type = bool
        .help = XXX
      display_phi_psi_correction = False
        .type = bool
        .help = XXX
      exclude_d_peptides = False
        .type = bool
        .style = hidden
        .help = Attempts to exclude D-peptide using the large CBD
      }
"""
  datatypes = ['model','phil']
  data_manager_options = ['model_skip_expand_with_mtrix']
  known_article_ids = ['molprobity']

  def validate(self):
    self.data_manager.has_models(raise_sorry=True)

  def get_results(self):
    return self.results

  def get_results_as_JSON(self):
    return self.results.as_JSON(self.info_json)

  def run(self):
    hierarchy = self.data_manager.get_model().get_hierarchy()
    self.info_json = {"model_name":self.data_manager.get_default_model_name(),
                      "time_analyzed": str(datetime.now())}
    self.results = cbetadev(
      pdb_hierarchy=hierarchy,
      outliers_only=self.params.outliers_only,
      apply_phi_psi_correction=self.params.cbetadev.apply_phi_psi_correction,
      display_phi_psi_correction=self.params.cbetadev.display_phi_psi_correction,
      exclude_d_peptides=self.params.cbetadev.exclude_d_peptides,
      out=self.logger,
      quiet=False)
    if self.params.cbetadev.output == "kin":
      self.logger.write(self.results.as_kinemage())
    elif self.params.cbetadev.output == "bullseye":
      filebase = os.path.basename(self.data_manager.get_model_names()[0])
      self.logger.write(self.results.as_bullseye_kinemage(pdbid=filebase))
    elif self.params.json:
      print(self.get_results_as_JSON())
    elif self.params.verbose:
      #pdb_file_str = os.path.basename(self.params.model)[:-4]
      #get input file name from data manager, strip file extension
      pdb_file_str = os.path.basename(self.data_manager.get_model_names()[0])[:-4]
      self.results.show_old_output(out=self.logger, prefix=pdb_file_str, verbose=True)
