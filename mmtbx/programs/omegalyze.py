"""Analyze omega angles in a model, report on cis peptides"""
from __future__ import absolute_import, division, print_function

import os
from mmtbx.validation.omegalyze import omegalyze
from libtbx.program_template import ProgramTemplate
from libtbx.utils import Sorry
from datetime import datetime

class Program(ProgramTemplate):
  prog = os.getenv('LIBTBX_DISPATCHER_NAME')
  description="""
%(prog)s file.pdb [params.eff] [options ...]

Options:

  model=input_file      input PDB or mmCIF file
  nontrans_only=True    only print nontrans residues (does not affect kinemage)
  text=True          verbose colon-delimited text output (default)
  kinemage=False        Create kinemage markup (overrides text output)
  json=False            Outputs results as JSON compatible dictionary
  help = False          Prints this help message if true

  text output is colon-delimited and follows the format:
    residue:type:omega:conformation
      'residue' is a unique residue identifier
      'type' is either proline or general case
      'omega' is the calculated omega dihedral for the peptide between this
        residue and the preceeding residue
      'conformation' is: cis for omega within 30 degrees of planar cis
                         trans for omega within 30 degrees of planar trans
                         twisted for omega not within 30 degrees of planar

  SUMMARY statistics provide:
    counts of cis prolines and twisted prolines relative to total prolines with
      measurable omega dihedrals across all chains
    counts of non-proline cis and twisted peptides relative to total non-proline
      peptides with measurable omega dihedrals across all chains

  Cis Prolines occur in ~5%% of prolines (1 in 20) at high resolution
  Non-Proline Cis residues occur in ~0.05%% of residues (1 in 2000) and require
    clear support from experimental data or homology.
  Twisted peptides are even less frequent and are highly suspect without
    high-resolution data.

Example:

  %(prog)s model=1ubq.pdb kinemage=True
""" % locals()

  master_phil_str = """
    nontrans_only = True
      .type = bool
      .help = "Controls whether trans peptides are stored and printed"
    text = True
      .type = bool
      .help = "Prints verbose, colon-delimited text output and summary"
    kinemage = False
      .type = bool
      .help = "Prints kinemage markup for cis-peptides"
    json = False
      .type = bool
      .help = "Prints results as JSON format dictionary"
    oneline = False
      .type = bool
      .help = "Prints oneline-style summary statistics"
    result_file = None
      .type = path
      .help = "Path for output file"
"""
  datatypes = ['model','phil']
  data_manager_options = ['model_skip_expand_with_mtrix']

  def validate(self):
    self.data_manager.has_models(raise_sorry=True)

  def run(self):
    f = None
    if self.params.result_file is not None:
      try:
        f = open(self.params.result_file, 'w')
        self.logger.register('file', f)
      except IOError:
        raise Sorry("The output file could not be opened")
    hierarchy = self.data_manager.get_model().get_hierarchy()
    self.info_json = {"model_name":self.data_manager.get_default_model_name(),
                      "time_analyzed": str(datetime.now())}

    self.results = omegalyze(
      pdb_hierarchy=hierarchy,
      nontrans_only=self.params.nontrans_only,
      out=self.logger,
      quiet=False)
    if self.params.kinemage:
      print(self.results.as_kinemage(), file=self.logger)
    elif self.params.oneline:
      self.results.summary_only(out=self.logger, pdbid=self.data_manager.get_default_model_name())#params.model)
    elif self.params.json:
      print(self.results.as_JSON(addon_json=self.info_json), file=self.logger)
    elif self.params.text:
      self.results.show_old_output(out=self.logger, verbose=True)
    if f:
      try:
        f.close()
      except Exception:
        raise Sorry("Could not close output file")

  def get_results(self):
    return self.results

  def get_results_as_JSON(self):
    return self.results.as_JSON(self.info_json)
