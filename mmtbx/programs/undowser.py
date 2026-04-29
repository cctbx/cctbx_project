"""Identify unlikely waters in a model based on clashes"""
from __future__ import absolute_import, division, print_function

import os
from mmtbx.validation.undowser import undowserlyze
from libtbx.program_template import ProgramTemplate
from datetime import datetime

class Program(ProgramTemplate):
  prog = os.getenv('LIBTBX_DISPATCHER_NAME')
  description="""\
%(prog)s file.pdb [params.eff] [options ...]

Options:

  model=input_file        input PDB file
  outliers_only=False   only print outliers
  keep_hydrogens=False      keep input hydrogen atoms if True, regenerate if False
  nuclear=False             use nuclear x-H distances and vdW radii
  json=False            Outputs results as JSON compatible dictionary
  verbose=False         verbose text output

Example:

  %(prog)s model=1ubq.pdb outliers_only=True
""" % locals()

  master_phil_str = """
  include scope mmtbx.validation.molprobity_cmdline_phil_str
  show_errors = False
    .type = bool
    .help = '''Print out errors'''
  keep_hydrogens = False
    .type = bool
    .help = '''Keep hydrogens in input file'''
  nuclear = False
    .type = bool
    .help = '''Use nuclear hydrogen positions'''
  time_limit = 120
    .type = int
    .help = '''Time limit (sec) for Reduce optimization'''
  json = False
    .type = bool
    .help = "Prints results as JSON format dictionary"
  use_parent = False
    .type = bool
  """
  datatypes = ['model','phil']
  data_manager_options = ['model_skip_expand_with_mtrix']
  known_article_ids = ['molprobity']

  def validate(self):
    self.data_manager.has_models(raise_sorry=True)

  def run(self):
    hierarchy = self.data_manager.get_model().get_hierarchy()
    self.info_json = {"model_name":self.data_manager.get_default_model_name(),
                      "time_analyzed": str(datetime.now())}
    self.results = undowserlyze(
      pdb_hierarchy=hierarchy,
      keep_hydrogens=self.params.keep_hydrogens,
      nuclear=self.params.nuclear,
      outliers_only=self.params.outliers_only,)
    if self.params.json:
      print(self.results.as_JSON(), file=self.logger)
    else:
      print(self.results.as_HTML())

  def get_results(self):
    return self.results

  def get_results_as_JSON(self):
    return self.results.as_JSON(self.info_json)
