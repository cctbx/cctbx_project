"""Validate RNA"""
from __future__ import absolute_import, division, print_function

import os
from mmtbx.model import manager
from mmtbx.validation.rna_validate import rna_validation
from libtbx.program_template import ProgramTemplate
from libtbx.utils import null_out
from datetime import datetime

class Program(ProgramTemplate):
  prog = os.getenv('LIBTBX_DISPATCHER_NAME')
  description="""\
%(prog)s file.pdb [params.eff] [options ...]

Options:

  model=input_file        input PDB file
  outliers_only=False   only print outliers
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
  json = False
    .type = bool
    .help = "Prints results as JSON format dictionary"
  use_parent = False
    .type = bool
  rna_sugar_pucker_analysis
    .short_caption = RNA sugar pucker analysis
    .style = box noauto auto_align menu_item parent_submenu:advanced
  {
    include scope mmtbx.monomer_library.rna_sugar_pucker_analysis.master_phil
  }
  """
  datatypes = ['model','phil']
  data_manager_options = ['model_skip_expand_with_mtrix']
  known_article_ids = ['molprobity']

  def validate(self):
    self.data_manager.has_models(raise_sorry=True)

  def run(self):
    model = self.data_manager.get_model()
    model.set_stop_for_unknowns(False)
    hierarchy = model.get_hierarchy()
    self.info_json = {"model_name":self.data_manager.get_default_model_name(),
                      "time_analyzed": str(datetime.now())}
    p = manager.get_default_pdb_interpretation_params()
    ##print(dir(p.pdb_interpretation))
    p.pdb_interpretation.allow_polymer_cross_special_position=True
    p.pdb_interpretation.flip_symmetric_amino_acids=False
    p.pdb_interpretation.clash_guard.nonbonded_distance_threshold = None
    model.set_log(log = null_out())
    model.process(make_restraints=True, pdb_interpretation_params=p)
    geometry = model.get_restraints_manager().geometry

    self.results = rna_validation(
      pdb_hierarchy=hierarchy,
      geometry_restraints_manager=geometry,
      params=self.params,
      outliers_only=self.params.outliers_only)
    if self.params.json:
      print(self.results.as_JSON(), file=self.logger)
    elif self.params.verbose:
      self.results.show(out=self.logger, verbose=True)

  def get_results(self):
    return self.results

  def get_results_as_JSON(self):
    return self.results.as_JSON(self.info_json)
