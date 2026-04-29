"""Check bond and angle values in a model"""
from __future__ import absolute_import, division, print_function

import os
from mmtbx.model import manager
from libtbx.program_template import ProgramTemplate
from libtbx.utils import null_out
import json
from mmtbx.validation.mp_validate_bonds import mp_validate_bonds
from datetime import datetime

class Program(ProgramTemplate):
  prog = os.getenv('LIBTBX_DISPATCHER_NAME')
  description="""\
%(prog)s file.pdb [params.eff] [options ...]

Options:

  model=input_file        input PDB file
  outliers_only=False   only print outliers
  json=False            Outputs results as JSON compatible dictionary
  use_cdl=True          Use the Conformational Dependent Library (cdl) for reference values
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
  use_cdl = True
    .type = bool
    .help = "Use conformational dependent library for reference values"
  use_parent = False
    .type = bool
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
                      "time_analyzed": str(datetime.now()),
                      "params": {"use_cdl":self.params.use_cdl,
                                 "outliers_only":self.params.outliers_only,
                                 }}
    p = manager.get_default_pdb_interpretation_params()
    ##print(dir(p.pdb_interpretation))
    p.pdb_interpretation.allow_polymer_cross_special_position=True
    p.pdb_interpretation.flip_symmetric_amino_acids=False
    p.pdb_interpretation.clash_guard.nonbonded_distance_threshold = None
    p.pdb_interpretation.restraints_library.cdl=self.params.use_cdl
    model.set_log(log = null_out())
    model.process(make_restraints=True, pdb_interpretation_params=p)
    geometry = model.get_restraints_manager().geometry
    atoms = hierarchy.atoms()
    self.results = mp_validate_bonds(
      pdb_hierarchy=hierarchy,
      geometry_restraints_manager=geometry,
      outliers_only=self.params.outliers_only)
    if self.params.json:
      print(self.results.as_JSON(addon_json=self.info_json), file=self.logger)
    elif self.params.verbose:
      self.results.show(out=self.logger, verbose=True)

  def get_results(self):
    return self.results

  def get_results_as_JSON(self):
    return self.results.as_JSON(self.info_json)
