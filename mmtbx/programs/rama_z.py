# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function

from libtbx.program_template import ProgramTemplate

from mmtbx.validation import rama_z
from libtbx.utils import Sorry

# =============================================================================

class Program(ProgramTemplate):

  description = '''
mmtbx.rama_z: Tool to calculate Rama-Z score. Validation of Ramachandran plot.

Usage examples:
  mmtbx.rama_z model1.pdb
  '''

  datatypes = ['model', 'phil']

  master_phil_str = """\
  write_HSL_models = False
    .type = bool
  write_HSL_plot = False
    .type = bool
  write_HSL_general_only = True
    .type = bool
  write_whole_plot = False
    .type = bool
  write_whole_general_only = True
    .type = bool
"""

  # ---------------------------------------------------------------------------
  def validate(self):
    print('Validating inputs', file=self.logger)
    # print(dir(self.data_manager))
    # print("def output", self.data_manager.get_default_output_filename())
    # print("def output", self.data_manager.get_default_output_model_filename())
    # STOP()
    self.data_manager.has_models(expected_n=1, exact_count=True, raise_sorry=True)
    m = self.data_manager.get_model()
    if m.get_hierarchy().models_size() != 1:
      raise Sorry("Multi-model files are not supported.")

  # ---------------------------------------------------------------------------
  def run(self):

    model = self.data_manager.get_model()
    # Temporarly disable
    # self.params.output.prefix=self.data_manager.get_default_model_name().split('.')[0]
    # self.params.output.suffix='_helix'
    # self.data_manager.set_default_output_filename(self.get_default_output_filename())

    self.rama_z = rama_z.rama_z(
        model = model,
        log = self.logger)

    # self.data_manager.write_model_file(model)
    # self.params.output.suffix='_sheet'
    # self.data_manager.set_default_output_filename(self.get_default_output_filename())
    # self.data_manager.write_model_file(model)

    result = self.get_results()
    if result is None:
      print("Calculation of z-score failed for some reason", file=self.logger)
    else:
      for k in ["whole", "helix", "sheet", "loop"]:
        rc = k[0].upper()
        v = result.get(rc, None)
        if v is None:
          print("z-score %-5s: None, residues: %d" % (k, result['residue_counts'][rc]), file=self.logger)
        else:
          print("z-score %-5s: %6.3f (%5.3f), residues: %d" % (k, v[0], v[1], result['residue_counts'][rc]), file=self.logger)

  # ---------------------------------------------------------------------------
  def get_results(self):
    r = self.rama_z.get_z_scores()
    r['residue_counts'] = self.rama_z.get_residue_counts()
    return r
