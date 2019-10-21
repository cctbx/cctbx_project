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
    include scope mmtbx.validation.rama_z.master_phil_str
"""

  # ---------------------------------------------------------------------------
  def validate(self):
    print('Validating inputs', file=self.logger)
    self.data_manager.has_models(expected_n=1, exact_count=True, raise_sorry=True)
    m = self.data_manager.get_model()
    if m.get_hierarchy().models_size() != 1:
      raise Sorry("Multi-model files are not supported.")

  # ---------------------------------------------------------------------------
  def run(self):
    model = self.data_manager.get_model()

    self.rama_z = rama_z.rama_z(
        model = model,
        params = self.params.rama_z,
        log = self.logger)

  # ---------------------------------------------------------------------------
  def get_results(self):
    r = self.rama_z.get_z_scores()
    r['residue_counts'] = self.rama_z.get_residue_counts()
    return r

