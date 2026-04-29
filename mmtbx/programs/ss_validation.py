"""Tool for validation of secondary structure annotations"""
# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function

from libtbx.program_template import ProgramTemplate

from mmtbx.secondary_structure import ss_validation

# =============================================================================

class Program(ProgramTemplate):

  description = '''
phenix.secondary_structure_validation: tool for validation of secondary
  structure annotations.

Usage examples:
  phenix.secondary_structure_validation model.pdb
  phenix.secondary_structure_validation model.cif
  phenix.secondary_structure_validation model.pdb nproc=7
  '''

  datatypes = ['model', 'phil']

  master_phil_str = ss_validation.master_phil_str

  # ---------------------------------------------------------------------------
  def validate(self):
    print('Validating inputs', file=self.logger)
    self.data_manager.has_models(raise_sorry=True)

  # ---------------------------------------------------------------------------
  def run(self):
    # I'm guessing self.data_manager, self.params and self.logger
    # are already defined here...
    print('Using model: %s' % self.data_manager.get_default_model_name(), file=self.logger)

    # this must be mmtbx.model.manager?
    model = self.data_manager.get_model()

    self.val_obj = ss_validation.validate(
        model=model,
        params = self.params.ss_validation,
        log = self.logger)

  # ---------------------------------------------------------------------------
  def get_results(self):
    return self.val_obj.get_results()

# =============================================================================
# end
