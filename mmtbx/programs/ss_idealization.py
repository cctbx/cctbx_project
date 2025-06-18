"""Tool for idealization of secondary structure"""
# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function

from libtbx.program_template import ProgramTemplate
from mmtbx.secondary_structure.build import ss_idealization
from libtbx import Auto
import os



# =============================================================================

class Program(ProgramTemplate):

  description = '''
phenix.ss_idealization: tool for idealization of secondary structure.
  Uses SS annotation in the input model file.

Usage examples:
  phenix.ss_idealization model.pdb
  phenix.ss_idealization model.cif
  phenix.ss_idealization model.pdb nproc=7
  '''

  datatypes = ['model', 'phil']

  master_phil_str = """
include scope mmtbx.secondary_structure.build.ss_idealization.ss_idealization_master_phil_str
  """

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
    self.model = self.data_manager.get_model()

    # Done in place!
    self.params.ss_idealization.enabled=True
    sss = ss_idealization.substitute_ss(
        model = self.model,
        params = self.params.ss_idealization,
        log=self.logger)
    sss.run()
    inp_fn = os.path.basename(self.data_manager.get_default_model_name())[:-4]
    fn = "%s" % self.get_default_output_filename(
                prefix='%s_' % inp_fn,
                suffix='ss_idealized',
                serial=Auto)
    actual_fn = self.data_manager.write_model_file(self.model, filename=fn)

  # ---------------------------------------------------------------------------
  def get_results(self):
    return self.model

# =============================================================================
# end
