"""Tool for sampling different conformations in attempt to
  fix Cablam outliers"""
# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function

from libtbx.program_template import ProgramTemplate

from mmtbx.building import cablam_idealization

import os

# =============================================================================

class Program(ProgramTemplate):

  description = '''
phenix.cablam_idealization: tool for sampling different conformations in attempt to
  fix Cablam outliers.

Usage examples:
  phenix.cablam_idealization model.pdb
  phenix.cablam_idealization model.cif
  '''

  datatypes = ['model', 'phil']

  master_phil_str = """
include scope mmtbx.building.cablam_idealization.master_phil_str
output {
  suffix = _cablam_fixed
    .type = str
}
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
    model = self.data_manager.get_model()

    self.output_fname_base = os.path.splitext(
        self.data_manager.get_default_model_name())[0] + self.params.output.suffix
    fo = open(self.output_fname_base+'.log', 'w')
    self.logger.register(label='logfile', file_object=fo)

    self.cablam_id = cablam_idealization.cablam_idealization(
        model=model,
        params = self.params.cablam_idealization,
        log = self.logger)

    results = self.cablam_id.get_results()
    print("Total number of tried outliers: %d" % results.n_tried_residues, file=self.logger)
    print("Number of rotated outliers: %d" % results.n_rotated_residues, file=self.logger)
    # I believe this should go to data_manager. Also not clear how output of
    # two files would affect data_manager.
    for m, fname_base in [
        (results.model, self.output_fname_base),
        (results.model_minimized, self.output_fname_base+"_minimized")]:
      if m is not None:
        self.final_file_name = self.data_manager.write_model_file(
            m, fname_base)
        print("Model written to '%s'" % self.final_file_name)

  # ---------------------------------------------------------------------------
  def get_results(self):
    return self.cablam_id.get_results()

# =============================================================================
# end
