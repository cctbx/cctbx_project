# -*- coding: utf-8 -*-
from __future__ import division, print_function

from libtbx.program_template import ProgramTemplate

from mmtbx.validation import evalurama

# =============================================================================

class Program(ProgramTemplate):

  description = '''
phenix.evalurama: tool for checking distribution of points inside favored
  ramachandran region

Usage examples:
  phenix.evalurama model1.pdb
  '''

  datatypes = ['model', 'phil']

  master_phil_str = """\
    include scope mmtbx.validation.evalurama.master_phil_str
"""

  # ---------------------------------------------------------------------------
  def validate(self):
    print('Validating inputs', file=self.logger)
    self.data_manager.has_models(expected_n=1, exact_count=True, raise_sorry=True)
    m = self.data_manager.get_model()
    assert m.get_hierarchy().models_size() == 1

  # ---------------------------------------------------------------------------
  def run(self):
    model = self.data_manager.get_model()

    self.evalurama = evalurama.eval(
        model = model,
        params = self.params.evalurama,
        log = self.logger)

    results = self.evalurama.get_results()
    print("Grand total results:")
    for k, v in results.iteritems():
      if v is not None:
        print(k,v, file=self.logger)

  # ---------------------------------------------------------------------------
  def get_results(self):
    return self.evalurama.get_results()

