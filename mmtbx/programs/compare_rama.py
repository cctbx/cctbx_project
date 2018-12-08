# -*- coding: utf-8 -*-
from __future__ import division, print_function

from libtbx.program_template import ProgramTemplate

from mmtbx.validation import compare_rama

# =============================================================================

class Program(ProgramTemplate):

  description = '''
phenix.compare_rama: tool for compare Ramachandran plots, e.g. before-after
  refinement.

Usage examples:
  phenix.compare_rama model1.pdb model2.pdb
  phenix.compare_rama model1.cif model2.cif
  '''

  datatypes = ['model', 'phil']

  master_phil_str = compare_rama.master_phil_str

  # ---------------------------------------------------------------------------
  def validate(self):
    print('Validating inputs', file=self.logger)
    self.data_manager.has_models(expected_n=2, exact_count=True, raise_sorry=True)
    model_1, model_2 = self._get_models()
    assert model_1.get_hierarchy().is_similar_hierarchy(model_2.get_hierarchy())
    for m in [model_1, model_2]:
      assert m.get_hierarchy().models_size() == 1

  # ---------------------------------------------------------------------------
  def run(self):
    # I'm guessing self.data_manager, self.params and self.logger
    # are already defined here...
    # print('Using model: %s' % self.data_manager.get_default_model_name(), file=self.logger)
    print(dir(self.data_manager))
    print(self.data_manager.get_model_names())
    # this must be mmtbx.model.manager?
    model_1, model_2 = self._get_models()

    self.rama_comp = compare_rama.rcompare(
        model1 = model_1,
        model2 = model_2,
        params = self.params.compare_rama,
        log = self.logger)

  # ---------------------------------------------------------------------------
  def get_results(self):
    return self.rama_comp.get_results()

  def _get_models(self):
    m_names = self.data_manager.get_model_names()
    model_1 = self.data_manager.get_model(filename=m_names[0])
    model_2 = self.data_manager.get_model(filename=m_names[1])
    return model_1, model_2

# =============================================================================
# end
