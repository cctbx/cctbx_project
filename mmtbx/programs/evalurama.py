# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function

from libtbx.program_template import ProgramTemplate

from mmtbx.validation import evalurama
from mmtbx.validation import ramalyze
from collections import OrderedDict

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
    self.data_manager.has_models(expected_n=1, exact_count=False, raise_sorry=True)
    m = self.data_manager.get_model()
    assert m.get_hierarchy().models_size() == 1

  # ---------------------------------------------------------------------------
  def run(self):
    models = []
    for name in self.data_manager.get_model_names():
      models.append(self.data_manager.get_model(name))

    self.evalurama = evalurama.eval(
        models = models,
        params = self.params.evalurama,
        log = self.logger)

    results = self.evalurama.get_results()
    res_sorted = OrderedDict(sorted(results.items()))
    print("", file=self.logger)
    print("Grand total results:", file=self.logger)
    print("Rama type     Peak coordinates              CC   Number of residues", file=self.logger)
    for k, v in res_sorted.iteritems():
      if v is not None:
        print("%10s %-30s %.3f %5d" % (ramalyze.res_type_labels[k[0]],k[1], v[0], v[1]), file=self.logger)

  # ---------------------------------------------------------------------------
  def get_results(self):
    return self.evalurama.get_results()

