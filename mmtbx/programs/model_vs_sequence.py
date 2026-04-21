"""Verify the sequence of each chain in a model file against a sequence file"""
from __future__ import absolute_import, division, print_function
try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  from libtbx.program_template import ProgramTemplate
from libtbx.utils import Sorry
import mmtbx.validation.sequence

master_phil_str = """
include scope mmtbx.validation.sequence.master_phil
"""

# =============================================================================

class Program(ProgramTemplate):

  description = '''
Verify the sequence of each chain in a model file to detect residue mismatches
and other inconsistencies (similar to validation upon PDB deposition).

phenix.model_vs_sequence model.pdb sequence.fa
'''

  datatypes = ['model', 'sequence', 'phil']

  master_phil_str = master_phil_str

  # ---------------------------------------------------------------------------

  def validate(self):
    print('Validating inputs...\n', file=self.logger)
    self.data_manager.has_models(
      raise_sorry=True,
      expected_n=1,
      exact_count=True)
    if not self.data_manager.has_sequences():
      raise Sorry('No sequence file supplied.')

  # ---------------------------------------------------------------------------

  def run(self):
    model_fn = self.data_manager.get_default_model_name()
    seq_fn = self.data_manager.get_default_sequence_name()
    print('Using model:    %s' % model_fn, file=self.logger)
    print('Using sequence: %s' % seq_fn, file=self.logger)

    model = self.data_manager.get_model()
    sequences = self.data_manager.get_sequence()

    if len(sequences) == 0:
      raise Sorry("There don't appear to be any valid sequences in %s!" % seq_fn)

    self.result = mmtbx.validation.sequence.validation(
      pdb_hierarchy=model.get_hierarchy(),
      sequences=sequences,
      params=self.params,
      log=self.logger)
    self.result.show(out=self.logger)

  def get_results(self):
    return self.result
