from __future__ import absolute_import, division, print_function

import os

from libtbx import group_args
from libtbx.program_template import ProgramTemplate

# =============================================================================
class Program(ProgramTemplate):

  description = '''
Program for preparing model and data files for depostion into the Proten
Data Bank.

Minimum required data:
  Model file
  Sequence file

The sequence file should have a sequence for each chain in the model file.

Currently, this program only combines the model and sequence into a single
mmCIF file. If the input model is in mmCIF format, extra loops in the
file that are not modified will be kept.

Adding sequences will populate the entity_poly, entity_poly_seq,
struct_ref, and struct_ref_seq loops. The struct_ref and struct_ref_seq
loops are works in progress and will require external information. Also,
the canonical sequence containing all residues for the structure should
be provided. If there are non-standard residues in the model, they will
be automatically changed to the 3 letter residue name surrounded by
parentheses.

More functionality is planned.
'''

  datatypes = ['model', 'phil', 'restraint', 'sequence']

  master_phil_str = '''
mmtbx.validation.sequence.sequence_alignment {
  include scope mmtbx.validation.sequence.master_phil
}
custom_residues = None
  .type = strings
  .help = Space-separated list of three letter codes to be included in \
    the entity_poly.pdbx_one_letter_code mmCIF loop \
    (e.g. custom_residues='SUI PYL')
keep_original_loops = True
  .type = bool
  .help = Preserves mmCIF data from the input model file (if available) that \
    is not overwritten by other input
output {
  output_suffix = '.deposit.cif'
    .type = str
}
'''
  # ---------------------------------------------------------------------------
  def validate(self):
    self.data_manager.has_models(expected_n=1, exact_count=True, raise_sorry=True)
    self.data_manager.has_sequences(raise_sorry=True)

  # ---------------------------------------------------------------------------
  def run(self):
    print('Using model: {model_name}'.format(
      model_name=self.data_manager.get_default_model_name()))
    print('Using sequence(s): {sequence_names}'.format(
      sequence_names=', '.join(self.data_manager.get_sequence_names())))

    model = self.data_manager.get_model()
    model.set_log(self.logger)
    model.get_restraints_manager()

    # add sequences
    sequences = list()
    for sequence_name in self.data_manager.get_sequence_names():
      sequences.extend(self.data_manager.get_sequence(sequence_name))

    # match input sequences with model chains
    print(file=self.logger)
    print('Matching sequences to chains', file=self.logger)
    print('----------------------------', file=self.logger)
    # No error trapping
    model.set_sequences(
      sequences,
      custom_residues=self.params.custom_residues,
      similarity_matrix=self.params.mmtbx.validation.sequence.sequence_alignment.similarity_matrix,
      min_allowable_identity=self.params.mmtbx.validation.sequence.sequence_alignment.min_allowable_identity)
    model._sequence_validation.show(out=self.logger)
    print(file=self.logger)

    self.cif_model = model.model_as_mmcif(
      keep_original_loops=self.params.keep_original_loops)

    # write output file
    self.output_file = os.path.splitext(
      self.data_manager.get_default_model_name())[0] + \
      self.params.output.output_suffix
    print('Writing mmCIF', file=self.logger)
    print('-------------', file=self.logger)
    print ('  Output file = %s' % self.output_file, file=self.logger)
    self.data_manager.write_model_file(
      self.output_file, self.cif_model, self.params.output.overwrite)

  # ---------------------------------------------------------------------------
  def get_results(self):
    return group_args(output_file=self.output_file,
                      cif_model=self.cif_model)

# =============================================================================
# end
