"""Program for preparing model and data files for depostion into the Protein
Data Bank
"""
from __future__ import absolute_import, division, print_function

import os

from six.moves import cStringIO as StringIO

from libtbx import group_args
from libtbx.program_template import ProgramTemplate
from libtbx.utils import Sorry

# =============================================================================
class Program(ProgramTemplate):

  description = '''
Program for preparing model and data files for depostion into the Protein
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
align_columns = False
  .type = bool
  .help = When set to True, the columns are aligned. This will take longer \
    because the column widths need to be deteremined before outputting.
include scope mmtbx.monomer_library.pdb_interpretation.grand_master_phil_str
include scope mmtbx.geometry_restraints.torsion_restraints.reference_model.reference_model_str
include scope mmtbx.geometry_restraints.external.external_energy_params_str
output {
  suffix = '.deposit'
    .type = str
    .help = Suffix string added to automatically generated output filenames
}

# PHIL parameters for current GUI
include scope libtbx.phil.interface.tracking_params
gui
  .help = "GUI-specific parameter required for output directory"
{
  output_dir = None
  .type = path
  .style = output_dir
}
'''
  # ---------------------------------------------------------------------------
  def custom_init(self):
    self.output_file = None

  # ---------------------------------------------------------------------------
  def validate(self):
    self.data_manager.has_models(expected_n=1, exact_count=True, raise_sorry=True)
    self.data_manager.has_sequences(raise_sorry=True)
    model = self.data_manager.get_model()
    if not model.input_model_format_cif():
      raise Sorry('Your input model must already be in mmCIF.')

  # ---------------------------------------------------------------------------
  def run(self):
    print('Using model: {model_name}'.format(
      model_name=self.data_manager.get_default_model_name()))
    print('Using sequence(s): {sequence_names}'.format(
      sequence_names=', '.join(self.data_manager.get_sequence_names())))

    # get model
    model = self.data_manager.get_model()
    model.set_log(self.logger)

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

    # When the input is already in mmCIF, just add sequence information
    found_cif_block = False
    if hasattr(model._model_input, 'cif_model'):
      for cif_block in model._model_input.cif_model.values():
        if '_atom_site' in cif_block:
          cif_block.update(model._sequence_validation.sequence_as_cif_block())
          out = StringIO()
          model._model_input.cif_model.show(
            out=out, align_columns=self.params.align_columns)
          self.cif_model = out.getvalue()
          found_cif_block = True
          break
    # otherwise, generate mmCIF
    if not found_cif_block:
      model.get_restraints_manager()
      self.cif_model = model.model_as_mmcif(
        align_columns=self.params.align_columns,
        keep_original_loops=self.params.keep_original_loops)

    # show sequence alignment
    model._sequence_validation.show(out=self.logger)
    print(file=self.logger)

    # write output file to current directory
    if self.params.output.prefix is None:
      self.params.output.prefix = os.path.split(os.path.splitext(
        self.data_manager.get_default_model_name())[0])[1]
    self.data_manager.set_default_output_filename(
      self.get_default_output_filename())
    self.output_file = self.data_manager.get_default_output_model_filename()
    print('Writing mmCIF', file=self.logger)
    print('-------------', file=self.logger)
    print ('  Output file = %s' % self.output_file, file=self.logger)
    self.data_manager.write_model_file(self.cif_model)

  # ---------------------------------------------------------------------------
  def get_results(self):
    return group_args(output_file=self.output_file,
                      cif_model=self.cif_model)

# =============================================================================
# end
