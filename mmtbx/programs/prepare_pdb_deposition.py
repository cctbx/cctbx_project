from __future__ import division, print_function

import iotbx.pdb
import iotbx.phil
import mmtbx.model
from libtbx import group_args
from libtbx.program_template import ProgramTemplate
from libtbx.utils import Sorry

'''
Code will be changed once new command-line parser is functioning
'''

# =============================================================================
description = '''
Program for preparing model and data files for depostion into the Proten Data
Bank

Minimum required data:
  Model file
  Sequence file

The sequence file should have a sequence for each chain in the model file.

Currently, this program only combines the model and sequence into a single mmCIF
file. More functionality is planned.
'''

# =============================================================================
master_params_str = '''
input {
  model_file = None
    .type = path
  sequence_file = None
    .type = path
}
sequence_alignment {
  include scope mmtbx.validation.sequence.master_phil
}
'''

master_params = iotbx.phil.parse(master_params_str, process_includes=True)

# =============================================================================
class Program(ProgramTemplate):

  def validate(self):
    print('Validating inputs', file=self.logger)
    if (len(self.data_manager.get_model_names()) == 0):
      raise Sorry('One model file is required.')
    if (len(self.data_manager.get_sequence_names()) == 0):
      raise Sorry('One sequence file is required.')

    print ('Input files:', file=self.logger)
    print ('  Model file = %s' % self.data_manager.get_default_model_name(),
           file=self.logger)
    print ('  Sequence file = %s' %
           self.data_manager.get_default_sequence_name(), file=self.logger)

  def run(self):
    model = self.data_manager.get_model()
    self.cif_blocks = list()

    # sequence block
    print ('Creating mmCIF block for sequence', file=self.logger)
    sequence = self.data_manager.get_sequence()
    hierarchy = model._pdb_hierarchy
    seq_block = hierarchy.as_cif_block_with_sequence(
      sequence, crystal_symmetry=model.crystal_symmetry(),
      alignment_params=self.params.sequence_alignment)
    self.cif_blocks.append(seq_block)

    # create additional cif blocks?

    # add cif blocks together
    print ('Creating complete mmCIF', file=self.logger)
    self.cif_model = model.model_as_mmcif(additional_blocks=self.cif_blocks)

    # write output file
    self.output_file = self.data_manager.get_default_model_name().split('.')[0]+\
                       '.deposit.cif'
    print ('Writing mmCIF', file=self.logger)
    print ('  Output file = %s' % self.output_file, file=self.logger)
    with open(self.output_file, 'wb') as f:
      f.write(self.cif_model)

    # update data manager for any downstream applications
    pdb_input = iotbx.pdb.input(self.output_file)
    model = mmtbx.model.manager(model_input=pdb_input, log=self.logger)
    self.data_manager.add_model(self.output_file, model)
    self.data_manager.set_default_model(self.output_file)

  def get_results(self):
    return group_args(output_file=self.output_file,
                      cif_model=self.cif_model)

# =============================================================================
