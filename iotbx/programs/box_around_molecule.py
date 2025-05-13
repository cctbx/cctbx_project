from __future__ import absolute_import, division, print_function
import os
from libtbx.program_template import ProgramTemplate
from cctbx.maptbx.box import shift_and_box_model

master_phil_str = '''
buffer_layer = 5.0
  .type = float
  .help = buffer around atoms, in Angstrom
shift_model = True
  .type = bool
  .help = shift model closer to the origin
output {
  suffix = _box
    .type = str
  format = *Auto pdb cif
    .type = choice(multi=False)
    .help = output model file type
}
'''

# ------------------------------------------------------------------------------

class Program(ProgramTemplate):
  description = '''
Script that creates a P1 box around a model.

Inputs:
  PDB or mmCIF file containing atomic model
  optional: buffer layer in Angstrom (default is 5.0 A)

Usage examples:
  1) Save model as PDB (input is mmCIF)
     iotbx.pdb.box_around_molecule.py test.cif format=pdb

  2) Save model as mmCIF (input is PDB)
    pdb.box_around_molecule.py test.pdb format=cif

  3) Default: input and output format are the same
     pdb.box_around_molecule.py test.pdb

  The program template accepts PDB and mmCIF input formats.
  The paramter output.format allows choosing the output format.
'''

  datatypes = ['model', 'phil']
  master_phil_str = master_phil_str

  # ----------------------------------------------------------------------------

  def validate(self):
    self.data_manager.has_models(expected_n  = 1,
                                 exact_count = True,
                                 raise_sorry = True)

  # ----------------------------------------------------------------------------

  def run(self):
    #
    self.model = self.data_manager.get_model()

    # shift_and_box_model creates a new model object, so the format needs to
    # be obtained here
    if self.params.output.format == 'Auto':
      if self.model.input_model_format_pdb():
        self.params.output.format = 'pdb'
      elif self.model.input_model_format_cif():
        self.params.output.format = 'cif'

    self.model = shift_and_box_model(model = self.model,
                                     box_cushion = self.params.buffer_layer,
                                     shift_model=self.params.shift_model)

    # Get output filename if not provided
    fn = self.data_manager.get_default_model_name()
    basename = os.path.splitext(os.path.basename(fn))[0]
    if self.params.output.prefix is None:
      self.params.output.prefix = basename
      self.data_manager.set_default_output_filename(self.get_default_output_filename())

    # save output
    output_fn = self.data_manager.write_model_file(self.model,
                                       format = self.params.output.format)
    print('Created file', output_fn, file=self.logger)
