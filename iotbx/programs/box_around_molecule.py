from __future__ import absolute_import, division, print_function
from libtbx.program_template import ProgramTemplate
#from libtbx.utils import date_and_time
import libtbx.load_env
from cctbx.maptbx.box import shift_and_box_model

master_phil_str = '''
buffer_layer = 5.0
  .type = float
  .help = buffer around atoms, in Angstrom
output
{
  format = *pdb cif
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
  1) Print PDB format string of model (input is mmCIF)
     iotbx.pdb.box_around_molecule.py test.cif format=pdb

  2) Print mmCIF format string of model (input is PDB format)
  pdb.box_around_molecule.py test.pdb format=cif

  3) Print and save model in PDB format (input is PDB format)
    pdb.box_around_molecule.py test.pdb format=pdb prefix=toto

  The program template accepts PDB and mmCIF input formats.
  The paramter output.format allows choosing the output format.
  If output filename or prefix is given, a file is saved; otherwise, the result
  is only printed.
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
    self.model = shift_and_box_model(model = self.model,
                                     box_cushion = self.params.buffer_layer
                                     )
    # PDB format
    if (self.params.output.format == 'pdb'):
      # If output file should contain REMARK with buffer_layer
      #model_str = 'REMARK %s --buffer-layer=%.6g %s\n' % (
      #            libtbx.env.dispatcher_name, self.params.buffer_layer, fn) + \
      #       'REMARK %s\n' % date_and_time() + \
      #       self.model.model_as_pdb()
      print(self.model.model_as_pdb(), file=self.logger)
    # mmCIF format
    if (self.params.output.format == 'cif'):
      print(self.model.model_as_mmcif(), file=self.logger)

    self.data_manager.write_model_file(self.model,
                                       format = self.params.output.format)
