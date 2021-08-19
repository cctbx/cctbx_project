from __future__ import absolute_import, division, print_function
import os
from libtbx.program_template import ProgramTemplate
#from libtbx.utils import null_out
from libtbx import group_args, phil
from libtbx.str_utils import make_sub_header
from libtbx.utils import Sorry
from mmtbx.hydrogens import reduce_hydrogen

master_phil_str = '''
use_neutron_distances = False
  .type = bool
  .help = Use neutron distances

output
  .style = menu_item auto_align
{
  file_name_prefix = None
    .type = path
    .short_caption = Prefix for file name
    .help = Prefix for file name
    .input_size = 400
}
'''

program_citations = phil.parse('''
citation {
  authors = Word, et. al.
  journal = J. Mol. Biol.
  volume = 285
  pages = 1711-1733
  year = 1999
  external = True
}
''')

# ------------------------------------------------------------------------------

class Program(ProgramTemplate):
  description = '''
Compute the MolProbity Probe score for a file, or a subset of the file.

Inputs:
  PDB or mmCIF file containing atomic model
  Ligand CIF file, if needed
'''
  datatypes = ['model', 'restraint', 'phil']
  master_phil_str = master_phil_str
  citations = program_citations
  epilog = '''
  For additional information and help, see http://molprobity.biochem.duke.edu
  '''

# ------------------------------------------------------------------------------

  def validate(self):
    self.data_manager.has_models(raise_sorry=True)
    if self.params.output.file_name_prefix is None:
      raise Sorry("Supply the prefix for an output file name using output.file_name_prefix=")

# ------------------------------------------------------------------------------

  def run(self):
    self.model = self.data_manager.get_model()
    #
    make_sub_header('Compute Probe Score', out=self.logger)
    #reduce_add_h_obj = reduce_hydrogen.place_hydrogens(
    #  model = self.model,
    #  use_neutron_distances = self.params.use_neutron_distances)
    #reduce_add_h_obj.run()
    #self.model = reduce_add_h_obj.get_model()
    #reduce_add_h_obj.show(log = self.logger)
    #
    base = self.params.output.file_name_prefix
    of = open("%s.txt"%base,"w")
    of.write("@todo")
    of.close()

# ------------------------------------------------------------------------------

  #def get_results(self):
  #  return group_args(model = self.model)
