from __future__ import absolute_import, division, print_function
import sys
import math
from datetime import datetime
from libtbx.program_template import ProgramTemplate
from libtbx import group_args, phil
from libtbx.str_utils import make_sub_header
from libtbx.utils import Sorry
import mmtbx
import mmtbx_probe_ext as probeExt
from mmtbx.probe import Helpers
from iotbx import pdb
from iotbx.pdb import common_residue_names_get_class

version = "0.1.0"

master_phil_str = '''
use_neutron_distances = False
  .type = bool
  .help = Use neutron distances (-nuclear in probe)

output
  .style = menu_item auto_align
{
  file_name = None
    .type = str
    .short_caption = Output file name
    .help = Output file name
}
''' + Helpers.probe_phil_parameters

program_citations = phil.parse('''
citation {
  authors = Word, et. al.
  journal = J. Mol. Biol.
  volume = 285
  pages = 1735-1747
  year = 1999
  external = True
}
''')

# ------------------------------------------------------------------------------

class Program(ProgramTemplate):
  description = '''
Reduce2 version {}
Add Hydrogens to a model and optimize their placement by adjusting movable groups and
flippable groups of atoms.

Inputs:
  PDB or mmCIF file containing atomic model
  Ligand CIF file, if needed
Output:
  PDB or mmCIF file with added hydrogens.

NOTES:
  Equivalent PHIL arguments for original Reduce command-line options:
     @todo
'''.format(version)
  datatypes = ['model', 'restraint', 'phil']
  master_phil_str = master_phil_str
  data_manager_options = ['model_skip_expand_with_mtrix']
  citations = program_citations
  epilog = '''
  For additional information and help, see http://kinemage.biochem.duke.edu/software/probe
  and http://molprobity.biochem.duke.edu
  '''

# ------------------------------------------------------------------------------

  def validate(self):
    self.data_manager.has_models(raise_sorry=True)
    if self.params.output.file_name is None:
      raise Sorry("Must specify output.file_name")

# ------------------------------------------------------------------------------

  def run(self):

    # @todo How to let it know to produce a PDB or CIF file as normal?  Is that by putting back in
    # get_results() below?

    # String describing the run that will be output to the specified file.
    outString = 'reduce2 v.{}, run {}\n'.format(version, datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    for a in sys.argv:
      outString += ' {}'.format(a)
    outString += '\n'

    make_sub_header('Interpret Model', out=self.logger)

    make_sub_header('Writing descriptive text', out=self.logger)

    # Write the output to the specified file.
    of = open(self.params.output.file_name,"w")
    of.write(outString)
    of.close()

# ------------------------------------------------------------------------------

  def Test(self):
    '''
      Run tests on the methods of the class.  Throw an assertion error if there is a problem with
      one of them and return normally if there is not a problem.
    '''

    #=====================================================================================
    # @todo Unit tests for other methods

# ------------------------------------------------------------------------------

  #def get_results(self):
  #  return group_args(model = self.model)

