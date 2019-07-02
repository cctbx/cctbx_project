from __future__ import absolute_import, division, print_function

import os
from mmtbx.validation.cbetadev import cbetadev
from libtbx.program_template import ProgramTemplate
#from libtbx.utils import Sorry

class Program(ProgramTemplate):
  prog = os.getenv('LIBTBX_DISPATCHER_NAME')
  description="""
  %(prog)s file.pdb [params.eff] [options ...]

Calculates ideal CB atom position based on mainchain geometry,
then displays deviation of modeled CB from that ideal positon.
Deviations of >= 0.25A are considered outliers

Options:

  model=input_file      input PDB file
  output=text or kin    select type of output
  outliers_only=False   suppress non-outlier results

Example:

  %(prog)s model=1ubq.pdb
""" % locals()

  master_phil_str = """
  include scope mmtbx.validation.molprobity_cmdline_phil_str
    cbetadev {
      output = *text kin
      .type = choice
      .help = '''choose output type'''
      }
"""
  datatypes = ['model','phil']
  known_article_ids = ['molprobity']

  def validate(self):
    self.data_manager.has_models(raise_sorry=True)

  def run(self):
    hierarchy = self.data_manager.get_model().get_hierarchy()
    hierarchy.atoms().reset_i_seq()
    result = cbetadev(
      pdb_hierarchy=hierarchy,
      outliers_only=self.params.outliers_only,
      out=self.logger,
      quiet=False)
    if self.params.cbetadev.output == "kin":
      self.logger.write(result.as_kinemage())
    elif self.params.verbose:
      #pdb_file_str = os.path.basename(self.params.model)[:-4]
      #get input file name from data manager, strip file extension
      pdb_file_str = os.path.basename(self.data_manager.get_model_names()[0])[:-4]
      result.show_old_output(out=self.logger, prefix=pdb_file_str, verbose=True)
