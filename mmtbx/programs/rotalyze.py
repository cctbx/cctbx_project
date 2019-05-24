from __future__ import absolute_import, division, print_function

import os
from mmtbx.validation.rotalyze import rotalyze
from libtbx.program_template import ProgramTemplate
from libtbx.utils import Sorry

class Program(ProgramTemplate):
  prog = os.getenv('LIBTBX_DISPATCHER_NAME')
  description="""\
%(prog)s file.pdb [params.eff] [options ...]

Options:

  model=input_file        input PDB file
  outliers_only=False   only print outliers
  verbose=False         verbose text output

Example:

  %(prog)s model=1ubq.pdb outliers_only=True
""" % locals()

  master_phil_str = """
  include scope mmtbx.validation.molprobity_cmdline_phil_str
  show_errors = False
    .type = bool
    .help = '''Print out errors'''
  wxplot = False
    .type = bool
    .help = Display interactive plots (requires wxPython and Matplotlib)
  """
  datatypes = ['model','phil']
  known_article_ids = ['molprobity']

  def validate(self):
    self.data_manager.has_models(raise_sorry=True)

  def run(self):
    hierarchy = self.data_manager.get_model().get_hierarchy()

    result = rotalyze(
      pdb_hierarchy=hierarchy,
      data_version="8000",#params.data_version,
      show_errors=self.params.show_errors,
      outliers_only=self.params.outliers_only,
      out=self.logger,
      quiet=False)
    if self.params.verbose:
      result.show_old_output(out=self.logger, verbose=True)
    if self.params.wxplot :
      try :
        import wxtbx.app
      except ImportError as e :
        raise Sorry("wxPython not available.")
      else :
        app = wxtbx.app.CCTBXApp(0)
        result.display_wx_plots()
        app.MainLoop()
