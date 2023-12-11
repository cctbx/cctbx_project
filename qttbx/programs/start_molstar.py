from __future__ import absolute_import, division, print_function
from phenix.program_template import ProgramTemplate
from libtbx import group_args
from qttbx.viewers.gui.apps.molstar_app import main


# =============================================================================

class Program(ProgramTemplate):

  description = """
  Star the Molstar viewer in a QT webview window
  """

  datatypes = ['phil','model', 'real_map']


  def validate(self):
    pass

  def run(self):
    viewer_name = self.params.viewer.choice
    main(viewer=viewer_name,dm=self.data_manager,log=self.logger)

  def get_results(self):
    return group_args()