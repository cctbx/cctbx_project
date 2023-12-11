from __future__ import absolute_import, division, print_function
from phenix.program_template import ProgramTemplate
from libtbx import group_args
from qttbx.viewers.gui.apps.demo import main


# =============================================================================

class Program(ProgramTemplate):

  description = """
  Demo program to visualize models in a QT gui connected to either Molstar or ChimeraX
  """

  datatypes = ['phil','model', 'real_map']


  master_phil_str = """

  viewer {
    choice = None
      .type = str
      .help = 'molstar' or 'chimerax'
      .short_caption = The viewer to use
    }

  """

  def validate(self):
    pass

  def run(self):
    viewer_name = self.params.viewer.choice
    main(viewer=viewer_name,dm=self.data_manager,log=self.logger)

  def get_results(self):
    return group_args()