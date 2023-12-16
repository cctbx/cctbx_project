from __future__ import absolute_import, division, print_function
from phenix.program_template import ProgramTemplate
from libtbx import group_args
from qttbx.viewers.gui.apps.main import main


# =============================================================================

class Program(ProgramTemplate):

  description = """
  Demo program to visualize models in a QT gui connected to either Molstar or ChimeraX
  """

  datatypes = ['phil','model','real_map']


  master_phil_str = """

  viewer_choice = None
    .type = str
    .help = 'molstar' or 'chimerax'
    .short_caption = The viewer to use. None for gui prompt

  show_tab = None
    .type = str
    .multiple = True
    .help = 'all' or None or 'atoms','cif', 'restraints', etc
    .short_caption = Names of the tabs to include in the viewer. \
    'viewer','selections' and 'files' are required and default
  """

  def validate(self):
    pass

  def run(self):
    viewer_name = self.params.viewer_choice
    main(dm=self.data_manager,params=self.params,log=self.logger)

  def get_results(self):
    return group_args()