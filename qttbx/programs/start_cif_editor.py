from __future__ import absolute_import, division, print_function
from phenix.program_template import ProgramTemplate
from libtbx import group_args
from qttbx.viewers.gui.apps.cif_editor import main


# =============================================================================

class Program(ProgramTemplate):

  description = """
  Demo program to visualize CIF formatted data and make edits
  """

  datatypes = ['phil','model','restraint']


  master_phil_str = """
  viewer_choice = 'molstar'
    .type = str
    .help = 'molstar' or 'chimerax'
    .short_caption = The viewer to use. None for gui prompt

  show_tab = None
    .type = str
    .multiple = True
    .help = 'all' or None or 'atoms','cif', 'restraints', etc
    .short_caption = Names of the tabs to include in the viewer. \
    'viewer','selections' and 'files' are required and default

  rest_server_port = 5000
    .type = int
    .help = "The port for http control"
  """

  def validate(self):
    pass

  def run(self):
    main(dm=self.data_manager,params=self.params,log=self.logger)

  def get_results(self):
    return group_args()
