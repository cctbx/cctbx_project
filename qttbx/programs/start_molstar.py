from __future__ import absolute_import, division, print_function
from phenix.program_template import ProgramTemplate
from libtbx import group_args
from mmtbx.monomer_library.pdb_interpretation import grand_master_phil_str
from qttbx.viewers.gui.apps.molstar_base_app import main
import mmtbx

# =============================================================================

class Program(ProgramTemplate):

  description = """
  Demo program to visualize models in a QT gui connected to either Molstar or ChimeraX
  """

  datatypes = ['phil','model','real_map']


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

    include scope mmtbx.monomer_library.pdb_interpretation.grand_master_phil_str

  """

  def validate(self):
    pass

  def run(self):
    viewer_name = self.params.viewer_choice
    main(dm=self.data_manager,params=self.params,log=self.logger)

  def get_results(self):
    return group_args()