from __future__ import absolute_import, division, print_function
from phenix.program_template import ProgramTemplate
import subprocess
from libtbx import group_args
from mmtbx.monomer_library.pdb_interpretation import grand_master_phil_str
import mmtbx

from qttbx.viewers.gui.view.apps.molstar_base_app import MolstarBaseAppView
from qttbx.viewers.gui.controller.apps.molstar_base_app import MolstarBaseAppController
from qttbx.viewers.gui.modelstate import State
from qttbx.viewers.gui.apps import check_program_access

QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)
# =============================================================================

class Program(ProgramTemplate):

  description = """
  Demo program to visualize models in a QT gui connected to either Molstar or ChimeraX
  """

  datatypes = ['phil','model']


  master_phil_str = """

  rest_server_port = 5000
    .type = int
    .help = "The port for http control"

    include scope mmtbx.monomer_library.pdb_interpretation.grand_master_phil_str

  """

  def validate(self):
    pass


  def get_results(self):
    return group_args()


  def check_program_access(self,programs):
    inaccessible_programs = []

    for program in programs:
        try:
            subprocess.run([program, '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        except FileNotFoundError:
            inaccessible_programs.append(program)
        except subprocess.CalledProcessError:
            pass

    return inaccessible_programs

  def run(self):


    # first check that the necessary programs are available
    programs_to_check = ['npm', 'http-server']
    inaccessible_programs = self.check_program_access(programs_to_check)

    if inaccessible_programs:
      print(f"The following required programs are inaccessible or not found: {', '.join(inaccessible_programs)}")
      sys.exit()


    # start app
    QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)
    app = QApplication(sys.argv)

    # get icon
    #icon_path =  Path(__file__).parent / '../view/assets/icons/phenix/icon.icns'
    #icon = QIcon(str(icon_path))
    #qapp.setWindowIcon(icon)
    
    # Core top level object initialization
    self.state = State(self.dm,params=self.params)
    self.view = MolstarBaseAppView()
    self.controller = MolstarBaseAppController(parent=state,view=view)

    # Start
    self.controller.view.show()
    sys.exit(app.exec_())
