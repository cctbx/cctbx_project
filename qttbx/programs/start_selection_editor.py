from __future__ import absolute_import, division, print_function
import subprocess
import sys
from pathlib import Path

from qttbx.programs.start_molstar_base import MolstarBaseProgram
from libtbx import group_args
import mmtbx
from mmtbx.monomer_library.pdb_interpretation import grand_master_phil_str
from qttbx.viewers.gui.view.apps.selection_editor import SelectionEditorAppView
from qttbx.viewers.gui.controller.apps.selection_editor import SelectionEditorAppController
from qttbx.viewers.gui.model.state import State

from PySide2.QtCore import Qt
from PySide2.QtGui import QIcon
from PySide2.QtWidgets import QApplication
QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)
# =============================================================================

master_phil_str = """

    include scope mmtbx.monomer_library.pdb_interpretation.grand_master_phil_str
    include scope qttbx.programs.start_molstar_base.master_phil_str

  """

class SelectionProgram(MolstarBaseProgram):

  description = """
  Program to view selections
  """

  datatypes = ['phil','model']

  master_phil_str = master_phil_str

  def run(self):

    # start app
    QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)
    app = QApplication(sys.argv)

    # get icon
    icon_path =  Path(__file__).parent / '../viewers/gui/view/assets/icons/phenix.icns'
    icon = QIcon(str(icon_path))
    app.setWindowIcon(icon)
    
    # Core top level object initialization
    self.state = State(self.data_manager,params=self.params)
    self.view = SelectionEditorAppView(params=self.params)
    self.controller = SelectionEditorAppController(parent=self.state,view=self.view)

    # Reach into the Console tab to make variables accessible
    self.view.python_console.jupyter_widget.kernel_manager.kernel.shell.push({'self': self})

    # Start
    self.controller.view.show()
     
    sys.exit(app.exec_())