
import time
import sys
from pathlib import Path
#import argparse

from PySide2.QtGui import QIcon
from PySide2.QtWidgets import QApplication, QWidget, QApplication, QWidget, QPushButton, QVBoxLayout
from PySide2.QtWidgets import QApplication, QDialog, QVBoxLayout, QPushButton, QMessageBox, QMainWindow

from PySide2.QtCore import QObject, QEvent, Qt,  QEvent, QSize
from PySide2.QtSvg import QSvgRenderer

from iotbx.data_manager import DataManager
from ..view.apps.demo import DemoView
from ..controller.apps.demo import DemoController
from ..state.state import State
from ...last.selection_utils import Selection, SelectionQuery
from . import ViewerChoiceDialog, check_program_access

QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)



class DemoApp:
  def __init__(self,state,view,controller):
    self.controller = controller
    self.view = view
    self.state = state
  



def main(viewer=None,dm=None,log=None):

  # first check that the necessary programs are available
  programs_to_check = ['npm', 'http-server']
  inaccessible_programs = check_program_access(programs_to_check)

  if inaccessible_programs:
      print(f"The following required programs are inaccessible or not found: {', '.join(inaccessible_programs)}")
      sys.exit()
  else:
      print("All programs are accessible.")


  choice = viewer # 'molstar' or 'chimerax'
  

  # start app
  QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)
  qapp = QApplication(sys.argv)

  # get icon
  icon_path =  Path(__file__).parent / '../view/assets/icons/phenix/icon.icns'
  icon = QIcon(str(icon_path))
  qapp.setWindowIcon(icon)

  # choose viewer
  if not choice:
    choice_dialog = ViewerChoiceDialog()
    res = choice_dialog.exec_()

    # If a choice was made, show the main window
    if res != QDialog.Accepted:
      QMessageBox.warning(None, "No Choice", "No option was selected. Exiting application.")
      qapp.quit()
      sys.exit()
    
    
    
    choice = choice_dialog.choice
  
  # Set up a data manager
  if not dm:
    dm = DataManager()

  
  # DEBUG: load some data automatically
  #dm.process_model_file("/Users/user/software/phenix/modules/cctbx_project/qttbx/data/1yjp.pdb")
  #dm.process_real_map_file("/Users/user/software/phenix/modules/cctbx_project/qttbx/data/1yjp_calc.mrc")
  #dm.process_model_file("/Users/user/software/phenix/modules/cctbx_project/qttbx/data/1yjp.cif")

  # Core top level object initialization
  state = State(dm)  
  view = DemoView(viewer_choice=choice)
  controller = DemoController(parent=state,view=view,viewer_choice=choice,log=log)
  app = DemoApp(state,view,controller)

  # DEBUG: Sync references for test data
  state.signals.references_change.emit()
  

  # # Reach into the Console tab to make variables accessible
  # try:
  #   import qtconsole
  #   app.view.python_console.jupyter_widget.kernel_manager.kernel.shell.push({'app': app})
  #   #include Selection dataclasses to build querys in console
  #   app.view.python_console.jupyter_widget.kernel_manager.kernel.shell.push({'Selection':Selection})
  #   app.view.python_console.jupyter_widget.kernel_manager.kernel.shell.push({'SelectionQuery':SelectionQuery})

  # except:
  #   print("No qtconsole found")


  # # Create an instance of the event filter
  # globalEventFilter = GlobalEventFilter()

  # # Install the event filter on the QApplication instance
  # qapp.installEventFilter(globalEventFilter)

  controller.view.show()

  sys.exit(qapp.exec_())

if __name__ == '__main__':
  main()