
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
QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)

class GlobalEventFilter(QObject):
  def eventFilter(self, watched, event):
    if event.type() == QEvent.Enter:
      if isinstance(watched, QWidget):
        print(f"Mouse entered: {watched.objectName()} ({type(watched).__name__})")
    elif event.type() == QEvent.Leave:
      if isinstance(watched, QWidget):
        print(f"Mouse left: {watched.objectName()} ({type(watched).__name__})")
    return super().eventFilter(watched, event)


class DemoApp:
  def __init__(self,state,view,controller):
    self.controller = controller
    self.view = view
    self.state = state
  

class ViewerChoiceDialog(QDialog):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Choose Viewer")
        self.choice = None
        self.initUI()

    def initUI(self):
        self.setMinimumSize(QSize(300, 100))  # Example size: 300x


        layout = QVBoxLayout()
        

        # Buttons for choices
        btnOption1 = QPushButton("Molstar", self)
        btnOption1.clicked.connect(self.select_molstar)

        btnOption2 = QPushButton("ChimeraX", self)
        btnOption2.clicked.connect(self.select_chimerax)

        layout.addWidget(btnOption1)
        layout.addWidget(btnOption2)

        self.setLayout(layout)

    def select_molstar(self):
        self.choice = 'molstar'
        self.accept()

    def select_chimerax(self):
        self.choice = 'chimerax'
        self.accept()

def main(viewer=None,dm=None,log=None):
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
  

  # Reach into the Console tab to make variables accessible
  app.view.python_console.jupyter_widget.kernel_manager.kernel.shell.push({'app': app})
  #include Selection dataclasses to build querys in console
  app.view.python_console.jupyter_widget.kernel_manager.kernel.shell.push({'Selection':Selection})
  app.view.python_console.jupyter_widget.kernel_manager.kernel.shell.push({'SelectionQuery':SelectionQuery})


  # # Create an instance of the event filter
  # globalEventFilter = GlobalEventFilter()

  # # Install the event filter on the QApplication instance
  # qapp.installEventFilter(globalEventFilter)

  controller.view.show()

  sys.exit(qapp.exec_())

# if __name__ == '__main__':
#   # arguments
#   parser = argparse.ArgumentParser(description='Phenix Viewer Demo')
#   parser.add_argument('--viewer', type=str, help="Either 'molstar' or 'chimerax'", required=False)
#   args = parser.parse_args()
#   choice = args.viewer # viewer choice
#   main(viewer=choice)