
import time
import sys
from pathlib import Path

from PySide2.QtGui import QIcon
from PySide2.QtWidgets import QApplication, QWidget, QApplication, QWidget
from PySide2.QtCore import QObject, QEvent, Qt,  QEvent
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
  


if __name__ == '__main__':

  # start app
  QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)
  qapp = QApplication(sys.argv)
  #qapp.setStyleSheet("QWidget { border: 1px solid black; }")

  # get icon
  icon_path =  Path(__file__).parent / '../views/assets/icons/phenix/icon.icns'
  icon = QIcon(str(icon_path))
  qapp.setWindowIcon(icon)


  dm = DataManager()
  # DEBUG: load some data automatically
  #dm.process_model_file("/Users/user/software/phenix/modules/cctbx_project/qttbx/data/1yjp.pdb")
  dm.process_real_map_file("/Users/user/software/phenix/modules/cctbx_project/qttbx/data/1yjp_calc.mrc")

  dm.process_model_file("/Users/user/software/phenix/modules/cctbx_project/qttbx/data/1yjp.cif")

  # Core top level object initialization
  state = State(dm)  
  view = DemoView()
  controller = DemoController(parent=state,view=view)
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