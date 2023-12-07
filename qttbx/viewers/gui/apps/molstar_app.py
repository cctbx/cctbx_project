
import time
import sys
from pathlib import Path

from PySide2.QtGui import QIcon
from PySide2.QtWidgets import QApplication, QWidget, QApplication, QWidget
from PySide2.QtCore import QObject, QEvent, Qt,  QEvent
from PySide2.QtSvg import QSvgRenderer

from iotbx.data_manager import DataManager
from ..view.apps.molstar_app import MolstarAppView
from ..controller.apps.molstar_app import MolstarAppController
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


class MolstarApp:
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
  icon_path =  Path(__file__).parent / '../view/assets/icons/phenix/icon.icns'
  #print("DEBUG: icon path for phenix: ",icon_path)
  icon = QIcon(str(icon_path))
  qapp.setWindowIcon(icon)


  dm = DataManager()

  # Core top level object initialization
  state = State(dm)  
  view = MolstarAppView()
  controller = MolstarAppController(parent=state,view=view)
  app = MolstarApp(state,view,controller)

  # DEBUG: Sync references for test data
  state.signals.references_change.emit()
  

  # Reach into the Console tab to make variables accessible
  app.view.python_console.jupyter_widget.kernel_manager.kernel.shell.push({'app': app})
  #include Selection dataclasses to build querys in console
  app.view.python_console.jupyter_widget.kernel_manager.kernel.shell.push({'Selection':Selection})
  app.view.python_console.jupyter_widget.kernel_manager.kernel.shell.push({'SelectionQuery':SelectionQuery})

  controller.view.show()

  sys.exit(qapp.exec_())