import sys
from pathlib import Path

from qtconsole.rich_jupyter_widget import RichJupyterWidget
from qtconsole.inprocess import QtInProcessKernelManager

from PySide2.QtCore import Slot
from PySide2.QtWidgets import (
    QLineEdit,
    QPlainTextEdit,
    QVBoxLayout
)
from .tab import GUITab



class JupyterTabWidget(GUITab):
  def __init__(self,parent=None):
    super().__init__(parent=parent)


    layout = QVBoxLayout()

    kernel_manager = QtInProcessKernelManager()


    kernel_manager.start_kernel()
    kernel_manager.kernel.shell.banner1 = f"""
    Execute python code here.
    Access the program with 'self'.

    For example:
      {'self.view':<30} : {'the Qt main window':<30}
      {'self.controller.molstar.graphics':<30} : {'the python interface for the Mol* viewer':<30}
      {'self.state.data_manager':<30} : {'the cctbx data manager':<30}
    """
    kernel_manager.kernel.gui = 'qt'

    kernel_client = kernel_manager.client()
    kernel_client.start_channels()

    self.jupyter_widget = RichJupyterWidget()
    self.jupyter_widget.kernel_manager = kernel_manager
    self.jupyter_widget.kernel_client = kernel_client



    layout.addWidget(self.jupyter_widget)
    self.setLayout(layout)



  def on_first_visit(self):
    pass