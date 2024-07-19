from pathlib import Path

from PySide2.QtGui import QIcon
from PySide2.QtWidgets import (
    QComboBox,
    QFrame,
    QHBoxLayout,
    QPushButton,
    QVBoxLayout
)

from ..table import PandasTableView
from ..cif import CifBrowserTabView
from ..widgets.tab import GUITab



class RestraintTableView(PandasTableView):

  def __init__(self, parent=None):
    super().__init__(parent=parent)

  


class RestraintBrowserTabView(CifBrowserTabView):
  """
  View restraints cif file structure
  """
  def __init__(self,parent=None):
    super().__init__(parent=parent,load_button=False,next_button=True) # Don't let users manually load restraint files

    self.combobox_data_block.hide() # Hide the data_block dropdown for simplicity
    #self.data_block_layout.hide()
    self.combo_button_layout.removeItem(self.data_block_layout)

  
