from PySide2.QtWidgets import QApplication
from PySide2.QtCore import QEvent

from qttbx.viewers.gui.controller.selection import SelectionTabController
from qttbx.viewers.gui.controller.molstar_selection import MolstarSelectionController
from qttbx.viewers.gui.controller.apps.molstar_base_app import BaseAppController


class SelectionEditorAppController(BaseAppController):
  """
  This is the top level Controller instance for the selection App
  """
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)

    self.molstar = MolstarSelectionController(parent=self,view=self.view.viewer_tab_view)
    self.selection = SelectionTabController(parent=self,view=self.view.selection_tab_view)