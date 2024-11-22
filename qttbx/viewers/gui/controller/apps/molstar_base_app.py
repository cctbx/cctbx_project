"""
The top level Controller for for the molstar base app, which implements a Molstar viewer running 
in a QtWebView with very rudimentary GUI controls. 
"""
from PySide2.QtWidgets import QApplication
from PySide2.QtCore import QEvent

from qttbx.viewers.gui.controller.molstar_base import MolstarBaseController
from qttbx.viewers.gui.controller.apps.base_app import BaseAppController


class MolstarBaseAppController(BaseAppController):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)
    self.molstar = MolstarBaseController(parent=self,view=self.view.viewer_tab_view)
