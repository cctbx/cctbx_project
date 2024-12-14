"""
The top level View for for the molstar base app, which implements a Molstar viewer running 
  in a QtWebView with very rudimentary GUI controls. 
"""

from PySide2.QtCore import Signal
from PySide2.QtGui import QGuiApplication
from PySide2.QtWidgets import QMainWindow

from qttbx.viewers.gui.view.molstar_selection import MolstarSelectionTabView
from qttbx.viewers.gui.view.widgets.console import JupyterTabWidget
from qttbx.viewers.gui.view.selection import SelectionsTabView, MolstarSelectionTabView
from qttbx.viewers.gui.view.apps.base_app import BaseAppView

class SelectionEditorAppView(BaseAppView):

  def __init__(self,parent=None,params=None):
    super().__init__(parent=parent)

    self.viewer_tab_view = MolstarSelectionTabView(parent=self)
    self.viewer_tab_view.order_index=0
    self.tabs.insertTab(0,self.viewer_tab_view, "Viewer")

    self.selection_tab_view = SelectionsTabView(parent=self)
    self.selection_tab_view.order_index=1
    self.tabs.insertTab(1,self.selection_tab_view, "Selections")

    # Pop out viewer by default
    self.tabs.simulate_drag_out(0)
    self.tabs.widget(0).set_focus_on() # focus on first tab for right side