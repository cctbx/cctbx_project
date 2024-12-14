"""
The top level View for for the molstar base app, which implements a Molstar viewer running 
  in a QtWebView with very rudimentary GUI controls. 
"""
from .base_app import BaseAppView
from qttbx.viewers.gui.view.molstar_base import MolstarTabViewBase

class MolstarBaseAppView(BaseAppView):
  def __init__(self,parent=None,params=None):
    super().__init__(parent=parent)

    # Viewer
    self.viewer_tab_view = MolstarTabViewBase(parent=self)
    self.viewer_tab_view.order_index=0
    self.tabs.insertTab(0,self.viewer_tab_view, "Viewer")
    self.tabs.widget(0).set_focus_on() # focus on first tab for right side

    