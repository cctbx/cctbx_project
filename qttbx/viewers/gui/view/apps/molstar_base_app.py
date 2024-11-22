"""
The top level View for for the molstar base app, which implements a Molstar viewer running 
  in a QtWebView with very rudimentary GUI controls. 
"""
from .base_app import BaseAppView
from qttbx.viewers.gui.view.molstar_base import MolstarTabView

class MolstarBaseAppView(BaseAppView):
  def __init__(self,parent=None,params=None):
    super().__init__(parent=parent)

    # Viewer
    self.viewer_tab_view = MolstarTabView(parent=self)
    self.viewer_tab_view.order_index=0
    self.tabs.insertTab(0,self.viewer_tab_view, "Viewer")
    self.tabs.findTabByName("Viewer").set_focus_on()

    