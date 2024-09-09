"""
The top level Controller for for the molstar base app, which implements a Molstar viewer running 
in a QtWebView with very rudimentary GUI controls. 
"""
from PySide2.QtWidgets import QApplication
from PySide2.QtCore import QEvent

from qttbx.viewers.gui.controller.molstar import MolstarController
from qttbx.viewers.gui.controller import Controller


class MolstarBaseAppController(Controller):
  """
  This is the top level Controller instance for the base molstar app
  """
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)


    # Main Level Components
    self.molstar = MolstarController(parent=self,view=self.view.viewer_tab_view)

    # signals
    self.view.signal_close.connect(self.close_event)
    self.state.signals.tab_change.connect(self.change_tab_to)

    # Finish
    self.state.start()


  def close_application(self):
    # manually call this function to close
    closeEvent = QEvent(QEvent.Close)
    QApplication.sendEvent(self.view, closeEvent)

  def close_event(self):
    # slot for built in close event
    if hasattr(self,"molstar"):
      self.molstar.close_viewer()
    if hasattr(self,"chimerax"):
      self.chimerax.close_viewer()

  def change_tab_to(self,name):
    """
    Change the active tab by a name
    """
    # TODO: Clean up tab change behavior for Selection and Viewer as child
    #if self.view._has_child_window:
    #self.log("Change to tab: ",name)
    tab_widget = self.view.tabs
    for index in range(tab_widget.count()):
      if tab_widget.tabText(index) == name:
        tab_widget.setCurrentIndex(index)
        return
