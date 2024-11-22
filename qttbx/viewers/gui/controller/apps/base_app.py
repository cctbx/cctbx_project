from PySide2.QtWidgets import QApplication
from PySide2.QtCore import QEvent

from qttbx.viewers.gui.controller.molstar_base import MolstarBaseController
from qttbx.viewers.gui.controller import Controller


class BaseAppController(Controller):
  """
  This is the top level Controller instance for the base molstar app
  """
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)

    # signals
    self.view.signal_close.connect(self.close_event)
    self.state.signals.tab_change.connect(self.change_tab_to)


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
    tab_widget = self.view.tabs
    for index in range(tab_widget.count()):
      if tab_widget.tabText(index) == name:
        tab_widget.setCurrentIndex(index)
        return
