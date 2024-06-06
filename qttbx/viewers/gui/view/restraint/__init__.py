
from PySide2.QtCore import Signal

from ..widgets.tab import GUITabWidget
from .restraint_files import RestraintFileTabView
from .restraint_browser import RestraintBrowserTabView


class RestraintTabView(GUITabWidget): # Top view
  trigger_processing = Signal(object)
  def __init__(self,parent=None):
    super().__init__(parent=parent)


    # Files
    self.files = RestraintFileTabView(parent=self)
    self.addTab(self.files, "Files")

    # Browser
    self.browser = RestraintBrowserTabView(parent=self)
    self.addTab(self.browser, "Browser")
