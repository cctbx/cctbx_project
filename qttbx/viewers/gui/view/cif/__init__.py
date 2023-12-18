
from PySide2.QtCore import Signal

from ..widgets.tab import GUITabWidget
from ..restraints.bonds import BondsTabView
from ..restraints.angles import AnglesTabView
from .cif_files import CifFileTabView
from .cif_browser import CifBrowserTabView


class CifTabView(GUITabWidget): # Top view
  trigger_processing = Signal(object)
  def __init__(self,parent=None):
    super().__init__(parent=parent)


    # Files
    self.files = CifFileTabView(parent=self)
    self.addTab(self.files, "Files")

    # Browser
    self.browser = CifBrowserTabView(parent=self)
    self.addTab(self.browser, "Browser")
