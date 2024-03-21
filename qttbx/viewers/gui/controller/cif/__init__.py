
from PySide2.QtWidgets import QApplication, QMessageBox

from PySide2.QtCore import QUrl, QThread, Signal, Slot, QObject, QThreadPool, QRunnable
from ..scroll_entry import ScrollEntryController
from ...view.widgets.scroll_entry import ScrollEntryView
from ..scroll_list import ScrollableListController
from ..controller import Controller
from ...state.restraints import Restraints
from ..data import DataTabController
from .cif_files import CifFileListController
from .cif_browser import CifBrowserController

import time

class CifTabController(Controller):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)

    self.files = CifFileListController(parent=self,view=self.view.files.list_view)
    self.browser = CifBrowserController(parent=self,view=self.view.browser)
