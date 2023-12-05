import time

from PySide2.QtCore import Signal
from PySide2.QtWidgets import QMainWindow
from PySide2.QtGui import QGuiApplication

from ..tabs.viewer_chimerax import ChimeraXTabView
from ..tabs.viewer_molstar import ViewerTabView
from ..tabs.console import JupyterTabWidget, JSConsoleTab
from ..tabs.selection import SelectionsTabView
from ..tabs.data import DataTabView
from ..tabs.sites import SitesTabView
from ..cif import CifTabView
#from ..tabs.restraints import RestraintsTopTabView
from ..tabs.restraints_table import RestraintsTableTopTabView
from ..tabs.qscore import QscoreTab
from ..widgets.tab import GUITabWidget



class DemoView(QMainWindow):
  signal_close = Signal()

  def __init__(self,parent=None):
    super().__init__(parent=parent)
    self._has_child_window = False
    

    # set title
    self.setWindowTitle("Phenix Viewer")
    #self.setWindowFlags(Qt.FramelessWindowHint)
    # titlebar = SnapTitleBar(self)
    # self.setMenuWidget(titlebar)

    # Retrieve screen size
    screen = self.screen().availableGeometry()
    screen_width = screen.width()
    screen_height = screen.height()

    # Set window size based on screen size (e.g., 80% of screen dimensions)
    w,h = screen_width * 0.5, screen_height * 1.0
    self.setGeometry(screen_width // 2, screen_height,w,h)

    # Main Level Components
    self.tabs = GUITabWidget(parent=self)
    self.tabs.tabBar().childSignal.connect(self.child_window_handler) 
    self.setCentralWidget(self.tabs)


    # Viewer
    self.viewer_tab_view = ViewerTabView(parent=self)
    self.tabs.insertTab(1,self.viewer_tab_view, "Molstar")

    # ChimeraX Viewer
    self.chimerax_tab_view = ChimeraXTabView(parent=self)
    self.tabs.insertTab(0,self.chimerax_tab_view, "ChimeraX")

    # Selections
    self.selection_tab_view = SelectionsTabView(parent=self)
    self.tabs.addTab(self.selection_tab_view, "Selections")

    # Data
    # recieves datamanager from outside
    self.data_tab_view = DataTabView(parent=self)
    self.tabs.addTab(self.data_tab_view, "Files")

    # Sites
    self.sites_tab_view = SitesTabView(parent=self)
    self.tabs.addTab(self.sites_tab_view, "Atoms")

    # Cif
    self.cif_tab_view = CifTabView(parent=self)
    self.tabs.addTab(self.cif_tab_view, "CIF")

    # # Restraints
    # self.restraints_tab_view = RestraintsTopTabView(parent=self)
    # self.tabs.addTab(self.restraints_tab_view, "Restraints")

    # Restraints Table
    self.restraints_table_tab_view = RestraintsTableTopTabView(parent=self)
    self.tabs.addTab(self.restraints_table_tab_view, "Restraints")


    # Qscore
    # self.qscore_tab_view = QscoreTab(parent=self)
    # self.tabs.addTab(self.qscore_tab_view, "Qscore")

    # Consoles
    self.consoles = GUITabWidget(parent=self)
     # Python console subtab
    self.python_console = JupyterTabWidget(parent=self.consoles)
    self.consoles.addTab(self.python_console, "Python")

    # javascript console subtab
    self.javascript_console = JSConsoleTab(parent=self.consoles,web_view=self.viewer_tab_view.web_view)
    self.consoles.addTab(self.javascript_console, "Javascript")

    self.tabs.addTab(self.consoles,"Console")

  
  def child_window_handler(self,event):
    self._has_child_window = True
    pass #TODO: Setting window to half of screen causes it to freeze. Why?
    # if event == "created":
    #   QTimer.singleShot(3, lambda: self.setHalfWindow())

  def setHalfWindow(self):
    # make window the full right half of screen
    # Get screen geometry
    desktop = QGuiApplication.primaryScreen()
    screen_rect = desktop.availableGeometry()

    # Set main window geometry to take up the right half of the screen
    self.setGeometry(screen_rect.width() // 2, screen_rect.y(),
                  screen_rect.width() // 2, screen_rect.height())
    self.update()
    self.repaint()
    self.setEnabled(True)
    self.setFocus()


  def closeEvent(self, event):
    self.signal_close.emit()
    event.accept()


