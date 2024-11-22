
from PySide2.QtCore import Signal
from PySide2.QtGui import QGuiApplication
from PySide2.QtWidgets import QMainWindow

from qttbx.viewers.gui.view.molstar_base import MolstarTabView
from qttbx.viewers.gui.view.widgets.console import JupyterTabWidget
from qttbx.viewers.gui.view.widgets.tab import GUITabWidget

class BaseAppView(QMainWindow):
  signal_close = Signal()

  def __init__(self,parent=None,params=None):
    super().__init__(parent=parent)
    self.params = params
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


    self.tabs = GUITabWidget(parent=self)
    self.tabs.tabBar().childSignal.connect(self.child_window_handler)
    self.setCentralWidget(self.tabs)

    # Consoles
    self.consoles = GUITabWidget(parent=self)
    # Python console subtab
    self.python_console = JupyterTabWidget(parent=self.consoles)
    self.consoles.addTab(self.python_console, "Python")
    self.tabs.addTab(self.consoles,"Console")

    # Not visible by default
    if self.params:
      if not self.params.console:
        self.tabs.toggle_tab_visible("Console",show=False)



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
