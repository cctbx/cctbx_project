"""
The view for the molstar web app. 
"""

from PySide2.QtWebEngineWidgets import QWebEngineView, QWebEnginePage
from PySide2.QtCore import Signal, QEventLoop, QObject, QUrl
from PySide2.QtGui import QDragEnterEvent, QDropEvent
from PySide2.QtWebEngineWidgets import QWebEngineView
from PySide2.QtWidgets import (
    QSizePolicy,
    QVBoxLayout
)

from qttbx.viewers.gui.view.widgets.tab import GUITab
from qttbx.viewers.gui.view.viewer_controls_base import ViewerControlsBaseView

class Signals(QObject):
  sync_signal = Signal()


class WebEnginePage(QWebEnginePage):
  """
  Subclassing this allows access to javascript console output
  """
  console_message = Signal(str, int, str)
  def __init__(self, *args, **kwargs):
      super().__init__(*args, **kwargs)
      self.debug = True

  def log(self,*args):
    if self.debug:
        print(*args)

  def javaScriptConsoleMessage(self, level, msg, line, sourceID):
    self.log(f"JS console message: {msg}, line: {line}, sourceID: {sourceID}") #get as stdout
    self.console_message.emit(msg, line, sourceID) # get for a gui tab



class MolstarWebEngineView(QWebEngineView):
    """
    Subclassing this disables the "File drop" behavior of QWebEngineView
    """
    def __init__(self,*args,**kwargs):
       super().__init__(*args,**kwargs)
       self.signals = Signals()
    
    def set_url(self,url):
        self.setUrl(QUrl(url))
    def dragEnterEvent(self, event: QDragEnterEvent):
        event.ignore()

    def dropEvent(self, event: QDropEvent):
        event.ignore()

    def runJavaScript(self, script, custom_callback=None):
        if custom_callback is None:
            custom_callback = lambda x: None
        return self.page().runJavaScript(script, 0, custom_callback)

    def runJavaScriptSync(self, script, custom_callback=None):
        """
        Run javascript synchronously, meaning Python waits for it to finish.
        """
        self.js_result = None  # Initialize a placeholder for the JS execution result
        self.custom_callback = custom_callback

        # Define a temporary callback that stores the result and optionally calls a custom callback
        def temp_callback(result):
            self.js_result = result  # Store the JavaScript execution result
            if self.custom_callback:
                self.custom_callback(result)
            self.loop.quit()  # Quit the event loop

        self.loop = QEventLoop()  # Create an event loop
        self.page().runJavaScript(script, 0, temp_callback)
        self.loop.exec_()  # Block until loop.quit() is called in temp_callback

        self.signals.sync_signal.emit()
        return self.js_result  # Return the JavaScript execution result

    def onJavaScriptResult(self, result):
        if self.custom_callback:
            self.custom_callback(result)  # Run the custom callback with the result
        self.loop.quit()  # Quit the event loop

class MolstarTabView(GUITab):
  """
  The QT GUI Tab for the viewer
  """
  def __init__(self,parent=None):
    super().__init__(parent)


    # Start GUI Components
    self.layout = QVBoxLayout()
    self.layout.setSpacing(5)
    self.layout.setContentsMargins(5,5,5,5)

    # Webview setup

    self.web_view = MolstarWebEngineView()
    self.web_page = WebEnginePage(self.web_view)
    self.web_view.setPage(self.web_page)

    self.web_view.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
    self.layout.addWidget(self.web_view)

    self.viewer_controls = ViewerControlsBaseView()

    # Set the layout for the whole viewer
    self.setLayout(self.layout)


