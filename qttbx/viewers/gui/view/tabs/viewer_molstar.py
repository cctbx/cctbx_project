from PySide2.QtWidgets import QSizePolicy, QVBoxLayout
from PySide2.QtWebEngineWidgets import QWebEngineView, QWebEnginePage
from PySide2.QtCore import Signal, QEventLoop
from PySide2.QtGui import QDragEnterEvent, QDropEvent


from ..widgets.tab import GUITab
from ..widgets.selection_controls import SelectionControlsView

from PySide2.QtWebEngineWidgets import QWebEngineView


class WebEnginePage(QWebEnginePage):
  """
  Subclassing this allows access to javascript console output
  """
  console_message = Signal(str, int, str)
  def __init__(self, *args, **kwargs):
      super().__init__(*args, **kwargs)


  def javaScriptConsoleMessage(self, level, msg, line, sourceID):
    print(f"JS console message: {msg}, line: {line}, sourceID: {sourceID}") #get as stdout
    self.console_message.emit(msg, line, sourceID) # get for a gui tab




class MolstarWebEngineView(QWebEngineView):
  """
  Subclassing this disables the "File drop" behavior of QWebEngineView
  """

  def dragEnterEvent(self, event: QDragEnterEvent):
    event.ignore()

  def dropEvent(self, event: QDropEvent):
    event.ignore()

  def runJavaScript(self,script,custom_callback=None):
    return self.page().runJavaScript(script,0,custom_callback)

  def runJavaScriptSync(self, script, custom_callback=None):
    """
    Run javascript syncronously, meaning Python waits for it to finish.
    Not using this should be preferred.
    """
    self.loop = QEventLoop()  # Create an event loop
    self.custom_callback = custom_callback
    self.page().runJavaScript(script, 0, self.onJavaScriptResult)
    self.loop.exec_()  # Block until loop.quit() is called

  def onJavaScriptResult(self, result):
    if self.custom_callback:
        self.custom_callback(result)  # Run the custom callback with the result
    self.loop.quit()  # Quit the event loop to unblock

class ViewerTabView(GUITab):
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

    self.selection_controls = SelectionControlsView()
    self.layout.addWidget(self.selection_controls)


    # Set the layout for the whole viewer
    self.setLayout(self.layout)
