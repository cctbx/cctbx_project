from PySide2.QtWidgets import QSizePolicy, QVBoxLayout, QPushButton, QSpacerItem
from PySide2.QtWebEngineWidgets import QWebEngineView, QWebEnginePage
from PySide2.QtCore import Signal
from PySide2.QtGui import QDragEnterEvent, QDropEvent


from ..widgets.tab import GUITab
from ..widgets.selection_controls import SelectionControlsView



class ChimeraXTabView(GUITab):
  """
  The QT GUI Tab for the viewer
  """
  def __init__(self,parent=None):
    super().__init__(parent)


    self.layout = QVBoxLayout()
    self.setLayout(self.layout)


    spacer = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)
    self.layout.addSpacerItem(spacer)

    self.button_start = QPushButton("Start ChimeraX")
    self.button_start.setMaximumWidth(128)
    self.layout.addWidget(self.button_start)

    self.selection_controls = SelectionControlsView()
    self.layout.addWidget(self.selection_controls)

