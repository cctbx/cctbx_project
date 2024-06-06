from PySide2.QtWidgets import QSizePolicy, QVBoxLayout, QPushButton, QSpacerItem, QLabel, QLineEdit


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

    # Label
    self.label = QLabel("Executable path", self)
    self.layout.addWidget(self.label)

    # Text Input Widget
    self.textInput = QLineEdit(self)
    self.textInput.setMaximumWidth(256)
    self.layout.addWidget(self.textInput)

    self.button_start = QPushButton("Start ChimeraX")
    self.button_start.setMaximumWidth(128)
    self.layout.addWidget(self.button_start)


    self.selection_controls = SelectionControlsView()
    self.layout.addWidget(self.selection_controls)
