from PySide2.QtCore import Qt
from PySide2.QtGui import QPalette, QColor
from PySide2.QtWidgets import (
    QScrollArea,
    QVBoxLayout,
    QWidget
)

class CustomScrollArea(QScrollArea):
  def __init__(self, widget, parent=None):
    super().__init__(parent)
    self.setWidget(widget)
    self.setWidgetResizable(True)
    viewport_palette = QPalette()
    viewport_palette.setColor(QPalette.Background, Qt.transparent)
    self.viewport().setPalette(viewport_palette)
    self.viewport().setAutoFillBackground(True)




class EntryContainer(QWidget):
  """
  A container for ALL entries
  """
  def __init__(self, parent=None):
    super().__init__(parent)
    self.layout = QVBoxLayout(self)
    palette = QPalette()
    palette.setColor(QPalette.Background, QColor(255, 255, 255))
    self.setPalette(palette)
    self.setAutoFillBackground(True)
    self.layout.setAlignment(Qt.AlignTop)
    self.layout.setContentsMargins(5, 5, 5, 5)
    self.setLayout(self.layout)



class ScrollableListView(QWidget):
  def __init__(self, parent=None):
    super().__init__(parent)
    self.parent_explicit = parent
    self.layout = QVBoxLayout(self)
    self.container = EntryContainer(parent=self)
    scroll = CustomScrollArea(self.container)
    self.layout.addWidget(scroll)
