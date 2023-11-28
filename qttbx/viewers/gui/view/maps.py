from pathlib import Path

from PySide2.QtWidgets import QLabel, QPushButton, QHBoxLayout
from PySide2.QtGui import QIcon

from .widgets.scroll_list import ScrollableListView
from .widgets.scroll_entry import ScrollEntryView
from .widgets import ISOWidget, OpacityWidget


class MapEntryView(ScrollEntryView):
  def __init__(self,parent=None):
    super().__init__(parent=parent)
    self.iso_label = QLabel("ISO: ")
    self.iso_label.setFixedWidth(80)
    self.layout.addWidget(self.iso_label)

    
    self.iso_widget = ISOWidget(parent=self)
    self.layout.addWidget(self.iso_widget)

    # Opacity
   # Color picking widget
    self.button_opacity = QPushButton()
    icon_path = Path(__file__).parent / 'assets/icons/material/opacity.svg'
    icon = QIcon(str(icon_path))
    self.button_opacity.setIcon(icon)
    self.button_opacity.setToolTip("Opacity")
    self.button_opacity.setFixedSize(self._all_button_width,self._all_button_height)
    #button_color.setContentsMargins(0,0,0,0)
    #button_color.setMaximumSize(QSize(maxs2,maxs2))
    self.layout.addWidget(self.button_opacity)
    self.button_opacity.clicked.connect(self.show_opacity)

    self.opacity_widget = OpacityWidget(parent=self)

    # Close
    self.button_close = QPushButton()
    icon_path = Path(__file__).parent / 'assets/icons/material/close.svg'
    icon = QIcon(str(icon_path))
    self.button_close.setIcon(icon)
    self.button_close.setToolTip("Remove")
    self.button_close.setFixedSize(self._all_button_width,self._all_button_height)
    self.layout.addWidget(self.button_close)


  def show_opacity(self):
      global_pos = self.button_opacity.mapToGlobal(self.button_opacity.pos())
      self.opacity_widget.move(global_pos)  # Position the widget
      self.opacity_widget.show()  # Show the widget








class MapListView(ScrollableListView):
  def __init__(self,parent=None,title="Maps"):
    super().__init__(parent=parent)
    header_layout = QHBoxLayout()
    label = QLabel(title)
    current_font = label.font()
    current_font.setPointSize(16)
    current_font.setBold(False)
    label.setFont(current_font)
    
    self.load_button = QPushButton()
    icon_path = Path(__file__).parent / './assets/icons/material/plus.svg'
    load_icon = QIcon(str(icon_path))
    self.load_button.setIcon(load_icon)
    self.load_button.setMaximumSize(50, 50)
    self.load_button.setContentsMargins(10, 10, 0, 0) 
    header_layout.addWidget(label)
    header_layout.addWidget(self.load_button)
    
    self.layout.insertLayout(0, header_layout)
    


