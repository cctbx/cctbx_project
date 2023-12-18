from pathlib import Path

from PySide2.QtWidgets import QHBoxLayout, QVBoxLayout, QLabel, QPushButton
from PySide2.QtGui import QIcon

from ..widgets.tab import GUITabWidget, GUITab


class PreprocessTab(GUITab):
  def __init__(self,parent=None,title="Preprocess",action="Process",tool_tip="Process to continue"):
    super().__init__(parent=parent)
    self.layout = QVBoxLayout()
    self.setLayout(self.layout)
    self.title = title
    header_layout = QHBoxLayout()
    label = QLabel(title)
    current_font = label.font()
    current_font.setPointSize(16)
    current_font.setBold(False)
    label.setFont(current_font)
    header_layout.addWidget(label)
    header_layout.addStretch(1)
    self.header_layout = header_layout

    # Process

    # Process text
    label_process = QLabel(action)
    label_process.setFont(current_font)
    header_layout.addWidget(label_process)

    # Button
    self.process_button = QPushButton()
    icon_path = Path(__file__).parent / '../assets/icons/material/gears.svg'
    process_icon = QIcon(str(icon_path))
    self.process_button.setIcon(process_icon)
    self.process_button.setMaximumSize(50, 50)
    self.process_button.setContentsMargins(10, 10, 0, 0)
    self.process_button.setToolTip(tool_tip)
    header_layout.addWidget(self.process_button)


    # # Load button
    # self.load_button = QPushButton()
    # icon_path = Path(__file__).parent / '../assets/icons/material/plus.svg'
    # load_icon = QIcon(str(icon_path))
    # self.load_button.setIcon(load_icon)
    # self.load_button.setMaximumSize(50, 50)
    # self.load_button.setContentsMargins(10, 10, 0, 0)
    # header_layout.addWidget(self.load_button)

    self.layout.insertLayout(0, header_layout)
