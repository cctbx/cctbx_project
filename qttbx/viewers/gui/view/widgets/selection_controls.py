from pathlib import Path

from PySide2.QtGui import QIcon
from PySide2.QtCore import Qt
from PySide2.QtWidgets import (
    QFrame,
    QHBoxLayout,
    QPushButton,
    QSizePolicy,
    QVBoxLayout,
    QWidget
)

from .history_line_edit import HistoryLineEdit
from .widgets import NoCheckComboBox
from .toggles import ToggleIconButton


class SelectionControlsView(QWidget):
  def __init__(self,parent=None):
    super().__init__(parent=parent)

    #Start settings/selection area
    self.settings_vbox = QVBoxLayout()
    policy = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
    self.setSizePolicy(policy)
    self.setFixedHeight(40)
    self.settings_vbox.setSpacing(3)
    self.settings_vbox.setContentsMargins(3,3,3,3)
    self.settings_vbox.setAlignment(Qt.AlignVCenter)

    # Create the first QHBoxLayout for selection
    self.selection_hbox = QHBoxLayout()
    #selection_label = QLabel("Selection:")
    #selection_label.setToolTip("Enter a Phenix selection string")
    self.selection_edit = HistoryLineEdit()
    self.selection_edit.setPlaceholderText("Enter selection string...")
    self.selection_edit.setToolTip("Enter a Phenix selection string")

    # Add widgets to the selection_hbox
    #self.selection_hbox.addWidget(selection_label)
    self.selection_hbox.addWidget(self.selection_edit)


    icon_path_checked = Path(__file__).parent / '../assets/icons/material/target_searching.svg'
    icon_path_unchecked = Path(__file__).parent / '../assets/icons/material/target_stop_searching.svg'
    self.selector_toggle = ToggleIconButton(str(icon_path_checked),str(icon_path_unchecked))
    self.selection_hbox.addWidget(self.selector_toggle)

    # # start selecting
    # self.start_selecting = QPushButton()
    # icon_path = Path(__file__).parent / '../assets/icons/material/target_searching.svg'
    # icon = QIcon(str(icon_path))
    # self.start_selecting.setIcon(icon)
    # self.start_selecting.setToolTip("Start selecting")
    # self.selection_hbox.addWidget(self.start_selecting)


    # # stop selecting
    # self.button_deselect = QPushButton()
    # icon_path = Path(__file__).parent / '../assets/icons/material/target_stop_searching.svg'
    # icon = QIcon(str(icon_path))
    # self.button_deselect.setIcon(icon)
    # self.button_deselect.setToolTip("Stop selecting, clear selections")
    # # self.selection_hbox.addWidget(self.checkbox_persist)
    # self.selection_hbox.addWidget(self.button_deselect)

    # Add selection button
    self.load_button = QPushButton()
    icon_path = Path(__file__).parent / '../assets/icons/material/plus.svg'
    load_icon = QIcon(str(icon_path))
    self.load_button.setIcon(load_icon)
    self.load_button.setToolTip("Keep selection")
    self.selection_hbox.addWidget(self.load_button)


    # Create a combobox for picking
    self.combo_box = NoCheckComboBox()

    # view = NoCheckListView()
    # combo_box.setView(view)

    # Adding items with icons to the QComboBox
    icon_path = Path(__file__).parent / '../assets/icons/other/WaterFull.svg'
    icon = QIcon(str(icon_path))
    self.combo_box.addItem(icon, "")
    icon_path = Path(__file__).parent / '../assets/icons/other/WaterSingle.svg'
    icon = QIcon(str(icon_path))
    self.combo_box.addItem(icon, "")
    self.combo_box.setToolTip("Pick atoms or residues")
    # Adjust icon size
    # current_size = self.combo_box.iconSize()
    # scale_factor = 1.5
    # scaled_size = QSize(current_size.width() * scale_factor, current_size.height() * scale_factor)
    # self.combo_box.setIconSize(scaled_size)


    # Add the QComboBox to the layout
    self.selection_hbox.addWidget(self.combo_box)


    # Create the vertical separator
    separator = QFrame()
    separator.setFrameShape(QFrame.VLine)
    separator.setFrameShadow(QFrame.Sunken)

    self.selection_hbox.addWidget(separator)

    # Clear
    self.button_clear = QPushButton()
    icon_path = Path(__file__).parent / '../assets/icons/material/reset.svg'
    icon = QIcon(str(icon_path))
    self.button_clear.setIcon(icon)
    self.button_clear.setToolTip("Reset. Close everything, start fresh.")
    # self.selection_hbox.addWidget(self.checkbox_persist)
    self.selection_hbox.addWidget(self.button_clear)


    # Add the selection_hbox to the settings_vbox
    self.settings_vbox.addLayout(self.selection_hbox)

    # Add the settings_vbox to the top-level layout
    self.setLayout(self.settings_vbox)
