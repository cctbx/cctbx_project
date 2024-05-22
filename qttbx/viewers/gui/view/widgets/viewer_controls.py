from PySide2.QtWidgets import QWidget, QFrame, QSizePolicy, QVBoxLayout,QHBoxLayout, QPushButton
from PySide2.QtGui import QIcon
from PySide2.QtCore import Qt


from .history_line_edit import HistoryLineEdit
from .widgets import NoCheckComboBox
from .toggles import ToggleIconButton

from ..widgets.scroll_list import ScrollableListView
from ..widgets.scroll_entry import ScrollEntryView
from ..widgets.representation_select import RepresentationSelect
from ..widgets.toggles import ToggleIconButton
from ..widgets.checkbox import ConditionalCheckBox
from ..models import GeoCheckBox

from pathlib import Path

class ViewerControlsView(QWidget):
  def __init__(self,parent=None):
    super().__init__(parent=parent)
    self._all_button_width = 48
    self._all_button_height = 48

    #Start settings/selection area
    self.settings_vbox = QVBoxLayout()
    policy = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
    self.setSizePolicy(policy)
    #self.setFixedHeight(40)
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

    #############
    ## Model Hbox
    #############
    self.model_hbox = QHBoxLayout()

    # Open in folder
    self.button_files = QPushButton()
    icon_path = Path(__file__).parent / '../assets/icons/material/folder.svg'
    icon = QIcon(str(icon_path))
    self.button_files.setIcon(icon)
    self.button_files.setToolTip("Open containing folder")
    self.button_files.setMaximumWidth(self._all_button_width)
    self.model_hbox.addWidget(self.button_files)

   # vertical separator
    separator = QFrame()
    separator.setFrameShape(QFrame.VLine)
    separator.setFrameShadow(QFrame.Sunken)
    self.model_hbox.addWidget(separator)

    self.button_geo = QPushButton()
    icon_path = Path(__file__).parent / '../assets/icons/material/gears.svg'
    icon = QIcon(str(icon_path))
    self.button_geo.setIcon(icon)
    self.button_geo.setToolTip("Manage geometry")
    # self.selection_hbox.addWidget(self.checkbox_persist)
    self.model_hbox.addWidget(self.button_geo)
    self.button_geo.setMaximumWidth(self._all_button_width)

    self.geo_checkbox = GeoCheckBox("Geo", self)
    #self.geo_checkbox.stateChanged.connect(self.showState)
    self.model_hbox.addWidget(self.geo_checkbox)


    # vertical separator
    separator = QFrame()
    separator.setFrameShape(QFrame.VLine)
    separator.setFrameShadow(QFrame.Sunken)
    self.model_hbox.addWidget(separator)

    # Visibility
    on_icon_path = Path(__file__).parent / '../assets/icons/material/eye_open.svg'
    off_icon_path = Path(__file__).parent / '../assets/icons/material/eye_closed.svg'
    self.button_viz = ToggleIconButton(on_icon_path, off_icon_path, parent=self)
    self.button_viz.setToolTip("Toggle visibility")
    self.button_viz.setMaximumWidth(self._all_button_width)
    #button_color.setContentsMargins(0,0,0,0)
    self.model_hbox.addWidget(self.button_viz)



    # Color theme widget # TODO: fix this
    # self.button_theme = ColorThemeButton(parent=self)
    # self.button_theme.button.setFixedSize(self._all_button_width,self._all_button_height)
    # self.layout.addWidget(self.button_theme)

    # Representations
    self.button_rep = RepresentationSelect(parent=self)
    self.button_rep.button.setMaximumWidth(self._all_button_width)
    self.model_hbox.addWidget(self.button_rep)

    # Color picking widget
    self.button_color = QPushButton()
    icon_path = Path(__file__).parent / '../assets/icons/material/paint_bucket.svg'
    icon = QIcon(str(icon_path))
    self.button_color.setIcon(icon)
    self.button_color.setToolTip("Color fill")
    self.button_color.setMaximumWidth(self._all_button_width)
    #button_color.setContentsMargins(0,0,0,0)
    #button_color.setMaximumSize(QSize(maxs2,maxs2))
    self.model_hbox.addWidget(self.button_color)


    # As restraint
    self.button_restraint = QPushButton()
    icon_path = Path(__file__).parent / '../assets/icons/material/measure.svg'
    icon = QIcon(str(icon_path))
    self.button_restraint.setIcon(icon)
    self.button_restraint.setToolTip("Selection as geometry edit")
    self.button_restraint.setMaximumWidth(self._all_button_width)
    #button_color.setContentsMargins(0,0,0,0)
    #button_color.setMaximumSize(QSize(maxs2,maxs2))
    self.model_hbox.addWidget(self.button_restraint)

    self.settings_vbox.addLayout(self.model_hbox)

    #self.row_2_hbox = QHBoxLayout()
    # from ..models import ModelEntryView
    # scroll_entry_view = ModelEntryView()
    #self.row_2_hbox.addWidget(self.selection_hbox)
    # # Visibility
    # on_icon_path = Path(__file__).parent / 'assets/icons/material/eye_open.svg'
    # off_icon_path = Path(__file__).parent / 'assets/icons/material/eye_closed.svg'
    # self.button_viz = ToggleIconButton(on_icon_path, off_icon_path, parent=self)
    # self.button_viz.setToolTip("Toggle visibility")
    # self.button_viz.setFixedSize(self._all_button_width,self._all_button_height)
    # #button_color.setContentsMargins(0,0,0,0)
    # self.row_2_hbox.addWidget(self.button_viz)



    # # Color theme widget # TODO: fix this
    # # self.button_theme = ColorThemeButton(parent=self)
    # # self.button_theme.button.setFixedSize(self._all_button_width,self._all_button_height)
    # # self.layout.addWidget(self.button_theme)

    # # Representations
    # self.button_rep = RepresentationSelect(parent=self)
    # self.button_rep.button.setFixedSize(self._all_button_width,self._all_button_height)
    # self.row_2_hbox.addWidget(self.button_rep)

    # # Color picking widget
    # self.button_color = QPushButton()
    # icon_path = Path(__file__).parent / 'assets/icons/material/paint_bucket.svg'
    # icon = QIcon(str(icon_path))
    # self.button_color.setIcon(icon)
    # self.button_color.setToolTip("Color fill")
    # self.button_color.setFixedSize(self._all_button_width,self._all_button_height)
    # #button_color.setContentsMargins(0,0,0,0)
    # #button_color.setMaximumSize(QSize(maxs2,maxs2))
    # self.row_2_hbox.addWidget(self.button_color)

    # # Create the second vertical separator
    # separator2 = QFrame()
    # separator2.setFrameShape(QFrame.VLine)
    # separator2.setFrameShadow(QFrame.Sunken)
    # self.row_2_hbox.addWidget(separator2)

    # # Open in folder
    # self.button_files = QPushButton()
    # icon_path = Path(__file__).parent / 'assets/icons/material/folder.svg'
    # icon = QIcon(str(icon_path))
    # self.button_files.setIcon(icon)
    # self.button_files.setToolTip("Open containing folder")
    # #button_color.setMaximumSize(50, 50)
    # self.button_files.setFixedSize(self._all_button_width,self._all_button_height)
    # self.row_2_hbox.addWidget(self.button_files)

    # # vertical separator
    # separator = QFrame()
    # separator.setFrameShape(QFrame.VLine)
    # separator.setFrameShadow(QFrame.Sunken)
    # self.row_2_hbox.addWidget(separator)

    # # Geometry
    # self.button_restraints = QPushButton()
    # icon_path = Path(__file__).parent / 'assets/icons/material/gears.svg'
    # icon = QIcon(str(icon_path))
    # self.button_restraints.setIcon(icon)
    # self.button_restraints.setToolTip("Manage restraints")
    # #button_color.setMaximumSize(50, 50)
    # self.button_restraints.setFixedSize(self._all_button_width,self._all_button_height)
    # self.row_2_hbox.addWidget(self.button_restraints)

    # self.geo_checkbox = GeoCheckBox("Geo", self)
    # #self.geo_checkbox.stateChanged.connect(self.showState)
    # self.row_2_hbox.addWidget(self.geo_checkbox)


    # # vertical separator
    # separator = QFrame()
    # separator.setFrameShape(QFrame.VLine)
    # separator.setFrameShadow(QFrame.Sunken)
    # self.row_2_hbox.addWidget(separator)

    # # Close
    # self.button_close = QPushButton()
    # icon_path = Path(__file__).parent / 'assets/icons/material/close.svg'
    # icon = QIcon(str(icon_path))
    # self.button_close.setIcon(icon)
    # self.button_close.setToolTip("Remove")
    # self.button_close.setFixedSize(self._all_button_width,self._all_button_height)
    # self.row_2_hbox.addWidget(self.button_close)



    # Add the settings_vbox to the top-level layout
    self.setLayout(self.settings_vbox)
