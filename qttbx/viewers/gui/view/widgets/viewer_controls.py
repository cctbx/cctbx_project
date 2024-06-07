from pathlib import Path

import qtawesome as qta
from PySide2.QtWidgets import QComboBox, QDialog, QFrame, QHBoxLayout, QLabel, QPushButton, QSizePolicy, QVBoxLayout, QWidget
from PySide2.QtWidgets import QSpacerItem,QLayoutItem, QApplication, QWidget, QVBoxLayout, QListWidget, QListWidgetItem, QLabel, QScrollArea
from PySide2.QtGui import QIcon
from PySide2.QtCore import Qt, QEvent
from PySide2.QtCore import QObject, Signal
from PySide2.QtGui import QIcon

from .history_line_edit import HistoryLineEdit

# Enable high DPI scaling
QApplication.setAttribute(Qt.AA_EnableHighDpiScaling)
QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps)


class ViewerControlsView(QWidget):
  def __init__(self, parent=None):
    super().__init__(parent=parent)
    self._all_button_width = 48
    self._all_button_height = 48
    self._all_icon_size = 16



    # Model label
    self.active_model_label = QComboBox()



    # Geometry
    #################
    self.button_geo = QPushButton()
    self.button_geo.setCheckable(True)
    icon = qta.icon("mdi.cogs")
    self.button_geo.setIcon(icon)
    self.button_geo.setToolTip("Manage geometry")
    self.button_geo.setMaximumWidth(self._all_button_width)
    #self.geo_checkbox = GeoCheckBox("", self)
    #geo_hbox = QHBoxLayout()
    #geo_hbox.addWidget(self.button_geo)
    #geo_hbox.addWidget(self.geo_checkbox)

    # Use as restraint edit
    #####################
    self.button_edits = QPushButton()
    # icon_path = Path(__file__).parent / '../assets/icons/material/edit.svg'
    # icon = QIcon(str(icon_path))
    icon = qta.icon("mdi.pencil")
    self.button_edits.setIcon(icon)
    self.button_edits.setToolTip("Selection as geometry edit")
    self.button_edits.setMaximumWidth(self._all_button_width)

    # Load restraints
    self.button_restraints = QPushButton()
    self.button_restraints.setCheckable(True)

    #icon_path = Path(__file__).parent / '../assets/icons/material/measure.svg'
    #icon = QIcon(str(icon_path))
    icon = qta.icon("mdi.set-square")

    self.button_restraints.setIcon(icon)
    self.button_restraints.setToolTip("Load restraints from library")
    self.button_restraints.setMaximumWidth(self._all_button_width)


    def spacer_max():
      return QSpacerItem(5, 5, QSizePolicy.Expanding, QSizePolicy.Maximum)
    def spacer_min():
      return QSpacerItem(20, 20, QSizePolicy.Expanding, QSizePolicy.Minimum)
    # Selection Search
    #####################
    self.button_search = QPushButton()
    # icon_path = Path(__file__).parent / '../assets/icons/material/search.svg'
    # icon = QIcon(str(icon_path))
    icon = qta.icon("mdi.magnify")
    self.button_search.setIcon(icon)
    self.button_search.setToolTip("Specify a selection")
    self.button_search.setMaximumWidth(self._all_button_width)    

    # search box goes in popup
    self.selection_edit = HistoryLineEdit()

    # Selection pick
    #################
    # Version 1: two distinct buttons
    # self.button_select = QPushButton()
    # icon_path = Path(__file__).parent / '../assets/icons/material/target_searching.svg'
    # icon = QIcon(str(icon_path))
    # self.button_select.setIcon(icon)
    # self.button_select.setToolTip("Start selecting")
    # self.button_select.setMaximumWidth(self._all_button_width)    
 

    # Version 2: Checkable buttons
    self.button_select = QPushButton()
    self.button_select.setCheckable(True)
    #icon_path = Path(__file__).parent / '../assets/icons/material/target_searching.svg'
    #icon = QIcon(str(icon_path))
    icon = qta.icon("mdi.crosshairs")
    self.button_select.setIcon(icon)
    self.button_select.setToolTip("Toggle selecting")
    self.button_select.setMaximumWidth(self._all_button_width)   

    # New button creation func
    self.button_clear = QPushButton()
    #icon_path = Path(__file__).parent / '../assets/icons/material/cancel.svg'
    #icon = QIcon(str(icon_path))
    icon = qta.icon("mdi.close-circle-outline")
    self.button_clear.setIcon(icon)
    self.button_clear.setToolTip("Clear all selection")
    self.button_clear.setMaximumWidth(self._all_button_width)   

    # Selection Add
    #################
    self.add_selection_button = QPushButton()
    #icon_path = Path(__file__).parent / '../assets/icons/material/plus.svg'
    #load_icon = QIcon(str(icon_path))
    icon = qta.icon("mdi.plus")
    self.add_selection_button.setIcon(icon)
    self.add_selection_button.setToolTip("Add selection to list")
    self.add_selection_button.setMaximumWidth(self._all_button_width)    


    # Picking level
    #################
    self.picking_level = QComboBox()
    icon_path = Path(__file__).parent / '../assets/icons/other/WaterFull.svg'
    icon = QIcon(str(icon_path))
    self.picking_level.addItem("Residue")
    icon_path = Path(__file__).parent / '../assets/icons/other/WaterSingle.svg'
    icon = QIcon(str(icon_path))
    self.picking_level.addItem("Atom")
    self.picking_level.setToolTip("Pick atoms or residues")

   

    def separator():
      separator = QFrame()
      separator.setFrameShape(QFrame.VLine)
      separator.setFrameShadow(QFrame.Sunken)
      return separator

    active_model_components = {
      #"":spacer_max(),
      "Filename:": self.active_model_label,
      #"File":self.button_files,
      "Geo": self.button_geo,
      "Edits": self.button_edits,
      "Restraints": self.button_restraints,
      "":spacer_min(),
    }
    selection_components = {
      #"":spacer_max(),
      "Search": self.button_search,
      "Select": self.button_select,
      "Clear": self.button_clear,
      "Add": self.add_selection_button,
      "Picking": self.picking_level,
      "":spacer_min(),
    }
    self.settings_vbox = QVBoxLayout()
    self.settings_vbox.setSpacing(3)
    self.settings_vbox.setContentsMargins(3, 3, 3, 3)
    self.settings_vbox.setAlignment(Qt.AlignVCenter)
    main_layout = QHBoxLayout()


    def add_component(layout,component):
      if isinstance(component,(QHBoxLayout,QVBoxLayout)):
        layout.addLayout(component)
      elif isinstance(component,QLayoutItem):
        layout.addItem(component)
      else:
        layout.addWidget(component)

    def form_panel(components, component_text):
      comp_box = QVBoxLayout()
      label = QLabel(component_text)
      label.setStyleSheet("font-weight: bold;")
      comp_box.addWidget(label)
      hbox = QHBoxLayout()
      for label_text, widget in components.items():
        vbox = QVBoxLayout()
        label = QLabel(label_text)
        add_component(vbox,label)
        add_component(vbox,widget)
        hbox.addLayout(vbox)
      
      comp_box.addLayout(hbox)
      return comp_box

    self.model_panel = form_panel(active_model_components, "Active model:")
    self.selection_panel = form_panel(selection_components, "Selection:")

    main_layout.addLayout(self.model_panel)
    main_layout.addWidget(separator())
    main_layout.addLayout(self.selection_panel)
    self.settings_vbox.addLayout(main_layout)
    self.setLayout(self.settings_vbox)

    # hierarchy select widgets
    self.search_select_dialog = None # set by controller


  
class ScrollableHierarchyWidget(QWidget):
  selected_change = Signal(str) # the value that was changed to

  def __init__(self,parent=None):
    super().__init__(parent=parent)
    self.setFocusPolicy(Qt.StrongFocus)

    # Create the main layout
    self.layout = QVBoxLayout()

    # Create a QListWidget instance
    self.list_widget = QListWidget()
    self.list_widget.keyPressFilter = KeyPressFilter(self.list_widget)
    self.list_widget.installEventFilter(self.list_widget.keyPressFilter)

    # Create a scroll area and set the list widget as its widget
    self.scroll_area = QScrollArea()
    self.scroll_area.setWidgetResizable(True)
    self.scroll_area.setWidget(self.list_widget)

    # Add the scroll area to the main layout
    self.layout.addWidget(self.scroll_area)

    # Add a label to show the current selection
    self.label = QLabel("Selected: None")
    self.layout.addWidget(self.label)

    # Set the layout for the main window
    self.setLayout(self.layout)



  def update_list_widget(self, items):
    # Clear the current items
    self.list_widget.clear()
    # Add new items to the list widget
    for item in items:
      self.list_widget.addItem(QListWidgetItem(item))
  



# End edits dialog
class SearchSelectDialog(QDialog):
  def __init__(self,title="Selection search",parent=None):
    super().__init__(parent)
    self.setWindowTitle(title)
    self.list_width = 150
    self.chain_scroller = ScrollableHierarchyWidget(self)
    self.comp_scroller = ScrollableHierarchyWidget(self)
    self.seq_scroller = ScrollableHierarchyWidget(self)
    self.atom_scroller = ScrollableHierarchyWidget(self)
    layout = QVBoxLayout(self)
    # navigate selection
    hierarchy_label = QLabel("Hierarchy filters:")
    layout.addWidget(hierarchy_label)
    hierarchy_layout  = QHBoxLayout()
    # Chains
    chain_layout = QVBoxLayout()
    chain_label = QLabel("    Chain")
    chain_list = self.chain_scroller
    chain_list.setFixedWidth(self.list_width)
    chain_layout.addWidget(chain_label)
    chain_layout.addWidget(chain_list)
    hierarchy_layout.addLayout(chain_layout)

    # Residue components
    comp_layout = QVBoxLayout()
    comp_label = QLabel("   Residue")
    comp_list = self.comp_scroller
    comp_list.setFixedWidth(self.list_width)

    comp_layout.addWidget(comp_label)
    comp_layout.addWidget(comp_list)
    hierarchy_layout.addLayout(comp_layout)


    # Residues
    seq_layout = QVBoxLayout()
    seq_label = QLabel("   Sequence")
    seq_list = self.seq_scroller
    seq_list.setFixedWidth(self.list_width)

    seq_layout.addWidget(seq_label)
    seq_layout.addWidget(seq_list)
    hierarchy_layout.addLayout(seq_layout)

    # Atoms
    atom_layout = QVBoxLayout()
    atom_label = QLabel("   Atom")
    atom_list = self.atom_scroller
    atom_list.setFixedWidth(self.list_width)

    atom_layout.addWidget(atom_label)
    atom_layout.addWidget(atom_list)
    hierarchy_layout.addLayout(atom_layout)

    layout.addLayout(hierarchy_layout)



    # String Selection
    selection_edit = parent.selection_edit
    
    search_label = QLabel("Phenix selection string:")
    selection_edit.setPlaceholderText("Enter selection string...")
    selection_edit.setToolTip("Enter a Phenix selection string")
    layout.addWidget(search_label)
    layout.addWidget(selection_edit)

  

class KeyPressFilter(QObject):
  space_pressed = Signal()
  def __init__(self, parent=None):
    super().__init__(parent)

  def eventFilter(self, obj, event):
    if event.type() == QEvent.KeyPress:
      if event.key() == Qt.Key_Space:
        self.space_pressed.emit()
        return True  # Event is handled
    return super().eventFilter(obj, event)
