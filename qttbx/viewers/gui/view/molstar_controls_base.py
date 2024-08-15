"""
The view which manages the minimum base level controls present in the molstar GUI tab
"""
from pathlib import Path

import qtawesome as qta
from PySide2.QtGui import QIcon
from PySide2.QtCore import Qt, QEvent
from PySide2.QtCore import QObject, Signal
from PySide2.QtGui import QIcon
from PySide2.QtWidgets import (
    QApplication,
    QComboBox,
    QDialog,
    QFrame,
    QHBoxLayout,
    QLabel,
    QLayoutItem,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QSpacerItem,
    QVBoxLayout,
    QWidget
)

from qttbx.viewers.gui.view.widgets.history_line_edit import HistoryLineEdit
from qttbx.viewers.gui.view.selection_controls import ScrollableHierarchyWidget, SearchSelectDialog

# Enable high DPI scaling
QApplication.setAttribute(Qt.AA_EnableHighDpiScaling)
QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps)


def spacer_max():
  return QSpacerItem(5, 5, QSizePolicy.Expanding, QSizePolicy.Maximum)
def spacer_min():
  return QSpacerItem(20, 20, QSizePolicy.Expanding, QSizePolicy.Minimum)

class MolstarControlsBaseView(QWidget):
  def __init__(self, parent=None):
    super().__init__(parent=parent)
    self._all_button_width = 48
    self._all_button_height = 48
    self._all_icon_size = 16



    # Model label
    self.active_model_label = QLabel()


    # Selection Search
    #####################
    self.button_search = QPushButton()
    icon = qta.icon("mdi.magnify")
    self.button_search.setIcon(icon)
    self.button_search.setToolTip("Specify a selection")
    self.button_search.setMaximumWidth(self._all_button_width)    

    # Clear selection
    self.button_cancel = QPushButton()
    icon = qta.icon("mdi.cancel")
    self.button_cancel.setIcon(icon)
    self.button_cancel.setToolTip("Specify a selection")
    self.button_cancel.setMaximumWidth(self._all_button_width)  

    # search box goes in popup
    self.selection_edit = HistoryLineEdit()

    
    def separator():
      separator = QFrame()
      separator.setFrameShape(QFrame.VLine)
      separator.setFrameShadow(QFrame.Sunken)
      return separator

    active_model_components = [
      (None,self.active_model_label),
      ("",spacer_min())
    ]
    selection_components = [
      ("Search", self.button_search),
      ("Clear", self.button_cancel),
      ("",spacer_min()),
    ]
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
      for label_text, widget in components:
        vbox = QVBoxLayout()
        if label_text:
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

