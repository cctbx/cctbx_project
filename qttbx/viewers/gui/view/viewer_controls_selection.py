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
from qttbx.viewers.gui.view.viewer_controls_base import *
from qttbx.viewers.gui.view.widgets.toggles import ToggleIconButton


class ViewerControlsSelectionView(ViewerControlsBaseView):
  def __init__(self, parent=None):
    super().__init__(parent=parent)
    print("Init ViewerControlsSelectionView")
    # Picking level
    self.picking_level = QComboBox()
    self.picking_level.addItem("Atom")
    self.picking_level.addItem("Residue")
    self.picking_level.setToolTip("Pick atoms or residues")


    # Start manual selection
    self.selector_toggle = ToggleIconButton(qta.icon("mdi.crosshairs"),qta.icon("mdi.crosshairs-off"))
    self.selector_toggle.setToolTip("Toggle selecting")
    self.selector_toggle.setMaximumWidth(self._all_button_width)

    # Clear selection
    self.button_cancel = QPushButton()
    icon = qta.icon("mdi.cancel")
    self.button_cancel.setIcon(icon)
    self.button_cancel.setToolTip("Clear selection")
    self.button_cancel.setMaximumWidth(self._all_button_width)  

    # Selection Add
    self.add_selection_button = QPushButton()
    icon = qta.icon("mdi.plus")
    self.add_selection_button.setIcon(icon)
    self.add_selection_button.setToolTip("Add selection to list")
    self.add_selection_button.setMaximumWidth(self._all_button_width)    
    

    selection_components = [
      ("Select",self.selector_toggle),
      ("Clear", self.button_cancel),
      ("Picking",self.picking_level),
      ("Add", self.add_selection_button),
      ("",spacer_min()),
    ]
    
    self.selection_panel = form_panel(selection_components, "Selection:")
    self.main_layout.addWidget(separator())
    self.main_layout.addLayout(self.selection_panel)

