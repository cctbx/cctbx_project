from pathlib import Path
import platform
import subprocess
import time

from PySide2.QtCore import Slot
from PySide2.QtCore import Qt
from PySide2.QtWidgets import (
    QFileDialog,
    QMenu,
    QPushButton
)

import iotbx
from mmtbx.monomer_library.pdb_interpretation import grand_master_phil_str

from .controller import Controller

class ViewerControlsBaseController(Controller):
  """
  A controller which manages the minimum base level controls present in the molstar GUI tab
  """
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)

    self.state.signals.active_change.connect(self.active_model_change)

    # Picking level
    self.picking_level = 'atom' # internal state
    self.view.picking_level.currentTextChanged.connect(self.picking_level_change) # user initiates new level
    self.state.signals.picking_level.connect(self.picking_level_change) # update view to reflect new level

    # Enable return key to execute selection
    self.view.selector_toggle.clicked.connect(self.start_selecting)
    self.view.button_cancel.clicked.connect(self.viewer.deselect_all)




  @property
  def viewer(self):
    return self.parent


  @property
  def picking_level_displayed(self):
    current_text = self.view.picking_level.currentText().lower()
    return current_text

  def picking_level_change(self, picking_level):
    picking_level = picking_level.lower()
    assert picking_level in ["atom","residue"]
    if picking_level != self.picking_level: # only do stuff if needed
      # update internal state
      self.picking_level = picking_level

      # Update widget, emit signal
      if picking_level == "residue":
        self.view.picking_level.setCurrentIndex(1)
        self.state.signals.picking_level.emit("residue")
      elif picking_level == "atom":
        self.view.picking_level.setCurrentIndex(0)
        self.state.signals.picking_level.emit("atom")
    
    assert self.picking_level_displayed == picking_level, f"picking level displayed: {self.picking_level_displayed}, picking level set: {picking_level}"
    


  def start_selecting(self,*args):
    # selection button was clicked. is_checked is NEW state
    if self.view.selector_toggle.is_on: # opposite of intuition, 'on'
      self.viewer.toggle_selection_mode(True)
    else:
      self.viewer.toggle_selection_mode(False)

  def clear_viewer(self):
    self.parent.clear_viewer("Clearing the viewer.")


  def active_model_change(self,ref):
    if ref and isinstance(ref,ModelRef):
      if model_ref.data.filename:
        label = model_ref.data.filename
        self.view.active_model_label.setText(label)

