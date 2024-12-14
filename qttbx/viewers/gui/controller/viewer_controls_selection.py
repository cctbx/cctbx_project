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

from qttbx.viewers.gui.controller.viewer_controls_base import ViewerControlsBaseController

class ViewerControlsSelectionController(ViewerControlsBaseController):
  """
  A controller which manages the selection controls
  """
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)


    self.picking_level = 'atom' # internal state
    self.state.signals.picking_level.connect(self.set_picking_level)
    self.view.picking_level.currentTextChanged.connect(self.picking_level_change) # user initiates new level
    self.state.signals.picking_level.connect(self.picking_level_change) # update view to reflect new level
    self.view.selector_toggle.clicked.connect(self.start_selecting)
    self.view.button_cancel.clicked.connect(self.viewer.select_none)
    self.view.add_selection_button.clicked.connect(self.add_selection)
    
  @property
  def graphics(self):
    return self.parent.graphics

  def set_picking_level(self,picking_level):
    if 'atom' in picking_level:
      self.set_granularity("element")
    elif "residue" in picking_level:
      self.set_granularity("residue")
    else:
      pass

  def set_granularity(self,granularity="residue"):
    assert granularity in ['element','residue'], 'Provide one of the implemented picking levels'
    self._picking_granularity = granularity
    self.graphics.set_granularity(granularity=granularity)


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
    
  def toggle_selection_mode(self,value):
    if value == True:
      self.graphics.selection_mode_on()
    else:
      self.graphics.selection_mode_off()



  def start_selecting(self,*args):
    if self.view.selector_toggle.is_on: # opposite of intuition, 'on'
      self.toggle_selection_mode(True)
    else:
      self.toggle_selection_mode(False)


  def add_selection(self):
    synced = self.viewer.sync_selection()
    ref = self.state.active_selection_ref
    ref.show = True
    self.state.signals.selection_added.emit(ref) # Update tables
    if synced:
      self.state.signals.tab_change.emit("Selections") # show selection tab
