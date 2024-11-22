from pathlib import Path
import platform
import subprocess
import time
from functools import partial

from PySide2.QtCore import Slot
from PySide2.QtCore import Qt
from PySide2.QtWidgets import (
    QFileDialog,
    QMenu,
    QPushButton
)

import iotbx
from mmtbx.monomer_library.pdb_interpretation import grand_master_phil_str
from qttbx.viewers.gui.controller import Controller
from qttbx.viewers.gui.model.refs import ModelRef

class ViewerControlsBaseController(Controller):
  """
  A controller which manages the minimum base level controls present in the molstar GUI tab
  """
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)

    self.state.signals.active_change.connect(self.active_model_change)
    self.representation_toggles = {'ball-and-stick':False,"cartoon":False}
    for key,action in self.view.button_rep.actions.items():
      action.triggered.connect(partial(self.representation_selected, action))

  @property
  def viewer(self):
    return self.parent

  def clear_viewer(self):
    self.parent.clear_viewer("Clearing the viewer.")

  def representation_selected(self,action):

    if action:
      key = action.text()
      if key.startswith("Ball"):
        key = 'ball-and-stick'
      else:
        key = 'cartoon'
      self.representation_toggles[key] = not self.representation_toggles[key]
      show = self.representation_toggles[key]
      print("Need to change representation. Temporarily disabled...")
      self.parent.select_all()
      self.parent.graphics.add_representation(key)
      

  def active_model_change(self,ref):
    if ref and isinstance(ref,ModelRef):
      if ref.filename:
        label = Path(ref.filename).name
        index = self.view.active_model_label.findText(label)
        if index == -1:
          self.view.active_model_label.addItem(label)
