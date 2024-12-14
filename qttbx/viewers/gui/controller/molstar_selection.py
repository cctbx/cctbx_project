import math

import numpy as np
from PySide2.QtWidgets import QMenu

from ..model.refs import SelectionRef
from .models import ModelLikeEntryController
from ..view.selection import SelectionEntryView
from ..view.widgets.dialogs import InfoDialog
from .scroll_list import ScrollableListController
from .controller import Controller
from .viewer_controls_selection import ViewerControlsSelectionController
from .molstar_base import MolstarBaseController
from qttbx.viewers.selection import Selection


class MolstarSelectionController(MolstarBaseController):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)

    # Selections
    self.state.signals.active_change.connect(self.change_selection)
    self.state.signals.select_all.connect(self.select_all)
    self.state.signals.select_none.connect(self.select_none)

    self.graphics_controls = ViewerControlsSelectionController(parent=self,view=self.view.viewer_controls)

  def change_selection(self,selection_ref):
    if isinstance(selection_ref,SelectionRef):
      current_selection = self.poll_selection()
      if current_selection != selection_ref.selection:
        print(selection_ref.selection.string)
        if len(selection_ref.selection)>0:
          self.graphics.select(selection_ref.selection.string)

      else:
        print("Skipping selection because same as polled")
        
    

  def sync_selection(self):
    """
    Poll selection from molstar viewer and set state. 
    Returns True if it produced a valid non-empty selection
    """
    selection = self.poll_selection()
    if selection.is_empty:
      self.log("Skipping add selection due to empty selection")
      self.state.active_selection_ref = None
      #self.state.active_selection = None
      return False
    else:
      sel_ref = SelectionRef(selection,model=self.state.active_model_ref,show=False)
      self.state.add_ref(sel_ref)
      self.state.active_selection_ref = sel_ref
      #self.state.active_selection = selection
      return True

  def poll_selection(self):
    """
    Poll selection from molstar viewer
    Returns selection object
    """
    atom_records = self.graphics.poll_selection()
    selection = Selection.from_atom_records(atom_records,self.state.active_model)
    return selection