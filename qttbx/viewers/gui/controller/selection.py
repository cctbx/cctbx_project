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

class SelectionEntryController(ModelLikeEntryController):
  def __init__(self,parent=None,view=None,ref=None):
    super().__init__(parent=parent,view=view,ref=ref)
    self.state.signals.selection_activated.connect(self.handle_selection_change)
    self.view.button_info.clicked.connect(self.display_info)


  def handle_selection_change(self,new_selection_ref):
    # Just disable the toggle if self is not active
    if (new_selection_ref is  None  or
     new_selection_ref.uuid != self.ref.uuid):
      self.view.active_toggle.is_checked = False

  def toggle_active_func(self,is_checked):
    if is_checked:
      self.state.active_selection_ref = self.ref

    else:
      if self.state.active_selection_ref == self.ref:
        self.state.active_selection_ref = None


  def display_info(self):
    text = f"""
    Model Reference label: {self.ref.model_ref.label}
    Number of atoms: {self.ref.model.get_number_of_atoms()}

    Phenix string: {self.ref.selection.phenix_string}
    """
    dialog = InfoDialog(text, title="Selection Info", parent=self.view)
    default_width = dialog.width()
    new_width = default_width * 6
    dialog.setMinimumWidth(new_width)

    dialog.exec_()

class SelectionListController(ScrollableListController):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)
    self.next_selection_number = 1 # for labeling new selections


    self.state.signals.selection_added.connect(self.update)


  def add_entry_from_ref(self,ref: SelectionRef,force_show=True):
    if ref not in self.refs:
      if not force_show and not ref.show:
        return 
      entry_view = SelectionEntryView()
      entry_controller = SelectionEntryController(parent=self,view=entry_view,ref=ref)
      is_checked = True
      label = ref.selection.phenix_string
      if label.strip() == "all":
        label = f"(all) {ref.model_ref.label}"
        is_checked = False # don't select full model at first
      elif len(label)<30:
        pass # keep phenix selection

      elif len(label)< 100:
        label = label[:30]+"..." # truncate
      else:
        # convert to number of atoms
        label = f"{ref.number_of_atoms} atoms selected"
      ref.label = label
      entry_controller.view.active_toggle.is_checked = is_checked
      entry_controller.view.label_name.setText(ref.label)
      self.add_entry(entry_controller)
      self.next_selection_number+=1


  def update(self):
    for ref in self.state.references_selection:
      self.add_entry_from_ref(ref,force_show=False)



class SelectionTabController(Controller):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)

    self.selection_list_controller = SelectionListController(parent=self,view=self.view.selection_list_view)

