import os
from pathlib import Path

from PySide2.QtWidgets import QFileDialog, QColorDialog

from ..models import ModelEntryController, ModelListController
from ..maps import MapEntryController, MapListController
from ...view.restraints.restraints_staging import RestraintStagingEntryView
from ...view.cif.cif_files import CifFileEntryView
from ...view.maps import MapEntryView, MapListView
from ..scroll_entry import ScrollEntryController
from ..selection import SelectionEntryController, SelectionListController
from ..scroll_list import ScrollableListController
from ..controller import Controller
from ...state.cif import CifFileData
from ...state.ref import CifFileRef

class RestraintStagingEntryController(SelectionEntryController):
  name = "Restraint"
  n_atoms = None

  @classmethod
  def init_and_check_natoms(cls,parent=None,view=None,ref=None):
    self = cls(parent=parent,view=view,ref=ref)
    sel_atoms = len(self.state.mol.sites.select_from_query(ref.query))
    if self.n_atoms is not None:
      if sel_atoms != self.n_atoms:
        self.state.notify(f"Cannot form restraint of type {self.name} from number of atoms: {sel_atoms}")
        return None
    return self

  def __init__(self,parent=None,view=None,ref=None):
    super().__init__(parent=parent,view=view,ref=ref)



class BondStagingEntryController(RestraintStagingEntryController):
  name = "Bond"
  counter = 1
  n_atoms = 2
  def __init__(self,parent=None,view=None,ref=None):
    super().__init__(parent=parent,view=view,ref=ref)

class AngleStagingEntryController(RestraintStagingEntryController):
  name = "Angle"
  counter = 1
  n_atoms = 3
  def __init__(self,parent=None,view=None,ref=None):
    super().__init__(parent=parent,view=view,ref=ref)

class DihedralStagingEntryController(RestraintStagingEntryController):
  name = "Dihedral"
  counter = 1
  n_atoms = 4
  def __init__(self,parent=None,view=None,ref=None):
    super().__init__(parent=parent,view=view,ref=ref)

class ChiralStagingEntryController(RestraintStagingEntryController):
  name = "Chiral"
  counter = 1
  n_atmos = 4
  def __init__(self,parent=None,view=None,ref=None):
    super().__init__(parent=parent,view=view,ref=ref)

class PlaneStagingEntryController(RestraintStagingEntryController):
  name = "Plane"
  counter = 1
  n_atms = None
  def __init__(self,parent=None,view=None,ref=None):
    super().__init__(parent=parent,view=view,ref=ref)

# End entries



class RestraintStagingListController(SelectionListController):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)

    self.state.signals.stage_restraint.connect(self.add_selection_ref)

    self.entry_controller_classes = {
      "bond":BondStagingEntryController,
      "angle":AngleStagingEntryController,
      "dihedral":DihedralStagingEntryController,
      "chiral":ChiralStagingEntryController,
      "plane":PlaneStagingEntryController,
    }

  def update(self):
    return None

  def add_selection_ref(self,ref,type):

    if ref not in self.refs:
      if ref.show_in_list:

        entry_view = RestraintStagingEntryView()
        controller_class = self.entry_controller_classes[type]
        entry_controller = controller_class.init_and_check_natoms(parent=self,view=entry_view,ref=ref)
        if entry_controller:
          entry_controller.view.active_toggle.is_checked = True
          self.add_entry(entry_controller)
          ref.label = f"{controller_class.name} {controller_class.counter}"
          entry_controller.view.label_name.setText(ref.label)
          controller_class.counter+=1
          # switch tab
          self.state.signals.tab_change.emit("Restraints") # show selection tab


