from PySide2.QtWidgets import QApplication, QPushButton, QMenu, QMainWindow, QVBoxLayout, QWidget
import numpy as np
import math

from ..state.ref import SelectionRef
from .models import ModelLikeEntryController
from ..view.tabs.selection import SelectionEntryView
from ..view.widgets import InfoDialog
from .scroll_list import ScrollableListController
from .controller import Controller
from ..state.base import ObjectFrame
from ..state.geometry import (
  BondGeometry,
  AngleGeometry,
  DihedralGeometry,
  ChiralGeometry,
  PlaneGeometry
)
from ..state.edits import (
  BondEdit,
  AngleEdit,
  DihedralEdit
)

from ..view.widgets import (
  BondEditDialog,
  AngleEditDialog,
  DihedralEditDialog,
  ChiralEditDialog,
  PlaneEditDialog
)
from ..state.ref import (
  BondEditsRef,
  AngleEditsRef,
  DihedralEditsRef

)

class SelectionEntryController(ModelLikeEntryController):
  def __init__(self,parent=None,view=None,ref=None):
    super().__init__(parent=parent,view=view,ref=ref)
    self.state.signals.selection_activated.connect(self.handle_selection_change)
    self.view.button_info.clicked.connect(self.display_info)


  def handle_selection_change(self,new_selection_ref):
    # Just disable the toggle if self is not active
    if (new_selection_ref is  None  or
     new_selection_ref.id != self.ref.id):
      self.view.active_toggle.is_checked = False

  def toggle_active_func(self,is_checked):
    # TODO: Move this to data tab?
    if is_checked:
      self.state.active_selection_ref = self.ref

    else:
      #print("The entry is unchecked.")
      if self.state.active_selection_ref == self.ref:
        self.state.active_selection_ref = None


  def showContextMenu(self, position):
    # Override the model-like-entry context menu
    contextMenu = QMenu(self.view)

    # Add actions to the context menu
    action1 = contextMenu.addAction("Edit as Bond restraint")
    action2 = contextMenu.addAction("Edit as Angle restraint")
    action3 = contextMenu.addAction("Edit as Dihedral restraint")
    action4 = contextMenu.addAction("Edit as Chiral restraint")
    action5 = contextMenu.addAction("Edit as Plane restraint")

    # Connect actions to functions/slots
    action1.triggered.connect(self.stage_as_bond_restraint)
    action2.triggered.connect(self.stage_as_angle_restraint)
    action3.triggered.connect(self.stage_as_dihedral_restraint)
    action4.triggered.connect(self.stage_as_chiral_restraint)
    action5.triggered.connect(self.stage_as_plane_restraint)

    # Show the context menu at the button's position
    contextMenu.exec_(self.view.button_restraints.mapToGlobal(position))

  def stage_as_bond_restraint(self):
    # NOTE: Commented out old way of modifying restraint data object
    # restraints = self.ref.model_ref.geometry_ref.data
    # sel = self.state.mol.sites.select_from_query(self.ref.query)
    # i_seqs = list(sel.index.values)
    # new_restraint = BondGeometry(i_seqs=i_seqs)
    # restraints.add_bond_restraint(new_restraint)
    # self.state.signals.geometry_change.emit(self.ref.model_ref.geometry_ref)

    # New way is to create an edit
    mol = self.state.active_model_ref.mol
    sel = mol.sites.select_from_query(self.ref.query)
    assert len(sel)==2, "Bond restraint must only have two atoms"
    i_seqs = list(sel.index)
    xyz = sel.xyz.astype(float)
    d = round(np.linalg.norm(xyz[1]-xyz[0]),3)
    defaults_dict = {
      "action":"add",
      "ideal":d,
      }

    other_defaults = BondEdit.defaults(BondEdit)
    for k,v in other_defaults.items():
      if k not in defaults_dict:
        if k.endswith("_new"):
          k = k.replace("_new","")
          defaults_dict[k] = v

    dialog = BondEditDialog(defaults_dict=defaults_dict)
    if dialog.exec_():
      print("Making edit")
      row = BondEdit(
        i_seqs = i_seqs,
        sel_strings = [mol.sites.select_query_from_i_seqs([i]).phenix_string for i in i_seqs],
        ideal_old=d,
        ideal_new = dialog.collectInputValues()["ideal"],
        sigma_old=None,
        sigma_new = dialog.collectInputValues()["sigma"],
        action="add"
        )

      edit_ref = None
      for ref_id,ref in self.state.references.items():
        if isinstance(ref,BondEditsRef):
          if ref.geometry_ref == self.state.active_model_ref.geometry_ref:
            edit_ref = ref
            edit_ref.data.rows.append(row)
      if edit_ref is None:
        objframe = ObjectFrame.from_rows([row])
        edit_ref = BondEditsRef(data=objframe,geometry_ref=self.state.active_model_ref.geometry_ref)
      self.state.add_ref(edit_ref)
    else:
      print("Dialog Cancelled")

  def stage_as_angle_restraint(self):
    #self.state.signals.stage_restraint.emit(self.ref,"angle")
    mol = self.state.active_model_ref.mol
    sel = mol.sites.select_from_query(self.ref.query)
    assert len(sel)==3, "Angle restraint must only have two atoms"
    i_seqs = list(sel.index)
    xyz = sel.xyz.astype(float)
    a,b,c = xyz

    def calculate_angle(A, B, C):
      # TODO: use centralized code
      # Calculate the lengths of the sides
      AB = math.sqrt((B[0] - A[0])**2 + (B[1] - A[1])**2)
      BC = math.sqrt((C[0] - B[0])**2 + (C[1] - B[1])**2)
      CA = math.sqrt((C[0] - A[0])**2 + (C[1] - A[1])**2)
      # Calculate the angle using the Law of Cosines
      cos_theta = (AB**2 + BC**2 - CA**2) / (2 * AB * BC)
      theta = math.acos(cos_theta)  # This gives the angle in radians
      
      # Convert to degrees, if desired
      angle_in_degrees = math.degrees(theta)
      
      return angle_in_degrees

    angle = round(calculate_angle(a,b,c),2)

    defaults_dict = {
      "action":"add",
      "ideal":angle,
      }

    other_defaults = AngleEdit.defaults(AngleEdit)
    for k,v in other_defaults.items():
      if k not in defaults_dict:
        if k.endswith("_new"):
          k = k.replace("_new","")
          defaults_dict[k] = v

    dialog = AngleEditDialog(defaults_dict=defaults_dict)
    if dialog.exec_():
      row = AngleEdit(
        i_seqs = i_seqs,
        sel_strings = [mol.sites.select_query_from_i_seqs([i]).phenix_string for i in i_seqs],
        ideal_old=angle,
        ideal_new = dialog.collectInputValues()["ideal"],
        sigma_old=None,
        sigma_new = dialog.collectInputValues()["sigma"],
        action="add"
        )

      edit_ref = None
      for ref_id,ref in self.state.references.items():
        if isinstance(ref,AngleEditsRef):
          #if ref.geometry_ref == self.state.active_model_ref.geometry_ref:
          edit_ref = ref
          edit_ref.data.rows.append(row)
      if edit_ref is None:
        objframe = ObjectFrame.from_rows([row])
        edit_ref = AngleEditsRef(data=objframe,geometry_ref=self.state.active_model_ref.geometry_ref)
      self.state.add_ref(edit_ref)
    else:
      print("Dialog Cancelled")

  def stage_as_dihedral_restraint(self):
    self.state.signals.stage_restraint.emit(self.ref,"dihedral")

  def stage_as_chiral_restraint(self):
    self.state.signals.stage_restraint.emit(self.ref,"chiral")

  def stage_as_plane_restraint(self):
    restraints = self.ref.model_ref.geometry_ref.data
    sel = self.state.mol.sites.select_from_query(self.ref.query)
    i_seqs = list(sel.index.values)
    new_restraint = PlaneGeometry(i_seqs=i_seqs,
                                   weights=[PlaneGeometry.default_weight for i in range(len(sel))])
    restraints.add_plane_restraint(new_restraint)
    self.state.signals.geometry_change.emit(self.ref.model_ref.geometry_ref)

  def display_info(self):
    text = f"""
    Model Reference label: {self.ref.model_ref.label}
    Number of atoms: {self.ref.number_of_atoms}

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

  def set_focus_on(self):
    # only set focus if populated
    if len(self.selection_list_controller.entries)>0:
      self.view.set_focus_on()