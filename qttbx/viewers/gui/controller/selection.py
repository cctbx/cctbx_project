from PySide2.QtWidgets import QApplication, QPushButton, QMenu, QMainWindow, QVBoxLayout, QWidget


from .models import ModelLikeEntryController
from ..view.tabs.selection import SelectionEntryView
from ..view.widgets import InfoDialog
from .scroll_list import ScrollableListController
from .controller import Controller
from ..state.restraints import (
  BondRestraint,
  AngleRestraint,
  DihedralRestraint,
  ChiralRestraint,
  PlaneRestraint
)

class SelectionEntryController(ModelLikeEntryController):
  def __init__(self,parent=None,view=None,ref=None):
    super().__init__(parent=parent,view=view,ref=ref)
    self.state.signals.selection_change.connect(self.handle_selection_change)
    self.view.button_info.clicked.connect(self.display_info)


  def handle_selection_change(self):
    # Just disable the toggle if self is not active
    if (self.state.active_selection_ref is  None  or
     self.state.active_selection_ref.id != self.ref.id):
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
    action1 = contextMenu.addAction("Stage as Bond restraint")
    action2 = contextMenu.addAction("Stage as Angle restraint")
    action3 = contextMenu.addAction("Stage as Dihedral restraint")
    action4 = contextMenu.addAction("Stage as Chiral restraint")
    action5 = contextMenu.addAction("Stage as Plane restraint")

    # Connect actions to functions/slots
    action1.triggered.connect(self.stage_as_bond_restraint)
    action2.triggered.connect(self.stage_as_angle_restraint)
    action3.triggered.connect(self.stage_as_dihedral_restraint)
    action4.triggered.connect(self.stage_as_chiral_restraint)
    action5.triggered.connect(self.stage_as_plane_restraint)

    # Show the context menu at the button's position
    contextMenu.exec_(self.view.button_restraints.mapToGlobal(position))

  def stage_as_bond_restraint(self):
    restraints = self.ref.model_ref.restraints_ref.data
    sel = self.state.mol.sites.select_from_query(self.ref.query)
    i_seqs = list(sel.index.values)
    new_restraint = BondRestraint(i_seqs=i_seqs)
    restraints.add_bond_restraint(new_restraint)
    self.state.signals.restraints_change.emit(self.ref.model_ref.restraints_ref)

  def stage_as_angle_restraint(self):
    self.state.signals.stage_restraint.emit(self.ref,"angle")

  def stage_as_dihedral_restraint(self):
    self.state.signals.stage_restraint.emit(self.ref,"dihedral")

  def stage_as_chiral_restraint(self):
    self.state.signals.stage_restraint.emit(self.ref,"chiral")

  def stage_as_plane_restraint(self):
    restraints = self.ref.model_ref.restraints_ref.data
    sel = self.state.mol.sites.select_from_query(self.ref.query)
    i_seqs = list(sel.index.values)
    new_restraint = PlaneRestraint(i_seqs=i_seqs,
                                   weights=[PlaneRestraint.default_weight for i in range(len(sel))])
    restraints.add_plane_restraint(new_restraint)
    self.state.signals.restraints_change.emit(self.ref.model_ref.restraints_ref)

  def display_info(self):
    # TODO: this is a view, should move to the view directory
    text = f"""
    Reference id: {self.ref.id}
    Model Reference id: {self.ref.model_ref.id}
    External ids:
    {self.ref.model_ref.external_ids}

    Phenix string: {self.ref.query.phenix_string}
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

    self.state.signals.selection_change.connect(self.update)


  def update(self):
    selection_list = self
    for ref in self.state.references_selection:
      if ref not in selection_list.refs:
        if ref.show_in_list:
          entry_view = SelectionEntryView()
          entry_controller = SelectionEntryController(parent=self,view=entry_view,ref=ref)
          entry_controller.view.active_toggle.is_checked = True
          selection_list.add_entry(entry_controller)
          ref.label = f"Selection {selection_list.next_selection_number}"
          entry_controller.view.label_name.setText(ref.label)
          selection_list.next_selection_number+=1
          # make new selection active



class SelectionTabController(Controller):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)

    self.selection_list_controller = SelectionListController(parent=self,view=self.view.selection_list_view)
