from PySide2.QtCore import Slot
from pathlib import Path
import platform
import subprocess
import numpy as np
from .controller import Controller
from ..state.ref import SelectionRef
from ..view.widgets.viewer_controls import SearchSelectDialog
from PySide2.QtCore import Slot
from PySide2.QtWidgets import QFileDialog
from PySide2.QtWidgets import QMenu, QPushButton

from PySide2.QtCore import Qt
from ..state.ref import GeometryRef, RestraintsRef, SelectionRef
from ..state.geometry import Geometry
from ..state.restraint import Restraints
from ..state.base import ObjectFrame

from .geometry.tables import GeometryTableTabController
from ..view.widgets import (
  BondEditDialog
)
from ..state.edits import (
  BondEdit
)
from ..state.ref import (
  BondEditsRef

)

class ViewerControlsController(Controller):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)


    ## Signals

    # Picking level
    self.view.picking_level.currentIndexChanged.connect(self.on_picking_change)
    self.state.signals.picking_level.connect(self.on_picking_change)

    # Edits
    self.view.button_edits.setContextMenuPolicy(Qt.CustomContextMenu)
    self.view.button_edits.customContextMenuRequested.connect(self.showContextMenuEdits)
    self.view.button_edits.mousePressEvent = self.buttonMousePressEventEdits  # Override mousePressEvent

    # Restraints
    self.view.button_restraints.clicked.connect(self.load_restraints)

    # Enable return key to execute selection
    self.view.selection_edit.returnPressed.connect(self.execute_selection)
    self.view.button_select.clicked.connect(self.start_selecting)
    self.view.button_clear.clicked.connect(self.clear_selection)

    self.view.add_selection_button.clicked.connect(self.add_selection)
    # # self.view.start_selecting.clicked.connect(self.start_selecting)
    # # self.view.button_deselect.clicked.connect(self.deselect)

    self.state.signals.model_change.connect(self.active_model_change)
    self.view.button_search.clicked.connect(self.search_select_dialog)
    # self.view.combo_box.currentIndexChanged.connect(self.on_picking_change)
    # self.view.button_clear.clicked.connect(self.clear_viewer) # TODO: Move out of selection controls
    # self.state.signals.picking_level.connect(self.on_picking_change)

   # Geometry
    self.view.button_geo.setContextMenuPolicy(Qt.CustomContextMenu)
    self.view.button_geo.customContextMenuRequested.connect(self.showContextMenuGeo)
    self.view.button_geo.mousePressEvent = self.buttonMousePressEventGeo  # Override mousePressEvent


    # Hierarchy Selector
    self.search_select_dialog_controller = None

  



  @property
  def viewer(self):
    return self.parent

  def clear_viewer(self):
    self.parent.clear_viewer("Clearing the viewer.")


  def active_model_change(self,model_ref):
    self.view.active_model_label.clear() #TODO: maintain a list instead
    label = model_ref.data.filename
    self.view.active_model_label.addItem(label)


  def open_file_explorer(self,path):
    if platform.system() == 'Windows':
      subprocess.run(['explorer', path])
    elif platform.system() == 'Darwin':
      subprocess.run(['open', path])
    elif platform.system() == 'Linux':
      subprocess.run(['xdg-open', path])

  @property
  def file_folder(self):

    folder = Path.home()
    if self.state.active_model_ref is not None:
      folder = Path(self.state.active_model_ref.data.filepath).parent

    if folder is None or not folder.exists():
      folder = Path.home()

    return folder

  @Slot()
  def on_picking_change(self,choice):
    current_index = self.view.picking_level.currentIndex()
    if 'residue' in choice:
      if current_index != 0:
        self.view.picking_level.setCurrentIndex(0)
        self.state.signals.picking_level.emit("residue")

    elif 'atom' in choice:
      if current_index != 1:
        self.view.picking_level.setCurrentIndex(1)
        self.state.signals.picking_level.emit("atom")

    


  def buttonMousePressEventEdits(self, event):
    if event.button() == Qt.LeftButton:
        self.showContextMenuEdits(event.pos())
    else:
        QPushButton.mousePressEvent(self.view.button_edits, event)  # Call the original mousePressEvent


  def showContextMenuEdits(self, position):
    # Override the default context menu
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
    contextMenu.exec_(self.view.button_edits.mapToGlobal(position))
    #contextMenu.exec_()

  def load_restraints(self,is_checked):
    # Restraints button clicked, is_checked is the NEW state
    if not self.view.button_restraints.isChecked():
      # clear restraints
      ref = self.state.active_model_ref.restraints_ref
      for sub_ref in ref.children:
        self.state.remove_ref(sub_ref)
      self.state.remove_ref(ref)

    else:
      print("loading, not checked")
      # add restraints
      restraints = Restraints.from_sites(self.state.sites)
      ref = RestraintsRef(data=restraints,model_ref=self.state.active_model_ref,show=True)
      self.state.add_ref(ref) # will emit signals.restraints_change
      # switch to tab
      self.state.signals.tab_change.emit("Restraints")


  def stage_as_bond_restraint(self):

    # Create an 'add' edit for a new restraint
    sites = self.state.mol.sites
    if not self.state.active_selection_ref:
      selection = self.viewer.poll_selection()
      self.state.add_ref_from_selection(selection,make_active=True,show=False)
    ref = self.state.active_selection_ref
    selection = ref.selection
    sel_sites = sites.select_from_selection(selection)
    n_atoms = len(sel_sites)
    assert n_atoms==2, "Bond restraint must only have two atoms"
    i_seqs = list(sel_sites.index)
    xyz = sel_sites.xyz.astype(float)
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
      labels_compositional = GeometryTableTabController.get_labels_compositional_from_iseqs(sites,i_seqs)
      row = BondEdit(
        i_seqs = i_seqs,
        sel_strings = GeometryTableTabController.get_sel_strings_from_iseqs(self.state.mol,i_seqs),
        labels_compositional = labels_compositional,
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
    list(sel.index)
    xyz = sel.xyz.astype(float)
    a,b,c = xyz

  def stage_as_dihedral_restraint(self):
    self.state.signals.stage_restraint.emit(self.ref,"dihedral")

  def stage_as_chiral_restraint(self):
    self.state.signals.stage_restraint.emit(self.ref,"chiral")

  def stage_as_plane_restraint(self):
    geometry = self.ref.model_ref.geometry_ref.data
    sel = self.state.mol.sites.select_from_query(self.ref.query)
    i_seqs = list(sel.index.values)
    new_restraint = PlaneGeometry(i_seqs=i_seqs,
                                   weights=[PlaneGeometry.default_weight for i in range(len(sel))])
    geometry.add_plane_restraint(new_restraint)
    self.state.signals.geometry_change.emit(self.ref.model_ref.geometry_ref)


  @Slot()
  def select_active_selection(self):
    #self.highlight_persist(None,value=True) # set the highlight persist automatically when selecting
    if self.state.active_selection_ref is None:
      self.viewer.deselect_all()
    else:
      self.viewer.select_from_query(self.state.active_selection_ref.query)

  def clear_selection(self):
    self.state.active_selection = None
    self.viewer.deselect_all()

  @Slot()
  def start_selecting(self,is_checked):
    # selection button was clicked. is_checked is NEW state
    if is_checked:
      self.viewer.toggle_selection_mode(True)
    else:
      self.viewer.toggle_selection_mode(False)


  @Slot()
  def execute_selection(self):
    """
    This is a selection from the text box
    """
    text = self.view.selection_edit.text()
    if text.startswith("select"):
      text = text[7:]
    elif text.startswith("sel "):
      text = text[4:]
    
    try:
      self.viewer.select_from_phenix_string(text)
      self.viewer.focus_selected()

      self.save_text_to_history()
    except:
      raise
      self.view.selection_edit.clear()
      self.view.selection_edit.setPlaceholderText(f"Unable to interpret selection: {text}")
  
  def buttonMousePressEventGeo(self, event):
    self.view.button_geo.setChecked(True)
    if event.button() == Qt.LeftButton:
        self.showContextMenuGeo(event.pos())
    else:
        QPushButton.mousePressEvent(self.view.button_geo, event)  # Call the original mousePressEvent

  def showContextMenuGeo(self, position):
    contextMenu = QMenu(self.view)

    # Add actions to the context menu
    action1 = contextMenu.addAction("Generate geometry")
    action2 = contextMenu.addAction("Open geometry file")

    # Connect actions to functions/slots
    action1.triggered.connect(self.process_and_make_geometry)
    action2.triggered.connect(self.load_geometry)

    # Show the context menu at the button's position
    contextMenu.exec_(self.view.button_geo.mapToGlobal(position))
    #self.update_geo_checkbox()

  def process_and_make_geometry(self):
    ref = self.state.active_model_ref
    if ref.has_geometry:
      self.state.notify("Already have geometry loaded. New processing will replace existing geometry.")
    try:
      geodata = Geometry.from_mmtbx_model(ref.model)
      georef = GeometryRef(geodata,ref)
      self.state.add_ref(georef)
    except:
      self.state.notify("Failed to process and make geometry.")
      raise

  def load_geometry(self):
    ref = self.state.active_model_ref
    if ref.has_geometry:
      self.state.notify("Already have geometry loaded. Will now replace existing geometry.")
    self.open_geometry_file_dialog()
  
  def open_geometry_file_dialog(self):

    self.openFileDialog = QFileDialog(self.view)
    self.openFileDialog.setFileMode(QFileDialog.AnyFile)
    if self.openFileDialog.exec_():
        file_path = self.openFileDialog.selectedFiles()[0]
        filepath = str(Path(file_path).absolute())
        print(f"Geometry file selected: {filepath}")
        data = Geometry.from_geo_file(filepath)
        ref = GeometryRef(data,self.state.active_model_ref)
        self.state.add_ref(ref)

  def save_text_to_history(self):
    text = self.view.selection_edit.text()
    if text:
      self.view.selection_edit.history.insert(0, text)
      self.view.selection_edit.index = -1
      self.view.selection_edit.clear()

  def add_selection(self):
    """
    This is when the 'add selection'
    """
    selection = self.viewer.poll_selection()
    
    sel = self.state.active_mol.select_from_selection(selection)
    # query = SelectionQuery.from_json(selection_query_json)
    # #assert len(query_dict)==1, "Multi structure queries not yet supported"
    # ref_id = query.params.refId
    # if ref_id not in self.state.references:
    #   ref_id = self.state.active_model_ref.id
    # ref = self.state.references[ref_id]
    # query.params.keymap = self.state.mmcif_column_map
    # sel_sites = self.state.active_mol.sites.select_from_query(query)
    if len(sel)>0:
      sel_ref = SelectionRef(selection,model_ref=self.state.active_model_ref,show=True)
      self.state.add_ref(sel_ref)
      self.state.active_selection_ref = sel_ref
      self.state.signals.tab_change.emit("Selections") # show selection tab
    else:
      print("Skipping add selection due to empty selection")


  def search_select_dialog(self):
    dialog = SearchSelectDialog(title="Selection search", parent=self.view)
    self.view.search_select_dialog = dialog
    self.search_select_dialog_controller = SearchSelectDialogController(parent=self,view=self.view.search_select_dialog)
    default_width = dialog.width()
    new_width = default_width * 6
    dialog.setMinimumWidth(new_width)
    dialog.exec_()

class SearchSelectDialogController(Controller):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)
    mol = self.state.active_model_ref.mol
    chain_options = sorted([str(e) for e in list(mol.sites["asym_id"].unique())])
    self.chain_scroller_controller = ScrollableHierarchyController(parent=self,view=self.view.chain_scroller)
    self.chain_scroller_controller.update_list_widget(chain_options)
    self.chain_scroller_controller.view.selected_change.connect(self.chain_change)



    residue_options = sorted([e for e in list(mol.sites["comp_id"].unique())])
    self.comp_scroller_controller = ScrollableHierarchyController(parent=self,view=self.view.comp_scroller)
    self.comp_scroller_controller.update_list_widget(residue_options)
    self.comp_scroller_controller.view.selected_change.connect(self.comp_change)

    residue_options = sorted([int(e) for e in list(mol.sites["seq_id"].unique())])
    residue_options = [str(e) for e in residue_options]
    self.seq_scroller_controller = ScrollableHierarchyController(parent=self,view=self.view.seq_scroller)
    self.seq_scroller_controller.update_list_widget(residue_options)
    self.seq_scroller_controller.view.selected_change.connect(self.seq_change)


    atom_options = sorted([str(e) for e in list(mol.sites["atom_id"].unique())])
    self.atom_scroller_controller = ScrollableHierarchyController(parent=self,view=self.view.atom_scroller)
    self.atom_scroller_controller.update_list_widget(atom_options)
    self.atom_scroller_controller.view.selected_change.connect(self.atom_change)

  def update_options(self,current,controller, mol, filters, group_by, sort_as_int=False):
    sites = mol.sites
    sites_sel = sites
    for key, value in filters.items():

      sites_sel = sites_sel[sites_sel[key] == value]

    options = sorted([str(e) for e in list(sites_sel[group_by].unique())])
    if sort_as_int:
      options = sorted([int(e) for e in options])
      options = [str(e) for e in options]

    controller.update_list_widget(options)
    if current:
      controller.set_selected_item(current)
    return options

  def is_specific(self,value):
    if value not in ["all",None]:
      return True
    else:
      return False

  def chain_change(self, selected):
    self.comp_scroller_controller.reset()
    self.seq_scroller_controller.reset()
    self.atom_scroller_controller.reset()

    mol = self.state.active_model_ref.mol
    self.chain_scroller_controller.selected
    comp = self.comp_scroller_controller.selected
    seq = self.seq_scroller_controller.selected 
    atom = self.atom_scroller_controller.selected
    # Update residues and atoms
    filters = {}
    sel_str_parts = []
    if self.is_specific(selected):
      sel_str_parts.append(f"chain {selected}")
      filters["asym_id"] = selected
    if self.is_specific(seq):
      sel_str_parts.append(f"resseq {seq}")
      filters["seq_id"] = seq
    if self.is_specific(comp):
      sel_str_parts.append(f"resname {comp}")
      filters["comp_id"] = comp
    if self.is_specific(atom): 
      sel_str_parts.append(f"name {atom}")
      filters["atom_id"] = atom
    if len(sel_str_parts)==0:
      sel_str = "all"
    else:
      sel_str = " and ".join(sel_str_parts)
    print("Sel str: ",sel_str)
    comp_options = self.update_options(comp,self.comp_scroller_controller, mol, filters, "comp_id")
    seq_options = self.update_options(seq,self.seq_scroller_controller, mol, filters, "seq_id", sort_as_int=True)
    atom_options = self.update_options(atom,self.atom_scroller_controller, mol,filters, "atom_id")
    all_options = comp_options +  seq_options + atom_options
    if len(all_options)==0:
      self.reset()
    else:
      self.parent.viewer.select_from_phenix_string(sel_str)
      self.parent.viewer.focus_selected()

  def comp_change(self, selected):
    self.seq_scroller_controller.reset()
    self.atom_scroller_controller.reset()
    mol = self.state.active_model_ref.mol
    chain = self.chain_scroller_controller.selected
    seq = self.seq_scroller_controller.selected 
    atom = self.atom_scroller_controller.selected
    comp = selected
    filters = {}
    
    sel_str_parts = []
    if self.is_specific(comp):
      sel_str_parts.append(f"resname {comp}")
      filters["comp_id"] = comp
    if self.is_specific(chain):
      sel_str_parts.append(f"chain {chain}")
      filters["asym_id"] = chain
    if self.is_specific(seq):
      sel_str_parts.append(f"resseq {seq}")
      filters["seq_id"] = seq
    if self.is_specific(atom): 
      sel_str_parts.append(f"name {atom}")
      filters["atom_id"] = atom

    if len(sel_str_parts)==0:
      sel_str = "all"
    else:
      sel_str = " and ".join(sel_str_parts)
    print("Sel str: ",sel_str)
    seq_options = self.update_options(seq,self.seq_scroller_controller, mol, filters, "seq_id", sort_as_int=True)
    atom_options = self.update_options(atom,self.atom_scroller_controller, mol, filters, "atom_id")
    all_options = seq_options + atom_options
    if len(all_options)==0:
      self.reset()
    else:
      self.parent.viewer.select_from_phenix_string(sel_str)
      self.parent.viewer.focus_selected()

  def seq_change(self, selected):
    mol = self.state.active_model_ref.mol
    chain = self.chain_scroller_controller.selected
    comp = self.comp_scroller_controller.selected
    atom = self.atom_scroller_controller.selected
    
    if selected != "all":
      selected = int(selected)

    seq = selected
    filters = {}
      
    #comp_options = self.update_options(self.comp_scroller_controller, mol, {"asym_id": chain, "seq_id":selected}, "comp_id", sort_as_int=False)
    sel_str_parts = []
    if self.is_specific(seq):
      sel_str_parts.append(f"resseq {seq}")
      filters["seq_id"] = seq
    if self.is_specific(chain):
      sel_str_parts.append(f"chain {chain}")
      filters["asym_id"] = chain
    if self.is_specific(comp):
      sel_str_parts.append(f"resname {comp}")
      filters["comp_id"] = comp

    if len(sel_str_parts)==0:
      sel_str = "all"
    else:
      sel_str = " and ".join(sel_str_parts)
    print("Sel str: ",sel_str)
    atom_options = self.update_options(atom,self.atom_scroller_controller, mol, filters, "atom_id")
    all_options = atom_options
    if len(all_options)==0:
      self.reset()
    else:
      self.parent.viewer.select_from_phenix_string(sel_str)
      self.parent.viewer.focus_selected()

  def atom_change(self,selected):
    self.parent.viewer.set_granularity(value="element")
    self.state.active_model_ref.mol
    chain = self.chain_scroller_controller.selected
    seq = self.seq_scroller_controller.selected
    comp = self.comp_scroller_controller.selected

    sel_str_parts = []
    if selected != "all":
      sel_str_parts.append(f"name {selected}")
    if self.is_specific(chain):
      sel_str_parts.append(f"chain {chain}")
    if self.is_specific(seq):
      sel_str_parts.append(f"resseq {seq}")
    if self.is_specific(comp):
      sel_str_parts.append(f"resname {comp}")

    if len(sel_str_parts)==0:
      sel_str = "all"
    else:
      sel_str = " and ".join(sel_str_parts)
    print("Sel str: ",sel_str)
    self.parent.viewer.select_from_phenix_string(sel_str)
    print("emitting atom")
    self.state.signals.picking_level.emit("atom")
    self.parent.viewer.focus_selected()

  def reset(self):
    self.chain_scroller_controller.reset()
    self.comp_scroller_controller.reset()
    self.seq_scroller_controller.reset()
    self.atom_scroller_controller.reset()
    self.parent.viewer.deselect_all()

class ScrollableHierarchyController(Controller):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)
    # Connect the list widget item selection event to a handler function
    self.view.list_widget.itemSelectionChanged.connect(self.on_item_selection_changed)
    
    self._selected = None
    self.last_selected_time = -1
    self.view.list_widget.keyPressFilter.space_pressed.connect(self.on_space_pressed)
    self._initial_list = None

  @property
  def selected(self):
    return self._selected

  @selected.setter
  def selected(self,value):
    self._selected = value

  def on_space_pressed(self):
    self.advance_one()

  def reset(self):
    self.selected = None
    self.update_list_widget(self._initial_list)

  def advance_one(self):
    # Get the current selected items
    selected_items = self.view.list_widget.selectedItems()
    if not selected_items:
      return  # No item is selected

    # Get the index of the current selected item
    current_index = self.view.list_widget.row(selected_items[0])
    # Calculate the next index
    next_index = current_index + 1

    # Check if the next index is within the valid range
    if next_index < self.view.list_widget.count():
      # Select the next item
      self.view.list_widget.setCurrentRow(next_index)
      # Emit the signal with the new selected value
      new_value = self.view.list_widget.item(next_index).text()
      self.view.selected_change.emit(new_value)


  def update_list_widget(self, items):
    items = list(items)
    if "all" not in items:
      items = ["all"]+ items
    if self._initial_list is None:
      self._initial_list = items
    self.view.update_list_widget(items)

  def set_selected_item(self, text):
    items = self.view.list_widget.findItems(text, Qt.MatchExactly)
    if items:
      self.view.list_widget.setCurrentItem(items[0])

  def on_item_selection_changed(self):
    # Get the selected items
    selected_items = self.view.list_widget.selectedItems()

    # Update the label text based on the selection
    if selected_items:
      selected_texts = [item.text() for item in selected_items]
      selected = ', '.join(selected_texts)
      self.view.label.setText(f"Selected: {selected}")
      self.selected = selected
      self.view.selected_change.emit(selected)
    else:
      self.view.label.setText("Selected: None")
      self.selected = None

    