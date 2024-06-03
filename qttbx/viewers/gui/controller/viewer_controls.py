from PySide2.QtCore import Slot
import json
import time
from .controller import Controller
from ...core.selection_utils import SelectionQuery
from ..state.ref import SelectionRef
from ..view.widgets.viewer_controls import SearchSelectDialog
from PySide2.QtCore import QObject, QTimer, Signal, Slot
from PySide2.QtWidgets import QApplication
from PySide2.QtWidgets import QFileDialog, QColorDialog
from PySide2.QtWidgets import QApplication, QPushButton, QMenu, QMainWindow, QVBoxLayout, QWidget

from PySide2.QtCore import Qt
from .scroll_list import ScrollableListController
from ..state.ref import SelectionRef, ModelRef, GeometryRef
from ..state.geometry import Geometry
from .widgets import InputDialog

class ViewerControlsController(Controller):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)
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
    self.view.button_geo.customContextMenuRequested.connect(self.showContextMenu)
    self.view.button_geo.mousePressEvent = self.buttonMousePressEvent  # Override mousePressEvent
    

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
    self.viewer.toggle_selection_mode(False)

  @Slot()
  def start_selecting(self):
    # selection button was clicked
    self.viewer.toggle_selection_mode(True)
    #self.execute_selection()


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
  
  def buttonMousePressEvent(self, event):
    if event.button() == Qt.LeftButton:
        self.showContextMenu(event.pos())
    else:
        QPushButton.mousePressEvent(self.view.button_geo, event)  # Call the original mousePressEvent

  def showContextMenu(self, position):
    contextMenu = QMenu(self.view)

    # Add actions to the context menu
    action1 = contextMenu.addAction("Generate geometry")
    action2 = contextMenu.addAction("Open geometry file")

    # Connect actions to functions/slots
    action1.triggered.connect(self.process_and_make_restraints)
    action2.triggered.connect(self.load_restraints)

    # Show the context menu at the button's position
    contextMenu.exec_(self.view.button_geo.mapToGlobal(position))

  def process_and_make_restraints(self):
      ref = self.state.active_model_ref
      if ref.has_restraints:
        self.state.notify("Already have restraints loaded. New processing will replace existing restraints.")
      try:
        restraints = Geometry.from_model_ref(ref)
        georef = GeometryRef(restraints,ref)
        self.state.add_ref(georef)
      except:
        self.state.notify("Failed to process and make restraints.")
        raise

  def load_restraints(self):
      if self.ref.has_restraints:
        self.state.notify("Already have restraints loaded. Will now replace existing restraints.")
      self.open_restraints_file_dialog()

  def open_restraints_file_dialog(self):

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

  @Slot()
  def on_picking_change(self,index):
    if index == 0:
      self.viewer.set_granularity(value='residue')
    elif index == 1:
      self.viewer.set_granularity(value='element')

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
    chain = self.chain_scroller_controller.selected
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
    mol = self.state.active_model_ref.mol
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

    