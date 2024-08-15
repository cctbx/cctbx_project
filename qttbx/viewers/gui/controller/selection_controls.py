"""
The controller for managing the selection controls
"""
from itertools import chain

from PySide2.QtCore import Slot
from PySide2.QtCore import Qt

from qttbx.viewers.gui.controller import Controller

class SearchSelectDialogController(Controller):
  """
  A controller which manages the hierarchical selection dialog window
  """
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)
    model = self.state.active_model
    chain_options, resname_options, resseq_options, atom_options = self.update_options("all")
    self.chain_scroller_controller = ScrollableHierarchyController(parent=self,view=self.view.chain_scroller)
    self.chain_scroller_controller.update_list_widget(chain_options)
    self.chain_scroller_controller.view.selected_change.connect(self.chain_change)

    self.comp_scroller_controller = ScrollableHierarchyController(parent=self,view=self.view.comp_scroller)
    self.comp_scroller_controller.update_list_widget(resname_options)
    self.comp_scroller_controller.view.selected_change.connect(self.comp_change)

    self.seq_scroller_controller = ScrollableHierarchyController(parent=self,view=self.view.seq_scroller)
    self.seq_scroller_controller.update_list_widget(resseq_options)
    self.seq_scroller_controller.view.selected_change.connect(self.seq_change)

    self.atom_scroller_controller = ScrollableHierarchyController(parent=self,view=self.view.atom_scroller)
    self.atom_scroller_controller.update_list_widget(atom_options)
    self.atom_scroller_controller.view.selected_change.connect(self.atom_change)

  def update_options(self,selection_string):
    model = self.state.active_model
    model = model.select(model.selection(selection_string))
    
    chain_options = sorted(list(set(model.chain_ids())))
    resname_options = sorted({
      resname.strip()
      for chain in model.chains()
      for residue_group in chain.residue_groups()
      for resname in residue_group.unique_resnames()
    })
    resseq_options = sorted({
      residue_group.resseq_as_int()
      for chain in model.chains()
      for residue_group in chain.residue_groups()
    })
    resseq_options = sorted([int(e) for e in resseq_options ])
    resseq_options = [str(e) for e in resseq_options]
    atom_options = sorted(list(set([a.name.strip() for a in model.get_atoms()])))
    alts_options = sorted(list(set([a.parent().altloc for a in model.get_atoms()])))

    return chain_options, resname_options, resseq_options, atom_options


  def is_specific(self,value):
    if value not in ["all",None]:
      return True
    else:
      return False

  def chain_change(self, selected):
    self.comp_scroller_controller.reset()
    self.seq_scroller_controller.reset()
    self.atom_scroller_controller.reset()

    #sel = Selection(model=self.state.active_model)
    self.chain_scroller_controller.selected
    comp = self.comp_scroller_controller.selected
    seq = self.seq_scroller_controller.selected 
    atom = self.atom_scroller_controller.selected
    # Update residues and atoms
    #filters = {}
    sel_str_parts = []
    if self.is_specific(selected):
      sel_str_parts.append(f"chain {selected}")
      #filters["chain"] = selected
    if self.is_specific(seq):
      sel_str_parts.append(f"resseq {seq}")
      #filters["resseq"] = seq
    if self.is_specific(comp):
      sel_str_parts.append(f"resname {comp}")
      #filters["resname"] = comp
    if self.is_specific(atom): 
      sel_str_parts.append(f"name {atom}")
      #filters["name"] = atom
    if len(sel_str_parts)==0:
      sel_str = "all"
    else:
      sel_str = " and ".join(sel_str_parts)
    self.log("Sel str: ",sel_str)
    chain_options, resname_options, resseq_options, atom_options = self.update_options(sel_str)

    self.comp_scroller_controller.update_list_widget(resname_options)
    self.comp_scroller_controller.set_selected_item(comp)

    self.seq_scroller_controller.update_list_widget(resseq_options)
    self.seq_scroller_controller.set_selected_item(seq)

    self.atom_scroller_controller.update_list_widget(atom_options)
    self.atom_scroller_controller.set_selected_item(atom)

    all_options = resname_options +  resseq_options + atom_options
    if len(all_options)==0:
      self.reset()
    else:
      self.parent.viewer.select_from_selection_string(sel_str)
      self.parent.viewer.focus_selected()

  def comp_change(self, selected):
    self.seq_scroller_controller.reset()
    self.atom_scroller_controller.reset()
    chain = self.chain_scroller_controller.selected
    seq = self.seq_scroller_controller.selected 
    atom = self.atom_scroller_controller.selected
    comp = selected
    
    sel_str_parts = []
    if self.is_specific(comp):
      sel_str_parts.append(f"resname {comp}")
    if self.is_specific(chain):
      sel_str_parts.append(f"chain {chain}")
    if self.is_specific(seq):
      sel_str_parts.append(f"resseq {seq}")
    if self.is_specific(atom): 
      sel_str_parts.append(f"name {atom}")

    if len(sel_str_parts)==0:
      sel_str = "all"
    else:
      sel_str = " and ".join(sel_str_parts)
    self.log("Sel str: ",sel_str)
    chain_options, resname_options, resseq_options, atom_options = self.update_options(sel_str)

    self.seq_scroller_controller.update_list_widget(resseq_options)
    self.seq_scroller_controller.set_selected_item(seq)

    self.atom_scroller_controller.update_list_widget(atom_options)
    self.atom_scroller_controller.set_selected_item(atom)

    all_options = resseq_options + atom_options
    if len(all_options)==0:
      self.reset()
    else:
      self.parent.viewer.select_from_selection_string(sel_str)
      self.parent.viewer.focus_selected()

  def seq_change(self, selected):
    chain = self.chain_scroller_controller.selected
    comp = self.comp_scroller_controller.selected
    atom = self.atom_scroller_controller.selected
    
    if selected != "all":
      selected = int(selected)

    seq = selected
      
    sel_str_parts = []
    if self.is_specific(seq):
      sel_str_parts.append(f"resseq {seq}")
    if self.is_specific(chain):
      sel_str_parts.append(f"chain {chain}")
    if self.is_specific(comp):
      sel_str_parts.append(f"resname {comp}")

    if len(sel_str_parts)==0:
      sel_str = "all"
    else:
      sel_str = " and ".join(sel_str_parts)
    self.log("Sel str: ",sel_str)
    chain_options, resname_options, resseq_options, atom_options = self.update_options(sel_str)

    self.atom_scroller_controller.update_list_widget(atom_options)
    self.atom_scroller_controller.set_selected_item(atom)

    all_options = atom_options
    if len(all_options)==0:
      self.reset()
    else:
      self.parent.viewer.select_from_selection_string(sel_str)
      self.parent.viewer.focus_selected()

  def atom_change(self,selected):
    self.parent.viewer.set_granularity(value="element")
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
    self.log("Sel str: ",sel_str)
    self.parent.viewer.select_from_selection_string(sel_str)
    self.log("emitting atom")
    self.state.signals.picking_level.emit("atom")
    self.parent.viewer.focus_selected()

  def reset(self):
    self.chain_scroller_controller.reset()
    self.comp_scroller_controller.reset()
    self.seq_scroller_controller.reset()
    self.atom_scroller_controller.reset()
    self.parent.viewer.deselect_all()

class ScrollableHierarchyController(Controller):
  """
  A controller which manages a single scroll window for a category
    of hierarchy (chain, resname, resseq, etc)
  """
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

    