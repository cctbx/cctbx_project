from PySide2.QtCore import Slot

from .controller import Controller
from ...last.selection_utils import SelectionQuery
from ..state.ref import SelectionRef

class SelectionControlsController(Controller):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)
    # Enable return key to execute selection
    self.view.selection_edit.returnPressed.connect(self.execute_selection)
    self.view.selector_toggle.clicked.connect(self.toggle_selection)
    # self.view.start_selecting.clicked.connect(self.start_selecting)
    # self.view.button_deselect.clicked.connect(self.deselect)


    self.view.load_button.clicked.connect(self.add_selection)
    self.view.combo_box.currentIndexChanged.connect(self.on_picking_change)
    self.view.button_clear.clicked.connect(self.clear_viewer) # TODO: Move out of selection controls
    self.state.signals.picking_level.connect(self.on_picking_change)

  @property
  def viewer(self):
    return self.parent
  
  def clear_viewer(self):
    self.parent.clear_viewer("Clearing the viewer.")

  def toggle_selection(self):
    if not self.view.selector_toggle.is_on:
      self.start_selecting()
    else:
      self.deselect()


  @Slot()
  def select_active_selection(self):
    #self.highlight_persist(None,value=True) # set the highlight persist automatically when selecting
    if self.state.active_selection_ref is None:
      self.viewer.deselect_all()
    else:
      self.viewer.select_from_json(self.state.active_selection_ref.query.to_json())

  def deselect(self):
    self.state.active_selection = None
    self.viewer.deselect_all()
    self.viewer.toggle_selection_mode(False)

  @Slot()
  def start_selecting(self):
    # selection button was clicked
    print("clicked")
    self.viewer.toggle_selection_mode(True)
    self.execute_selection()
  
  @Slot()
  def execute_selection(self):
    """
    This is a selection from the text box
    """    
    text = self.view.selection_edit.text()
    if text:
      try:
        if text.startswith("select"):
          text = text[7:]
        elif text.startswith("sel "):
          text = text[4:]
        self.viewer.select_from_phenix(text)
        self.save_text_to_history()
      except:
        raise
        self.view.selection_edit.clear()
        self.view.selection_edit.setPlaceholderText(f"Unable to interpret selection: {text}")

  def exec_selection_box(self):
    sel_str = self.view.selection_edit.text()
    if len(sel_str.strip())>0:
      query = self.state.mol.atom_sites.select_as_query(sel_str)
      return query

  
  
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
    print("controls add_selection()")
    self.viewer.poll_selection(callback=self._add_selection)
    
  def _add_selection(self,query_dict):
    print(f"controls _add_selection() callback for query_dict: {query_dict}")
    assert len(query_dict)==1, "Multi structure queries not yet supported"
    ref_id = list(query_dict.keys())[0]
    query = query_dict[ref_id]
    if ref_id not in self.state.references:
      ref_id = self.state.active_model_ref.id
    ref = self.state.references[ref_id]
    sel_sites = self.state.active_mol.atom_sites.select_from_query(query)
    if len(sel_sites)>0:
      sel_ref = SelectionRef(data=query,model_ref=ref,show=True)
      self.state.add_ref(sel_ref)
      self.state.signals.tab_change.emit("Selections") # show selection tab
    else:
      print("Skipping add selection due to empty selection")
      

