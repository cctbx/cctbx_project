"""
The controller for the molstar web app
"""
import time

from libtbx.test_utils import approx_equal
from libtbx.utils import null_out, Sorry
from PySide2.QtCore import QTimer
from PySide2.QtWidgets import QMessageBox

from .controller import Controller
from ...molstar.molstar import MolstarGraphics
from .molstar_controls_base import MolstarControlsController
from ..controller.molstar_controls_base import MolstarControlsController
from ..model.selection import Selection
from ..model.ref import ModelRef, SelectionRef

class sync_manager:
  """
  Use this class to run commands which require syncronization with the typescript app
  """
  def __init__(self, molstar_controller):
    self.molstar_controller = molstar_controller
    
  def __enter__(self):
    pass

  def __exit__(self, exc_type, exc_val, exc_tb):
    # Start timer when complete
    self.molstar_controller.sync_timer.start()


class MolstarController(Controller):
  """
  A controller which manages external graphics for Molstar. However it does not actually
    produce any JavaScript code. That is delegated to its 'graphics' member variable, and 
    instance of the MolstarGraphics class. 

  So an example program flow would be
    1. This class responds to a signal or has a method directly called
    2. This class calls the corresponding method on self.graphics (MolstarGraphics instance)
    3. The MolstarGraphics instance runs Javascript code in the QtWebView to affect the Molstar webapp

  In this way there is a separation between the Qt GUI infrastructure (signals and slots) and 
    the Python-to-Javascript api present in MolstarGraphics. No other class besides this 
    one (MolstarController) should interact directly with the MolstarGraphics instance.

  """
  def __init__(self,parent=None,view=None,web_view=None):
    super().__init__(parent=parent,view=view)
    self._web_view = web_view
    
    self.graphics = MolstarGraphics(self.web_view)
    self.graphics.state = self.state
    if self.view:
      self.graphics_controls = MolstarControlsController(parent=self,view=self.view.viewer_controls)


    self._blocking_commands = True
    self._picking_granularity = "residue"
    self.references_remote_map = {} # Local ref_id : molstar_ref_id
    #self.sync_manager = SyncManager(self)


    # Signals
    self.graphics.web_view.loadStarted.connect(self._on_load_started)
    self.graphics.web_view.loadFinished.connect(self._on_load_finished_pre_sync)
    self.state.signals.has_synced.connect(self._on_load_finished_post_sync)

    # Maps/Models
    self.state.signals.model_change.connect(self.load_model_from_ref)

    # Selections
    self.state.signals.select_all.connect(self.select_all)
    self.state.signals.selection_activated.connect(self.select_from_ref)
    self.state.signals.deselect_all.connect(self.deselect_all)
    self.state.signals.selection_deactivated.connect(self.deselect_from_ref)
    self.state.signals.selection_focus.connect(self.focus_selected)


    #timer for update
    self.sync_timer = QTimer()
    self.sync_timer.setInterval(2000)
    self.sync_timer.timeout.connect(self._update_state_from_remote)
    self.timer_accumulator = 0
    self.timer_max_retries = 10

    # Start by default
    self.start_viewer()

  @property
  def web_view(self):
    if self._web_view:
      return self._web_view
    else:
      return self.view.web_view

  def _on_load_started(self):
    self._blocking_commands = True
    self.graphics._blocking_commands = True
    if self.view:
      self.view.parent_explicit.setEnabled(False)

  def _on_load_finished_pre_sync(self, ok):
    if ok:
      self.log("Page loaded successfully. Waiting for sync")
      self.sync_timer.start()
    else:
      self.log("An error occurred while loading web app.")
      self._blocking_commands = True

  def _on_load_finished_post_sync(self, ok):
    if self.view:
      self.view.parent_explicit.setEnabled(True) # enable gui
    self._blocking_commands = False
    self.graphics._blocking_commands = False
    self._load_all_from_ref()

  def _update_state_from_remote(self):
    self.state.molstarState = self.graphics.sync_remote()
    if self.state.molstarState and self.state.molstarState.has_synced:
      self.sync_timer.stop()
    self.timer_accumulator+=1
    if self.timer_accumulator>self.timer_max_retries:
      self.sync_timer.stop() # give up

  # API
  def start_viewer(self):
    self.graphics.start_viewer()

  def log_message(self,message):
    self.graphics.log_message(message)


  def _load_all_from_ref(self,*args):
    # Molstar must load inital models here
    for ref in self.state.references_model:
      self.load_model_from_ref(ref)


  def load_model_from_ref(self,ref,label=None,format='pdb'):
    if not isinstance(ref,ModelRef):
      return None
    if self._blocking_commands:
      print("Skipping load_model_from_ref due to self._blocking_commands=True")
      return
    if ref is not None:
      # Check if already loaded
      if "molstar" in ref.identifiers and "molstar":
        if ref.identifiers["molstar"] is not None and ref.identifiers["molstar"] in self.state.external_loaded["molstar"]:
          self.log("Not loading model because already loaded into molstar")
          self.log(self.state.external_loaded)
          self.log(ref.identifiers["molstar"])
        return

      # Ref needs to be loaded
      with sync_manager(self):
        #self.graphics._set_sync_state(self.state.to_json())

        self.graphics.load_model_from_mmtbx(
          model=ref.data.model, # mmtbx model
          format=format,
          ref_id=ref.uuid,
          label=label,
          callback=None
        )
        #self.sync_manager.has_synced = False


  # Selection
  def sync_selection(self):
    """
    Poll selection from molstar viewer and set state. 
    Returns True if it produced a valid non-empty selection
    """
    selection = self.poll_selection()
    if selection.is_empty:
      self.state.active_selection_ref = None
      return False
    else:
      sel = self.state.active_mol.select_from_selection(selection)
      if len(sel)>0:
        sel_ref = SelectionRef(selection,model_ref=self.state.active_model_ref,show=False)
        self.state.add_ref(sel_ref)
        self.state.active_selection_ref = sel_ref
        return True
      else:
        self.log("Skipping add selection due to empty selection")
        return False

  def poll_selection(self):
    """
    Poll selection from molstar viewer
    Returns selection object
    """
    atom_records = self.graphics.poll_selection()
    selection = Selection.from_atom_records(atom_records,self.state.active_model)
    return selection
  
  def test_selection_string(self,phenix_string):
    # Do raw selection
    model = self.state.active_model
    xyz1 = model.get_hierarchy().atoms().extract_xyz()

    # Select in molstar and poll back
    self.select_from_selection_string(phenix_string)
    selection = self.poll_selection()
    model2 = model.select(model.selection(selection.string))
    xyz2 = model2.get_hierarchy().atoms().extract_xyz()
    if approx_equal(xyz1,xyz1,eps=3):
      return True
    else:
      return False
    



  def select_from_ref(self,ref: SelectionRef):
    return self.select_from_selection(ref.data)
  def deselect_from_ref(self,ref:SelectionRef):
    return self.deselect_all()
    
  def select_from_selection(self,selection: Selection):
    return self.graphics.select_from_selection(selection)

  def select_from_selection_string(self,selection_string: str):
    selection = Selection.from_selection_string(selection_string,model=self.state.active_model)
    return self.graphics.select_from_selection(selection)

  def focus_selected(self):
    return self.graphics.focus_selected()

  def select_all(self):
    self.graphics.select_all()

  def deselect_all(self):
    self.graphics.deselect_all()

  def close_viewer(self):
    self.graphics.close_viewer()

  def picking_level(self,picking_level):
    if 'atom' in picking_level:
      self.set_granularity("element")
    elif "residue" in picking_level:
      self.set_granularity("residue")
    else:
      pass

  def set_granularity(self,value="residue"):
    assert value in ['element','residue'], 'Provide one of the implemented picking levels'
    self._picking_granularity = value
    self.viewer._set_granularity(value=value)