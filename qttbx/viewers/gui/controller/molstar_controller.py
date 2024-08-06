import time

from PySide2.QtCore import QTimer
from PySide2.QtWidgets import QMessageBox

from .controller import Controller
from ...molstar.molstar import MolstarViewer
from .style import ModelStyleController, MapStyleController
from ..controller.viewer_controls import ViewerControlsController
from ...core.selection import Selection
from ..state.ref import SelectionRef

class sync_manager:
  """
  Use this to run commands which require syncronization with the typescript app
  """
  def __init__(self, molstar_controller):
      self.molstar_controller = molstar_controller
      self.state = self.molstar_controller.state
      self.debug=True

  def log(self,*args):
    if self.debug:
      print(*args)
    
  def __enter__(self):
      self.log("Entering sync manager")
      #self.state.has_synced = False


  def __exit__(self, exc_type, exc_val, exc_tb):
      self.molstar_controller.sync_timer.start()


class MolstarController(Controller):
  api_function_names = [
    "load_model",
    "load_model_from_string",
    "load_map",
    "log_message",
    "select_from_phenix_string",
  ]
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)

    self.viewer = MolstarViewer(self.view.web_view)
    self.viewer.state = self.state
    self.model_style_controller = ModelStyleController(parent=self,view=None)
    self.map_style_controller = MapStyleController(parent=self,view=None)
    self.viewer_controls = ViewerControlsController(parent=self,view=self.view.viewer_controls)


    self._blocking_commands = True
    self._picking_granularity = "residue"
    self.references_remote_map = {} # Local ref_id : molstar_ref_id
    #self.sync_manager = SyncManager(self)


    # Signals
    self.viewer.web_view.loadStarted.connect(self._on_load_started)
    self.viewer.web_view.loadFinished.connect(self._on_load_finished_pre_sync)
    self.state.signals.has_synced.connect(self._on_load_finished_post_sync)


    self.state.signals.model_change.connect(self.load_model_from_ref)
    self.state.signals.map_change.connect(self.load_map_from_ref)

    # Selections
    self.state.signals.select_all.connect(self.select_all)
    self.state.signals.selection_activated.connect(self.select_from_ref)


    self.state.signals.deselect_all.connect(self.deselect_all)
    self.state.signals.selection_deactivated.connect(self.deselect_from_ref)
    
    self.state.signals.selection_hide.connect(self.selection_hide)
    self.state.signals.selection_show.connect(self.selection_show)
    self.state.signals.selection_rep_show.connect(self.selection_show)

    self.state.signals.selection_focus.connect(self.focus_selected)
    self.state.signals.clear.connect(self.clear_viewer)
    self.state.signals.picking_level.connect(self.picking_level)

    # Styles
    self.state.signals.color_selection.connect(self.color_selection)
    self.state.signals.representation_selection.connect(self.representation_selection)

    #timer for update
    self.sync_timer = QTimer()
    self.sync_timer.setInterval(2000)
    self.sync_timer.timeout.connect(self._update_state_from_remote)
    self.timer_accumulator = 0
    self.timer_max_retries = 10

    # Start by default
    self.start_viewer()

  def _on_load_started(self):
    self._blocking_commands = True
    self.viewer._blocking_commands = True
    self.view.parent_explicit.setEnabled(False)

  def _on_load_finished_pre_sync(self, ok):
    if ok:
      self.log("Page loaded successfully. Waiting for sync")
      self.sync_timer.start()
    else:
      self.log("An error occurred while loading web app.")
      self._blocking_commands = True

  def _on_load_finished_post_sync(self, ok):
    self.view.parent_explicit.setEnabled(True) # enable gui
    self._blocking_commands = False
    self.viewer._blocking_commands = False
    self._load_all_from_ref()

  def _update_state_from_remote(self):
    self.state.phenixState = self.viewer.sync_remote()
    if self.state.phenixState and self.state.phenixState.has_synced:
      self.sync_timer.stop()
    self.timer_accumulator+=1
    if self.timer_accumulator>self.timer_max_retries:
      self.sync_timer.stop() # give up

  # API
  def start_viewer(self):
    self.viewer.start_viewer()

  def log_message(self,message):
    self.viewer.log_message(message)



  def _load_all_from_ref(self,*args):
    # Molstar must load inital models here
    for ref in self.state.references_model:
      self.load_model_from_ref(ref)
    for ref in self.state.references_map:
      self.load_map_from_ref(ref)



  def load_model_from_ref(self,ref,label=None,format='pdb'):
    if self._blocking_commands:
      print("Skipping load_model_from_ref due to self._blocking_commands=True")
      return
    if ref is not None:
      if ref.id_molstar is not None and ref.id_molstar in self.state.external_loaded["molstar"]:
        self.log("Not loading model because already loaded into molstar")
        self.log(self.state.external_loaded)
        self.log(ref.id_molstar)

        return
      # Ref needs to be loaded
      with sync_manager(self):
        #self.viewer._set_sync_state(self.state.to_json())

        self.viewer.load_model_from_mmtbx(
          model=ref.model,
          format=format,
          ref_id=ref.id,
          label=label,
          callback=None
        )
        #self.sync_manager.has_synced = False



  def load_model_from_string(self,model_str=None,label=None,format='pdb'):
    if model_str is not None:
      # make a ref first
      ref = self.state.add_ref_from_model_string(model_str,label=label,format=format)
      self.load_model_from_ref(ref,label=label,format=format)

  def load_model(self,filename=None,label=None,format='pdb'):
    if filename is not None:
      if label is None:
        label = filename
      # make a ref first
      ref = self.state.add_ref_from_model_file(filename=filename,label=label,format=format)
      self.load_model_from_ref(ref,label=label,format=format)

  # Maps
  def load_map(self,filename=None,volume_id=None,model_id=None,label=None):
    if volume_id is None:
      volume_id = filename
    ref = self.state.add_ref_from_map_file(filename=filename,volume_id=volume_id,model_id=model_id,label=label)
    if ref.model_ref is None:
      if self.state.active_model_ref is None:
        if len(self.state.references_model)>0:
          self.state.active_model_ref = self.state.references_model[0]
    if ref.model_ref is None and self.state.active_model_ref is not None:
      ref.model_ref = self.state.active_model_ref
    if ref.model_ref is not None:
      self.load_map_from_ref(map_ref=ref)


  def load_map_from_ref(self,map_ref=None,model_ref=None):

    if map_ref is not None and map_ref.id not in self.state.external_loaded["molstar"]:
       with sync_manager(self):
          if model_ref is None and map_ref.model_ref is None:
            map_ref.model_ref = self.state.active_model_ref

          self.viewer.load_map(filename=map_ref.data.filepath,volume_id=map_ref.id,model_id=map_ref.model_ref.id)


  # Selection
  @property
  def n_atoms_selected(self):
    synced = self.sync_selection()
    if not synced:
      return 0
    ref = self.state.active_selection_ref
    sel_sites = self.state.mol.sites.select_from_selection(ref.selection)
    return len(sel_sites)
    
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
    return self.viewer.poll_selection()

  def _poll_selection_callback(self,callback,selection_json):
    #self.log("Calling MolstarController._poll_selection_callback(callback, selection_json)")
    #self.log(f"where\n\tcallback:{callback}\n\n\tselection_json:{selection_json}")
    assert isinstance(selection_json,str) and len(selection_json.strip())>0, "Failure to recieve selection json"
    query_atoms = SelectionQuery.from_json(selection_json)

    # get the relevant model ref (to get mol)
    ref_id = query_atoms.params.refId
    model_ref = self.state.references[ref_id]
    query_atoms.params.keymap = self.state.mmcif_column_map

    # select sites
    sites_sel = model_ref.mol.atom_sites.select_from_query(query_atoms)

    # make condensed query
    query_condensed = model_ref.mol.atom_sites._convert_sites_to_query(sites_sel)
    query_condensed.params.refId = ref_id

    query_dict = {ref_id:query_condensed}
    if callback is not None:
      callback(query_dict)

  def select_from_phenix_string(self,phenix_string):
    return self.viewer.select_from_phenix_string(phenix_string)

  def select_from_ref(self,ref: SelectionRef):
    return self.select_from_selection(ref.data)
  def deselect_from_ref(self,ref:SelectionRef):
    return self.deselect_all()
    
  def select_from_selection(self,selection: Selection):
    return self.viewer.select_from_selection(selection)

  def focus_selected(self):
    return self.viewer.focus_selected()

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

  def toggle_selection_mode(self,value):
    self.viewer._toggle_selection_mode(value)

  def select_all(self):
    self.viewer.select_all()

  def deselect_all(self):
    self.viewer.deselect_all()

  # Other

  def clear_viewer(self,msg):

    # clear local state
    self.state.external_loaded["molstar"] = []

    # Un-toggle any entry
    for ref in self.state.references.values():
      ref.entry = None

    # clear viewer
    self.viewer.clear_viewer()

    # reset active objects
    self.state.active_model_ref = None
    self.state.active_map_ref = None

    # just using the above causes problems within molstar, but only for web view

    self.view.web_view.reload() # so refresh the whole page. Not ideal, large visual disruption.
    QMessageBox.information(self.view, 'Notice', msg)



  def close_viewer(self):
    self.viewer.close_viewer()

  # Style
  def set_iso(self,ref,value):
    # This one is unusual in that remote uses local id
    self.viewer.set_iso(ref.id,value)



  def representation_selection(self,selection_ref,representation_name,show):
    if self.state.active_selection_ref != selection_ref:
      self.state.active_selection_ref = selection_ref
    if show:
      transparency = 0.0
    else:
      transparency = 1.0
    self.viewer.transparency_selection(selection_ref.selection,"all",representation_name,transparency)
  
  def selection_hide(self,selection_ref,representation_name=None):

    # Make sure it is selected (it probably already is)
    if self.state.active_selection_ref != selection_ref:
      self.state.active_selection_ref = selection_ref

    # Tell molstar to just hide what it has selected
    #self.viewer.selection_hide() # No arguments is all representations
    self.viewer.transparency_selection(selection_ref.selection,"all","all",1.0)

  
  def selection_show(self,selection_ref,representation_name=None):

    # Make sure it is selected (it probably already is)
    if self.state.active_selection_ref != selection_ref:
      self.state.active_selection_ref = selection_ref

    # Tell molstar to just hide what it has selected
    self.viewer.transparency_selection(selection_ref.selection,"all","all",0.0)
    #self.viewer.selection_show(representation_name=representation_name) # No arguments is all representations

  def color_selection(self,selection_ref,color):
    # Make sure it is selected (it probably already is)
    if self.state.active_selection_ref != selection_ref:
      self.state.active_selection_ref = selection_ref
    self.viewer.color_selection(selection_ref.selection,color)
