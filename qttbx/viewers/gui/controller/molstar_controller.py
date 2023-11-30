from pathlib import Path
import json
import re
import copy
import time
import sys
from functools import partial
from typing import Optional

import pandas as pd
from PySide2.QtCore import QUrl, Signal, QObject, QTimer
from PySide2.QtWidgets import QFileDialog, QMessageBox


from .controller import Controller
from ...molstar.molstar import MolstarViewer
from .style import ModelStyleController, MapStyleController
from ..controller.selection_controls import SelectionControlsController
from ...last.selection_utils import Selection, SelectionQuery

class MolstarController(Controller):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)

    self.viewer = MolstarViewer(self.view.web_view)
    self.model_style_controller = ModelStyleController(parent=self,view=None)
    self.map_style_controller = MapStyleController(parent=self,view=None)
    self.selection_controls = SelectionControlsController(parent=self,view=self.view.selection_controls)

    # local state
    self.loaded_map_refs = []
    self.loaded_model_refs = []
    self._blocking_commands = False
    self._picking_granularity = "residue"
    self.references_remote_map = {} # Local ref_id : molstar_ref_id
    self._has_synced = None
    self._sync_timer = QTimer()
    self._sync_timer.timeout.connect(self.sync)  # Connect the timeout signal to my_function
    self._sync_interval = 500
    self._sync_count = 0
    self._max_sync_count = 20
    self._sync_failure = False

    # Signals
    self.viewer.web_view.loadStarted.connect(self._on_load_started)
    self.viewer.web_view.loadFinished.connect(self._on_load_finished)

    self.state.signals.model_change.connect(self.load_model_from_ref)
    self.state.signals.map_change.connect(self.load_active_map)

    self.state.signals.select.connect(self.select_from_ref)
    self.state.signals.clear.connect(self.clear_viewer)




    # # self.state.signals.selection_change.connect(self.selection_controls.select_active_selection)
    # # self.state.signals.color_change.connect(self.set_color)
    # # self.state.signals.repr_change.connect(self.set_representation)
    # # self.state.signals.viz_change.connect(self.set_visibility)
    # self.state.signals.data_change.connect(self._data_change)

    # Start by default
    self.start_viewer()

  def _on_load_started(self):
    self._blocking_commands = True
    self.viewer._blocking_commands = True
    self.view.parent_explicit.setEnabled(False)

  def _on_load_finished(self, ok):
    
    self.view.parent_explicit.setEnabled(True)
    if ok:
      print("Page loaded successfully. Accepting commands")
      self._blocking_commands = False
      self.viewer._blocking_commands = False
      # Debugging
      # Set the loaded model to be visible
      # self.parent.data.model_list_controller.entries[0].toggle_active_func(True)
    else:
      print("An error occurred while loading web app.")
      self._blocking_commands = True


  # API



  def start_viewer(self):
    self.viewer.start_viewer()


  # def _data_change(self):
  #   self.volume_streamer.pack_data_manager(self.state.data_manager)



  def load_model_from_ref(self,ref,label=None,format='pdb'):
    if ref is not None and ref not in self.loaded_model_refs:
      print("molstar controller adding model from ref: ",ref.id)
      self.viewer._load_model_from_mmtbx(
        model=ref.model,
        format=format,
        ref_id=ref.id,
        label=label,
        callback=None
      )

      # important! add ref to local state
      self.loaded_model_refs.append(ref)
      
      self.has_synced = False # flag that we do not know the remote ref yet
    


  @property
  def has_synced(self):
    return self._has_synced

  @has_synced.setter
  def has_synced(self,value):
    print("Setting has_synced")
    if self._has_synced is False and value is True:
      self._has_synced = value
      self._sync_timer.stop()
    elif ((self._has_synced is True and value is False) or self._has_synced is None):
      self._has_synced = value
      self._sync_timer.start(self._sync_interval)
    elif self._has_synced is True and value is True:
      assert False, "Error, should not reach here"


  # Maps

  def load_active_map(self,ref):
    assert ref is None or ref is self.state.active_map_ref
    if ref is not None and ref not in self.loaded_map_refs:
      self.load_map_from_ref(ref)


  def load_map_from_ref(self,map_ref):
    if map_ref.model_ref is None:
      map_ref.model_ref = self.state.active_model_ref
    self.viewer.load_map(filename=map_ref.data.filepath,volume_id=map_ref.id,model_id=map_ref.model_ref.id)
    self.loaded_map_refs.append(map_ref)
    self.has_synced = False


  # Selection

  def poll_selection(self,callback=None):


    self.viewer.poll_selection(callback=partial(self._poll_selection,callback))
     #self.viewer.poll_selection(callback=None)

  def _poll_selection(self,callback,selection_json):
    if not isinstance(selection_json,str) or len(selection_json.strip())==0:
      return None
    query_atoms = SelectionQuery.from_json(selection_json)
    
    # get the relevant model ref (to get mol)
    ref_id = query_atoms.params.refId
    model_ref = self.state.references[ref_id]

    # select sites
    sites_sel = model_ref.mol.atom_sites.select_from_query(query_atoms)

    # make condensed query
    query_condensed = model_ref.mol.atom_sites._convert_sites_to_query(sites_sel)
    query_condensed.params.refId = ref_id

    query_dict = {ref_id:query_condensed}
    if callback is not None:
      callback(query_dict)

  def select_from_phenix_string(self,selection_phenix: str):
    # convert to query
    query = self.state.mol.sites._select_query_from_str_phenix(selection_phenix)
    query.params.refId = self.state.active_model_ref.id
    self.viewer.select_from_query(query.to_json())

  def select_from_ref(self,ref):
    self.viewer.select_from_query(ref.query.to_json())

  def set_granularity(self,value="residue"):
    assert value in ['element','residue'], 'Provide one of the implemented picking levels'
    self._picking_granularity = value
    self.viewer._set_granularity(value=value)

  def toggle_selection_mode(self,value):
    self.viewer._toggle_selection_mode(value)

  def deselect_all(self):
    self.viewer.deselect_all()

  # Other

  def clear_viewer(self,msg):

    # clear local state
    self.loaded_map_refs = []
    self.loaded_model_refs = []
    
    # Un-toggle any entry
    for ref in self.state.references.values():
      if ref.entry is not None:
        ref.entry.is_active = False

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

  def sync(self):
    print(f"sync called with count: {self._sync_count} at time: {time.time()}, has synced: {self.has_synced}")
    if self.has_synced:
      pass
    elif not self.has_synced and self._sync_count < self._max_sync_count:
      self._sync_count+=1
      self.viewer._sync(callback=self._sync_callback)
      
    else:
      
      print("DEBUG: warning, have not synced")
      self._sync_failure = True
      msg = QMessageBox.warning(self.view,"Warning", "The GUI and the molstar viewer are out of sync, program will exit.")
      self.parent.close_application()
      sys.exit()

  def _sync_callback(self,result):
    print("sync callback result:",result,", type: ",type(result), ",time: ",time.time())
    if isinstance(result,str) and len(result.strip())>0:
      result = json.loads(result)
      print("sync result: ",result)

      has_synced = result["hasSynced"]
      if has_synced:
        self.has_synced = has_synced
        self._sync_count = 0
        model_result = json.loads(result["refMapping"])
        for remote_id,local_id in model_result.items():
          assert local_id in self.state.references, "Recieved local ref_id from remote that is not known"
          if local_id not in self.references_remote_map:
            self.references_remote_map[local_id] = remote_id
        volume_result = result["volumeEntries"]
        if isinstance(volume_result,str):
          volume_result = json.loads(volume_result)
        print("volume_result_from_sync_callback:")
        print(volume_result)

      # sync references
      for ref_id,ref in self.state.references.items():
        if ref.id in self.references_remote_map:
          ref.external_ids["molstar"] = self.references_remote_map[ref_id]

      # catch failure
      if self._sync_failure:
        msg = QMessageBox.warning(self.view,"Warning", "The GUI and the molstar viewer are out of sync, program will exit.")
        self.parent.close_application()
        sys.exit()
        
  # Style
  def set_iso(self,ref,value):
    # This one is unusual in that remote uses local id
    self.viewer.set_iso(ref.id,value)


  def show_ref(self,ref,representation: Optional[str] = None):
    model_id = ref.model_ref.external_ids["molstar"]
    if representation is None:
      representation = ref.style.representation
    else:
      representation = [representation]

    for rep_name in representation:
      self.viewer.show_query(model_id,ref.query.to_json(),rep_name)

  def hide_ref(self,ref,representation: Optional[str] = None):
    model_id = ref.model_ref.external_ids["molstar"]
    if representation is None:
      representation = ref.style.representation
    else:
      representation = [representation]

    for rep_name in representation:
      self.viewer.hide_query(model_id,ref.query.to_json(),rep_name)

  def color_ref(self,ref,color):
    model_id = ref.model_ref.external_ids["molstar"]
    self.viewer.color_query(model_id,ref.query.to_json(),color)

    