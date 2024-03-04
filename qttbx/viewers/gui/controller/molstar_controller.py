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
from ...last.python_utils import DotDict

class SyncManager:
  """
  Manages the synchronization of state for the molstar typescript app
  """
  def __init__(self,controller):
    self.controller= controller
    self._has_synced = None
    self._sync_timer = QTimer()
    self._sync_timer.timeout.connect(self.try_sync)  # Connect the timeout signal to my_function
    self._sync_interval = 500
    self._sync_count = 0
    self._max_sync_count = 20
    self._sync_failure = False
    self._received_first_sync = False


  @property
  def has_synced(self):
    return self._has_synced

  @has_synced.setter
  def has_synced(self,value):
    """
    Change the sync state
    """
    if (value is True and (self._has_synced is False or self._has_synced is None)):
      self._has_synced = value
      self.controller.state.has_synced = value
      self._sync_timer.stop()
    elif (value is False and (self._has_synced is True or self._has_synced is None)):
      self._has_synced = value
      self.controller.state.has_synced = value
      self._sync_timer.start(self._sync_interval)
    elif self._has_synced is True and value is True:
      assert False, "Error, should not reach here"


  # def check_alive(self):
  #   self.controller.viewer._get_sync_state(self.controller.state.to_json(),callback=self._check_alive_callback)

  # def _check_alive_callback(self,result):
  #   if isinstance(result,str) and len(result.strip())>0:
  #     result = json.loads(result)
  #   if "hasSynced" in result:
  #     self._is_alive = True
  #   else:
  #     self._is_alive = False


  def try_sync(self,callback=None):
    self._has_synced = False
    # callback is run if successful sync
    #print(f"try_sync() ... sync called with count: {self._sync_count} at time: {time.time()}, has synced: {self.has_synced}")
    if self.has_synced:
      print("case1: has synced")
      pass
    elif not self.has_synced and self._sync_count < self._max_sync_count:
      print("case 2: keep trying to sync")
      self._sync_count+=1
      self.controller.viewer._get_sync_state(callback=self._try_sync_callback)
    else:
      print("case 3: sync failure")
      self._sync_failure = True
      msg = QMessageBox.warning(self.controller.view,"Warning", "The GUI and the molstar viewer are out of sync, program will exit.")
      self.controller.parent.close_application()
      sys.exit()

  def _try_sync_callback(self,result,second_callback=None): #  get_state function in ts
    received_result = False
    print("sync callback result:",result,", type: ",type(result), ",time: ",time.time())
    if isinstance(result,str) and len(result.strip())>0:
      result = json.loads(result)
      print("sync result: ",json.dumps(result,indent=2))
      if 'hasSynced' in result:
        received_result = True
        if not self._received_first_sync:
          self._received_first_sync = True
          result["hasSynced"] = True # The very first sync is True
          print("Setting result['hasSynced'] = True because it is the first sync state")
          print("sync callback result:",result,", type: ",type(result), ",time: ",time.time())


    # update molstar ids
    if received_result:
      if 'references' in result: # if not, not well formed response
        for ref_id, ref_dict in result["references"].items():
          for external_key, external_id in ref_dict["external_ids"].items():
            ref = self.controller.state.references[ref_id]
            ref.external_ids[external_key] = external_id
            self.controller.state.external_loaded["molstar"].append(ref.id)

    if received_result:
      if result["hasSynced"]:
        self.has_synced = True
        self._sync_count = 0
        self.controller.viewer._set_sync_state(self.controller.state.to_json())

    # catch failure
    if self._sync_failure:
      msg = QMessageBox.warning(self.view,"Warning", "The GUI and the molstar viewer are out of sync, program will exit.")
      self.controller.parent.close_application()
      sys.exit()

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
    self.model_style_controller = ModelStyleController(parent=self,view=None)
    self.map_style_controller = MapStyleController(parent=self,view=None)
    self.selection_controls = SelectionControlsController(parent=self,view=self.view.selection_controls)

    self._blocking_commands = False
    self._picking_granularity = "residue"
    self.references_remote_map = {} # Local ref_id : molstar_ref_id
    self.sync_manager = SyncManager(self)


    # Signals
    self.viewer.web_view.loadStarted.connect(self._on_load_started)
    self.viewer.web_view.loadFinished.connect(self._on_load_finished)

    self.state.signals.model_change.connect(self.load_model_from_ref)
    self.state.signals.map_change.connect(self.load_map_from_ref)

    self.state.signals.select.connect(self.select_from_ref)
    self.state.signals.clear.connect(self.clear_viewer)


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
      time.sleep(0.5)
      self._load_all_from_ref()

    else:
      print("An error occurred while loading web app.")
      self._blocking_commands = True




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
    if ref is not None and ref.id not in self.state.external_loaded["molstar"]:
        self.viewer._set_sync_state(self.state.to_json())

        self.viewer.load_model_from_mmtbx(
          model=ref.model,
          format=format,
          ref_id=ref.id,
          label=label,
          callback=None
        )
        self.sync_manager.has_synced = False


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
        self.viewer._set_sync_state(self.state.to_json())
        if model_ref is None and map_ref.model_ref is None:
          map_ref.model_ref = self.state.active_model_ref

        self.viewer.load_map(filename=map_ref.data.filepath,volume_id=map_ref.id,model_id=map_ref.model_ref.id)
        self.sync_manager.has_synced = False


  # Selection

  def poll_selection(self,callback=None):
    self.viewer.poll_selection(callback=partial(self._poll_selection_callback,callback))

  def _poll_selection_callback(self,callback,selection_json):
    #print("Calling MolstarController._poll_selection_callback(callback, selection_json)")
    #print(f"where\n\tcallback:{callback}\n\n\tselection_json:{selection_json}")
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

  def select_from_phenix_string(self,selection_phenix=None):
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
