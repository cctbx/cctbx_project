from pathlib import Path
import json
import re
import copy
import time
import sys
from functools import partial
from typing import Optional
from contextlib import contextmanager

import pandas as pd
from PySide2.QtCore import QUrl, Signal, QObject, QTimer
from PySide2.QtWidgets import QFileDialog, QMessageBox


from .controller import Controller
from ...molstar.molstar import MolstarViewer
from .style import ModelStyleController, MapStyleController
#from ..controller.selection_controls import SelectionControlsController
from ..controller.viewer_controls import ViewerControlsController
from ...core.selection_utils import Selection, SelectionQuery
from ...core.python_utils import DotDict

class sync_manager:
  """
  Use this to run commands which require syncronization with the typescript app
  """
  def __init__(self, molstar_controller):
      self.molstar_controller = molstar_controller
      self.state = self.molstar_controller.state

  def __enter__(self):
      print("Entering sync manager")
      self.state.has_synced = False


  def __exit__(self, exc_type, exc_val, exc_tb):
      self.molstar_controller.timer.start()


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
    #self.selection_controls = SelectionControlsController(parent=self,view=self.view.selection_controls)
    self.viewer_controls = ViewerControlsController(parent=self,view=self.view.viewer_controls)

    self._blocking_commands = False
    self._picking_granularity = "residue"
    self.references_remote_map = {} # Local ref_id : molstar_ref_id
    #self.sync_manager = SyncManager(self)


    # Signals
    self.viewer.web_view.loadStarted.connect(self._on_load_started)
    self.viewer.web_view.loadFinished.connect(self._on_load_finished)

    self.state.signals.model_change.connect(self.load_model_from_ref)
    self.state.signals.map_change.connect(self.load_map_from_ref)

    self.state.signals.select.connect(self.select_from_ref)
    self.state.signals.clear.connect(self.clear_viewer)
    self.state.signals.selection_change.connect(self.select_from_ref)

    #timer for update
    self.timer = QTimer()
    self.timer.setInterval(2000)
    self.timer.timeout.connect(self._update_state_from_remote)
    self.timer.start()
    self.timer_accumulator = 0
    self.timer_max_retries = 10

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



  def _update_state_from_remote(self):
    self.state.phenixState = self.viewer.sync_remote()
    if self.state.has_synced:
      self.timer.stop()
    self.timer_accumulator+=1
    if self.timer_accumulator>self.timer_max_retries:
      self.timer.stop() # give up

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

    if ref is not None:
      if ref.id_molstar is not None and ref.id_molstar in self.state.external_loaded["molstar"]:
        print("Not loading model because already loaded into molstar")
        print(self.state.external_loaded)
        print(ref.id_molstar)

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

  def poll_selection(self,callback=None):
    return self.viewer.poll_selection(callback=callback)

  def _poll_selection_callback(self,callback,selection_json):
    #print("Calling MolstarController._poll_selection_callback(callback, selection_json)")
    #print(f"where\n\tcallback:{callback}\n\n\tselection_json:{selection_json}")
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

  def select_from_phenix_string(self,selection_phenix=None):
    # convert to query
    query = self.state.mol.sites._select_query_from_str_phenix(selection_phenix)
    query.params.refId = self.state.active_model_ref.id
    self.viewer.select_from_query(query)

  def select_from_ref(self,ref):
    self.viewer.select_from_query(ref.query)

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
    phenix_ref = self.state.phenixState.references[ref.id]
    model_id = phenix_ref.id_molstar
    representations = phenix_ref.representations
    for rep_name in representations:
      self.viewer.show_query(model_id,ref.query.to_json(),rep_name)


  def hide_ref(self,ref,representation: Optional[str] = None):
    phenix_ref = self.state.phenixState.references[ref.id]
    model_id = phenix_ref.id_molstar
    representations = phenix_ref.representations
    for rep_name in representations:
      self.viewer.hide_query(model_id,ref.query.to_json(),rep_name)

  def color_ref(self,ref,color):
    model_id = ref.model_ref.id_molstar
    self.viewer.color_query(model_id,ref.query.to_json(),color)
