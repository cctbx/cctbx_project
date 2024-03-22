from pathlib import Path
import json
import re
import copy
from typing import Optional

import pandas as pd

from .controller import Controller
from ..state.ref import ModelRef, MapRef
from ...chimerax.chimerax import ChimeraXViewer
from .style import ModelStyleController, MapStyleController
from ..controller.selection_controls import SelectionControlsController
from ...core.selection_utils import Selection, SelectionQuery

class ChimeraXController(Controller):
  def __init__(self,parent=None,view=None):
    super().__init__(parent=parent,view=view)

    self.viewer = ChimeraXViewer() # Not formally a controller, but similar role
    self.viewer.controller = self
    self.selection_controls = SelectionControlsController(parent=self,view=self.view.selection_controls)

    # Local state
    self.loaded_map_refs = []
    self.loaded_model_refs = []
    #self.references_remote_map = {} # key: local ref_id, value: remote ref_id (ChimeraX model spec)
    self._picking_granularity = "residue"


    # Signals
    self.view.button_start.clicked.connect(self.start_viewer)
    self.state.signals.model_change.connect(self.load_model_from_ref)
    self.state.signals.map_change.connect(self.load_map_from_ref)
    self.state.signals.select.connect(self.select_from_ref)
    self.state.signals.clear.connect(self.clear_viewer)

    # # self.state.signals.selection_change.connect(self.selection_controls.select_active_selection)
    # # self.state.signals.color_change.connect(self.set_color)
    # # self.state.signals.repr_change.connect(self.set_representation)
    # # self.state.signals.viz_change.connect(self.set_visibility)
    # self.state.signals.data_change.connect(self._data_change)
    self.view.textInput.setPlaceholderText("Find automatically")


  def start_viewer(self):
    self.model_style_controller = ModelStyleController(parent=self,view=None)
    self.map_style_controller = MapStyleController(parent=self,view=None)
    if self.view.textInput.text() == 'Find automatically':
      self.view.textInput.setText(self.viewer.find_viewer())
    self.viewer.command = self.view.textInput.text()
    self.viewer.start_viewer(json_response=True)
    self.viewer.send_command("set bgColor white")

  # Parsing the output of a response is kept in the controller for now because it is messy
  def _interpret_response_text(self,response_text,debug=False):
    if debug:
      print("0: Response text:",type(response_text))
      print(response_text)
      print("1: Response json:",type(json.loads(response_text)))
      response_json = json.loads(response_text)
      print(response_json)
      print("2: 'json values' key:",type(response_json["json values"]))
      json_values = response_json["json values"]
      print("3: json value 0:",type(json_values[0]))
      print(json_values[0])
      if json_values[0] is None:
        return None
      print("4: json value 0 as json:",type(json.loads(json_values[0])))
      json_values_0 = json.loads(json_values[0])
      print(json_values_0)
      if not isinstance(json_values_0,list):
        json_values_0 = [json_values_0]
      print("5: json value 00:",type(json_values_0[0]))
      print(json_values_0[0])

    objects = [json.loads(e) for e in json.loads(response_text)["json values"]][0]
    return objects

  def _process_remote_ref_response_model(self,response):
    response_objects = self._interpret_response_text(response.text,debug=True)
    assert 'model specs' in response_objects, f"Unexpected response objects: {response_objects}"
    model_specs = response_objects['model specs']
    assert len(model_specs) ==1, f'Unsure how to interpret multiple model specs in response: {model_specs}'
    remote_ref_id = model_specs[0]
    return remote_ref_id

  # Models
  def _load_active_model(self,ref):
    assert ref is None or ref is self.state.active_model_ref
    if ref is not None and ref not in self.loaded_model_refs:
      self.load_model_from_ref(ref)

  def load_model_from_ref(self,ref,format='pdb',callback=None):
    if self.viewer._connected:
      if ref.data.filepath is not None and Path(ref.data.filepath).exists():
        if ref not in self.loaded_model_refs:
          response = self.viewer.load_model(
            filename=ref.data.filepath,
            format=format,
            ref_id=ref.id,
            label=ref.label,
            callback=callback
          )
        else:
          response = self.viewer.load_model_from_mmtbx(
            model=ref.model,
            format=format,
            ref_id=ref.id,
            label=ref.label,
            callback=callback
          )

        if response is not None and response.status_code == 200:
          remote_ref_id = self._process_remote_ref_response_model(response)
          #print("ADDING CHIMERAX REMOTE MODEL ID FOR MODEL: ",remote_ref_id)
          ref.external_ids["chimerax"] = remote_ref_id
          self.state.external_loaded["chimerax"].append(ref.id)
        else:
          print(f"Not adding remote ref from response: {response}")
        self.loaded_model_refs.append(ref)
        return response

  # Maps

  def _process_remote_ref_response_map(self,response):
    response_json = json.loads(response.text)
    note = response_json['log messages']['note'][0]
    pattern = r"as (.*?),"
    match = re.search(pattern, note)
    if match:
        found_text = match.group(1)  # The text between "as " and ","
        return found_text
    else:
        print("No match found for volume spec value")

  def _load_active_map(self,ref):
    assert ref is None or ref is self.state.active_map_ref
    if ref is not None and ref not in self.loaded_map_refs:
      self.load_map_from_ref(ref)



  def load_map_from_ref(self,ref,callback=None):
    if self.viewer._connected:
      if ref.data.filepath is not None and Path(ref.data.filepath).exists():
        if ref not in self.loaded_map_refs:
          response = self.viewer.load_map(
            filename=ref.data.filepath,
            ref_id_map=ref.id,
            ref_id_model = self.state.active_model_ref.id, # Should be stored on map_ref?
            label=ref.label,
            callback=callback
          )
        else:
          response = self.viewer.load_map_from_mmtbx(
            map_manager=ref.model,
            ref_id_map=ref.id,
            ref_id_model = self.state.active_model_ref.id,
            label=ref.label,
            callback=callback
          )

        if response is not None and response.status_code == 200:
          remote_ref_id = self._process_remote_ref_response_map(response)
          if remote_ref_id is not None:
            #print("ADDING CHIMERAX REMOTE MODEL ID FOR MAP: ",remote_ref_id)
            ref.external_ids["chimerax"] = remote_ref_id
        else:
          print(f"Not adding remote ref from response: {response}")
        self.loaded_map_refs.append(ref)
        self.viewer.set_transparency_map(ref.external_ids['chimerax'],0.7)
        return response

  # Selection

  def _process_remote_ref_response_selection(self,response,debug=False):
    response_text = response.text
    if debug:
      print("0: Response text:",type(response_text))
      print(response_text)
      print("1: Response json:",type(json.loads(response_text)))
    response_json = json.loads(response_text)
    if debug:
      print(response_json)
      print("2: 'json values' key:",type(response_json["json values"]))
    json_values = response_json["json values"]
    if debug:
      print("3: json value 0:",type(json_values[0]))
      print(json_values[0])
    if json_values[0] is None:
      return None
    if debug:
      print("4: json value 0 as json:",type(json.loads(json_values[0])))
    json_values_0 = json.loads(json_values[0])
    return json_values_0


  def poll_selection(self,callback=None):
    if not self.viewer._connected:
      return


    # get selected atoms
    if self._picking_granularity == "residue":
      self.viewer._select_up_residues()
    response = self.viewer.poll_selection()

    json_objects = self._process_remote_ref_response_selection(response)

    #print("chimerax controller poll_selection() json_objects:")
    if json_objects is None:
      return None


    # Translate chimerax response to per-atom dictionaries
    # Updated pattern to correctly handle the value after '@'
    pattern = r"(?:([^/]+))?/([^:]*)(?::([^@]*))?(?:@(.*))?"

    atoms = []
    for obj in json_objects:  # atoms
        spec = obj['spec']

        match = re.match(pattern, spec)
        if match:
          atom_dict = {}
          ref_id, asym_id, seq_id, atom_id = match.groups()
          if ref_id is None:
            ref_id = self.state.active_model_ref.id

          # print(spec)
          # print(f"Value before '/': {value0}")
          # print(f"Value after '/': {value1}")
          # print(f"Value after ':': {value2}")
          # print(f"Value after '@': {value3}")

          atom_dict["ref_id"] = ref_id
          atom_dict["asym_id"] = asym_id
          atom_dict["seq_id"] = seq_id
          atom_dict["atom_id"] = atom_id
          atoms.append(atom_dict)


        else:
            print("No match found")
            assert False, f"Must be able to translate atom spec... {spec}"


    #get df of all selected atoms
    df = pd.DataFrame.from_records(atoms)

    # group selections by reference
    grouped_dfs = {ref_id: group for ref_id, group in df.groupby('ref_id')}

    output = {}
    for ref_id,df in grouped_dfs.items():
      # build per-atom selections
      selections = []
      for key,df in grouped_dfs.items():
        records = df.to_dict("records")
        for record in records:
          #make selection
          d = copy.deepcopy(Selection.default_select_all)
          d['asym_id']["ops"][0]['value'] = record["asym_id"]
          d['seq_id']["ops"][0]['value'] = int(record["seq_id"])
          d['atom_id']["ops"][0]['value'] = record["atom_id"]
          selections.append(Selection(d,params=copy.deepcopy(Selection.default_params)))



      # make per-atom query
      query = SelectionQuery(selections=selections)

      # get the relevant mol
      model_ref = self.state.references[ref_id]

      # select sites
      sites_sel = model_ref.mol.atom_sites.select_from_query(query)

      # make condensed query
      query = model_ref.mol.atom_sites._convert_sites_to_query(sites_sel)
      query.params.refId = ref_id
      output[ref_id] = query

    if callback is not None:
      callback(output)
    return output
  # end poll_selection

  def select_from_phenix_string(self,selection_phenix: str):
    if self.viewer._connected:
      ref= self.state.active_model_ref
      model_id = ref.external_ids["chimerax"]
      return self.viewer.select_from_phenix_string(model_id,selection_phenix)

  def select_from_ref(self,ref):
    if self.viewer._connected:
      model_id = ref.model_ref.external_ids["chimerax"]
      phenix_string = ref.query.phenix_string
      self.viewer.select_from_phenix_string(model_id,phenix_string)

  def set_granularity(self,value="residue"):
    assert value in ['element','residue'], 'Provide one of the implemented picking levels'
    self._picking_granularity = value

  def toggle_selection_mode(self,value):
    pass

  def deselect_all(self):
    self.viewer.deselect_all()


  # Other

  def clear_viewer(self):
    self.viewer.clear_viewer()

  def close_viewer(self):
    if self.viewer._connected:
      self.viewer.close_viewer()


  def reset_camera(self):
    self.viewer.reset_camera()


  def sync_ref_mapping(self):
    pass



  # Style

  # Style
  def set_iso(self,ref,value):
    if self.viewer._connected:
      model_id = ref.external_ids["chimerax"]
      self.viewer.set_iso(model_id,value)


  def show_ref(self,ref,representation: Optional[str] = None):
    if self.viewer._connected:
      model_id = ref.model_ref.external_ids["chimerax"]
      if representation is None:
        representation = ref.style.representation
      else:
        representation = [representation]

      for rep_name in representation:
        # if isinstance(ref,(ModelRef,MapRef)):
        #   self.viewer.show_model(model_id)
        # else:
        #   self.viewer.show_query(model_id,ref.query.to_json(),rep_name)
        self.viewer.show_query(model_id,ref.query.to_json(),rep_name)

  def hide_ref(self,ref,representation: Optional[str] = None):
    if self.viewer._connected:
      model_id = ref.model_ref.external_ids["chimerax"]

      if representation is None:
        representation = ref.style.representation
      else:
        representation = [representation]

      for rep_name in representation:
        # if isinstance(ref,(ModelRef,MapRef)):
        #   self.viewer.hide_model(model_id)
        # else:
        #   self.viewer.hide_query(model_id,ref.query.to_json(),rep_name)
        self.viewer.hide_query(model_id,ref.query.to_json(),rep_name)

  def color_ref(self,ref,color):
    if self.viewer._connected:
      model_id = ref.model_ref.external_ids["chimerax"]
      self.viewer.color_query(model_id,ref.query.to_json(),color)
