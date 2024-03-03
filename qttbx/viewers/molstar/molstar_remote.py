from __future__ import absolute_import, division, print_function


import requests
from pathlib import Path
import sys
import time
import json
import tempfile
import glob
import os
import random
import requests
import threading
from typing import Optional
from requests.exceptions import RequestException

import subprocess
import sys
import time

try:
  from urllib.parse import urlencode
except ImportError:
  from urllib import urlencode

from libtbx.utils import Sorry



from PySide2.QtCore import QObject, Signal
from PySide2.QtCore import QUrl, Signal, QObject
from qttbx.viewers import ModelViewer
from .server_utils import find_open_port


# =============================================================================
class PhenixMolstarViewer(ModelViewer):

  viewer_name = 'Phenix Molstar '

  def __init__(self):
    super().__init__()

    # local state
    self.loaded_model_refs = []
    self.loaded_map_refs = []
    self.log_molstar = []
    self._port = None
    self._url = None

  @property
  def url(self):
    if self._url is None:
      return f"http://127.0.0.1:{self.port}/process_json"
    else:
      return self._url
  @url.setter
  def url(self,value):
    self._url = value
  @property
  def port(self):
    if self._port is None:
      # randomly select port
      self._port = find_open_port()
    return self._port
  @port.setter
  def port(self,value):
    self._port = value
  # ---------------------------------------------------------------------------
  # Start API
  # ---------------------------------------------------------------------------


  # ---------------------------------------------------------------------------
  # Status

  def is_available(self):
    '''
    Function for determining if Phenix Molstar  is available

    Parameters
    ----------
      Nothing

    Returns
    -------
      True if Phenix Molstar  is available
    '''
    self.find_viewer()
    if self.command:
      return True
    else:
      return False

  def find_viewer(self):
    '''
    Function for finding Phenix Molstar

    Parameters
    ----------
      Nothing

    Returns
    -------
      Command for running Phenix Molstar
    '''

    self.command = str(Path(__file__).parent / Path("../../command_line/start_viewer.py"))
    return self.command

  def start_viewer(self, timeout=60, show_tabs=["all"]):
    '''
    Function for starting the Phenix Molstar  REST server

    Parameters
    ----------
      timeout: int
        The number of seconds to wait for the REST server to become
        available before raising a Sorry

      json_response: bool
        Pass json true to Phenix Molstar  remotecontrol. Returns the log as
        a json object for inspection.

    Returns
    -------
      Nothing
    '''

    # append some standard search locations for each platform
    # the existing $PATH is searched first
    if self.command is None:
      self.find_viewer()


    self.run_basic_checks()

    # construct command
    cmd = (["python", self.command] +
            [f"show_tab={tab_name}" for tab_name in show_tabs]+
            ["rest_server_port={}".format(self.port)]
    )

    # clean environment for launching
    env = os.environ.copy()
    # for v in ['PYTHONPATH', 'LD_LIBRARY_PATH', 'DYLD_LIBRARY_PATH', 'DYLD_FALLBACK_LIBRARY_PATH']:
    #   env.pop(v, None)

    #  Start the subprocess detached from the parent process
    proccess = subprocess.Popen(cmd, start_new_session=True)



    # # Launch a thread to handle stdout
    # self.stdout_thread = threading.Thread(target=read_output, args=(self.process.stdout,))
    # self.stdout_thread.start()



    print()
    print('-'*79)
    print('Starting Phenix Molstar REST server')
    print(" ".join(cmd))
    print(self.url)
    counter = 0
    while counter<timeout:
      output = self._check_status()
      if self._connected:
        break
      counter += 1
      time.sleep(1)
    if not self._connected:
      raise Sorry('The Phenix Molstar  REST server is not reachable at {} after '
                  '{} seconds.'.format(self.url, counter))
    print('Phenix Molstar  is ready')
    print('-'*79)
    print()

  # ---------------------------------------------------------------------------
  def _check_status(self):
    '''
    Check if the REST server is available
    '''
    output = None
    try:
      data = {"status_ok":None}
      data = json.dumps(data)
      headers = {'Content-Type': 'application/json'}
      output = requests.post(self.url, headers=headers, json=data)
      if output.status_code == 200:
        self._connected = True
    except requests.exceptions.ConnectionError:
      self._connected = False
    return output

  def close_viewer(self):
    print('='*79)
    if self._connected:
      print('Shutting down Phenix Molstar ')
      params = {'command': 'close_viewer'}
      self._run_command(params)
    else:
      print('Phenix Molstar  already shut down')
    rc = self.process.returncode
    stdout, stderr = self.process.communicate()
    # if self.process:
    #   self.process.wait()
    #self.stdout_thread.join()
    # print('-'*79)
    # print(stdout)
    # print('-'*79)
    # print(stderr)
    print('='*79)

  def log_message(self,message=None):
    command = {
      "command":"log_message",
      "kwargs":{
        "message":message
      }
    }
    self.send_command(command)

  # ---------------------------------------------------------------------------
  # Remote communication

  def _run_command(self, data):
    print("DEBUG: _run_command() in molstar remote")

    try:
      # Template data
      #
      # data = {"command":"load_model_from_string",
      #         "kwargs":{
      #             "model_str":debug_model_str,
      #             "label":"mylabel",
      #             "format":"pdb",
      #           }
      #         }
      data = json.dumps(data,indent=2)
      print("Data being sent: ")
      print(data)
      response = requests.post(self.url, json=data)
      return response

    except requests.exceptions.ConnectionError:
      return None


  def send_command(self, data):
    print("DEBUG: send_command() in molstar remote")
    if self._connected == True:
      return self._run_command(data)


  # ---------------------------------------------------------------------------
  # Models

  def load_model(self, filename=None,format='pdb',label=None):
    """
    Load a model directly from file.
    """
    print("Debug: load_model() in molstar remote")
    filename = str(filename)

    command = {
      "command":"load_model",
      "kwargs":{
        "filename":filename,
        "format":format,
        "label":label,
      }
    }
    print("command to send:")
    self.send_command(command)


  def load_model_from_mmtbx(self,model,format='pdb',label=None):
    # This function does connect to identical one in MolstarViewer, because need
    # string serialization
    assert format in ['pdb','mmcif'], "Use one of the supported format names"
    if format == 'pdb':
      model_str = model.model_as_pdb()
    elif format == "mmcif":
      model_str = model.model_as_mmcif()

    return self.load_model_from_string(model_str,format=format,label=label)

  def load_model_from_string(self, model_str,format='pdb',label=None,callback=None):
    command = {
      "command":"load_model_from_string",
      "kwargs":{
        "model_str":model_str,
        "format":format,
        "label":label,
      }
    }
    self.send_command(command)

  # ---------------------------------------------------------------------------
  # Maps

  def load_map(self,filename=None,volume_id=None,model_id=None,label=None):
    """
    Load a map from disk.
    """
    filename = str(filename)
    command = {
      "command":"load_map",
      "kwargs":{
        "filename":filename,
        "volume_id":volume_id,
        "model_id":model_id,
        'label':label,
      }
    }
    self.send_command(command)


  def load_map_from_map_manager(self,map_manager=None,volume_id=None,model_id=None):
    with tempfile.NamedTemporaryFile(mode='w+t', suffix=".mrc",delete=False) as temp_file:
      map_manager.write_map(temp_file.name)
      return self.load_map(filename=temp_file.name,volume_id=volume_id,model_id=model_id)


  # ---------------------------------------------------------------------------
  # Selection


  def select_from_query(self,query_json):
    command = f"result = await {self.plugin_prefix}.phenix.select("+query_json+");"
    self.send_command(command,sync=True)


  def poll_selection(self,callback=None,queue=False):
    #print("Callback manager callback: ",self.command_queue.callback_manager.callback)
    # get selection
    command = f"""
    {self.plugin_prefix}.phenix.pollSelection();
    """
    self.send_command(command,callback=callback,queue=queue,wrap_async=False)


  def deselect_all(self,queue=False):
    command = f"{self.plugin_prefix}.phenix.deselectAll()"
    self.send_command(command,queue=queue)


  # ---------------------------------------------------------------------------
  # Other

  def clear_viewer(self,queue=False):
    # Remove all objects from the viewer
    command = f"{self.plugin_prefix}.plugin.clear()"
    self.send_command(command,queue=queue)

  def reset_camera(self,queue=False):
    command = f"{self.plugin_prefix}.plugin.managers.camera.reset();"
    self.send_command(command,queue=queue)

  def _toggle_selection_mode(self,value):
    if value == True:
      value = 'true'
    else:
      value = 'false'
    command = f"""
    {self.plugin_prefix}.phenix.toggleSelectionMode({value});
    """
    self.send_command(command)

  def _set_granularity(self,value="residue"):
    assert value in ['element','residue'], 'Provide one of the implemented picking levels'
    if value == "element":
      command = f"{self.plugin_prefix}.plugin.managers.interactivity.setProps({{ granularity: 'element' }})"
    elif value == "residue":
      command = f"{self.plugin_prefix}.plugin.managers.interactivity.setProps({{ granularity: 'residue' }})"
    self.send_command(command)

  def _get_sync_state(self,callback=None,verbose=False):
    # get the remote: local -> reference mapping from the web app
    command = f"""
    {self.plugin_prefix}.phenix.getState();
    """
    if verbose:
      assert callback is None, "Cannot use custom callback and verbose together"
      def callback(x):
        print(json.dumps(json.loads(x),indent=2))
    return self.send_command(command,callback=callback,wrap_async=False,queue=False,sync=True)


  def _set_sync_state(self,state_json):
    command = f"""
    {self.plugin_prefix}.phenix.setState('{state_json}')
    """
    self.send_command(command)

  # ---------------------------------------------------------------------------
  # Style


  # Volume ISO
  def set_iso(self,volume_id,value):
    # volume_id here is local id
    params = '{"entry":{"name":"'+volume_id+'","params":{"view":{"name":"camera-target","params":{"radius":5,"selectionDetailLevel":0,"isSelection":false,"bottomLeft":[6.000999927520752,6.004000186920166,6.011000156402588],"topRight":[30.89900016784668,11.321000099182129,18.44099998474121]}},"detailLevel":0,"channels":{"em":{"isoValue":{"kind":"absolute","absoluteValue":'+str(value)+'},"color":6524815,"wireframe":false,"opacity":0.3}}}}}'
    command = f"""
    const volumeEntry = {self.plugin_prefix}.phenix.getVolumeEntry('{volume_id}');
    volumeEntry.source.params.isoValue.absoluteValue = {value};
    const volumeStreamingRef = {self.plugin_prefix}.refMapping_volume['{volume_id}']
    await {self.plugin_prefix}.plugin.build().to(volumeStreamingRef).update({params}).commit();
    """
    self.send_command(command)

  # Transparency Volume
  def set_transparency_map(self,model_id,value):
    raise NotImplementedError

  # Show/Hide model/selection, works by setting transparency

  def show_model(self,model_id,representation_name: Optional[str] = None):
    raise NotImplementedError

  def show_selected(self,representation_name: Optional[str] = None):
    raise NotImplementedError

  def show_query(self,model_id: str, query_json: str, representation_name: str):
    self.transparency_query(model_id=model_id, query_json=query_json, representation_name=representation_name,value=0.0)


  def hide_model(self,model_id,representation_name: Optional[str] = None):
    raise NotImplementedError

  def hide_selected(self,representation_name: Optional[str] = None):
    raise NotImplementedError

  def hide_query(self,model_id: str, query_json: str, representation_name: str):
    self.transparency_query(model_id=model_id, query_json=query_json, representation_name=representation_name,value=1.0)


  # Transparency for model/selection

  def transparency_model(self,model_id: str, representation_name: Optional[str] = None):
    raise NotImplementedError

  def transparency_selected(self,representation_name: Optional[str] = None):
    raise NotImplementedError

  def transparency_query(self,model_id: str, query_json: str, representation_name: str, value: float):
    command = f"""
    const query = {self.plugin_prefix}.phenix.queryFromJSON('{query_json}');
    query.params.refId = '{model_id}'
    await {self.plugin_prefix}.phenix.setTransparencyFromQuery(query, '{representation_name}', '{value}')
    """
    self.send_command(command)


  # Color for model/selection

  def color_model(self,model_id: str, color: str):
    raise NotImplementedError

  def color_selected(self, color: str):
    raise NotImplementedError

  def color_query(self,model_id: str, query_json: str, color: str):
    command = f"""
    const query = {self.plugin_prefix}.phenix.queryFromJSON('{query_json}');
    query.params.refId = '{model_id}'
    await {self.plugin_prefix}.phenix.setQueryColor(query,'{color}')
    this.viewer.phenix.deselectAll();
    """
    self.send_command(command)

  # Representation for model/selection

  def representation_model(self,model_id: str, representation_name: str):
    query = SelectionQuery.from_model_ref()

  def representation_selected(self, representation_name: str):
    raise NotImplementedError


  def representation_query(self,model_id: str, query_json: str, representation_name: str):
    assert representation_name in ['ball-and-stick','ribbon']

    # add representation
    command = f"""
    const query = {self.plugin_prefix}.phenix.queryFromJSON('{query_json}');
    query.params.refId = '{model_id}'
    await {self.plugin_prefix}.phenix.addRepr(query,'{representation_name}')
    """
    self.send_command(command)

  def _get_representation_names(self,model_id: str, query_json: str):
    # add representation
    command = f"""
    const query = {self.plugin_prefix}.phenix.queryFromJSON('{query_json}');
    query.params.refId = '{model_id}'
    {self.plugin_prefix}.phenix.getRepresentationNames(query)
    """
    result = self.send_command(command)
