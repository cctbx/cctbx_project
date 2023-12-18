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


from ..last.selection_utils import SelectionQuery
from .sel_convert_chimera import (
  translate_phenix_selection_string,
  convert_selection_str_to_int,
  convert_selection_int_to_str
)




# =============================================================================
class ChimeraXViewer(ModelViewer):

  viewer_name = 'ChimeraX'
  
  def __init__(self):
    super().__init__()

    # local state
    self.loaded_model_refs = []
    self.loaded_map_refs = []


  # ---------------------------------------------------------------------------
  # Start API
  # ---------------------------------------------------------------------------


  # ---------------------------------------------------------------------------
  # Status

  def is_available(self):
    '''
    Function for determining if ChimeraX is available

    Parameters
    ----------
      Nothing

    Returns
    -------
      True if ChimeraX is available
    '''
    self.find_viewer()
    if self.command:
      return True
    else:
      return False

  def find_viewer(self):
    '''
    Function for finding ChimeraX

    Parameters
    ----------
      Nothing

    Returns
    -------
      Command for running ChimeraX
    '''

    if self.command is not None:
      return self.command  # no work to do

    search_paths = os.getenv('PATH')
    search_command = self.viewer_name
    if search_paths is not None:
      if sys.platform == 'win32':
        search_paths = search_paths.split(';')
      else:
        search_paths = search_paths.split(':')

    # /Applications/ChimeraX<version>.app
    if sys.platform == 'darwin':
      known_paths = glob.glob('/Applications/ChimeraX*.app/Contents/MacOS')
      known_paths.sort(reverse=True)
      search_paths += known_paths
    # /usr/bin/chimerax
    elif sys.platform.startswith('linux'):
      search_paths += ['/usr/bin']
      search_command = 'chimerax'
    # C:\Program Files\ChimeraX <version>\bin
    elif sys.platform == 'win32':
      known_paths = glob.glob('C:\\Program Files\\ChimeraX*\\bin')
      known_paths.sort(reverse=True)
      search_paths += known_paths
      search_command = 'ChimeraX.exe'

    if len(search_paths) > 0:
      self.command = self.find_command(cmd=search_command, path=search_paths)

    return self.command

  def start_viewer(self, timeout=60,json_response=False):
    '''
    Function for starting the ChimeraX REST server

    Parameters
    ----------
      timeout: int
        The number of seconds to wait for the REST server to become
        available before raising a Sorry

      json_response: bool
        Pass json true to ChimeraX remotecontrol. Returns the log as 
        a json object for inspection.

    Returns
    -------
      Nothing
    '''

    # append some standard search locations for each platform
    # the existing $PATH is searched first
    if self.command is None:
      self.find_viewer()

    # randomly select port
    if self.port is None:
      self.port = random.randint(49152, 65535)

    # set REST server information
    self.flags = ['--cmd', 'remotecontrol rest start port {}'.format(self.port)]
    
    if json_response:
      self.flags[-1]+=" json true"

    self.url = "http://127.0.0.1:{}/".format(self.port)

    self.run_basic_checks()

    # construct ChimeraX command
    # ChimeraX --cmd "remotecontrol rest start port <port>"
    cmd = [self.command] + self.flags

    # clean environment for launching ChimeraX
    env = os.environ.copy()
    for v in ['PYTHONPATH', 'LD_LIBRARY_PATH', 'DYLD_LIBRARY_PATH', 'DYLD_FALLBACK_LIBRARY_PATH']:
      env.pop(v, None)

    # start ChimeraX server and wait until it is ready
    self.process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env)

    print()
    print('-'*79)
    print('Starting ChimeraX REST server')
    print(self.command)
    print(self.url)
    counter = 0
    while counter<timeout:
      output = self._check_status()
      if self._connected:
        break
      counter += 1
      time.sleep(1)
    if not self._connected:
      raise Sorry('The ChimeraX REST server is not reachable at {} after '
                  '{} seconds.'.format(self.url, counter))
    print('ChimeraX is ready')
    print('-'*79)
    print()

  # ---------------------------------------------------------------------------
  def _check_status(self):
    '''
    Check if the REST server is available
    '''
    output = None
    try:
      output = requests.get(url=self.url + 'cmdline.html')
      if output.status_code == 200:
        self._connected = True
    except requests.exceptions.ConnectionError:
      self._connected = False
    return output

  def close_viewer(self):
    print('='*79)
    if self._connected:
      print('Shutting down ChimeraX')
      params = {'command': 'quit'}
      self._run_command(params)
    else:
      print('ChimeraX already shut down')
    rc = self.process.returncode
    stdout, stderr = self.process.communicate()
    # print('-'*79)
    # print(stdout)
    # print('-'*79)
    # print(stderr)
    print('='*79)

  # ---------------------------------------------------------------------------
  # Remote communication

  def _run_command(self, params):
    '''
    Make requests call to REST server
    https://www.cgl.ucsf.edu/chimerax/docs/user/commands/remotecontrol.html
    '''
    try:
      # Encode query parameters using urllib.parse.urlencode
      encoded_params = urlencode(params)

      # Append the encoded query parameters to the URL
      url = self.url+"run"
      url_with_params = url+"?"+encoded_params

      # Send a POST request with form data as query parameters in the URL
      response = requests.post(url_with_params)

      return response

    except requests.exceptions.ConnectionError:
      return None


  def send_command(self, cmds):
    
    # print("="*79)
    # print(f"ChimeraX command: (connected={self._connected})")
    # print()
    # if len(cmds)==1:
    #   print(cmds[0])
    # else:
    #   print(cmds)
    # print()
    # print("="*79)

    if self._connected == True:
      if not isinstance(cmds,list):
          cmds = [cmds]
      params = [("command",c) for c in cmds]
      return self._run_command(params)

  
  # ---------------------------------------------------------------------------
  # Models

  def _load_file(self, filename=None):
    command = f'open {filename}'
    return self.send_command(command)

  def load_model(self, filename=None,format=None,label=None,ref_id=None,callback=None):
    return self._load_file(filename)

  def _load_model_from_string(self, model_str,format='pdb',label=None,ref_id=None,callback=None):

    with tempfile.NamedTemporaryFile(mode='w+t', suffix="."+format,delete=False) as temp_file:
      temp_file.write(model_str)
    
    response = self.load_model(filename=temp_file.name,label=label,ref_id=ref_id,callback=callback)
    os.remove(temp_file.name)
    return response

  def _load_model_from_mmtbx(self,model,format='pdb',label=None,ref_id=None,callback=None):
    assert format in ['pdb','mmcif'], "Use one of the supported format names"
    if format == 'pdb':
      model_str = model.model_as_pdb()
    elif format == "mmcif":
      model_str = model.model_as_mmcif()

    return self._load_model_from_string(model_str,format=format,label=label,ref_id=ref_id,callback=callback)

  # ---------------------------------------------------------------------------
  # Maps

  # def load_map_from_mmtbx(self,map_manager,filename=None,label=None,ref_id_map=None,ref_id_model=None,callback=None):

  #   with tempfile.NamedTemporaryFile(mode='wb', suffix=".mrc",delete=False) as temp_file:
  #     temp_filename = temp_file.name
    
  #   map_manager.write_map(temp_filename)
  #   response = self.load_map(filename=temp_filename)
  #   os.remove(temp_filename)
  #   return response

  def load_map(self, filename=None,label=None,ref_id_map=None,ref_id_model=None,callback=None):
    """
    Load a map from disk.
    """
    filename = str(filename)
    return self._load_file(filename=filename)
    


  # ---------------------------------------------------------------------------
  # Selection

  def poll_selection(self):
    response = self.send_command('info sel')
    return response

  def _select_up_residues(self):
    response = self.send_command("select sel residues true")


  def deselect_all(self):
    response = self.send_command('~select')
  
  def select_from_query(self,model_id: str, query_json: str):
    """
    This takes a json string, converts to query object, and then translates. 
    The idea is to maintain a looser coupling with upstream by never using 
    query objects in the api.
    """
    query = SelectionQuery.from_json(query_json)
    self.select_from_phenix_string(model_id,query.phenix_string)

  def select_from_phenix_string(self,model_id: str, phenix_selection: str, focus=True):
    # focus by default for consistency with molstar
    #print(f"DEBUG: select_from_phenix_string({model_id},{phenix_selection})")
    model_number = int(model_id.replace("#",""))
    # if phenix_selection in ["","all"]:
    #   chimerax_selection = model_id
    # else:
    chimerax_selection = translate_phenix_selection_string(phenix_selection,model_number=model_number)
    print("DEBUG: chimerax_selection: ",chimerax_selection)
    command = f"sel {chimerax_selection}"
    self.send_command(command)

    if focus:
      command = f"view {chimerax_selection}"
      self.send_command(command)



  # Other
  def clear_viewer(self):
    # Remove all objects from the viewer
    response = self.send_command('close')

  def reset_camera(self,queue=False):
    response = self.send_command("view")


  # ---------------------------------------------------------------------------
  # Style

  #  Volume ISO
  def set_iso(self,volume_id,value):
    # ref_id is remote ref_id here, the ChimeraX spec, like '#1'
    command = f"volume {volume_id} level {value}"
    self.send_command(command)

  # Volume Transparency
  def set_transparency_map(self,model_id,value):
    command = f"volume {model_id} transparency  {value}"
    self.send_command(command)

# Show/Hide model/selection, works by setting transparency

  def _get_level_from_representation(self,representation_name):
    # Convert from upstream/molstar representation names and chimerax level keyword
    if representation_name is None:
      level = 'models'
    elif representation_name == "ball-and-stick":
      level = 'atoms'
    elif representation_name == 'cartoon':
      level = 'cartoons'
    else:
      assert False, f"Unrecognized representation name: {representation_name}"
    return level

  def show_model(self,model_id,representation_name: Optional[str] = None):
    level = self._get_level_from_representation(representation_name)
    command = f"show {model_id} {level}"
    self.send_command(command)

  def show_selected(self,representation_name: Optional[str] = None):
    level = self._get_level_from_representation(representation_name)
    command = f"show sel {level}"
    self.send_command(command)

  def show_query(self,model_id: str, query_json: str, representation_name: Optional[str] = None ):
    self.deselect_all()
    self.select_from_query(model_id,query_json)
    self.show_selected(representation_name=representation_name)
    self.deselect_all()

  def hide_model(self,model_id,representation_name: Optional[str] = None):
    level = self._get_level_from_representation(representation_name)
    command = f"hide {model_id} {level}"
    self.send_command(command)

  def hide_selected(self,representation_name: Optional[str] = None):
    level = self._get_level_from_representation(representation_name)
    command = f"hide sel {level}"
    self.send_command(command)

  def hide_query(self,model_id: str, query_json: str,representation_name: Optional[str] = None ):
    self.deselect_all()
    self.select_from_query(model_id,query_json)
    self.hide_selected(representation_name=representation_name)
    self.deselect_all()


  # Transparency for model/selection

  def transparency_model(self,model_id: str, representation_name: Optional[str] = None):
    raise NotImplementedError

  def transparency_selected(self,representation_name: Optional[str] = None):
    raise NotImplementedError

  def transparency_query(self,model_id: str, query_json: str, representation_name: str, value: float):
    raise NotImplementedError


  # Color for model/selection

  def color_model(self,model_id: str,color: str):
    command = f"color {model_id} {color}"
    self.send_command(command)

  def color_selected(self,color: str):
    command = f"color sel {color}"
    self.send_command(command)

  def color_query(self,model_id: str, query_json: str, color: str):
    self.deselect_all()
    self.select_from_query(model_id,query_json)
    self.color_selected(color) # currently works on all reps
    self.deselect_all()

  # Representation for model/selection

  def representation_model(self,model_id: str, representation_name: str):
    self.show_model(model_id,representation_name)

  def representation_selected(self,representation_name: str):
    self.show_selected(representation_name)

  def representation_query(self,model_id: str, query_json: str, representation_name: str):
    self.show_query(model_id,query_json,representation_name)

    # add representation
    command = f"""
    const query = {self.plugin_prefix}.queryFromJSON('{query_json}');
    query.params.refId = '{model_id}'
    await {self.plugin_prefix}.visual.addRepr(query,'{representation_name}')
    """
    self.send_command(command)

  def _get_representation_names(self,model_id: str, query_json: str):
    # add representation
    command = f"""
    const query = {self.plugin_prefix}.queryFromJSON('{query_json}');
    query.params.refId = '{model_id}'
    {self.plugin_prefix}.visual.getRepresentationNames(query)
    """
    result = self.send_command(command)
    print("RESULT")
    print(result)

  # ---------------------------------------------------------------------------
  # functions for wxPython GUI
  def is_alive(self):
    self._check_status()
    return self._connected

  def quit(self):
    return self.close_viewer()

# =============================================================================

