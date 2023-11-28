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
try:
  from qttbx.viewers import ModelViewer
except:
  from viewers import ModelViewer



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
    
    print("="*79)
    print(f"ChimeraX command: (connected={self._connected})")
    print()
    if len(cmds)==1:
      print(cmds[0])
    else:
      print(cmds)
    print()
    print("="*79)

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

  def load_model_from_string(self, model_str,format='pdb',label=None,ref_id=None,callback=None):

    with tempfile.NamedTemporaryFile(mode='w+t', suffix="."+format,delete=False) as temp_file:
      temp_file.write(model_str)
    
    response = self.load_model(filename=temp_file.name,label=label,ref_id=ref_id,callback=callback)
    os.remove(temp_file.name)
    return response

  def load_model_from_mmtbx(self,model,format='pdb',label=None,ref_id=None,callback=None):
    assert format in ['pdb','mmcif'], "Use one of the supported format names"
    if format == 'pdb':
      model_str = model.model_as_pdb()
    elif format == "mmcif":
      model_str = model.model_as_mmcif()

    return self.load_model_from_string(model_str,format=format,label=label,ref_id=ref_id,callback=callback)

  # ---------------------------------------------------------------------------
  # Maps

  def load_map_from_mmtbx(self,map_manager,filename=None,label=None,ref_id_map=None,ref_id_model=None,callback=None):

    with tempfile.NamedTemporaryFile(mode='wb', suffix=".mrc",delete=False) as temp_file:
      temp_filename = temp_file.name
    
    map_manager.write_map(temp_filename)
    response = self.load_map(filename=temp_filename)
    os.remove(temp_filename)
    return response

  def load_map(self, filename=None,label=None,ref_id_map=None,ref_id_model=None,callback=None):
    """
    Load a map from disk.
    """
    filename = str(filename)
    return self._load_file(filename=filename)
    


  # ---------------------------------------------------------------------------
  # Selection

  def interpret_response_text(self,response_text,debug=True):
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

  def poll_selection(self):
    response = self.send_command('info sel')
    return response

  def _select_up_residues(self):
    response = self.send_command("select sel residues true")


  def deselect_all(self):
    response = self.send_command('~select')

  # Other
  def clear_viewer(self):
    # Remove all objects from the viewer
    response = self.send_command('close')

  def close_viewer(self):
    print('='*79)
    if self._connected:
      print('Shutting down ChimeraX')
      params = {'command': 'quit'}
      self._run_command(params)
    else:
      print('ChimeraX already shut down')
    if self.process is not None:
      rc = self.process.returncode
      stdout, stderr = self.process.communicate()
    # print('-'*79)
    # print(stdout)
    # print('-'*79)
    # print(stderr)
    print('='*79)


  def reset_camera(self,queue=False):
    response = self.send_command("view")


  # ---------------------------------------------------------------------------
  # Style

  def set_iso(self,ref_id,value):
    # ref_id is remote ref_id here, the ChimeraX spec, like '#1'
    command = f"volume {ref_id} level {value}"
    self.send_command(command)

  def set_transparency(self,model_id,value):
    command = f"volume {model_id} transparency  {value}"
    self.send_command(command)



  def hide_model(self,remote_ref_id):
    command = f"hide {remote_ref_id} models"
    self.send_command(command)

  def show_model(self,remote_ref_id):
    command = f"show {remote_ref_id} models"
    self.send_command(command)

  def show_selection(self):
    command = f"show sel atoms"
    self.send_command(command)

  def hide_selection(self):
    command = f"hide sel atoms"
    self.send_command(command)

  def color_selection(self,value):
    command = f"color sel {value}"
    self.send_command(command)

  def color_model(self,ref_id,value):
    command = f"color {ref_id} {value}"
    self.send_command(command)
  
  def show_ribbon_model(self,ref_id):
    command = f"show {ref_id} cartoons"
    self.send_command(command)

  def show_ribbon_selection(self):
    command = f"show sel cartoons"
    self.send_command(command)

  def hide_ribbon_model(self,ref_id):
    command = f"hide {ref_id} cartoons"
    self.send_command(command)

  def hide_ribbon_selection(self):
    command = f"hide sel cartoons"
    self.send_command(command)

  # def select_from_phenix(self,selection_string):
  #   # convert to json
  #   query = self.state.mol.sites._select_query_from_str_phenix(selection_string)
  #   self.select_from_query(query)

  # def select_from_json(self,query_json):
  #   # json representation of SelectionQuery class
  #   query = SelectionQuery.from_json(query_json)
  #   self.select_from_query(query)

  # def select_from_query(self,query):
  #   if query.params.refId in ["",None]:
  #     query.params.refId = self.state.active_model_ref.id
  #   query_json = query.to_json()
  #   command = f"result = await {self.plugin_prefix}.select("+query_json+");"
  #   self.send_command(command)

  # def poll_selection(self,callback=None,queue=False):
  #   # return the selection on the viewer side, 
  #   # which may not be the active selection stored in self.state
  #   #print("(#1) molstar.poll_selection")
  #   command = f"{self.plugin_prefix}.pollSelection()"
  #   self.send_command(command,callback=callback,queue=queue,wrap_async=False)

  # # def clear_selection(self,queue=False):
  # #   # This version removes coloring. Use deselect_all to retain color.
  # #   command = f"{self.plugin_prefix}.clearSelection()"
  # #   self.send_command(command,queue=queue)

  # def deselect_all(self,queue=False):
  #   command = f"{self.plugin_prefix}.visual.deselectAll()"
  #   self.send_command(command,queue=queue)

  # # def highlight_selection(self):
  # #   # Highlighting usually comes with selection
  # #   raise NotImplementedError

  # # def focus_selection(self):
  # #   raise NotImplementedError

  # # def select_and_focus(self,selection_string):
  # #   # Select, highlight, and focus
  # #   # Implementing this first because that is 
  # #   # the existing behavior
  # #   print("select_and_focus()")


  # # ---------------------------------------------------------------------------
  # # Other

  # def clear_viewer(self,queue=False):
  #   # Remove all objects from the viewer
  #   self.loaded_map_refs = []
  #   self.loaded_model_refs = []
  #   command = f"{self.plugin_prefix}.plugin.clear()"
  #   self.send_command(command,queue=queue)
  #   for ref in self.state.references.values():
  #     if ref.entry is not None:
  #       ref.entry.is_active = False

   
  # def reset_camera(self,queue=False):
  #   command = f"{self.plugin_prefix}.plugin.managers.camera.reset();"
  #   self.send_command(command,queue=queue)

  # def toggle_selection_mode(self,value):
  #   if value == True:
  #     value = 'true'
  #   else:
  #     value = 'false'
  #   command = f"""
  #   {self.plugin_prefix}.toggleSelectionMode({value});
  #   """
  #   self.send_command(command)

  # def set_granularity(self,value="residue"):
  #   assert value in ['element','residue'], 'Provide one of the implemented picking levels'
  #   if value == "element":
  #     command = f"{self.plugin_prefix}.plugin.managers.interactivity.setProps({{ granularity: 'element' }})"
  #   elif value == "residue":
  #     command = f"{self.plugin_prefix}.plugin.managers.interactivity.setProps({{ granularity: 'residue' }})"
  #   self.send_command(command)

  # # ---------------------------------------------------------------------------
  # # Style

  # def set_iso(self,ref,value):
  #   #print(f"settings iso for ref: {ref.id} to {value}")
  #   command = f"""
  #   {self.plugin_prefix}.volumeRefInfo.params.values.entries[0].source.params.isoValue.absoluteValue = {value};
  #   {self.plugin_prefix}.plugin.build().to({self.plugin_prefix}.volumeStreamingRef).update().commit();
  #   """
  #   self.send_command(command)

  # def set_color(self,ref,color,theme='uniform',queue=False):
  #   assert theme in ['uniform'], "Provide predefined color theme"
  #   if theme == 'uniform':
  #     command = f"""
  #     var query = {self.plugin_prefix}.queryFromString('{ref.query.to_json()}');
  #     {self.plugin_prefix}.visual.colorSelection(query, '{color}');
  #     """
  #     self.send_command(command,queue=queue)


  # def set_representation(self,ref,representation_name,queue=False):
  #   # select
  #   self.select_from_json(ref.query.to_json())
  #   # hide
  #   self.hide_selection()

  #   # add representation
  #   command = f"""
  #   var query = {self.plugin_prefix}.queryFromString('{ref.query.to_json()}');
  #   await {self.plugin_prefix}.visual.addRepr(query,'{representation_name}')
  #   """
  #   self.send_command(command,queue=queue)

  # def set_visibility(self,ref_id,is_visible,queue=False):
  #   ref = self.state.references[ref_id]
  #   # select
  #   self.select_from_json(ref.query.to_json())
  #   if not is_visible:
  #     # hide
  #     self.hide_selection(queue=queue)
  #   else:
  #     # show
  #     self.show_selection(queue=queue)
    

  # def hide_selection(self,queue=False):
  #   command = f"{self.plugin_prefix}.hideSelection();"
  #   self.send_command(command,queue=queue)
  #   #self.state.emitter.signal_repr_change.emit(ref_id,[])

  # def show_selection(self,queue=False):
  #   command = f"{self.plugin_prefix}.showSelection();"
  #   self.send_command(command,queue=queue)
