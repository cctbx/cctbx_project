"""
This file defines the API for communication with the Phenix implementation of the molstar web app
"""

from pathlib import Path
import time
import json
from typing import Optional

import requests
import urllib.parse


from PySide2.QtCore import QUrl
from PySide2.QtWebEngineWidgets import QWebEnginePage, QWebEngineSettings
from qttbx.viewers import ModelViewer

from libtbx.utils import Sorry
from libtbx import group_args

from ..gui.model.molstar import MolstarState
from .server_utils import  NodeHttpServer
from ..gui.model.selection import Selection


# =============================================================================

class MolstarGraphics(ModelViewer):
  """
  The Python interface for the molstar viewer. A specific Molstar plugin is written to pair
  with this class.
  """
  viewer_name = 'Molstar'
  def __init__(self,web_view,use_web_view=True,dm=None,config_json_file=None):
    if config_json_file is None:
      config_json_file = Path(__file__).parent / Path("config.json")
    ModelViewer.__init__(self)

    self.web_view = web_view
    self.state = group_args(molstarState=None)
    self.dm = dm


    self.use_web_view = use_web_view
    self.selenium_driver = None
    self.log_list = []
    self.debug = True
    self._initial_sync_done = False # Set to True the first time communication is established with js viewer




    # get config variables
    with open(config_json_file,"r") as fh:
      config = json.load(fh)
      self.config = group_args(**config)
      expected_config_params = [
        "node_js_path",
        "volume_server_relative_path",
        "pack_script_relative_path",
        "molstar_app_relative_path",
      ]
      for name in expected_config_params:
        assert name in config, f"Missing necessary config variable: '{name}'"

      for key,value in config.items():
        if key == 'node_js_path' and value == '':
          value = 'npm'
        else:
          value = str(Path(__file__).parent / Path(value))
        setattr(self.config,key, value)

    # Flags
    self._blocking_commands = False

  def log(self,*args):
    if self.debug:
      print(*args)

  # ---------------------------------------------------------------------------
  # Start API
  # ---------------------------------------------------------------------------


  # ---------------------------------------------------------------------------
  # Status
  def is_available(self):
    '''
    Function for determining if Molstar is available

    Parameters
    ----------
      Nothing

    Returns
    -------
      True if available
    '''
    self.find_viewer()
    if self.command:
      return True
    else:
      return False

  def find_viewer(self):
    '''
    Function for finding Molstar

    Parameters
    ----------
      Nothing

    Returns
    -------
      Command for running Molstar
    '''
    self.app_root_dir = Path(self.config.molstar_app_relative_path)

    self.node_js_path = self.config.node_js_path
    self.command = ['http-server',str(self.app_root_dir)]
    return self.command


  def start_viewer(self,plugin_prefix="this.viewer",volume_streaming=False,timeout=60):
    '''
    Function for starting Molstar. Sequence of events:
      1. Start web server for molstar app
      2. Start volume server for volume streaming

    Parameters
    ----------

    Returns
    -------
      Nothing
    '''

    self.plugin_prefix = plugin_prefix


    # Start volume streaming server TODO: This should only occur if volumes will be used
    # if volume_streaming:
    #   self.log()
    #   self.log('-'*79)
    #   self.log('Starting volume streaming server')
    #   self.volume_streamer = VolumeStreamingManager(
    #             node_js_path = self.config.node_js_path,
    #             volume_server_relative_path = self.config.volume_server_relative_path,
    #             pack_script_relative_path = self.config.pack_script_relative_path,
    #             default_server_port=1336,
    #             debug = True
    #   )
    #   self.volume_streamer.start_server()
    #   self.log(self.volume_streamer.url)
    #   self.log('-'*79)
    #   self.log()


    # Start node http-server
    self.log()
    self.log('-'*79)
    self.log('Starting HTTP server for Molstar')
    self.command = self.find_viewer()
    self.view_server = NodeHttpServer(self.command)
    self.view_server.start()
    self.command = self.view_server.command
    self.port = self.view_server.port
    self.url = self.view_server.url


    # Set url on web view
    self.web_view.setUrl(QUrl(self.url))


    counter = 0
    while counter<timeout:
      self._check_status()
      if self._connected:
        break
      counter += 1
      time.sleep(1)
    if not self._connected:
      raise Sorry('The Molstar on the QT web view is not reachable at {} after '
                  '{} seconds.'.format(self.url, counter))
    self.log('Molstar is ready')
    #self.send_command(f"{self.plugin_prefix}.postInit()")
    self.log('-'*79)
    self.log()


  def _check_status(self):
    '''
    Check if the server is available
    '''
    output = None
    try:
      output = requests.get(url=self.url)
      if output.status_code == 200:
        self._connected = True
    except requests.exceptions.ConnectionError:
      self._connected = False
    return self._connected



  def close_viewer(self):
    self.view_server.stop()
    if hasattr(self,"volume_streamer"):
      self.volume_streamer.stop_server()
    self.log('='*79)

  def log_message(self,message):
    self.log_list.append(message)

  # ---------------------------------------------------------------------------
  # Remote communication

  def send_command(self, js_command,callback=print,sync=False,log_js=True):
    if log_js:
      self.log("JavaScript command:")
      lines = js_command.split("\n")
      self.log("Total Lines: ",len(lines))
      if len(lines)>40:
        lines = lines[:20]+["","...",""]+ lines[-20:]
        js_command_print = "\n".join(lines)
      else:
        js_command_print = js_command
      self.log(js_command_print)
      self.log("Callback:",callback)
    if not sync:
      js_command= f"""
      (async () => {{
          {js_command}
      }})();
      """
      result =  self.web_view.runJavaScript(js_command,custom_callback=callback)
    else:
      result =  self.web_view.runJavaScriptSync(js_command,custom_callback=callback)


    return result


  # ---------------------------------------------------------------------------
  # Models

  def _load_model_build_js(self,model_str,format='pdb',label=None,ref_id=None):
    assert ref_id is not None, 'Cannot load into molstar without defining a Python identifier (ref_id)'
    assert label is not None, 'Cannot load into molstar without label'
    js_str = f"""
    var model_str = `{model_str}`
    {self.plugin_prefix}.phenix.loadStructureFromPdbString(model_str,'{format}', '{label}', '{ref_id}')
    """
    return js_str


  def load_model(self,filename=None,ref_id=None):
    model = self.dm.get_model(filename=filename)
    format='pdb'
    label="model"
    callback=None

    func_name = f"model_as_{format}"
    func = getattr(model,func_name)
    model_str = func()

    command =  self._load_model_build_js(model_str,format=format,label=label,ref_id=ref_id)
    self.send_command(command,sync=False)


  # ---------------------------------------------------------------------------
  # Selection

  def _set_query_string(self,selection: Selection):
    # Returns js string to set 'query' from a Selection object
    molstar_syntax = selection.molstar_syntax
    js_str = f"""
    const MS = {self.plugin_prefix}.MS
    const sel = MS.struct.generator.atomGroups({{
              'atom-test': {molstar_syntax}}})
    const query = {self.plugin_prefix}.StructureSelectionQuery('Phenix Query',sel)
    """
    return js_str
    

  def select_from_selection(self,selection: Selection):
    """
    Make a selection from a Selection object. All other selection functions lead here
    """
    self.deselect_all()
    command = f"""
    {self._set_query_string(selection)}
    {self.plugin_prefix}.phenix.selectFromQuery(query);
    """
    self.send_command(command)

  
  # def select_from_selection_string(self,phenix_string):
  #   """
  #   Make a selection from a Phenix selection string. Is converted to molstar
  #     syntax via a Selection instance
  #   """
  #   selection = Selection.from_selection_string(phenix_string)
  #   self.select_from_selection(selection)

  def poll_selection(self,callback=None):
    """
    Get the current selection as a new selection object.
    """
    command = f"""
    {self.plugin_prefix}.phenix.pollSelection();
    """
    result_str = self.send_command(command,callback=callback,sync=True)
    try:
      atom_records = json.loads(result_str)
      selection = Selection.from_atom_records(atom_records)
      return selection
    except:
      self.log(result_str)
      raise
      return None

  def focus_selected(self):
    """
    Focus on the selected region
    """
    command = f"""
    {self.plugin_prefix}.phenix.focusSelected();
    """
    self.send_command(command,sync=True)

  def select_all(self):
    command = f"{self.plugin_prefix}.phenix.selectAll()"
    self.send_command(command)

  def deselect_all(self):
    command = f"{self.plugin_prefix}.phenix.deselectAll()"
    self.send_command(command)


  # ---------------------------------------------------------------------------
  # Other

  def clear_viewer(self):
    # Remove all objects from the viewer
    command = f"{self.plugin_prefix}.plugin.clear()"
    self.send_command(command)

  def reset_camera(self):
    command = f"{self.plugin_prefix}.plugin.managers.camera.reset();"
    self.send_command(command)

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


  # ---------------------------------------------------------------------------
  # Synchronization

  def sync_remote(self,callback=None,verbose=True):
    # get the remote: local -> reference mapping from the web app
    command = f"""
    {self.plugin_prefix}.phenix.getState();
    """
    # Inline sync check function
    def _validate_sync_output(sync_output):
      return_val = None
      if isinstance(sync_output,(str,dict)):
        try:
          output = json.loads(sync_output)
          #assert isinstance(output,dict)
          return_val = output
        except:
          #raise
          pass
      return return_val

    if verbose:
      assert callback is None, "Cannot use custom callback and verbose together"

      # Inline print callback
      def callback(x):
        self.log("Verbose sync ouput: ")
        output =  _validate_sync_output(x)
        if output is not None:
          self.log(json.dumps(output,indent=2))

    #Run and get output
    output = self.send_command(command,callback=callback,sync=True,log_js=True)
    self.log("Returned sync output:")
    self.log(type(output))
    self.log(output)
    output = _validate_sync_output(output)
    self.log("Building phenix state")
    self.log(type(output))
    self.log(output)
    if output is not None:
      output = MolstarState.from_dict(output)
      if output.has_synced:
        self._initial_sync_done = True
      elif not output.has_synced and not self._initial_sync_done:
        self._set_sync_state()

    
    return output


  def _set_sync_state(self):
    command = f"""
    {self.plugin_prefix}.phenix.setState()
    """
    self.send_command(command)

  