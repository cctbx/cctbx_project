import requests
from pathlib import Path
import time
import json
from functools import partial
from typing import Optional, List

from PySide2.QtCore import QObject, Signal
from PySide2.QtCore import QUrl, Signal, QObject, QTimer
try:
  from qttbx.viewers import ModelViewer
except:
  from viewers import ModelViewer

try:
  from selenium import webdriver
  from webdriver_manager.chrome import ChromeDriverManager
  from selenium.webdriver.chrome.service import Service
except:
  pass
from iotbx.data_manager import DataManager
from libtbx.utils import Sorry


from .volume_streaming import VolumeStreamingManager
from ..gui.controller.controller import Controller
from ..gui.controller.selection_controls import SelectionControlsController
from .server_utils import  NodeHttpServer
from .volume_streaming import VolumeStreamingManager
from ..last.selection_utils import SelectionQuery
from ..gui.controller.style import ModelStyleController, MapStyleController
from ..last.python_utils import DotDict

class CallbackManager:
  """
  Intermediate layer for callback functions to prevent undesired
  repeated calls from subsequent events emitting the same signal.
  """
  def __init__(self):
    super().__init__()
    self.callback = None
    self.cmd = None

  def call(self, *args, **kwargs):
  
    if self.callback is not None:
      print(f"Callback manager calling: {self.callback}")
      func = self.callback
      self.callback = None
      self.cmd = None
      func(*args, **kwargs)


  def add_callback(self,func,cmd=None):
    if func is not None:
      assert self.callback is None, f"A second callback was added ({func}) before the first was executed ({self.callback})"
      self.callback = func
      self.cmd = cmd

class CommandQueue(QObject):
  """
  Enable batching commands. This is crucial to avoid async bugs, where the order of javascript commands given
  is not necessarily the order they finish.
  """
  commandRequested = Signal(str, object,object,object,bool) # (cmd,web_view,selenium_driver,callback,sync)
  #jsResultReady = Signal(str) # signal that the webview js command return value is ready
  def __init__(self):
    super().__init__()
    self.commands = [] # list of commands
    self.callback = None
    self.callback_manager = CallbackManager()

    
    self.commandRequested.connect(self._execute_command)



  def add(self,cmd,callback=None):
    if callback is not None:
      if self.callback is None:
        self.callback = callback
      else:
        assert self.callback == callback, "Cannot form a command queue with multiple callbacks"
    self.commands.append(cmd)
  
  def run(self,web_view=None,selenium_driver=None,wrap_async=True,sync=False):
    command = "\n".join(self.commands)
    if wrap_async:
      command= f"""
      (async () => {{
          {command}
      }})();
      """
    self.commandRequested.emit(command,web_view,selenium_driver,self.callback,sync)
    self.clear()


  def clear(self):
    self.commands = []
    self.callback = None

 
  def _execute_command(self,cmd,web_view,selenium_driver,callback=None,sync=False):
    """
    Actually execute a command in a web view. Requires a QT Web view to be
    set up, and this function conencted to self.emitter.commandRequested
    """
    #print(f"Executing command at time: {time.time()}")
    #print(cmd)

    self.callback_manager.add_callback(callback)

    # web view
    if web_view is not None:
      if not sync:
        web_view.runJavaScript(cmd,self._handle_command_result)
      else:
        web_view.runJavaScriptSync(cmd,self._handle_command_result)

    # selenium
    if selenium_driver is not None:
      result = selenium_driver.execute_script(cmd)
  
      print(f"Result going to selenium callback: {result}")
      if callback is not None:
        callback(result)

  def _handle_command_result(self,result):
    """
    The standing callback for all javascript run by _execute_command. It's only
    function is to emit a signal that a Javascript result is ready. The callback manager
    should be the one that listens to this signal, and runs actually usefull callbacks.
    """
    #(f"Handling callback for result at time: {time.time()}")
    #print(result)
    #print("Deliberately not running callback")
    self.callback_manager.call(result)
    #print(self.callback_manager.cmd)
    #self.jsResultReady.emit(result) # Emit signal to run function specific callback after result is ready


# =============================================================================
class MolstarViewer(ModelViewer):
  """
  The Python interface for the molstar viewer. The specific Molstar plugin written to pair
  with this class. It functions as a controller for the 'viewer' tab in the GUI.
  """
  viewer_name = 'Molstar'
  def __init__(self,web_view,use_web_view=True,use_selenium=False,config_json_file=None):
    if config_json_file is None:
      config_json_file = Path(__file__).parent / Path("config.json")
    ModelViewer.__init__(self)
    
    self.web_view = web_view
    self.command_queue = CommandQueue()
    #self.references_remote_map={} # Keys: remote (molstar) ref_ids, Values: local ref_ids

    # Use selenium to get interactive app in chrome
    self.use_selenium = use_selenium
    self.use_web_view = use_web_view
    self.selenium_driver = None
    
    
    # get config variables
    with open(config_json_file,"r") as fh:
      config = json.load(fh)
      self.config = DotDict()
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
        self.config[key] = value

      if self.config["node_js_path"] == "":
        self.config["node_js_path"] = 'npm'

    # Flags
    self._blocking_commands = False


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


  def start_viewer(self,plugin_prefix="this.viewer",timeout=60):
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

    # Start the HTTP server using built in server
    # self.http_server_thread = HttpServerThread()
    # self.http_server_thread.start()
    # self.ip = self.http_server_thread.ip
    # self.port = self.http_server_thread.port
    # self.url = self.http_server_thread.url
    self.plugin_prefix = plugin_prefix


    # Start volume streaming server TODO: This should only occur if volumes will be used
    print()
    print('-'*79)
    print('Starting volume streaming server')
    self.volume_streamer = VolumeStreamingManager(
              node_js_path = self.config.node_js_path,
              volume_server_relative_path = self.config.volume_server_relative_path,
              pack_script_relative_path = self.config.pack_script_relative_path,
              default_server_port=1336,
              debug = True
    )
    self.volume_streamer.start_server()
    print(self.volume_streamer.url)
    print('-'*79)
    print()


    # Start node http-server
    print()
    print('-'*79)
    print('Starting HTTP server for Molstar')
    self.command = self.find_viewer()
    self.view_server = NodeHttpServer(self.command)
    self.view_server.start()
    self.command = self.view_server.command
    self.port = self.view_server.port
    self.url = self.view_server.url
    # Set url on web view
    self.web_view.setUrl(QUrl(self.url))
    if self.use_selenium:
      #service = Service(executable_path='/Users/user/software/phenix/modules/cctbx_project/qttbx/chromedriver-mac-arm64')
      #self.selenium_driver = webdriver.Chrome(service=service)
      self.selenium_driver = webdriver.Chrome()
      #
      #self.selenium_driver.set_script_timeout(30)
      #self.selenium_driver = webdriver.Chrome(ChromeDriverManager().install())

      self.selenium_driver.get(self.url)


      
    counter = 0
    while counter<timeout:
      output = self._check_status()
      if self._connected:
        break
      counter += 1
      time.sleep(1)
    if not self._connected:
      raise Sorry('The Molstar on the QT web view is not reachable at {} after '
                  '{} seconds.'.format(self.url, counter))
    print('Molstar is ready')
    self.send_command(f"{self.plugin_prefix}.postInit()")
    print('-'*79)
    print()




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
    self.volume_streamer.stop_server()
    print('='*79)

  # ---------------------------------------------------------------------------
  # Remote communication

  def send_command(self, cmds,callback=None,queue=False,wrap_async=True,sync=False):
    """
    This function takes a list of commands or a single command. 
    Each command is a string of Javascript 
    code that will be executed. 
    """
    print("="*79)
    print(f"Molstar plugin command: (connected={self._connected})")
    print()
    if len(cmds)==1:
      print(cmds[0])
    else:
      print(cmds)
    print()
    print(f"wrap_async={wrap_async}")
    print("="*79)
    out = "fish"
    if self._connected:
      if self._blocking_commands:
        print(f"Commands ignored due to not accepting commands: {cmds}")
        return False
      if not isinstance(cmds,list):
        assert isinstance(cmds,str), "send_command takes a string or list of strings"
        cmds = [cmds]
      if not isinstance(callback,list):
        callback= [callback for cmd in cmds]
      for cmd,callback in zip(cmds,callback):
        #self._run_command(cmd,callback=callback)
        self.command_queue.add(cmd,callback)
      if not queue:
        if self.use_web_view:
          web_view = self.web_view
        else:
          web_view = None
        out = self.command_queue.run(web_view,self.selenium_driver,wrap_async=wrap_async,sync=sync)
    return out
  def execute_command_queue(self,wrap_async=True,sync=False):
    if self.use_web_view:
      web_view = self.web_view
    else:
      web_view = None
    self.command_queue.run(web_view,self.selenium_driver,wrap_async=wrap_async,sync=sync)
  
  # ---------------------------------------------------------------------------
  # Models

  def _load_model_build_js(self,model_str,format='pdb',label=None,ref_id=None):
    assert ref_id is not None, 'Cannot load into molstar without preexisting Python ref'
    assert label is not None, 'Cannot load into molstar without label'
    js_str = f"""
    var model_str = `{model_str}`
    {self.plugin_prefix}.phenix.loadStructureFromPdbString(model_str,'{format}', '{label}', '{ref_id}')
    """
    return js_str


  def _load_model_from_mmtbx(self,model,format='pdb',label=None,ref_id=None,callback=None):

    func_name = f"model_as_{format}"
    func = getattr(model,func_name)
    model_string = func()
    return self._load_model_from_string(model_string,format=format,label=label,ref_id=ref_id,callback=callback)


  def _load_model_from_string(self, model_str,format='pdb',label=None,ref_id=None,callback=None,queue=False):
    """
    Load a model using a raw string. ref_id is important
    """
    if ref_id is None:
      ref_id = 'null'
    if label is None:
      label = 'model'
  
    command =  self._load_model_build_js(model_str,format=format,label=label,ref_id=ref_id)
    self.send_command(command,callback=callback,queue=queue)


  def load_model(self, filename=None,format='pdb',label=None,ref_id=None,callback=None):
    """
    Load a model directly from file. Note if using upstream gui, it will not automatically be aware.
    """
    filename = str(filename)
    
    with open(filename,"r") as fh:
      contents = fh.read()
    self.load_model_from_string(contents,format=format,label=label,ref_id=ref_id,callback=callback)


  # ---------------------------------------------------------------------------
  # Maps

  # def load_map_from_mmtbx(self,map_manager,filename=None,ref_id_map=None,ref_id_model=None,callback=None):
  #   assert filename is not None, "Provide map filename explicitly"
  #   self.volume_streamer.pack_volume_path(filename=filename,ref_id_map=ref_id_map)
  #   js_str = self._load_map_build_js(ref_id_map=ref_id_map,ref_id_model=ref_id_model)
  #   self.send_command(js_str,callback=callback)

  def _load_map_build_js(self,volume_id=None,model_id=None):
    iso = 0.5 # TODO: sane default
    assert volume_id is not None, "At this state, maps must already have a ref id in the volume streamer"
    assert model_id is not None, "Maps cannot be loaded without an accompanying model"
    url_server =self.volume_streamer.url
    
    js_str = f"""
    {self.plugin_prefix}.volumeServerURL = '{url_server}';
    await {self.plugin_prefix}.phenix.loadMap('{model_id}','{volume_id}');
    """
    return js_str
  
  def load_map(self,filename,volume_id,model_id):
    """
    Load a map from disk.
    """
    filename = str(filename)
    
    self.volume_streamer.pack_volume_path(volume_path=filename,volume_id=volume_id)
    js_str = self._load_map_build_js(volume_id=volume_id,model_id=model_id)
    self.send_command(js_str)


  # ---------------------------------------------------------------------------
  # Selection


  def select_from_query(self,query_json):
    command = f"result = await {self.plugin_prefix}.phenix.select("+query_json+");"
    self.send_command(command,sync=True)


  def poll_selection(self,callback=None,queue=False):
    print("Callback manager callback: ",self.command_queue.callback_manager.callback)
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

  def _get_sync_state(self,state_json,callback=None):
    # get the remote: local -> reference mapping from the web app
    command = f"""
    {self.plugin_prefix}.phenix.getState('{state_json}');
    """
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
    print("RESULT")
    print(result)