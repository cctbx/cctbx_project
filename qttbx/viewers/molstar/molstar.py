

from pathlib import Path
import time
import json
from typing import Optional

import requests

from PySide2.QtCore import QUrl
try:
  from qttbx.viewers import ModelViewer
except:
  from viewers import ModelViewer

from libtbx.utils import Sorry
from libtbx import group_args

from .volume_streaming import VolumeStreamingManager
from ..gui.state.state import PhenixState
from ..gui.state.color import Color
from .server_utils import  NodeHttpServer
from .volume_streaming import VolumeStreamingManager
from ..core.selection import Selection
from ..core.python_utils import DotDict



# =============================================================================

class MolstarViewer(ModelViewer):
  """
  The Python interface for the molstar viewer. The specific Molstar plugin written to pair
  with this class. It functions as a controller for the 'viewer' tab in the GUI.
  """
  viewer_name = 'Molstar'
  def __init__(self,web_view,use_web_view=True,config_json_file=None):
    if config_json_file is None:
      config_json_file = Path(__file__).parent / Path("config.json")
    ModelViewer.__init__(self)

    self.web_view = web_view
    self.state = group_args(phenixState=None)
    #self.command_queue = CommandQueue()
    #self.references_remote_map={} # Keys: remote (molstar) ref_ids, Values: local ref_ids


    self.use_web_view = use_web_view
    self.selenium_driver = None
    self.log_list = []




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
    print('Molstar is ready')
    #self.send_command(f"{self.plugin_prefix}.postInit()")
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

  def log_message(self,message):
    self.log_list.append(message)

  # ---------------------------------------------------------------------------
  # Remote communication

  def send_command(self, js_command,callback=print,sync=False,log_js=True):
    # queue needs to be refactored out

    if log_js:
      print("JavaScript command:")
      lines = js_command.split("\n")
      print("Total Lines: ",len(lines))
      if len(lines)>40:
        lines = lines[:20]+["","...",""]+ lines[-20:]
        js_command_print = "\n".join(lines)
      else:
        js_command_print = js_command
      print(js_command_print)
      print("Callback:",callback)
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


  # def send_command(self, cmds,callback=None,queue=False,wrap_async=True,sync=False):
  #   """
  #   This function takes a list of commands or a single command.
  #   Each command is a string of Javascript
  #   code that will be executed.
  #   """
  #   # print("="*79)
  #   # print(f"Molstar plugin command: (connected={self._connected})")
  #   # print()
  #   # if len(cmds)==1:
  #   #   print(cmds[0])
  #   # else:
  #   #   print(cmds)
  #   # print()
  #   # print(f"wrap_async={wrap_async}")
  #   # print("="*79)
  #   out = None
  #   if self._connected:
  #     if self._blocking_commands:
  #       print(f"Commands ignored due to not accepting commands: {cmds}")
  #       return False
  #     if not isinstance(cmds,list):
  #       assert isinstance(cmds,str), "send_command takes a string or list of strings"
  #       cmds = [cmds]
  #     if not isinstance(callback,list):
  #       callback= [callback for cmd in cmds]
  #     for cmd,callback in zip(cmds,callback):
  #       #self._run_command(cmd,callback=callback)
  #       self.command_queue.add(cmd,callback)
  #     if not queue:
  #       if self.use_web_view:
  #         web_view = self.web_view
  #       else:
  #         web_view = None
  #       out = self.command_queue.run(web_view,self.selenium_driver,wrap_async=wrap_async,sync=sync)
  #   return out
  # def execute_command_queue(self,wrap_async=True,sync=False):
  #   if self.use_web_view:
  #     web_view = self.web_view
  #   else:
  #     web_view = None
  #   self.command_queue.run(web_view,self.selenium_driver,wrap_async=wrap_async,sync=sync)

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


  def load_model_from_mmtbx(self,model,format='pdb',label=None,ref_id=None,callback=None):

    func_name = f"model_as_{format}"
    func = getattr(model,func_name)
    model_string = func()
    return self.load_model_from_string(model_string,format=format,label=label,ref_id=ref_id,callback=callback)


  def load_model_from_string(self, model_str,format='pdb',label=None,ref_id=None,callback=None):
    """
    Load a model using a raw string. ref_id is important
    """
    if ref_id is None:
      ref_id = 'null'
    if label is None:
      label = 'model'

    command =  self._load_model_build_js(model_str,format=format,label=label,ref_id=ref_id)
    self.send_command(command,sync=True)


  # don't use this, it is confusing, but it loads without a ref
  # def load_model(self, filename=None,format='pdb',label=None,ref_id=None,callback=None):
  #   """
  #   Load a model directly from file. Note if using upstream gui, it will not automatically be aware.
  #   """
  #   filename = str(filename)

  #   dm = DataManager()
  #   dm.process_model_file(filename)
  #   model = dm.get_model()

  #   format="pdb" # only one supported for now
  #   self.load_model_from_string(model.model_as_pdb(),format=format,label=label,ref_id=ref_id,callback=callback)


  # ---------------------------------------------------------------------------
  # Maps

  # def load_map_from_mmtbx(self,map_manager,filename=None,ref_id_map=None,ref_id_model=None,callback=None):
  #   assert filename is not None, "Provide map filename explicitly"
  #   self.volume_streamer.pack_volume_path(filename=filename,ref_id_map=ref_id_map)
  #   js_str = self._load_map_build_js(ref_id_map=ref_id_map,ref_id_model=ref_id_model)
  #   self.send_command(js_str,callback=callback)

  def _load_map_build_js(self,volume_id=None,model_id=None):
    assert volume_id is not None, "At this state, maps must already have a ref id in the volume streamer"
    assert model_id is not None, "Maps cannot be loaded without an accompanying model"
    url_server =self.volume_streamer.url

    js_str = f"""
    {self.plugin_prefix}.volumeServerURL = '{url_server}';
    await {self.plugin_prefix}.phenix.loadMap('{model_id}','{volume_id}');
    """
    return js_str

  def load_map(self,filename=None,volume_id=None,model_id=None):
    """
    Load a map from disk.
    """
    filename = str(filename)

    self.volume_streamer.pack_volume_path(volume_path=filename,volume_id=volume_id)
    js_str = self._load_map_build_js(volume_id=volume_id,model_id=model_id)
    self.send_command(js_str)


  # ---------------------------------------------------------------------------
  # Selection

  def select_from_selection(self,selection: Selection):
    self.deselect_all()
    molstar_syntax = selection.molstar_syntax
    command = f"""
    const MS = {self.plugin_prefix}.MS
    const sel = MS.struct.generator.atomGroups({{
              'atom-test': {molstar_syntax}}})
    {self.plugin_prefix}.phenix.selectFromSel(sel);
    """
    self.send_command(command)

  # def _select_from_json(self,selectionJSON):
  #   self.deselect_all()
  #   selection_dict = json.loads(selectionJSON)
  #   molstar_syntax = selection_dict["molstar_syntax"]
  #   command = f"""
  #   const MS = {self.plugin_prefix}.MS
  #   const.StructureSelectionQuery = {self.plugin_prefix}.StructureSelectionQuery
  #   const selectionQuery = StructureSelectionQuery('Custom Query',
  #           MS.struct.generator.atomGroups({{
  #             'atom-test': {molstar_syntax}}}),
  #       )
  #   {self.plugin_prefix}.phenix.selectFromQuery(selectionQuery);
  #   """
  #   self.send_command(command)
  #   # self.send_command(command)
  #   # command = f"""
  #   # const molstarSelection = '{}'
  #   # {self.plugin_prefix}.phenix.selectFromJSON('{selectionJSON}');
  #   # """

  def select_from_phenix_string(self,phenix_string):
    selection = Selection.from_phenix_string(phenix_string)
    self.select_from_selection(selection)

  def poll_selection(self,callback=None):
    # get selection
    command = f"""
    {self.plugin_prefix}.phenix.pollSelection();
    """
    result_str = self.send_command(command,callback=callback,sync=True)
    try:
      atom_records = json.loads(result_str)
      selection = Selection.from_atom_records(atom_records)
      return selection
    except:
      print(result_str)
      raise
      return None
  def focus_selected(self):

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
        print("Verbose sync ouput: ")#,x,", type: ",type(x))
        output =  _validate_sync_output(x)
        if output is not None:
          print(json.dumps(output,indent=2))

    #Run and get output
    output = self.send_command(command,callback=callback,sync=True,log_js=True)
    #print("Returned sync output:")
    #print(type(output))
    #print(output)
    output = _validate_sync_output(output)
    #print("Building phenix state")
    #print(type(output))
    #print(output)
    if output is not None:
      output = PhenixState.from_dict(output)
    return output


  def _set_sync_state(self,state_json):
    assert False
    # command = f"""
    # {self.plugin_prefix}.phenix.setState('{state_json}')
    # """
    # self.send_command(command)

  # ---------------------------------------------------------------------------
  # Style


  # Volume ISO
  def set_iso(self,volume_id,value):
    # volume_id here is local id
    params = '{"entry":{"name":"'+volume_id+'","params":{"view":{"name":"camera-target","params":{"radius":5,"selectionDetailLevel":0,"isSelection":false,"bottomLeft":[6.000999927520752,6.004000186920166,6.011000156402588],"topRight":[30.89900016784668,11.321000099182129,18.44099998474121]}},"detailLevel":0,"channels":{"em":{"isoValue":{"kind":"absolute","absoluteValue":'+str(value)+'},"color":6524815,"wireframe":false,"opacity":0.5}}}}}'
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

  # def show_query(self,model_id: str, query_json: str, representation_name: str):
  #   self.transparency_query(model_id=model_id,
  #                           query_json=query_json,
  #                           representation_name=representation_name,
  #                           component_key="all",
  #                           value=0.0)


  def hide_model(self,model_id,representation_name: Optional[str] = None):
    raise NotImplementedError

  def selection_hide(self,component_name: Optional[str] = None, representation_name: Optional[str] = None):
    self.transparency_selected(component_name=component_name, representation_name=representation_name,value=1.0)


  def selection_show(self,component_name: Optional[str] = None, representation_name: Optional[str] = None):
    self.transparency_selected(component_name=component_name, representation_name=representation_name,value=0.0)

  def representation_add(self,representation_name):
    command = f"""
    await {self.plugin_prefix}.phenix.addRepresentationSelected('{representation_name}')
    """
    self.send_command(command)

  # Transparency for model/selection

  def transparency_model(self,model_id: str, representation_name: Optional[str] = None):
    raise NotImplementedError

  
  def transparency_selected(self,component_name: Optional[str] = None, representation_name: Optional[str] = None, value: Optional[float] = None):
    if component_name is None:
      component_name = 'undefined'
    else:
      component_name = f"'{component_name}'"
    if representation_name is None:
      representation_name = 'undefined'
    else:
      representation_name = f"'{representation_name}'"
    if value is None:
      value = 0.0
    command = f"""
    await {self.plugin_prefix}.phenix.setTransparencySelected({component_name},{representation_name},'{value}')
    """
    self.send_command(command)

  def transparency_query(self,selection: str, model_id: str, representation_name: str,component_key: str, value: float):
    """
    Set transparency. The primary method to show/hide molecules.
    Params:
      model_id: the Python side identifier for what to operate on. (refId)
      selection: What to operate on. Is The json format of a Selection type.
      representation_name: Specify what representation to operate on (ball-and-stick, cartoon, etc)
      component_key: Spcify what component to operate on (structure-component-static-polymer)
      opacity: A float 0-1 (visible-invisible)
    """
    command = f"""
    const selection = {self.plugin_prefix}.phenix.selectionFromJSON('{query_json}');
    await {self.plugin_prefix}.phenix.setTranfsparencyFromQuery(selection, {model_id}, '{representation_name}', '{component_key}', '{value}')
    """
    self.send_command(command)


  # Color for model/selection


  def color_selected(self, color: str):
    color = Color(color)

    command = f"""
    await {self.plugin_prefix}.phenix.setColorSelected({color.R},{color.G},{color.B})
    """
    self.send_command(command)
    

  # def color_query(self,model_id: str, query_json: str, color: str):
  #   command = f"""
  #   const query = {self.plugin_prefix}.phenix.queryFromJSON('{query_json}');
  #   query.params.refId = '{model_id}'
  #   await {self.plugin_prefix}.phenix.setQueryColor(query,'{color}')
  #   this.viewer.phenix.deselectAll();
  #   """
  #   self.send_command(command)

  # Representation for model/selection

  # def representation_model(self,model_id: str, representation_name: str):
  #   query = SelectionQuery.from_model_ref()

  def representation_selected(self, representation_name: str):
    # add representation
    command = f"""

    await {self.plugin_prefix}.phenix.addRepresentationSelected('{representation_name}')
    """
    self.send_command(command)


  # def representation_query(self,model_id: str, query_json: str, representation_name: str):
  #   assert representation_name in ['ball-and-stick','ribbon']

  #   # add representation
  #   command = f"""
  #   const query = {self.plugin_prefix}.phenix.queryFromJSON('{query_json}');
  #   query.params.refId = '{model_id}'
  #   await {self.plugin_prefix}.phenix.addRepr(query,'{representation_name}')
  #   """
  #   self.send_command(command)

  def _get_representation_names(self,model_id: str):
    command = f"""
    {self.plugin_prefix}.phenix.getRepresentationNames('{model_id}')
    """
    return self.send_command(command)
