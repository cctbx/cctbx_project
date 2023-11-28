from pathlib import Path
import socket
import tempfile
import subprocess
import json


class VolumeStreamingManager:
  """
  The order of operations to display a volume in molstar:
  
  1. Initiate a temporary directory for .mdb files 
        (database files for the volume server)
  2. Start the volume server configured with the temp directory
        (molstar/lib/commonjs/servers/volume/server.js)
  3. Generate an .mdb file in the temp directory
  4. Subsequent volumes just need an .mdb created in the right place

  This class is to manage all of this from Python. 
  

  NOTE: For now, 'em' is hardcoded as a map id type
  """
  def __init__(self,
               server_js_path="/Users/user/Desktop/CambridgeTrip/Molstar/molstar/lib/commonjs/servers/volume/server.js",
               pack_js_path = "/Users/user/Desktop/CambridgeTrip/Molstar/molstar/lib/commonjs/servers/volume/pack.js",
               node_js_path = "/opt/homebrew/bin/node",
               default_server_port=1336,
               debug = True
              ):
    self.server_process = None
    self.node_js_path = Path(node_js_path).absolute()
    self.temp_dir = tempfile.TemporaryDirectory()
    self.server_js_path = Path(server_js_path).absolute()
    self.pack_js_path = Path(pack_js_path).absolute()
    
    # get any free port of default populated
    if not self.check_port_free(default_server_port):
      default_server_port = self.find_open_port()
    self.server_port = default_server_port
    
    
    self.data = {} # keys are data_manager keys (dm.get_real_map_names())
                   # values are {"mdb":<filepath to mdb file>,
                   #             "map":<filepath to map>,
                   #             "source":<volume server source,"em" or "xray">,
                   #             "id": <volume server id, Path(mdb).stem
    self.debug = debug
    if debug:
      print(self.mdb_path)


  def __str__(self):
    # data
    def json_serializable(obj):
      return str(obj)
    data = json.dumps(self.data,default=json_serializable,indent=2)
    pid = None
    running = self.server_process is not None
    if running: 
      pid = self.server_process.pid
    s = f'Volume Server:\n\tURL: {self.url}\n\tRunning: {running}\n\tPID: {pid}\n\tdata: {data}'
    return s


  def check_port_free(self,port,ip='localhost'):
    try:
      # Create a new socket using the AF_INET address family (IPv4)
      # and the SOCK_STREAM socket type (TCP)
      s = socket.create_connection((ip, port), timeout=1)
      s.close()
      return True
    except ConnectionRefusedError:
      return False
    
  def find_open_port(self):
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
      s.bind(("", 0))
      s.listen(1)
      port = s.getsockname()[1]
    return port
    
  @property
  def mdb_path(self):
    return Path(self.temp_dir.name).absolute()

  @property
  def url(self):
    return f"http://localhost:{self.server_port}/VolumeServer"

  # def pack_data_manager(self,dm,state):
  #   """
  #   Pack a data manager into .mdb files. This goes through
  #   a map manager intermediate
  #   """
  #   names = dm.get_real_map_names()
  #   for name in names:
      
  #     mm = dm.get_real_map(filename=name)
  #     self.pack_map_manager(mm,filename=None,key=name)
  
  # def pack_map_manager(self,map_manager,filename=None,ref_id_map=None):
  #   """
  #   Pack a map manager into a .mdb file. A filename to the map
  #   must be specified or present on mm.file_name.

  #   Key is how this map_manager will be accessed in this class, ideally
  #   the same as the key in the datamanager, which is also a filename.
  #   """
  #   key = ref_id_map 
  #   if map_manager.file_name is not None:
  #     filename = map_manager.file_name
  #   else:
  #     assert filename is not None, "Error: Provide a file name for this map manager."
    
  #   if key is None:
  #     key = Path(filename).stem
  #   self.pack_volume_path(filename,key=key)

  def pack_volume_path(self,volume_path,volume_id):
    """
    Run the VolumeServer pack.js on a volume file.
    Store the resulting .mdb in a temporary directory
    Update self.data which keeps track of all packed volumes
    """
    volume_path = Path(volume_path).absolute()
    volume_name = volume_path.stem

      
    mdb_path = Path(self.mdb_path,"em",f"{volume_id}.mdb")
    command = [str(self.node_js_path),
               str(self.pack_js_path),
               "em",
                str(volume_path),
                str(mdb_path)]
    print('Pack command: ',command)
    result = subprocess.run(
      command,
      env={'TEMP_DIR': self.mdb_path},
      stdout=subprocess.PIPE,
      stderr=subprocess.PIPE,
      text=True
    )
    if self.debug:
      print("stdout:", result.stdout)
      print("stderr:", result.stderr)
      
    self.data[volume_id] = {
      "path_mdb":mdb_path,
      "path_map":volume_path,
      "source":"em",
      'label':volume_name,
      "id":volume_id
    }

  
  @property
  def available_volumes(self):
    files = list(self.mdb_path.glob("**/*"))
    return files
  


  def start_server(self):
    if self.server_process:
      print("Volume Server is already running.")
    command = [str(self.node_js_path),
               str(self.server_js_path),
               "--idMap",
               "em",
               str(self.mdb_path) + "/em/${id}.mdb",
               "--defaultPort",
               str(self.server_port)]
    
    print("Starting Volume Server with command: "," ".join(command))
    self.server_process = subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    print(f"Volume Server started with PID: {self.server_process.pid} on port: {self.server_port}")

  def stop_server(self):
    print("Cleanup temp dir...")

    self.temp_dir.cleanup()
    if not self.server_process:
      print("Volume Server is not running.")
      return
  
    self.server_process.terminate()
    self.server_process.wait()
    print(f"Volume Server with PID {self.server_process.pid} terminated.")
    self.server_process = None