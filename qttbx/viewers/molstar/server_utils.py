"""
Utilities for running the molstar web server
"""
import http.server
import socketserver
import socket
import subprocess

from PySide2.QtCore import QThread,Signal, Slot


class NodeHttpServer:
  def __init__(self,command,default_port=8080):
    assert isinstance(command,list), "Provide command as a list of strings"
    if not self.check_port_free(default_port):
      default_port = self.find_open_port()
    self.port = default_port
    self.url = f"http://localhost:{self.port}"
    self.process = None
    self.command_list = command+['--port',str(self.port)]
    self.command = ' '.join(self.command_list)
    self.debug = True

  def log(self,*args):
    if self.debug:
      print(*args)
  @staticmethod
  def find_open_port():
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
      s.bind(("", 0))
      s.listen(1)
      port = s.getsockname()[1]
    return port


  @staticmethod
  def check_port_free(port, ip='localhost'):
    try:
      with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind((ip, port))
        # If we get here, it means the bind was successful,
        # indicating the port is free.
        return True
    except OSError:
      # If an OSError is caught, it likely means the port is in use.
      return False


  def start(self):
    if self.process is None:
      print(f"Starting HTTP server at: {self.url}")
      print(f"Command used: {self.command}")
      print("Command list: ",self.command_list)

      self.process = subprocess.Popen(self.command_list,stdout=None,stderr=None)
    else:
      print("HTTP server is already running.")

  def stop(self):
    if self.process:
      print("Stopping HTTP server...")
      self.process.terminate()
      self.process = None
    else:
      print("HTTP server is not running.")


class CustomHandler(http.server.SimpleHTTPRequestHandler):
  """
  This class allows custom directory path for server root
  """
  allow_reuse_address = True
  def __init__(self, *args, **kwargs):
    super().__init__(*args, directory=__file__ , **kwargs)


class HttpServerThread(QThread):
  stop_signal = Signal()

  def __init__(self,ip='localhost',default_port=8080,check_free=False):
    self.ip = ip
    if check_free:
      if not self.check_port_free(default_port,ip=ip):
        new_port = self.find_open_port()
        print(f"Default port {default_port} not free, chose {new_port}")
        default_port = new_port
    self.port = default_port
    self.url = f"http://{self.ip}:{self.port}"
    super(HttpServerThread, self).__init__()
    self.stop_signal.connect(self.stop)
    self.is_running = True



  def find_open_port(self):
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
      s.bind(("", 0))
      s.listen(1)
      port = s.getsockname()[1]
    return port

  def check_port_free(self,port,ip='localhost'):
    try:
      # Create a new socket using the AF_INET address family (IPv4)
      # and the SOCK_STREAM socket type (TCP)
      s = socket.create_connection((ip, port), timeout=1)
      s.close()
      return True
    except ConnectionRefusedError:
      return False

  def run(self):
    with socketserver.TCPServer((self.ip, self.port), CustomHandler) as httpd:
      httpd.timeout = 1  # timeout in seconds

      print(f"HTTP Server running at: {self.url}")
      while self.is_running:
        httpd.handle_request()  # This should now be non-blocking because of the timeout
      print("HTTP Server stopped.")

  @Slot()
  def stop(self):
    print("Stopping HTTP Server...")
    self.is_running = False
