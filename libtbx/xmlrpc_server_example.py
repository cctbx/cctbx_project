from __future__ import absolute_import, division, print_function

# This is an example of how a 3rd-party program with Python embedded, such
# as Coot or PyMOL, can be interfaced with CCTBX-based software.  Something
# much like this is used for the Phenix GUI extensions to those programs.
# I haven't tried this with any other software, but anything with a reasonably
# recent version of Python and support for either persistent Python threads
# or some sort of timer callback should be able to use it.

DEFAULT_PORT = 40000

import os, sys, string, signal
import xmlrpclib

try :
  from SimpleXMLRPCServer import SimpleXMLRPCServer
  class external_xmlrpc_server(SimpleXMLRPCServer):
    def __init__(self, addr, cctbx_interface):
      self.cctbx_interface = cctbx_interface
      SimpleXMLRPCServer.__init__(self, addr, logRequests=0)

    def _dispatch(self, method, params):
      if not self.cctbx_interface.enable_xmlrpc :
        return -1
      result = -1
      func = getattr(self.cctbx_interface, method, None)
      if not callable(func):
        print("%s is not a callable object!" % method)
      else :
        result = func(*params)
        if result is None :
          result = -1
      return result

  class external_xmlrpc_interface(object):
    def __init__(self, program_id, auto_start=True, verbose=False):
      self.enable_xmlrpc = True
      self.xmlrpc_server = None
      self.cctbx_server = None
      self.verbose = verbose
      self.timeout = string.atoi(os.environ.get("CCTBX_XMLRPC_TIMEOUT", "250"))
      self.program_id = program_id
      self.supported_modules = []
      self.setup_modules()
      self.setup_server()
      if auto_start :
        self.start_server()

    def setup_modules(self):
      pass

    def add_module(self, module_object=None, module_path=None):
      if module_object is not None :
        self.supported_modules.append(module_object)
      elif module_path is not None :
        module_object = __import__(module_path)
        self.supported_modules.append(module_object)

    def setup_server(self):
      port = os.environ.get("CCTBX_%s_PORT" % self.program_id, DEFAULT_PORT)
      if port is not None :
        self.port = int(port)
        self.xmlrpc_server = external_xmlrpc_server(("127.0.0.1", self.port),
                                                    self)
        if self.verbose :
          print("Listening on port %s" % port)
      cctbx_port = os.environ.get("CCTBX_XMLRPC_PORT", None)
      if cctbx_port is not None :
        uri = "http://localhost:%s/RPC2" % cctbx_port
        self.cctbx_server = xmlrpclib.ServerProxy(uri=uri)
        if self.verbose :
          print("Connecting to XML-RPC server on port %s" % cctbx_port)

    def start_server(self):
      if self.xmlrpc_server is not None :
        print("XML-RPC server started on port %d" % self.port)
        self.xmlrpc_server.serve_forever()

    def start_server_in_separate_thread(self):
      import threading
      t = threading.Thread(target=self.start_server)
      t.setDaemon(1)
      t.start()

    def set_socket_timeout(self, timeout):
      if self.xmlrpc_server is not None :
        self.xmlrpc_server.socket.settimeout(timeout)

    def timeout_func(self, *args):
      if self.xmlrpc_server is not None :
        self.xmlrpc_server.handle_request()
      return True

    def is_alive(self):
      return True

    # XXX: this should be replaced by the proper quit function for the program
    # being extended - e.g. cmd.quit() in PyMOL.
    def quit(self):
      print("quitting")
      sys.stdout.flush()
      os.kill(os.getpid(), signal.SIGKILL)

    def __getattr__(self, name):
      for module_object in self.supported_modules :
        if hasattr(module_object, name):
          return getattr(module_object, name)
      return None

except KeyboardInterrupt :
  raise
except ImportError :
  def external_xmlrpc_server(*args, **kwds):
    raise Exception("SimpleXMLRPCServer not available on this platform.")

  def external_cctbx_interface(*args, **kwds):
    raise Exception("SimpleXMLRPCServer not available on this platform.")

def test_server():
  class test_module(object):
    def echo_test(self):
      print("hello, world!")
      sys.stdout.flush()
      return True

#  os.environ["CCTBX_TEST_PORT"] = "48000"
  test_server = external_xmlrpc_interface("TEST", auto_start=False,
                                          verbose=False)
  module_object = test_module()
  test_server.add_module(module_object)
  test_server.start_server()

def coot_server():
  server = external_xmlrpc_interface("COOT",
    auto_start=False,
    verbose=True)
  server.set_socket_timeout(0.01)
  import coot
  import gobject
  server.add_module(coot)
  gobject.timeout_add(200, server.timeout_func)

if __name__ == "__main__" :
  #test_server()
  coot_server()

#---end
