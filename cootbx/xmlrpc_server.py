
"""
Prototype for Coot XML-RPC plugin used by the Phenix GUI.  This is a fully
functional implementation, but the code used in Phenix contains many additional
methods in the coot_interface class which may also be called over XML-RPC.  A
simple example is update_model() in the coot_interface class.

libtbx.xmlrpc_utils contains tools used to launch and communicate with Coot
from another Python process, including a fault-sensitive replacement for the
ServerProxy class in xmlrpclib.  However, the client code could just as easily
be written in an entirely different language.
"""

from __future__ import division, print_function
import SimpleXMLRPCServer
import traceback
import os

from libtbx.utils import to_str

try :
  import coot_python
except Exception as e :
  print("Could not import coot_python module!")
  print("Coot GUI extensions will be disabled.")
  class empty(object):
    def main_menubar(self):
      return None
    def main_toolbar(self):
      return None
  coot_python = empty()

def coot_startup():
  gui = coot_interface()

class coot_interface(object):
  """
  Class for managing communications, including a toolbar button to toggle the
  XML-RPC connection (which doesn't actually change the socket connection,
  but determines whether requests are ignored or processed).  Any method
  defined for this class (or a subclass, if desired) will become callable
  over XMl-RPC.
  """
  def __init__(self):
    self.flag_enable_xmlrpc = True
    self.xmlrpc_server = None
    self._server_toggle_btn = None
    self._current_model = None
    toolbar = coot_python.main_toolbar()
    port = 40000
    if ("COOT_XMLRPC_PORT" in os.environ):
      port = string.atoi(os.environ["COOT_XMLRPC_PORT"])
    try :
      self.xmlrpc_server = coot_xmlrpc_server(
        interface=self,
        addr=("127.0.0.1", port))
      self.xmlrpc_server.socket.settimeout(0.01)
    except Exception as e :
      print("Error starting XML-RPC server:")
      print(str(e))
    else :
      import gobject
      import gtk
      print("xml-rpc server running on port %d" % port)
      # timeout used to be set to whatever the Phenix preferences have it as,
      # but 250ms seems to be a universally good choice, and much shorter
      # intervals screw up the Coot GUI (at least on Windows)
      gobject.timeout_add(250, self.timeout_func)
      if (toolbar is not None):
        self._server_toggle_btn = gtk.ToggleToolButton()
        self._server_toggle_btn.set_label("Connected to PHENIX")
        self._server_toggle_btn.set_is_important(True)
        toolbar.insert(self._server_toggle_btn, -1)
        self._server_toggle_btn.connect("clicked", self.toggle_update)
        self._server_toggle_btn.set_active(True)
        self._server_toggle_btn.show()

  def timeout_func(self, *args):
    if (self.xmlrpc_server is not None):
      self.xmlrpc_server.handle_request()
    return True

  def toggle_update(self, *args):
    if (self._server_toggle_btn is not None):
      if self._server_toggle_btn.get_active():
        self._server_toggle_btn.set_label("XML-RPC enabled")
        self.flag_enable_xmlrpc = True
      else :
        self._server_toggle_btn.set_label("XML-RPC disabled")
        self.flag_enable_xmlrpc = False

  def update_model(self, pdb_file):
    """
    Example of a user-defined function which will automatically be made
    accessible via XML-RPC.
    """
    import coot
    if (self.current_imol is None):
      self.current_imol = read_pdb(to_str(pdb_file))
    else :
      clear_and_update_molecule_from_file(self.current_imol,
        pdb_file)

class coot_xmlrpc_server(SimpleXMLRPCServer.SimpleXMLRPCServer):
  """
  Replacement for the SimpleXMLRPCServer class, which can call methods
  defined in multiple places.
  """
  def __init__(self, interface, **kwds):
    self._interface = interface
    SimpleXMLRPCServer.SimpleXMLRPCServer.__init__(self, **kwds)

  def _dispatch(self, method, params):
    """
    Handle a method request for either a) a method of the interface object, or
    b) wrapped C++ function exposed by the 'coot' module; both options will
    be tried (in that order) if necessary.  Note that Python's None instance
    is problematic, so it is replaced with -1.
    """
    if not self._interface.flag_enable_xmlrpc :
      return -1
    import coot
    result = -1
    func = None
    if hasattr(self._interface, method):
      func = getattr(self._interface, method)
    elif hasattr(coot, method):
      func = getattr(coot, method)
    if not hasattr(func, "__call__"):
      print("%s is not a callable object!" % method)
    else :
      try :
        result = func(*params)
      except Exception as e :
        traceback_str = "\n".join(traceback.format_tb(sys.exc_info()[2]))
        raise Exception("%s\nOriginal traceback:%s" % (str(e), traceback_str))
      else :
        if result is None :
          result = -1
    return result

if (__name__ == "__main__"):
  coot_startup()
