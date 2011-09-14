import threading
import pickle
import SocketServer
import sys
import traceback
import os
import socket
import time

import FileClient
import FileSocket

# Socket Server

class FileServer:
  def __init__(self, port='', handler=''):
    sys.setcheckinterval(0)

    local = socket.getfqdn()
    server_address = (local, port)

    if handler:
      ddbs = _FileServer(server_address, handler)
    else:
      ddbs = _FileServer(server_address, FileRequestHandler)

    # now serve requests forever
    ddbs.keep_alive = 1
    while ddbs.keep_alive:
      ddbs.handle_request()


class _FileServer(SocketServer.ThreadingTCPServer):

  def server_bind(self):
    """Override server_bind to store the server name."""
    SocketServer.ThreadingTCPServer.server_bind(self)
    host, port = self.socket.getsockname()
    if not host or host == '0.0.0.0':
      host = socket.getfqdn()

##############################################################################
# Attributes accessable using self.server.server_name for example
##############################################################################
    self.server_name = host
    self.server_port = port
    self.main_thread = threading.enumerate()[0]
    self.ServerLock = threading.Lock()
    self.pid = os.getpid()

    self.FileLockDictionary = {}

class FileRequestHandler(SocketServer.StreamRequestHandler):

#------------------------------------------------------------------------------
# File functions
#------------------------------------------------------------------------------

  def _LockFile(self, file, id):
    DICT = self.server.FileLockDictionary
    rc = 1
    self.server.ServerLock.acquire()
    try:
      if id in DICT.keys():
        rc = 0
      if file in DICT.values():
        rc = 0
      if rc:
        DICT[id] = file
    finally:
      self.server.ServerLock.release()
    return rc

  def LockFile(self, file, id):
    rc = self._LockFile(file, id)
    while not rc:
      time.sleep(.01)
      rc = self._LockFile(file, id)

  def ReadPickleFile(self, file):
    try:
      f = open(file, 'rb')
      obj = pickle.load(f)
      f.close()
    except Exception:
      print 'failed to load',file
      obj = None
    return obj

  def WritePickleFile(self, obj, file):
    try:
      f = open(file, 'wb')
      pickle.dump(obj, f)
      f.close()
    except Exception:
      print 'failed to dump',file

  def ReadFile(self, file):
    try:
      f = open(file, 'rb')
      lines = f.read()
      f.close()
    except Exception:
      print 'failed to read',file
      lines = ''
    return lines

  def WriteFile(self, lines, file):
    try:
      f = open(file, 'wb')
      f.write(lines)
      f.close()
    except Exception:
      print 'failed to write',file

  def _UnlockFile(self, file, id):
    DICT = self.server.FileLockDictionary
    rc = 1
    self.server.ServerLock.acquire()
    try:
      if not id in DICT.keys():
        rc = 0
      if not file in DICT.values():
        rc = 0
      if rc:
        del DICT[id]
    finally:
      self.server.ServerLock.release()
    return rc

  def UnlockFile(self, file, id):
    rc = self._UnlockFile(file, id)
    while not rc:
      time.sleep(.01)
      rc = self._UnlockFile(file, id)

  def FileExists(self, file):
    return os.path.exists(file)

#------------------------------------------------------------------------------
# Server information
#------------------------------------------------------------------------------
  def GetServerDetails(self):
    #user = self.GetUser()
    server_sc = self.GetSocketConnection()
    #return user, server_sc
    return server_sc

  #def GetUser(self):
  #  return GetUser()

  def GetSocketConnection(self):
    """
    Returns the socket connection data
    """
    s = FileSocket.SocketConnection(self.server.server_name,
                                      os.getcwd(),
                                      self.server.server_port,
                                      os.getpid()
                                      )
    return s

  def GetHost(self):
    return self.server.server_name

  def GetPort(self):
    return self.server.server_port

#------------------------------------------------------------------------------
# Client tester
#------------------------------------------------------------------------------
  def _tester(self, host, port):
    client = FileClient.FileClient(host = host,
                                     port = port)
    try:
      client.tester()
      self.event.set()
    except Exception:
      pass

  def ClientTest(self, host, port, timeout=2):
    self.event = threading.Event()
    self.event.clear()

    t = threading.Thread(target=self._tester,
                         args=(host,port))
    t.setDaemon(1)
    t.start()

    self.event.wait(timeout) # timed wait may be too long

    if self.event.isSet():
      return 1
    return None

#------------------------------------------------------------------------------
# Misc
#------------------------------------------------------------------------------
  def getpid(self):
    return os.getpid()

  def shutdown(self):
    self.server.ServerLock.acquire()
    try:
      self.server.keep_alive = 0
    finally:
      self.server.ServerLock.release()

  def tester(self):
    """
    Called by client to test if server is alive
    """
    return " Hello from\n\t%s\n on %s %s" \
           % (str(self), str(self.server.server_name), sys.platform)

#------------------------------------------------------------------------------
# handler
#------------------------------------------------------------------------------
  def handle(self):
    """
    This method handles the arguments of a remote call
    """
    try:

      while self.server.keep_alive: # still alive?
        # get method name from client
        try:
          method = pickle.load(self.rfile) # use pickle to avoid blocking
        except Exception:
          #print 'failed to load method'
          #self.wfile.flush() #doesn't fixed socket problem!!!!!!!!!!
          break
        #except EOFError, socket.error: # socket closed?
        #  break

        if method == 'shutdown':
          self.server.keep_alive = 0

        # get arguments from client
        try:
          args = pickle.load(self.rfile)
        except Exception:
          args = ()
          print 'Method',method
          print 'error in args ',args
        try:
          kw = pickle.load(self.rfile)
        except Exception:
          kw = {}
          print 'error in kw ',kw

        #print method
        #print args
        #print kw

        # get method pointer
        meth_obj = getattr(self,method)
        #print meth_obj

        # call method
        try:
          result = apply(meth_obj,args,kw)
          pickle.dump(result,self.wfile,1) # binary by default
        except Exception:
          try:
            pickle.dump('0\n',self.wfile,1)
          except Exception:
            pass
          print 'Method',meth_obj
          print 'Args  ',args
          print 'Kw    ',kw
          raise
        try:
          self.wfile.flush()
        except Exception:
          pass

    except Exception:
      print self.__class__.__name__+'.handle: '
      traceback.print_exc()
      pass

    sys.stdout.flush()

#-----------------------------------------------------------------------------

def StartServer():
  s = FileSocket.GetSocket()
  port = (s.getsockname()[1])
  t = threading.Thread(target=FileServer,
                       args=(port,))
  #t.setDaemon(1)
  return (t,port)

def BindClient(host, port):
  nobind = 1
  while nobind:
    try:
      client = FileClient.FileClient(host = host,
                                     port = port
                                     )
      client.tester()
      nobind = 0
    except Exception:
      time.sleep(.2)
      nobind += 1
      if nobind > 10:
        return None

  return client

def ReadServerFile():
  host = None
  port = None
  try:
    f = open('server.port', 'rb')
    host, port = pickle.load(f)
    f.close()
  except Exception:
    print 'failed to read server port file'
  return host, port

def WriteServerFile(port):
  try:
    f = open('server.port', 'wb')
    pickle.dump((socket.getfqdn(),
                 port),
                f
                )
    f.close()
  except Exception:
    print 'failed to write server port file'

def GetServerClient():
  host, port = ReadServerFile()
  client = None
  if host:
    client = BindClient(host, port)

  return client
