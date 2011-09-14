import pickle
import socket
import sys
import os

class FileClient:

  def __init__(self,
               host='',
               port=8000):
    if host=='127.0.0.1' or host=='localhost':
      host = ''
    if not host: # or host[0].isalpha():
      host = socket.getfqdn()
    self.host=host
    if self.host.find('\n')>-1:
      self.host = self.host.split('\n')[0]
    self.IP = socket.gethostbyname(self.host)
    self.port=int(port)
    self.sock=None

  def GetServerDetails(self,*arg,**kw):
    return self._remote_call('GetServerDetails',arg,kw)

  def GetHost(self,*arg,**kw):
    return self._remote_call('GetHost',arg,kw)

  def GetPort(self,*arg,**kw):
    return self._remote_call('GetPort',arg,kw)

  #def GetUser(self,*arg,**kw):
  #  return self._remote_call('GetUser',arg,kw)

##############################################################################
# File functions
##############################################################################

  def LockFile(self,*arg,**kw):
    return self._remote_call('LockFile',arg,kw)

  def ReadFile(self,*arg,**kw):
    return self._remote_call('ReadFile',arg,kw)

  def WriteFile(self,*arg,**kw):
    return self._remote_call('WriteFile',arg,kw)

  def ReadPickleFile(self,*arg,**kw):
    return self._remote_call('ReadPickleFile',arg,kw)

  def WritePickleFile(self,*arg,**kw):
    return self._remote_call('WritePickleFile',arg,kw)

  def UnlockFile(self,*arg,**kw):
    return self._remote_call('UnlockFile',arg,kw)

  def FileExists(self,*arg,**kw):
    return self._remote_call('FileExists',arg,kw)

##############################################################################
# Process commands
##############################################################################
  def getpid(self,*arg,**kw):
    return self._remote_call('getpid',arg,kw)

  def tester(self,*arg,**kw):
    return self._remote_call('tester',arg,kw)

##############################################################################
# Client functions
##############################################################################

  def _remote_call(self,
                   meth,
                   args=(),
                   kwds={}):
    result = None
    if self.sock is None:
      self.sock=socket.socket(socket.AF_INET,socket.SOCK_STREAM)
      self.sock.connect((self.host, self.port))
      self.send = self.sock.makefile('wb')
      self.recv = self.sock.makefile('rb')
    pickle.dump(meth,self.send,1) # binary by default
    pickle.dump(args,self.send,1)
    pickle.dump(kwds,self.send,1)

    # this flush somethings causes a crash WHY? maybe try/except
    self.send.flush()
    try:
      result = pickle.load(self.recv)
    except Exception:
      result = None
    #if result.__class__.__name__ == "PropagateExceptionClass":
    #  cmd = "raise %s, \"%s\"" % (result.type, result.message)
    #  exec cmd

    # Closing the socket connection increases the speed of transfer but ruins
    # the possibility of contacting the clients.
    self.send.close()
    self.recv.close()
    self.sock.close()
    self.sock=None
    return result

  def __repr__(self):
    outl = '\nClass '+self.__class__.__name__
    outl += '\n\thost '+str(self.host)
    outl += '\n\tIP   '+str(self.IP)
    outl += '\n\tport '+str(self.port)
    outl += '\n\tsock '+str(self.sock)
    try:
      outl += '\n\tsend '+str(self.send)
      outl += '\n\trecv '+str(self.recv)
    except Exception:
      pass
    return outl

  def shutdown(self):
    try:
      # multiple connections are sometimes required...
      self._remote_call('shutdown',(),{})
      self.sock=None
      self._remote_call('shutdown',(),{})
      self.sock=None
      self._remote_call('shutdown',(),{})
      self.sock=None
      self._remote_call('shutdown',(),{})
    except Exception:
      pass
    return None

#-----------------------------------------------------------------------------

def LockReadProcessWriteUnlock(client, file, id, func):
  """
  This method will lock a file, read the contents and pass it to
  the user defined function.  This function should return the
  return item and the new contents of the file.  The contents is
  writing and the file unlocked.
  """
  client.LockFile(file, id)
  try:
    lines = client.ReadFile(file)

    return_item = None
    if lines:
      return_item, lines = func(lines)
      client.WriteFile(lines, file)

  finally:
    client.UnlockFile(file, id)

  return return_item

def LockReadPickleProcessWritePickleUnlock(client, file, id, func):
  """
  Same as above except a pickled obj is used.
  """
  client.LockFile(file, id)
  try:
    obj = client.ReadPickleFile(file)

    return_item = None
    if obj:
      return_item, obj = func(obj)
      client.WritePickleFile(obj, file)

  finally:
    client.UnlockFile(file, id)

  return return_item

if __name__=="__main__":

  import FileServer
  import time
  import libtbx.load_env

  cmd = libtbx.env.under_dist("libtbx", "libtbx/server/StartServer.py")
  python_path = sys.executable

  if sys.platform == 'win32':
    os.spawnv(os.P_NOWAIT, python_path,
              (python_path, cmd)
              )
  else:
    os.spawnvp(os.P_NOWAIT, python_path,
              (python_path, cmd)
               )

  time.sleep(5)

  client = FileServer.GetServerClient()

  print client
  if client:
    print client.tester()
    #client.shutdown()
