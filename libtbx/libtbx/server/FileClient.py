import pickle
import cPickle
import socket
import SocketServer
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
    return self._remote_call('ReadFile',arg,kw)

  def WritePickleFile(self,*arg,**kw):
    return self._remote_call('WriteFile',arg,kw)

  def UnlockFile(self,*arg,**kw):
    return self._remote_call('UnlockFile',arg,kw)

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
    except:
      result = None
    if result.__class__.__name__ == "PropagateExceptionClass":
      cmd = "raise %s, \"%s\"" % (result.type, result.message)
      exec cmd

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
    except:
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
    except:
      pass
    return None

if __name__=="__main__":

  import FileServer
  import time

  cmd = 'StartServer.py'
  python_path = sys.executable

  if sys.platform == 'win32':
    os.spawnv(os.P_NOWAIT, python_path,
              (python_path, cmd)
              )
  else:
    os.spawnvp(os.P_NOWAIT, python_path,
              (python_path, cmd)
               )

  time.sleep(2)

  client = FileServer.GetServerClient()

  print client
  if client:
    print client.tester()
    #client.shutdown()
