from __future__ import absolute_import, division, print_function
import socket
import random
import os, sys
from string import *

def GetSocket(HOST=""):
  unbind=1
  lower_limit=20000
  upper_limit=65535
  if not HOST:
    HOST = socket.getfqdn(socket.gethostname())
  i = 0
  while unbind:
    portnumber = int(random.random()*(upper_limit-lower_limit)+lower_limit)
    try:
      s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
      s.bind((HOST, portnumber))
      s.listen(5)
      portnumber = (s.getsockname()[1])
      unbind=0
      return s
    except Exception:
      s.close()
      del s

    if i > 5:
      HOST = 'localhost'
    if i > 10:
      print('Port assignment failed.', HOST, portnumber)
      sys.exit()

    i += 1


class SocketConnection:
  def __init__(self, host, dir='', port=0, pid=0, user=''):
    """
    Contains all the information about a socket connection and more.

    Arguments.

      host - hostname of machine #IP address of the host

      dir - Working directory on host

      port - Port on host

      pid - PID of process on host

      name - hostname determined from host

      user - process owner

    Attributes

      secret - secret key for passing sensitive data
    """
    self.SetHost(host)
    self.IP   = socket.gethostbyname(self.host)
    self.port = int(port)

    try:
      dir = self._deconPath(dir,[])
    except Exception:
      dir = ' '
    self.dir=dir

    self.pid=int(pid)
    self.user = user

  def __repr__(self):
    try:
      return '\n '+'-'*10+' '+self.host+' ('+self.IP+') '+ \
             '\n '+'-'*10+' '+os.path.join(*tuple(self.dir))+ \
             '\n '+'-'*10+' Port '+str(self.port)+' PID '+str(self.pid)
    except Exception:
      return 'socket connection '+str(self.__dict__)

  def _deconPath(self, path, plist):
    path = os.path.split(path)
    tmp_list = plist[:]
    if path[1]:
      tmp_list.append(path[1])
      tmp_list = self._deconPath(path[0], tmp_list)
    else:
      tmp_list.append(path[0])
      tmp_list.reverse()
    return tmp_list

  def GetStr(self):
    return self.__repr__()

  def GetPort(self):
    return int(self.port)

  def GetHost(self):
    return self.host

  def GetHostname(self):
    return self.host

  def GetIP(self):
    return self.IP

  def GetDir(self):
    dir = os.path.join(*self.dir)
    dir = replace(dir,'\\','/')
    return dir

  def GetPID(self):
    return self.pid

  def GetUser(self):
    return self.user

  def SetPort(self, port):
    self.port=int(port)

  def SetHost(self, host):
    for i in host:
      if not i in letters:
        self.host = socket.gethostbyaddr(host)[0]
        break
    else:
      self.host=socket.getfqdn(host)

  def SetHostOLD(self, host):
    for i in host:
      if i in letters:
        self.host = socket.gethostbyname(host)
        break
    else:
      self.host=host

  def SetDir(self, dir):
    self.dir=self._deconPath(dir, [])

  def SetPID(self, pid):
    self.pid=pid

  def SetUser(self, user):
    self.user=user
