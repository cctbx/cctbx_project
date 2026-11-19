"""With nothing listening on a localhost port inside the range the OS uses for
automatically chosen local ports (32768-60999 on Linux, which holds the
viewer ports of libtbx.xmlrpc_utils.external_program_server), a retried
connection can be given the server's own port and connect to itself (TCP
simultaneous open). The request is then read back as the reply and the GUI
reported "BadStatusLine: POST /RPC2 HTTP/1.1" while Coot was starting or
after it had closed. Such a connection must count as refused, which keeps
the request queued.
"""
import errno
import socket

from libtbx import xmlrpc_utils

class FakeSocket(object):
  """A connected socket; self_connected mimics a connection to itself."""
  def __init__(self, self_connected):
    self.self_connected = self_connected
    self.closed = False
    self.sent = []
  def getsockname(self):
    if self.self_connected:
      return ("127.0.0.1", 40500)
    return ("127.0.0.1", 51234)
  def getpeername(self):
    return ("127.0.0.1", 40500)
  def setsockopt(self, *args):
    pass
  def sendall(self, data):
    self.sent.append(data)
  def close(self):
    self.closed = True

def exercise():
  # http.client takes socket.create_connection when a connection is made
  sockets = []
  def fake_create_connection(self_connected):
    def create_connection(*args, **kwds):
      sockets.append(FakeSocket(self_connected))
      return sockets[-1]
    return create_connection
  connection_class = xmlrpc_utils.NoSelfConnectHTTPConnection
  saved = socket.create_connection
  try:
    socket.create_connection = fake_create_connection(True)
    connection = connection_class("127.0.0.1", 40500)
    try:
      connection.connect()
    except (IOError, OSError) as e:
      assert e.errno == errno.ECONNREFUSED, e.errno
    else:
      raise AssertionError("a connection to itself was accepted")
    assert sockets[-1].closed and (sockets[-1].sent == [])
    # through the GUI's caching proxy: no error, and the request stays
    # queued as for a refused connection
    proxy = xmlrpc_utils.ServerProxy(uri="http://127.0.0.1:40500/RPC2")
    assert proxy.is_alive() is None
    assert [ name for (name, params) in proxy._pending ] == ["is_alive"], \
      proxy._pending
    assert proxy.get_error_messages() != []
    # an ordinary connection is kept
    socket.create_connection = fake_create_connection(False)
    connection = connection_class("127.0.0.1", 40500)
    connection.connect()
    assert (connection.sock is sockets[-1]) and (not sockets[-1].closed)
  finally:
    socket.create_connection = saved

if __name__ == "__main__":
  exercise()
  print("OK")
