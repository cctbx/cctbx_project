import pickle
import socket
from FileServer import StartServer, BindClient
from FileClient import FileClient

host = None
client = None
try:
  f = open('server.port', 'rb')
  host, port = pickle.load(f)
  f.close()
except:
  print 'failed to read server port file'

if host:
  client = BindClient(host, port)

if not client:
  print 'Failed to find running server'
  thread, port = StartServer()
  try:
    f = open('server.port', 'wb')
    pickle.dump((socket.getfqdn(),
                 port),
                f
                )
    f.close()
  except:
    print 'failed to write server port file'

  thread.run()
