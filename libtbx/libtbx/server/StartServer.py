import pickle
import socket
from FileServer import StartServer, GetServerClient, WriteServerFile

client = GetServerClient()

if not client:
  print 'Failed to find running server'
  thread, port = StartServer()
  WriteServerFile(port)
  thread.start()
