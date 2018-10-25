from __future__ import absolute_import, division, print_function
from .FileServer import StartServer, GetServerClient, WriteServerFile

def run():
  client = GetServerClient()
  if not client:
    print('Failed to find running server')
    thread, port = StartServer()
    WriteServerFile(port)
    thread.start()

if (__name__ == "__main__"):
  run()
