from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
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
