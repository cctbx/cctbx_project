from BaseHTTPServer import HTTPServer
import cgi, sys
from multiprocessing import Process, current_process

from urlparse import urlparse
#backward compatibility with Python 2.5
try: from urlparse import parse_qs
except Exception: from cgi import parse_qs

def note(format, *args):
    sys.stderr.write('[%s]\t%s\n' % (current_process().name, format%args))

from spotfinder.servers.spotfinder_server_read_file import image_request_handler as irhbase
from spotfinder.servers.spotfinder_server_read_file import common_parameters

class image_request_handler(irhbase):

  def log_message(self, format, *args):
    note(format, *args)

def serve_forever(server):
    note('starting server')
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass


def runpool(address, number_of_processes,handler):
    # create a single server object -- children will each inherit a copy
    server = HTTPServer(address, handler)

    # create child processes to act as workers
    for i in range(number_of_processes-1):
        Process(target=serve_forever, args=(server,)).start()

    # main process also acts as a worker
    serve_forever(server)


if __name__=="__main__":
  import sys
  outer_resolution = None
  minimum_spot_area = None
  minimum_signal_height = None
  try:
    port = int(sys.argv[1])
    NUMBER_OF_PROCESSES = int(sys.argv[2])
    if len(sys.argv)>3:
      outer_resolution = float(sys.argv[3])
    if len(sys.argv)>4:
      minimum_spot_area = int(sys.argv[4])
    if len(sys.argv)>5:
      minimum_signal_height = float(sys.argv[5])
  except Exception:
    print """
Usage:  libtbx.python mp_spotfinder_server_read_file.py <port number> <number of processes> [<outer resolution> [<minimum spot area> [<minimum signal height>]]]
"""
  common_parameters(outer_resolution,minimum_spot_area,minimum_signal_height)

  server_address = ('', port)

  image_request_handler.protocol_version = "HTTP/1.0"

  print "Serving %d-process HTTP on"%NUMBER_OF_PROCESSES, server_address[0], "port", server_address[1], "..."
  print "To exit press Ctrl-C"
  runpool(server_address,NUMBER_OF_PROCESSES,image_request_handler)
