from BaseHTTPServer import BaseHTTPRequestHandler,HTTPServer
from iotbx.detectors import ImageFactory
from scitbx.array_family import flex
from libtbx.development.timers import Timer,Profiler
import StringIO, cgi, sys
from spotfinder.applications.stats_distl import optionally_add_saturation_webice,key_adaptor
from multiprocessing import Process, current_process, freeze_support

from urlparse import urlparse
#backward compatibility with Python 2.5
try: from urlparse import parse_qs
except: from cgi import parse_qs

def note(format, *args):
    sys.stderr.write('[%s]\t%s\n' % (current_process().name, format%args))

from spotfinder.servers.spotfinder_server_read_file import module_safe_items, module_image_stats
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
  try:
    port = int(sys.argv[1])
    NUMBER_OF_PROCESSES = int(sys.argv[2])
    if len(sys.argv)>3:
      outer_resolution = float(sys.argv[3])
    if len(sys.argv)>4:
      minimum_spot_area = int(sys.argv[4])
  except:
    print """
Usage:  libtbx.python mp_spotfinder_server_read_file.py <port number> <number of processes> [<outer resolution> [<minimum spot area.]]
"""
  common_parameters(outer_resolution,minimum_spot_area)

  server_address = ('', port)

  image_request_handler.protocol_version = "HTTP/1.0"

  print "Serving %d-process HTTP on"%NUMBER_OF_PROCESSES, server_address[0], "port", server_address[1], "..."
  print "To exit press Ctrl-C"
  runpool(server_address,NUMBER_OF_PROCESSES,image_request_handler)
