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
from spotfinder.servers.spotfinder_server_read_file import generate_common_parameters # import dependency

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
