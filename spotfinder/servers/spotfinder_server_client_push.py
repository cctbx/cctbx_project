from SimpleHTTPServer import SimpleHTTPRequestHandler
from BaseHTTPServer import BaseHTTPRequestHandler,HTTPServer
from iotbx.detectors import ImageFactory
from scitbx.array_family import flex
from libtbx.development.timers import Timer,Profiler
import StringIO, cgi, sys

from urlparse import urlparse
#backward compatibility with Python 2.5
try: from urlparse import parse_qs
except: from cgi import parse_qs

class image_request_handler(BaseHTTPRequestHandler):

  def do_POST(self):
    T = Timer("do_POST")
    parsed = urlparse(self.path)
    qs = parse_qs(parsed.query)

    # Get arguments by reading body of request.
    # We read this in chunks to avoid straining
    # socket.read(); around the 10 or 15Mb mark, some platforms
    # begin to have problems (bug #792570).
    max_chunk_size = 10*1024*1024
    size_remaining = int(self.headers["content-length"])
    L = []
    while size_remaining:
        chunk_size = min(size_remaining, max_chunk_size)
        L.append(self.rfile.read(chunk_size))
        size_remaining -= len(L[-1])
    data = ''.join(L)
    post_data = StringIO.StringIO(data)

    # Parse the multipart/form-data
    contentTypeHeader = self.headers.getheaders('content-type').pop()

    # Extract the boundary parameter in the content-type header
    headerParameters = contentTypeHeader.split(";")
    boundary = headerParameters[1].split("=")
    boundary = boundary[1].strip()

    parts = cgi.parse_multipart(post_data,
      {"boundary":boundary,
       "content-disposition":self.headers.getheaders('content-disposition')
      })
    print "*****************************"
    for item in parts.keys():
      print item, len(parts[item][0])
    print "*****************************"

    from iotbx.detectors.image_from_http_request import module_or_slice_from_http_request
    imgobj = module_or_slice_from_http_request(parts)
    imgobj.read()
    print "Final image object:"
    imgobj.show_header()

    from spotfinder.diffraction.imagefiles import image_files, file_names
    from spotfinder.diffraction.imagefiles import Spotspickle_argument_module

    from spotfinder.applications.overall_procedure import spotfinder_no_pickle

    class server_imagefiles(image_files):
      def __init__(self): pass

    Files = server_imagefiles()
    Files.filenames = file_names(Spotspickle_argument_module(imgobj.filename))
    Files.images = [imgobj]

    S = spotfinder_no_pickle(Files, s3_passthru = "-s3 4",
                             spot_convention = 0)

    frames = Files.frames(1)

    logfile = StringIO.StringIO()
    sys.stdout = logfile

    from spotfinder.applications.stats_distl import pretty_image_stats,notes
    for frame in frames:
      pretty_image_stats(S,frame)
      notes(S,frames[0])

    sys.stdout = sys.__stdout__
    log = logfile.getvalue()
    print log

    ctype = 'text/plain'
    self.send_response(200)
    self.send_header("Content-type", ctype)
    self.send_header("Content-length",len(log))
    self.end_headers()
    self.wfile.write(log)
    self.opt_logging()

  def opt_logging(self):
    pass

if __name__=="__main__":
  import sys
  try:
    port = int(sys.argv[1])
  except:
    print """
Usage:  libtbx.python adsc_server.py <port number>
"""
  server_address = ('', port)

  image_request_handler.protocol_version = "HTTP/1.0"
  httpd = HTTPServer(server_address, image_request_handler)

  sa = httpd.socket.getsockname()
  print "Serving HTTP on", sa[0], "port", sa[1], "..."
  httpd.serve_forever()
