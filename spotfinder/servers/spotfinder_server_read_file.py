from __future__ import absolute_import, division, print_function
from BaseHTTPServer import BaseHTTPRequestHandler
from scitbx.array_family import flex
from libtbx.development.timers import Timer
from six.moves import StringIO
import cgi, sys, copy
from spotfinder.applications.stats_distl import optionally_add_saturation_webice,key_adaptor

from urlparse import urlparse
#backward compatibility with Python 2.5
try: from urlparse import parse_qs
except Exception: from cgi import parse_qs

def module_safe_items(image):
  return [
      ("%7d","Spot Total",key_adaptor(image,'N_spots_total')),
      ("%7d","Method-2 Resolution Total",key_adaptor(image,'N_spots_resolution')),
      ("%7d","Good Bragg Candidates",key_adaptor(image,'N_spots_inlier')),
      ("%7.0f","Total integrated signal, pixel-ADC units above local background",
                 key_adaptor(image,'ad_hoc_signal1')),
      ("%7d","Ice Rings",key_adaptor(image,'ice-ring_impact')),

      ("%7.2f","Method 1 Resolution",key_adaptor(image,'distl_resolution')),
      ("%7.2f","Method 2 Resolution",key_adaptor(image,'resolution')),

      ("%7.1f","Maximum unit cell",key_adaptor(image,'maxcel')),
  ]

def module_image_stats(S,key):
    # List of spots between specified high- and low-resolution limits
    image = S.images[key]
    spots = image.__getitem__('spots_inlier')

    integrated = 0.0

    for i,spot in enumerate(spots):
     integrated += flex.sum(spot.wts)
    image["ad_hoc_signal1"]=integrated

    canonical_info = []
    canonical_info.extend(module_safe_items(image))
    optionally_add_saturation_webice(canonical_info,image)

    for item in canonical_info:
      if item[2]==None:
        print("%63s : None"%item[1])
      else:
        print("%63s : %s"%(item[1],item[0]%item[2]))

class image_request_handler(BaseHTTPRequestHandler):

  def shutdown(self):
      def my_shutdown(arg1):
        print("IN SHUTDOWN THREAD")
        arg1.server.shutdown()
      import thread
      thread.start_new_thread(my_shutdown,(self,))
      #must be called in a different thread or deadlock.
      log = ""
      self.send_response(200)
      self.send_header("Content-type", 'text/plain')
      self.send_header("Content-length",len(log))
      self.end_headers()
      self.wfile.write(log)

  def do_POST(self):
    T = Timer("do_POST")
    parsed = urlparse(self.path)
    qs = parse_qs(parsed.query)

    expect = self.headers.getheaders("Expect")
    if len(expect)>=1:
      if True in [item.find("200")>=0 for item in expect]:
        self.send_response(200) # untested; has no apparent affect on libcurl
        return

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
    post_data = StringIO(data)

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
    print("*****************************")
    for item in parts.keys():
      if len(parts[item][0])< 1000:
        print(item, parts[item])
    print("*****************************")

    if parts["filename"][0].find("EXIT")>=0:
      self.shutdown()
      return

    from spotfinder.diffraction.imagefiles import spotfinder_image_files as ImageFiles
    from spotfinder.diffraction.imagefiles import Spotspickle_argument_module
    response_params = copy.deepcopy(common_parameters_singleton).extract()

    Files = ImageFiles(Spotspickle_argument_module(parts["filename"][0]),response_params)

    print("Final image object:")
    Files.images[0].show_header()
    print("beam_center_convention",Files.images[0].beam_center_convention)
    print("beam_center_reference_frame",Files.images[0].beam_center_reference_frame)

    logfile = StringIO()
    if response_params.distl.bins.verbose: sys.stdout = logfile

    from spotfinder.applications.wrappers import spotfinder_factory
    S = spotfinder_factory(None, Files, response_params)
    print()
    sys.stdout = sys.__stdout__

    frames = Files.frames()

    sys.stdout = logfile

    print("Image: %s"%parts["filename"][0])
    from spotfinder.applications.stats_distl import pretty_image_stats,notes
    for frame in frames:
      #pretty_image_stats(S,frame)
      #notes(S,frames[0])
      module_image_stats(S,frame)

    sys.stdout = sys.__stdout__
    log = logfile.getvalue()
    print(log)

    ctype = 'text/plain'
    self.send_response(200)
    self.send_header("Content-type", ctype)
    self.send_header("Content-length",len(log))
    self.end_headers()
    self.wfile.write(log)
    self.opt_logging()

  def opt_logging(self):
    pass

  def do_GET(self):
    T = Timer("do_GET")
    parsed = urlparse(self.path)
    qs = parse_qs(parsed.query)

    expect = self.headers.getheaders("Expect")
    if len(expect)>=1:
      if True in [item.find("200")>=0 for item in expect]:
        self.send_response(200) # untested; has no apparent affect on libcurl
        return

    log = self.do_GET_run(qs)

    ctype = 'text/plain'
    self.send_response(200)
    self.send_header("Content-type", ctype)
    self.send_header("Content-length",len(log))
    self.end_headers()
    self.wfile.write(log)
    self.opt_logging()

  def do_GET_run(self,qs): #similar to the run() function in apache.py module
    from libtbx.utils import Sorry
    import os
    from spotfinder.servers import LoggingFramework

    base_params = copy.deepcopy(common_parameters_singleton)
    argument_interpreter = base_params.command_line_argument_interpreter()
    phil_objects = []

    for key in qs.keys():
      arg = "%s=%s"%(key,qs.get(key,"")[0])
      try: command_line_params = argument_interpreter.process(arg=arg)
      except Exception: return str(Sorry("Unknown file or keyword: %s" % arg))
      else: phil_objects.append(command_line_params)

    working_params = base_params.fetch(sources=phil_objects)
    params = working_params.extract()
    #working_params.show()
    if not os.path.isfile(params.distl.image):
      return  str(Sorry("%s is not a readable file" % params.distl.image))

    print("Image: %s"%params.distl.image)

    logfile = LoggingFramework()
    from spotfinder.applications import signal_strength
    try:
      signal_strength.run_signal_strength(params)
    except Exception:
      import traceback
      logger = StringIO()
      logger.write(
      "Sorry, can't process %s.  Please contact authors.\n"% params.distl.image)
      traceback.print_exc(file=logger)
      return str(Sorry( logger.getvalue() ))

    return logfile.getvalue()

common_parameters_singleton = None

def generate_common_parameters(input_parameters):
  global common_parameters_singleton
  common_parameters_singleton = input_parameters
