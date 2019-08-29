from __future__ import absolute_import, division, print_function
from mod_python import apache
from six.moves import StringIO

def handler(req):
  req.content_type = "text/plain"
  from mod_python.util import FieldStorage
  FS = FieldStorage(req)
  logfile = StringIO()

  if req.filename.find("distl.signal_strength_bcsb")>=0:
    from spotfinder.servers.apache_bcsb import run as run_command
  elif req.filename.find("distl.signal_strength")>0:
    run_command = run
  elif req.filename.find("dials.find_spots")>0:
    from spotfinder.servers.apache_dials import run as run_command

  logfile.write(run_command(args=FS))
  log = logfile.getvalue()
  req.set_content_length(len(log))
  req.write(log)
  return apache.OK

def run(args, verbose=False):
  from libtbx.utils import Sorry
  import os
  from spotfinder.command_line.signal_strength import master_params
  from spotfinder.servers import LoggingFramework

  #For the Apache server version, do not allow site, user, or dataset preferences
  #all parameters are to be passed in through the http: query line

  logfile = LoggingFramework()

  phil_objects = []
  argument_interpreter = master_params.command_line_argument_interpreter(
    home_scope="distl")

  for key in args.keys():
      arg = "%s=%s"%(key,args.get(key,""))
      try: command_line_params = argument_interpreter.process(arg=arg)
      except Exception: return str(Sorry("Unknown file or keyword: %s" % arg))
      else: phil_objects.append(command_line_params)

  working_params = master_params.fetch(sources=phil_objects)
  params = working_params.extract()
  #working_params.show()

  if not os.path.isfile(params.distl.image):
    return str(Sorry("%s is not a readable file" % params.distl.image))

  print("Image: %s"%params.distl.image)

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
