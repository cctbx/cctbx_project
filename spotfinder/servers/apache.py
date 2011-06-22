from mod_python import apache
import StringIO,sys

def handler(req):
  req.content_type = "text/plain"
  from mod_python.util import FieldStorage
  FS = FieldStorage(req)
  logfile = StringIO.StringIO()
  logfile.write(run(args=FS))
  log = logfile.getvalue()
  req.set_content_length(len(log))
  req.write(log)
  return apache.OK

class LoggingFramework:
  def __init__(self):
    self.k = StringIO.StringIO()
    self.current_out = sys.stdout
    self.current_err = sys.stderr
    sys.stdout = self.k
    sys.stderr = self.k

  def __del__(self):
    sys.stdout = self.current_out
    sys.stderr = self.current_err
    self.k.flush()
    self.k.close()

  def getvalue(self): return self.k.getvalue()

def run(args, verbose=False):
  from libtbx.utils import Sorry
  import os
  from labelit.preferences import labelit_phil
  from spotfinder.command_line.signal_strength import master_params

  #For the Apache server version, do not allow site, user, or dataset preferences
  #all parameters are to be passed in through the http: query line
  labelit_phil.rollback_dataset_preferences()

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

  print "Image: %s"%params.distl.image

  from spotfinder.applications import signal_strength
  try:
    signal_strength.run_signal_strength(params)
  except Exception:
    import traceback
    logger = StringIO.StringIO()
    logger.write(
    "Sorry, can't process %s.  Please contact authors.\n"% params.distl.image)
    traceback.print_exc(file=logger)
    return str(Sorry( logger.getvalue() ))

  return logfile.getvalue()
