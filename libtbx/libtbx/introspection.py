import sys

def frame_object(frames_back=0):
  try: raise Exception
  except:
    t = sys.exc_info()[2]
    f = t.tb_frame
    for i in xrange(frames_back+1):
      if (f.f_back is None): break
      f = f.f_back
    return f

def varnames(frames_back=0):
  return frame_object(frames_back=frames_back+1).f_code.co_varnames

def adopt_init_args(exclusions=[], prefix="", frames_back=0):
  frame = frame_object(frames_back=frames_back+1)
  varnames = frame.f_code.co_varnames
  exclusions.append(varnames[0]) # self
  init_locals = frame.f_locals
  self = init_locals[varnames[0]]
  for varname in varnames:
    if (varname not in exclusions):
      setattr(self, prefix+varname, init_locals[varname])
  if ("__init__varnames__" not in exclusions):
    setattr(self, "__init__varnames__", varnames)

class caller_location:

  def __init__(self, frames_back=0):
    f = frame_object(frames_back=frames_back+1)
    self.file_name = f.f_code.co_filename
    self.line_number = f.f_lineno

  def __str__(self):
    return "%s(%d)" % (self.file_name, self.line_number)

def check_point(frames_back=0):
  print caller_location(frames_back=frames_back+1)
  sys.stdout.flush()
