from __future__ import division

def load_cxi_phil(path, args=[]):
  import os
  from labelit.phil_preferences import iotbx_defs, libtbx_defs
  from iotbx import phil
  from libtbx.phil.command_line import argument_interpreter
  from libtbx.utils import Sorry

  exts = ["", ".params", ".phil"]
  foundIt = False
  for ext in exts:
    if os.path.exists(path + ext):
      foundIt = True
      path += ext
      break
  if not foundIt:
    raise Sorry("Target not found: " + path)

  master_phil = phil.parse(input_string=iotbx_defs + libtbx_defs,
                           process_includes=True)

  # Temporarily set the working directory to the directory of the
  # input file.  This emulates cpp(1)-like behaviour of "include"
  # statements, except relative paths of any nested includes will be
  # expanded with respect to this directory as well.
  cwd = os.getcwd()
  (head, tail) = os.path.split(path)
  if os.path.isdir(head):
    os.chdir(head)

  stream = open(tail)
  horizons_phil = master_phil.fetch(
    sources=[phil.parse(stream.read(), process_includes=True)])
  stream.close()
  os.chdir(cwd)

  argument_interpreter = argument_interpreter(
    master_phil=master_phil
  )
  consume = []
  for arg in args:
    try:
      command_line_params = argument_interpreter.process(
        arg=arg
      )
      horizons_phil = horizons_phil.fetch(sources=[command_line_params,])
      consume.append(arg)

    except Sorry,e:
      pass

  for item in consume:
    args.remove(item)

  return horizons_phil.extract()
