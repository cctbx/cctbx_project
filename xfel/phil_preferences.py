from __future__ import division

def load_cxi_phil(path, args=[]):
  from labelit.phil_preferences import iotbx_defs, libtbx_defs
  from iotbx import phil

  cxi_phil_str = iotbx_defs + libtbx_defs
  stream = open(path)
  master_phil = phil.parse(input_string=cxi_phil_str, process_includes=True)
  horizons_phil = master_phil.fetch(sources=[phil.parse(stream.read())])
  stream.close()

  from libtbx.phil.command_line import argument_interpreter
  from libtbx.utils import Sorry

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
