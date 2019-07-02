from __future__ import absolute_import, division, print_function

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

  horizons_phil = master_phil.fetch(
    sources=[phil.parse(file_name=path, process_includes=True)])

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

    except Sorry as e:
      pass

  for item in consume:
    args.remove(item)

  if len(args) > 0:
    raise Sorry("Not all arguments processed")

  params = horizons_phil.extract()
  if params.distl.tile_translations is not None and params.distl.quad_translations is not None:
    return params

  from spotfinder.applications.xfel.cxi_phil import cxi_versioned_extract
  args = ["distl.detector_format_version=%s"%params.distl.detector_format_version]

  versioned_extract = cxi_versioned_extract(args).persist.commands

  if params.distl.quad_translations is None:
    params.distl.quad_translations = versioned_extract.distl.quad_translations

  if params.distl.tile_translations is None:
    params.distl.tile_translations = versioned_extract.distl.tile_translations

  return params
