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

    except Sorry,e:
      pass

  for item in consume:
    args.remove(item)

  if len(args) > 0:
    raise Sorry("Not all arguments processed")

  params = horizons_phil.extract()
  assert params.distl.detector_format_version is not None
  message = "deprecated: %s and detector_format_version should not be specified in your phil file. Remove them from phil files of the form cxi-7.1.phil or xpp-7.1.phil, included from your xtal_target phil file"
  if params.distl.tile_translations is not None:
    raise Sorry(message%"tile_translations")
  if params.distl.quad_translations is not None:
    raise Sorry(message%"quad_translations")

  from spotfinder.applications.xfel.cxi_phil import cxi_versioned_extract
  args = ["distl.detector_format_version=%s"%params.distl.detector_format_version]

  versioned_extract = cxi_versioned_extract(args).persist.commands

  params.distl.quad_translations = versioned_extract.distl.quad_translations
  params.distl.tile_translations = versioned_extract.distl.tile_translations

  return params
