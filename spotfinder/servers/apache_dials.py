from __future__ import division
from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
from builtins import str
import io

def run(args, verbose=False):
  from libtbx.utils import Sorry
  try:
    from dials.array_family import flex
  except ImportError:
    return str(Sorry("DIALS is not configured"))

  from iotbx.phil import parse
  import os
  from spotfinder.servers import LoggingFramework
  from dials.array_family import flex
  from dxtbx.datablock import DataBlockFactory
  phil_scope = parse("""
  file_name = None
    .type = str
  frame_number = None
    .type = int
  stats = True
    .type = bool
  include scope dials.algorithms.spot_finding.factory.phil_scope
  """, process_includes=True)

  #For the Apache server version, do not allow site, user, or dataset preferences
  #all parameters are to be passed in through the http: query line

  logfile = LoggingFramework()

  phil_objects = []

  for key in list(args.keys()):
    arg = "%s=%s"%(key,args.get(key,""))
    try: phil_objects.append(parse(arg))
    except Exception: return str(Sorry("Unknown file or keyword: %s" % arg))

  working_params = phil_scope.fetch(sources=phil_objects)
  params = working_params.extract()
  #working_params.show()

  if not os.path.isfile(params.file_name):
    return str(Sorry("%s is not a readable file" % params.file_name))

  print("Image: %s\n"%params.file_name)

  try:
    datablock = DataBlockFactory.from_filenames([params.file_name])[0]
    imageset = datablock.extract_imagesets()[0]
    if datablock.num_images() > 0 and params.frame_number is not None:
      print("Frame number", params.frame_number)
      imageset = imageset[params.frame_number:params.frame_number+1]
      datablock = DataBlockFactory.from_imageset(imageset)[0]
    reflections = flex.reflection_table.from_observations(datablock, params)

    if params.stats:
      from dials.algorithms.spot_finding.per_image_analysis import  stats_single_image
      print(stats_single_image(imageset, reflections, i=None, resolution_analysis=True, plot=False))

  except Exception:
    import traceback
    logger = io.StringIO()
    logger.write(
    "Sorry, can't process %s.  Please contact authors.\n"% params.file_name)
    traceback.print_exc(file=logger)
    return str(Sorry( logger.getvalue() )) + logfile.getvalue()

  print("Found %d strong reflections"%len(reflections))

  return logfile.getvalue()
