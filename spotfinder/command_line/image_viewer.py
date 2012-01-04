# LIBTBX_SET_DISPATCHER_NAME distl.image_viewer
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1
import sys,os
from spotfinder.command_line.signal_strength import master_params

def run(args, command_name="distl.image_viewer"):
  help_str="""Same as distl.signal_strength (type that command for help) except that
  the image_viewer starts up an interactive GUI to visualize spotfinder spots.
"""

  if (len(args) == 0 or args[0] in ["H","h","-H","-h","help","--help","-help"]):
    print "usage:   %s image_filename [parameter=value ...]" % command_name
    print "example: %s lysozyme_001.img distl.res.outer=2.0 distl.res.inner=6.0 distl.minimum_spot_area=8"%command_name
    master_params.show(attributes_level=1,expert_level=1)
    print help_str
    return

  print "%s: characterization of candidate Bragg spots"%command_name

  phil_objects = []
  argument_interpreter = master_params.command_line_argument_interpreter(
    home_scope="distl")
  image_file_name = None
  moving_pdb_file_name = None
  for arg in args:
    if (os.path.isfile(arg)):
      if (image_file_name is None): image_file_name = arg
      else: raise Sorry("Too many file names.")
    else:
      try: command_line_params = argument_interpreter.process(arg=arg)
      except KeyboardInterrupt: raise
      except Exception: raise Sorry("Unknown file or keyword: %s" % arg)
      else: phil_objects.append(command_line_params)

  working_params = master_params.fetch(sources=phil_objects)
  params = working_params.extract()

  def raise_missing(what):
      raise Sorry("""\
Missing file name for %(what)s structure:
  Please add
    %(what)s=file_name
  to the command line to specify the %(what)s structure.""" % vars())

  if (image_file_name is None):
    if (params.distl.image is None): raise_missing("file name")
  else:
    params.distl.image = image_file_name

  working_params = master_params.format(python_object=params)
  #working_params.show(expert_level=1)

  #Now actually run the program logic
  from spotfinder.applications import image_viewer as app_image_viewer
  modeler = app_image_viewer.run_signal_strength_class(params=working_params.extract())
  modeler.view()

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
