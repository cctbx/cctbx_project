# LIBTBX_SET_DISPATCHER_NAME distl.mp_spotfinder_server_read_file
import libtbx.phil
import libtbx.phil.command_line
from libtbx.utils import Sorry
from spotfinder.command_line.signal_strength import master_params

def run(args, command_name="distl.mp_spotfinder_server_read_file"):
  help_str="""Multiprocessing server to find Bragg spots & quantify signal strength.
Full documentation: http://cci.lbl.gov/publications/download/ccn_jul2010_page18.pdf
Allowed parameters:
  distl.port=8125
  distl.processors=1
  distl.minimum_signal_strength=2.5
  distl.minimum_spot_area=5
  distl.res.outer=3.0 [outer resolution limit]
"""

  if (len(args)>=1 and args[0] in ["H","h","-H","-h","help","--help","-help"]):
    print "usage:   %s [parameter=value ...]" % command_name
    print "example: %s distl.port=8125 distl.processors=8"%command_name
    print help_str
    return
  else:  print help_str

  phil_objects = []
  argument_interpreter = libtbx.phil.command_line.argument_interpreter(
    master_phil=master_params, home_scope="distl")
  for arg in args:
      try: command_line_params = argument_interpreter.process(arg=arg)
      except KeyboardInterrupt: raise
      except: raise Sorry("Unknown file or keyword: %s" % arg)
      else: phil_objects.append(command_line_params)

  working_params = master_params.fetch(sources=phil_objects)
  params = working_params.extract()

  working_params = master_params.format(python_object=params)
  working_params.show()

  #Now actually run the program logic
  from spotfinder.servers import mp_spotfinder_server_read_file as srv

  NUMBER_OF_PROCESSES = params.distl.processors

  srv.common_parameters(outer_resolution=params.distl.res.outer,
                    minimum_spot_area=params.distl.minimum_spot_area,
                    minimum_signal_height=params.distl.minimum_signal_height)

  server_address = ('', params.distl.port)

  srv.image_request_handler.protocol_version = "HTTP/1.0"

  print "Serving %d-process HTTP on"%NUMBER_OF_PROCESSES, server_address[0], "port", server_address[1], "..."
  print "To exit press Ctrl-C"
  srv.runpool(server_address,NUMBER_OF_PROCESSES,srv.image_request_handler)

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
