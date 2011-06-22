# LIBTBX_SET_DISPATCHER_NAME distl.mp_spotfinder_server_read_file
from libtbx.utils import Sorry
from spotfinder.command_line.signal_strength import master_params

def run(args, command_name="distl.mp_spotfinder_server_read_file"):
  help_str="""Multiprocessing server to find Bragg spots & quantify signal strength.
Full documentation: http://cci.lbl.gov/publications/download/ccn_jul2010_page18.pdf
Allowed parameters:
"""

  if (len(args)>=1 and args[0] in ["H","h","-H","-h","help","--help","-help"]):
    print "usage:   %s [parameter=value ...]" % command_name
    print "example: %s distl.port=8125 distl.processors=8"%command_name
    print help_str
    return
  else:  print help_str

  phil_objects = []
  argument_interpreter = master_params.command_line_argument_interpreter(
    home_scope="distl")
  for arg in args:
      try: command_line_params = argument_interpreter.process(arg=arg)
      except KeyboardInterrupt: raise
      except Exception: raise Sorry("Unknown file or keyword: %s" % arg)
      else: phil_objects.append(command_line_params)

  working_params = master_params.fetch(sources=phil_objects)
  params = working_params.extract()

  working_params = master_params.format(python_object=params)
  #working_params.show()

  screen = 100
  for D in working_params.all_definitions():
    fp = D.object.full_path()
    if fp in ["distl.image", "distl.verbose", "distl.pdf_output"]: continue
    name = "  %s=%s"%(D.object.full_path(),D.object.extract())
    help = D.object.help
    if len(name) + len(help) < screen and len(name) < 36:
        print "%-36s"%name,"[%s]"%help
    else:
      print "%-36s"%name,
      tokens = ("[%s]"%help).split()
      reserve = min(screen-36,screen-len(name))
      while len(tokens)>0:
        reserve -= len(tokens[0])
        if len(tokens[0])>(screen-36):break
        while reserve<0:
          print
          print " "*36,
          reserve = screen-36
        print tokens.pop(0),
      print

  #Now actually run the program logic
  from spotfinder.servers import mp_spotfinder_server_read_file as srv

  NUMBER_OF_PROCESSES = params.distl.processors

  from labelit.preferences import labelit_commands
  labelit_commands.distl = params.distl
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
