# LIBTBX_SET_DISPATCHER_NAME iotbx.remediator

import libtbx.phil
from libtbx.utils import Usage
import os
import sys

master_phil = libtbx.phil.parse("""
remediator {
  file_name = None
    .type = path
    .short_caption = PDB file
    .help = '''input coordinate file (pdb)'''
  output_file = None
    .type = path
    .optional = True
    .help = '''Enter a .pdb output name'''
  version = *3.2 2.3
    .type = choice(multi=False)
  dict = None
    .type = path
    .optional = True
    .help = '''custom definition file'''
}
""", process_includes=True)

def run(args):
  from iotbx.remediation import remediator
  from iotbx import file_reader
  interpreter = master_phil.command_line_argument_interpreter()
  pdb_file = None
  sources = []
  for arg in args :
    if os.path.isfile(arg) :
      input_file = file_reader.any_file(arg)
      if (input_file.file_type == "pdb") :
        pdb_file = input_file
        sources.append(interpreter.process(arg="file_name=\"%s\"" % arg))
    else :
      arg_phil = interpreter.process(arg=arg)
      sources.append(arg_phil)
  work_phil = master_phil.fetch(sources=sources)
  work_params = work_phil.extract()
  if (work_params.remediator.file_name is None) :
    if (pdb_file is None) :
      summary = remediator.get_summary()
      raise Usage(summary)
    else :
      work_params.remediator.file_name = pdb_file.file_name
  params = work_params.remediator
  remediator.remediator(params)

if(__name__ == "__main__"):
  run(sys.argv[1:])
