# LIBTBX_SET_DISPATCHER_NAME iotbx.pdb_remediator

import libtbx.phil
from libtbx.utils import Usage, Sorry
import os
import sys

master_phil = libtbx.phil.parse("""
remediator
  .short_caption = PDB remediation
  .caption = This utility converts between version 2.3 and version 3.2 of \
    the PDB format, which differ primarily in the treatment of hydrogen \
    atoms.  PHENIX always uses v3.2 atom names, as do recent \
    versions of Coot, Refmac, and other programs using the CCP4 monomer \
    library, but it can also interpret previous formats.  You may however \
    need to convert between formats for certain applications.  (The PDB will \
    always convert to v3.2 upon submission.)
  .style = box auto_align caption_img:icons/custom/phenix.pdbtools.png
{
  file_name = None
    .type = path
    .short_caption = PDB file
    .help = '''input coordinate file (pdb)'''
    .style = bold file_type:pdb input_file
  output_file = None
    .type = path
    .optional = True
    .help = '''Enter a .pdb output name'''
    .style = new_file file_type:pdb
  version = *3.2 2.3
    .type = choice(multi=False)
    .short_caption = Format version
  dict = None
    .type = path
    .optional = True
    .short_caption = Custom definition file
    .help = '''custom definition file'''
}
""", process_includes=True)

def run(args=(), params=None, out=sys.stdout):
  from iotbx.pdb.remediation import remediator
  from iotbx import file_reader
  if (params is None) :
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
  else : # XXX for phenix GUI
    work_params = params
    if (work_params.remediator.output_file is None) :
      base, ext = os.path.splitext(work_params.remediator.file_name)
      work_params.remediator.output_file = base + "_remediated.pdb"
  if (work_params.remediator.file_name is None) :
    if (pdb_file is None) :
      summary = remediator.get_summary()
      raise Usage(summary)
    else :
      work_params.remediator.file_name = pdb_file.file_name
  params = work_params.remediator
  remediator.remediator(params)
  return work_params.remediator.output_file

def validate_params (params) :
  if (params.remediator.file_name is None) :
    raise Sorry("Please specify a PDB file.")
  return True

if(__name__ == "__main__"):
  run(sys.argv[1:])
