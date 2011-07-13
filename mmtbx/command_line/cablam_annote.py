# LIBTBX_SET_DISPATCHER_NAME phenix.cablam_annote

from iotbx import pdb
import libtbx.phil
from libtbx.utils import Usage
import sys, os

master_phil = libtbx.phil.parse("""
cablam {
  file_name = None
    .type = path
    .multiple = true
    .short_caption = PDB file
    .help = '''input coordinate file (pdb)'''
  kin = False
    .type = bool
    .short_caption = generate kinemage output
    .help = '''generate kinemage output'''
}
""", process_includes=True)

def run(args):
  from mmtbx.cablam import cablam_annote, cablam_math, cablam_res
  from iotbx import file_reader
  interpreter = master_phil.command_line_argument_interpreter()
  pdb_files = []
  sources = []
  for arg in args :
    if os.path.isfile(arg) :
      input_file = file_reader.any_file(arg)
      if (input_file.file_type == "pdb") :
        pdb_files.append(input_file)
        sources.append(interpreter.process(arg="file_name=\"%s\"" % arg))
    else :
      arg_phil = interpreter.process(arg=arg)
      sources.append(arg_phil)
  work_phil = master_phil.fetch(sources=sources)
  work_params = work_phil.extract()
  if len(work_params.cablam.file_name) == 0 :
    if len(pdb_files) == 0 :
      summary = cablam_annote.get_summary()
      raise Usage(summary)
    else :
      for pdb_file in pdb_files:
        work_params.cablam.file_name.append(pdb_file.file_name)
  params = work_params.cablam
  #parser = argparse.ArgumentParser()
  #parser.add_argument('file_or_dir',
  #  help='the path to the file or directory to be operated on')
  #parser.add_argument('--kin', action='store_true',
  #  help='flag to get kinemage output')

  #args = parser.parse_args()

  #if os.path.isdir(args.file_or_dir):
  #  fileset = os.listdir(args.file_or_dir)
  #  dirpath = args.file_or_dir
  #elif os.path.isfile(args.file_or_dir):
  #  fileset = [args.file_or_dir]
  #  dirpath = None
  #else:
  #  sys.stderr.write("Could not identify valid target file or dir.\n")
  #  ## print the help section
  #  sys.exit()

  #for filename in fileset:
  #  if dirpath: #must add the path if using the listed contents of a dir
  #    filename = os.path.join(dirpath,filename)
  #  else:
  #    pass
  for filename in params.file_name:
    pdb_io = pdb.input(filename)
    pdbid = os.path.basename(filename)

    sys.stderr.write(pdbid+'\n')
    hierarchy = pdb_io.construct_hierarchy()
    resdata = cablam_res.construct_linked_residues(
      hierarchy, targetatoms=["CA","O","C","N"], pdbid=pdbid)

    cablam_math.CApseudos(resdata, dodihedrals=True, doangles=True)

    dsspcodes = ['H','G','I','E','B','S','T','X']
    for dsspcode in dsspcodes:
      cablam_annote.annote_dssp_3d(resdata,dsspcode)

    reskeys = resdata.keys()
    reskeys.sort()

    if resdata:
      if params.kin:
        cablam_annote.print_annoted_kin(resdata)
      else:
        cablam_annote.print_annoted_text_human(resdata)

#{{{ __main__
if __name__ == "__main__":
  run(sys.argv[1:])
#}}}
