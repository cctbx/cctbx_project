from __future__ import absolute_import, division, print_function
from libtbx.utils import Sorry, Usage
import libtbx.phil.command_line
import sys

master_phil = libtbx.phil.parse("""
resname = None
  .type = str
d_max = None
  .type = float
protein_only = False
  .type = bool
xray_only = True
  .type = bool
data_only = False
  .type = bool
quiet = False
  .type = bool
""")

def run(args, out=sys.stdout):
  if (len(args) == 0) or ("--help" in args):
    raise Usage("""mmtbx.find_residue_in_pdb RESNAME [options]

Use the RCSB web services to retrieve a list of PDB structures containing the
specified chemical ID.

Full parameters:
%s
""" % master_phil.as_str(prefix="  "))
  sources = []
  def process_unknown(arg):
    if (1 <= len(arg) <= 3) and (arg.isalnum()):
      return libtbx.phil.parse("resname=%s" % arg)
  cai = libtbx.phil.command_line.argument_interpreter(master_phil=master_phil)
  working_phil = cai.process_and_fetch(args=args,
    custom_processor=process_unknown)
  params = working_phil.extract()
  if (params.resname is None):
    raise Sorry("No residue ID specified.")
  from mmtbx.wwpdb import rcsb_web_services
  pdb_ids = rcsb_web_services.chemical_id_search(
    resname=params.resname,
    d_max=params.d_max,
    protein_only=params.protein_only,
    xray_only=params.xray_only,
    data_only=params.data_only,
    )
  pdb_ids = [ id.lower() for id in pdb_ids ]
  if (len(pdb_ids) == 0):
    raise Sorry("No structures found matching the specified criteria.")
  else :
    if (not params.quiet):
      print("%d PDB IDs retrieved:" % len(pdb_ids), file=out)
      i = 0
      while (i < len(pdb_ids)):
        print("  %s" % " ".join(pdb_ids[i:i+16]), file=out)
        i += 16
    else :
      print("%d PDB IDs matching" % len(pdb_ids), file=out)

if (__name__ == "__main__"):
  run(sys.argv[1:])
