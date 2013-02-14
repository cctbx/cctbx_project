
from __future__ import division
from libtbx.utils import Sorry, Usage
import libtbx.phil
import os
import sys

master_phil = libtbx.phil.parse("""
resname = None
  .type = str
d_max = None
  .type = float
polymeric_type = *Any Free Polymeric
  .type = choice
xray_only = True
  .type = bool
data_only = False
  .type = bool
""")

def run (args, out=sys.stdout) :
  if (len(args) == 0) or ("--help" in args) :
    raise Usage("""mmtbx.find_residue_in_pdb RESNAME [options]

Use the RCSB web services to retrieve a list of PDB structures containing the
specified chemical ID.

Full parameters:
%s
""" % master_phil.as_str(prefix="  "))
  sources = []
  for arg in args :
    if (1 <= len(arg) <= 3) and (arg.isalnum()) :
      sources.append(libtbx.phil.parse("resname=%s" % arg))
    elif os.path.isfile(arg) :
      try :
        sources.append(libtbx.phil.parse(file_name=arg))
      except RuntimeError, e :
        raise Sorry("Can't parse %s as a parameter file:\n%s" % (arg, str(e)))
    else :
      try :
        sources.append(libtbx.phil.parse(arg))
      except RuntimeError, e :
        raise Sorry("Unrecognized argument '%s'." % arg)
  params = master_phil.fetch(sources=sources).extract()
  if (params.resname is None) :
    raise Sorry("No residue ID specified.")
  from mmtbx.wwpdb import rcsb_web_services
  pdb_ids = rcsb_web_services.chemical_id_search(
    resname=params.resname,
    d_max=params.d_max,
    polymeric_type=params.polymeric_type,
    xray_only=params.xray_only,
    data_only=params.data_only)
  pdb_ids = [ id.lower() for id in pdb_ids ]
  if (len(pdb_ids) == 0) :
    raise Sorry("No structures found matching the specified criteria.")
  else :
    print >> out, "%d PDB IDs retrieved:" % len(pdb_ids)
    i = 0
    while (i < len(pdb_ids)) :
      print >> out, "  %s" % " ".join(pdb_ids[i:i+16])
      i += 16

if (__name__ == "__main__") :
  run(sys.argv[1:])
