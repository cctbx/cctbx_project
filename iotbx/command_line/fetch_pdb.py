# LIBTBX_SET_DISPATCHER_NAME phenix.fetch_pdb

# XXX most of this code is unused when run from the command line, but the
# PHENIX GUI includes a simple frontend that uses the phil interface.

import libtbx.phil
from libtbx.utils import Sorry
import sys

master_phil = libtbx.phil.parse("""
fetch_pdb
  .caption = PHENIX can automatically retrieve data from the PDB via the RCSB \
    web server.  If you intend to re-refine or re-build the structure we \
    recommend creating a new project, but this is not required.  This utility \
    only downloads the data and does no further processing, but the \
    "Import CIF structure factors" module under "Reflection tools" will \
    generate an MTZ file, and can also access the PDB directly.
  .style = auto_align caption_img:icons/custom/pdb_import64.png \
    caption_width:400
{
  pdb_ids = None
    .type = strings
    .short_caption = PDB ID(s)
    .input_size = 400
    .style = bold
  action = *pdb_only all_data all_plus_mtz
    .type = choice
    .caption = Download_PDB_file(s) Download_all_data Download_all_data_and_convert_CIF_to_MTZ
    .style = bold
  site = *rcsb pdbe
    .type = choice
    .caption = RCSB PDBe
    .short_caption = Mirror site
    .style = bold
}""")

def run (args=(), params=None, out=sys.stdout) :
  assert (params is not None)
  import iotbx.pdb.fetch
  output_files = []
  errors = []
  mirror = "--mirror=%s" % params.fetch_pdb.site
  for id in params.fetch_pdb.pdb_ids :
    args = [mirror, id]
    if (params.fetch_pdb.action in ["all_data","all_plus_mtz"]) :
      args.insert(0, "--all")
      if (params.fetch_pdb.action == "all_plus_mtz") :
        args.insert(1, "--mtz")
      try :
        data_files = iotbx.pdb.fetch.run(args=args)
        print "\n".join(data_files)
        output_files.extend(data_files)
      except Exception, e :
        errors.append(str(e))
    else :
      pdb_file = iotbx.pdb.fetch.run(args=[mirror,id])
      print pdb_file
      output_files.append(pdb_file)
  return output_files, errors

def validate_params (params) :
  if (params.fetch_pdb.pdb_ids is None) or (len(params.fetch_pdb.pdb_ids)==0) :
    raise Sorry("No PDB IDs specified!")
  return True

if __name__ == "__main__" :
  import iotbx.pdb.fetch
  iotbx.pdb.fetch.run(sys.argv[1:])
