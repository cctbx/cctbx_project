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
  action = *pdb_only all_data
    .type = choice
    .caption = Download_PDB_file(s) Download_all_data
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
    pdb_file = iotbx.pdb.fetch.run(args=[id])
    output_files.append(pdb_file)
    if (params.fetch_pdb.action == "all_data") :
      for file_flag in ["-x", "-f"] :
        try :
          data_file = iotbx.pdb.fetch.run(args=[file_flag,id,mirror])
          print data_file
          output_files.append(data_file)
        except Exception, e :
          errors.append(str(e))
  return output_files, errors

def validate_params (params) :
  if (params.fetch_pdb.pdb_ids is None) or (len(params.fetch_pdb.pdb_ids)==0) :
    raise Sorry("No PDB IDs specified!")
  return True

if __name__ == "__main__" :
  import iotbx.pdb.fetch
  iotbx.pdb.fetch.run(sys.argv[1:])
