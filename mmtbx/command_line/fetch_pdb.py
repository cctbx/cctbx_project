"""
Provides a command-line utility for fetching PDB files and their associated
reflection data.
"""

from __future__ import absolute_import, division, print_function

# LIBTBX_SET_DISPATCHER_NAME phenix.fetch_pdb
# LIBTBX_SET_DISPATCHER_NAME iotbx.fetch_pdb

# XXX most of this code is unused when run from the command line, but the
# PHENIX GUI includes a simple frontend that uses the phil interface.

import libtbx.phil
from libtbx.utils import Sorry, Usage
from libtbx import easy_run
from iotbx.pdb.fetch import get_pdb
import sys
import os

master_phil = libtbx.phil.parse("""
fetch_pdb
  .caption = PHENIX can automatically retrieve data from the PDB via the RCSB \
    web server.  If you intend to re-refine or re-build the structure we \
    recommend creating a new project, but this is not required.  Note that \
    you may also use this tool to generate an MTZ file from the mmCIF \
    structure factors (if available), but the options \
    are more limited than what is available in the phenix.cif_as_mtz.
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
  site = *rcsb pdbe pdbj
    .type = choice
    .caption = RCSB PDBe PDBj
    .short_caption = Mirror site
    .style = bold
}""")

# XXX for use in the PHENIX GUI only - run2 is used from the command line
def run(args=(), params=None, out=sys.stdout):
  """
  For use in PHENIX GUI only, fetches pdb filesand/or
  reflection data from the PDB.

  Parameters
  ----------
  args : list of str, optional
  params : libtbx.phil.scope_extract, optional
  out : file, optional

  Returns
  -------
  output_files : list of str
  errors : list of str
  """
  assert (params is not None)
  output_files = []
  errors = []
  mirror = "--mirror=%s" % params.fetch_pdb.site
  for id in params.fetch_pdb.pdb_ids :
    args = [mirror, id]
    if (params.fetch_pdb.action in ["all_data","all_plus_mtz"]):
      args.insert(0, "--all")
      if (params.fetch_pdb.action == "all_plus_mtz"):
        args.insert(1, "--mtz")
      try :
        data_files = run2(args=args, log=out)
        print("\n".join(data_files), file=out)
        output_files.extend(data_files)
      except Exception as e:
        errors.append(str(e))
    else :
      pdb_file = run2(args=[mirror,id], log=out)
      print(pdb_file, file=out)
      output_files.append(pdb_file)
  return output_files, errors

def run2(args, log=sys.stdout):
  """
  Fetches pdb files and/or reflection data from the PDB.

  Parameters
  ----------
  args : list of str
  log : file, optional

  Returns
  -------
  str or list of str
      List of file names that were downloaded.
  """
  if len(args) < 1 :
    raise Usage("""\
phenix.fetch_pdb [-x|-f|--all] [--mtz] [-q] ID1 [ID2, ...]

Command-line options:
  -x      Get structure factors (mmCIF file)
  -c      Get model file in mmCIF format
  -f      Get sequence (FASTA file)
  --all   Download all available data
  --mtz   Download structure factors and PDB file, and generate MTZ
  -q      suppress printed output
""")
  quiet = False
  convert_to_mtz  = False
  data_type = "pdb"
  format = "pdb"
  mirror = "rcsb"
  ids = []
  for arg in args :
    if (arg == "--all"):
      data_type = "all"
    elif (arg == "-x"):
      data_type = "xray"
    elif (arg == "-f"):
      data_type = "fasta"
    elif (arg == "-q"):
      quiet = True
    elif (arg == "--mtz"):
      convert_to_mtz = True
      data_type = "all"
    elif (arg == "-c"):
      format = "cif"
    elif (arg.startswith("--mirror=")):
      mirror = arg.split("=")[1]
      if (not mirror in ["rcsb", "pdbe", "pdbj", "pdb-redo"]):
        raise Sorry(
          "Unrecognized mirror site '%s' (choices: rcsb, pdbe, pdbj, pdb-redo)" %
          mirror)
    else :
      ids.append(arg)
  if (len(ids) == 0):
    raise Sorry("No PDB IDs specified.")
  if (data_type != "all"):
    #mirror = "rcsb"
    files = []
    for id in ids :
      try:
        files.append(get_pdb(id, data_type, mirror, log, format=format))
      except Sorry as e:
        if data_type == "pdb" and format == "pdb":
          # try with mmCIF
          print("  PDB format is not available, trying to get mmCIF", file=log)
          files.append(get_pdb(id, data_type, mirror, log, format="cif"))
    if (len(files) == 1):
      return files[0]
    return files
  else :
    files = []
    for id in ids :
      for data_type_, data_format in [("pdb", "pdb"), ("fasta", "pdb"),
                                      ("xray", "pdb"), ("pdb", "cif")] :
        try:
          files.append(get_pdb(id, data_type_, mirror, log, format=data_format))
        except Sorry as e:
          print (str(e))
      if (convert_to_mtz):
        misc_args = ["--merge", "--map_to_asu", "--extend_flags",
                     "--ignore_bad_sigmas"]
        easy_run.call("phenix.cif_as_mtz %s-sf.cif %s" %
          (id, " ".join(misc_args)))
        if os.path.isfile("%s-sf.mtz" % id):
          os.rename("%s-sf.mtz" % id, "%s.mtz" % id)
          print("Converted structure factors saved to %s.mtz" % id, file=log)
        #  os.remove("%s-sf.cif" % id)
        files[-1] = os.path.abspath("%s.mtz" % id)
        if (not os.path.isfile("%s.mtz" % id)):
          raise Sorry("MTZ conversion failed - try running phenix.cif_as_mtz "+
            "manually (and check %s-sf.cif for format errors)." % id)
        from iotbx.file_reader import any_file
        mtz_in = any_file("%s.mtz" % id)
        mtz_in.assert_file_type("hkl")
        for array in mtz_in.file_server.miller_arrays :
          if (array.anomalous_flag()):
            print("  %s is anomalous" % array.info().label_string(), file=log)
    return files

def validate_params(params):
  """
  Validates that the input parameters specify a PDB file to download.

  Parameters
  ----------
  params : libtbx.phil.scope_extract

  Returns
  -------
  bool

  Raises
  ------
  Sorry
      Raised if a pdb file is not specified to be fetched within params.
  """
  if (params.fetch_pdb.pdb_ids is None) or (len(params.fetch_pdb.pdb_ids)==0):
    raise Sorry("No PDB IDs specified!")
  return True

if __name__ == "__main__" :
  run2(sys.argv[1:])
