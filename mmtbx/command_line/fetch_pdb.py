from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.fetch_pdb

# XXX most of this code is unused when run from the command line, but the
# PHENIX GUI includes a simple frontend that uses the phil interface.

import libtbx.phil
from libtbx.utils import Sorry, Usage
from libtbx import easy_run
import sys
import os

master_phil = libtbx.phil.parse("""
fetch_pdb
  .caption = PHENIX can automatically retrieve data from the PDB via the RCSB \
    web server.  If you intend to re-refine or re-build the structure we \
    recommend creating a new project, but this is not required.  Note that \
    you may also use this tool to generate an MTZ file from the mmCIF \
    structure factors (if available) and even calculate maps, but the options \
    are more limited than what is available in the phenix.cif_as_mtz and \
    phenix.maps GUIs.
  .style = auto_align caption_img:icons/custom/pdb_import64.png \
    caption_width:400
{
  pdb_ids = None
    .type = strings
    .short_caption = PDB ID(s)
    .input_size = 400
    .style = bold
  action = *pdb_only all_data all_plus_mtz all_plus_maps
    .type = choice
    .caption = Download_PDB_file(s) Download_all_data Download_all_data_and_convert_CIF_to_MTZ Download_all_data_and_create_maps
    .style = bold
  site = *rcsb pdbe pdbj
    .type = choice
    .caption = RCSB PDBe PDBj
    .short_caption = Mirror site
    .style = bold
}""")

# XXX for use in the PHENIX GUI only - run2 is used from the command line
def run (args=(), params=None, out=sys.stdout) :
  assert (params is not None)
  output_files = []
  errors = []
  mirror = "--mirror=%s" % params.fetch_pdb.site
  for id in params.fetch_pdb.pdb_ids :
    args = [mirror, id]
    if (params.fetch_pdb.action in ["all_data","all_plus_mtz",
        "all_plus_maps"]) :
      args.insert(0, "--all")
      if (params.fetch_pdb.action == "all_plus_mtz") :
        args.insert(1, "--mtz")
      elif (params.fetch_pdb.action == "all_plus_maps") :
        args.insert(1, "--maps")
      try :
        data_files = run2(args=args, log=out)
        print >> out, "\n".join(data_files)
        output_files.extend(data_files)
      except Exception, e :
        errors.append(str(e))
    else :
      pdb_file = run2(args=[mirror,id], log=out)
      print >> out, pdb_file
      output_files.append(pdb_file)
  return output_files, errors

def run2 (args, log=sys.stdout) :
  if len(args) < 1 :
    raise Usage("""\
phenix.fetch_pdb [-x|-f|--all] [--mtz] [--maps] [-q] ID1 [ID2, ...]

Command-line options:
  -x      Get structure factors (mmCIF file)
  -c      Get model file in mmCIF format
  -f      Get sequence (FASTA file)
  --all   Download all available data
  --mtz   Download structure factors and PDB file, and generate MTZ
  --maps  As for --mtz, plus create 2mFo-DFc and mFo-DFc map coefficients (and
          anomalous map, if possible)
  --fill  Fill missing F(obs) with F(calc) in map coefficients
  -q      suppress printed output
""")
  from iotbx.pdb.fetch import get_pdb
  quiet = False
  convert_to_mtz = maps = fill_maps = False
  data_type = "pdb"
  format = "pdb"
  mirror = "rcsb"
  ids = []
  for arg in args :
    if (arg == "--all") :
      data_type = "all"
    elif (arg == "-x") :
      data_type = "xray"
    elif (arg == "-f") :
      data_type = "fasta"
    elif (arg == "-q") :
      quiet = True
    elif (arg == "--mtz") :
      convert_to_mtz = True
      data_type = "all"
    elif (arg == "--maps") :
      convert_to_mtz = True
      data_type = "all"
      maps = True
    elif (arg == "-c") :
      format = "cif"
    elif (arg.startswith("--mirror=")) :
      mirror = arg.split("=")[1]
      if (not mirror in ["rcsb", "pdbe", "pdbj"]) :
        raise Sorry(
          "Unrecognized mirror site '%s' (choices: rcsb, pdbe, pdbj)" %
          mirror)
    else :
      ids.append(arg)
  if (len(ids) == 0) :
    raise Sorry("No PDB IDs specified.")
  if (data_type != "all") :
    #mirror = "rcsb"
    files = []
    for id in ids :
      files.append(get_pdb(id, data_type, mirror, log, format=format))
    if (len(files) == 1) :
      return files[0]
    return files
  else :
    files = []
    for id in ids :
      for data_type_, data_format in [("pdb", "pdb"), ("fasta", "pdb"),
                                      ("xray", "pdb"), ("pdb", "cif")] :
        files.append(get_pdb(id, data_type_, mirror, log, format=data_format))
      if (convert_to_mtz) :
        misc_args = ["--merge", "--map_to_asu", "--extend_flags",
                     "--ignore_bad_sigmas"]
        easy_run.call("phenix.cif_as_mtz %s-sf.cif --symmetry=%s.pdb %s" %
          (id,id, " ".join(misc_args)))
        if os.path.isfile("%s-sf.mtz" % id) :
          os.rename("%s-sf.mtz" % id, "%s.mtz" % id)
          print >> log, "Converted structure factors saved to %s.mtz" % id
        #  os.remove("%s-sf.cif" % id)
        files[-1] = os.path.abspath("%s.mtz" % id)
        if (not os.path.isfile("%s.mtz" % id)) :
          raise Sorry("MTZ conversion failed - try running phenix.cif_as_mtz "+
            "manually (and check %s-sf.cif for format errors)." % id)
        from iotbx.file_reader import any_file
        mtz_in = any_file("%s.mtz" % id)
        mtz_in.assert_file_type("hkl")
        for array in mtz_in.file_server.miller_arrays :
          if (array.anomalous_flag()) :
            print >> log, "  %s is anomalous" % array.info().label_string()
      if (maps) :
        assert os.path.isfile("%s.mtz" % id)
        from mmtbx.maps.utils import create_map_from_pdb_and_mtz
        create_map_from_pdb_and_mtz(
          pdb_file="%s.pdb" % id,
          mtz_file="%s.mtz" % id,
          output_file="%s_maps.mtz" % id,
          fill=fill_maps,
          out=log,
          remove_unknown_scatterering_type=True,
          assume_pdb_data=True)
        assert os.path.isfile("%s_maps.mtz" % id)
        print >> log, "Map coefficients saved to %s_maps.mtz" % id
        files.append(os.path.abspath("%s_maps.mtz" % id))
    return files

def validate_params (params) :
  if (params.fetch_pdb.pdb_ids is None) or (len(params.fetch_pdb.pdb_ids)==0) :
    raise Sorry("No PDB IDs specified!")
  return True

if __name__ == "__main__" :
  run2(sys.argv[1:])
