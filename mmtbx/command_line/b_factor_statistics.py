# LIBTBX_SET_DISPATCHER_NAME phenix.b_factor_statistics

# see also wxtbx/adp_statistics.py

from __future__ import division
from libtbx.utils import Sorry
from libtbx.str_utils import make_sub_header
import os
import sys
import iotbx.cif
import iotbx.pdb

master_phil = """
adp_statistics
  .short_caption = ADP statistics
  .caption = This tool displays tables and plots of statistics for isotropic \
    and anisotropic ADPs in a PDB file.  You can get the same information on \
    the command line by running phenix.pdbtools model_statistics=True \
    model.pdb.
  .style = box auto_align caption_img:icons/custom/phenix.pdbtools.png
{
  pdb_file = None
    .type = path
    .optional = False
    .short_caption = PDB file
    .style = bold file_type:pdb input_file
  cif_file = None
    .type = path
    .multiple = True
    .optional = True
    .short_caption = Restraints (CIF) file
    .style = file_type:cif input_file
  selection = all
    .type = atom_selection
    .short_caption = Atom selection
    # FIXME
    .style = no_view
}
"""

def run (args=(), params=None, out=None) :
  if (out is None) :
    out = sys.stdout
  if (params is None) :
    import iotbx.phil
    cmdline = iotbx.phil.process_command_line_with_files(
      args=args,
      master_phil_string=master_phil,
      pdb_file_def="adp_statistics.pdb_file",
      cif_file_def="adp_statistics.cif_file",
      usage_string="""\
phenix.b_factor_statistics model.pdb [restraints.cif] [selection=...]

Show statistics for atomic displacement parameters (ADPs) or B-factors,
including TLS contribution if present.""")
    params = cmdline.work.extract()
  validate_params(params)
  import mmtbx.model
  m_input = iotbx.pdb.input(file_name=params.adp_statistics.pdb_file)
  restraint_objects = []
  for fn in params.adp_statistics.cif_file:
    cif_object = iotbx.cif.reader(file_path=fn).model()
    restraint_objects.append((fn, cif_object))
  model = mmtbx.model.manager(
      model_input = m_input,
      restraint_objects = restraint_objects,
      build_grm = True,
      log = out)
  make_sub_header("Analyzing model B-factors", out=out)
  if (params.adp_statistics.selection is not None) :
    selection = model.selection(params.adp_statistics.selection)
    n_sel = selection.count(True)
    if (n_sel == 0) :
      raise Sorry("No atoms in selection!")
    else :
      model = model.select(selection)
      print >> out, "Extracted %d atoms in selection:" % n_sel
      print >> out, "  %s" % params.adp_statistics.selection
      print >> out, ""
  stats = model.adp_statistics()
  stats.file_name = params.adp_statistics.pdb_file
  stats.selection = params.adp_statistics.selection
  stats.show_1(out=out)
  return stats

def validate_params (params) :
  if (params.adp_statistics.pdb_file is None) :
    raise Sorry("A PDB file is required to run this program.")
  elif (not os.path.isfile(params.adp_statistics.pdb_file)) :
    raise Sorry("The path %s is not a file or does not exist.")

if (__name__ == "__main__") :
  run(sys.argv[1:])
