
import iotbx.phil
from libtbx.utils import Sorry, Usage
from libtbx import group_args
import sys

master_phil = iotbx.phil.parse("""
symmetry_search
  .short_caption = PDB symmetry search
  .caption = This utility allows you to search the PDB for structures with \
    similar unit cell parameters.  Crystallization artifacts due to \
    impurities in the protein solution can often be detected this way, if the \
    protein which actually crystallized has been solved before.  Note that \
    a large number of false positives are usually expected for genuinely \
    novel structures, so the presence of similar unit cells is not \
    necessarily a bad sign.
  .style = box auto_align caption_img:icons/custom/pdb_import64.png
{
  file_name = None
    .type = path
    .short_caption = File name (PDB or MTZ)
  unit_cell = None
    .type = unit_cell
  space_group = None
    .type = space_group
  max_rmsd = None
    .type = float
    .short_caption = Max. RMSD to consider
  max_hits_to_display = 50
    .type = int
    .short_caption = Max. number of hits
}
""")

def run (args=(), params=None, out=sys.stdout) :
  if (len(args) > 0) :
    cmdline = iotbx.phil.process_command_line_with_files(
      args=args,
      master_phil=master_phil,
      pdb_file_def="symmetry_search.file_name",
      reflection_file_def="symmetry_search.file_name")
    params = cmdline.work.extract().symmetry_search
  elif (params is None) :
    raise Usage("""mmtbx.search_pdb_symmetry [file] [space_group] [unit_cell]
  Utility for finding similar unit cells deposited in the PDB.
""")
  else :
    params = params.symmetry_search
  from mmtbx import pdb_symmetry
  from iotbx import crystal_symmetry_from_any
  from cctbx import crystal
  db = pdb_symmetry.load_db()
  if (params.file_name is not None) :
    symm = crystal_symmetry_from_any.extract_from(file_name=params.file_name)
    if (symm is None) :
      raise Sorry("The file %s does not include symmetry information." %
        params.file_name)
    elif (symm.space_group() is None) or (symm.unit_cell() is None) :
      raise Sorry("Incomplete symmetry information in %s." % params.file_name)
  else :
    symm = crystal.symmetry(
      unit_cell=params.unit_cell,
      space_group_info=params.space_group)
  print >> out, ""
  print >> out, "Input symmetry:"
  symm.show_summary()
  scores = pdb_symmetry.symmetry_search(db, symm, max_rmsd=params.max_rmsd)
  niggli_cell = symm.niggli_cell().unit_cell().parameters()
  print >> out, ""
  print >> out, "Top %d matches (sorted by RMSD):"
  results = []
  for scored in scores[:params.max_hits_to_display] :
    print >> out, "%s (rmsd = %.3f, volume ratio = %.2f)" % \
      (scored.entry.pdb_id, scored.rmsd, scored.volume_ratio)
    print >> out, "    Unit cell: %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f" % \
      scored.entry.crystal_symmetry.unit_cell().parameters()
    print >> out, "  Niggli cell: %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f" % \
      scored.entry.niggli_cell.unit_cell().parameters()
    print >> out, "  Target cell: %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f" % \
      niggli_cell
    print >> out, ""
    results.append(group_args(
      pdb_id=scored.entry.pdb_id,
      rmsd=scored.rmsd,
      volume_ratio=scored.volume_ratio,
      pdb_symmetry=scored.entry.crystal_symmetry))
  return group_args(
    crystal_symmetry=symm,
    hits=results)

def validate_params (params) :
  params = params.symmetry_search
  have_symm = (not None in [params.unit_cell, params.space_group])
  if (not have_symm) and (params.file_name is None) :
    raise Sorry("Missing or incomplete symmetry information.")
  return True

if (__name__ == "__main__") :
  run(sys.argv[1:])
