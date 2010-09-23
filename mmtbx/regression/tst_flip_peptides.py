
import libtbx.load_env
from libtbx import easy_run
from libtbx.test_utils import contains_lines, approx_equal
import os

def exercise () :
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/3ifk.pdb",
    test=os.path.isfile)
  mtz_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/reflection_files/3ifk.mtz",
    test=os.path.isfile)
  if (pdb_file is None) or (mtz_file is None) :
    print "Skipping tst_flip_peptides: input file(s) not available."
    return
  cmd = " ".join([ "mmtbx.flip_peptides", pdb_file, mtz_file ])
  output = easy_run.fully_buffered(cmd).stdout_lines
  output_str = "\n".join(output)
  print output_str
  assert contains_lines(output_str, "wrote new PDB file to 3ifk_new.pdb")
  residues = set([])
  for line in output :
    if line.startswith("  \"GLN B   3 \"") :
      assert (line.endswith("deg") or line.endswith("<<<"))
    elif line.startswith("  \"LEU A   4 \"") :
      assert (line.endswith("deg") or not line.endswith("<<<"))
      residues.add("LEU A   4")
    elif line.startswith("  \"LEU B   4 \"") :
      assert (line.endswith("deg") or not line.endswith("<<<"))
      residues.add("LEU B   4")
  assert (residues == set(["LEU A   4", "LEU B   4"]))
  from iotbx import file_reader
  pdb_old = file_reader.any_file(pdb_file).file_object
  pdb_new = file_reader.any_file("3ifk_new.pdb").file_object
  xrs = pdb_old.xray_structure_simple()
  xrs_shifted = pdb_new.xray_structure_simple()
  distances = xrs.distances(xrs_shifted)
  mmm = distances.min_max_mean()
  assert approx_equal(mmm.max, 3.42, 0.1)
  hierarchy_old = pdb_old.construct_hierarchy()
  cache = hierarchy_old.atom_selection_cache()
  isel = cache.iselection("chain A and resseq 4 and name O")
  assert approx_equal(distances[12], 3.42, 0.1)

if __name__ == "__main__" :
  exercise()
  print "OK"
