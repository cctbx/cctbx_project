from libtbx import test_utils
import libtbx.load_env

def run():
  tst_list = (
  "$D/iotbx/tst_wildcard.py",
  "$D/iotbx/tst_simple_parser.py",
  "$D/iotbx/tst_phil.py",
  "$D/iotbx/tst_crystal_symmetry_from_any.py",
  "$D/iotbx/kriber/tst_strudat.py",
  ["$D/iotbx/pdb/hybrid_36.py", "exercise"],
  "$D/include/iotbx/pdb/tst_hierarchy.py",
  "$D/include/iotbx/pdb/tst_ext.py",
  "$D/iotbx/pdb/tst_atom_selection.py",
  "$D/include/iotbx/pdb/tst_rna_dna_atom_names.py",
  "$D/iotbx/pdb/tst_atom_name_interpretation.py",
  "$D/iotbx/pdb/tst_extract_rfactors_resolutions_sigma.py",
  "$D/iotbx/pdb/tst_pdb.py",
  "$D/iotbx/examples/pdb_truncate_to_ala/tst.py",
  "$D/iotbx/cns/space_group_symbols.py",
  "$D/iotbx/cns/tst_cns.py",
  ["$D/iotbx/scalepack/tst_merge.py", "P31"],
  "$D/iotbx/scalepack/no_merge_original_index.py",
  "$D/include/iotbx/mtz/tst_ext.py",
  "$D/iotbx/mtz/extract_from_symop_lib.py",
  ["$D/iotbx/mtz/tst.py", "P31"],
  "$D/iotbx/examples/tst_mtz_free_flipper.py",
  "$D/iotbx/tst_reflection_file_utils.py",
  "$D/iotbx/detectors/tst_adsc.py",
  "$D/iotbx/detectors/tst_debug_write.py",
  "$D/iotbx/xplor/tst_xplormap.py",
  ["$D/iotbx/tst_phases.py", "P31"],
  ["$D/iotbx/regression/tst_reflection_statistics.py", "Fdd2 P31m"]
  )

  build_dir = libtbx.env.under_build("iotbx")
  dist_dir = libtbx.env.dist_path("iotbx")

  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run()
