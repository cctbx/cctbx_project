from libtbx import test_utils
import libtbx.load_env

def run():
  tst_list = (
  "$D/tst_wildcard.py",
  "$D/tst_simple_parser.py",
  "$D/tst_phil.py",
  "$D/ranges.py",
  "$D/tst_crystal_symmetry_from_any.py",
  "$D/kriber/tst_strudat.py",
  "$D/cif/tests/tst_geometry.py",
  "$D/cif/tests/tst_lex_parse_build.py",
  "$D/cif/tests/tst_model.py",
  "$D/cif/tests/tst_restraints.py",
  "$D/cif/tests/tst_validation.py",
  "$D/shelx/tst_lex_parse_build.py",
  "$D/shelx/tst_hklf.py",
  "$D/shelx/tst_writer.py",
  "$D/shelx/tst_fvar_encoding.py",
  ["$D/pdb/hybrid_36.py", "exercise"],
  "$B/pdb/hybrid_36_fem",
  "$D/pdb/tst_hierarchy.py",
  "$D/pdb/tst_ext.py",
  "$D/pdb/tst_atom_selection.py",
  "$D/pdb/tst_rna_dna_atom_names.py",
  "$D/pdb/tst_atom_name_interpretation.py",
  "$D/pdb/tst_extract_rfactors_resolutions_sigma.py",
  "$D/pdb/tst_pdb.py",
  "$D/pdb/secondary_structure.py",
  "$D/examples/pdb_to_map_simple.py",
  "$D/examples/pdb_truncate_to_ala/tst.py",
  "$D/examples/pdb_tardy_conf_sampling_simple.py",
  "$D/cns/space_group_symbols.py",
  "$D/cns/tst_cns.py",
  ["$D/scalepack/tst_merge.py", "P31"],
  "$D/scalepack/no_merge_original_index.py",
  "$D/ccp4_map/tst.py",
  "$D/mtz/tst_ext.py",
  "$D/mtz/tst_extract_from_symmetry_lib.py",
  ["$D/mtz/tst.py", "P31"],
  "$D/examples/tst_mtz_free_flipper.py",
  "$D/tst_reflection_file_utils.py",
  "$D/detectors/tst_adsc.py",
  "$D/detectors/tst_debug_write.py",
  "$D/xplor/tst_xplormap.py",
  ["$D/tst_phases.py", "P31"],
  "$D/regression/tst_lattice_symmetry.py",
  ["$D/regression/tst_reflection_statistics.py", "Fdd2 P31m"],
  "$D/tst_data_plots.py",
  "$D/tst_csv_utils.py",
  "$D/tst_file_reader.py",
  "$D/tst_bioinformatics.py",
  "$D/regression/tst_add_conformations.py",
  "$D/regression/tst_symmetry.py",
  )

  build_dir = libtbx.env.under_build("iotbx")
  dist_dir = libtbx.env.dist_path("iotbx")

  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run()
