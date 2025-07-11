"""
List of tests to run in regression tests
"""

from __future__ import absolute_import, division, print_function
from libtbx import test_utils
import sys
import libtbx.load_env

tst_list_base = [
  "$D/regression/tst_mrc_io.py",
  "$D/regression/tst_wildcard.py",
  "$D/gui_tools/tst.py",
  "$D/regression/tst_simple_parser.py",
  "$D/regression/tst_phil.py",
  "$D/regression/tst_pdb_cif_inputs.py",
  "$D/regression/tst_pdb_cif_cells.py",
  "$D/regression/tst_data_manager.py",
  "$D/regression/tst_map_manager_wrapping.py",
  "$D/regression/tst_map_manager.py",
  "$D/regression/tst_map_manager_2.py",
  "$D/regression/tst_map_model_manager.py",
  "$D/regression/tst_map_model_manager_2.py",
  "$D/regression/tst_map_model_manager_3.py",
  "$D/regression/tst_map_model_manager_4.py",
  "$D/regression/tst_map_model_manager_model_sharpening_5.py",
  "$D/regression/tst_map_model_manager_model_sharpening_5_cif.py",
  "$D/regression/tst_map_model_manager_call_consistency.py",
  "$D/regression/tst_map_model_manager_external_sharpening_7.py",
  "$D/regression/tst_map_model_manager_half_map_sharpening_6.py",
  "$D/regression/tst_map_model_manager_tls_from_map_8.py",
  "$D/regression/tst_map_model_manager_cif.py",
  "$D/regression/tst_map_model_manager_9_remove_origin_shift_and_unit_cell_crystal_symmetry.py",
  "$D/regression/tst_map_model_manager_9_remove_origin_shift_and_unit_cell_crystal_symmetry_cif.py",
  "$D/regression/tst_map_model_manager_local_resolution_10.py",
  "$D/regression/tst_map_tools.py",
  "$D/regression/tst_patterson.py",
  "$D/regression/tst_restraints_merge.py",
  "$D/regression/tst_atom_selections_10k.py",
  "$D/ranges.py",
  "$D/regression/tst_crystal_symmetry_from_any.py",
  "$D/regression/tst_poscar.py",
  "$D/kriber/tst_strudat.py",
  "$D/cif/tests/tst_geometry.py",
  "$D/cif/tests/tst_crystal_symmetry_builder.py",
  "$D/cif/tests/tst_lex_parse_build.py",
  "$D/cif/tests/tst_model.py",
  "$D/cif/tests/tst_restraints.py",
  "$D/cif/tests/tst_validation.py",
  "$D/cif/tests/tst_ucif_examples_compilation.py",
  "$D/cif/tests/tst_parser.py",
  "$D/cif/tests/tst_citations.py",
  "$D/cif/tests/tst_model_builder.py",
  "$D/shelx/tst_lex_parse_build.py",
  "$D/shelx/tst_hklf.py",
  "$D/shelx/tst_writer.py",
  "$D/shelx/tst_fvar_encoding.py",
  "$D/pdb/tst_pdb.py",
  "$D/pdb/tst_mmcif.py",
  "$D/pdb/tst_mmcif_hierarchy.py",
  "$D/pdb/tst_mmcif_hierarchy_2.py",
  "$D/pdb/tst_tls.py",
  ["$D/pdb/hybrid_36.py", "exercise"],
  "$B/pdb/hybrid_36_fem",
  "$D/pdb/tst_hierarchy.py",
  "$D/pdb/tst_hierarchy_atom_sort.py",
  "$D/pdb/tst_hierarchy_flip_symmetric.py",
  "$D/regression/tst_selected_hierarchy_flip.py",
  "$D/pdb/tst_ext.py",
  "$D/pdb/tst_atom_selection.py",
  "$D/pdb/tst_rna_dna_atom_names.py",
  "$D/pdb/tst_atom_name_interpretation.py",
  "$D/pdb/tst_extract_rfactors_resolutions_sigma.py",
  "$D/pdb/modified_aa_names.py",
  "$D/pdb/modified_rna_dna_names.py",
  "$D/regression/secondary_structure/tst_sheet.py",
  "$D/regression/secondary_structure/tst_annotation.py",
  "$D/regression/secondary_structure/tst_annotation_long.py",
  "$D/pdb/secondary_structure.py",
  "$D/pdb/tst_atom_selection_string.py",
  "$D/pdb/tst_secondary_structure.py",
  "$D/pdb/tst_utils.py",
  "$D/pdb/tst_secondary_structure_2.py",
  "$D/pdb/remediation/tst_remediator.py",
  "$D/examples/pdb_to_map_simple.py",
  "$D/examples/pdb_truncate_to_ala/tst.py",
  "$D/examples/pdb_tardy_conf_sampling_simple.py",
  "$D/regression/tst_examples.py",
  "$D/cns/space_group_symbols.py",
  "$D/cns/tst_cns.py",
  ["$D/scalepack/tst_merge.py", "P31"],
  "$D/scalepack/no_merge_original_index.py",
  "$D/ccp4_map/tst.py",
  "$D/mtz/tst_ext.py",
  "$D/mtz/tst_extract_from_symmetry_lib.py",
  "$D/mtz/tst_dano.py",
  "$D/mtz/tst_miller_dict.py",
  "$D/mtz/tst_unmerged.py",
  ["$D/mtz/tst.py", "P31"],
  "$D/examples/tst_mtz_free_flipper.py",
  "$D/regression/tst_reflection_file_utils.py",
  "$D/detectors/tst_adsc.py",
  "$D/detectors/tst_detectors.py",
  "$D/xplor/tst_xplormap.py",
  ["$D/regression/tst_phases.py", "P31"],
  "$D/regression/tst_pdbx_mmcif_tutorial.py",
  "$D/regression/tst_lattice_symmetry.py",
  ["$D/regression/tst_reflection_statistics.py", "Fdd2 P31m"],
  "$D/regression/tst_data_plots.py",
  "$D/regression/tst_csv_utils.py",
  "$D/regression/tst_file_reader.py",
  "$D/regression/tst_bioinformatics.py",
  "$D/regression/tst_box_around_molecule.py",
  "$D/regression/tst_mmcif_segids.py",
  "$D/regression/tst_mmcif_input.py",
  "$D/regression/tst_hierarchy_forward_compatible_pdb.py",
  "$D/regression/tst_mmcif_multimodel.py",
  "$D/regression/tst_add_conformations.py",
  "$D/regression/tst_symmetry.py",
  "$D/regression/tst_reindex.py",
  "$D/regression/tst_reflection_file_editor.py",
  "$D/regression/tst_split_models.py",
  "$D/regression/tst_pdb_as_fasta.py",
  "$D/regression/tst_pdb_link_records.py",
  "$D/regression/tst_merging_statistics.py",
  "$D/regression/tst_simple_map_coefficients.py",
  "$D/regression/tst_sort_atoms.py",
  "$D/xds/tests/tst_xparm.py",
  "$D/xds/tests/tst_xds_inp.py",
  "$D/xds/tests/tst_integrate_hkl.py",
  "$D/xds/tests/tst_spots_xds.py",
  "$D/regression/tst_pdb_as_cif.py",
  "$D/scalepack/tst_no_merge_original_index.py",
  "$D/regression/tst_export_scalepack_unmerged.py",
  ["$D/dsn6/tst.py", "P31"],
  "$D/regression/ncs/tst_mtrix_biomt_cmdl.py",
  "$D/regression/ncs/tst_mmcif_biomt_reduction_output.py",
  "$D/regression/ncs/tst_ncs_search_ligs.py",
  "$D/regression/ncs/tst_ncs_search_broken_chain.py",
  "$D/regression/ncs/tst_ncs_search_shortcut_1.py",
  "$D/regression/ncs/tst_ncs_groups_preprocessing.py",
  "$D/regression/ncs/tst_ncs_input.py",
  "$D/regression/ncs/tst_ncs_user_selections.py",
  "$D/regression/ncs/tst_ncs.py",
  "$D/regression/ncs/tst_ncs_without_validation.py",
  "$D/pdb/tst_read_mtrix_records_from_cif.py",
  "$D/regression/tst_show_systematic_absences.py",
  "$D/regression/tst_miller_sort_asu.py",
  "$D/regression/tst_reflection_file_reader.py",
  "$D/regression/tst_xray_scale.py",
  "$D/bioinformatics/test/tst_alignment_as_hsearch.py",
  "$D/bioinformatics/test/tst_ebi_wu_blast_xml.py",
  "$D/bioinformatics/test/tst_ncbi_blast_xml.py",
  "$D/bioinformatics/pdb_info.py",
  "$D/regression/tst_cif_as_pdb_1atom.py",
  "$D/regression/tst_cif_1.py",
  "$D/regression/tst_split_data_cif.py",
  "$D/regression/tst_all_chain_ids.py",
  "$D/regression/tst_extract_xtal_data.py",
  "$D/regression/tst_cli_parser.py",
  "$D/regression/tst_mtz_as_cif.py",
  "$D/regression/tst_group_rounding.py",
  "$D/regression/tst_hierarchy_occupancies_rounding.py",
  "$D/regression/tst_hierarchy_merge_atoms_at_end_to_residues.py",
  "$D/regression/tst_hierarchy_long_chain_ids_1.py",
  "$D/regression/tst_hierarchy_long_resname_1.py",
  "$D/regression/tst_hierarchy_long_resname_2.py",
  "$D/regression/tst_hierarchy_long_resname_3.py",
  "$D/regression/tst_hierarchy_long_resname_4.py",
  "$D/regression/tst_hierarchy_copy_select.py",
  "$D/regression/tst_hierarchy_id_str.py",
  "$D/regression/tst_hierarchy_altlocs.py",
  "$D/regression/tst_fetch.py",
  ]

# failing tests on Windows, Python 2.7
tst_list_windows_fail = [
  "$D/detectors/tst_debug_write.py",
]

tst_list_fail = [
  "$D/regression/ncs/tst_ncs_reordered_chains.py",
  "$D/regression/tst_mmcif_to_from_hierarchy.py",
  ]
if sys.platform == 'win32':
  tst_list_fail += tst_list_windows_fail
else:
  tst_list_base += tst_list_windows_fail

tst_list_py3_unstable = []
tst_list_unstable = list()
if sys.version_info > (3, 0):
  tst_list_unstable += tst_list_py3_unstable
else:
  tst_list_base += tst_list_py3_unstable

# final lists
tst_list = tst_list_base
tst_list_expected_failures = tst_list_fail
tst_list_expected_unstable = tst_list_unstable

def run():
  build_dir = libtbx.env.under_build("iotbx")
  dist_dir = libtbx.env.dist_path("iotbx")

  test_utils.run_tests(build_dir, dist_dir, tst_list)

if (__name__ == "__main__"):
  run()
