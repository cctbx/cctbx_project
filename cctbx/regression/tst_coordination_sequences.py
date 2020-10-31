from __future__ import absolute_import, division, print_function
from iotbx.kriber import strudat
from cctbx import crystal
import cctbx.crystal.coordination_sequences
from iotbx.option_parser import iotbx_option_parser
import libtbx.load_env
import sys, os
from six.moves import zip

def exercise_simple(structure, distance_cutoff, max_shell, verbose):
  asu_mappings = structure.asu_mappings(
    buffer_thickness=distance_cutoff)
  pair_asu_table = cctbx.crystal.pair_asu_table(
    asu_mappings=asu_mappings)
  pair_asu_table.add_all_pairs(
    distance_cutoff=distance_cutoff)
  term_table_slow = crystal.coordination_sequences.simple_and_slow(
    pair_asu_table=pair_asu_table,
    max_shell=max_shell)
  term_table_simple_asu = cctbx.crystal.coordination_sequences.simple(
    pair_asu_table=pair_asu_table,
    max_shell=max_shell)
  site_symmetry_table = structure.site_symmetry_table()
  full_pair_sym_table = pair_asu_table.extract_pair_sym_table() \
    .full_connectivity(site_symmetry_table=site_symmetry_table)
  term_table_simple_sym = cctbx.crystal.coordination_sequences.simple_sym(
    full_pair_sym_table=full_pair_sym_table,
    site_symmetry_table=site_symmetry_table,
    max_shell=max_shell)
  if (verbose):
    print("term_table_slow:")
    cctbx.crystal.coordination_sequences.show_terms(
      structure=structure,
      term_table=term_table_slow)
    print()
    print("term_table_simple_asu:")
    cctbx.crystal.coordination_sequences.show_terms(
      structure=structure,
      term_table=term_table_simple_asu)
    print()
    print("term_table_simple_sym:")
    cctbx.crystal.coordination_sequences.show_terms(
      structure=structure,
      term_table=term_table_simple_sym)
    print()
  for term_table_simple in [term_table_simple_asu, term_table_simple_sym]:
    for terms_slow,terms_simple in zip(term_table_slow, term_table_simple):
      assert terms_slow == list(terms_simple)

def exercise_shell_asu_tables(structure, verbose):
  asu_mappings = structure.asu_mappings(buffer_thickness=10)
  bond_asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
  bond_asu_table.add_all_pairs(distance_cutoff=3.5)
  shell_asu_tables = crystal.coordination_sequences.shell_asu_tables(
    pair_asu_table=bond_asu_table,
    max_shell=3)
  site_symmetry_table = structure.site_symmetry_table()
  full_pair_sym_table = bond_asu_table.extract_pair_sym_table() \
    .full_connectivity(site_symmetry_table=site_symmetry_table)
  shell_sym_tables_orig = crystal.coordination_sequences.shell_sym_tables(
    full_pair_sym_table=full_pair_sym_table,
    site_symmetry_table=site_symmetry_table,
    max_shell=3)
  shell_sym_tables = [_.tidy(site_symmetry_table=site_symmetry_table)
    for _ in shell_sym_tables_orig]
  have_redundancies = False
  for o_pst,t_pst in zip(shell_sym_tables_orig[1:], shell_sym_tables[1:]):
    for o_pair_sym_dict,t_pair_sym_dict in zip(o_pst, t_pst):
      assert list(o_pair_sym_dict.keys()) == list(t_pair_sym_dict.keys())
      if (verbose and not have_redundancies):
        for j_seq,o_sym_ops in o_pair_sym_dict.items():
          t_sym_ops = t_pair_sym_dict[j_seq]
          if (len(t_sym_ops) != len(o_sym_ops)):
            have_redundancies = True
  if (have_redundancies):
    print("crystal.coordination_sequences.shell_sym_tables redundancies:")
    print("original:")
    o_pst.show()
    print()
    print("tidy:")
    t_pst.show()
    print()
  for shell_asu_table,shell_sym_table in zip(
        shell_asu_tables, shell_sym_tables):
    if (0 or verbose):
      pairs_1 = structure.show_distances(pair_asu_table=shell_asu_table) \
        .distances_info
      print(list(pairs_1.pair_counts))
      assert pairs_1.pair_counts == shell_asu_table.pair_counts()
      print()
    sym_table = shell_asu_table.extract_pair_sym_table(
      skip_j_seq_less_than_i_seq=False)
    asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
    asu_table.add_pair_sym_table(sym_table=sym_table)
    if (0 or verbose):
      pairs_2 = structure.show_distances(pair_asu_table=asu_table) \
        .distances_info
      print(list(pairs_2.pair_counts))
      assert pairs_2.pair_counts == asu_table.pair_counts()
      print()
      assert pairs_1.pair_counts.all_eq(pairs_2.pair_counts)
    assert asu_table == shell_asu_table
    shell_sym_from_asu_table = shell_asu_table.extract_pair_sym_table() \
      .tidy(site_symmetry_table=site_symmetry_table)
    from six.moves import cStringIO as StringIO
    sio_sym = StringIO()
    shell_sym_table.show(f=sio_sym)
    sio_asu = StringIO()
    shell_sym_from_asu_table.show(f=sio_asu)
    from libtbx.test_utils import show_diff
    assert not show_diff(sio_sym.getvalue(), sio_asu.getvalue())

def exercise(args, distance_cutoff=3.5, max_shell=5):
  command_line = (iotbx_option_parser()
    .option(None, "--tag",
      action="store",
      type="string")
    .option(None, "--full",
      action="store_true")
    .option(None, "--verbose",
      action="store_true")
  ).process(args=args)
  co = command_line.options
  atlas_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/misc/strudat_zeolite_atlas",
    test=os.path.isfile)
  if (atlas_file is None):
    print("Skipping exercise(): input file not available")
    return
  with open(atlas_file) as f:
    all_entries = strudat.read_all_entries(f)
  for i,entry in enumerate(all_entries.entries):
    structure = entry.as_xray_structure()
    if (co.tag is not None):
      if (co.tag != entry.tag):
        continue
    elif (not (co.full or i % 20 == 0)):
      continue
    if (co.verbose):
      print("strudat tag:", entry.tag)
    exercise_simple(
      structure, distance_cutoff, max_shell, co.verbose)
    exercise_shell_asu_tables(structure, co.verbose)

def run():
  exercise(sys.argv[1:])
  print("OK")

if (__name__ == "__main__"):
  run()
