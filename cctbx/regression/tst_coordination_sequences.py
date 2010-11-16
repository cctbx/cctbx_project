from iotbx.kriber import strudat
from cctbx import crystal
import cctbx.crystal.coordination_sequences
from iotbx.option_parser import iotbx_option_parser
import libtbx.load_env
import sys, os

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
  term_table_simple = cctbx.crystal.coordination_sequences.simple(
    pair_asu_table=pair_asu_table,
    max_shell=max_shell)
  if (verbose):
    print "term_table_slow:"
    cctbx.crystal.coordination_sequences.show_terms(
      structure=structure,
      term_table=term_table_slow)
    print
    print "term_table_simple:"
    cctbx.crystal.coordination_sequences.show_terms(
      structure=structure,
      term_table=term_table_simple)
    print
  for terms_slow,terms_simple in zip(term_table_slow, term_table_simple):
    assert terms_slow == list(terms_simple)

def exercise_shell_asu_tables(structure, verbose):
  asu_mappings = structure.asu_mappings(buffer_thickness=10)
  bond_asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
  bond_asu_table.add_all_pairs(distance_cutoff=3.5)
  shell_asu_tables = crystal.coordination_sequences.shell_asu_tables(
    pair_asu_table=bond_asu_table,
    max_shell=3)
  for shell_asu_table in shell_asu_tables:
    if (0 or verbose):
      pairs_1 = structure.show_distances(pair_asu_table=shell_asu_table)
      print list(pairs_1.pair_counts)
      assert pairs_1.pair_counts == shell_asu_table.pair_counts()
      print
    sym_table = shell_asu_table.extract_pair_sym_table(
      skip_j_seq_less_than_i_seq=False)
    asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
    asu_table.add_pair_sym_table(sym_table=sym_table)
    if (0 or verbose):
      pairs_2 = structure.show_distances(pair_asu_table=asu_table)
      print list(pairs_2.pair_counts)
      assert pairs_2.pair_counts == asu_table.pair_counts()
      print
      assert pairs_1.pair_counts.all_eq(pairs_2.pair_counts)
    assert asu_table == shell_asu_table

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
  atlas_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/misc/strudat_zeolite_atlas",
    test=os.path.isfile)
  if (atlas_file is None):
    print "Skipping exercise(): input file not available"
    return
  all_entries = strudat.read_all_entries(open(atlas_file))
  for i,entry in enumerate(all_entries.entries):
    structure = entry.as_xray_structure()
    if (command_line.options.tag is not None):
      if (command_line.options.tag != entry.tag):
        continue
    elif (not (command_line.options.full or i % 20 == 0)):
      continue
    exercise_simple(
      structure, distance_cutoff, max_shell, command_line.options.verbose)
    exercise_shell_asu_tables(structure, command_line.options.verbose)

def run():
  exercise(sys.argv[1:])
  print "OK"

if (__name__ == "__main__"):
  run()
