from iotbx.kriber import strudat
from cctbx import crystal
import cctbx.crystal.coordination_sequences
import sys, os

def exercise(distance_cutoff=3.5, n_shells=5):
  verbose = "--Verbose" in sys.argv[1:]
  atlas_file = os.path.join(os.environ["LIBTBX_DIST_ROOT"],
    "regression", "misc", "strudat_zeolite_atlas")
  if (not os.path.isfile(atlas_file)): return
  all_entries = strudat.read_all_entries(open(atlas_file))
  for i,entry in enumerate(all_entries.entries):
    structure = entry.as_xray_structure()
    if ("--Full" in sys.argv[1:] or i % 20 == 0):
      asu_mappings = structure.asu_mappings(
        buffer_thickness=distance_cutoff)
      pair_asu_table = cctbx.crystal.pair_asu_table(
        asu_mappings=asu_mappings)
      pair_asu_table.add_all_pairs(
        distance_cutoff=distance_cutoff)
      term_table_slow = crystal.coordination_sequences.simple_and_slow(
        pair_asu_table=pair_asu_table,
        n_shells=n_shells)
      term_table_simple = cctbx.crystal.coordination_sequences_simple(
        asu_mappings=pair_asu_table.asu_mappings(),
        pair_asu_table_table=pair_asu_table.table(),
        n_shells=n_shells)
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

def run():
  exercise()
  print "OK"

if (__name__ == "__main__"):
  run()
