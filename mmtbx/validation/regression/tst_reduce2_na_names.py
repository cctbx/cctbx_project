from __future__ import absolute_import, division, print_function
import time
import mmtbx.model
import iotbx.pdb
from mmtbx.hydrogens import reduce_hydrogen
from libtbx.utils import null_out

# ------------------------------------------------------------------------------
# Minimal regression test: reduce2 must emit PDB v3 (prime) atom names for the
# hydrogens it adds to RNA/DNA, not the old monomer-library (asterisk)
# convention.
#
# The monomer library names ribose/deoxyribose and 2'-OH hydrogens with
# asterisks (H4*, H3*, H2*, H2*1, H2*2, H1*, H5*1, H5*2, HO2*).  PDB v3
# requires prime notation (H4', H3', H2', H2'', H1', H5', H5'', HO2').
# Without the translation step in place_hydrogens(), reduce2 passes the
# old asterisk names straight through.  This test fails on such a build
# (e.g. plain master) and passes once the names are corrected.
#
# A single RNA residue (A) plus a single DNA residue (DA) is enough to
# exercise every distinct sugar/backbone hydrogen, including the RNA-only
# 2'-hydroxyl (HO2') and the geminal deoxyribose pair (H2'/H2'').  This is
# the small, self-contained companion to tst_reduce2_atom_names.py (which
# checks all eight standard bases against the remediation dictionary).
# ------------------------------------------------------------------------------

pdb_str = """\
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1
ATOM      1  P     A A   1      10.000  10.000  10.000  1.00 20.00           P
ATOM      2  OP1   A A   1      10.800  10.900  10.800  1.00 20.00           O
ATOM      3  OP2   A A   1       8.700  10.600  10.300  1.00 20.00           O
ATOM      4  O5'   A A   1      10.400   8.600   9.600  1.00 20.00           O
ATOM      5  C5'   A A   1      11.700   8.100   9.900  1.00 20.00           C
ATOM      6  C4'   A A   1      11.600   6.600   9.700  1.00 20.00           C
ATOM      7  O4'   A A   1      10.400   6.100  10.300  1.00 20.00           O
ATOM      8  C3'   A A   1      12.800   5.900  10.300  1.00 20.00           C
ATOM      9  O3'   A A   1      12.500   4.500  10.400  1.00 20.00           O
ATOM     10  C2'   A A   1      11.600   5.500  11.200  1.00 20.00           C
ATOM     11  O2'   A A   1      12.000   5.500  12.600  1.00 20.00           O
ATOM     12  C1'   A A   1      10.500   6.500  11.700  1.00 20.00           C
ATOM     13  N9    A A   1       9.100   6.000  11.800  1.00 20.00           N
ATOM     14  C8    A A   1       8.800   4.700  12.100  1.00 20.00           C
ATOM     15  N7    A A   1       7.500   4.500  12.100  1.00 20.00           N
ATOM     16  C5    A A   1       7.000   5.700  11.800  1.00 20.00           C
ATOM     17  C6    A A   1       5.600   6.100  11.600  1.00 20.00           C
ATOM     18  N6    A A   1       4.600   5.200  11.700  1.00 20.00           N
ATOM     19  N1    A A   1       5.400   7.400  11.300  1.00 20.00           N
ATOM     20  C2    A A   1       6.500   8.200  11.200  1.00 20.00           C
ATOM     21  N3    A A   1       7.800   7.800  11.400  1.00 20.00           N
ATOM     22  C4    A A   1       7.900   6.500  11.700  1.00 20.00           C
ATOM     23  P    DA B   1      20.000  10.000  10.000  1.00 20.00           P
ATOM     24  OP1  DA B   1      20.800  10.900  10.800  1.00 20.00           O
ATOM     25  OP2  DA B   1      18.700  10.600  10.300  1.00 20.00           O
ATOM     26  O5'  DA B   1      20.400   8.600   9.600  1.00 20.00           O
ATOM     27  C5'  DA B   1      21.700   8.100   9.900  1.00 20.00           C
ATOM     28  C4'  DA B   1      21.600   6.600   9.700  1.00 20.00           C
ATOM     29  O4'  DA B   1      20.400   6.100  10.300  1.00 20.00           O
ATOM     30  C3'  DA B   1      22.800   5.900  10.300  1.00 20.00           C
ATOM     31  O3'  DA B   1      22.500   4.500  10.400  1.00 20.00           O
ATOM     32  C2'  DA B   1      21.600   5.500  11.200  1.00 20.00           C
ATOM     33  C1'  DA B   1      20.500   6.500  11.700  1.00 20.00           C
ATOM     34  N9   DA B   1      19.100   6.000  11.800  1.00 20.00           N
ATOM     35  C8   DA B   1      18.800   4.700  12.100  1.00 20.00           C
ATOM     36  N7   DA B   1      17.500   4.500  12.100  1.00 20.00           N
ATOM     37  C5   DA B   1      17.000   5.700  11.800  1.00 20.00           C
ATOM     38  C6   DA B   1      15.600   6.100  11.600  1.00 20.00           C
ATOM     39  N6   DA B   1      14.600   5.200  11.700  1.00 20.00           N
ATOM     40  N1   DA B   1      15.400   7.400  11.300  1.00 20.00           N
ATOM     41  C2   DA B   1      16.500   8.200  11.200  1.00 20.00           C
ATOM     42  N3   DA B   1      17.800   7.800  11.400  1.00 20.00           N
ATOM     43  C4   DA B   1      17.900   6.500  11.700  1.00 20.00           C
END
"""

# Sugar/backbone hydrogens that reduce2 should add, per residue, in PDB v3
# (prime) convention.  The base hydrogens (H2, H8, H61, H62) are unaffected
# by the asterisk/prime issue and are not asserted here.
expected_rna_A = set(["H1'", "H2'", "H3'", "H4'", "H5'", "H5''", "HO2'"])
expected_dna_DA = set(["H1'", "H2'", "H2''", "H3'", "H4'", "H5'", "H5''"])

# ------------------------------------------------------------------------------

def get_added_h_names():
  '''
  Run reduce2 hydrogen placement on pdb_str (heavy atoms only) and return a
  dict mapping residue name -> set of newly added atom names.
  '''
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  model = mmtbx.model.manager(
    model_input       = pdb_inp,
    stop_for_unknowns = False,
    log               = null_out())
  model.add_crystal_symmetry_if_necessary()
  original = set(atom.name.strip() for atom in model.get_hierarchy().atoms())

  reduce_add_h_obj = reduce_hydrogen.place_hydrogens(
    model                = model,
    use_neutron_distances = False,
    n_terminal_charge    = "residue_one",
    exclude_water        = True,
    stop_for_unknowns    = False,
    keep_existing_H      = False)
  reduce_add_h_obj.run()
  model_h = reduce_add_h_obj.get_model()

  added = {}
  for atom in model_h.get_hierarchy().atoms():
    name = atom.name.strip()
    if name not in original:
      resname = atom.parent().resname.strip()
      added.setdefault(resname, set()).add(name)
  return added

# ------------------------------------------------------------------------------

def exercise():
  added = get_added_h_names()

  # something was actually placed
  all_added = set()
  for names in added.values():
    all_added.update(names)
  assert len(all_added) > 0, "reduce2 added no hydrogens"

  # the bug: no added hydrogen may carry an old-convention asterisk name
  asterisk = sorted(n for n in all_added if "*" in n)
  assert not asterisk, \
    "reduce2 added old-convention (asterisk) NA hydrogen names: %s" % asterisk

  # positive check: the expected PDB v3 (prime) names are present
  assert "A" in added, "RNA residue A was not hydrogenated"
  assert "DA" in added, "DNA residue DA was not hydrogenated"
  missing_rna = expected_rna_A - added["A"]
  assert not missing_rna, "missing RNA prime H names: %s" % sorted(missing_rna)
  missing_dna = expected_dna_DA - added["DA"]
  assert not missing_dna, "missing DNA prime H names: %s" % sorted(missing_dna)

# ------------------------------------------------------------------------------

def run():
  exercise()

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f"%(time.time()-t0))
