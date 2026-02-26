"""Test that reduce2 adds only new-convention PDB atom names to nucleic acids.

Constructs a PDB with one each of the RNA bases (A, G, C, U) and DNA bases
(DA, DG, DC, DT), runs reduce2 hydrogen placement, and verifies that all
newly-added atom names use the modern PDB convention (primes instead of
asterisks, proper digit placement, etc).

The old-name set is built from iotbx/pdb/remediation/remediation.dict.
"""

from __future__ import absolute_import, division, print_function
import os
import sys
import mmtbx.model
import mmtbx.hydrogens.reduce_hydrogen as reduce_hydrogen
import iotbx.pdb
import iotbx.phil
import libtbx.load_env
from libtbx.utils import null_out

# Minimal nucleic acid PDB: RNA strand (A, G, C, U) + DNA strand (DA, DG, DC, DT).
# Heavy atoms only, no hydrogens.  Coordinates are approximate but chemically
# reasonable so that reduce2 can place hydrogens.
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
ATOM     23  P     G A   2      13.100   3.600   9.200  1.00 20.00           P
ATOM     24  OP1   G A   2      14.500   3.900   8.900  1.00 20.00           O
ATOM     25  OP2   G A   2      12.200   3.900   8.100  1.00 20.00           O
ATOM     26  O5'   G A   2      12.800   2.100   9.500  1.00 20.00           O
ATOM     27  C5'   G A   2      13.500   1.100  10.200  1.00 20.00           C
ATOM     28  C4'   G A   2      12.600  -0.100  10.300  1.00 20.00           C
ATOM     29  O4'   G A   2      11.300   0.300  10.800  1.00 20.00           O
ATOM     30  C3'   G A   2      13.100  -1.200  11.200  1.00 20.00           C
ATOM     31  O3'   G A   2      13.000  -2.500  10.600  1.00 20.00           O
ATOM     32  C2'   G A   2      12.100  -1.100  12.300  1.00 20.00           C
ATOM     33  O2'   G A   2      12.700  -1.300  13.600  1.00 20.00           O
ATOM     34  C1'   G A   2      11.100   0.000  11.900  1.00 20.00           C
ATOM     35  N9    G A   2       9.700  -0.400  12.100  1.00 20.00           N
ATOM     36  C8    G A   2       9.200  -1.700  12.000  1.00 20.00           C
ATOM     37  N7    G A   2       7.900  -1.700  12.200  1.00 20.00           N
ATOM     38  C5    G A   2       7.600  -0.400  12.400  1.00 20.00           C
ATOM     39  C6    G A   2       6.300   0.200  12.700  1.00 20.00           C
ATOM     40  O6    G A   2       5.200  -0.400  12.800  1.00 20.00           O
ATOM     41  N1    G A   2       6.400   1.600  12.800  1.00 20.00           N
ATOM     42  C2    G A   2       7.600   2.300  12.700  1.00 20.00           C
ATOM     43  N2    G A   2       7.500   3.600  12.800  1.00 20.00           N
ATOM     44  N3    G A   2       8.800   1.800  12.400  1.00 20.00           N
ATOM     45  C4    G A   2       8.700   0.400  12.300  1.00 20.00           C
ATOM     46  P     C A   3      14.100  -3.600  11.100  1.00 20.00           P
ATOM     47  OP1   C A   3      15.300  -3.100  11.700  1.00 20.00           O
ATOM     48  OP2   C A   3      13.400  -4.800  11.600  1.00 20.00           O
ATOM     49  O5'   C A   3      14.400  -3.800   9.500  1.00 20.00           O
ATOM     50  C5'   C A   3      15.200  -4.900   9.100  1.00 20.00           C
ATOM     51  C4'   C A   3      14.400  -6.100   8.700  1.00 20.00           C
ATOM     52  O4'   C A   3      13.200  -5.700   8.000  1.00 20.00           O
ATOM     53  C3'   C A   3      15.100  -7.100   7.800  1.00 20.00           C
ATOM     54  O3'   C A   3      15.200  -8.400   8.400  1.00 20.00           O
ATOM     55  C2'   C A   3      14.100  -7.100   6.700  1.00 20.00           C
ATOM     56  O2'   C A   3      14.700  -7.500   5.500  1.00 20.00           O
ATOM     57  C1'   C A   3      13.100  -6.000   6.600  1.00 20.00           C
ATOM     58  N1    C A   3      11.700  -6.400   6.400  1.00 20.00           N
ATOM     59  C2    C A   3      11.300  -7.100   5.200  1.00 20.00           C
ATOM     60  O2    C A   3      12.100  -7.400   4.300  1.00 20.00           O
ATOM     61  N3    C A   3       9.900  -7.400   5.100  1.00 20.00           N
ATOM     62  C4    C A   3       9.100  -7.000   6.100  1.00 20.00           C
ATOM     63  N4    C A   3       7.800  -7.300   6.000  1.00 20.00           N
ATOM     64  C5    C A   3       9.500  -6.200   7.200  1.00 20.00           C
ATOM     65  C6    C A   3      10.800  -6.000   7.300  1.00 20.00           C
ATOM     66  P     U A   4      16.200  -9.500   7.700  1.00 20.00           P
ATOM     67  OP1   U A   4      17.500  -9.000   8.200  1.00 20.00           O
ATOM     68  OP2   U A   4      15.600 -10.800   8.000  1.00 20.00           O
ATOM     69  O5'   U A   4      16.300  -9.300   6.100  1.00 20.00           O
ATOM     70  C5'   U A   4      17.300 -10.000   5.400  1.00 20.00           C
ATOM     71  C4'   U A   4      16.700 -11.100   4.500  1.00 20.00           C
ATOM     72  O4'   U A   4      15.500 -10.600   3.800  1.00 20.00           O
ATOM     73  C3'   U A   4      17.500 -11.800   3.400  1.00 20.00           C
ATOM     74  O3'   U A   4      17.600 -13.200   3.600  1.00 20.00           O
ATOM     75  C2'   U A   4      16.600 -11.500   2.200  1.00 20.00           C
ATOM     76  O2'   U A   4      17.300 -11.600   1.000  1.00 20.00           O
ATOM     77  C1'   U A   4      15.500 -10.500   2.500  1.00 20.00           C
ATOM     78  N1    U A   4      14.200 -11.100   2.200  1.00 20.00           N
ATOM     79  C2    U A   4      14.000 -11.800   1.000  1.00 20.00           C
ATOM     80  O2    U A   4      14.900 -12.000   0.200  1.00 20.00           O
ATOM     81  N3    U A   4      12.700 -12.200   0.800  1.00 20.00           N
ATOM     82  C4    U A   4      11.600 -12.000   1.600  1.00 20.00           C
ATOM     83  O4    U A   4      10.500 -12.400   1.200  1.00 20.00           O
ATOM     84  C5    U A   4      11.900 -11.200   2.800  1.00 20.00           C
ATOM     85  C6    U A   4      13.100 -10.800   3.000  1.00 20.00           C
ATOM     86  P    DA B   1      20.000  10.000  10.000  1.00 20.00           P
ATOM     87  OP1  DA B   1      20.800  10.900  10.800  1.00 20.00           O
ATOM     88  OP2  DA B   1      18.700  10.600  10.300  1.00 20.00           O
ATOM     89  O5'  DA B   1      20.400   8.600   9.600  1.00 20.00           O
ATOM     90  C5'  DA B   1      21.700   8.100   9.900  1.00 20.00           C
ATOM     91  C4'  DA B   1      21.600   6.600   9.700  1.00 20.00           C
ATOM     92  O4'  DA B   1      20.400   6.100  10.300  1.00 20.00           O
ATOM     93  C3'  DA B   1      22.800   5.900  10.300  1.00 20.00           C
ATOM     94  O3'  DA B   1      22.500   4.500  10.400  1.00 20.00           O
ATOM     95  C2'  DA B   1      21.600   5.500  11.200  1.00 20.00           C
ATOM     96  C1'  DA B   1      20.500   6.500  11.700  1.00 20.00           C
ATOM     97  N9   DA B   1      19.100   6.000  11.800  1.00 20.00           N
ATOM     98  C8   DA B   1      18.800   4.700  12.100  1.00 20.00           C
ATOM     99  N7   DA B   1      17.500   4.500  12.100  1.00 20.00           N
ATOM    100  C5   DA B   1      17.000   5.700  11.800  1.00 20.00           C
ATOM    101  C6   DA B   1      15.600   6.100  11.600  1.00 20.00           C
ATOM    102  N6   DA B   1      14.600   5.200  11.700  1.00 20.00           N
ATOM    103  N1   DA B   1      15.400   7.400  11.300  1.00 20.00           N
ATOM    104  C2   DA B   1      16.500   8.200  11.200  1.00 20.00           C
ATOM    105  N3   DA B   1      17.800   7.800  11.400  1.00 20.00           N
ATOM    106  C4   DA B   1      17.900   6.500  11.700  1.00 20.00           C
ATOM    107  P    DG B   2      23.100   3.600   9.200  1.00 20.00           P
ATOM    108  OP1  DG B   2      24.500   3.900   8.900  1.00 20.00           O
ATOM    109  OP2  DG B   2      22.200   3.900   8.100  1.00 20.00           O
ATOM    110  O5'  DG B   2      22.800   2.100   9.500  1.00 20.00           O
ATOM    111  C5'  DG B   2      23.500   1.100  10.200  1.00 20.00           C
ATOM    112  C4'  DG B   2      22.600  -0.100  10.300  1.00 20.00           C
ATOM    113  O4'  DG B   2      21.300   0.300  10.800  1.00 20.00           O
ATOM    114  C3'  DG B   2      23.100  -1.200  11.200  1.00 20.00           C
ATOM    115  O3'  DG B   2      23.000  -2.500  10.600  1.00 20.00           O
ATOM    116  C2'  DG B   2      22.100  -1.100  12.300  1.00 20.00           C
ATOM    117  C1'  DG B   2      21.100   0.000  11.900  1.00 20.00           C
ATOM    118  N9   DG B   2      19.700  -0.400  12.100  1.00 20.00           N
ATOM    119  C8   DG B   2      19.200  -1.700  12.000  1.00 20.00           C
ATOM    120  N7   DG B   2      17.900  -1.700  12.200  1.00 20.00           N
ATOM    121  C5   DG B   2      17.600  -0.400  12.400  1.00 20.00           C
ATOM    122  C6   DG B   2      16.300   0.200  12.700  1.00 20.00           C
ATOM    123  O6   DG B   2      15.200  -0.400  12.800  1.00 20.00           O
ATOM    124  N1   DG B   2      16.400   1.600  12.800  1.00 20.00           N
ATOM    125  C2   DG B   2      17.600   2.300  12.700  1.00 20.00           C
ATOM    126  N2   DG B   2      17.500   3.600  12.800  1.00 20.00           N
ATOM    127  N3   DG B   2      18.800   1.800  12.400  1.00 20.00           N
ATOM    128  C4   DG B   2      18.700   0.400  12.300  1.00 20.00           C
ATOM    129  P    DC B   3      24.100  -3.600  11.100  1.00 20.00           P
ATOM    130  OP1  DC B   3      25.300  -3.100  11.700  1.00 20.00           O
ATOM    131  OP2  DC B   3      23.400  -4.800  11.600  1.00 20.00           O
ATOM    132  O5'  DC B   3      24.400  -3.800   9.500  1.00 20.00           O
ATOM    133  C5'  DC B   3      25.200  -4.900   9.100  1.00 20.00           C
ATOM    134  C4'  DC B   3      24.400  -6.100   8.700  1.00 20.00           C
ATOM    135  O4'  DC B   3      23.200  -5.700   8.000  1.00 20.00           O
ATOM    136  C3'  DC B   3      25.100  -7.100   7.800  1.00 20.00           C
ATOM    137  O3'  DC B   3      25.200  -8.400   8.400  1.00 20.00           O
ATOM    138  C2'  DC B   3      24.100  -7.100   6.700  1.00 20.00           C
ATOM    139  C1'  DC B   3      23.100  -6.000   6.600  1.00 20.00           C
ATOM    140  N1   DC B   3      21.700  -6.400   6.400  1.00 20.00           N
ATOM    141  C2   DC B   3      21.300  -7.100   5.200  1.00 20.00           C
ATOM    142  O2   DC B   3      22.100  -7.400   4.300  1.00 20.00           O
ATOM    143  N3   DC B   3      19.900  -7.400   5.100  1.00 20.00           N
ATOM    144  C4   DC B   3      19.100  -7.000   6.100  1.00 20.00           C
ATOM    145  N4   DC B   3      17.800  -7.300   6.000  1.00 20.00           N
ATOM    146  C5   DC B   3      19.500  -6.200   7.200  1.00 20.00           C
ATOM    147  C6   DC B   3      20.800  -6.000   7.300  1.00 20.00           C
ATOM    148  P    DT B   4      26.200  -9.500   7.700  1.00 20.00           P
ATOM    149  OP1  DT B   4      27.500  -9.000   8.200  1.00 20.00           O
ATOM    150  OP2  DT B   4      25.600 -10.800   8.000  1.00 20.00           O
ATOM    151  O5'  DT B   4      26.300  -9.300   6.100  1.00 20.00           O
ATOM    152  C5'  DT B   4      27.300 -10.000   5.400  1.00 20.00           C
ATOM    153  C4'  DT B   4      26.700 -11.100   4.500  1.00 20.00           C
ATOM    154  O4'  DT B   4      25.500 -10.600   3.800  1.00 20.00           O
ATOM    155  C3'  DT B   4      27.500 -11.800   3.400  1.00 20.00           C
ATOM    156  O3'  DT B   4      27.600 -13.200   3.600  1.00 20.00           O
ATOM    157  C2'  DT B   4      26.600 -11.500   2.200  1.00 20.00           C
ATOM    158  C1'  DT B   4      25.500 -10.500   2.500  1.00 20.00           C
ATOM    159  N1   DT B   4      24.200 -11.100   2.200  1.00 20.00           N
ATOM    160  C2   DT B   4      24.000 -11.800   1.000  1.00 20.00           C
ATOM    161  O2   DT B   4      24.900 -12.000   0.200  1.00 20.00           O
ATOM    162  N3   DT B   4      22.700 -12.200   0.800  1.00 20.00           N
ATOM    163  C4   DT B   4      21.600 -12.000   1.600  1.00 20.00           C
ATOM    164  O4   DT B   4      20.500 -12.400   1.200  1.00 20.00           O
ATOM    165  C5   DT B   4      21.900 -11.200   2.800  1.00 20.00           C
ATOM    166  C7   DT B   4      20.800 -10.900   3.800  1.00 20.00           C
ATOM    167  C6   DT B   4      23.100 -10.800   3.000  1.00 20.00           C
END
"""

def parse_remediation_dict():
  """Parse remediation.dict to get old atom names per nucleic acid residue."""
  dict_path = libtbx.env.find_in_repositories(
    relative_path="iotbx/pdb/remediation/remediation.dict",
    test=os.path.isfile)
  assert dict_path is not None, "Cannot find remediation.dict"

  na_resnames = {'A', 'G', 'C', 'U', 'DA', 'DG', 'DC', 'DT'}
  # old_name -> set of resnames it applies to
  old_names = {}
  with open(dict_path) as f:
    for line in f:
      line = line.strip()
      if not line or line.startswith('#'):
        continue
      parts = line.split(':')
      if len(parts) != 2:
        continue
      new_side = parts[0].strip().split()
      old_side = parts[1].strip().split()
      if len(new_side) >= 2 and len(old_side) >= 2:
        new_name, new_res = new_side[0], new_side[1]
        old_name, old_res = old_side[0], old_side[1]
        if new_res in na_resnames:
          key = (old_name, old_res)
          old_names[key] = new_name
  return old_names

def exercise():
  """Run reduce2 on nucleic acid PDB and check for old-style atom names."""
  # Build old-name lookup: (atom_name, resname) -> new_name
  old_name_map = parse_remediation_dict()
  # Also build a set of all old atom names regardless of residue, for a
  # broader check
  all_old_names = set()
  for (old_name, old_res), new_name in old_name_map.items():
    all_old_names.add(old_name)

  # Build model
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  model = mmtbx.model.manager(
    model_input=pdb_inp,
    stop_for_unknowns=False,
    log=null_out())
  model.add_crystal_symmetry_if_necessary()

  # Record original atom names
  original_names = set()
  for atom in model.get_hierarchy().atoms():
    original_names.add(atom.name.strip())

  # Run reduce2 hydrogen placement
  reduce_obj = reduce_hydrogen.place_hydrogens(
    model=model,
    use_neutron_distances=False,
    n_terminal_charge="residue_one",
    exclude_water=True,
    stop_for_unknowns=False,
    keep_existing_H=False,
  )
  reduce_obj.run()
  h_model = reduce_obj.get_model()

  # Collect newly-added atom names per residue
  new_atoms = []
  for atom in h_model.get_hierarchy().atoms():
    name = atom.name.strip()
    if name not in original_names:
      rg = atom.parent()
      resname = rg.resname.strip() if rg else "UNK"
      new_atoms.append((name, resname))

  assert len(new_atoms) > 0, "No new atoms were added by reduce2"

  # Check each newly-added atom name against old-name mapping
  bad_atoms = []
  for atom_name, resname in new_atoms:
    # Check residue-specific old name
    if (atom_name, resname) in old_name_map:
      new_name = old_name_map[(atom_name, resname)]
      bad_atoms.append(
        "%s %s (should be %s)" % (atom_name, resname, new_name))
    # Also check for the telltale patterns of old naming:
    # asterisk instead of prime
    elif '*' in atom_name:
      bad_atoms.append("%s %s (contains asterisk)" % (atom_name, resname))

  if bad_atoms:
    print("FAIL: reduce2 added old-convention atom names:", file=sys.stderr)
    for b in bad_atoms:
      print("  %s" % b, file=sys.stderr)
    raise AssertionError(
      "reduce2 added %d old-convention atom name(s):\n  %s"
      % (len(bad_atoms), "\n  ".join(bad_atoms)))

  # Verify we got hydrogens with new-convention names
  h_count = 0
  prime_count = 0
  for atom_name, resname in new_atoms:
    if atom_name.startswith('H') or atom_name.startswith('D'):
      h_count += 1
    if "'" in atom_name:
      prime_count += 1

  print("OK: reduce2 added %d new atoms (%d hydrogens, %d with primes)"
        % (len(new_atoms), h_count, prime_count))
  assert h_count > 0, "No hydrogen atoms were added"
  assert prime_count > 0, "No prime-notation atoms found in added hydrogens"

if __name__ == "__main__":
  if not libtbx.env.has_module(name="reduce"):
    print("Skipping: reduce module not available")
  else:
    exercise()
    print("OK")
