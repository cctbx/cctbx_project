from __future__ import absolute_import, division, print_function
import mmtbx.clashes
from libtbx.utils import null_out

pdb_good = """
CRYST1   23.341   28.568   19.164  90.00  90.00  90.00 P 1
ATOM      1  N   ARG A  58       9.158  17.337   8.990  1.00  7.73           N
ATOM      2  CA  ARG A  58      10.275  17.092   9.895  1.00  8.65           C
ATOM      3  C   ARG A  58       9.831  16.274  11.105  1.00  9.84           C
ATOM      4  O   ARG A  58       9.166  16.793  12.002  1.00  8.58           O
ATOM      5  CB  ARG A  58      10.895  18.414  10.352  1.00 20.00           C
ATOM      6  CG  ARG A  58      11.449  19.263   9.219  1.00 20.00           C
ATOM      7  CD  ARG A  58      12.052  20.557   9.743  1.00 20.00           C
ATOM      8  NE  ARG A  58      12.585  21.387   8.667  1.00 20.00           N
ATOM      9  CZ  ARG A  58      13.157  22.572   8.852  1.00 20.00           C
ATOM     10  NH1 ARG A  58      13.273  23.072  10.075  1.00 20.00           N
ATOM     11  NH2 ARG A  58      13.614  23.259   7.813  1.00 20.00           N
ATOM     12  N   GLU A  59      10.199  14.996  11.127  1.00  9.39           N
ATOM     13  CA  GLU A  59      10.987  14.395  10.056  1.00 11.89           C
ATOM     14  C   GLU A  59      10.157  13.393   9.261  1.00  9.81           C
ATOM     15  O   GLU A  59       9.832  12.313   9.753  1.00  8.97           O
ATOM     16  CB  GLU A  59      12.233  13.714  10.624  1.00 20.00           C
ATOM     17  CG  GLU A  59      13.173  14.654  11.361  1.00 20.00           C
ATOM     18  CD  GLU A  59      14.393  13.943  11.915  1.00 20.00           C
ATOM     19  OE1 GLU A  59      14.461  12.701  11.801  1.00 20.00           O
ATOM     20  OE2 GLU A  59      15.283  14.626  12.464  1.00 20.00           O
TER
"""

pdb_poor = """
CRYST1   23.341   28.568   19.164  90.00  90.00  90.00 P 1
ATOM      1  N   ARG A  58       9.158  17.337   8.990  1.00  7.73           N
ATOM      2  CA  ARG A  58      10.275  17.092   9.895  1.00  8.65           C
ATOM      3  C   ARG A  58       9.831  16.274  11.105  1.00  9.84           C
ATOM      4  O   ARG A  58       9.166  16.793  12.002  1.00  8.58           O
ATOM      5  CB  ARG A  58      10.895  18.414  10.352  1.00 20.00           C
ATOM      6  CG  ARG A  58      12.359  18.576   9.974  1.00 20.00           C
ATOM      7  CD  ARG A  58      13.136  17.290  10.213  1.00 20.00           C
ATOM      8  NE  ARG A  58      14.545  17.429   9.859  1.00 20.00           N
ATOM      9  CZ  ARG A  58      15.444  16.459   9.982  1.00 20.00           C
ATOM     10  NH1 ARG A  58      15.084  15.272  10.451  1.00 20.00           N
ATOM     11  NH2 ARG A  58      16.707  16.675   9.635  1.00 20.00           N
ATOM     12  N   GLU A  59      10.199  14.996  11.127  1.00  9.39           N
ATOM     13  CA  GLU A  59      10.987  14.395  10.056  1.00 11.89           C
ATOM     14  C   GLU A  59      10.157  13.393   9.261  1.00  9.81           C
ATOM     15  O   GLU A  59       9.832  12.313   9.753  1.00  8.97           O
ATOM     16  CB  GLU A  59      12.233  13.714  10.624  1.00 20.00           C
ATOM     17  CG  GLU A  59      13.155  14.647  11.392  1.00 20.00           C
ATOM     18  CD  GLU A  59      14.055  15.459  10.480  1.00 20.00           C
ATOM     19  OE1 GLU A  59      14.651  16.447  10.957  1.00 20.00           O
ATOM     20  OE2 GLU A  59      14.167  15.108   9.286  1.00 20.00           O
TER      21      GLU A  59
END
"""

def tst_01():
  o = mmtbx.clashes.from_pdb(pdb_str=pdb_good, clash_threshold=2.3)
  o.show()
  assert o.clashing_pairs() == [(9, 10), (18, 19), (3, 11)]
  print()
  o = mmtbx.clashes.from_pdb(pdb_str=pdb_poor, clash_threshold=1.5)
  o.show()
  expected_pairs = [(9, 18), (7, 19), (8, 19), (9, 19), (9, 17)]
  found_pairs = o.clashing_pairs()
  expected_pairs = sorted(expected_pairs)
  found_pairs = sorted(found_pairs)
  assert expected_pairs == found_pairs

def tst_02():
  # Test for remove_clashes
  import mmtbx.model
  from mmtbx.clashes import remove_clashes
  import iotbx.pdb
  import sys
  pdb_inp = iotbx.pdb.input(lines=pdb_poor.splitlines(),source_info='None')
  model= mmtbx.model.manager(model_input=pdb_inp)
  model.process(make_restraints=True)

  model.set_log(log = null_out())

  print("\n","-"*79)
  print(" Summary of input model statistics ")
  print("-"*79)
  model.get_restraints_manager()
  geometry = model.geometry_statistics()
  geometry.show(log = sys.stdout)


  rc=remove_clashes(model=model)

  print("\n","-"*79)
  print("Starting residues: %d " % (
     rc.model.get_hierarchy().overall_counts().n_residues))
  print("Side-chains removed: %d    Residues removed: %d" %(
     rc.side_chains_removed,
     rc.residues_removed))
  print("Final residues: %d " % (
     rc.new_model.get_hierarchy().overall_counts().n_residues))

  rc.new_model.set_log(log = null_out())
  rc.new_model.get_restraints_manager()
  new_geometry = rc.new_model.geometry_statistics()
  new_geometry.show(log = sys.stdout)
  assert rc.side_chains_removed==1
  assert rc.residues_removed==0


if (__name__ == "__main__"):
  tst_01()
  tst_02()
