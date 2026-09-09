from __future__ import absolute_import, division, print_function
import math

from wxtbx.regression._warnings import install_wx_deprecation_filters
install_wx_deprecation_filters()

# Two residues, a numbering gap, then an alt-conf residue and a
# partial-occupancy residue: exercises every branch of analyze.make_plots.
pdb_str = """\
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1
ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00 10.00           N
ATOM      2  CA  GLY A   1       1.500   0.000   0.000  1.00 10.00           C
ATOM      3  C   GLY A   1       2.000   1.400   0.000  1.00 10.00           C
ATOM      4  O   GLY A   1       1.300   2.400   0.000  1.00 10.00           O
ATOM      5  N   GLY A   2       3.300   1.500   0.000  1.00 20.00           N
ATOM      6  CA  GLY A   2       4.000   2.800   0.000  1.00 20.00           C
ATOM      7  C   GLY A   2       5.500   2.800   0.000  1.00 20.00           C
ATOM      8  O   GLY A   2       6.200   3.800   0.000  1.00 20.00           O
ATOM      9  N  AALA A   5       7.000   0.000   0.000  0.50 30.00           N
ATOM     10  CA AALA A   5       8.500   0.000   0.000  0.50 30.00           C
ATOM     11  C  AALA A   5       9.000   1.400   0.000  0.50 30.00           C
ATOM     12  O  AALA A   5       8.300   2.400   0.000  0.50 30.00           O
ATOM     13  CB AALA A   5       9.000  -1.000   1.000  0.50 30.00           C
ATOM     14  N  BALA A   5       7.000   0.100   0.000  0.50 32.00           N
ATOM     15  CA BALA A   5       8.500   0.100   0.000  0.50 32.00           C
ATOM     16  C  BALA A   5       9.000   1.500   0.000  0.50 32.00           C
ATOM     17  O  BALA A   5       8.300   2.500   0.000  0.50 32.00           O
ATOM     18  CB BALA A   5       9.000  -0.900   1.000  0.50 32.00           C
ATOM     19  N   GLY A   6      10.300   1.500   0.000  0.80 40.00           N
ATOM     20  CA  GLY A   6      11.000   2.800   0.000  0.80 40.00           C
ATOM     21  C   GLY A   6      12.500   2.800   0.000  0.80 40.00           C
ATOM     22  O   GLY A   6      13.200   3.800   0.000  0.80 40.00           O
END
"""

def exercise():
  import iotbx.pdb
  from libtbx.utils import null_out
  from wxtbx import b_plot
  pdb_in = iotbx.pdb.input(source_info=None, lines=pdb_str)
  hierarchy = pdb_in.construct_hierarchy()
  hierarchy.atoms().reset_i_seq()
  xrs = pdb_in.xray_structure_simple()
  params = b_plot.master_phil.extract()
  result = b_plot.analyze(pdb_hierarchy=hierarchy, xray_structure=xrs,
    params=params.b_plot, out=null_out())
  plots = result.make_plots()
  assert len(plots) == 1
  label, b_vals, is_altconf, is_partocc, labels = plots[0]
  assert label == "Chain 'A' (1 - 6)", label
  # residues 1, 2, [gap of 3: 3, 4, and the re-counted 5], 5, 6
  assert labels == ["1", "2", None, None, None, "5", "6"], labels
  assert len(b_vals) == len(is_altconf) == len(is_partocc) == 7
  # gap entries are NaN, real residues carry their mean B
  assert [math.isnan(v) for v in b_vals] == [False, False, True, True, True,
                                             False, False]
  assert abs(b_vals[0] - 10) < 1e-6 and abs(b_vals[1] - 20) < 1e-6
  assert abs(b_vals[5] - 31) < 1e-6 and abs(b_vals[6] - 40) < 1e-6
  # only residue 5 has an alt conf; only residue 6 has a partial occupancy
  assert [math.isnan(v) for v in is_altconf] == [True, True, True, True, True,
                                                 False, True]
  assert [math.isnan(v) for v in is_partocc] == [True, True, True, True, True,
                                                 True, False]
  assert abs(is_altconf[5] - 8) < 1e-6   # max(min(avg_b) - 2, 0)
  assert abs(is_partocc[6] - 8) < 1e-6

if __name__ == "__main__":
  exercise()
  print("OK")
