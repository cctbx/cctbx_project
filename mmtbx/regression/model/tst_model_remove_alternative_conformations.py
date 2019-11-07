from __future__ import absolute_import, division, print_function
import mmtbx.model
import iotbx.pdb
from libtbx.test_utils import show_diff

tst_pdb_str = """\
ATOM  18920  N   PHE L   2     201.672 235.270 272.436  1.00 33.60           N
ATOM  18921  CA  PHE L   2     201.753 236.596 271.840  1.00 33.60           C
ATOM  18922  C   PHE L   2     200.597 237.520 272.192  1.00 33.60           C
ATOM  18923  O   PHE L   2     200.835 238.721 272.367  1.00 33.60           O
ATOM  18924  CB  PHE L   2     201.883 236.467 270.321  1.00 32.90           C
ATOM  18925  CG  PHE L   2     201.850 237.775 269.616  1.00 32.90           C
ATOM  18926  CD1 PHE L   2     202.865 238.689 269.797  1.00 32.90           C
ATOM  18927  CD2 PHE L   2     200.820 238.083 268.748  1.00 32.90           C
ATOM  18928  CE1 PHE L   2     202.829 239.908 269.166  1.00 32.90           C
ATOM  18929  CE2 PHE L   2     200.799 239.290 268.089  1.00 32.90           C
ATOM  18930  CZ  PHE L   2     201.803 240.201 268.300  1.00 32.90           C
ATOM  18931  N   SER L   3     199.357 236.979 272.328  1.00 34.28           N
ATOM  18932  C   SER L   3     197.907 238.857 271.631  1.00 34.28           C
ATOM  18933  O   SER L   3     198.129 240.024 271.974  1.00 34.28           O
ATOM  18934  CA ASER L   3     198.124 237.714 272.618  0.50 34.28           C
ATOM  18935  CB ASER L   3     198.121 238.241 274.053  0.50 35.83           C
ATOM  18936  OG ASER L   3     198.875 239.432 274.158  0.50 35.83           O
ATOM  18937  CA BSER L   3     198.130 237.722 272.621  0.50 34.28           C
ATOM  18938  CB BSER L   3     198.142 238.287 274.041  0.50 35.83           C
ATOM  18939  OG BSER L   3     197.037 239.151 274.236  0.50 35.83           O
ATOM  18940  N   PRO L   4     197.482 238.564 270.395  1.00 33.41           N
ATOM  18941  CA  PRO L   4     197.440 239.596 269.355  1.00 33.41           C
ATOM  18942  C   PRO L   4     196.353 240.646 269.519  1.00 33.41           C
ATOM  18943  O   PRO L   4     196.212 241.494 268.633  1.00 33.41           O
ATOM  18944  CB  PRO L   4     197.209 238.776 268.079  1.00 33.14           C
ATOM  18945  CG  PRO L   4     196.474 237.597 268.538  1.00 33.14           C
ATOM  18946  CD  PRO L   4     197.000 237.268 269.901  1.00 33.14           C
"""

def exercise_1(prefix="tst_model_remove_alternative_conformations_1"):
  """ Make sure that CA in SER3 gets to the right position. """
  inp = iotbx.pdb.input(lines=tst_pdb_str, source_info=None)
  model = mmtbx.model.manager(
      model_input = inp)
  model.remove_alternative_conformations(always_keep_one_conformer=False)
  assert not show_diff(model.model_as_pdb(), """\
ATOM      1  N   PHE L   2     201.672 235.270 272.436  1.00 33.60           N
ATOM      2  CA  PHE L   2     201.753 236.596 271.840  1.00 33.60           C
ATOM      3  C   PHE L   2     200.597 237.520 272.192  1.00 33.60           C
ATOM      4  O   PHE L   2     200.835 238.721 272.367  1.00 33.60           O
ATOM      5  CB  PHE L   2     201.883 236.467 270.321  1.00 32.90           C
ATOM      6  CG  PHE L   2     201.850 237.775 269.616  1.00 32.90           C
ATOM      7  CD1 PHE L   2     202.865 238.689 269.797  1.00 32.90           C
ATOM      8  CD2 PHE L   2     200.820 238.083 268.748  1.00 32.90           C
ATOM      9  CE1 PHE L   2     202.829 239.908 269.166  1.00 32.90           C
ATOM     10  CE2 PHE L   2     200.799 239.290 268.089  1.00 32.90           C
ATOM     11  CZ  PHE L   2     201.803 240.201 268.300  1.00 32.90           C
ATOM     12  N   SER L   3     199.357 236.979 272.328  1.00 34.28           N
ATOM     13  CA  SER L   3     198.124 237.714 272.618  1.00 34.28           C
ATOM     14  C   SER L   3     197.907 238.857 271.631  1.00 34.28           C
ATOM     15  O   SER L   3     198.129 240.024 271.974  1.00 34.28           O
ATOM     16  CB  SER L   3     198.121 238.241 274.053  1.00 35.83           C
ATOM     17  OG  SER L   3     198.875 239.432 274.158  1.00 35.83           O
ATOM     18  N   PRO L   4     197.482 238.564 270.395  1.00 33.41           N
ATOM     19  CA  PRO L   4     197.440 239.596 269.355  1.00 33.41           C
ATOM     20  C   PRO L   4     196.353 240.646 269.519  1.00 33.41           C
ATOM     21  O   PRO L   4     196.212 241.494 268.633  1.00 33.41           O
ATOM     22  CB  PRO L   4     197.209 238.776 268.079  1.00 33.14           C
ATOM     23  CG  PRO L   4     196.474 237.597 268.538  1.00 33.14           C
ATOM     24  CD  PRO L   4     197.000 237.268 269.901  1.00 33.14           C
TER
END
""")

if __name__ == '__main__':
  exercise_1()
