from __future__ import absolute_import, division, print_function
"""Tests for DNA handling in rna_validate.rna_puckers.

Until recently rna_puckers skipped every DNA residue, via `if (not ra1.is_rna)`, where
is_rna is defined as "has an O2' atom". These cover the lifted gate and the two guards
that lifting it turned out to need.

The fragment below is DNA chain B residues 12-15 of 6bdb, chosen because residue 14 has
a delta of 120.4 that lands between the two pucker windows, the intermediate case that
motivates the DNA rule.
"""

from mmtbx.validation import rna_validate
from iotbx.data_manager import DataManager
from libtbx.test_utils import approx_equal
import json
import sys

# 6bdb chain B, residues 12-15. Real coordinates; residue 14 (DC) has delta 120.4.
dna_fragment = """\
ATOM      1  P    DT B  12       7.703  84.467  28.653  1.00 28.49           P
ATOM      2  OP1  DT B  12       6.906  84.630  27.419  1.00 42.85           O
ATOM      3  OP2  DT B  12       7.902  85.593  29.585  1.00 35.28           O
ATOM      4  O5'  DT B  12       7.092  83.234  29.440  1.00 35.19           O
ATOM      5  C5'  DT B  12       6.559  82.196  28.682  1.00 34.76           C
ATOM      6  C4'  DT B  12       7.139  80.850  29.070  1.00 23.05           C
ATOM      7  O4'  DT B  12       8.569  80.882  29.224  1.00 23.66           O
ATOM      8  C3'  DT B  12       6.603  80.238  30.355  1.00 22.66           C
ATOM      9  O3'  DT B  12       6.158  78.993  29.992  1.00 20.03           O
ATOM     10  C2'  DT B  12       7.847  80.135  31.265  1.00 20.20           C
ATOM     11  C1'  DT B  12       8.912  79.941  30.210  1.00 19.16           C
ATOM     12  N1   DT B  12      10.317  80.209  30.588  1.00 20.50           N
ATOM     13  C2   DT B  12      11.183  79.164  30.726  1.00 21.62           C
ATOM     14  O2   DT B  12      10.836  78.005  30.659  1.00 17.26           O
ATOM     15  N3   DT B  12      12.475  79.514  30.981  1.00 18.84           N
ATOM     16  C4   DT B  12      12.990  80.794  31.070  1.00 21.42           C
ATOM     17  O4   DT B  12      14.175  81.004  31.291  1.00 23.24           O
ATOM     18  C5   DT B  12      12.035  81.849  30.888  1.00 21.59           C
ATOM     19  C7   DT B  12      12.473  83.282  30.980  1.00 24.39           C
ATOM     20  C6   DT B  12      10.764  81.517  30.643  1.00 20.26           C
ATOM     21  P    DT B  13       4.827  78.389  30.617  1.00 17.89           P
ATOM     22  OP1  DT B  13       4.292  77.461  29.590  1.00 21.53           O
ATOM     23  OP2  DT B  13       3.986  79.469  31.160  1.00 19.28           O
ATOM     24  O5'  DT B  13       5.389  77.539  31.843  1.00 16.21           O
ATOM     25  C5'  DT B  13       6.073  76.325  31.578  1.00 16.27           C
ATOM     26  C4'  DT B  13       6.815  75.858  32.811  1.00 14.84           C
ATOM     27  O4'  DT B  13       7.989  76.695  33.029  1.00 15.23           O
ATOM     28  C3'  DT B  13       5.998  75.942  34.107  1.00 12.73           C
ATOM     29  O3'  DT B  13       6.161  74.720  34.859  1.00 15.24           O
ATOM     30  C2'  DT B  13       6.616  77.169  34.820  1.00 11.70           C
ATOM     31  C1'  DT B  13       8.058  77.014  34.404  1.00 15.74           C
ATOM     32  N1   DT B  13       8.910  78.235  34.542  1.00 14.99           N
ATOM     33  C2   DT B  13      10.263  78.063  34.459  1.00 15.60           C
ATOM     34  O2   DT B  13      10.771  76.978  34.308  1.00 17.63           O
ATOM     35  N3   DT B  13      11.007  79.203  34.566  1.00 14.36           N
ATOM     36  C4   DT B  13      10.547  80.479  34.747  1.00 16.78           C
ATOM     37  O4   DT B  13      11.321  81.445  34.821  1.00 19.23           O
ATOM     38  C5   DT B  13       9.101  80.615  34.823  1.00 14.78           C
ATOM     39  C7   DT B  13       8.489  81.985  35.006  1.00 16.40           C
ATOM     40  C6   DT B  13       8.347  79.496  34.721  1.00 16.03           C
ATOM     41  P    DC B  14       5.053  74.294  35.934  1.00 12.44           P
ATOM     42  OP1  DC B  14       5.467  72.903  36.266  1.00 12.83           O
ATOM     43  OP2  DC B  14       3.707  74.497  35.385  1.00 12.20           O
ATOM     44  O5'  DC B  14       5.289  75.246  37.187  1.00 12.38           O
ATOM     45  C5'  DC B  14       6.193  74.879  38.236  1.00 12.48           C
ATOM     46  C4'  DC B  14       7.639  75.078  37.818  1.00 12.02           C
ATOM     47  O4'  DC B  14       7.915  76.504  37.748  1.00 12.04           O
ATOM     48  C3'  DC B  14       8.687  74.506  38.815  1.00 14.31           C
ATOM     49  O3'  DC B  14       9.518  73.519  38.201  1.00 16.27           O
ATOM     50  C2'  DC B  14       9.529  75.730  39.201  1.00 15.61           C
ATOM     51  C1'  DC B  14       9.277  76.633  37.997  1.00 12.72           C
ATOM     52  N1   DC B  14       9.645  78.045  38.204  1.00 13.61           N
ATOM     53  C2   DC B  14      10.956  78.401  37.953  1.00 16.45           C
ATOM     54  O2   DC B  14      11.733  77.528  37.556  1.00 16.25           O
ATOM     55  N3   DC B  14      11.346  79.683  38.145  1.00 14.70           N
ATOM     56  C4   DC B  14      10.476  80.590  38.579  1.00 15.62           C
ATOM     57  N4   DC B  14      10.930  81.847  38.744  1.00 16.34           N
ATOM     58  C5   DC B  14       9.115  80.251  38.871  1.00 14.53           C
ATOM     59  C6   DC B  14       8.737  78.971  38.670  1.00 13.16           C
ATOM     60  P    DG B  15       8.924  72.242  37.437  1.00 14.46           P
ATOM     61  OP1  DG B  15       8.738  72.612  36.023  1.00 16.47           O
ATOM     62  OP2  DG B  15       7.767  71.735  38.190  1.00 14.01           O
ATOM     63  O5'  DG B  15      10.119  71.199  37.443  1.00 16.57           O
ATOM     64  C5'  DG B  15      10.491  70.495  38.645  1.00 14.32           C
ATOM     65  C4'  DG B  15      11.999  70.524  38.796  1.00 13.50           C
ATOM     66  O4'  DG B  15      12.429  71.870  39.141  1.00 13.53           O
ATOM     67  C3'  DG B  15      12.585  69.639  39.870  1.00 13.83           C
ATOM     68  O3'  DG B  15      13.946  69.369  39.471  1.00 15.90           O
ATOM     69  C2'  DG B  15      12.504  70.561  41.103  1.00 13.83           C
ATOM     70  C1'  DG B  15      12.827  71.922  40.502  1.00 13.31           C
ATOM     71  N9   DG B  15      12.114  73.045  41.090  1.00 13.97           N
ATOM     72  C8   DG B  15      10.803  73.075  41.512  1.00 13.35           C
ATOM     73  N7   DG B  15      10.416  74.260  41.911  1.00 12.35           N
ATOM     74  C5   DG B  15      11.527  75.080  41.703  1.00 14.04           C
ATOM     75  C6   DG B  15      11.701  76.475  41.921  1.00 13.47           C
ATOM     76  O6   DG B  15      10.886  77.302  42.365  1.00 13.44           O
ATOM     77  N1   DG B  15      12.977  76.894  41.559  1.00 13.02           N
ATOM     78  C2   DG B  15      13.964  76.070  41.059  1.00 15.03           C
ATOM     79  N2   DG B  15      15.137  76.653  40.786  1.00 16.73           N
ATOM     80  N3   DG B  15      13.814  74.766  40.846  1.00 14.08           N
ATOM     81  C4   DG B  15      12.567  74.349  41.172  1.00 14.45           C
TER
END
"""


def hierarchy_from(pdb_str):
  dm = DataManager()
  dm.process_model_str("1", pdb_str)
  return dm.get_model("1").get_hierarchy()


def exercise_is_dna_resname():
  assert rna_validate.is_dna_resname("DA")
  assert rna_validate.is_dna_resname(" DT ")
  assert rna_validate.is_dna_resname("dc")
  assert rna_validate.is_dna_resname("DI")
  # RNA residues and amino acids must not pass
  assert not rna_validate.is_dna_resname("A")
  assert not rna_validate.is_dna_resname("U")
  assert not rna_validate.is_dna_resname("ALA")
  assert not rna_validate.is_dna_resname("")
  print("    exercise_is_dna_resname: OK")


def exercise_dna_is_evaluated():
  """DNA residues reach the pucker analysis instead of being skipped.

  Before the gate was lifted this returned nothing at all for DNA, which is why
  MolProbity reported no pucker validation for DNA-only structures.
  """
  h = hierarchy_from(dna_fragment)
  puckers = rna_validate.rna_puckers(pdb_hierarchy=h, outliers_only=False)
  assert puckers.n_total > 0, "DNA residues are being skipped; the is_rna gate is back"
  # All four residues get a delta, since delta is C5'-C4'-C3'-O3' and needs nothing
  # from the following residue. Residue 15 is last, so it has no P-perp and therefore
  # no verdict, but it is still counted.
  assert puckers.n_total == 4, "expected 4 evaluable residues, got %d" % puckers.n_total
  by_seq = dict((r.get("resseq").strip(), r)
                for r in json.loads(puckers.as_JSON())["flat_results"])
  assert approx_equal(by_seq["12"]["delta_angle"], 125.294, eps=1.e-2)
  assert approx_equal(by_seq["14"]["delta_angle"], 120.397, eps=1.e-2)
  assert by_seq["15"]["outlier"] is None, \
    "the last residue has no following phosphate, so it cannot have a verdict"
  print("    exercise_dna_is_evaluated: OK")


def exercise_include_dna_switch():
  """include_dna=False restores the pre-change behaviour exactly."""
  h = hierarchy_from(dna_fragment)
  off = rna_validate.rna_puckers(pdb_hierarchy=h, outliers_only=False,
                                 include_dna=False)
  assert off.n_total == 0, \
    "include_dna=False must skip DNA, got %d results" % off.n_total
  on = rna_validate.rna_puckers(pdb_hierarchy=h, outliers_only=False)
  assert on.n_total == 4, "include_dna should default to True"
  print("    exercise_include_dna_switch: OK")


def exercise_intermediate_delta_not_an_outlier():
  """6bdb B/14 has delta 120.4, between the windows, and must not be flagged.

  This is the end-to-end version of the rule tested in
  tst_rna_sugar_pucker_analysis_dna: it confirms is_dna is actually threaded through
  from rna_puckers, not merely available on evaluate().
  """
  h = hierarchy_from(dna_fragment)
  puckers = rna_validate.rna_puckers(pdb_hierarchy=h, outliers_only=False)
  assert puckers.n_outliers == 0, \
    "intermediate-delta DNA should not be flagged, got %d outliers" % puckers.n_outliers
  print("    exercise_intermediate_delta_not_an_outlier: OK")


def exercise_incomplete_sugar_does_not_crash():
  """A DNA residue missing C1' must be skipped, not raise.

  evaluate() indexes the sugar atom dict directly, so a missing atom raises KeyError.
  The RNA path never reached this because is_rna already required an O2'. Admitting
  DNA exposed it, and 24 of 420 benchmark structures crashed before the guard.
  """
  stripped = "\n".join(line for line in dna_fragment.splitlines()
                       if not (line.startswith("ATOM") and
                               line[12:16].strip() == "C1'" and
                               line[22:26].strip() == "14"))
  h = hierarchy_from(stripped + "\n")
  puckers = rna_validate.rna_puckers(pdb_hierarchy=h, outliers_only=False)
  # residue 14 is dropped, the other three still evaluate
  assert puckers.n_total == 3, \
    "expected the incomplete residue to be skipped, got %d" % puckers.n_total
  seqs = set(r.get("resseq").strip()
             for r in json.loads(puckers.as_JSON())["flat_results"])
  assert "14" not in seqs, "the residue missing C1' should be absent, not guessed at"
  print("    exercise_incomplete_sugar_does_not_crash: OK")


def exercise_rna_unaffected():
  """The RNA path must be untouched by any of this.

  A ribose residue still fails is_dna_resname, so it can never pick up DNA
  thresholds even though it now travels through the same branch.
  """
  rna_like = dna_fragment.replace(" DT B", "  U B").replace(" DC B", "  C B") \
                         .replace(" DG B", "  G B")
  h = hierarchy_from(rna_like)
  puckers = rna_validate.rna_puckers(pdb_hierarchy=h, outliers_only=False)
  # no O2' present, so these are not real RNA and residue_analysis rejects them;
  # the point is that they are not silently rescored under DNA thresholds
  assert puckers.n_total == 0
  print("    exercise_rna_unaffected: OK")


def exercise():
  exercise_is_dna_resname()
  exercise_dna_is_evaluated()
  exercise_include_dna_switch()
  exercise_intermediate_delta_not_an_outlier()
  exercise_incomplete_sugar_does_not_crash()
  exercise_rna_unaffected()
  print("OK")


if (__name__ == "__main__"):
  assert len(sys.argv[1:]) == 0
  exercise()
