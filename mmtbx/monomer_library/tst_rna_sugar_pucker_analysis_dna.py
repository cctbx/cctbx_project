from __future__ import absolute_import, division, print_function
"""Tests for the DNA path through rna_sugar_pucker_analysis.evaluate.

Companion to tst_rna_sugar_pucker_analysis.py, which covers the RNA path. The two
guarantees worth pinning down here are that is_dna=False changes nothing, and that
is_dna=True applies a genuinely different RULE rather than merely different constants.

Coordinates are lifted from real deposited structures rather than invented, so each
case is a geometry that actually occurs. The PDB id and residue are named in each
docstring so a failure can be traced back to something inspectable.
"""

from mmtbx.monomer_library import rna_sugar_pucker_analysis
from libtbx.test_utils import approx_equal
import sys
from six.moves import zip


class atom(object):
  def __init__(O, xyz):
    O.xyz = xyz


def atom_dict(xyz_list):
  assert len(xyz_list) == 5
  return dict([(key, atom(xyz))
    for key, xyz in zip(["C1'", "O3'", "C3'", "C4'", "C5'"], xyz_list)])


# 3rmp E/13 DA. P-perp 1.84 implies C2'-endo, delta 80.4 sits squarely in the
# C3'-endo window. A flat contradiction, and far enough from the P-perp cutoff that
# the confidence margin does not suppress it.
CONTRADICTION = (
  [(11.329, -1.856, -47.031), (9.567, -5.032, -47.521), (10.601, -4.032, -47.803),
   (11.827, -4.093, -46.893), (12.854, -5.154, -47.253)],
  (11.680, -0.659, -47.797), (8.974, -5.471, -46.078))

# 6bdb B/14 DC. delta 120.4 falls in the gap between the two DNA windows, belonging
# to neither. This is the case the DNA rule exists for: the RNA thresholds call it an
# outlier, but a residue whose delta is merely intermediate is not evidence of an
# error, and roughly a third of DNA sugars live in that region.
DELTA_IN_GAP = (
  [(9.277, 76.633, 37.997), (9.518, 73.519, 38.201), (8.687, 74.506, 38.815),
   (7.639, 75.078, 37.818), (6.193, 74.879, 38.236)],
  (9.645, 78.045, 38.204), (8.924, 72.242, 37.437))

# 4nqa F/512 DG. A real contradiction (delta 104.7 against a C2'-endo P-perp) but
# P-perp is 2.899, a thousandth of an Angstrom from the 2.9 cutoff, so the pucker
# call it contradicts is not trustworthy in the first place.
INSIDE_MARGIN = (
  [(0.191, 81.271, -11.295), (-1.759, 78.794, -12.455), (-0.409, 79.242, -12.438),
   (0.198, 78.985, -11.063), (1.157, 77.811, -11.003)],
  (1.076, 82.427, -11.409), (-2.586, 78.769, -13.833))


def evaluate(params, case, is_dna):
  sugar, outbound, next_p = case
  return rna_sugar_pucker_analysis.evaluate(
    params=params,
    residue_1_deoxy_ribo_atom_dict=atom_dict(sugar),
    residue_1_c1p_outbound_atom=atom(outbound),
    residue_2_p_atom=atom(next_p),
    is_dna=is_dna)


def exercise_defaults_unchanged():
  """is_dna must default to False, and the RNA constants must be untouched.

  This is the backward-compatibility guarantee. pdb_interpretation builds restraints
  through this same function, so a DNA threshold leaking into the default path would
  silently change refinement behaviour for every RNA structure.
  """
  params = rna_sugar_pucker_analysis.master_phil.fetch().extract()
  assert approx_equal(params.delta_range_2p_min, 129.0)
  assert approx_equal(params.delta_range_2p_max, 162.0)
  assert approx_equal(params.delta_range_3p_min, 65.0)
  assert approx_equal(params.delta_range_3p_max, 104.0)
  assert approx_equal(params.o3p_distance_c1p_outbound_line_2p_max, 2.4)
  assert approx_equal(params.p_distance_c1p_outbound_line_2p_max, 2.9)
  # calling without is_dna must equal calling with is_dna=False
  implicit = rna_sugar_pucker_analysis.evaluate(
    params=params,
    residue_1_deoxy_ribo_atom_dict=atom_dict(DELTA_IN_GAP[0]),
    residue_1_c1p_outbound_atom=atom(DELTA_IN_GAP[1]),
    residue_2_p_atom=atom(DELTA_IN_GAP[2]))
  explicit = evaluate(params, DELTA_IN_GAP, is_dna=False)
  assert implicit.is_delta_outlier == explicit.is_delta_outlier
  assert approx_equal(implicit.delta, explicit.delta)
  print("    exercise_defaults_unchanged: OK")


def exercise_dna_thresholds_present():
  params = rna_sugar_pucker_analysis.master_phil.fetch().extract()
  assert approx_equal(params.dna_delta_range_2p_min, 123.0)
  assert approx_equal(params.dna_delta_range_2p_max, 162.0)
  assert approx_equal(params.dna_delta_range_3p_min, 62.0)
  assert approx_equal(params.dna_delta_range_3p_max, 110.0)
  assert approx_equal(params.dna_o3p_distance_c1p_outbound_line_2p_max, 1.95)
  assert approx_equal(params.dna_pucker_confidence_margin, 0.8)
  # the two windows must not meet, or a residue could satisfy both at once
  assert params.dna_delta_range_3p_max < params.dna_delta_range_2p_min
  print("    exercise_dna_thresholds_present: OK")


def exercise_contradiction_is_flagged():
  """A delta squarely in the opposite window is an outlier under either path."""
  params = rna_sugar_pucker_analysis.master_phil.fetch().extract()
  ana = evaluate(params, CONTRADICTION, is_dna=True)
  assert approx_equal(ana.delta, 80.432, eps=1.e-2)
  assert approx_equal(ana.p_distance_c1p_outbound_line, 1.835, eps=1.e-2)
  assert ana.is_delta_outlier is True
  # and the RNA path agrees here, so this case does not by itself justify is_dna
  assert evaluate(params, CONTRADICTION, is_dna=False).is_delta_outlier is True
  print("    exercise_contradiction_is_flagged: OK")


def exercise_intermediate_delta_not_flagged():
  """The case the DNA rule exists for: RNA thresholds flag it, DNA must not.

  delta 120.4 belongs to neither DNA window. Under the RNA path that counts as
  disagreeing with the P-perp call and is reported as an outlier. Under the DNA path
  it is merely uninformative, which is the correct reading, because deoxyribose
  populates the region between the two puckers far more heavily than ribose does.
  """
  params = rna_sugar_pucker_analysis.master_phil.fetch().extract()
  ana = evaluate(params, DELTA_IN_GAP, is_dna=True)
  assert approx_equal(ana.delta, 120.397, eps=1.e-2)
  assert ana.is_delta_outlier is False
  assert evaluate(params, DELTA_IN_GAP, is_dna=False).is_delta_outlier is True
  print("    exercise_intermediate_delta_not_flagged: OK")


def exercise_confidence_margin():
  """A contradiction is suppressed when P-perp is too close to its own cutoff."""
  params = rna_sugar_pucker_analysis.master_phil.fetch().extract()
  ana = evaluate(params, INSIDE_MARGIN, is_dna=True)
  assert approx_equal(ana.p_distance_c1p_outbound_line, 2.899, eps=1.e-2)
  assert ana.is_delta_outlier is False, \
    "P-perp is 0.001 A from the cutoff; the pucker call it contradicts is not trustworthy"
  # drop the margin and the same residue becomes an outlier, which shows the margin
  # is what suppressed it rather than the geometry being unremarkable
  params.dna_pucker_confidence_margin = 0.0
  assert evaluate(params, INSIDE_MARGIN, is_dna=True).is_delta_outlier is True
  print("    exercise_confidence_margin: OK")


def exercise_threshold_fallback():
  """A DNA threshold set to None falls back to the RNA value, it does not disable.

  The dna_* parameters are overrides layered on the RNA defaults, so clearing one
  means "use the general value here", not "skip this test". Pinned down because the
  two readings are both plausible and the difference is silent: disabling would make
  is_2p_delta None, whereas falling back leaves it a real boolean.
  """
  params = rna_sugar_pucker_analysis.master_phil.fetch().extract()
  params.dna_delta_range_2p_min = None
  params.dna_delta_range_2p_max = None
  ana = evaluate(params, DELTA_IN_GAP, is_dna=True)
  assert ana.delta is not None
  assert ana.is_2p_delta is not None, "cleared DNA threshold should fall back, not disable"
  # delta 120.4 is outside the RNA 2p window (129-162), which is what it now uses
  assert ana.is_2p_delta is False
  # clearing BOTH the DNA and RNA values is what actually disables the test
  params.delta_range_2p_min = None
  params.delta_range_2p_max = None
  assert evaluate(params, DELTA_IN_GAP, is_dna=True).is_2p_delta is None
  print("    exercise_threshold_fallback: OK")


def exercise():
  exercise_defaults_unchanged()
  exercise_dna_thresholds_present()
  exercise_contradiction_is_flagged()
  exercise_intermediate_delta_not_flagged()
  exercise_confidence_margin()
  exercise_threshold_fallback()
  print("OK")


if (__name__ == "__main__"):
  assert len(sys.argv[1:]) == 0
  exercise()
