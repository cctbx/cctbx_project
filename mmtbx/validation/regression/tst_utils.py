from __future__ import absolute_import, division, print_function
from libtbx.test_utils import approx_equal
from mmtbx.validation.utils import (
  _clash_severity,
  _cbeta_severity,
  _bond_angle_severity,
  _omega_twist_severity,
  _rna_suite_severity,
  _rna_pucker_severity,
  calculate_overall_residue_quality_score,
)


def exercise_clash_severity():
  # Single clash at key overlaps (non-linear power curve)
  assert approx_equal(_clash_severity(0.4, 1), 1.39, eps=0.05)
  assert approx_equal(_clash_severity(0.5, 1), 2.02, eps=0.05)
  assert approx_equal(_clash_severity(0.7, 1), 3.42, eps=0.05)
  assert approx_equal(_clash_severity(0.9, 1), 4.98, eps=0.05)
  assert approx_equal(_clash_severity(1.5, 1), 10.30, eps=0.05)

  # At threshold (0.1 A offset), severity should be zero
  assert _clash_severity(0.1, 1) == 0.0

  # Multi-clash log2 bonus
  assert approx_equal(_clash_severity(0.5, 2), 3.02, eps=0.05)   # +1.0
  assert approx_equal(_clash_severity(0.5, 5), 4.34, eps=0.05)   # +2.3
  assert approx_equal(_clash_severity(0.5, 16), 6.02, eps=0.05)  # +4.0 (cap)

  # Cap: 16 and 50 clashes give the same bonus
  assert approx_equal(
    _clash_severity(0.5, 50),
    _clash_severity(0.5, 16),
    eps=0.01)

  print("  exercise_clash_severity: OK")


def exercise_cbeta_severity():
  # Logarithmic, anchored at the two values the old linear form documented
  assert approx_equal(_cbeta_severity(0.25), 1.0, eps=0.01)
  assert approx_equal(_cbeta_severity(0.50), 4.0, eps=0.01)
  assert approx_equal(_cbeta_severity(1.00), 7.0, eps=0.01)
  assert approx_equal(_cbeta_severity(2.00), 10.0, eps=0.01)

  # At or below the outlier threshold, the floor. The caller only reaches this
  # function for residues already flagged, so 0 would be the wrong answer.
  assert _cbeta_severity(0.13) == 1.0
  assert _cbeta_severity(0.10) == 1.0
  assert _cbeta_severity(None) == 1.0

  # Monotonic, and bounded in practice. The largest deviation in 139,418 measured
  # C-beta evaluations was 2.93 A; the old linear form scored that 33.6, more than
  # twice the ceiling of every discrete tier in the metric.
  vals = [_cbeta_severity(d / 100.0) for d in range(25, 300, 5)]
  assert all(b >= a for a, b in zip(vals, vals[1:])), "cbeta severity must be monotonic"
  assert _cbeta_severity(2.93) < 15.0, "worst observed C-beta must stay under the ceiling"

  print("  exercise_cbeta_severity: OK")


def exercise_bond_angle_severity():
  # Logarithmic, anchored at the two values the old linear form documented
  assert approx_equal(_bond_angle_severity(1, 4.0), 1.0, eps=0.01)
  assert approx_equal(_bond_angle_severity(1, 10.0), 4.0, eps=0.01)

  # Sub-threshold sigma cannot drop below the floor
  assert approx_equal(_bond_angle_severity(1, 2.0), 1.0, eps=0.01)
  # Sign of the deviation is irrelevant; only the magnitude counts
  assert approx_equal(_bond_angle_severity(1, -10.0), _bond_angle_severity(1, 10.0))

  # Count bonus is log2 capped at 4.0, matching _clash_severity. It used to be
  # +0.5 per extra outlier with no cap, the same unbounded growth in a second term.
  assert approx_equal(_bond_angle_severity(2, 4.0), 2.0, eps=0.01)
  assert approx_equal(_bond_angle_severity(16, 4.0), 5.0, eps=0.01)
  assert approx_equal(_bond_angle_severity(64, 4.0), _bond_angle_severity(16, 4.0)), \
      "count bonus must be capped"

  # Fallback without sigma: the floor, plus the count bonus
  assert approx_equal(_bond_angle_severity(1, None), 1.0, eps=0.01)
  assert approx_equal(_bond_angle_severity(4, None), 3.0, eps=0.01)

  # The whole point of the change: the far tail must stay under the metric's ceiling.
  # The worst bond in 1.2M restraints was 164 sigma, which the old linear form scored 81.
  assert _bond_angle_severity(1, 164.0) < 15.0, \
      "worst observed bond deviation must stay under the twisted-peptide ceiling"
  assert _bond_angle_severity(1, 31.1) < _omega_twist_severity(90.0), \
      "a single ligand-link bond outlier must not outrank a perpendicular peptide"

  print("  exercise_bond_angle_severity: OK")


def exercise_cbeta_chirality_suppression():
  """A handedness swap and its C-beta deviation are one fact, not two findings.

  98.9% of C-beta outliers past 1.8 A already carry a chirality flag, so scoring both
  gave those residues 10.0 for the swap plus 25-33 on top. Only the typed handedness
  count suppresses: tetrahedral geometry is a different problem, and the untyped
  fallback cannot tell the two apart.
  """
  base = {'ramalyze_type': 'not_evaluated', 'rotalyze_category': 'not_evaluated',
          'cablam_outlier_type': 'not_evaluated', 'omega_type': 'not_evaluated',
          'is_cbeta_outlier': True, 'cbeta_deviation': 2.2}

  d = dict(base)
  alone = calculate_overall_residue_quality_score(d)

  d = dict(base, num_chiral_handedness_res=1)
  suppressed = calculate_overall_residue_quality_score(d)
  assert approx_equal(suppressed, 10.0, eps=0.01), \
      "handedness must suppress C-beta entirely, got %s" % suppressed
  assert suppressed < alone + 10.0, "C-beta must not be added on top of handedness"

  # Tetrahedral does NOT suppress: different problem, both should count
  d = dict(base, num_chiral_tetrahedral_res=1)
  tetra = calculate_overall_residue_quality_score(d)
  assert tetra > 5.0, "tetrahedral must not suppress C-beta, got %s" % tetra

  # The untyped fallback does not suppress either
  d = dict(base, num_chiral_outliers_res=1)
  untyped = calculate_overall_residue_quality_score(d)
  assert untyped > 5.0, "untyped chirality must not suppress C-beta, got %s" % untyped

  # A C-beta outlier with no chirality flag is untouched.
  # eps allows for the score being rounded to 1 decimal place by the caller.
  assert approx_equal(alone, _cbeta_severity(2.2), eps=0.06)

  print("  exercise_cbeta_chirality_suppression: OK")


def _make_residue(**overrides):
  """Build a residue_data dict with sensible defaults."""
  data = {
    'ramalyze_type': 'general',
    'ramalyze_category': 'favored',
    'rotalyze_category': 'favored',
    'is_glycine': False,
    'is_cbeta_outlier': False,
    'cbeta_deviation': 0.0,
    'cablam_outlier_type': 'favored',
    'omega_type': 'trans',
    'is_proline': False,
    'num_bad_clashes_res': 0,
    'worst_clash_overlap': 0.0,
    'num_bond_outliers_res': 0,
    'worst_bond_deviation': None,
    'num_angle_outliers_res': 0,
    'worst_angle_deviation': None,
    'num_chiral_handedness_res': 0,
    'num_chiral_tetrahedral_res': 0,
    'num_chiral_pseudochiral_res': 0,
    'num_chiral_outliers_res': 0,
    'is_rna_residue': False,
  }
  data.update(overrides)
  return data


def exercise_residue_quality_score():
  score = calculate_overall_residue_quality_score

  # Clean residue: all metrics present, no outliers -> 0.0
  assert score(_make_residue()) == 0.0

  # No applicable metrics -> None (glycine with everything not_evaluated)
  assert score({
    'ramalyze_type': 'not_applicable',
    'rotalyze_category': 'not_evaluated',
    'is_glycine': True,
    'cablam_outlier_type': 'not_evaluated',
    'omega_type': 'not_evaluated',
  }) is None

  # Twist only -> 15.0
  assert score(_make_residue(omega_type='twisted')) == 15.0

  # Rama outlier only -> 5.0
  assert score(_make_residue(ramalyze_category='outlier')) == 5.0

  # Rotamer outlier only -> 3.0
  assert score(_make_residue(rotalyze_category='outlier')) == 3.0

  # Single 0.4 A clash -> 1.4
  assert score(_make_residue(
    num_bad_clashes_res=1, worst_clash_overlap=-0.4)) == 1.4

  # Cis non-proline -> 8.0
  assert score(_make_residue(omega_type='cis', is_proline=False)) == 8.0

  # Cis proline is not an outlier -> 0.0
  assert score(_make_residue(omega_type='cis', is_proline=True)) == 0.0

  # Chiral handedness swap -> 10.0
  assert score(_make_residue(num_chiral_handedness_res=1)) == 10.0

  # Multi-issue: twist + rama + rota -> 15 + 0.25*(5+3) = 17.0
  assert score(_make_residue(
    omega_type='twisted',
    ramalyze_category='outlier',
    rotalyze_category='outlier')) == 17.0

  # 15 clashes at 1.54 A -> ~14.6
  assert score(_make_residue(
    num_bad_clashes_res=15, worst_clash_overlap=-1.54)) == 14.6

  print("  exercise_residue_quality_score: OK")


def exercise_ranking_invariants():
  score = calculate_overall_residue_quality_score

  twist_only = score(_make_residue(omega_type='twisted'))
  big_clash = score(_make_residue(
    num_bad_clashes_res=1, worst_clash_overlap=-1.5))
  rama_only = score(_make_residue(ramalyze_category='outlier'))
  rota_only = score(_make_residue(rotalyze_category='outlier'))
  small_clash = score(_make_residue(
    num_bad_clashes_res=1, worst_clash_overlap=-0.4))

  # Twisted peptide > single large clash > rama > rotamer > small clash
  assert twist_only > big_clash, \
    "twist (%s) should outrank single 1.5A clash (%s)" % (twist_only, big_clash)
  assert big_clash > rama_only, \
    "1.5A clash (%s) should outrank rama outlier (%s)" % (big_clash, rama_only)
  assert rama_only > rota_only, \
    "rama (%s) should outrank rotamer (%s)" % (rama_only, rota_only)
  assert rota_only > small_clash, \
    "rotamer (%s) should outrank small clash (%s)" % (rota_only, small_clash)

  # Twist + other issues should outrank clash-only with many clashes
  twist_plus = score(_make_residue(
    omega_type='twisted',
    ramalyze_category='outlier',
    rotalyze_category='outlier'))
  many_clashes = score(_make_residue(
    num_bad_clashes_res=18, worst_clash_overlap=-1.54))
  assert twist_plus > many_clashes, \
    "twist+rama+rota (%s) should outrank 18 clashes at 1.54A (%s)" % (
      twist_plus, many_clashes)

  print("  exercise_ranking_invariants: OK")



def exercise_omega_twist_severity():
  # The anchors: MolProbity's Twisted border, and a perpendicular peptide.
  assert approx_equal(_omega_twist_severity(150.0), 3.0), \
    "a peptide 30 deg off trans sits at the border and takes the floor"
  assert approx_equal(_omega_twist_severity(90.0), 15.0), \
    "perpendicular keeps the old flat value"
  # sin^2 in between. 135 deg is 45 deg of twist, exactly half the conjugation lost.
  assert approx_equal(_omega_twist_severity(135.0), 7.0, eps=0.01)
  assert approx_equal(_omega_twist_severity(120.0), 11.0, eps=0.01)
  # Symmetric about perpendicular, and about the sign of omega: the same twist either
  # side of trans, or on the cis side, is the same amount of lost planarity.
  for a, b in ((150.0, -150.0), (150.0, 30.0), (120.0, 60.0), (135.0, -45.0)):
    assert approx_equal(_omega_twist_severity(a), _omega_twist_severity(b)), \
      "twist %s and %s describe the same non-planarity" % (a, b)
  # Monotonic from border to perpendicular.
  vals = [_omega_twist_severity(180.0 - t) for t in range(30, 91, 5)]
  assert vals == sorted(vals), "severity must not fall as the twist grows"
  # Never below the floor, never above the old flat tier, whatever comes in.
  for om in (179.0, 151.0, 150.0, 90.0, 1.0, 0.0, -90.0, -179.0):
    v = _omega_twist_severity(om)
    assert 3.0 <= v <= 15.0, "%s gave %s, outside [3, 15]" % (om, v)
  # No angle available: fall back to the old behaviour rather than under-reporting.
  assert approx_equal(_omega_twist_severity(None), 15.0), \
    "a caller that cannot supply omega must be no worse off than before"
  # A borderline twist must stay above the high-triage cut of 2.0, which is the whole
  # reason the floor is 3.0 and not 0.0.
  assert _omega_twist_severity(150.0) > 2.0
  print("  exercise_omega_twist_severity: OK")


def exercise_rna_suite_severity():
  # A suite matching no named cluster is the serious case.
  assert approx_equal(_rna_suite_severity(True, None), 4.0)
  assert approx_equal(_rna_suite_severity(True, 0.9), 4.0), \
    "an outlier stays an outlier however good its suiteness looks"
  # Assigned to a cluster but fitting it poorly: the "wannabe" tier. This is the
  # reason for grading at all, since a binary flag cannot express it.
  assert approx_equal(_rna_suite_severity(False, 0.1), 1.5)
  assert approx_equal(_rna_suite_severity(False, 0.29), 1.5)
  # A good fit, or no suiteness reported, costs nothing.
  assert approx_equal(_rna_suite_severity(False, 0.3), 0.0), \
    "0.3 is the boundary and must count as an acceptable fit"
  assert approx_equal(_rna_suite_severity(False, 0.95), 0.0)
  assert approx_equal(_rna_suite_severity(False, None), 0.0)
  print("  exercise_rna_suite_severity: OK")


def exercise_rna_pucker_severity():
  # A delta outlier means P-perp and delta disagree about the pucker: a real,
  # rare, chemically well-defined error, at the Ramachandran tier.
  assert approx_equal(_rna_pucker_severity(True, False), 5.0)
  assert approx_equal(_rna_pucker_severity(True, True), 5.0), \
    "delta dominates when both are set"
  # Epsilon out of range while the pucker itself is consistent is milder.
  assert approx_equal(_rna_pucker_severity(False, True), 3.0)
  assert approx_equal(_rna_pucker_severity(False, False), 0.0)
  # None must behave as "not an outlier". The validator reports None for a residue
  # it could not judge, and an unjudged residue must not be scored as a problem.
  assert approx_equal(_rna_pucker_severity(None, None), 0.0)
  # The untyped fallback is treated as the delta tier.
  assert approx_equal(_rna_pucker_severity(False, False, True), 5.0)
  print("  exercise_rna_pucker_severity: OK")


def exercise_rna_severity_ordering():
  """The tiers must stay ordered, since the ranking depends on it.

  Pucker outliers are rarer than suite outliers, which is why they outrank them.
  If someone retunes one constant without the other, this is what catches it.
  """
  pucker_delta = _rna_pucker_severity(True, False)
  pucker_eps = _rna_pucker_severity(False, True)
  suite_out = _rna_suite_severity(True, None)
  suite_wannabe = _rna_suite_severity(False, 0.1)
  assert pucker_delta > suite_out > pucker_eps > suite_wannabe > 0.0, \
    "severity tiers out of order: %s" % [pucker_delta, suite_out,
                                         pucker_eps, suite_wannabe]
  print("  exercise_rna_severity_ordering: OK")


def exercise():
  exercise_clash_severity()
  exercise_cbeta_severity()
  exercise_bond_angle_severity()
  exercise_cbeta_chirality_suppression()
  exercise_residue_quality_score()
  exercise_ranking_invariants()
  exercise_omega_twist_severity()
  exercise_rna_suite_severity()
  exercise_rna_pucker_severity()
  exercise_rna_severity_ordering()


if __name__ == "__main__":
  exercise()
  print("OK")
