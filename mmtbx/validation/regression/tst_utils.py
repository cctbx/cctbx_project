from __future__ import absolute_import, division, print_function
from libtbx.test_utils import approx_equal
from mmtbx.validation.utils import (
  _clash_severity,
  _cbeta_severity,
  _bond_angle_severity,
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
  assert approx_equal(_cbeta_severity(0.25), 1.44, eps=0.01)
  assert approx_equal(_cbeta_severity(0.50), 4.44, eps=0.01)
  # At baseline (0.13), severity is zero
  assert _cbeta_severity(0.13) == 0.0
  # Below baseline returns zero
  assert _cbeta_severity(0.10) == 0.0

  print("  exercise_cbeta_severity: OK")


def exercise_bond_angle_severity():
  # With sigma values
  assert approx_equal(_bond_angle_severity(1, 4.0), 1.0, eps=0.01)
  assert approx_equal(_bond_angle_severity(1, 6.0), 2.0, eps=0.01)
  # Multiple outliers add 0.5 each
  assert approx_equal(_bond_angle_severity(3, 6.0), 3.0, eps=0.01)

  # Fallback without sigma
  assert approx_equal(_bond_angle_severity(1, None), 1.0, eps=0.01)
  assert approx_equal(_bond_angle_severity(3, None), 2.0, eps=0.01)

  print("  exercise_bond_angle_severity: OK")


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


def exercise():
  exercise_clash_severity()
  exercise_cbeta_severity()
  exercise_bond_angle_severity()
  exercise_residue_quality_score()
  exercise_ranking_invariants()


if __name__ == "__main__":
  exercise()
  print("OK")
