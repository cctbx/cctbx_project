"""Bond-valence-sum oxidation-state estimate from the vendored bvparm2020.cif
[35] (plan 6.4).

For a metal with a set of candidate oxidation states and its observed
first-shell donor distances, this computes the Brown-Altermatt bond-valence sum
V_q = sum_i exp((R0(q) - R_i) / B) [27] under each candidate q and returns the q
whose V_q is closest to q itself: the most self-consistent oxidation state for
the observed geometry, plus the per-state V map. This is the same estimate as
the calibration BVS in ~/qmi/calibration/prepare_qm_calcs.py.

It is a CROSS-CHECK, advisory only: BVS uses the observed crystal bond lengths,
so it is reliable only at high resolution and for an undamaged metal, and the
crystallographic state may be X-ray-photoreduced away from the biological one
[42]. qmi reports the BVS estimate alongside the ligand-field/geometry guess in
electronic.py; it does not override it.

Parameters come from bvparm2020.cif [35] (no mmtbx.ions), parsed with iotbx.cif,
keyed on (cation, oxidation, donor element), first/most-reliable row per combo
(adapted from the calibration parser). Caveat noted there: Fe(II)/Fe(III) R0 are
close and S-donor rows sparse, so BVS barely separates oxidation states at
thiolate sites.
"""

import os
from math import exp
from functools import lru_cache

import iotbx.cif
from libtbx import group_args

from mmtbx.geometry_restraints.qmi import metals


_CIF_PATH = os.path.join(os.path.dirname(__file__), "bvparm2020.cif")

# Below this separation in |V - q| between the best and second-best candidate,
# the BVS verdict is not stable against re-refinement of the same data and is
# reported as a close call. Calibrated on the benchmark's five Fe sites, where
# margins of 0.10-0.32 flip between refinements of one crystal and only the
# 0.73 case is stable [84].
CLOSE_CALL_MARGIN = 0.35


@lru_cache(maxsize=1)
def _table():
  """Parse bvparm2020.cif [35] into {(cation, oxidation, donor): (R0, B)}.

  Rows are listed in decreasing reliability, so the FIRST row for a combo wins,
  with one exception: the file's own header warns that "parameters marked at
  'unchecked' are unreliable", and 28 combos lead with an unchecked row while a
  checked alternative follows. Those prefer the first CHECKED row. Two of the
  28 are load-bearing for qmi: Fe(2+)-S (2.12 unchecked vs 2.125 checked, which
  sets the rubredoxin-type estimate) and Cu(1+)-O (1.610 unchecked vs 1.601,
  and a 1.504 row from transition-metal complexes).
  """
  model = iotbx.cif.reader(file_path=_CIF_PATH).model()
  if model is None:
    raise RuntimeError(f"could not parse bvparm CIF: {_CIF_PATH}")
  blk = list(model.values())[0]
  details = blk.get("_valence_param_details", None)
  rows = zip(blk["_valence_param_atom_1"], blk["_valence_param_atom_1_valence"],
             blk["_valence_param_atom_2"], blk["_valence_param_Ro"],
             blk["_valence_param_B"],
             details if details is not None
             else [""] * len(blk["_valence_param_atom_1"]))
  ordered = {}  # key -> [(params, unchecked), ...] in file order
  for cation, valence, anion, ro, b, det in rows:
    try:
      key = (metals.normalize_element(cation), int(valence),
             metals.normalize_element(anion))
      params = (float(ro), float(b))
    except ValueError:
      continue
    ordered.setdefault(key, []).append(
      (params, "unchecked" in str(det).lower()))
  table = {}
  for key, cands in ordered.items():
    checked = [p for p, unchecked in cands if not unchecked]
    table[key] = checked[0] if checked else cands[0][0]
  return table


def valence_params(cation, oxidation, donor):
  """(R0, B) for a (cation, oxidation, donor-element) bond, or None if the
  vendored table has no entry."""
  return _table().get((metals.normalize_element(cation), oxidation,
                       metals.normalize_element(donor)))


def valence_sum(cation, oxidation, contacts):
  """Brown-Altermatt bond-valence sum V = sum_i exp((R0 - R_i)/B) [27] over the
  first-shell *contacts* (``(donor_element, distance_A)`` tuples) under one
  oxidation state. Returns ``(V, n_used)`` (n_used = donors with a parameter);
  V is None if none matched."""
  sv = [exp((p[0] - dist) / p[1])
        for donor, dist in contacts
        for p in [valence_params(cation, oxidation, donor)]
        if p is not None]
  if not sv:
    return None, 0
  return sum(sv), len(sv)


def estimate_oxidation(cation, ox_states, contacts):
  """BVS estimate of the most likely oxidation state for the observed geometry
  (plan 6.4), over the candidate *ox_states*.

  For each candidate q, computes V_q = valence_sum(cation, q, contacts) and
  returns the q that minimises ``|V_q - q|`` (the most self-consistent state),
  with the per-state V map. Returns group_args(best, valences={q: V_q}, n_used),
  or None if no candidate had vendored parameters. Advisory cross-check
  (observed bond lengths, reliable only at high resolution / undamaged metal
  [42]); it does not override the electronic guess.
  """
  valences, counts = {}, {}
  for q in ox_states:
    V, n = valence_sum(cation, q, contacts)
    if V is not None:
      valences[q] = V
      counts[q] = n
  if not valences:
    return None
  n_used = max(counts.values())

  # Asymmetric parameter coverage: if the candidate states were summed over
  # DIFFERENT numbers of bonds, the comparison is meaningless, because a state
  # missing a parameter silently loses that bond's valence and is pushed down.
  # Live case: Mn(2+) has S parameters and Mn(3+) does not, so every Mn thiolate
  # site would be forced to Mn(II). Refuse rather than report a biased winner.
  if len(set(counts.values())) > 1:
    return group_args(best=None, valences=valences, n_used=n_used,
                      n_by_state=counts, margin=None,
                      coverage_asymmetric=True, close_call=False)

  best = min(valences, key=lambda q: abs(valences[q] - q))
  # Margin between the winner and the runner-up, in |V - q|. The discrimination
  # here is ~0.01-0.06 A of mean bond length, smaller than the spread between
  # re-refinements of the same crystal [84], so a small margin means the verdict
  # is a property of the refinement rather than of the site.
  devs = sorted(abs(valences[q] - q) for q in valences)
  margin = (devs[1] - devs[0]) if len(devs) > 1 else None
  return group_args(best=best, valences=valences, n_used=n_used,
                    n_by_state=counts, margin=margin,
                    coverage_asymmetric=False,
                    close_call=(margin is not None
                                and margin < CLOSE_CALL_MARGIN))
