"""Bond-valence-sum oxidation-state estimate from the vendored bvparm2020.cif
[35] (plan 6.4).

For a metal with a set of candidate oxidation states and its observed first-shell
donor distances, this computes the Brown-Altermatt bond-valence sum
V_q = sum_i exp((R0(q) - R_i) / B) [27] under each candidate q and returns the q
whose V_q is closest to q itself: the most self-consistent oxidation state for
the observed geometry, plus the per-state V map. This is the same estimate as the
calibration BVS in ~/qmi/calibration/prepare_qm_calcs.py.

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

from __future__ import absolute_import, division, print_function

import os
from math import exp
from functools import lru_cache

import iotbx.cif
from libtbx import group_args

from mmtbx.geometry_restraints.qmi import metals


_CIF_PATH = os.path.join(os.path.dirname(__file__), "bvparm2020.cif")


@lru_cache(maxsize=1)
def _table():
  """Parse bvparm2020.cif [35] into {(cation, oxidation, donor): (R0, B)},
  keeping the first (most reliable) row per combo, via iotbx.cif."""
  model = iotbx.cif.reader(file_path=_CIF_PATH).model()
  if model is None:
    raise RuntimeError(f"could not parse bvparm CIF: {_CIF_PATH}")
  blk = list(model.values())[0]
  table = {}
  for cation, valence, anion, ro, b in zip(
      blk["_valence_param_atom_1"], blk["_valence_param_atom_1_valence"],
      blk["_valence_param_atom_2"], blk["_valence_param_Ro"],
      blk["_valence_param_B"]):
    try:
      key = (metals.normalize_element(cation), int(valence),
             metals.normalize_element(anion))
      params = (float(ro), float(b))
    except ValueError:
      continue
    if key not in table:                        # first row = most reliable
      table[key] = params
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
  valences, n_used = {}, 0
  for q in ox_states:
    V, n = valence_sum(cation, q, contacts)
    if V is not None:
      valences[q] = V
      n_used = max(n_used, n)
  if not valences:
    return None
  best = min(valences, key=lambda q: abs(valences[q] - q))
  return group_args(best=best, valences=valences, n_used=n_used)
