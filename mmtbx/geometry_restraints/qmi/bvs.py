"""Bond-valence sum from the vendored bvparm2020.cif [35], as a CONFIRM-ONLY,
resolution-gated oxidation cross-check (plan 6.4).

BVS is deliberately NOT used to deduce the oxidation state: bond lengths are too
resolution-sensitive for that, and (for redox metals) the crystallographic state
is often X-ray-photoreduced away from the biological one [10]. The primary
determination is the ligand-field/geometry guess in electronic.py. Here BVS only
CONFIRMS a supplied oxidation state, and only at sufficiently high resolution:
given the metal, a candidate oxidation state, and its first-shell donor
distances, it computes V = sum_i exp((R0 - R_i)/B) [27] and flags agreement
(|V - ox| <= tol) or disagreement. Below the resolution gate it returns None.

The CIF parse and valence sum are adapted from the calibration BVS in
~/qmi/calibration/prepare_qm_calcs.py (parses bvparm2020.cif with iotbx.cif,
keyed on (cation, oxidation, donor element), first/most-reliable row per combo).
That helper *deduced* the best self-consistent state; this module keeps the same
parse + sum but exposes only a confirm-only wrapper. Caveat noted there:
Fe(II)/Fe(III) R0 are close and S-donor rows sparse, so BVS barely separates
oxidation states at thiolate sites.
"""

from __future__ import absolute_import, division, print_function

import os
from math import exp
from functools import lru_cache

import iotbx.cif
from libtbx import group_args


_CIF_PATH = os.path.join(os.path.dirname(__file__), "bvparm2020.cif")

# BVS agreement tolerance in valence units (CheckMyMetal acceptable band [10]).
DEFAULT_TOLERANCE = 0.5
# Bond lengths are only trustworthy for BVS at high resolution. Conservative
# default policy threshold (Angstrom), NOT a measured constant; the plan calls
# for "very high resolution", so a caller may tighten it.
DEFAULT_RESOLUTION_GATE_A = 1.5


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
      key = (cation.capitalize(), int(valence), anion.capitalize())
    except ValueError:
      continue
    if key not in table:                        # first row = most reliable
      table[key] = (float(ro), float(b))
  return table


def valence_params(cation, oxidation, donor):
  """(R0, B) for a (cation, oxidation, donor-element) bond, or None if the
  vendored table has no entry."""
  return _table().get((cation.strip().capitalize(), oxidation,
                       donor.strip().capitalize()))


def valence_sum(cation, oxidation, contacts):
  """Brown-Altermatt bond-valence sum V = sum_i exp((R0 - R_i)/B) [27] over the
  first-shell *contacts* (``(donor_element, distance_A)`` tuples), under the
  given oxidation state. Returns ``(V, n_used)`` where n_used is the number of
  donors that had a parameter; ``V`` is None if none matched (a partial sum
  would be meaningless)."""
  sv = [exp((p[0] - dist) / p[1])
        for donor, dist in contacts
        for p in [valence_params(cation, oxidation, donor)]
        if p is not None]
  if not sv:
    return None, 0
  return sum(sv), len(sv)


def confirm_oxidation(cation, oxidation, contacts, resolution,
                      tolerance=DEFAULT_TOLERANCE,
                      resolution_gate_A=DEFAULT_RESOLUTION_GATE_A):
  """CONFIRM (never determine) a supplied oxidation state via BVS (plan 6.4).

  Returns None when *resolution* is unknown or coarser than the gate (BVS not
  trustworthy), or when the shell has no vendored parameters. Otherwise returns
  group_args(valence_sum, n_used, nearest_int, agrees, tolerance, resolution),
  with ``agrees = |V - oxidation| <= tolerance``. The caller FLAGS a
  disagreement; it never changes the oxidation state on this basis. Note BVS
  confirms the *crystallographic* state, which for a redox metal may be
  photoreduced away from the biological one [10].
  """
  if resolution is None or resolution > resolution_gate_A:
    return None
  V, n_used = valence_sum(cation, oxidation, contacts)
  if V is None:
    return None
  return group_args(valence_sum=V, n_used=n_used, nearest_int=int(round(V)),
                    agrees=abs(V - oxidation) <= tolerance,
                    tolerance=tolerance, resolution=resolution)
