"""Sanity-bounds QC for emitted metal-ligand distances (plan Sec 9.4).

Every emitted M-donor bond is cross-checked against two independent references
and *flagged, not fixed* (the QM target always stands):

- an empirical target window: the consensus of the Harding [54] and Bazayeva
  et al. [55] protein/CSD mean M-donor coordination distances for the metal and
  donor element, widened for xTB tolerance.  Outside it => the distance is
  atypical for that metal-donor pair.
- a bond-valence plausibility window from the vendored bond-valence parameters
  [35] (the same table qmi.bvs uses) over the metal's accessible oxidation
  states.  Outside it => the distance is physically implausible (a collapsed or
  dissociating contact) -- a stronger signal of a wrong state, bad cluster or
  SCF failure.

Coarse by design: the emitted value is a QM-optimised distance carrying xTB
error, so both windows are per (metal, donor element) with generous tolerance --
a plausibility check, not a crystallographic-accuracy check.
"""

from __future__ import absolute_import, division, print_function

from math import log

from libtbx import group_args

from mmtbx.geometry_restraints.qmi import bvs, electronic, metals


# Consensus mean M-donor coordination distances (Angstrom) as (mean_lo, mean_hi),
# per donor ELEMENT: the range spanned by the Harding (2006) Table 1 [54] PDB/CSD
# means and the Bazayeva et al. (2024) Table S2 [55] short (coordinating) peaks,
# over every donor kind of that element (water / carboxylate / carbonyl /
# hydroxyl O; His N; Cys/Met S).  The 2024 long secondary peaks (far O of
# monodentate carboxylate, second shell) are excluded.
TARGET_MEAN_RANGE_A = {
  # Zn  O 1.99-2.09 [54]; N His 2.01-2.05 [54,55]; S Cys 2.28-2.34 [54,55]
  ("Zn", "O"): (1.99, 2.09), ("Zn", "N"): (2.01, 2.05), ("Zn", "S"): (2.28, 2.34),
  # Mg  O 2.00-2.30 [54,55]; N His 2.10-2.20 [55] (absent in [54])
  ("Mg", "O"): (2.00, 2.30), ("Mg", "N"): (2.10, 2.20),
  # Mn  O 2.10-2.22 [54,55]; N His 2.16-2.25 [54,55]; S Cys 2.35 [54]
  ("Mn", "O"): (2.10, 2.22), ("Mn", "N"): (2.16, 2.25), ("Mn", "S"): (2.35, 2.35),
  # Fe  O 1.99-2.17 [54,55]; N His 2.00-2.20 [54,55]; S Cys/Met 2.28-2.32 [54,55]
  ("Fe", "O"): (1.99, 2.17), ("Fe", "N"): (2.00, 2.20), ("Fe", "S"): (2.28, 2.32),
  # Cu  O 1.96-2.13 [54]; N His 2.00-2.03 [54,55]; S Cys 2.15-2.2 / Met 2.5 [54,55]
  ("Cu", "O"): (1.96, 2.13), ("Cu", "N"): (2.00, 2.03), ("Cu", "S"): (2.15, 2.50),
}

# Widen the observed mean-range for the xTB M-L error plus the natural spread of
# the reference distributions.  An engineering tolerance, not a database value.
TARGET_PAD_A = 0.25

# Bond-valence plausibility band: an individual coordination bond carries a
# valence s = exp((R0 - R)/B) [27] between s_max (a collapsed, near-covalent
# contact) and s_min (a dissociating one).  Inverting, R = R0 - B*ln(s).  Wide by
# design, so only a gross error trips it.
BVS_VALENCE_BAND = (0.05, 1.5)   # (s_min, s_max)


def target_window(metal, donor_element):
  """(min, max) plausible distance (Angstrom) from the empirical target range
  for a (metal, donor element), or None when no target is tabulated."""
  r = TARGET_MEAN_RANGE_A.get((metals.normalize_element(metal),
                               metals.normalize_element(donor_element)))
  if r is None:
    return None
  return (r[0] - TARGET_PAD_A, r[1] + TARGET_PAD_A)


def bvs_window(metal, donor_element):
  """(min, max) physically-plausible distance (Angstrom) from the bond-valence
  params [35] over the metal's accessible oxidation states, or None when no
  parameter matches the (metal, donor element)."""
  s_min, s_max = BVS_VALENCE_BAND
  ox_states = electronic.accessible_oxidation_states(metal)
  if not ox_states:                    # None (no prior) or empty -> no window
    return None
  lo = hi = None
  for q in ox_states:
    params = bvs.valence_params(metal, q, donor_element)
    if params is None:
      continue
    r0, b = params
    r_short = r0 - b * log(s_max)   # highest valence -> shortest bond
    r_long = r0 - b * log(s_min)    # lowest valence -> longest bond
    lo = r_short if lo is None else min(lo, r_short)
    hi = r_long if hi is None else max(hi, r_long)
  if lo is None:
    return None
  return (lo, hi)


def check_bond(metal, donor_element, distance):
  """Flags for one emitted M-donor distance (Angstrom).

  Returns a list with any of ``'atypical'`` (outside the empirical target
  window) and ``'implausible'`` (outside the bond-valence window); empty when
  the distance passes every window for which a reference exists.
  """
  flags = []
  tw = target_window(metal, donor_element)
  if tw is not None and not (tw[0] <= distance <= tw[1]):
    flags.append("atypical")
  bw = bvs_window(metal, donor_element)
  if bw is not None and not (bw[0] <= distance <= bw[1]):
    flags.append("implausible")
  return flags


def check_site(metal, bonds):
  """Sanity-check every emitted bond of a site.

  *bonds* is an iterable of ``(donor_element, distance)`` pairs.  Returns
  ``group_args(n_atypical, n_implausible, details)`` where *details* is a list of
  short strings (one per flagged bond) for the log and provenance header.
  """
  n_atypical = n_implausible = 0
  details = []
  m = metals.normalize_element(metal)
  for donor_element, distance in bonds:
    flags = check_bond(metal, donor_element, distance)
    if not flags:
      continue
    if "atypical" in flags:
      n_atypical += 1
    if "implausible" in flags:
      n_implausible += 1
    details.append(f"{m}-{metals.normalize_element(donor_element)} "
                   f"{distance:.2f} A ({'/'.join(flags)})")
  return group_args(n_atypical=n_atypical, n_implausible=n_implausible,
                    details=details)
