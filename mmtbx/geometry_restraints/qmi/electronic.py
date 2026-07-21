"""Educated-guess electronic-state engine (plan 6): oxidation/spin to a
(d-count, unpaired-electron, multiplicity) state, from ligand-field rules.

This module is the spin-rule core (plan 6.2). Its numbers are not empirical
data: the unpaired-electron counts are *derived* by Hund's-rule filling of the
d-orbitals split by the coordination geometry (ligand-field theory, references
[12]), so per the companion spec's provenance carve-out ("ligand-field spin
counts ... appear concretely") they live in code, not a checksummed data file.
The d-count is a definition (group number minus oxidation state) and the group
numbers are the periodic table.

This module also carries the COMMITTING predictor (plan Sec 6). Given a metal,
its CShM geometry class, coordination number, and a donor summary (from
:mod:`mmtbx.geometry_restraints.qmi.donors`, all refinement-robust features),
:func:`predict_state` ALWAYS commits a single (oxidation, multiplicity) with a
confidence label -- it never defers to the user and never widens a sigma;
uncertainty is reported as the confidence, a different object than a restraint
sigma. Oxidation is narrowed by donor hardness + CN + geometry (plan Sec 6.3),
never by bond lengths. Spin is high-spin by default and flips to low spin only
on a determinable strong-field signature -- in biological macromolecules that is
a heme with a strong axial-donor pair (bis-His / His-Met) [45, 46, 48], not a
pi-acid diatomic (essentially exogenous here); non-heme biological Fe is
high-spin at every geometry [41].

Note: the ``oxidation_state`` here is the metal's *formal* charge. The QM
cluster's total charge is not this value: the builder assembles it from the
metal oxidation state plus the ligand formal charges (via
``quantum_interface.electrons``).
"""

from __future__ import absolute_import, division, print_function

from libtbx import group_args

from mmtbx.geometry_restraints.qmi import metals
from mmtbx.geometry_restraints.qmi.util import (
  LINEAR, TRIGONAL, TETRAHEDRAL, SQUARE_PLANAR, TRIGONAL_BIPYRAMIDAL,
  SQUARE_PYRAMIDAL, OCTAHEDRAL)


# Periodic-table group numbers (1-18) for the metals in qmi scope. Definitional
# (not empirical): for a transition metal the group number is the count of
# (n-1)d + ns valence electrons, so d-count = group - oxidation. Main-group
# metals (Na, Mg, K, Ca) have no d electrons in these states; the same
# arithmetic yields d-count 0 (closed shell), which is correct.
GROUP = {
  "Na": 1, "K": 1,
  "Mg": 2, "Ca": 2,
  "Sc": 3, "Ti": 4, "V": 5, "Cr": 6, "Mn": 7, "Fe": 8, "Co": 9, "Ni": 10,
  "Cu": 11, "Zn": 12,
}


def d_electron_count(element, oxidation):
  """d-electron count = group number - oxidation state (a definition).

  Raises KeyError for a metal outside :data:`GROUP` and ValueError if the
  result is outside 0..10.
  """
  el = metals.normalize_element(element)
  d = GROUP[el] - oxidation
  if not 0 <= d <= 10:
    raise ValueError(f"{el}{oxidation:+d} gives d-count {d}, outside 0..10")
  return d


# Accessible oxidation states and the maximum-likelihood default, per metal, for
# MONONUCLEAR biological sites (plan 6.1/6.3/6.6). Empirical reference data,
# cited [36]: the accessible states are standard inorganic chemistry (Zn is
# redox-inert d10, only +2; Mn spans +2/+3 in mononuclear sites, with +4 only in
# high-valent/polynuclear contexts such as the OEC, out of Phase-1 scope; Mg is
# main-group +2). ``default`` is the dominant biological resting state, used
# ONLY when the user supplies no oxidation state (plan 6.6); a user metal_state
# override is authoritative and bypasses these priors.
#
OXIDATION_PRIORS = {
  "Zn": group_args(accessible=(2,), default=2),
  "Mn": group_args(accessible=(2, 3), default=2),
  "Mg": group_args(accessible=(2,), default=2),
  "Fe": group_args(accessible=(2, 3), default=2),
  "Cu": group_args(accessible=(1, 2), default=2),
}
# ``default`` is the fallback oxidation the predictor commits when no robust
# donor/geometry signature resolves the state (plan Sec 6.6). Every metal here
# gets one: qmi always commits (it never defers). Mn(II) [43] and Cu(II) [39, 40]
# dominate their surveys; Fe defaults to +2 (ferrous), the lean of a photoreduced
# crystal [42], but predict_oxidation() overrides it toward Fe(III) on a ferric
# donor signature (tyrosinate [50], hard-anionic-O charge <= -2 [37, 44]) and
# back toward Fe(II) on a ferrous one (2-His-1-carboxylate facial triad [41], an
# all-thiolate Fe(S-Cys)n set [49]) -- with the confidence set accordingly.
#
# SPIN INVARIANT (load-bearing): the multiplicity commit is HIGH-SPIN by default
# and flips to low spin ONLY in predict_spin's heme branch (octahedral d4-d7 with
# a strong axial pair [45]). Non-heme biological Fe/Mn are high-spin at every
# geometry [41, 47]. Adding a metal with an accessible biological low-spin state
# OUTSIDE heme (e.g. Co, Ni) requires extending predict_spin's strong-field
# detection first, or the multiplicity will be silently wrong.


def accessible_oxidation_states(element):
  """Accessible oxidation states for a metal in mononuclear biological sites
  (cited [36]) as a tuple, or None if the metal is outside the vetted priors."""
  p = OXIDATION_PRIORS.get(metals.normalize_element(element))
  return p.accessible if p is not None else None


def default_oxidation(element):
  """Dominant biological resting-state oxidation to commit for a metal when
  structure leaves more than one candidate (plan 6.6), or None when the metal
  is outside the vetted priors OR has no dominant state (Fe)."""
  p = OXIDATION_PRIORS.get(metals.normalize_element(element))
  return p.default if p is not None else None


# Confidence labels for a committed prediction. Reported to the user (plan Sec 6,
# 9.3), NEVER folded into a restraint sigma.
CONF_HIGH, CONF_MEDIUM, CONF_LOW = "high", "medium", "low"
_CONF_RANK = {CONF_HIGH: 2, CONF_MEDIUM: 1, CONF_LOW: 0}


def min_confidence(*labels):
  """The weakest of several confidence labels (a commit is only as confident as
  its least-confident part)."""
  return min(labels, key=lambda c: _CONF_RANK[c])


def predict_oxidation(element, geometry, cn=None, summary=None):
  """Commit ONE oxidation state (+ confidence + reason) from refinement-robust
  features (plan Sec 6.3): donor identity/charge (via *summary* from
  :func:`mmtbx.geometry_restraints.qmi.donors.summarize`), coordination number,
  and CShM geometry class. NEVER bond lengths, energy, or diffraction.

  The rule per branch and its literature basis:

  - Single accessible state (Zn, Mg): committed, high confidence.
  - Cu [39, 40]: low CN (2-3) / linear-trigonal -> Cu(I); CN>=5 / square-planar/
    pyramidal/bipyramidal/octahedral -> Cu(II); CN 4 / tetrahedral is ambiguous
    -> Cu(II) default, low confidence. (Jahn-Teller 4+2 is a bond-length feature
    and is excluded.)
  - Fe: donor hardness commits ferrous vs ferric, Fe(III) being the harder acid
    [37, 44]; a coordinating tyrosinate is a ferric diagnostic [50]; the 2-His-
    1-carboxylate facial triad is the ferrous oxygenase motif [41]; an all-
    thiolate Fe(S-Cys)n set stays ferrous on a soft-S guard [49]; otherwise the
    ferrous photoreduction lean [42].
  - Mn: no robust II/III discriminator survives (Jahn-Teller excluded), so the
    dominant Mn(II) default [43], low confidence.

  Returns group_args(oxidation, confidence, reason), or None if the metal is
  outside the vetted priors (:data:`OXIDATION_PRIORS`).
  """
  el = metals.normalize_element(element)
  accessible = accessible_oxidation_states(el)
  if accessible is None:
    return None
  if len(accessible) == 1:
    return group_args(oxidation=accessible[0], confidence=CONF_HIGH,
                      reason="only one accessible oxidation state")

  if el == "Cu":                                                        # [39, 40]
    if (cn is not None and cn <= 3) or geometry in (LINEAR, TRIGONAL):
      return group_args(oxidation=1, confidence=CONF_HIGH,
                        reason="low CN / linear-trigonal geometry")
    if (cn is not None and cn >= 5) or geometry in (
        SQUARE_PLANAR, SQUARE_PYRAMIDAL, TRIGONAL_BIPYRAMIDAL, OCTAHEDRAL):
      return group_args(oxidation=2, confidence=CONF_HIGH,
                        reason="CN>=5 / square-planar-pyramidal-octahedral")
    return group_args(oxidation=2, confidence=CONF_LOW,
                      reason="CN 4 / tetrahedral is ambiguous (Cu(II) default)")

  if el == "Fe":
    if summary is not None:
      if summary.n_tyrosinate >= 1:                                     # [50]
        return group_args(oxidation=3, confidence=CONF_HIGH,
                          reason="coordinating tyrosinate")
      if (summary.n_his >= 2 and summary.n_carboxylate == 1             # [41]
          and summary.n_tyrosinate == 0):
        return group_args(oxidation=2, confidence=CONF_HIGH,
                          reason="2-His-1-carboxylate facial triad")
      if summary.n_cys_s >= 2 and summary.q_hard == 0:                  # [49]
        return group_args(oxidation=2, confidence=CONF_LOW,
                          reason="all-thiolate soft-S donor set")
      if summary.q_hard <= -2:                                         # [37, 44]
        return group_args(oxidation=3, confidence=CONF_MEDIUM,
                          reason="hard anionic O charge <= -2")
    return group_args(oxidation=2, confidence=CONF_LOW,                 # [42]
                      reason="no ferric signature (ferrous lean)")

  # Mn: dominant Mn(II) default, no robust II/III discriminator.            [43]
  default = default_oxidation(el)
  return group_args(oxidation=default if default is not None else accessible[0],
                    confidence=CONF_LOW,
                    reason="no robust discriminator (dominant default)")


def _aufbau_unpaired(n, levels):
  """Unpaired electrons when *n* electrons fill orbital energy *levels* (a tuple
  of orbital counts per level, low energy first), Hund's rule within each level.

  A level of ``size`` orbitals holding ``k`` electrons has ``k`` unpaired for
  ``k <= size`` (one each) and ``2*size - k`` above that (the rest paired).
  Weak field (high spin) is one level of all five d-orbitals, ``levels=(5,)``;
  strong-field octahedral fills t2g then eg, ``levels=(3, 2)``.
  """
  unpaired = 0
  for size in levels:
    k = min(n, 2 * size)
    unpaired += k if k <= size else 2 * size - k
    n -= k
    if n == 0:
      break
  return unpaired


# Orbital-splitting levels (low energy first) that are strong-field / low-spin
# distinct from the weak-field (5,) filling. Only octahedral differs.
_LOW_SPIN_LEVELS = {OCTAHEDRAL: (3, 2)}   # t2g, eg


def is_encoded(geometry, d_count):
  """True when the committed spin for (geometry, d_count) rests on an explicit
  rule rather than the weak-field fallback.

  Octahedral and tetrahedral are encoded for all d-counts. Square-planar is
  encoded for the strong-field d8 (diamagnetic) and the d9 doublet (Cu(II)). For
  any other geometry the spin is encoded only where the d-count is field-
  INDEPENDENT (0-3, 9, 10 give the same unpaired count in any field); the
  remaining field-dependent d-counts (4-7) at 5-coordinate / trigonal-prism
  geometries fall back to the weak-field high-spin count -- predict_spin still
  commits them high-spin, but flags the lower certainty via ``encoded``."""
  if geometry in (OCTAHEDRAL, TETRAHEDRAL):
    return True
  if geometry == SQUARE_PLANAR:
    return d_count in (8, 9)
  return d_count in (0, 1, 2, 3, 9, 10)


def is_field_dependent(geometry, d_count):
  """True where high-spin vs low-spin is a real choice (octahedral d4-d7);
  elsewhere the spin is fixed by geometry + d-count."""
  return geometry == OCTAHEDRAL and 4 <= d_count <= 7


def unpaired_electrons(geometry, d_count, low_spin=False):
  """Unpaired d electrons for a coordination geometry and d-count.

  Tetrahedral is always high spin (small Delta_t). Square-planar d8 is
  diamagnetic (strong-field). Octahedral is high spin unless ``low_spin`` and
  the d-count is field-dependent. Un-encoded (geometry, d_count) fall back to
  the weak-field high-spin count (the caller should flag; see
  :func:`is_encoded`).
  """
  if not 0 <= d_count <= 10:
    raise ValueError(f"d-count {d_count} outside 0..10")
  if geometry == SQUARE_PLANAR and d_count == 8:
    return 0
  if low_spin and geometry in _LOW_SPIN_LEVELS:
    return _aufbau_unpaired(d_count, _LOW_SPIN_LEVELS[geometry])
  return _aufbau_unpaired(d_count, (5,))   # weak-field / high-spin default


def multiplicity(unpaired):
  """Spin multiplicity 2S+1 = unpaired electrons + 1."""
  return unpaired + 1


def spin_state(element, oxidation, geometry, low_spin=False):
  """Convenience: (element, oxidation, geometry) -> the spin state as
  group_args(d_count, unpaired, multiplicity, low_spin, field_dependent,
  encoded). Ties d_electron_count + unpaired_electrons together for the caller.
  """
  d = d_electron_count(element, oxidation)
  u = unpaired_electrons(geometry, d, low_spin=low_spin)
  return group_args(
    d_count=d, unpaired=u, multiplicity=u + 1, low_spin=low_spin,
    field_dependent=is_field_dependent(geometry, d),
    encoded=is_encoded(geometry, d))


def predict_spin(element, oxidation, geometry, summary=None):
  """Commit low_spin (+ confidence + reason) for (element, oxidation, geometry)
  from robust features (plan Sec 6.2/6.6). High-spin by default; low spin only in
  the heme branch. The rule per branch and its basis:

  - Not field-dependent (all but octahedral d4-d7): spin fixed by geometry +
    d-count, high confidence.
  - Heme Fe (FeN4 topology in *summary*), octahedral d4-d7:
      6-coordinate with a strong axial pair (His/His or His/Met) -> low spin [45];
      6-coordinate with a weak axial (aquo/hydroxo met) -> high spin [48];
      otherwise (single axial / 5-coordinate) -> high spin [46].
  - Non-heme biological field-dependent site: weak-field donor set -> high spin,
    Fe [41], Mn [47, 12].

  Returns group_args(low_spin, confidence, reason).
  """
  d = d_electron_count(element, oxidation)      # may raise for a bad override
  if not is_field_dependent(geometry, d):
    return group_args(low_spin=False, confidence=CONF_HIGH,
                      reason="spin fixed by geometry and d-count")
  el = metals.normalize_element(element)
  if el == "Fe" and summary is not None and summary.heme:
    if summary.n_axial >= 2 and summary.axial_all_strong:              # [45]
      return group_args(low_spin=True, confidence=CONF_HIGH,
                        reason="6-coordinate heme, strong axial pair "
                               "(His/His or His/Met)")
    if summary.n_axial >= 2:                                           # [48]
      return group_args(low_spin=False, confidence=CONF_HIGH,
                        reason="6-coordinate heme, weak axial (aquo/hydroxo)")
    return group_args(low_spin=False, confidence=CONF_HIGH,            # [46]
                      reason="5-coordinate heme (single axial)")
  # Non-heme biological weak field -> high spin.                       [41, 47, 12]
  return group_args(low_spin=False, confidence=CONF_HIGH,
                    reason="non-heme biological weak-field donor set")


def predict_state(element, geometry, cn=None, summary=None):
  """Committing educated-guess electronic state (plan Sec 6): ALWAYS commits a
  single (oxidation_state, multiplicity) with a confidence label, from
  refinement-robust features only (donor *summary* + *cn* + CShM *geometry*
  class), never bond lengths / energy / diffraction.

  *summary* is a :func:`mmtbx.geometry_restraints.qmi.donors.summarize`
  group_args (donor hardness + heme topology), or None (oxidation/spin then rest
  on geometry/CN priors alone).

  Returns group_args(oxidation_state, multiplicity, low_spin, d_count,
  confidence, oxidation_confidence, spin_confidence, oxidation_reason,
  spin_reason, encoded), or None if the metal is outside the vetted priors and
  so cannot be predicted (the caller then requires a user override).
  """
  ox = predict_oxidation(element, geometry, cn, summary)
  if ox is None:
    return None
  try:
    d = d_electron_count(element, ox.oxidation)
    sp = predict_spin(element, ox.oxidation, geometry, summary)
    unpaired = unpaired_electrons(geometry, d, low_spin=sp.low_spin)
  except (KeyError, ValueError):
    return None
  return group_args(
    oxidation_state=ox.oxidation, multiplicity=unpaired + 1,
    low_spin=sp.low_spin, d_count=d,
    confidence=min_confidence(ox.confidence, sp.confidence),
    oxidation_confidence=ox.confidence, spin_confidence=sp.confidence,
    oxidation_reason=ox.reason, spin_reason=sp.reason,
    encoded=is_encoded(geometry, d))
