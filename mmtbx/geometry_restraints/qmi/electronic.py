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
confidence label: it never defers to the user and never widens a sigma;
uncertainty is reported as the confidence, a different object than a restraint
sigma. Oxidation is narrowed by donor hardness + CN + geometry (plan Sec 6.3);
the only bond-length input is the carboxylate binding mode that sets the CN, not
a fit of the target distances. Spin is high-spin by default and flips to low
spin only on a determinable strong-field signature; in biological macromolecules
that is a heme with a strong axial-donor pair (bis-His / His-Met) [45, 46, 48],
not a pi-acid diatomic (essentially exogenous here); non-heme biological Fe is
high-spin at every geometry [41].

Note: the ``oxidation_state`` here is the metal's *formal* charge. The QM
cluster's total charge is not this value: the builder assembles it from the
metal oxidation state plus the ligand formal charges (via
``quantum_interface.electrons``).
"""

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
# gets one: qmi always commits (it never defers). Mn(II) [43] and Cu(II) [39,
# 40] dominate their surveys; Fe defaults to +2 (ferrous), the lean of a
# photoreduced crystal [42], but predict_oxidation() overrides it toward Fe(III)
# on a ferric donor signature (tyrosinate [50], hard-anionic-O charge <= -2 [37,
# 44]) and back toward Fe(II) on a ferrous one (2-His-1-carboxylate facial triad
# [41], an rubredoxin-type Fe(S-Cys)4 set, now read as FERRIC [76, 77, 78]),
# with the confidence set accordingly.
#
# SPIN INVARIANT (load-bearing): the multiplicity commit is HIGH-SPIN by default
# and flips to low spin ONLY in predict_spin's heme branch (octahedral d4-d7
# with a strong axial pair [45]). Non-heme biological Fe/Mn are high-spin at
# every geometry (Fe [41]; Mn(II)/Mn(III) [47, 61, 62], weak-field O/N per LFT
# [12]). Adding a metal with an accessible biological low-spin state OUTSIDE
# heme (e.g. Co, Ni) requires extending predict_spin's strong-field detection
# first, or the multiplicity will be silently wrong.


def accessible_oxidation_states(element):
  """Accessible oxidation states for a metal in mononuclear biological sites
  (cited [36]) as a tuple, or None if the metal is outside the vetted priors."""
  p = OXIDATION_PRIORS.get(metals.normalize_element(element))
  return p.accessible if p is not None else None


def default_oxidation(element):
  """Dominant biological resting-state oxidation to commit for a metal when
  structure leaves more than one candidate (plan 6.6), or None when the metal
  is outside the vetted priors."""
  p = OXIDATION_PRIORS.get(metals.normalize_element(element))
  return p.default if p is not None else None


# Confidence labels for a committed prediction. Reported to the user (plan Sec
# 6, 9.3), NEVER folded into a restraint sigma.
CONF_HIGH, CONF_MEDIUM, CONF_LOW = "high", "medium", "low"
_CONF_RANK = {CONF_HIGH: 2, CONF_MEDIUM: 1, CONF_LOW: 0}


def min_confidence(*labels):
  """The weakest of several confidence labels (a commit is only as confident as
  its least-confident part)."""
  return min(labels, key=lambda c: _CONF_RANK[c])


def _heme_oxidation(summary):
  """Commit a heme Fe oxidation from the axial ligand alone (plan Sec 6.3, Phase
  5). The porphyrin N4 is redox-invariant, so the sixth/axial ligand is the
  diagnostic; heme oxidation topology never supports high confidence.

  Ferric on any diagnostic ferric axial ligand, else a LOW-confidence ferrous
  lean. Ordering is first-match; a phenolate coexisting with a thiolate/water
  still yields ferric, so anionic O leads.
  """
  if summary.q_hard <= -1:  # [52, 50, 44]
    # Anionic axial O (Tyr phenolate; catalase-type ferric). The carboxylate
    # limb of q_hard is HSAB-extrapolated: no characterised axial-carboxylate
    # heme, but the direction (anionic O favours the harder Fe(III)) holds.
    #
    # This limb fires on a SINGLE anionic group, unlike the non-heme -2
    # threshold, so a charge the bond graph supplied is load-bearing here on its
    # own and the reason must say so rather than asserting a phenolate.
    why = "heme, anionic axial O (phenolate/catalase-type ferric)"
    n_bg = getattr(summary, "n_anionic_unplaced", 0)
    if n_bg:
      why += f", {n_bg} anionic O group(s) read from the ligand bond graph"
    return group_args(oxidation=3, confidence=CONF_MEDIUM, reason=why)
  if summary.n_cys_s >= 1:                                             # [51]
    # Coordinating Cys thiolate (charge < 0; neutral Met thioether is not
    # counted in n_cys_s). Heme-c peripheral thioether Cys sit off the metal and
    # are dropped by the first-shell cutoff, so cytochrome c is not misread
    # here.
    return group_args(
      oxidation=3, confidence=CONF_MEDIUM,
      reason="thiolate-ligated heme (P450/CPO/NOS ferric resting)")
  if summary.n_water_o >= 1:                                           # [48]
    # Bound axial water: aquo/hydroxo met (ferric). Deoxy ferrous is
    # 5-coordinate, so a sixth aqua ligand implies met. Keyed strictly on
    # water_o (HOH), so a diatomic OXY/CMO (het O/C) never trips it; a
    # non-coordinating distal water inside the cutoff stays a modelling risk.
    return group_args(oxidation=3, confidence=CONF_MEDIUM,
                      reason="aquo/hydroxo met-heme (bound axial water)")
  # Single-His deoxy, bis-His / His-Met carriers: topology cannot decide, and
  # this branch systematically leans the wrong way for the 5-coordinate
  # proximal-His peroxidases (HRP/CcP/APX, ferric resting [53]) and oxidised
  # bis-His/His-Met carriers ([45]). Flagged LOW; the sole honest lean is
  # photoreduction.
  return group_args(oxidation=2, confidence=CONF_LOW,  # [42, 46, 53, 45]
                    reason="single-His / weak-field heme; oxidation "
                           "indeterminate (globin-ferrous vs "
                           "peroxidase-ferric), photoreduction ferrous lean")


def predict_oxidation(element, geometry, cn=None, summary=None):
  """Commit ONE oxidation state (+ confidence + reason) from refinement-robust
  features (plan Sec 6.3): donor identity/charge (via *summary* from
  :func:`mmtbx.geometry_restraints.qmi.donors.summarize`), coordination number,
  and CShM geometry class. NEVER bond lengths, energy, or diffraction.

  Caveat for future rules: *geometry* is the CShM class of the EFFECTIVE-vertex
  set (a symmetric bidentate carboxylate merged to one vertex), while *cn* is
  the COORDINATIVE count (that chelate counted as its two O). So for a chelating
  site the shape's vertex count is below *cn* by the number of symmetric
  chelates; a rule that combines the two must not assume the shape's vertex
  count equals *cn*. The current rules stay consistent (the Cu CN cutoffs are
  monotone; the Mn-SOD gate requires TBP, which pins the chelate count to zero
  at cn==5).

  The rule per branch and its literature basis:

  - Single accessible state (Zn, Mg): committed, high confidence.
  - Cu [39, 40]: low CN (2-3) / linear-trigonal -> Cu(I); CN>=5 / square-planar/
    pyramidal/bipyramidal/octahedral -> Cu(II); CN 4 / tetrahedral is ambiguous
    -> Cu(II) default, low confidence. (Jahn-Teller 4+2 is a bond-length feature
    and is excluded.)
  - Fe heme (FeN4): the AXIAL ligand commits ferrous vs ferric (Sec 6.3, Phase
    5), the porphyrin N4 being redox-invariant: an anionic axial O (phenolate,
    catalase-type [52, 50]), a Cys thiolate (P450/CPO/NOS resting [51]), or a
    bound axial water (aquo/hydroxo met [48]) commits Fe(III); a single-His /
    bis-His / His-Met heme is indeterminate and leans ferrous at LOW confidence
    ([42, 46], with the known 5-coordinate-His peroxidase miss [53, 45]).
  - Fe non-heme: donor hardness commits ferrous vs ferric, Fe(III) being the
    harder acid [37, 44]; a coordinating tyrosinate is a ferric diagnostic [50];
    the 2-His-1-carboxylate facial triad is the ferrous oxygenase motif [41]; an
    a rubredoxin-type Fe(S-Cys)4 set commits FERRIC at LOW confidence
    [76, 77, 78], the ferrous reservoir case and photoreduction being carried in
    the confidence rather than the state [84, 87]; otherwise
    the ferrous photoreduction lean [42].
  - Mn: dominant Mn(II) default [56, 57], overridden to Mn(III) by the SOD-fold
    signature (5-coordinate trigonal-bipyramidal 3-His/carboxylate/solvent,
    resting Mn(III)-hydroxo [58, 59, 62, 63, 64], MEDIUM) or a coordinating
    tyrosinate ([60], LOW). No hardness / CN / Jahn-Teller rule: none separates
    Mn(II) from Mn(III) refinement-robustly.

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

  if el == "Cu":  # [39, 40]
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
      if summary.heme:
        return _heme_oxidation(summary)
      if summary.n_tyrosinate >= 1:                                     # [50]
        return group_args(oxidation=3, confidence=CONF_HIGH,
                          reason="coordinating tyrosinate")
      if (summary.n_his >= 2 and summary.n_carboxylate == 1             # [41]
          and summary.n_tyrosinate == 0):
        return group_args(oxidation=2, confidence=CONF_HIGH,
                          reason="2-His-1-carboxylate facial triad")
      # Rubredoxin-type Fe(S-Cys)4: FERRIC as isolated. Rubredoxin is purified
      # and crystallised as the red oxidised protein, is auto-oxidisable in air
      # [76], is high-spin Fe(III) by Mossbauer [77], and its couple sits near
      # -50 to 0 mV, far below O2/H2O, so a protein handled without a reductant
      # is ferric [78]. This REPLACES an inverted guard that read the same donor
      # set as ferrous; its citation [49] documents a site crystallised oxidised
      # and reduced only by adding dithionite, and never supported the ferrous
      # direction.
      #
      # Deliberately narrow: exactly four donors, all Cys thiolate, nothing else
      # (cn == 4 with four Cys residues admits no fifth donor). Broader forms
      # mis-fire. Nitrile hydratase (3 Cys-derived S + backbone N, CN 5) is a
      # LOW-SPIN Fe(III) [86] that this must not touch, and a >=2-Cys form also
      # catches every [2Fe-2S] iron, which presents 2 Cys S + 2 inorganic S.
      #
      # LOW confidence, and the reason says why: the site is ferrous in
      # reservoir and carrier roles (e.g. the rubredoxin-like centre of
      # 3-hydroxyanthranilate-3,4-dioxygenase [87], deposited ferrous many times
      # over), and a crystal at routine dose is partly photoreduced regardless
      # [84], so no single structure's Fe-S distances settle it.
      if cn == 4 and summary.n_cys_s == 4:                     # [76, 77, 78]
        return group_args(oxidation=3, confidence=CONF_LOW,
                          reason="rubredoxin-type Fe(S-Cys)4, ferric as "
                                 "isolated; ferrous in reservoir roles and "
                                 "photoreduced at routine dose")
      if summary.q_hard <= -2:  # [37, 44]
        why = "hard anionic O charge <= -2"
        n_bg = getattr(summary, "n_anionic_unplaced", 0)
        if n_bg:
          why += f", {n_bg} anionic O group(s) read from the ligand bond graph"
        return group_args(oxidation=3, confidence=CONF_MEDIUM, reason=why)
    return group_args(oxidation=2, confidence=CONF_LOW,                 # [42]
                      reason="no ferric signature (ferrous lean)")

  if el == "Mn":
    if summary is not None:
      # SOD-fold resting/oxidised Mn(III): a 5-coordinate trigonal-bipyramidal
      # 3-His / 1-carboxylate / 1-solvent site, resting Mn(III)-hydroxo by QM
      # and neutron structure [58, 59, 62, 64]. MEDIUM only, the reduced
      # Mn(II)-SOD state is isostructural [57].
      if (cn == 5 and geometry == TRIGONAL_BIPYRAMIDAL
          and summary.n_his == 3 and summary.n_carboxylate == 1
          and summary.n_water_o >= 1 and summary.n_tyrosinate == 0
          and summary.n_cys_s == 0 and not summary.heme):
        return group_args(oxidation=3, confidence=CONF_MEDIUM,
                          reason="SOD-fold 5-coordinate TBP "
                          "(3-His/carboxylate/solvent)")
      # A coordinating Tyr phenolate leans Mn(III) (transferrin-type [60]); rare
      # and indirectly supported, so LOW (cf. the HIGH-confidence Fe rule [50]).
      if summary.n_tyrosinate >= 1:
        return group_args(oxidation=3, confidence=CONF_LOW,
                          reason="coordinating tyrosinate")
    # Dominant Mn(II) default; no robust II/III discriminator otherwise [56,
    # 57].
    return group_args(oxidation=2, confidence=CONF_LOW,
                      reason="no ferric signature (Mn(II) default)")

  # Any other multi-state metal: dominant-default fallback.
  default = default_oxidation(el)
  return group_args(oxidation=default if default is not None else accessible[0],
                    confidence=CONF_LOW,
                    reason="no robust discriminator (dominant default)")


def oxidation_uses_donors(element):
  """Whether :func:`predict_oxidation`'s branch for *element* reads the donor
  summary at all.

  Only the Fe and Mn branches dereference ``summary``; Zn and Mg short-circuit
  on a single accessible state, and Cu decides from CN and CShM class alone.
  Callers use this to tell "the donor evidence was empty" (which should lower
  confidence) from "the donor evidence was irrelevant" (which should not):
  capping a Zn whose oxidation is definitional would contradict its own reason
  string.

  MUST track predict_oxidation: if a branch starts or stops reading *summary*,
  update this. Covered by tst_qmi.
  """
  return metals.normalize_element(element) in ("Fe", "Mn")


def oxidation_uses_donor_charge(element):
  """Whether :func:`predict_oxidation`'s branch for *element* reads
  ``summary.q_hard``, as opposed to the motif counts.

  Only Fe does, in both its heme and non-heme limbs. Mn reads the summary but
  keys on motifs (n_his / n_carboxylate / n_tyrosinate) and never on q_hard, so
  an anionic charge recovered from the bond graph is not evidence for a Mn site
  and must not lift its confidence cap.

  MUST track predict_oxidation, as :func:`oxidation_uses_donors` must.
  """
  return metals.normalize_element(element) == "Fe"


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
  geometries fall back to the weak-field high-spin count; predict_spin still
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
  from robust features (plan Sec 6.2/6.6). High-spin by default; low spin only
  in the heme branch. The rule per branch and its basis:

  - Not field-dependent (all but octahedral d4-d7): spin fixed by geometry +
    d-count, high confidence.
  - Heme Fe (FeN4 topology in *summary*), octahedral d4-d7:
      6-coordinate with a strong axial pair (His/His or His/Met) -> low spin
      [45]; 6-coordinate with a weak axial (aquo/hydroxo met) -> high spin [48];
      otherwise (single axial / 5-coordinate) -> high spin [46].
  - Non-heme biological field-dependent site: weak-field donor set -> high spin,
    Fe [41], Mn [47, 61, 62, 12].

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
  # Non-heme biological weak field -> high spin. [41, 47, 61, 62, 12]
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


def enumerate_states(element, geometry=None):
  """Accessible (oxidation, multiplicity) states for a multi-state scan.

  Each accessible oxidation is crossed with the ligand-field spin extremes:
  weak-field high-spin (the ``(5,)`` filling) and strong-field low-spin (the
  ``(3,2)`` t2g/eg filling). The two extremes coincide (a single state) for
  d0-d3, d8-d10 and closed shells.

  *geometry* selects how much of that space to emit:

  * ``None`` (the ``states=all`` scan): GEOMETRY-BLIND, both extremes for every
    oxidation, including low-spin at geometries the predictor would commit high-
    spin (e.g. tetrahedral). Brute-force everything, trust nothing.
  * a CShM class (the ``states=reasonable`` scan): only the spin(s) the ligand
    field permits at that geometry (tetrahedral -> high-spin only, square-planar
    d8 -> diamagnetic, octahedral d4-d7 -> both). An unclassified geometry
    cannot prune, so pass ``None`` and it degrades to the blind scan.

  This is a scan enumeration, NOT a committing prediction (predict_spin stays
  geometry-aware and commits a single state).

  Returns a list of ``(oxidation, multiplicity)``, or None if the metal is
  outside the vetted priors (:data:`OXIDATION_PRIORS`)."""
  accessible = accessible_oxidation_states(element)
  if accessible is None:
    return None
  out = []
  for ox in accessible:
    try:
      d = d_electron_count(element, ox)
    except (KeyError, ValueError):
      continue
    if geometry is None:                        # geometry-blind (states=all)
      hs = _aufbau_unpaired(d, (5,)) + 1  # weak field (one 5-orbital level)
      ls = _aufbau_unpaired(d, (3, 2)) + 1      # strong field (t2g then eg)
    else:  # ligand-field-allowed (reasonable)
      hs = unpaired_electrons(geometry, d, low_spin=False) + 1
      ls = unpaired_electrons(geometry, d, low_spin=True) + 1
    for mult in ([hs] if hs == ls else [hs, ls]):
      out.append((ox, mult))
  return out
