"""Educated-guess electronic-state engine (plan 6): oxidation/spin to a
(d-count, unpaired-electron, multiplicity) state, from ligand-field rules.

This module is the spin-rule core (plan 6.2). Its numbers are not empirical
data: the unpaired-electron counts are *derived* by Hund's-rule filling of the
d-orbitals split by the coordination geometry (ligand-field theory, references
[12]), so per the companion spec's provenance carve-out ("ligand-field spin
counts ... appear concretely") they live in code, not a checksummed data file.
The d-count is a definition (group number minus oxidation state) and the group
numbers are the periodic table.

What is NOT here, because it is genuinely empirical reference data (later
steps): which oxidation states a metal actually adopts, the geometry priors
that narrow oxidation (plan 6.3), the spectrochemical ordering, and the
strong-field signatures that flip octahedral d4-d7 to low spin (plan 6.6).

High-spin is the default (plan 6.6); the caller flips to low spin only on a
strong-field signature, and only where :func:`is_field_dependent` is True.

Note: the ``oxidation_state`` here is the metal's *formal* charge. The QM
cluster's total charge is not this value: the builder assembles it from the
metal oxidation state plus the ligand formal charges (via
``quantum_interface.electrons``).
"""

from __future__ import absolute_import, division, print_function

from libtbx import group_args


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
  el = element.strip().capitalize()
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
# INVARIANT (load-bearing): guess_state is unconditionally HIGH-SPIN. That is
# safe ONLY because every metal here is closed-shell (Zn d10, Mg d0) or
# biologically high-spin (Mn(II) d5, Mn(III) d4). Do NOT add a metal that has an
# accessible biological low-spin state (Fe, Co, Ni) without first implementing
# the low-spin flip (plan 6.6) + strong-field detection (spectrochemical series
# [13]), or the engine will silently emit the wrong multiplicity (e.g. low-spin
# d6 Fe(II) heme is S=0, not the S=2 the HS default would give).
OXIDATION_PRIORS = {
  "Zn": group_args(accessible=(2,), default=2),
  "Mn": group_args(accessible=(2, 3), default=2),
  "Mg": group_args(accessible=(2,), default=2),
  "Fe": group_args(accessible=(2, 3), default=None),
  "Cu": group_args(accessible=(1, 2), default=2),
}
# ``default`` is the DOMINANT biological resting oxidation state (a coarse
# survey-based prior, not a per-site likelihood) to commit when structure leaves
# more than one candidate. It is set where one state clearly dominates the
# surveys -- Mn(II) [43]; Cu(II), the common oxidised state [39, 40] (and Cu
# usually narrows by geometry anyway) -- and None for Fe, where Fe(II)/Fe(III)
# are both common and indistinguishable from geometry [41] and the crystal is
# often photoreduced [42]. So an ambiguous Fe site is NOT committed to a single
# guess: it awaits a user override or the plan Sec 7 run-both-and-inflate.
# candidate_oxidations() returns the real candidate set regardless.
# ``default`` is a weak max-likelihood fallback for the multi-state metals: Fe
# has no defensible structure-independent default [41] and photoreduced crystals
# lean ferrous [42] (so +2 is the fallback); Cu +2 is the common oxidised state.
# candidate_oxidations() carries the real candidate set, so the default is used
# only as the ordering/"if forced to pick one" seed, not to commit a state.


def accessible_oxidation_states(element):
  """Accessible oxidation states for a metal in mononuclear biological sites
  (cited [36]) as a tuple, or None if the metal is outside the vetted priors."""
  p = OXIDATION_PRIORS.get(element.strip().capitalize())
  return p.accessible if p is not None else None


def default_oxidation(element):
  """Dominant biological resting-state oxidation to commit for a metal when
  structure leaves more than one candidate (plan 6.6), or None when the metal
  is outside the vetted priors OR has no dominant state (Fe)."""
  p = OXIDATION_PRIORS.get(element.strip().capitalize())
  return p.default if p is not None else None


def candidate_oxidations(element, geometry, cn=None):
  """Plausible oxidation-state CANDIDATES, narrowed from the coordination
  geometry + number where structure allows (plan 6.3). NEVER uses bond lengths.

  Returns a list (subset of the accessible states), max-likelihood (default)
  first, or None if the metal has no priors. Structure narrows Cu well but not
  Fe/Mn, so those carry both:

  - Single accessible state (Zn, Mg): that state.
  - Cu [39, 40]: low CN (2-3) or linear/trigonal -> Cu(I); CN 5-6 or square-
    planar/pyramidal/octahedral -> Cu(II); CN 4 / tetrahedral is ambiguous ->
    both. (A type-1 blue-copper 2His/Cys/Met site is entatic and resists this;
    donor-aware low-confidence flagging is future scope.)
  - Fe [41, 42]: high-spin d6 and d5 share CN/geometry/donors and the
    crystallographic state is often photoreduced, so both {II, III} are carried.
  - Mn: carries both (a mild Mn(III) Jahn-Teller up-weight needs axial detection
    and is deferred).

  Consuming a multi-candidate result (run each, inflate the restraint) is the
  plan Sec 7 job; this function does not commit a single state.
  """
  el = element.strip().capitalize()
  accessible = accessible_oxidation_states(el)
  if accessible is None:
    return None
  if len(accessible) == 1:
    return list(accessible)

  if el == "Cu":
    if (cn is not None and cn <= 3) or geometry in ("linear", "trigonal"):
      cands = [1]
    elif (cn is not None and cn >= 5) or geometry in (
        "square_planar", "square_pyramidal", "octahedral"):
      cands = [2]
    else:
      cands = list(accessible)          # CN 4 / tetrahedral -> ambiguous
  else:
    cands = list(accessible)            # Fe, Mn -> carry both

  default = default_oxidation(el)       # order max-likelihood first
  if default in cands and len(cands) > 1:
    cands = [default] + [o for o in cands if o != default]
  return cands


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
_LOW_SPIN_LEVELS = {"octahedral": (3, 2)}   # t2g, eg


def is_encoded(geometry, d_count):
  """True when the spin rule for (geometry, d_count) is well defined here.

  Octahedral and tetrahedral are encoded for all d-counts; square-planar only
  for the strong-field d8 (diamagnetic). Other geometries (5-coordinate, etc.)
  are not encoded in Phase 1: the caller should flag and use the high-spin
  default (plan 6.6)."""
  if geometry in ("octahedral", "tetrahedral"):
    return True
  if geometry == "square_planar":
    return d_count == 8
  return False


def is_field_dependent(geometry, d_count):
  """True where high-spin vs low-spin is a real choice (octahedral d4-d7);
  elsewhere the spin is fixed by geometry + d-count."""
  return geometry == "octahedral" and 4 <= d_count <= 7


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
  if geometry == "square_planar" and d_count == 8:
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


def guess_state(element, geometry):
  """Educated-guess electronic state (plan 6) for a metal with no user override:
  the biological-default oxidation (6.6) and the high-spin ligand-field spin for
  the coordination *geometry* (a shape name, or None if unclassified).

  Returns group_args(oxidation_state, multiplicity, d_count, unpaired,
  field_dependent, encoded), or None if the metal is outside the vetted priors
  (:data:`OXIDATION_PRIORS`) and so cannot be guessed. Always high-spin here
  (plan 6.6 default); the low-spin flip for field-dependent octahedral d4-d7,
  and oxidation narrowing (6.3), come later. When ``encoded`` is False the shape
  has no spin rule and the high-spin default is a flagged fallback.
  """
  oxidation = default_oxidation(element)
  if oxidation is None:
    return None
  st = spin_state(element, oxidation, geometry)
  return group_args(
    oxidation_state=oxidation, multiplicity=st.multiplicity,
    d_count=st.d_count, unpaired=st.unpaired,
    field_dependent=st.field_dependent, encoded=st.encoded)
