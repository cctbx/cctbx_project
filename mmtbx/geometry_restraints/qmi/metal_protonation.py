"""Protonation rules for metal-coordinating donors.

A program that protonates a structure without looking at its metals gets the
metal sites wrong in a predictable way: it builds every side chain in the
neutral form it would take in bulk water.  A cysteine on a zinc is a thiolate, a
histidine on an iron donates a lone pair and must be H-free on the coordinating
ring nitrogen, and the N-terminus of an ATCUN copper site is a deprotonated
amidate.  These rules say which hydrogens a metal site does not want.

The rules are kept separate from the pipeline that applies them so that a
protonation program can consume them and place the right hydrogens in the FIRST
place, rather than a downstream consumer removing hydrogens that were just
added.  Nothing here imports the rest of qmi.  That is why deciding and doing
are separate:

- :func:`decide` is a pure function from facts about one metal-donor pair to an
  action.  It touches no hierarchy and does no I/O, so it can be called from
  inside a placement loop, before the hydrogen in question exists.
- :func:`propose` walks a hierarchy and returns one change record per donor it
  would touch, naming the hydrogen, the rule and the reason.  It mutates
  nothing.
- :func:`apply_changes` is the reference applier, for a caller that already has a
  protonated structure and wants it corrected.

WHAT THIS MODULE DELIBERATELY DOES NOT DECIDE.  Two policy questions belong to
the caller, and are required arguments rather than defaults:

- *which atoms are metals worth scoring against*.  A consumer that models only
  some metals wants donors near the others left alone; a general protonation
  program wants every metal.  Neither is a fact about chemistry.
- *how close counts as coordinating*.  The right cutoff depends on what the
  caller does with a false positive.  Supplying a default here would let a
  Ser hydroxyl near a passing potassium quietly lose its proton.

KNOWN LIMITS.  Each donor is scored against its NEAREST metal only, so in a
di-metal site the other metal cannot induce the chemistry; the caller's own
site list is the place to handle that.  Rules keyed on residue name cannot see
the same chemistry on a ligand (a het-group thiol on a zinc is missed).

References: HSAB and metal-induced thiolate/phenolate, Pearson 1963 [37] and
Dudev & Lim 2003 [38]; ATCUN N-terminal deprotonation [71].
"""

from libtbx import group_args
from scitbx import matrix


# Anionic donors: (resname, donor heavy-atom) -> the acid H names on that donor.
COORD_DEPROT_RULES = {
  ("CYS", "SG"):  frozenset(("HG", "DG")),
  ("TYR", "OH"):  frozenset(("HH", "DH")),
  ("ASP", "OD1"): frozenset(("HD2", "DD2")),
  ("ASP", "OD2"): frozenset(("HD2", "DD2")),
  ("GLU", "OE1"): frozenset(("HE2", "DE2")),
  ("GLU", "OE2"): frozenset(("HE2", "DE2")),
}

# Thiolate and phenolate formation is metal-INDUCED: only a soft or borderline,
# Lewis-acidic cation supplies the covalency that stabilises the anion.  A
# carboxylate is anionic by pH whatever the metal, so it is NOT gated.
SOFT_DONOR_KEYS = frozenset((("CYS", "SG"), ("TYR", "OH")))

# The d-block transition metals + group 12: acidic/soft enough to induce a
# thiolate or phenolate.  Excludes the hard alkali/alkaline-earth ions, which
# rarely bind S or phenol O at all and do not induce the anion [37, 38].  The
# p-block metals are left out as a conservative choice, not a chemical claim: a
# Pb or Sn thiolate is real and needs the caller to extend this set.
ANION_INDUCING_METALS = frozenset((
  "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
  "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
  "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
))

# Metal-induced N-terminal deprotonation at Cu(II)/Ni(II), the ATCUN motif [71].
ATCUN_METALS = frozenset(("Cu", "Ni"))
BACKBONE_AMIDE_H = frozenset(("H", "HN", "D", "DN"))
NTERM_AMINE_H = frozenset(("H1", "H2", "H3", "D1", "D2", "D3"))

# Histidine stays NEUTRAL, never an imidazolate: the coordinating ring N donates
# its lone pair and must be H-free, and the partner ring N keeps a proton.
HIS_RING_NEIGHBOURS = {"ND1": ("CG", "CE1"), "NE2": ("CD2", "CE1")}
HIS_RING_H = {"ND1": ("HD1", "DD1"), "NE2": ("HE2", "DE2")}
HIS_PARTNER = {"ND1": "NE2", "NE2": "ND1"}

NH_BOND_A = 1.01     # placed N-H length, relaxed by any subsequent optimisation

# Actions a rule can ask for.
KEEP = "keep"
REMOVE_H = "remove_h"
MOVE_TO_PARTNER = "move_to_partner"


def default_policy(**kwargs):
  """The rule set, as data a caller may override.

  These are choices, not laws: ``anion_inducing_metals`` is a conservative
  element list with no oxidation state in it (it cannot tell a Cu(I) thioether
  from a Cu(II) thiolate), and ``deprotonate_observed_h`` decides what to do
  when a hydrogen is a measurement rather than a placement.
  """
  policy = group_args(
    anion_inducing_metals=ANION_INDUCING_METALS,
    atcun_metals=ATCUN_METALS,
    # A neutron or sub-angstrom structure OBSERVES its hydrogens.  Removing one
    # because a name table says a coordinating Cys must be a thiolate discards
    # data; the default refuses, and a caller that knows better can say so.
    deprotonate_observed_h=False)
  for key, value in kwargs.items():
    setattr(policy, key, value)
  return policy


def donor_facts(resname, atom_name, element, metal_element, distance,
                h_names=(), partner_h_names=(), is_standard_aa=False,
                h_observed=False):
  """One metal-donor pair, reduced to what :func:`decide` needs.

  *h_names* are the hydrogens presently on the donor heavy atom (empty is
  meaningful: it says the donor is already deprotonated).  *partner_h_names*
  matters only for an imidazole ring N, where the partner N's protonation
  distinguishes a doubly-protonated HIP from a wrong tautomer.
  """
  return group_args(
    resname=str(resname).strip().upper(),
    atom_name=str(atom_name).strip().upper(),
    element=str(element).strip().capitalize(),
    metal_element=str(metal_element).strip().capitalize(),
    distance=distance,
    h_names=tuple(str(h).strip().upper() for h in h_names),
    partner_h_names=tuple(str(h).strip().upper() for h in partner_h_names),
    is_standard_aa=bool(is_standard_aa),
    h_observed=bool(h_observed))


def _outcome(action, rule, reason, h_name=None, h_select=None,
             new_name=None, mark_acceptor=False):
  return group_args(action=action, rule=rule, reason=reason, h_name=h_name,
                    h_select=h_select, new_name=new_name,
                    mark_acceptor=mark_acceptor)


def decide(facts, policy=None):
  """What, if anything, to do about the hydrogens on one coordinating donor.

  Pure: no hierarchy, no geometry beyond the distance already in *facts*, no
  I/O.  Returns ``group_args(action, rule, reason, h_name, h_select, new_name,
  mark_acceptor)``, where *action* is :data:`KEEP`, :data:`REMOVE_H` or
  :data:`MOVE_TO_PARTNER`.

  ``h_name`` names the target hydrogen when the rule can; ``h_select`` instead
  names a geometric criterion the caller must resolve (only
  ``'nearest_to_metal'``, for an -NH3+ whose three hydrogens are equivalent
  until you look at where the metal is).  ``mark_acceptor`` says the donor
  becomes a hydrogen-bond acceptor once stripped, which a placement program
  needs in order to re-optimise the surrounding network.
  """
  if policy is None:
    policy = default_policy()
  key = (facts.resname, facts.atom_name)

  if facts.h_observed and not policy.deprotonate_observed_h:
    return _outcome(KEEP, "observed-h",
                    "hydrogen is observed, not placed; not removing data")

  # Cys / Tyr / Asp / Glu.
  if key in COORD_DEPROT_RULES:
    if not facts.h_names:
      return _outcome(KEEP, "anionic-donor", "already deprotonated")
    acid_h = [h for h in facts.h_names if h in COORD_DEPROT_RULES[key]]
    if not acid_h:
      return _outcome(KEEP, "anionic-donor", "no acid hydrogen present")
    if (key in SOFT_DONOR_KEYS
        and facts.metal_element not in policy.anion_inducing_metals):
      # A hard, low-charge ion does not induce a thiolate or phenolate, and
      # mostly does not bind these donors at all [37, 38].
      return _outcome(KEEP, "hsab-gate",
                      f"{facts.metal_element} does not induce an anionic "
                      f"{facts.resname} donor")
    return _outcome(REMOVE_H, "anionic-donor",
                    f"{facts.resname} {facts.atom_name} coordinating "
                    f"{facts.metal_element} is anionic",
                    h_name=acid_h[0], mark_acceptor=True)

  # ATCUN, at a coordinating Cu(II)/Ni(II) N-terminus [71].
  if facts.atom_name == "N" and facts.is_standard_aa:
    if facts.metal_element not in policy.atcun_metals:
      return _outcome(KEEP, "atcun-gate",
                      f"{facts.metal_element} is not an ATCUN-type ion")
    amine_h = [h for h in facts.h_names if h in NTERM_AMINE_H]
    # All three present means a free alpha-amino terminus (-NH3+), which binds
    # as -NH2.  Requiring three is what keeps this off a backbone amide, which
    # carries one, and off an already-neutral amine.
    if len(amine_h) >= 3:
      return _outcome(REMOVE_H, "atcun-amine",
                      "N-terminal -NH3+ binds an ATCUN ion as -NH2",
                      h_select="nearest_to_metal", mark_acceptor=True)
    amide_h = [h for h in facts.h_names if h in BACKBONE_AMIDE_H]
    if amide_h:
      return _outcome(REMOVE_H, "atcun-amidate",
                      "backbone amide coordinating an ATCUN ion is an amidate",
                      h_name=amide_h[0], mark_acceptor=True)
    return _outcome(KEEP, "atcun", "no removable N-terminal hydrogen")

  # Histidine: keep the ring NEUTRAL.  Never strip the last ring proton, which
  # would make an imidazolate.
  if facts.resname == "HIS" and facts.atom_name in HIS_RING_H:
    ring_h = [h for h in facts.h_names if h in HIS_RING_H[facts.atom_name]]
    if not ring_h:
      return _outcome(KEEP, "his-tautomer",
                      "coordinating ring N already H-free (neutral tautomer)")
    partner = HIS_PARTNER[facts.atom_name]
    if any(h in HIS_RING_H[partner] for h in facts.partner_h_names):
      # Both ring N protonated: an imidazolium.  Dropping the coordinating one
      # leaves a neutral ring.
      return _outcome(REMOVE_H, "his-hip",
                      "doubly-protonated ring; removing the coordinating "
                      "proton gives a neutral His",
                      h_name=ring_h[0], mark_acceptor=True)
    # Only the coordinating N carries a proton: the wrong tautomer.  Move it
    # across the ring rather than removing it.  The name encodes the tautomer
    # site, so it is renamed with it.
    isotope = "D" if ring_h[0].startswith("D") else "H"
    return _outcome(MOVE_TO_PARTNER, "his-tautomer",
                    "wrong tautomer; moving the proton to the partner ring N",
                    h_name=ring_h[0],
                    new_name=isotope + ("E2" if partner == "NE2" else "D1"),
                    mark_acceptor=True)

  return _outcome(KEEP, "no-rule", "no metal-protonation rule for this donor")


def imidazole_nh_xyz(atoms_by_name, n_name, nh_len=NH_BOND_A):
  """In-plane bisector position for an N-H on imidazole ring nitrogen *n_name*,
  or None when a ring atom needed for the construction is missing."""
  neighbours = HIS_RING_NEIGHBOURS.get(n_name)
  if neighbours is None:
    return None
  n_atom = atoms_by_name.get(n_name)
  a1 = atoms_by_name.get(neighbours[0])
  a2 = atoms_by_name.get(neighbours[1])
  if n_atom is None or a1 is None or a2 is None:
    return None
  npos = matrix.col(n_atom.xyz)
  v1 = (matrix.col(a1.xyz) - npos).normalize()
  v2 = (matrix.col(a2.xyz) - npos).normalize()
  bisector = v1 + v2
  if bisector.length() < 1e-6:
    return None
  return tuple(npos - bisector.normalize() * nh_len)   # outward, in ring plane


def _change(action, h_atom, donor_atom, metal_atom, distance, rule, reason,
          new_name=None, new_xyz=None, mark_acceptor=False):
  ag = h_atom.parent()
  rg = ag.parent() if ag is not None else None
  ch = rg.parent() if rg is not None else None
  return group_args(
    action=action, h_atom=h_atom, donor_atom=donor_atom, metal_atom=metal_atom,
    distance=distance, rule=rule, reason=reason,
    new_name=new_name, new_xyz=new_xyz, mark_acceptor=mark_acceptor,
    # Identity that survives a re-sort or a reprocessing, unlike i_seq.
    h_id=(ch.id.strip() if ch is not None else "",
          rg.resseq.strip() if rg is not None else "",
          rg.icode.strip() if rg is not None else "",
          ag.altloc.strip() if ag is not None else "",
          ag.resname.strip().upper() if ag is not None else "",
          h_atom.name.strip().upper()))


def _nearest(atom, metal_atoms):
  xyz = matrix.col(atom.xyz)
  return min(((m, (xyz - matrix.col(m.xyz)).length()) for m in metal_atoms),
             key=lambda md: md[1])


def propose(hierarchy, metal_atoms, coord_cutoff, policy=None,
            is_standard_aa=None, normalize_element=None, h_observed=False):
  """The metal-protonation changes *hierarchy* needs.  Mutates nothing.

  *metal_atoms* are the metals to score donors against and *coord_cutoff* is a
  callable ``(metal_element, donor_element) -> float``; both are required, see
  the module docstring on why neither gets a default here.

  *is_standard_aa* is a callable ``resname -> bool`` used to recognise a peptide
  backbone N (the ATCUN rules); without it those rules do not fire.
  *normalize_element* maps an element string to its canonical symbol.
  """
  if policy is None:
    policy = default_policy()
  if normalize_element is None:
    normalize_element = lambda e: str(e).strip().capitalize()
  if is_standard_aa is None:
    is_standard_aa = lambda resname: False
  if not metal_atoms:
    return []

  changes, claimed = [], set()

  def _emit(donor_atom, ag, facts_kwargs, metal, distance):
    """Decide for one donor and turn the outcome into an change."""
    by_name = {a.name.strip().upper(): a for a in ag.atoms()}
    outcome = decide(donor_facts(**facts_kwargs), policy=policy)
    if outcome.action == KEEP:
      return
    if outcome.h_select == "nearest_to_metal":
      m_xyz = matrix.col(metal.xyz)
      cands = [a for a in ag.atoms()
               if a.name.strip().upper() in NTERM_AMINE_H
               and a.memory_id() not in claimed]
      if not cands:
        return
      h_atom = min(cands,
                   key=lambda h: (matrix.col(h.xyz) - m_xyz).length())
    else:
      h_atom = by_name.get(outcome.h_name)
    if h_atom is None or h_atom.memory_id() in claimed:
      return
    new_xyz = None
    if outcome.action == MOVE_TO_PARTNER:
      new_xyz = imidazole_nh_xyz(by_name, HIS_PARTNER[donor_atom.name.strip()
                                                     .upper()])
      if new_xyz is None:
        # No ring geometry to place the partner proton.  Keeping the one we
        # have is the safe failure: removing it would leave an imidazolate.
        return
    claimed.add(h_atom.memory_id())
    changes.append(_change(outcome.action, h_atom, donor_atom, metal, distance,
                       outcome.rule, outcome.reason,
                       new_name=outcome.new_name, new_xyz=new_xyz,
                       mark_acceptor=outcome.mark_acceptor))

  soft_metals = [m for m in metal_atoms
                 if normalize_element(m.element) in policy.anion_inducing_metals]

  for atom in hierarchy.atoms():
    ag = atom.parent()
    if ag is None:
      continue
    resname = ag.resname.strip().upper()
    name = atom.name.strip().upper()
    key = (resname, name)
    scored = metal_atoms
    if key in COORD_DEPROT_RULES:
      if key in SOFT_DONOR_KEYS:
        scored = soft_metals
    elif name == "N" and is_standard_aa(resname):
      pass
    else:
      continue
    if not scored:
      continue
    metal, distance = _nearest(atom, scored)
    if distance > coord_cutoff(metal.element, atom.element):
      continue
    h_on_donor = _hydrogens_on(ag, resname, name)
    _emit(atom, ag, dict(
      resname=resname, atom_name=name, element=atom.element,
      metal_element=normalize_element(metal.element), distance=distance,
      h_names=h_on_donor, is_standard_aa=is_standard_aa(resname),
      h_observed=h_observed), metal, distance)

  # Histidine, scored per ring: the coordinating N is whichever ring N is
  # nearest a metal, so a bridging His still picks the right one.
  for ag in hierarchy.atom_groups():
    if ag.resname.strip().upper() != "HIS":
      continue
    by_name = {a.name.strip().upper(): a for a in ag.atoms()}
    best = None
    for ring_n in ("ND1", "NE2"):
      n_atom = by_name.get(ring_n)
      if n_atom is None:
        continue
      metal, distance = _nearest(n_atom, metal_atoms)
      if best is None or distance < best[2]:
        best = (ring_n, n_atom, distance, metal)
    if best is None:
      continue
    ring_n, n_atom, distance, metal = best
    if distance > coord_cutoff(metal.element, n_atom.element):
      continue
    _emit(n_atom, ag, dict(
      resname="HIS", atom_name=ring_n, element=n_atom.element,
      metal_element=normalize_element(metal.element), distance=distance,
      h_names=[h for h in HIS_RING_H[ring_n] if h in by_name],
      partner_h_names=[h for h in HIS_RING_H[HIS_PARTNER[ring_n]]
                       if h in by_name],
      h_observed=h_observed), metal, distance)
  return changes


def _hydrogens_on(ag, resname, donor_name):
  """The hydrogen names the rules care about for this donor, read off *ag*."""
  present = set(a.name.strip().upper() for a in ag.atoms())
  wanted = set()
  wanted |= COORD_DEPROT_RULES.get((resname, donor_name), frozenset())
  if donor_name == "N":
    wanted |= BACKBONE_AMIDE_H | NTERM_AMINE_H
  return sorted(present & wanted)


def apply_changes(hierarchy, changes):
  """Apply *changes* to *hierarchy* in place.  Reference applier.

  Does NOT re-sort the hierarchy: a caller may hold i_seq-indexed structures,
  and silently reordering their atoms would break them.  Sort afterwards if you
  need to.
  """
  n_removed = n_moved = 0
  for change in changes:
    if change.action == REMOVE_H:
      ag = change.h_atom.parent()
      if ag is not None:
        ag.remove_atom(change.h_atom)
        n_removed += 1
    elif change.action == MOVE_TO_PARTNER:
      if change.new_name:
        raw = change.h_atom.name
        lead = len(raw) - len(raw.lstrip())
        change.h_atom.name = (" " * lead + change.new_name
                            + " " * max(0, len(raw) - lead - len(change.new_name)))
      if change.new_xyz is not None:
        change.h_atom.set_xyz(change.new_xyz)
      n_moved += 1
  return group_args(n_removed=n_removed, n_moved=n_moved)
