"""Bond-graph perception of metal-coordinating donor atoms.

Answers, from the BOND GRAPH rather than from a residue-name whitelist, the
questions a caller must settle before it can reason about a metal's first shell:

- what formal charge does this donor's resonance group carry?  (an oxidation
  predictor needs the anionic charge donated to the metal)
- how many sigma connections does it have, and which of them are hydrogens?
  (a protonation program needs this to decide whether an H is removable)

Both consumers otherwise reach for a residue-name table, which reads a cofactor,
buffer or substrate ligand as nothing at all.  The bond graph does not care what
the residue is called.

RESONANCE GROUPS.  A Lewis-structure walk parks a delocalised charge on ONE
arbitrary member of the group: for a carboxylate the -1 lands on whichever
oxygen the monomer library made single-bonded, which is the metal-facing one
only about half the time.  Reading the coordinating ATOM's charge therefore
misses roughly half of all carboxylates, and is WORSE than the whitelist it
would replace.  :func:`resonance_group` collects the equivalent terminal atoms
so the charge can be read off the group, where it is well defined.  For the
twenty standard amino acids the group and the residue coincide (one carboxylate
per Asp/Glu, one phenolate per Tyr), so this generalises the per-residue
convention rather than replacing it; on a multi-arm ligand such as citrate, with
two coordinating carboxylate groups on ONE residue, the two differ and the group
is the chemically meaningful unit.

Inputs are plain data (atom elements, bond pairs, a charge map) rather than a
model or a hierarchy, so the rules are testable without constructing a
structure, and the module stays a leaf: it imports nothing from mmtbx.
"""


_HYDROGEN = frozenset(("H", "D"))


def adjacency(bond_pairs, n_atoms=None):
  """Undirected adjacency ``{i: set(j)}`` from an iterable of ``(i, j)`` bond
  index pairs.  Isolated atoms are absent unless *n_atoms* is given, in which
  case every index below it gets an entry."""
  adj = {}
  if n_atoms is not None:
    for i in range(n_atoms):
      adj[i] = set()
  for i, j in bond_pairs:
    adj.setdefault(i, set()).add(j)
    adj.setdefault(j, set()).add(i)
  return adj


class DonorPerception(object):
  """Bond-graph facts about individual atoms.

  Parameters
  ----------
  elements : sequence of str
      Element symbol per atom index (any capitalisation).
  bond_pairs : iterable of (int, int)
      The covalent bond graph.  Metal-ligand bonds must be EXCLUDED: a
      coordinating N that counts the metal as a neighbour looks saturated, and
      every sigma-count rule downstream then misfires.  Callers building the
      graph through pdb_interpretation get this by setting
      ``automatic_linking.link_metals = False``.
  charges : dict, optional
      ``{atom index: formal charge}``, sign as a chemist writes it (carboxylate
      O is -1).  Absent indices are neutral.
  metal_iseqs : iterable of int, optional
      Excluded from every neighbour list, so a metal that IS bonded in the
      supplied graph still cannot be mistaken for a substituent, and its own
      charge can never be absorbed into a donor group.
  unresolved_iseqs : iterable of int, optional
      Atoms whose charge the supplying Lewis walk never solved, typically a
      conjugated ring it gave up on.  Their charge is whichever partial value
      the walk stopped at, so :func:`group_is_resolved` reports it as untrusted
      rather than letting a caller read it as a result.
  """

  def __init__(self, elements, bond_pairs, charges=None, metal_iseqs=(),
               unresolved_iseqs=()):
    self.elements = [str(e).strip().capitalize() for e in elements]
    self.charges = dict(charges or {})
    self.metals = frozenset(metal_iseqs)
    self.unresolved = frozenset(unresolved_iseqs)
    self._adj = adjacency(bond_pairs, n_atoms=len(self.elements))

  def element(self, i):
    return self.elements[i]

  def is_hydrogen(self, i):
    return self.elements[i] in _HYDROGEN

  def neighbours(self, i):
    """Every bonded neighbour except a metal."""
    return sorted(j for j in self._adj.get(i, ()) if j not in self.metals)

  def heavy_neighbours(self, i):
    """Bonded non-hydrogen, non-metal neighbours."""
    return [j for j in self.neighbours(i) if not self.is_hydrogen(j)]

  def bonded_h(self, i):
    """Bonded hydrogens (or deuteriums)."""
    return [j for j in self.neighbours(i) if self.is_hydrogen(j)]

  def n_sigma(self, i):
    """Sigma connection count: heavy neighbours plus bonded H, metals excluded.

    This is the quantity that says whether a donor's valence is already
    satisfied without the metal, and so whether an H on it is chemically
    required or is an artefact of protonating in the metal's absence."""
    return len(self.neighbours(i))

  def resonance_group(self, i):
    """The atoms over which *i*'s formal charge is delocalised, as a sorted
    tuple including *i*.

    A TERMINAL donor (one heavy neighbour) is grouped with the same-element
    terminal siblings on that neighbour: the two oxygens of a carboxylate, the
    three of a phosphate.  Anything else (a hydroxyl bridging nothing, an ether
    or ring O, a thioether S) is its own group.  Deliberately shallow: it walks
    one bond, so it will not merge the two arms of a citrate, which are
    genuinely independent charges and must stay independent.
    """
    heavy = self.heavy_neighbours(i)
    if len(heavy) != 1:
      return (i,)
    parent = heavy[0]
    el = self.elements[i]
    group = [j for j in self.neighbours(parent)
             if self.elements[j] == el and len(self.heavy_neighbours(j)) == 1]
    return tuple(sorted(set(group) | {i}))

  def group_charge(self, i):
    """Net formal charge of *i*'s functional group: its resonance group plus the
    atom that group hangs off.

    This, not ``charges[i]``, is the charge a coordinating atom presents to a
    metal: see the module docstring on where a Lewis walk parks a delocalised
    charge.

    The PARENT has to be included, because a group's charge is not always
    entirely on its terminal atoms.  A nitro group is the case that shows it:
    N is formally +1 and one O is -1, so the two oxygens alone read -1 while
    the fragment is neutral.  Adding the parent leaves a carboxylate (-1),
    phosphate (-3) and sulfate (-2) untouched, since their parent atoms are
    formally neutral, and corrects nitro to 0 and a zwitterionic sulfate
    (S(2+) with four O(-)) from -4 to -2.
    """
    total = sum(self.charges.get(j, 0) for j in self.resonance_group(i))
    parent = self.group_parent(i)
    if parent is not None:
      total += self.charges.get(parent, 0)
    return total

  def group_parent(self, i):
    """The heavy atom *i*'s resonance group hangs off, or None for a group that
    is not terminal (and so has no single parent)."""
    heavy = self.heavy_neighbours(i)
    return heavy[0] if len(heavy) == 1 else None

  def group_is_resolved(self, i):
    """False when *i*'s resonance group, or the atom it hangs off, sits on a
    ring the Lewis walk never solved, so :func:`group_charge` there is a
    stopping point rather than an answer.

    The parent matters as much as the group: a terminal O on an unsolved
    aromatic ring is its own group and never intersects the ring, so testing the
    group alone would report an unsolved phenol-type charge as trustworthy."""
    touched = set(self.resonance_group(i))
    parent = self.group_parent(i)
    if parent is not None:
      touched.add(parent)
    return not (touched & self.unresolved)

  def group_charges(self, donor_iseqs):
    """``{group: charge}`` over the resonance groups of *donor_iseqs*, so a
    bidentate donor whose two oxygens belong to one carboxylate is counted
    ONCE.  Keys are the group tuples from :func:`resonance_group`."""
    out = {}
    for i in donor_iseqs:
      out[self.resonance_group(i)] = self.group_charge(i)
    return out


def from_electron_distribution(atoms, electron_distribution, metal_iseqs=()):
  """Build a :class:`DonorPerception` from an ``mmtbx.ligands.electrons``
  result.

  *atoms* is the parallel atom array (``hierarchy.atoms()``).  The bond graph is
  read from the distribution's tuple keys, and the charges from
  ``get_charged_atoms()``, whose value is a SURPLUS ELECTRON COUNT and so is
  negated here to give a formal charge.  Rings the walk never solved are carried
  through as untrusted (see :func:`DonorPerception.group_is_resolved`).
  """
  bonds = [k for k in electron_distribution if isinstance(k, tuple)]
  charges = {}
  for atom, n_electrons in electron_distribution.get_charged_atoms():
    charges[atom.i_seq] = -n_electrons
  unresolved = set()
  for cycle in getattr(electron_distribution, "unresolved_cycles", ()):
    unresolved.update(cycle)
  return DonorPerception([a.element for a in atoms], bonds,
                         charges=charges, metal_iseqs=metal_iseqs,
                         unresolved_iseqs=unresolved)
