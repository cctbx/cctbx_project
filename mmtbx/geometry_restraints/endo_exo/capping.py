"""Hydrogen cap placement at covalent boundary sites."""

from __future__ import absolute_import, division, print_function

import math

import iotbx.pdb
import iotbx.pdb.hierarchy
from cctbx.array_family import flex
from scitbx import matrix

#: Names a backbone nitrogen's added hydrogens may take.  The monomer library
#: bonds these to N; a name outside the set leaves the atom in the hierarchy
#: with no bond proxy at all, which is a free hydrogen radical rather than an
#: error, so the names here are restricted to ones it accepts.
AMINE_HYDROGEN_NAMES = ('H1', 'H2', 'H3')


def _free_resseq(chain):
  """Return a residue sequence number no residue in *chain*'s id is using.

  The LOWEST free number, not one past the highest.  A chain with a gap in
  it has unused numbers below its maximum by definition, so this stays
  inside the four characters a PDB file allows, where counting up from the
  maximum reaches hybrid-36 as soon as anything -- a water, typically -- is
  numbered 9999.  Hybrid-36 holds such a value in the hierarchy and mmCIF
  decodes it back to a plain integer, but it is a cctbx convention rather
  than a standard, so it is worth not needing.

  Taken over every chain object sharing the id, not just this one: a TER
  between polymer and solvent splits them, and numbering within one alone
  gives the repair an id the other already holds.
  """
  chain_id = chain.id
  root = chain.parent()
  if root is not None:
    groups = [rg for other in root.chains() if other.id == chain_id
              for rg in other.residue_groups()]
  else:
    groups = list(chain.residue_groups())
  taken = {rg.resseq_as_int() for rg in groups}
  candidate = 1
  while candidate in taken:
    candidate += 1
  # hy36 only past 9999, which needs every number below it to be occupied.
  return iotbx.pdb.hy36encode(width=4, value=candidate)


class HydrogenCapper:
  """Place hydrogen cap atoms at covalent boundary sites.

  Parameters
  ----------
  log : file-like or None, optional
      Destination for diagnostic messages.
  """

  def __init__(self, log=None):
    self.log = log

  def cap_atom(self, anchor, cap):
    """Move *cap* to a hydrogen position 1.1 A along the anchor->cap vector.

    Parameters
    ----------
    anchor : iotbx.pdb.hierarchy.atom or None
    cap : iotbx.pdb.hierarchy.atom or None
        Both may be ``None`` (no-op).
    """
    if anchor is None or cap is None:
      return
    v = flex.double(cap.xyz) - flex.double(anchor.xyz)
    v_norm = v.norm()
    assert v_norm > 1e-6, "anchor and cap must be distinct atoms"
    u = v / v_norm
    cap.set_element('H')
    # cap.set_name(('H' + cap.name.strip()).rjust(4))
    cap.set_xyz(tuple(flex.double(anchor.xyz) + u * 1.1))
    print('Capped atom:', file=self.log)
    print('  ' + anchor.format_atom_record().rstrip(), file=self.log)
    print('  ' + cap.format_atom_record().rstrip(), file=self.log)

  def complete_amine(self, nitrogen, neighbours):
    """Add hydrogens to *nitrogen* until it is a neutral primary amine.

    A backbone nitrogen whose preceding residue is absent from the model
    arrives short of a bond.  Capping cannot help: there is nothing to
    sever and nothing to replace, so the atom reaches the region carrying
    an unpaired electron and a QM code reads it as a radical, an amide
    anion, or whatever its charge heuristic decides.  Hydrogens are added
    until it carries three substituents, which repairs the valence and
    restores the charge the residue has in an unbroken chain.

    What it does not do is rescue a hydrogen bond: a nitrogen near enough
    to the metal to reach here keeps its own amide hydrogen either way, and
    one far enough out to be cut is cut before this runs.  It also stands
    the new hydrogen in the volume the missing carbonyl carbon occupied, so
    a contact that is covalent in the full structure can read as a hydrogen
    bond in the region.  Restoring the acyl group rather than a hydrogen
    would avoid both, and would reach the carbonyl half of the same break.

    As many as are missing, not one: the nitrogen may arrive with a single
    amide hydrogen or with none at all, and adding one to a bare nitrogen
    would leave exactly the two-coordinate centre this exists to remove.
    Counted over substituents rather than hydrogens, so proline's ring
    nitrogen, which already holds CA and CD, gets one and becomes a
    secondary amine.

    These ADD atoms rather than converting one, so they are not caps: they
    have no element to restore and must stay out of ``cap_iseqs``, which
    consumers round-trip through the original element to build restraints.

    Parameters
    ----------
    nitrogen : iotbx.pdb.hierarchy.atom
    neighbours : sequence of iotbx.pdb.hierarchy.atom
        The nitrogen's bonded atoms, used to point each new hydrogen away
        from everything already there.

    Returns
    -------
    list of iotbx.pdb.hierarchy.atom
        The atoms added, empty if the nitrogen needed none or no free name
        was left.
    """
    atom_group = nitrogen.parent()
    if atom_group is None:
      return []
    placed = list(neighbours)
    # Counted over EVERY substituent, not just the hydrogens: a neutral amine
    # nitrogen carries three.  Proline's ring nitrogen already holds CA and CD,
    # so it needs one and asking for two puts the second on top of the first.
    wanted = 3 - len(placed)
    added = []
    for _ in range(max(0, wanted)):
      taken = {a.name.strip().upper() for a in atom_group.atoms()}
      name = next((n for n in AMINE_HYDROGEN_NAMES if n not in taken), None)
      if name is None:
        break

      # Away from everything already on the nitrogen: the sum of the unit
      # vectors to them points into the crowd, so its negation points at the
      # vacancy.  Hydrogens placed here count, or the second would land on
      # the first.
      units = []
      for other in placed:
        v = flex.double(other.xyz) - flex.double(nitrogen.xyz)
        norm = v.norm()
        if norm > 1e-6:
          units.append(v / norm)
      direction = flex.double([0.0, 0.0, 0.0])
      for unit in units:
        direction -= unit
      norm = direction.norm()
      if norm < 1e-6:               # substituents cancel; any direction will do
        direction = flex.double([0.0, 0.0, 1.0])
        norm = 1.0
      direction /= norm

      if len(units) == 1:
        # Anti to a lone substituent is 180 degrees, a linear R-N-H.  Tilt off
        # the axis to the trigonal angle instead; the plane is arbitrary,
        # since nothing else is bonded yet to define one.
        axis = units[0]
        trial = flex.double([1.0, 0.0, 0.0])
        if abs(axis.dot(trial)) > 0.9:
          trial = flex.double([0.0, 1.0, 0.0])
        perpendicular = trial - axis * axis.dot(trial)
        perpendicular /= perpendicular.norm()
        direction = -axis * 0.5 + perpendicular * (3.0 ** 0.5 / 2.0)
        direction /= direction.norm()

      hydrogen = iotbx.pdb.hierarchy.atom()
      hydrogen.name = name.ljust(4)
      hydrogen.element = ' H'
      hydrogen.xyz = tuple(flex.double(nitrogen.xyz) + direction * 1.01)
      hydrogen.occ = nitrogen.occ
      hydrogen.b = nitrogen.b
      # Copied, not defaulted: a blank segid beside siblings carrying spaces
      # makes the residue non-unique and interpretation refuses the region.
      hydrogen.segid = nitrogen.segid
      if any(hydrogen.distance(other) < 0.5 for other in placed):
        break                       # nowhere left to put one; leave it short
      atom_group.append_atom(hydrogen)
      placed.append(hydrogen)
      added.append(hydrogen)
      print('Completed amine:', file=self.log)
      print('  ' + nitrogen.format_atom_record().rstrip(), file=self.log)
      print('  ' + hydrogen.format_atom_record().rstrip(), file=self.log)
    return added

  def complete_carbonyl(self, carbon, neighbours):
    """Cap *carbon* with an NH2 residue, restoring the amide.

    The mirror of :meth:`complete_amine` on the other half of a severed
    peptide bond.  Where the nitrogen half has to guess a torsion, this one
    guesses nothing: the carbon already holds its CA and its O, so the
    missing nitrogen is fixed by sp2 planarity and the two hydrogens by the
    amide plane.  Restoring the amide rather than adding a hydrogen matters
    here because a backbone carbonyl oxygen is often what coordinates the
    metal, and an aldehyde understates its donor strength.

    Placed as a separate residue group, since the parent residue already has
    an atom named N.  ``NH2`` is the monomer library's own C-terminal amide
    cap, so the restraints build without a bespoke definition.

    Parameters
    ----------
    carbon : iotbx.pdb.hierarchy.atom
    neighbours : sequence of iotbx.pdb.hierarchy.atom
        The carbon's bonded atoms; its CA and O fix the amide plane.

    Returns
    -------
    list of iotbx.pdb.hierarchy.atom
        The atoms added, empty if the geometry could not be determined.
    """
    residue_group = carbon.parent().parent()
    chain = residue_group.parent() if residue_group is not None else None
    if chain is None:
      return []
    alpha = next((a for a in neighbours
                  if a.name.strip().upper() == 'CA'), None)
    oxygen = next((a for a in neighbours
                   if a.element.strip().upper() == 'O'), None)
    if alpha is None or oxygen is None:
      return []

    origin = flex.double(carbon.xyz)
    to_alpha = flex.double(alpha.xyz) - origin
    to_oxygen = flex.double(oxygen.xyz) - origin
    if to_alpha.norm() < 1e-6 or to_oxygen.norm() < 1e-6:
      return []
    to_alpha /= to_alpha.norm()
    to_oxygen /= to_oxygen.norm()
    # Trigonal planar at the carbonyl carbon: the third substituent opposes
    # the other two.
    direction = -(to_alpha + to_oxygen)
    if direction.norm() < 1e-6:
      return []
    direction /= direction.norm()
    nitrogen_xyz = origin + direction * 1.335

    # The amide is planar, so both hydrogens lie in the CA/C/O plane at the
    # trigonal angle from the C-N bond.
    normal = flex.double(list(matrix.col(list(to_alpha)).cross(
      matrix.col(list(to_oxygen)))))
    if normal.norm() < 1e-6:
      return []
    normal /= normal.norm()
    back = -direction              # from N toward its carbonyl C
    in_plane = flex.double(list(matrix.col(list(normal)).cross(
      matrix.col(list(back)))))
    in_plane /= in_plane.norm()

    new_group = iotbx.pdb.hierarchy.residue_group(
      resseq=_free_resseq(chain))
    chain.append_residue_group(new_group)
    # The parent's altloc, not blank: a backbone split over two conformers
    # otherwise gets two NH2 groups that both read as whole residues, which
    # interpretation then bonds to each other and to both carbonyls.
    atom_group = iotbx.pdb.hierarchy.atom_group(
      resname='NH2', altloc=carbon.parent().altloc)
    new_group.append_atom_group(atom_group)

    added = []
    for name, xyz in (
        ('N', tuple(nitrogen_xyz)),
        ('HN1', tuple(nitrogen_xyz
                      + (back * -0.5 + in_plane * (3.0 ** 0.5 / 2.0)) * 1.01)),
        ('HN2', tuple(nitrogen_xyz
                      + (back * -0.5 - in_plane * (3.0 ** 0.5 / 2.0)) * 1.01))):
      atom = iotbx.pdb.hierarchy.atom()
      atom.name = name.ljust(4)
      atom.element = (' N' if name == 'N' else ' H')
      atom.xyz = xyz
      atom.occ = carbon.occ
      atom.b = carbon.b
      atom.segid = carbon.segid
      atom_group.append_atom(atom)
      added.append(atom)
    print('Completed carbonyl:', file=self.log)
    print('  ' + carbon.format_atom_record().rstrip(), file=self.log)
    for atom in added:
      print('  ' + atom.format_atom_record().rstrip(), file=self.log)
    return added

  def complete_carboxylate(self, carbon, neighbours):
    """Add OXT to *carbon*, completing a C-terminal carboxylate.

    The other reading of a carbonyl carbon that has lost a bond.  Where a
    later residue exists in the chain the gap is interior and the carbon is
    missing an amide (see :meth:`complete_carbonyl`); where none does, the
    chain genuinely ends there and what is missing is the second carboxylate
    oxygen, which a depositor often leaves unmodelled.

    Determined by the same geometry as the amide nitrogen: the third
    substituent on an sp2 carbon, fixed by the CA and O already present,
    with no free dihedral.  Unlike the amide this needs no new residue
    group -- OXT belongs to the residue's own -- and it gives the residue
    the -1 a C-terminus carries.

    Parameters
    ----------
    carbon : iotbx.pdb.hierarchy.atom
    neighbours : sequence of iotbx.pdb.hierarchy.atom

    Returns
    -------
    list of iotbx.pdb.hierarchy.atom
        The atom added, empty if the geometry could not be determined or an
        OXT is already there.
    """
    atom_group = carbon.parent()
    if atom_group is None:
      return []
    if any(a.name.strip().upper() == 'OXT' for a in atom_group.atoms()):
      return []
    alpha = next((a for a in neighbours
                  if a.name.strip().upper() == 'CA'), None)
    oxygen = next((a for a in neighbours
                   if a.element.strip().upper() == 'O'), None)
    if alpha is None or oxygen is None:
      return []

    origin = flex.double(carbon.xyz)
    to_alpha = flex.double(alpha.xyz) - origin
    to_oxygen = flex.double(oxygen.xyz) - origin
    if to_alpha.norm() < 1e-6 or to_oxygen.norm() < 1e-6:
      return []
    to_alpha /= to_alpha.norm()
    to_oxygen /= to_oxygen.norm()
    direction = -(to_alpha + to_oxygen)
    if direction.norm() < 1e-6:
      return []
    direction /= direction.norm()

    hydroxyl = iotbx.pdb.hierarchy.atom()
    hydroxyl.name = 'OXT '
    hydroxyl.element = ' O'
    hydroxyl.xyz = tuple(origin + direction * 1.249)
    hydroxyl.occ = carbon.occ
    hydroxyl.b = carbon.b
    hydroxyl.segid = carbon.segid
    atom_group.append_atom(hydroxyl)
    print('Completed carboxylate:', file=self.log)
    print('  ' + carbon.format_atom_record().rstrip(), file=self.log)
    print('  ' + hydroxyl.format_atom_record().rstrip(), file=self.log)
    return [hydroxyl]
