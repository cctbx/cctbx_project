"""Bond-cut heuristics (backbone / sidechain) for locating good
hydrogen-capping sites."""

from __future__ import absolute_import, division, print_function

from collections import deque

from mmtbx.geometry_restraints.endo_exo.util import _neighbour_iseqs


PREFERRED_CUTS = {
  'ALA': {'CA', 'CB'},
  'ARG': {'CD', 'CG'},
  'ASN': {'CA', 'CB'},
  'ASP': {'CA', 'CB'},
  'CYS': {'CA', 'CB'},
  'GLN': {'CB', 'CG'},
  'GLU': {'CB', 'CG'},
  'GLY': None,
  'HIS': {'CA', 'CB'},
  'ILE': {'CA', 'CB'},
  'LEU': {'CB', 'CG'},
  'LYS': {'CD', 'CE'},
  'MET': {'CB', 'CG'},
  'PHE': {'CA', 'CB'},
  # Proline's sidechain closes back onto its own N, so no single bond
  # detaches it: cutting one ring bond severs nothing and cutting two
  # strands the atoms between them.  An empty entry means the sidechain is
  # never cut, leaving the backbone rules to bound the residue.
  'PRO': frozenset(),
  'SER': {'CA', 'CB'},
  'THR': {'CA', 'CB'},
  'TRP': {'CA', 'CB'},
  'TYR': {'CA', 'CB'},
  'VAL': {'CA', 'CB'},
}


class BondCutDetector:
  """Identify covalent bonds that are good hydrogen-capping sites.

  Parameters
  ----------
  use_preferred_cuts : bool, optional
      When ``True`` (default) the ``PREFERRED_CUTS`` lookup table is
      consulted before falling back to geometric heuristics.
  log : file-like or None, optional
      Destination for diagnostic messages.
  """

  def __init__(self, use_preferred_cuts=True, log=None):
    self.use_preferred_cuts = use_preferred_cuts
    self.log = log

  # ------------------------------------------------------------------
  # Public interface
  # ------------------------------------------------------------------

  def is_cc_single_sp3_bond(self, resname, atom1, atom2, adjacency,
                             atoms_by_i_seq=None):
    """Return ``True`` if the bond atom1->atom2 is a suitable C-C capping site.

    The residue's ``PREFERRED_CUTS`` entry is consulted first: it names the
    bond that retains the group a coordinating residue is in the region
    for, and BFS arriving from that group reaches it first.  A bond the
    entry does not name falls through to a geometric heuristic, which cuts
    at any sp3 C-C: both atoms carbon, bonded, separated by 1.42-1.68 A,
    both of degree 4, and neither carbon looking unsaturated.  That is what
    trims a residue held only as scaffold, whose preferred bond BFS never
    reaches from the sidechain side.

    An entry may be empty, meaning the residue has no cuttable sidechain
    bond at all; the backbone rules then bound it on their own.  Setting
    ``use_preferred_cuts=False``, or a *resname* absent from the table,
    leaves only the geometric heuristic.

    Parameters
    ----------
    resname : str
        Three-letter residue name (upper-case).
    atom1 : iotbx.pdb.hierarchy.atom
    atom2 : iotbx.pdb.hierarchy.atom
    adjacency : dict
        Local covalent graph.
    atoms_by_i_seq : dict or None, optional
        Map ``{i_seq: atom}`` used for the unsaturation check.  Only
        consulted in the geometric branch.

    Returns
    -------
    bool
    """
    if self.use_preferred_cuts:
      preferred = PREFERRED_CUTS.get(resname)
      if preferred is not None:
        if (atom1.name.strip().upper() in preferred
            and atom2.name.strip().upper() in preferred):
          return True
        if not self._is_backbone_ward_of(preferred, atom1, atom2, adjacency):
          return False

    if atom1.element.strip().upper() != 'C':
      return False
    if atom2.element.strip().upper() != 'C':
      return False
    nbr1 = _neighbour_iseqs(adjacency, atom1.i_seq)
    # Defensive: the BFS only offers bonded pairs, but bad input (a missing
    # bond proxy) would otherwise reach the degree test below.
    if atom2.i_seq not in nbr1:
      return False

    cc_dist = atom1.distance(atom2)
    if not (1.42 <= cc_dist <= 1.68):
      return False

    if len(nbr1) != 4 or len(_neighbour_iseqs(adjacency, atom2.i_seq)) != 4:
      return False

    if atoms_by_i_seq is not None:
      if self._looks_unsaturated(atom1, adjacency, atoms_by_i_seq) or \
         self._looks_unsaturated(atom2, adjacency, atoms_by_i_seq):
        return False

    return True

  def is_ca_c_bond(self, atom1, atom2, adjacency):
    """Return ``True`` if atom1->atom2 is a backbone CA-C bond.

    Parameters
    ----------
    atom1 : iotbx.pdb.hierarchy.atom
    atom2 : iotbx.pdb.hierarchy.atom
    adjacency : dict

    Returns
    -------
    bool
    """
    if not (atom1.name.strip().upper() == 'CA' and
            atom2.name.strip().upper() == 'C'):
      return False
    return atom2.i_seq in _neighbour_iseqs(adjacency, atom1.i_seq)

  def is_ca_n_bond(self, atom1, atom2, adjacency):
    """Return ``True`` if atom1->atom2 is a backbone CA-N bond.

    Parameters
    ----------
    atom1 : iotbx.pdb.hierarchy.atom
    atom2 : iotbx.pdb.hierarchy.atom
    adjacency : dict

    Returns
    -------
    bool
    """
    if not (atom1.name.strip().upper() == 'CA' and
            atom2.name.strip().upper() == 'N'):
      return False
    return atom2.i_seq in _neighbour_iseqs(adjacency, atom1.i_seq)

  # ------------------------------------------------------------------
  # Private helpers
  # ------------------------------------------------------------------

  @staticmethod
  def _hops_to_ca(atom, adjacency):
    """Return ``{atom name: bonds from CA}`` within *atom*'s residue group.

    Empty if the residue has no CA, or if *atom* carries no hierarchy to
    ask (a standalone atom built outside a model).  Hydrogens are
    included; they are leaves and never affect a comparison between two
    heavy atoms.

    Parameters
    ----------
    atom : iotbx.pdb.hierarchy.atom
    adjacency : dict

    Returns
    -------
    dict
    """
    atom_group = atom.parent()
    residue_group = None if atom_group is None else atom_group.parent()
    if residue_group is None:
      return {}
    in_residue = {}
    for atom_group in residue_group.atom_groups():
      for residue_atom in atom_group.atoms():
        in_residue[residue_atom.i_seq] = residue_atom.name.strip().upper()
    ca_iseq = next((i for i, n in in_residue.items() if n == 'CA'), None)
    if ca_iseq is None:
      return {}

    hops = {ca_iseq: 0}
    queue = deque([ca_iseq])
    while queue:
      current = queue.popleft()
      for neighbour in _neighbour_iseqs(adjacency, current):
        if neighbour in hops or neighbour not in in_residue:
          continue
        hops[neighbour] = hops[current] + 1
        queue.append(neighbour)
    return {in_residue[i]: h for i, h in hops.items()}

  @classmethod
  def _is_backbone_ward_of(cls, preferred, atom1, atom2, adjacency):
    """Return ``True`` if this bond lies nearer the backbone than *preferred*.

    The table entry names the bond that retains a residue's functional
    group, so a bond further out along the sidechain must not be cut
    instead: doing so severs the group the entry exists to keep, and can
    strand a lone carbon between two cuts.  A bond nearer the backbone is
    fair game, and is what trims a residue BFS reached from the mainchain.

    Measured in bonds from CA, so it does not depend on which end of the
    bond BFS arrived from.  Residues with no CA fall through as ``True``,
    leaving the geometric heuristic in charge.

    Parameters
    ----------
    preferred : set of str
        The residue's ``PREFERRED_CUTS`` entry.
    atom1, atom2 : iotbx.pdb.hierarchy.atom
    adjacency : dict

    Returns
    -------
    bool
    """
    if not preferred:
      return False  # an empty entry means this sidechain is never cut
    hops = cls._hops_to_ca(atom1, adjacency)
    if not hops:
      return True
    entry_depth = min((hops[name] for name in preferred if name in hops),
                      default=None)
    if entry_depth is None:
      return True
    bond_depth = min(
      (hops[a.name.strip().upper()] for a in (atom1, atom2)
       if a.name.strip().upper() in hops), default=None)
    return bond_depth is not None and bond_depth < entry_depth

  @staticmethod
  def _looks_unsaturated(carbon_atom, adjacency, atoms_by_i_seq):
    """Return ``True`` if *carbon_atom* appears to be in an unsaturated environment.

    Parameters
    ----------
    carbon_atom : iotbx.pdb.hierarchy.atom
    adjacency : dict
    atoms_by_i_seq : dict

    Returns
    -------
    bool
    """
    for neighbor_i_seq in _neighbour_iseqs(adjacency, carbon_atom.i_seq):
      neighbor_atom = atoms_by_i_seq.get(neighbor_i_seq)
      if neighbor_atom is None:
        continue
      neighbor_element = neighbor_atom.element.strip().upper()
      dist = carbon_atom.distance(neighbor_atom)
      if neighbor_element in {'O', 'N', 'S'} and dist < 1.34:
        return True
      if neighbor_element == 'C' and dist < 1.42:
        return True
    return False
