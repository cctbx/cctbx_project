"""Bond-cut heuristics (backbone / sidechain) for locating good
hydrogen-capping sites."""

from __future__ import absolute_import, division, print_function

from mmtbx.geometry_restraints.endoexo.util import _neighbour_iseqs


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
  'PRO': None,
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

    Two modes:

    * If ``use_preferred_cuts`` is ``True`` *and* *resname* has an entry in
      ``PREFERRED_CUTS``, the bond is suitable iff both atom names appear in
      that entry (the geometric heuristic is not consulted).
    * Otherwise (``use_preferred_cuts=False`` or *resname* absent from
      ``PREFERRED_CUTS``) a geometric / valence heuristic is applied:
      both atoms must be carbon, bonded, separated by 1.42-1.68 A, both of
      degree 4, and neither carbon may look unsaturated.

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
        consulted in the heuristic branch.

    Returns
    -------
    bool
    """
    if self.use_preferred_cuts and PREFERRED_CUTS.get(resname) is not None:
      preferred = PREFERRED_CUTS[resname]
      if (atom1.name.strip().upper() in preferred and
          atom2.name.strip().upper() in preferred):
        # print(f'Found preferred cut atoms for {resname}:', file=self.log)
        # print(f'{atom1.format_atom_record().rstrip()}', file=self.log)
        # print(f'{atom2.format_atom_record().rstrip()}', file=self.log)
        #   f'format_atom_record().rstrip()'
        #   f'{atom1.name.strip().upper()} and {atom2.name.strip().upper()}',
        #   file=self.log,
        # )
        return True
      return False
    else:
      if atom1.element.strip().upper() != 'C':
        return False
      if atom2.element.strip().upper() != 'C':
        return False
      nbr1 = _neighbour_iseqs(adjacency, atom1.i_seq)
      # The check below is purely defensive since the BFS should only call this method on bonded pairs, but it
      # guards against bad input data (e.g. missing bond proxies) and prevents a KeyError in that case.
      if atom2.i_seq not in nbr1:
        return False

      cc_dist = atom1.distance(atom2)
      if not (1.42 <= cc_dist <= 1.68):
        return False

      deg1 = len(nbr1)
      deg2 = len(_neighbour_iseqs(adjacency, atom2.i_seq))
      if deg1 != 4 or deg2 != 4:
      # if deg1 < 2 or deg1 > 4 or deg2 < 2 or deg2 > 4:
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
