"""mmtbx.endoexo — QM region builder with BFS expansion and hydrogen capping.

Grows a QM region around each seed site (metal atoms by default, or a
user-supplied selection) by breadth-first traversal of the covalent graph,
optionally caps dangling bonds with hydrogen atoms, estimates the net charge
of the surrounding region, and writes a PDB file, an mmCIF file, and a
sidecar PHIL file per seed.

Usage (via dispatcher):

    mmtbx.endoexo model.pdb
    mmtbx.endoexo model.pdb selection="chain A and resseq 100" radius=5.0
    mmtbx.endoexo model.pdb selection="chain A and resseq 100" \\
                             selection="chain B and resseq 200"

Programmatic use::

    from mmtbx.programs.endoexo import Program
    from iotbx.cli_parser import run_program
    run_program(program_class=Program)
"""

from __future__ import absolute_import, division, print_function

import os
from collections import defaultdict, deque

try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  from libtbx.program_template import ProgramTemplate

import libtbx.phil
from libtbx.utils import Sorry
# from iotbx.cli_parser import run_program

from cctbx import sgtbx
from cctbx.array_family import flex
from iotbx.pdb import common_residue_names_get_class
from scitbx import matrix

from scipy.spatial import ConvexHull, KDTree

from mmtbx.qmi import metals as qmi_metals


def _canon_op(op):
  """Return a hash-stable copy of *op*.

  cctbx's ``sgtbx.rt_mx`` has an inconsistent ``__hash__`` vs
  ``__eq__`` after composition: two operations that compare equal
  can hash to different values, which breaks set / dict membership.
  Round-tripping through the canonical xyz string fixes it.
  """
  return sgtbx.rt_mx(op.as_xyz())


def _neighbour_iseqs(adjacency, i_seq):
  """Return the set of bare neighbour i_seqs for *i_seq* in the tagged
  adjacency, dropping the per-edge ``sym_op``.

  Used by bond-cut detection and degree counting where only covalent
  connectivity matters, not which symmetry image the neighbour belongs
  to.
  """
  return {j for (j, _op) in adjacency.get(i_seq, set())}

STANDARD_RESIDUE_CHARGES = {
  'ALA':  0,
  'ARG':  1,
  'ASN':  0,
  'ASP': -1,
  'CYS':  0,
  'GLN':  0,
  'GLU': -1,
  'GLY':  0,
  'HIS':  None,  # can be +1, 0 or -1 depending on protonation state
  'ILE':  0,
  'LEU':  0,
  'LYS':  1,
  'MET':  0,
  'PHE':  0,
  'PRO':  0,
  'SER':  0,
  'THR':  0,
  'TRP':  0,
  'TYR':  0,
  'VAL':  0,
}

# formal_charge is the residue's charge in its standard RCSB CCD reference
# form (ASP/GLU/CYS/TYR neutral; LYS/ARG/HIS+ protonated).  protonation_sites
# lists every canonical H name that, when present, contributes +1 to the
# reference charge -- one canonical name per chemically distinct H slot.
# Net sidechain charge = formal_charge - len(canonical_sites) + n_present.
CHARGED_SIDECHAINS = {
  'ASP': {
    'formal_charge': 0,
    'charged_heavy_atoms': {'OD1', 'OD2'},
    'protonation_sites': {'HD2', 'DD2'},
  },
  'CYS': {
    'formal_charge': 0,
    'charged_heavy_atoms': {'SG'},
    'protonation_sites': {'HG', 'DG'},
  },
  'GLU': {
    'formal_charge': 0,
    'charged_heavy_atoms': {'OE1', 'OE2'},
    'protonation_sites': {'HE2', 'DE2'},
  },
  'LYS': {
    'formal_charge': 1,
    'charged_heavy_atoms': {'NZ'},
    'protonation_sites': {'HZ1', 'HZ2', 'HZ3', 'DZ1', 'DZ2', 'DZ3'},
  },
  'ARG': {
    'formal_charge': 1,
    'charged_heavy_atoms': {'NE', 'NH1', 'NH2', 'CZ'},
    'protonation_sites': {'HH11', 'HH12', 'HH21', 'HH22', 'DH11', 'DH12', 'DH21', 'DH22'},
  },
  'HIS': {
    'formal_charge': 1,
    'charged_heavy_atoms': {'ND1', 'NE2'},
    'protonation_sites': {'HD1', 'HE2', 'DD1', 'DE2'},
  },
  'TYR': {
    'formal_charge': 0,
    'charged_heavy_atoms': {'OH'},
    'protonation_sites': {'HH', 'DH'},
  },
}

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


# =============================================================================
# PHIL parameter string
# =============================================================================

master_phil_str = """
selection = None
  .type = str
  .multiple = True
  .help = "Atom selection string(s) for the initial seed region(s). Each entry produces one QM region output file. \
May be specified multiple times (e.g. selection='chain A and resseq 100' selection='chain B and resseq 200'). \
If no selection is given, all metals in the structure are used as seeds (one output per metal)."
metal_element = None
  .type = str
  .multiple = True
  .help = "Restrict the metal-scan seed search to atoms of these element(s) \
(e.g. metal_element=Fe, or metal_element=Fe metal_element=Cu). Element symbols \
are case-insensitive. Only consulted when no `selection` is given; if any \
selection is provided, the explicit selection wins. Default (no values): \
seed on every metal found by mmtbx.qmi.metals."
altloc = auto
  .type = str
  .help = "Which altloc letter to retain per residue.  'auto' (default) picks the \
highest mean-occupancy non-blank altloc per residue.  A specific letter (e.g. 'A', 'B') \
keeps that letter, falling back to the highest-occupancy altloc with a warning if the \
letter is absent from a residue.  'all' disables altloc filtering."
radius = 5.0
  .type = float(value_min=0)
  .help = "Radius of the buffer region around the selected scatterer."
skip_radius_search = False
  .type = bool
  .help = "If True, the initial radius search is skipped and only the seed atoms themselves seed the QM region (BFS expansion still applies)."
metal_ligand_cutoff = 3.0
  .type = float(value_min=0)
  .help = "Distance cutoff in Angstrom used to add fallback metal-ligand edges when bond proxies do not include coordination bonds."
max_depth = 3
  .type = int(value_min=0)
  .help = "Maximum BFS depth from any atom within the QM region."
use_preferred_cuts = True
  .type = bool
  .help = "Whether to use preferred cut atoms for each residue type when identifying candidate bonds for capping, instead of relying on heuristics alone."
include_waters_in_convex_hull = True
  .type = bool
  .help = "Whether to check for water molecules inside the convex hull of the selected QM region and add them to the QM region if found."
do_capping = True
  .type = bool
  .help = "Whether to perform capping of boundary atoms based on heuristics. If False, the output QM region will have uncapped dangling bonds."
residues_to_include = None
  .type = str
  .help = "Selection string for residues to almost always include in the output (depending on chain specified in the selection string),
  regardless of the sidechain rules. For example: 'chain A and resseq 50-100'."
include_terminal_charges = False
  .type = bool
  .help = "If True, estimate free peptide termini charges inside the truncated QM region."
write_files = True
  .type = bool
  .help = "If True (default), write a PDB, an mmCIF, and a sidecar PHIL \
file per seed to the current working directory. Set to False when calling \
the program in-memory (e.g. via Program(...).run() + get_results()) and \
the per-seed Model objects are consumed directly without a disk round-trip."
n_terminus_charge = 1
  .type = int(value_min=0, value_max=1)
  .help = "Charge assigned to each detected free N-terminus when include_terminal_charges=True."
c_terminus_charge = -1
  .type = int(value_min=-1, value_max=0)
  .help = "Charge assigned to each detected free C-terminus when include_terminal_charges=True."
"""


# Sidecar PHIL emitted alongside each QM-region PDB/mmCIF.  Indices are
# 0-based positional indices into the output PDB/mmCIF atom list.
master_sidecar_phil_str = """
endoexo_region {
  cap_atoms = None
    .type = ints
    .help = "Indices of hydrogen cap atoms in the QM-region output file."
  cap_original_elements = None
    .type = strings
    .help = "Element each cap atom carried before being replaced by H, \
parallel to cap_atoms.  Downstream consumers that rebuild restraints on \
the sub-model need to restore these elements for the duration of \
pdb_interpretation (so the monomer-library lookup matches), then \
switch back to H for electron counting."
  seed_atoms = None
    .type = ints
    .help = "Indices of seed atoms in the QM-region output file (the metal \
atom in metal-scan mode, or all atoms matched by the selection string)."
  selection_string = None
    .type = str
    .help = "Original CCTBX selection string that seeded the region; \
None for metal-scan regions."
}
"""


# =============================================================================
# SeedFinder
# =============================================================================

class SeedFinder:
  """Locate seed atoms within a CCTBX model hierarchy.

  Seeds are either all metal atoms in the structure (default) or the atoms
  matched by one or more user-supplied CCTBX selection strings.
  """

  def find_metals(self, model, element_filter=None):
    """Return metal atom objects in *model*.

    Delegates to ``mmtbx.qmi.metals.metal_atoms``. When *element_filter*
    is given, its element symbols are passed straight through as the
    ``metals=`` argument, restricting the scan to those element(s).
    When ``None``, the canonical ``METALS`` recognition list is used.

    Parameters
    ----------
    model : mmtbx.model.manager
    element_filter : iterable of str or None, optional
        Element symbols (case-/whitespace-tolerant, matched against
        ``element.strip().capitalize()``).

    Returns
    -------
    list of iotbx.pdb.hierarchy.atom
    """
    if element_filter:
      wanted = {
        str(e).strip().capitalize()
        for e in element_filter if str(e).strip()
      }
      if wanted:
        return qmi_metals.metal_atoms(model, metals=wanted)
    return qmi_metals.metal_atoms(model)

  def find_by_selection(self, model, selection_str):
    """Return atom objects matched by *selection_str*.

    Parameters
    ----------
    model : mmtbx.model.manager
    selection_str : str
        CCTBX atom selection string.

    Returns
    -------
    list of iotbx.pdb.hierarchy.atom
    """
    mask = model.selection(selection_str)
    atoms = model.get_hierarchy().atoms()
    return [atoms[i] for i in mask.iselection()]

  def find(self, model, selection_strings=None, element_filter=None):
    """Return seed groups for the given model.

    Each group is a ``(label, atoms)`` tuple where *label* is the selection
    string that produced the group (or ``None`` for metal-scan groups) and
    *atoms* is the list of seed atoms for that group.

    When *selection_strings* is a non-empty list each entry produces one
    group, and *element_filter* is ignored.  When it is empty or ``None``
    every metal atom in the model becomes its own group, optionally
    restricted to the element(s) listed in *element_filter*.

    Parameters
    ----------
    model : mmtbx.model.manager
    selection_strings : list of str or None, optional
    element_filter : iterable of str or None, optional
        Element symbols (case-insensitive) restricting the metal-scan.
        Only consulted when *selection_strings* is empty.

    Returns
    -------
    list of tuple
        Each element is ``(str or None, list of iotbx.pdb.hierarchy.atom)``.
    """
    if selection_strings:
      return [
        (sel_str, self.find_by_selection(model, sel_str))
        for sel_str in selection_strings
      ]
    return [
      (None, [m])
      for m in self.find_metals(model, element_filter=element_filter)
    ]

  def count_metals(self, model):
    """Return the number of metal atoms in *model*.

    Parameters
    ----------
    model : mmtbx.model.manager

    Returns
    -------
    int
    """
    return qmi_metals.count_metals(model)


# =============================================================================
# AtomGraphBuilder
# =============================================================================

class AtomGraphBuilder:
  """Build a unified, symmetry-aware covalent adjacency graph from CCTBX
  restraint proxies.

  Each edge is tagged with the ``rt_mx_ji`` operation that brings atom
  *j* into contact with atom *i*. Intra-ASU bonds carry the identity
  operation; symmetry-crossing bonds carry the corresponding
  ``rt_mx_ji`` from ``bond_proxies_asu``.

  This unification lets the downstream BFS treat lattice translations
  and point-group symmetry uniformly.
  """

  def build_adjacency(self, bond_proxies_simple, bond_proxies_asu,
                      asu_mappings):
    """Build a tagged adjacency graph combining intra-ASU and
    symmetry-crossing bond proxies.

    Parameters
    ----------
    bond_proxies_simple : iterable
        Simple (intra-ASU) bond proxy objects with ``.i_seqs``.
    bond_proxies_asu : iterable
        ASU bond proxy objects with ``.i_seq``, ``.j_seq``, ``.j_sym``.
        The corresponding ``rt_mx_ji`` is resolved via *asu_mappings*.
    asu_mappings : cctbx.crystal.direct_space_asu.asu_mappings
        Maps each proxy's ``j_sym`` index to the ``rt_mx_ji`` operation
        that brings *j* into contact with *i*; obtained from
        ``grm.pair_proxies(sites_cart=...).bond_proxies.asu_mappings()``.

    Returns
    -------
    collections.defaultdict of set
        ``adj[i_seq]`` is a set of ``(j_seq, rt_mx_ji)`` tuples. Local
        (intra-ASU) edges carry the identity ``rt_mx``; the reverse
        edge ``j -> i`` is recorded with the inverse operation.
    """
    identity = _canon_op(sgtbx.rt_mx())
    adj = defaultdict(set)
    for proxy in bond_proxies_simple:
      i_seq, j_seq = proxy.i_seqs
      adj[i_seq].add((j_seq, identity))
      adj[j_seq].add((i_seq, identity))
    for proxy in bond_proxies_asu:
      i_seq, j_seq = proxy.i_seq, proxy.j_seq
      sym_op = _canon_op(asu_mappings.get_rt_mx_ji(proxy))
      adj[i_seq].add((j_seq, sym_op))
      adj[j_seq].add((i_seq, _canon_op(sym_op.inverse())))
    return adj

  def add_seed_contact_edges(self, seeds, model, adjacency, cutoff):
    """Add distance-based seed-ligand edges to *adjacency* in-place,
    including symmetry images.

    Uses ``cctbx.crystal.pair_asu_table`` to find pairs within
    *cutoff* of each seed including those reached by point-group or
    lattice-translation operations.  Each edge is tagged with the
    ``rt_mx_ji`` that maps the ligand atom into contact with the
    seed.  Hydrogen and deuterium atoms are skipped.  Edges that
    already exist (same ``(j_seq, sym_op)`` after canonicalisation)
    are not duplicated.

    Parameters
    ----------
    seeds : list of iotbx.pdb.hierarchy.atom
        Seed atom objects.
    model : mmtbx.model.manager
    adjacency : collections.defaultdict of set
        Tagged adjacency; modified in-place.
    cutoff : float
        Distance threshold in Angstrom.

    Returns
    -------
    int
        Number of new edges added.
    """
    atoms = model.get_hierarchy().atoms()
    xs = model.get_xray_structure()
    pat = xs.pair_asu_table(distance_cutoff=cutoff)
    # ``all_interactions_from_inside_asu=True`` is required to enumerate
    # every close contact of an ASU atom: the default collapses pairs
    # that are sym-equivalent under the site symmetry of the seed (e.g.
    # the three Asp 93 contacts to an Fe sitting on a 3-fold axis would
    # otherwise be returned as one representative).
    sym_table = pat.extract_pair_sym_table(
      skip_j_seq_less_than_i_seq=False,
      all_interactions_from_inside_asu=True)
    added_edges = 0

    for seed in seeds:
      i_seq = seed.i_seq
      for j_seq, rt_mx_ji_list in dict(sym_table[i_seq]).items():
        if j_seq == i_seq:
          continue
        atom = atoms[j_seq]
        if atom.element.strip().upper() in ('H', 'D'):
          continue
        for rt_mx_ji in rt_mx_ji_list:
          op = _canon_op(rt_mx_ji)
          edge = (j_seq, op)
          if edge in adjacency[i_seq]:
            continue
          adjacency[i_seq].add(edge)
          adjacency[j_seq].add((i_seq, _canon_op(rt_mx_ji.inverse())))
          added_edges += 1

    return added_edges

  def atoms_within_radius(self, center_atom, model, radius):
    """Return a boolean mask for atoms within *radius* of *center_atom*.

    Uses a KD-tree for O(n log n) neighbour lookup.  The tree is built
    once per *model* and reused across calls; this matters for
    structures with many seed atoms (e.g. 12 Fe sites in 8fum) where
    rebuilding the tree for every seed dominated the radius-search
    cost.

    Parameters
    ----------
    center_atom : iotbx.pdb.hierarchy.atom
    model : mmtbx.model.manager
    radius : float

    Returns
    -------
    cctbx.array_family.flex.bool
    """
    if getattr(self, '_kdtree_model_id', None) != id(model):
      self._kdtree = KDTree(model.get_sites_cart())
      self._kdtree_model_id = id(model)
      self._kdtree_n_atoms = model.get_number_of_atoms()
    indices = self._kdtree.query_ball_point(center_atom.xyz, radius)
    selected = flex.bool(self._kdtree_n_atoms, False)
    for idx in indices:
      selected[idx] = True
    return selected

  def atoms_within_radius_string(self, center_atom, model, radius):
    """Select atoms within *radius* using CCTBX selection syntax.

    Parameters
    ----------
    center_atom : iotbx.pdb.hierarchy.atom
    model : mmtbx.model.manager
    radius : float

    Returns
    -------
    cctbx.array_family.flex.bool
    """
    atom_group = center_atom.parent()
    residue_group = atom_group.parent()
    chain = residue_group.parent()

    center_terms = [
      f'chain {chain.id.strip()}',
      f'resseq {residue_group.resseq.strip()}',
      f'resname {atom_group.resname.strip()}',
      f'name {center_atom.name.strip()}',
    ]
    altloc = atom_group.altloc.strip()
    if altloc:
      center_terms.append(f'altloc "{altloc}"')

    center_selector = ' and '.join(center_terms)
    selection_str = f'within({radius}, ({center_selector}))'
    return model.selection(selection_str)

  def atoms_within_radius_best(self, center_atom, model, radius):
    """Return a boolean mask for atoms within *radius* of *center_atom*.

    Uses the cached KD-tree path; falls back to the CCTBX selection-string
    approach if the tree returns no atoms (e.g. the seed atom is outside
    the model's coordinate range).

    Parameters
    ----------
    center_atom : iotbx.pdb.hierarchy.atom
    model : mmtbx.model.manager
    radius : float

    Returns
    -------
    cctbx.array_family.flex.bool
    """
    mask = self.atoms_within_radius(center_atom, model, radius)
    if mask.count(True) == 0:
      mask = self.atoms_within_radius_string(center_atom, model, radius)
    return mask


# =============================================================================
# BondCutDetector  (backbone / sidechain cut heuristics)
# =============================================================================

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


# =============================================================================
# QMRegionGrower
# =============================================================================

class QMRegionGrower:
  """Grow a QM region by BFS traversal of a covalent graph.

  The grower stops BFS at "cuttable" bonds (sidechain C-C sp3 cuts and
  backbone CA-C / CA-N cuts) and records the far atom as a tentative cap.
  It does **not** actually place hydrogens; cap placement is performed by
  :class:`HydrogenCapper` at the calling site after sub-model selection.

  Parameters
  ----------
  bond_cut_detector : BondCutDetector
      Classifier for cuttable bonds.
  log : file-like or None, optional
      Destination for diagnostic messages.
  """

  def __init__(self, bond_cut_detector, log=None):
    self.bond_cut_detector = bond_cut_detector
    self.log = log

  # ------------------------------------------------------------------
  # Public interface
  # ------------------------------------------------------------------

  def grow_by_depth(self, seed_atoms, adjacency, model, max_depth=3):
    """Symmetry-aware BFS over ``(iseq, sym_op)`` nodes, recording cap
    candidates by pruning.

    BFS state (visited set, queue, cap_candidates) keys on nodes
    ``(iseq, sym_op)`` rather than bare iseqs. Each adjacency edge
    ``(j_seq, edge_op)`` composes with the current node's op to
    produce the neighbour node
    ``(j_seq, current_op.multiply(edge_op))``. Lattice translations
    and point-group symmetry are therefore handled uniformly.

    Whenever a cuttable bond ``current -> neighbour`` is encountered,
    *neighbour* is marked visited so BFS will not expand past it, and
    ``(neighbour, current)`` is stored as a tentative cap via
    :meth:`_try_mark_cap`.  Two situations promote a tentative cap back to
    interior:

    * **Re-encounter** -- *neighbour* is later reached from a different
      node.  The cap is then a regular interior node on at least two
      paths, so it is demoted via :meth:`_demote_cap_candidate`.
    * **Adjacent-cap conflict** -- the would-be cap is directly bonded to
      an existing cap candidate (detected at cap-creation time inside
      :meth:`_try_mark_cap`).  Both nodes are promoted to interior; this
      catches the edge case the re-encounter check cannot see, since
      cap candidates are never themselves enqueued.

    Together these guarantee that every surviving cap candidate has only
    its recorded anchor as a QM-region neighbour.

    Parameters
    ----------
    seed_atoms : set of int
        Seed i_seq values. Internally promoted to nodes carrying the
        identity ``rt_mx``.
    adjacency : collections.defaultdict of set
        Tagged adjacency: ``adj[i_seq]`` is a set of
        ``(j_seq, edge_op)`` tuples (see
        :meth:`AtomGraphBuilder.build_adjacency`). Read-only;
        not mutated.
    model : mmtbx.model.manager
    max_depth : int, optional
        Kept for API compatibility; currently unused.  Default is 3.

    Returns
    -------
    visited : set of (int, rt_mx)
        All nodes in the grown QM region (interior + surviving caps).
    cap_candidates : dict
        ``{(cap_iseq, cap_op): (anchor_iseq, anchor_op)}``; each cap is
        guaranteed to have only the recorded anchor as a QM-region
        neighbour.
    """
    atoms = model.get_hierarchy().atoms()
    identity = _canon_op(sgtbx.rt_mx())
    seed_nodes = {(iseq, identity) for iseq in seed_atoms}
    visited = set(seed_nodes)
    cap_candidates = {}
    queue = deque(seed_nodes)

    while queue:
      current = queue.popleft()
      curr_iseq, curr_op = current
      curr_resname = atoms[curr_iseq].parent().resname.strip().upper()

      for (nbr_iseq, edge_op) in list(adjacency[curr_iseq]):
        nbr_op = _canon_op(curr_op.multiply(edge_op))
        neighbour = (nbr_iseq, nbr_op)

        if neighbour in visited:
          if (neighbour in cap_candidates and
              cap_candidates[neighbour] != current):
            self._demote_cap_candidate(
              neighbour, cap_candidates, visited,
              seed_nodes, adjacency, atoms, queue,
              protected={current},
            )
          continue

        visited.add(neighbour)

        if self.bond_cut_detector.is_cc_single_sp3_bond(
            curr_resname, atoms[curr_iseq], atoms[nbr_iseq], adjacency):
          print('Found C-C single bond between:', file=self.log)
          print('  ' + atoms[curr_iseq].format_atom_record().rstrip(),
                file=self.log)
          print('  ' + atoms[nbr_iseq].format_atom_record().rstrip(),
                file=self.log)
          self._try_mark_cap(
            neighbour, current, cap_candidates, visited,
            seed_nodes, adjacency, atoms, queue,
          )
          continue

        if self.bond_cut_detector.is_ca_c_bond(
            atoms[curr_iseq], atoms[nbr_iseq], adjacency):
          if self._any_amide_of_current_residue_in_visited(
              curr_iseq, curr_op, visited, atoms):
            print('Found backbone CA-C bond between:', file=self.log)
            print('  ' + atoms[curr_iseq].format_atom_record().rstrip(),
                  file=self.log)
            print('  ' + atoms[nbr_iseq].format_atom_record().rstrip(),
                  file=self.log)
            self._try_mark_cap(
              neighbour, current, cap_candidates, visited,
              seed_nodes, adjacency, atoms, queue,
            )
            continue

        if self.bond_cut_detector.is_ca_n_bond(
            atoms[curr_iseq], atoms[nbr_iseq], adjacency):
          if self._any_amide_of_next_residue_in_visited(
              curr_iseq, curr_op, visited, atoms):
            print('Found backbone CA-N bond between:', file=self.log)
            print('  ' + atoms[curr_iseq].format_atom_record().rstrip(),
                  file=self.log)
            print('  ' + atoms[nbr_iseq].format_atom_record().rstrip(),
                  file=self.log)
            self._try_mark_cap(
              neighbour, current, cap_candidates, visited,
              seed_nodes, adjacency, atoms, queue,
            )
            continue

        queue.append(neighbour)

    print(
      f'Found {len(cap_candidates)} candidate atoms for hydrogen capping '
      f'based on heuristics.',
      file=self.log,
    )
    return visited, cap_candidates

  def _demote_cap_candidate(self, cap, cap_candidates, visited,
                            seed_nodes, adjacency, atoms, queue,
                            protected=set()):
    """Promote *cap* (a node) back to interior and re-open BFS through it.

    Removes *cap* from ``cap_candidates`` and re-enqueues it so the BFS
    can reach the neighbours that were previously hidden behind the cut.
    Every other covalent-neighbour node of *cap* is discarded from
    ``visited`` (and from ``cap_candidates``) so it can be rediscovered
    through the now-interior cap.  Three sets of nodes are protected
    from this discard: seeds, the cap's original anchor (a confirmed
    interior node whose subtree was already explored), and any
    nodes in *protected*.

    Parameters
    ----------
    cap : (int, rt_mx)
        Node of the cap candidate being demoted.
    cap_candidates : dict
        ``{cap_node: anchor_node}``; modified in-place.
    visited : set of (int, rt_mx)
        Modified in-place.
    seed_nodes : set of (int, rt_mx)
        Seed nodes; always protected from discard.
    adjacency : collections.defaultdict of set
        Read-only; tagged ``{i_seq: {(j_seq, edge_op), ...}}``.
    atoms : cctbx.array_family.flex array
        Atom objects indexed by i_seq, used for log messages.
    queue : collections.deque
        BFS queue; *cap* is appended for further expansion.
    protected : set of (int, rt_mx), optional
        Additional nodes to skip when discarding.
    """
    cap_iseq, cap_op = cap
    print('Demoting cap candidate:', file=self.log)
    print('  ' + atoms[cap_iseq].format_atom_record().rstrip(), file=self.log)
    original_anchor = cap_candidates.pop(cap)
    protected = set(protected) | {original_anchor}
    for (nb_iseq, edge_op) in adjacency[cap_iseq]:
      nb_node = (nb_iseq, _canon_op(cap_op.multiply(edge_op)))
      if nb_node in protected or nb_node in seed_nodes:
        continue
      visited.discard(nb_node)
      cap_candidates.pop(nb_node, None)
    queue.append(cap)

  def _try_mark_cap(self, candidate, anchor, cap_candidates, visited,
                     seed_nodes, adjacency, atoms, queue):
    """Record *candidate* (a node) as a cap with *anchor* (a node),
    unless that would create an adjacent-cap conflict.

    A conflict arises if *candidate* shares a covalent bond with an
    existing cap candidate: capping both would leave each with a
    non-anchor QM-region neighbour and place two Hs at chemically
    nonsensical positions.  In that case the conflicting cap(s) are
    demoted to interior and *candidate* is also promoted to interior
    (enqueued for further BFS expansion).

    Parameters
    ----------
    candidate : (int, rt_mx)
        Node of the would-be cap (the atom on the far side of the
        cuttable bond).
    anchor : (int, rt_mx)
        Node of the would-be anchor (the atom on the QM-region side
        of the cuttable bond).
    cap_candidates : dict
        ``{cap_node: anchor_node}``; modified in-place.
    visited : set of (int, rt_mx)
        Modified in-place via :meth:`_demote_cap_candidate` if a conflict
        is found.
    seed_nodes : set of (int, rt_mx)
        Passed through to demotion.
    adjacency : collections.defaultdict of set
        Read-only.
    atoms : cctbx.array_family.flex array
        Atom objects indexed by i_seq, used for log messages.
    queue : collections.deque
        BFS queue; *candidate* is appended to it if promoted to interior.
    """
    cand_iseq, cand_op = candidate
    conflicting = []
    for (nb_iseq, edge_op) in adjacency[cand_iseq]:
      nb_node = (nb_iseq, _canon_op(cand_op.multiply(edge_op)))
      if nb_node == anchor:
        continue
      if nb_node in cap_candidates:
        conflicting.append(nb_node)
    if conflicting:
      print('Adjacent-cap conflict; promoting to interior:', file=self.log)
      print('  ' + atoms[cand_iseq].format_atom_record().rstrip(),
            file=self.log)
      for c in conflicting:
        self._demote_cap_candidate(
          c, cap_candidates, visited, seed_nodes, adjacency,
          atoms, queue, protected={anchor, candidate},
        )
      queue.append(candidate)
      return
    cap_candidates[candidate] = anchor

  def grow_by_distance(self, start_atom_index, adjacency, model,
                       max_distance=5.0):
    """Expand from *start_atom_index* by BFS within *max_distance*.

    Only neighbours whose distance from the start atom is within
    *max_distance* are added to the visited set and enqueued.

    Parameters
    ----------
    start_atom_index : int
    adjacency : collections.defaultdict of set
    model : mmtbx.model.manager
    max_distance : float, optional
        Default is 5.0 Angstrom.

    Returns
    -------
    set of int
        i_seqs of atoms in the grown region.
    """
    atoms = model.get_hierarchy().atoms()
    start_atom = atoms[start_atom_index]
    visited = {start_atom_index}
    queue = deque([start_atom_index])

    while queue:
      curr_idx = queue.popleft()
      for neighbor in adjacency[curr_idx]:
        if neighbor in visited:
          continue
        if start_atom.distance(atoms[neighbor]) > max_distance:
          continue
        visited.add(neighbor)
        queue.append(neighbor)

    return visited

  # ------------------------------------------------------------------
  # Private helpers - amide-group checks
  # ------------------------------------------------------------------

  def _any_amide_in_residue_group_in_visited(self, residue_group, sym_op,
                                             visited):
    """Return ``True`` if any backbone amide atom of *residue_group*, under
    the given ``sym_op`` symmetry image, is in *visited*.

    "Amide atoms" here means the backbone N, C, and O atoms across all
    alternate conformations of the residue group. Each amide atom is
    looked up in *visited* as a node ``(i_seq, sym_op)`` -- so two
    different symmetry images of the same residue are treated as
    distinct.

    Parameters
    ----------
    residue_group : iotbx.pdb.hierarchy.residue_group
    sym_op : cctbx.sgtbx.rt_mx
        Symmetry image of the residue to check.
    visited : set of (int, rt_mx)

    Returns
    -------
    bool
    """
    amide_atom_names = {'N', 'C', 'O'}
    amide_i_seqs = [
      residue_atom.i_seq
      for atom_group in residue_group.atom_groups()
      for residue_atom in atom_group.atoms()
      if residue_atom.name.strip().upper() in amide_atom_names
    ]
    if not amide_i_seqs:
      return False
    return any((i_seq, sym_op) in visited for i_seq in amide_i_seqs)

  def _any_amide_of_current_residue_in_visited(self, atom_idx, sym_op,
                                               visited, atoms):
    """Return ``True`` if any amide atom of *atom_idx*'s residue, under
    *sym_op*, is in *visited*.

    Parameters
    ----------
    atom_idx : int
    sym_op : cctbx.sgtbx.rt_mx
        Symmetry image of the residue to check.
    visited : set of (int, rt_mx)
    atoms : flex array of iotbx.pdb.hierarchy.atom

    Returns
    -------
    bool
    """
    residue_group = atoms[atom_idx].parent().parent()
    return self._any_amide_in_residue_group_in_visited(
      residue_group, sym_op, visited)

  def _any_amide_of_next_residue_in_visited(self, atom_idx, sym_op,
                                            visited, atoms):
    """Return ``True`` if any amide atom of the next residue, under
    *sym_op*, is in *visited*.

    Uses positional order in the chain's ``residue_groups()`` list rather
    than arithmetic on resseq, so insertion codes and non-contiguous
    numbering are handled correctly.

    Parameters
    ----------
    atom_idx : int
    sym_op : cctbx.sgtbx.rt_mx
        Symmetry image of the residue to check.
    visited : set of (int, rt_mx)
    atoms : flex array of iotbx.pdb.hierarchy.atom

    Returns
    -------
    bool
    """
    residue_group = atoms[atom_idx].parent().parent()
    chain = residue_group.parent()

    rg_list = list(chain.residue_groups())
    try:
      current_pos = rg_list.index(residue_group)
    except ValueError:
      return False

    if current_pos + 1 >= len(rg_list):
      return False

    next_residue_group = rg_list[current_pos + 1]
    return self._any_amide_in_residue_group_in_visited(
      next_residue_group, sym_op, visited)


# =============================================================================
# HydrogenCapper
# =============================================================================

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


# =============================================================================
# ChargeEstimator
# =============================================================================

class ChargeEstimator:
  """Estimate the net amino-acid sidechain charge of a model.

  Walks the model's residues, looks up each standard amino acid's expected
  charge from ``CHARGED_SIDECHAINS`` / ``STANDARD_RESIDUE_CHARGES`` based on
  which protonation-site H atoms are actually present, and (optionally)
  adds a contribution for free peptide termini.  The model may be a full
  structure, a truncated QM region, or any subset thereof.

  Parameters
  ----------
  include_terminal_charges : bool, optional
      If ``True``, detect free peptide termini and add their charges.
      Default is ``False``.
  n_terminus_charge : int, optional
      Charge assigned to each detected free N-terminus.  Default is ``1``.
  c_terminus_charge : int, optional
      Charge assigned to each detected free C-terminus.  Default is ``-1``.
  log : file-like or None, optional
      Destination for diagnostic messages.
  """

  def __init__(self, include_terminal_charges=False,
               n_terminus_charge=1, c_terminus_charge=-1, log=None):
    self.include_terminal_charges = include_terminal_charges
    self.n_terminus_charge = n_terminus_charge
    self.c_terminus_charge = c_terminus_charge
    self.log = log

  # ------------------------------------------------------------------
  # Public interface
  # ------------------------------------------------------------------

  def calculate(self, model):
    """Estimate the net amino-acid sidechain charge of *model*.

    Parameters
    ----------
    model : mmtbx.model.manager

    Returns
    -------
    dict
        Keys: ``total_charge``, ``sidechain_charge``, ``terminal_charge``,
        ``contributors``, ``residue_counts``, ``unknown_residues``,
        ``residue_count_total``.
    """
    residue_data = self._collect_residue_data(model.get_hierarchy().atoms())

    total_sidechain_charge = 0
    total_terminal_charge = 0
    residue_counts = defaultdict(int)
    unknown_residues = set()
    residue_contributions = {}

    for key, data in residue_data.items():
      resname = data['resname']
      if resname not in STANDARD_RESIDUE_CHARGES:
        unknown_residues.add(resname)
        continue
      residue_counts[resname] += 1
      sidechain_charge = self._sidechain_charge(resname, data['atom_names'])
      total_sidechain_charge += sidechain_charge
      residue_contributions[key] = {
        'resname': resname,
        'sidechain_charge': sidechain_charge,
        'terminal_charge': 0,
      }

    if self.include_terminal_charges:
      total_terminal_charge = self._add_terminal_charges(
        residue_data, residue_contributions, model
      )

    contributors = self._build_contributors_list(residue_contributions)
    total_charge = total_sidechain_charge + total_terminal_charge

    return {
      'total_charge': total_charge,
      'sidechain_charge': total_sidechain_charge,
      'terminal_charge': total_terminal_charge,
      'contributors': contributors,
      'residue_counts': dict(sorted(residue_counts.items())),
      'unknown_residues': sorted(unknown_residues),
      'residue_count_total': len(residue_data),
    }

  def show(self, seed_index, charge_summary, out=None):
    """Write a human-readable charge summary to *out*.

    No-op when *out* is ``None``.

    Parameters
    ----------
    seed_index : int
    charge_summary : dict
        As returned by :meth:`calculate`.
    out : file-like or None, optional
    """
    if out is None:
      return
    fmt = self._fmt_signed
    print(
      f'Estimated amino-acid net charge (seed {seed_index}): '
      f'{fmt(charge_summary["total_charge"])}',
      file=out,
    )
    print(f'  sidechain contribution: {fmt(charge_summary["sidechain_charge"])}',
          file=out)
    if self.include_terminal_charges:
      print(f'  terminal contribution : {fmt(charge_summary["terminal_charge"])}',
            file=out)
    print(f'  residues counted      : {charge_summary["residue_count_total"]}',
          file=out)
    if charge_summary['residue_counts']:
      counts_str = ', '.join(
        f'{resname}:{count}'
        for resname, count in charge_summary['residue_counts'].items()
      )
      print(f'  residue composition   : {counts_str}', file=out)
    if charge_summary['unknown_residues']:
      print(
        '  skipped non-standard residues: ' +
        ', '.join(charge_summary['unknown_residues']),
        file=out,
      )
    if charge_summary['contributors']:
      print('  charge-contributing residues:', file=out)
      for entry in charge_summary['contributors']:
        residue_id = f'chain {entry["chain_id"]} resseq {entry["resseq"]}'
        if entry['icode']:
          residue_id += f' icode {entry["icode"]}'
        if entry['altloc']:
          residue_id += f' altloc {entry["altloc"]}'
        print(
          f'    {entry["resname"]} {residue_id}: '
          f'sidechain: {fmt(entry["sidechain_charge"])}, '
          f'terminal: {fmt(entry["terminal_charge"])}; '
          f'total: {fmt(entry["total_charge"])}',
          file=out,
        )

  # ------------------------------------------------------------------
  # Private helpers
  # ------------------------------------------------------------------

  @staticmethod
  def _fmt_signed(x):
    return f'{x:+}' if x else '0'

  @staticmethod
  def _collect_residue_data(atoms):
    """Group atoms by residue and collect atom names and backbone i_seqs.

    Parameters
    ----------
    atoms : flex array of iotbx.pdb.hierarchy.atom

    Returns
    -------
    dict
        Keyed by ``(chain_id, resseq, icode, altloc)`` tuples.
    """
    residue_data = {}
    for atom in atoms:
      i_seq = atom.i_seq
      atom_group = atom.parent()
      residue_group = atom_group.parent()
      chain = residue_group.parent()
      resname = atom_group.resname.strip().upper()

      key = (
        chain.id.strip(),
        residue_group.resseq.strip(),
        residue_group.icode.strip(),
        atom_group.altloc.strip(),
      )
      if key not in residue_data:
        residue_data[key] = {
          'resname': resname,
          'atom_names': set(),
          'n_i_seqs': set(),
          'c_i_seqs': set(),
        }

      atom_name = atom.name.strip().upper()
      residue_data[key]['atom_names'].add(atom_name)
      if atom_name == 'N':
        residue_data[key]['n_i_seqs'].add(i_seq)
      elif atom_name == 'C':
        residue_data[key]['c_i_seqs'].add(i_seq)

    return residue_data

  @staticmethod
  def _is_charged_sidechain_present(resname, atom_names):
    """Return ``True`` if the charge-bearing heavy atom(s) of *resname* are
    in *atom_names*.

    Residues without a ``charged_heavy_atoms`` entry in
    ``CHARGED_SIDECHAINS`` are considered always present (return ``True``);
    otherwise the residue's sidechain is treated as truncated/absent when
    none of its charged heavy atoms appear.

    Parameters
    ----------
    resname : str
        Three-letter residue name (upper-case).
    atom_names : set of str
        Upper-case atom names observed in the residue.

    Returns
    -------
    bool
    """
    charged_atoms = CHARGED_SIDECHAINS.get(resname, {}).get(
      'charged_heavy_atoms', set()
    )
    if not charged_atoms:
      return True
    return bool(charged_atoms.intersection(atom_names))

  @staticmethod
  def _calculate_side_chain_charge(resname, atom_names):
    """Compute the sidechain net charge from the present protonation-site H atoms.

    Starts from ``formal_charge`` (the residue's charge in its CCD
    reference protonation state) and adjusts for missing or extra
    protonation-site hydrogens.  Deuterium names (``D...``) are folded
    onto their canonical hydrogen names so D/H labelling does not affect
    the count.

    Parameters
    ----------
    resname : str
        Three-letter residue name (upper-case); must be a key in
        ``CHARGED_SIDECHAINS``.
    atom_names : set of str
        Upper-case atom names observed in the residue.

    Returns
    -------
    int
        Net sidechain charge.
    """
    entry = CHARGED_SIDECHAINS.get(resname, {})
    formal_charge = entry.get('formal_charge', 0)
    protonation_sites = entry.get('protonation_sites', set())

    def _canonical_h_site(name):
      return 'H' + name[1:] if name.startswith('D') else name

    sites = {_canonical_h_site(n) for n in protonation_sites}
    present = {_canonical_h_site(n) for n in atom_names}
    return formal_charge - len(sites) + len(sites.intersection(present))

  def _sidechain_charge(self, resname, atom_names):
    """Return the net sidechain charge for *resname* given its present atoms.

    Returns ``0`` if the residue's charge-bearing heavy atoms are absent
    (truncated sidechain).  Otherwise, residues listed in
    ``CHARGED_SIDECHAINS`` get a protonation-aware charge; remaining
    residues fall back to ``STANDARD_RESIDUE_CHARGES``.

    Parameters
    ----------
    resname : str
        Three-letter residue name (upper-case).
    atom_names : set of str
        Upper-case atom names observed in the residue.

    Returns
    -------
    int
    """
    if not self._is_charged_sidechain_present(resname, atom_names):
      return 0
    if resname in CHARGED_SIDECHAINS:
      return self._calculate_side_chain_charge(resname, atom_names)
    return STANDARD_RESIDUE_CHARGES.get(resname, 0)

  def _add_terminal_charges(self, residue_data, residue_contributions, model):
    """Detect free termini and accumulate terminal charges.

    A residue is treated as having a free N-terminus if it is the first
    residue group in its chain (within *model*), and likewise a free
    C-terminus if it is the last.

    Parameters
    ----------
    residue_data : dict
    residue_contributions : dict
        Modified in-place.
    model : mmtbx.model.manager

    Returns
    -------
    int
        Total terminal charge contribution.
    """
    first_rgs, last_rgs = self._chain_terminus_keys(model)
    total = 0
    for key, data in residue_data.items():
      if key not in residue_contributions:
        continue
      rg_key = key[:3]  # (chain_id, resseq, icode); altloc is dropped

      if data['n_i_seqs'] and rg_key in first_rgs:
        residue_contributions[key]['terminal_charge'] += self.n_terminus_charge
        total += self.n_terminus_charge
      if data['c_i_seqs'] and rg_key in last_rgs:
        residue_contributions[key]['terminal_charge'] += self.c_terminus_charge
        total += self.c_terminus_charge

    return total

  @staticmethod
  def _chain_terminus_keys(model):
    """Return ``(first, last)`` sets of residue-group keys for chain termini.

    Each key is ``(chain_id, resseq, icode)`` for the first / last residue
    group in its chain within *model*.

    Parameters
    ----------
    model : mmtbx.model.manager

    Returns
    -------
    first : set of tuple
    last : set of tuple
    """
    first = set()
    last = set()
    for chain in model.get_hierarchy().chains():
      rgs = list(chain.residue_groups())
      if not rgs:
        continue
      chain_id = chain.id.strip()
      first.add((chain_id, rgs[0].resseq.strip(), rgs[0].icode.strip()))
      last.add((chain_id, rgs[-1].resseq.strip(), rgs[-1].icode.strip()))
    return first, last

  @staticmethod
  def _build_contributors_list(residue_contributions):
    """Flatten *residue_contributions* into a sorted list of nonzero entries.

    Residues whose total charge (sidechain + terminal) is zero are
    omitted from the output.

    Parameters
    ----------
    residue_contributions : dict
        Keyed by ``(chain_id, resseq, icode, altloc)``; values carry
        ``resname``, ``sidechain_charge``, and ``terminal_charge``.

    Returns
    -------
    list of dict
        One entry per nonzero-charge residue, sorted by key.  Each entry
        has keys ``chain_id``, ``resseq``, ``icode``, ``altloc``,
        ``resname``, ``sidechain_charge``, ``terminal_charge``,
        ``total_charge``.
    """
    contributors = []
    for key, contribution in sorted(residue_contributions.items()):
      chain_id, resseq, icode, altloc = key
      sc = contribution['sidechain_charge']
      tc = contribution['terminal_charge']
      total = sc + tc
      if abs(total) == 0:
        continue
      contributors.append({
        'chain_id': chain_id,
        'resseq': resseq,
        'icode': icode,
        'altloc': altloc,
        'resname': contribution['resname'],
        'sidechain_charge': sc,
        'terminal_charge': tc,
        'total_charge': total,
      })
    return contributors


# =============================================================================
# Program  (ProgramTemplate entry point)
# =============================================================================

class Program(ProgramTemplate):
  """Extract QM regions with BFS expansion and hydrogen capping.

  Seeds the QM region either from all metals in the structure (default) or
  from a user-supplied CCTBX selection string (``selection`` parameter).
  """

  description = '''
  Grows a QM region around each seed site by BFS, optionally caps dangling
  bonds with hydrogen atoms, estimates the net charge, and writes a PDB
  file, an mmCIF file, and a sidecar PHIL file per seed.  Seeds are all
  metals in the structure unless a custom selection string is provided
  via the ``selection`` parameter.
  '''

  datatypes = ['model', 'phil']
  master_phil_str = master_phil_str

  # ------------------------------------------------------------------
  # ProgramTemplate interface
  # ------------------------------------------------------------------

  def validate(self):
    """Validate user inputs before :meth:`run` is called."""
    if not self.data_manager.has_models():
      raise Sorry('No model provided. Please supply a PDB or mmCIF file.')

    model = self.data_manager.get_model()
    selection_strings = [s for s in (self.params.selection or []) if s]
    for sel_str in selection_strings:
      try:
        model.selection(sel_str)
      except Exception as e:
        raise Sorry(
          f"Invalid selection string '{sel_str}': {e}"
        )

  def run(self):
    """Locate seeds, grow QM regions, and write a PDB, mmCIF, and sidecar
    PHIL file per seed."""
    model = self.data_manager.get_model()

    # model.add_hydrogens(
    #   element='H',
    #   neutron=True,
    #   occupancy=1.0,
    # )

    model.process(
      pdb_interpretation_params=model.get_current_pdb_interpretation_params(),
      make_restraints=True,
    )

    model = self._apply_altloc_filter(model)

    self._graph_builder = AtomGraphBuilder()
    self._capper = HydrogenCapper(log=self.logger)
    self._bond_cut_detector = BondCutDetector(
      use_preferred_cuts=self.params.use_preferred_cuts,
      log=self.logger,
    )
    self._region_grower = QMRegionGrower(
      self._bond_cut_detector,
      log=self.logger,
    )
    self._charge_estimator = ChargeEstimator(
      include_terminal_charges=self.params.include_terminal_charges,
      n_terminus_charge=self.params.n_terminus_charge,
      c_terminus_charge=self.params.c_terminus_charge,
      log=self.logger,
    )
    self._results = []

    seed_finder = SeedFinder()
    selection_strings = [s for s in (self.params.selection or []) if s]
    element_filter = [
      e for e in (self.params.metal_element or []) if e and e.strip()
    ]
    if element_filter and selection_strings:
      print(
        f"Note: ignoring metal_element={element_filter} because an "
        f"explicit selection was provided.",
        file=self.logger,
      )
    seed_groups = seed_finder.find(
      model,
      selection_strings=selection_strings,
      element_filter=element_filter,
    )

    if not seed_groups:
      raise Sorry('No seed atoms found in the model.')

    if selection_strings:
      for sel_str, seeds in seed_groups:
        if not seeds:
          raise Sorry(
            f"Selection '{sel_str}' matched no atoms in the model."
          )
        print(
          f"Selection '{sel_str}' matched {len(seeds)} atom(s):",
          file=self.logger,
        )
        self._print_seeds(seeds, label='selected atoms')
    else:
      seeds_all = [s for _, grp in seed_groups for s in grp]
      self._print_seeds(seeds_all, label='seed atoms')

    seeds_flat = [atom for _, grp in seed_groups for atom in grp]

    grm = model.get_restraints_manager().geometry
    sites_cart = model.get_sites_cart()
    bond_proxies_simple, bond_proxies_asu = grm.get_all_bond_proxies(
      sites_cart=sites_cart
    )
    asu_mappings = grm.pair_proxies(
      sites_cart=sites_cart).bond_proxies.asu_mappings()

    adjacency = self._graph_builder.build_adjacency(
      bond_proxies_simple, bond_proxies_asu, asu_mappings)

    cutoff = self.params.metal_ligand_cutoff
    if not self.params.skip_radius_search:
      added_edges = self._graph_builder.add_seed_contact_edges(
        seeds_flat, model, adjacency, cutoff=cutoff
      )
      print(
        f'Added {added_edges} distance-based seed-contact edges '
        f'(cutoff={cutoff:.2f} A)',
        file=self.logger,
      )
    else:
      print(
        'Skipping distance-based seed-contact edges '
        '(skip_radius_search=True).',
        file=self.logger,
      )
    print(
      f'Always-included seed-centered radius: {self.params.radius:.2f} A',
      file=self.logger,
    )

    for seed_index, (sel_str, seeds) in enumerate(seed_groups, start=1):
      result = self._process_seed(
        seed_index, seeds, model, adjacency, selection_str=sel_str
      )
      self._results.append(result)

  def get_results(self):
    """Return per-seed result dicts produced during :meth:`run`.

    Returns
    -------
    list of dict
        Each dict contains:

        file_name : str
            Output filename stem (no extension); the PDB and mmCIF copies
            are written as ``file_name + '.pdb'`` and ``file_name + '.cif'``.
        n_atoms : int
            Number of atoms in the QM region.
        charge_summary : dict
            As returned by :meth:`ChargeEstimator.calculate`.
        model : mmtbx.model.manager
            Truncated sub-model with caps applied and the parent's
            restraints manager attached. Suitable for direct in-memory
            consumption by downstream tools (e.g., ``mmtbx.qmi``) without
            the disk round-trip through the written PDB/mmCIF.
        seed_iseqs : list of int
            Sorted 0-based positional indices of the seed atoms inside
            ``model``.
        cap_iseqs : list of int
            Sorted 0-based positional indices of the cap atoms inside
            ``model`` (empty when capping is disabled).
        selection_string : str or None
            The CCTBX selection string that produced this seed group, or
            ``None`` for metal-scan groups.
    """
    return self._results

  # ------------------------------------------------------------------
  # Private helpers
  # ------------------------------------------------------------------

  def _apply_altloc_filter(self, model):
    """Return *model* with non-selected altlocs removed per residue.

    Behaviour is controlled by ``params.altloc``:

    * ``'all'`` -- no filtering; *model* is returned unchanged.
    * ``'auto'`` -- for each residue group containing non-blank altlocs,
      retain the altloc with the highest mean atom occupancy and drop
      the others.
    * any specific letter (e.g. ``'A'``) -- retain that letter where it
      is present; if a residue has non-blank altlocs but not the
      requested letter, fall back to the highest-occupancy altloc and
      emit a warning.

    Atoms with empty altloc (the shared backbone) are always kept.

    Parameters
    ----------
    model : mmtbx.model.manager

    Returns
    -------
    mmtbx.model.manager
        Either *model* itself (no atoms to drop) or a freshly selected
        sub-model whose restraints have been re-indexed by
        :meth:`mmtbx.model.manager.select`.
    """
    from libtbx import Auto
    altloc_param = self.params.altloc
    # PHIL turns the unquoted token `auto` into the libtbx.Auto
    # sentinel object, not the string "auto"; treat both as auto mode.
    if altloc_param is Auto or altloc_param is None:
      mode_auto = True
      letter = None
    else:
      stripped = altloc_param.strip()
      if stripped.lower() == 'all':
        return model
      if stripped.lower() == 'auto':
        mode_auto = True
        letter = None
      else:
        mode_auto = False
        letter = stripped

    keep = flex.bool(model.get_number_of_atoms(), True)
    drop_count = 0

    for chain in model.get_hierarchy().chains():
      for rg in chain.residue_groups():
        alt_ags = [
          ag for ag in rg.atom_groups() if ag.altloc.strip() != ''
        ]
        if not alt_ags:
          continue

        if mode_auto:
          chosen = max(alt_ags, key=self._mean_occupancy)
        else:
          matching = [
            ag for ag in alt_ags if ag.altloc.strip() == letter
          ]
          if matching:
            chosen = matching[0]
          else:
            chosen = max(alt_ags, key=self._mean_occupancy)
            print(
              f'Warning: altloc "{letter}" not found in '
              f'chain {chain.id.strip()} resseq {rg.resseq.strip()}; '
              f'falling back to altloc "{chosen.altloc.strip()}" '
              f'(highest mean occupancy).',
              file=self.logger,
            )

        for ag in alt_ags:
          if ag is chosen:
            continue
          for atom in ag.atoms():
            keep[atom.i_seq] = False
            drop_count += 1

    if drop_count == 0:
      return model

    filtered = model.select(keep)
    mode_str = 'auto' if mode_auto else letter
    print(
      f'Altloc filter (altloc={mode_str}): dropped {drop_count} '
      f'atom(s); {filtered.get_number_of_atoms()} atoms retained.',
      file=self.logger,
    )
    return filtered

  @staticmethod
  def _mean_occupancy(atom_group):
    """Return the mean atom occupancy of *atom_group*.

    Returns ``0.0`` for an empty atom group so :func:`max` never sees a
    missing value.

    Parameters
    ----------
    atom_group : iotbx.pdb.hierarchy.atom_group

    Returns
    -------
    float
    """
    atoms = atom_group.atoms()
    n = len(atoms)
    if n == 0:
      return 0.0
    return sum(a.occ for a in atoms) / n

  def _print_seeds(self, seeds, label='seeds'):
    """Print a labelled list of seed atoms to ``self.logger``.

    Parameters
    ----------
    seeds : list of iotbx.pdb.hierarchy.atom
    label : str, optional
        Human-readable label.  Default is ``'seeds'``.
    """
    print(f'Found {label}:', file=self.logger)
    for idx, atom in enumerate(seeds, start=1):
      print(f'  {idx}: {atom.format_atom_record().rstrip()}',
            file=self.logger)

  def _process_seed(self, seed_index, seeds, model, adjacency,
                    selection_str=None):
    """Run the full pipeline for one seed site.

    Parameters
    ----------
    seed_index : int
        1-based index used in output filenames and messages.
    seeds : list of iotbx.pdb.hierarchy.atom
        Seed atom(s).
    model : mmtbx.model.manager
    adjacency : collections.defaultdict of set
        Full covalent graph (not modified).
    selection_str : str or None, optional
        Original selection string used to derive the output filename.

    Returns
    -------
    dict
        Keys: ``file_name``, ``n_atoms``, ``charge_summary``, ``model``,
        ``seed_iseqs``, ``cap_iseqs``, ``selection_string``.
    """
    qm_atoms = self._seed_qm_region(seeds, model)
    visited_nodes, cap_nodes = self._region_grower.grow_by_depth(
      qm_atoms, adjacency, model, max_depth=self.params.max_depth
    )

    visited_nodes = self._add_hull_waters(model, visited_nodes)

    (model_sel, seed_indices, cap_indices,
     cap_original_elements) = self._materialize_qm_region(
      model, visited_nodes, cap_nodes, seeds)

    charge_summary = self._charge_estimator.calculate(model=model_sel)
    self._charge_estimator.show(seed_index, charge_summary, out=self.logger)

    file_name = self._make_output_filename(
      seed_index, seeds, selection_str=selection_str
    )
    if self.params.write_files:
      self._write_submodel(
        model_sel, model.crystal_symmetry(), file_name)
      self._write_sidecar(
        file_name, cap_indices, cap_original_elements,
        seed_indices, selection_str)

    return {
      'file_name': file_name,
      'n_atoms': model_sel.get_number_of_atoms(),
      'charge_summary': charge_summary,
      # In-memory hand-off for downstream consumers (e.g., mmtbx.qmi
      # predictor): the truncated sub-model with restraints attached,
      # and the positional indices of seed and cap atoms inside it.
      'model': model_sel,
      'seed_iseqs': seed_indices,
      'cap_iseqs': cap_indices,
      'cap_original_elements': cap_original_elements,
      'selection_string': selection_str,
    }

  # ------------------------------------------------------------------
  # Pipeline steps
  # ------------------------------------------------------------------

  def _seed_qm_region(self, seeds, model):
    """Return the initial i_seq set for the QM region.

    When *skip_radius_search* is False (default) all atoms within
    ``params.radius`` of every seed are included.  When *skip_radius_search*
    is True only the seed atoms themselves are added and BFS expansion
    (controlled by ``params.max_depth``) is relied upon to grow the region.

    Parameters
    ----------
    seeds : list of iotbx.pdb.hierarchy.atom
    model : mmtbx.model.manager

    Returns
    -------
    set of int
    """
    qm_atoms = set()
    if self.params.skip_radius_search:
      for seed in seeds:
        qm_atoms.add(seed.i_seq)
    else:
      for seed in seeds:
        mask = self._graph_builder.atoms_within_radius_best(
          seed, model, self.params.radius
        )
        qm_atoms |= set(mask.iselection())
    return qm_atoms

  def _add_hull_waters(self, model, visited_nodes):
    """Extend *visited_nodes* with water residue groups whose
    representative oxygen sits inside the convex hull of the visited
    set.

    The hull is built from the *materialized* Cartesian positions of
    the visited nodes, not from ASU coordinates -- so symmetry
    images that already participate in the QM region contribute to
    the bounding volume on equal footing with the parent atoms.

    For every water residue group in the model and every space-group
    operator, the operator is composed with the integer lattice shift
    that brings the water's representative atom closest to the hull
    centroid; if the resulting position lies inside the hull, nodes
    for every atom in that residue group (under that combined op)
    are added.  No-op when
    ``params.include_waters_in_convex_hull`` is ``False``.

    Parameters
    ----------
    model : mmtbx.model.manager
    visited_nodes : set of (int, sgtbx.rt_mx)

    Returns
    -------
    set of (int, sgtbx.rt_mx)
        Union of *visited_nodes* with the new water nodes.
    """
    if not self.params.include_waters_in_convex_hull:
      return visited_nodes

    cs = model.crystal_symmetry()
    uc = cs.unit_cell()
    sg = cs.space_group()
    sites_frac = uc.fractionalize(model.get_sites_cart())

    # Hull built from materialized node positions
    hull_cart = flex.vec3_double()
    for (iseq, op) in visited_nodes:
      r = matrix.sqr(op.r().as_double())
      t = matrix.col(op.t().as_double())
      f = matrix.col(sites_frac[iseq])
      hull_cart.append(uc.orthogonalize((r * f + t).elems))

    if hull_cart.size() < 4:
      return visited_nodes
    try:
      hull = ConvexHull(hull_cart)
    except Exception:
      return visited_nodes

    centroid_frac = uc.fractionalize(hull_cart.mean())

    visited_iseqs_by_op = {}
    for (iseq, op) in visited_nodes:
      visited_iseqs_by_op.setdefault(op.as_xyz(), set()).add(iseq)

    # Collect water residues once per *parent model*: ``_add_hull_waters``
    # is called per seed but the water set is the same every call, and on
    # large structures (e.g. ~12k waters across 12 seeds) the list-comp
    # over alt-confed atoms dominated the per-call cost.
    if (getattr(self, '_waters_cache', None) is None
        or self._waters_cache_model_id != id(model)):
      cached = []
      for chain in model.get_hierarchy().chains():
        for rg in chain.residue_groups():
          if common_residue_names_get_class(
              rg.atom_groups()[0].resname) != 'common_water':
            continue
          rep_iseq = rg.atom_groups()[0].atoms()[0].i_seq
          iseqs = [a.i_seq for ag in rg.atom_groups()
                   for a in ag.atoms()]
          cached.append((rep_iseq, iseqs))
      self._waters_cache = cached
      self._waters_cache_model_id = id(model)
      self._waters_cache_rep_frac = (
        sites_frac.select(flex.size_t(w[0] for w in cached))
        if cached else None)
    waters = self._waters_cache
    if not waters:
      return visited_nodes

    n_waters = len(waters)
    rep_frac_orig = self._waters_cache_rep_frac
    cx, cy, cz = centroid_frac

    centroid_tuple = (cx, cy, cz)
    new_nodes = set()
    for sg_op in sg.all_ops():
      # Bulk apply r * f + t over all water rep positions.  ``v * R``
      # in flex is row-vector * matrix; using ``R.transpose()`` makes
      # the result equivalent to the column-vector form ``R * v``.
      r_t = matrix.sqr(sg_op.r().as_double()).transpose()
      t_tuple = sg_op.t().as_double()
      rep_frac_sym = rep_frac_orig * r_t + t_tuple

      # Nearest-image shift toward centroid, vectorized per axis.
      diff = rep_frac_sym - centroid_tuple
      xs_d, ys_d, zs_d = diff.parts()
      shifts_x = xs_d.iround()
      shifts_y = ys_d.iround()
      shifts_z = zs_d.iround()
      shifts_vec = flex.vec3_double(
        shifts_x.as_double(), shifts_y.as_double(), shifts_z.as_double())
      rep_cart = uc.orthogonalize(rep_frac_sym - shifts_vec)

      # Batched hull test: one flex dot product per hull plane over all
      # waters, executed in C.
      inside = flex.bool(n_waters, True)
      for eq in hull.equations:
        d = rep_cart.dot((float(eq[0]), float(eq[1]), float(eq[2]))) \
            + float(eq[3])
        inside &= (d <= 1e-6)

      for k in inside.iselection():
        sx, sy, sz = shifts_x[k], shifts_y[k], shifts_z[k]
        shift_op = sgtbx.rt_mx(f"x{-sx:+d},y{-sy:+d},z{-sz:+d}")
        final_op = _canon_op(shift_op.multiply(sg_op))
        already = visited_iseqs_by_op.get(final_op.as_xyz(), set())
        for iseq in waters[k][1]:
          if iseq not in already:
            new_nodes.add((iseq, final_op))

    if new_nodes:
      print(
        f'Adding {len(new_nodes)} water atom nodes inside convex '
        f'hull to QM region.',
        file=self.logger,
      )
    return visited_nodes | new_nodes

  def _materialize_qm_region(self, model, visited_nodes, cap_nodes,
                              seeds):
    """Build a materialized QM sub-model from the BFS nodes.

    Nodes carrying the same ``i_seq`` but different ``sym_op`` are
    realised as separate atoms in the output: the parent atom's
    position is transformed by the ``sym_op``'s fractional rotation
    plus translation, and the result lives in its own chain block.
    Atoms whose materialized positions coincide (special-position
    case: multiple symmetry operations producing the same physical
    point) are deduplicated by Cartesian-position key.

    Chain IDs are assigned per ``(original_chain_id, sym_op)`` pair so
    each symmetry image is distinguishable in the output.

    Parameters
    ----------
    model : mmtbx.model.manager
        The full ASU model with restraints attached.
    visited_nodes : set of (int, sgtbx.rt_mx)
        BFS visited set; each tuple is ``(parent_iseq, sym_op)``.
    cap_nodes : dict
        ``{(cap_iseq, cap_op): (anchor_iseq, anchor_op)}``.
    seeds : list of iotbx.pdb.hierarchy.atom
        Seed atoms (always present at the identity image).

    Returns
    -------
    model_sel : mmtbx.model.manager
        Sub-model containing the materialized atoms.
    seed_indices : list of int
        Positional indices of seed atoms in *model_sel*.
    cap_indices : list of int
        Positional indices of cap atoms in *model_sel*.
    """
    import iotbx.pdb
    import mmtbx.model as mmtbx_model

    parent_hier = model.get_hierarchy()
    parent_atoms = parent_hier.atoms()
    cs = model.crystal_symmetry()
    unit_cell = cs.unit_cell()
    identity_xyz = sgtbx.rt_mx().as_xyz()

    # Group iseqs by sym_op (keyed by canonical xyz string)
    by_op = defaultdict(set)
    for (iseq, op) in visited_nodes:
      by_op[op.as_xyz()].add(iseq)
    op_keys = sorted(by_op.keys(),
                     key=lambda k: (k != identity_xyz, k))

    # Build chain ID map: (original_chain_id, sym_op_xyz) -> new_chain_id.
    # Stick to single-character chain IDs because PDB's chain field is
    # one character wide; multi-char IDs overflow into the residue-name
    # field and break the bond graph that ``make_restraints`` builds.
    orig_chain_ids = sorted({ch.id.strip() for ch in parent_hier.chains()})
    available_chain_ids = iter(c for c in
      "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
      if c not in orig_chain_ids)
    chain_id_for = {}
    for op_idx, op_xyz in enumerate(op_keys):
      for orig_id in orig_chain_ids:
        if op_xyz == identity_xyz:
          chain_id_for[(orig_id, op_xyz)] = orig_id
        else:
          try:
            chain_id_for[(orig_id, op_xyz)] = next(available_chain_ids)
          except StopIteration:
            raise RuntimeError(
              "exhausted single-character chain IDs while "
              "materialising symmetry images")

    # Build an empty target hierarchy and append per-sym_op sub-models
    out_root = iotbx.pdb.hierarchy.root()
    out_hier_model = iotbx.pdb.hierarchy.model(id="1")
    out_root.append_model(out_hier_model)

    # As we copy each per-sym_op sub-hierarchy into the master,
    # remember which parent node each fresh atom came from.
    # ``model.select(mask)`` preserves atom order, so the k-th atom in
    # the copied sub-hierarchy corresponds to the k-th parent i_seq in
    # ascending order from ``by_op[op_xyz]``.
    node_of_atom = {}  # output atom -> (parent_iseq, op_xyz)
    n_parent_atoms = parent_atoms.size()
    for op_xyz in op_keys:
      op = sgtbx.rt_mx(op_xyz)
      r_mat = matrix.sqr(op.r().as_double())
      t_vec = matrix.col(op.t().as_double())

      iseqs_in_order = sorted(by_op[op_xyz])
      mask = flex.bool(n_parent_atoms, False)
      for i in iseqs_in_order:
        mask[i] = True
      sub_model = model.select(mask)
      sub_hier = sub_model.get_hierarchy()

      # Transform coordinates from ASU positions to the image.
      if op_xyz != identity_xyz:
        new_sites = flex.vec3_double()
        for xyz in sub_model.get_sites_cart():
          frac = matrix.col(unit_cell.fractionalize(xyz))
          new_frac = r_mat * frac + t_vec
          new_sites.append(unit_cell.orthogonalize(new_frac.elems))
        sub_model.set_sites_cart(new_sites)
        # Rename chains so this image is distinguishable
        for ch in sub_hier.chains():
          ch.id = chain_id_for[(ch.id.strip(), op_xyz)]

      # Detached-copy each chain into the master and record the parent
      # node for every fresh atom object.
      k = 0
      for ch in sub_hier.chains():
        new_ch = ch.detached_copy()
        out_hier_model.append_chain(new_ch)
        for new_atom in new_ch.atoms():
          node_of_atom[new_atom] = (iseqs_in_order[k], op_xyz)
          k += 1

    # Single pass: deduplicate atoms occupying the same physical
    # position (special-position case: multiple sym_ops produce the
    # same point), drop atom groups / residue groups / chains left
    # empty by that removal, and index the survivors by parent node
    # so the cap pass below can find them.  Dedup is O(n^2) within
    # each (element, name) bucket; the 0.2 A tolerance absorbs the
    # fractionalize -> rotate/translate -> orthogonalize drift.
    POS_TOL = 0.2
    seen_by_key = defaultdict(list)  # (element, name) -> [atom, ...]
    atom_for_node = {}  # (parent_iseq, op_xyz) -> surviving atom
    for ch in list(out_hier_model.chains()):
      for rg in list(ch.residue_groups()):
        for ag in list(rg.atom_groups()):
          for atom in list(ag.atoms()):
            ename = (atom.element.strip(), atom.name.strip())
            if any(atom.distance(prev) < POS_TOL
                   for prev in seen_by_key[ename]):
              ag.remove_atom(atom)
            else:
              seen_by_key[ename].append(atom)
              atom_for_node[node_of_atom[atom]] = atom
          if len(list(ag.atoms())) == 0:
            rg.remove_atom_group(ag)
        if len(list(rg.atom_groups())) == 0:
          ch.remove_residue_group(rg)
      if len(list(ch.residue_groups())) == 0:
        out_hier_model.remove_chain(ch)
    out_root.atoms().reset_serial()

    def _lookup_node(iseq, op_xyz):
      return atom_for_node.get((iseq, op_xyz))

    # Apply hydrogen capping on the materialized atoms.  Remember the
    # cap atom's original element so downstream consumers that rebuild
    # restraints can temporarily restore it (pdb_interpretation matches
    # name + element against the monomer library; the H we stamp here
    # would clash with the library's expected element for that name
    # inside non-standard residues).
    orig_element_for_iseq = {}
    if self.params.do_capping:
      for cap_node, anchor_node in cap_nodes.items():
        cap_iseq_orig, cap_op = cap_node
        anchor_iseq_orig, anchor_op = anchor_node
        cap_atom = _lookup_node(cap_iseq_orig, cap_op.as_xyz())
        anchor_atom = _lookup_node(anchor_iseq_orig, anchor_op.as_xyz())
        if cap_atom is not None:
          orig_element_for_iseq[cap_atom] = cap_atom.element.strip()
        self._capper.cap_atom(anchor_atom, cap_atom)

    # Assemble the final mmtbx.model.manager
    model_sel = mmtbx_model.manager(
      model_input=None,
      pdb_hierarchy=out_root,
      crystal_symmetry=cs)

    # ``atom_by_key`` already holds the materialized atoms; their
    # ``.i_seq`` is the post-cap value after ``reset_serial`` above and
    # is not changed by ``mmtbx_model.manager`` wrapping the hierarchy.
    seed_indices = sorted(
      a.i_seq for a in (_lookup_node(s.i_seq, identity_xyz)
                        for s in seeds) if a is not None)
    cap_atoms_sorted = sorted(
      (a for a in (_lookup_node(iseq, op.as_xyz())
                   for (iseq, op) in cap_nodes) if a is not None),
      key=lambda a: a.i_seq)
    cap_indices = [a.i_seq for a in cap_atoms_sorted]
    cap_original_elements = [orig_element_for_iseq.get(a, 'C')
                             for a in cap_atoms_sorted]

    return model_sel, seed_indices, cap_indices, cap_original_elements

  def _write_submodel(self, model_sel, crystal_symmetry, file_name):
    """Write *model_sel* as both a PDB and an mmCIF file.

    mmCIF is emitted alongside the PDB because mmCIF stores the element in
    ``_atom_site.type_symbol`` as a first-class field, so capped atoms keep
    their hydrogen identity through downstream re-parsing instead of being
    silently re-classified from their atom names.

    Parameters
    ----------
    model_sel : mmtbx.model.manager
    crystal_symmetry : cctbx.crystal.symmetry
    file_name : str
        Filename stem (no extension); ``.pdb`` and ``.cif`` are appended here.
    """
    hierarchy = model_sel.get_hierarchy()
    hierarchy.write_pdb_file(
      crystal_symmetry=crystal_symmetry,
      file_name=file_name + '.pdb',
    )
    hierarchy.write_mmcif_file(
      crystal_symmetry=crystal_symmetry,
      file_name=file_name + '.cif',
      data_block_name='qm_region',
    )
    print(
      f'Wrote QM region to {file_name}.pdb and {file_name}.cif '
      f'({model_sel.get_number_of_atoms()} atoms)',
      file=self.logger,
    )

  def _write_sidecar(self, file_name, cap_indices, cap_original_elements,
                     seed_indices, selection_string):
    """Write the per-region sidecar PHIL file as ``<file_name>.phil``.

    The sidecar carries metadata needed by downstream tools: which atoms
    are hydrogen caps (so the QM driver can freeze them), which atoms
    seeded the region (so the charge/spin classifier knows the metal or
    user-selected anchor), and the original selection string for
    provenance.

    Parameters
    ----------
    file_name : str
        Filename stem (no extension); ``.phil`` is appended here.
    cap_indices : list of int
        Sorted 0-based positional indices of cap atoms in the output file.
    cap_original_elements : list of str
        Element symbol each cap atom carried before being replaced by H
        (parallel to *cap_indices*); needed by downstream consumers that
        rebuild restraints on the sub-model.
    seed_indices : list of int
        Sorted 0-based positional indices of seed atoms in the output file.
    selection_string : str or None
        The CCTBX selection string that seeded the region, or ``None`` for
        metal-scan regions.
    """
    sidecar_phil = libtbx.phil.parse(master_sidecar_phil_str)
    sidecar_phil_extract = sidecar_phil.extract()
    sidecar_phil_extract.endoexo_region.cap_atoms = list(cap_indices)
    sidecar_phil_extract.endoexo_region.cap_original_elements = list(
      cap_original_elements)
    sidecar_phil_extract.endoexo_region.seed_atoms = list(seed_indices)
    sidecar_phil_extract.endoexo_region.selection_string = selection_string
    sidecar_phil_path = file_name + '.phil'
    with open(sidecar_phil_path, 'w') as file:
      sidecar_phil.format(sidecar_phil_extract).show(out=file)
    print(f'Wrote sidecar to {sidecar_phil_path}', file=self.logger)

  def _make_output_filename(self, seed_index, seeds, selection_str=None):
    """Build the output filename stem from the model name and seed identity.

    When *selection_str* is given the stem is::

        {model_stem}_sel{seed_index:03d}_within{radius}A_depth{depth}

    When seeds are individual metal atoms the stem is::

        {model_stem}_{element}_chain{chain}_res{resseq}_within{radius}A_depth{depth}

    Parameters
    ----------
    seed_index : int
    seeds : list of iotbx.pdb.hierarchy.atom
    selection_str : str or None, optional

    Returns
    -------
    str
        Filename stem without extension.
    """
    model_name = self.data_manager.get_default_model_name()
    model_stem = os.path.splitext(os.path.basename(model_name))[0]

    common_suffix = (
      f'_within{self.params.radius:.2f}A'
      f'_depth{self.params.max_depth}'
    )

    if selection_str:
      seed_suffix = f'_sel{seed_index:03d}' + common_suffix
    else:
      seed = seeds[0]
      atom_group = seed.parent()
      residue_group = atom_group.parent()
      chain = residue_group.parent()
      element = seed.element.strip().upper()
      chain_id = chain.id.strip()
      resseq = residue_group.resseq.strip()
      seed_suffix = (
        f'_{element}_chain{chain_id}_res{resseq}' + common_suffix
      )

    return self.get_default_output_filename(
      prefix=model_stem,
      suffix=seed_suffix,
    )

