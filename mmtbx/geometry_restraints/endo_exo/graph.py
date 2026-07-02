"""Symmetry-aware covalent adjacency graph construction for the endo_exo
QM-region builder."""

from __future__ import absolute_import, division, print_function

from collections import defaultdict

from cctbx import sgtbx
from cctbx.array_family import flex
from scipy.spatial import KDTree

from mmtbx.geometry_restraints.endo_exo.util import _canon_op


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
    # that are sym-equivalent under the site symmetry of the seed.
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
    once per *model* and reused across calls.

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

  def seed_sym_nodes_within_radius(self, seeds, model, radius):
    """Return symmetry-image ``(iseq, op)`` nodes within *radius* of a seed.

    The KD-tree radius search (:meth:`atoms_within_radius_best`) only sees
    ASU atoms, so a metal on or near a special position misses the
    symmetry-related copies of its coordinating residues that are
    physically inside the buffer sphere.  This enumerates those copies via
    ``pair_asu_table`` -- the same mechanism :meth:`add_seed_contact_edges`
    uses for coordination edges -- and returns one ``(j_seq, rt_mx_ji)``
    node per atom-image within *radius* of any seed.

    Only **non-identity** images are returned; the identity image is
    already covered by the KD-tree search, so the caller unions the two.

    Parameters
    ----------
    seeds : list of iotbx.pdb.hierarchy.atom
        Seed atom objects.
    model : mmtbx.model.manager
    radius : float
        Distance threshold in Angstrom.

    Returns
    -------
    set of (int, sgtbx.rt_mx)
        Non-identity atom-image nodes within *radius* of a seed.
    """
    identity_xyz = _canon_op(sgtbx.rt_mx()).as_xyz()
    xs = model.get_xray_structure()
    pat = xs.pair_asu_table(distance_cutoff=radius)
    # ``all_interactions_from_inside_asu=True`` enumerates every close
    # contact of an ASU atom rather than collapsing site-symmetry-equivalent
    # pairs (see add_seed_contact_edges).
    sym_table = pat.extract_pair_sym_table(
      skip_j_seq_less_than_i_seq=False,
      all_interactions_from_inside_asu=True)
    nodes = set()
    for seed in seeds:
      for j_seq, rt_mx_ji_list in dict(sym_table[seed.i_seq]).items():
        for rt_mx_ji in rt_mx_ji_list:
          op = _canon_op(rt_mx_ji)
          if op.as_xyz() == identity_xyz:
            continue
          nodes.add((j_seq, op))
    return nodes
