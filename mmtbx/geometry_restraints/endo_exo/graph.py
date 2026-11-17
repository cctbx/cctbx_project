"""Symmetry-aware covalent adjacency graph construction for the endo_exo
QM-region builder."""

from collections import defaultdict

from cctbx import sgtbx
from cctbx.crystal import super_cell


def sites_within_radius(sites_cart, centre, radius):
  """Return the sites within *radius* of *centre*, as a selection.

  Parameters
  ----------
  sites_cart : flex.vec3_double
  centre : tuple of float
      Cartesian coordinates.
  radius : float
      In Angstrom.

  Returns
  -------
  flex.bool
  """
  return (sites_cart - centre).dot() <= radius * radius


def nodes_within_radius(sites_cart, crystal_symmetry, radius, selection,
                        symmetry_tables=None):
  """Return every site, in place or as a symmetry copy, within *radius* of
  a selected site.

  Site symmetry is taken from the positions alone: a site whose copy lies
  within 0.5 A of it counts as on the symmetry element, whatever its
  occupancy.

  Parameters
  ----------
  sites_cart : flex.vec3_double
  crystal_symmetry : cctbx.crystal.symmetry
  radius : float
      In Angstrom.
  selection : list of int
      The sites to search around.
  symmetry_tables : cctbx.crystal.super_cell.sym_equiv_tables, optional
      From :func:`cctbx.crystal.super_cell.get_sym_equiv_tables` for
      *sites_cart* and *radius*, to search several selections without
      rebuilding them.

  Returns
  -------
  set of (int, sgtbx.rt_mx)
      ``(j_seq, rt_mx)`` nodes, the identity for a site in place.
  """
  identity = sgtbx.rt_mx()
  nodes = set()
  for i_seq in selection:
    for j_seq in sites_within_radius(
        sites_cart, sites_cart[i_seq], radius).iselection():
      nodes.add((j_seq, identity))
  siiu, _ = super_cell.get_siiu(
    crystal_symmetry=crystal_symmetry, sites_cart=sites_cart,
    selection=selection, select_within_radius=radius, buffer=0,
    symmetry_tables=symmetry_tables)
  for j_seq, rt_mx_ji_list in siiu.items():
    for rt_mx_ji in rt_mx_ji_list:
      nodes.add((j_seq, rt_mx_ji))
  return nodes


class AtomGraphBuilder:
  """Build a unified, symmetry-aware covalent adjacency graph from CCTBX
  restraint proxies.

  Each edge is tagged with the ``rt_mx_ji`` operation that brings atom
  *j* into contact with atom *i*. Intra-ASU bonds carry the identity
  operation; symmetry-crossing bonds carry the corresponding
  ``rt_mx_ji`` from ``bond_proxies_asu``.  Traversal therefore visits
  ``(i_seq, rt_mx)`` nodes rather than bare atoms, composing each node's
  operation with the edge operation as it moves along a bond, and lattice
  translations and point-group symmetry are handled uniformly.

  The proxies come from a restraints manager built with metal linking
  disabled, so metal-donor coordination bonds are absent from the graph
  and the region does not grow through a metal centre.
  """

  def __init__(self):
    self._sym_tables = None
    self._sym_tables_model = None
    self._sym_tables_radius = None

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
    identity = sgtbx.rt_mx()
    adj = defaultdict(set)
    for proxy in bond_proxies_simple:
      i_seq, j_seq = proxy.i_seqs
      adj[i_seq].add((j_seq, identity))
      adj[j_seq].add((i_seq, identity))
    for proxy in bond_proxies_asu:
      i_seq, j_seq = proxy.i_seq, proxy.j_seq
      # rt_mx hashes its denominators, and these operations and the
      # products of multiply() come on reduced ones, so every operation in a
      # node is put on the denominators (1, 12) that rt_mx() and the
      # space-group operations use.
      sym_op = asu_mappings.get_rt_mx_ji(proxy).new_denominators(1, 12)
      adj[i_seq].add((j_seq, sym_op))
      adj[j_seq].add((i_seq, sym_op.inverse()))
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
        if atom.element_is_hydrogen():
          continue
        for rt_mx_ji in rt_mx_ji_list:
          op = rt_mx_ji.new_denominators(1, 12)
          edge = (j_seq, op)
          if edge in adjacency[i_seq]:
            continue
          adjacency[i_seq].add(edge)
          adjacency[j_seq].add(
            (i_seq, rt_mx_ji.inverse().new_denominators(1, 12)))
          added_edges += 1

    return added_edges

  def _sym_equiv_tables(self, model, radius):
    """Return the whole-structure symmetry tables for *radius*.

    They cover every close contact in the structure, so they depend on the
    model and the radius but not on which seed is being grown.  They are
    built once and reused across seeds.  Site symmetry is taken from the
    positions alone, as in :func:`nodes_within_radius`.

    Parameters
    ----------
    model : mmtbx.model.manager
    radius : float

    Returns
    -------
    cctbx.crystal.super_cell.sym_equiv_tables
    """
    if (self._sym_tables_model is not model
        or self._sym_tables_radius != radius):
      self._sym_tables = super_cell.get_sym_equiv_tables(
        sites_cart=model.get_sites_cart(),
        crystal_symmetry=model.crystal_symmetry(), radius=radius)
      self._sym_tables_model = model
      self._sym_tables_radius = radius
    return self._sym_tables

  def seed_nodes_within_radius(self, seeds, model, radius):
    """Return the ``(iseq, op)`` nodes within *radius* of a seed.

    Atoms in place and their symmetry copies both count, so a seed on or
    near a special position gets the copies of its neighbouring residues
    that lie inside the buffer sphere.

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
        Atom nodes within *radius* of a seed, the identity for an atom in
        place.
    """
    return nodes_within_radius(
      sites_cart=model.get_sites_cart(),
      crystal_symmetry=model.crystal_symmetry(), radius=radius,
      selection=[seed.i_seq for seed in seeds],
      symmetry_tables=self._sym_equiv_tables(model, radius))
