"""Breadth-first QM-region growth over the covalent graph, with cap-candidate
detection."""

from __future__ import absolute_import, division, print_function

from collections import deque

from cctbx import sgtbx

from mmtbx.geometry_restraints.endo_exo.util import _canon_op, _neighbour_iseqs
from mmtbx.geometry_restraints.endo_exo.cutting import (
  BondCutDetector, PREFERRED_CUTS)


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
    # Used only by the preferred-cut fallback: when a preferred cut is no
    # longer achievable we re-cut inward with the pure geometric C-C
    # heuristic, regardless of the residue's PREFERRED_CUTS entry.
    self._geom_detector = BondCutDetector(use_preferred_cuts=False, log=log)

  # ------------------------------------------------------------------
  # Public interface
  # ------------------------------------------------------------------

  def grow_by_depth(self, seed_atoms, adjacency, model, max_depth=3,
                    forced_cut_bonds=None, geometric_for_overgrowth=False):
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

    * **Re-encounter**: *neighbour* is later reached from a different
      node.  The cap is then a regular interior node on at least two
      paths, so it is demoted via :meth:`_demote_cap_candidate`.
    * **Adjacent-cap conflict**: the would-be cap is directly bonded to
      an existing cap candidate (detected at cap-creation time inside
      :meth:`_try_mark_cap`).  Both nodes are promoted to interior; this
      catches the edge case the re-encounter check cannot see, since
      cap candidates are never themselves enqueued.

    Together these guarantee that every surviving cap candidate has only
    its recorded anchor as a QM-region neighbour.

    Parameters
    ----------
    seed_atoms : set of int or (int, rt_mx)
        Seed nodes.  Bare i_seqs are promoted to nodes carrying the
        identity ``rt_mx``; ``(iseq, op)`` tuples (symmetry-image seeds
        from the symmetry-aware radius search) are used as-is.
    adjacency : collections.defaultdict of set
        Tagged adjacency: ``adj[i_seq]`` is a set of
        ``(j_seq, edge_op)`` tuples (see
        :meth:`AtomGraphBuilder.build_adjacency`). Read-only;
        not mutated.
    model : mmtbx.model.manager
    max_depth : int, optional
        Kept for API compatibility; currently unused.  Default is 3.
    forced_cut_bonds : set of frozenset, optional
        Bonds (each a ``frozenset({iseq_a, iseq_b})``) to treat as
        cuttable regardless of residue type or geometry.  Matched on the
        bare iseq pair, so every symmetry image of the bond is cut.  Used
        by :meth:`grow_region` (preferred-cut fallback) to re-cut residues
        whose preferred cut became unachievable.
    geometric_for_overgrowth : bool, optional
        When ``True``, residues with no seed (radius) atom (i.e. those
        reached purely as BFS backbone overgrowth) defer to the
        geometric C-C sp3 heuristic instead of their ``PREFERRED_CUTS``
        entry.  The heuristic cuts at the first sp3 C-C bond along the
        BFS path (typically CA-CB), trimming such residues to backbone.
        Residues carrying a seed atom keep their preferred cuts.

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
    forced_cut_bonds = forced_cut_bonds or set()
    identity = _canon_op(sgtbx.rt_mx())
    # Seeds may be bare iseqs (identity image) or pre-built (iseq, op)
    # nodes (symmetry images from the symmetry-aware radius search);
    # normalise both to nodes.
    seed_nodes = {s if isinstance(s, tuple) else (s, identity)
                  for s in seed_atoms}
    # iseqs belonging to a residue that contains at least one seed
    # (radius) atom.  Residues absent from this set are pure BFS
    # overgrowth and, when geometric_for_overgrowth is on, defer to the
    # geometric C-C heuristic.  Empty (and unused) otherwise.
    #
    # APPROXIMATION: keyed on bare iseq, so the sym_op is dropped.  With
    # symmetry-aware seeding the same residue can appear under several ops,
    # and "has a seed" is really a per-(residue, op) property.  Here, if
    # ANY image of a residue carries a seed, every image is treated as
    # seeded -> uses preferred cuts. This errs toward the preferred cut
    # (the pre-overgrowth-rule behaviour), so it can only over-include, not
    # break. It bites only when one image of a residue coordinates (seeded)
    # while another image of the SAME residue is pure overgrowth in the same
    # sphere (rare). To make it exact, key on (residue, op) instead.
    residue_has_seed_iseqs = set()
    if geometric_for_overgrowth:
      for (s_iseq, _s_op) in seed_nodes:
        residue_group = atoms[s_iseq].parent().parent()
        for atom_group in residue_group.atom_groups():
          for residue_atom in atom_group.atoms():
            residue_has_seed_iseqs.add(residue_atom.i_seq)
    visited = set(seed_nodes)
    cap_candidates = {}
    queue = deque(seed_nodes)

    while queue:
      current = queue.popleft()
      curr_iseq, curr_op = current
      curr_resname = atoms[curr_iseq].parent().resname.strip().upper()
      # Overgrowth residues (no seed atom) defer to the geometric C-C
      # heuristic; seeded residues keep their configured (preferred) cuts.
      if geometric_for_overgrowth and curr_iseq not in residue_has_seed_iseqs:
        curr_cut_detector = self._geom_detector
      else:
        curr_cut_detector = self.bond_cut_detector

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

        if frozenset((curr_iseq, nbr_iseq)) in forced_cut_bonds:
          print('Forced cut (preferred-cut fallback) between:', file=self.log)
          print('  ' + atoms[curr_iseq].format_atom_record().rstrip(),
                file=self.log)
          print('  ' + atoms[nbr_iseq].format_atom_record().rstrip(),
                file=self.log)
          self._try_mark_cap(
            neighbour, current, cap_candidates, visited,
            seed_nodes, adjacency, atoms, queue,
          )
          continue

        if curr_cut_detector.is_cc_single_sp3_bond(
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

  def grow_region(self, seed_atoms, adjacency, model, max_depth=3,
                  preferred_cut_fallback=False):
    """Grow the QM region, optionally re-cutting residues whose preferred
    cut became unachievable.

    Single public entry point.  With ``preferred_cut_fallback=False`` this
    is exactly one :meth:`grow_by_depth` pass.  With it ``True`` two things
    happen, both of which defer to the geometric C-C sp3 heuristic when the
    preferred cut is not the right boundary:

    * **Overgrowth residues**: residues with no seed (radius) atom use
      the geometric heuristic from the start (see
      ``geometric_for_overgrowth`` on :meth:`grow_by_depth`), so a residue
      reached only as backbone overgrowth is cut at its first sp3 C-C bond
      (typically CA-CB) rather than at its preferred site.
    * **Consumed preferred cuts**: the pass is repeated, accumulating
      geometric fallback cuts, until no new ones are found (below).

    A preferred cut is "consumed" when both of its endpoints end up
    interior to the region, typically because the radius search seeded
    atoms on both sides of it.  The canonical cut can then no longer be
    made, and BFS walks inward through the (now interior) sidechain into
    the backbone.  For each such residue this computes the outermost
    geometric C-C cut inward of the consumed bond, adds it to a
    ``forced_cut_bonds`` set, and re-runs :meth:`grow_by_depth`.  The trim
    follows from re-running: an atom survives only if it is still reachable
    from a seed without crossing a cut, so the backbone overgrowth (and
    anything dragged in through it) falls away on its own.  Seeds
    (including all radius atoms) are never trimmed, so the radius floor is
    preserved.

    The loop terminates because each round either adds a genuinely new
    forced cut or returns none: once a residue's fallback cut is already
    in ``forced_cut_bonds`` it is filtered out and the loop stops, even
    though the preferred bond stays nominally consumed.  In practice it is
    two passes (discover, then apply-and-confirm).

    Parameters
    ----------
    seed_atoms : set of int
    adjacency : collections.defaultdict of set
    model : mmtbx.model.manager
    max_depth : int, optional
    preferred_cut_fallback : bool, optional
        When ``False`` (default) a single :meth:`grow_by_depth` pass is
        returned unchanged.

    Returns
    -------
    visited : set of (int, rt_mx)
    cap_candidates : dict
        Same contract as :meth:`grow_by_depth`.
    """
    forced_cut_bonds = set()
    while True:
      visited, cap_candidates = self.grow_by_depth(
        seed_atoms, adjacency, model, max_depth=max_depth,
        forced_cut_bonds=forced_cut_bonds,
        geometric_for_overgrowth=preferred_cut_fallback)
      if not preferred_cut_fallback:
        return visited, cap_candidates
      new_cuts = self._fallback_cuts_for_consumed_preferred(
        visited, cap_candidates, adjacency, model, existing=forced_cut_bonds)
      if not new_cuts:
        return visited, cap_candidates
      forced_cut_bonds |= new_cuts

  def _fallback_cuts_for_consumed_preferred(self, visited, cap_candidates,
                                            adjacency, model, existing):
    """Return geometric fallback cuts for residues whose preferred cut is
    consumed (both endpoints interior and not actually cut).

    Parameters
    ----------
    visited : set of (int, rt_mx)
    cap_candidates : dict
        ``{cap_node: anchor_node}`` from the last BFS pass.
    adjacency : collections.defaultdict of set
    model : mmtbx.model.manager
    existing : set of frozenset
        Forced cuts already applied; results in this set are filtered out
        (this is what lets :meth:`grow_region` terminate).

    Returns
    -------
    set of frozenset
        New ``frozenset({iseq_a, iseq_b})`` bonds to force-cut.
    """
    atoms = model.get_hierarchy().atoms()
    new_cuts = set()
    handled = set()
    for (iseq, op) in visited:
      atom = atoms[iseq]
      preferred = PREFERRED_CUTS.get(atom.parent().resname.strip().upper())
      if not preferred:
        continue
      name = atom.name.strip().upper()
      if name not in preferred:
        continue
      partner = self._sibling_in_residue(atom, (preferred - {name}).pop())
      if partner is None or (partner.i_seq, op) not in visited:
        continue
      bond = frozenset((iseq, partner.i_seq))
      if bond in handled:
        continue
      handled.add(bond)
      # If one endpoint is the cap of the other, the preferred cut actually
      # fired: the bond is the boundary, not consumed.  Leave it.
      if (cap_candidates.get((iseq, op)) == (partner.i_seq, op) or
          cap_candidates.get((partner.i_seq, op)) == (iseq, op)):
        continue
      ca_iseq, dist = self._residue_ca_distances(iseq, adjacency, atoms)
      if ca_iseq is None:
        continue
      d_self, d_partner = dist.get(iseq), dist.get(partner.i_seq)
      if d_self is None or d_partner is None:
        continue
      # Cut from the backbone-side endpoint (covalently closer to CA),
      # keeping the tip (functional-group) side.
      if d_self <= d_partner:
        cut_iseq, keep_iseq = iseq, partner.i_seq
      else:
        cut_iseq, keep_iseq = partner.i_seq, iseq
      fallback = self._first_geom_cut_inward(
        cut_iseq, keep_iseq, dist, adjacency, atoms)
      if fallback is not None and fallback not in existing:
        new_cuts.add(fallback)
    return new_cuts

  def _first_geom_cut_inward(self, cut_iseq, keep_iseq, dist, adjacency,
                             atoms):
    """Walk from *cut_iseq* toward CA (away from *keep_iseq*) and return
    the first bond that passes the geometric C-C heuristic.

    Parameters
    ----------
    cut_iseq : int
        Backbone-side endpoint of the consumed preferred bond.
    keep_iseq : int
        Tip-side endpoint (the step direction is away from this).
    dist : dict
        ``{iseq: covalent-distance-to-CA}`` within the residue, from
        :meth:`_residue_ca_distances`.
    adjacency : collections.defaultdict of set
    atoms : flex array of iotbx.pdb.hierarchy.atom

    Returns
    -------
    frozenset or None
        ``frozenset({iseq_a, iseq_b})`` of the first geometric cut, or
        ``None`` if CA is reached without one (the existing backbone cut
        rules then apply on the next BFS pass).
    """
    cur, prev = cut_iseq, keep_iseq
    while True:
      # Heavy neighbour strictly closer to CA, within the residue, not
      # where we came from.  min() picks the closest if a branch forks.
      candidates = [
        nb for nb in _neighbour_iseqs(adjacency, cur)
        if nb != prev and nb in dist and dist[nb] < dist[cur]
        and not atoms[nb].element_is_hydrogen()]
      if not candidates:
        return None
      nxt = min(candidates, key=lambda n: dist[n])
      if self._geom_detector.is_cc_single_sp3_bond(
          None, atoms[cur], atoms[nxt], adjacency):
        return frozenset((cur, nxt))
      prev, cur = cur, nxt

  @staticmethod
  def _sibling_in_residue(atom, name):
    """Return the atom named *name* in *atom*'s residue group, or None."""
    for ag in atom.parent().parent().atom_groups():
      for sibling in ag.atoms():
        if sibling.name.strip().upper() == name:
          return sibling
    return None

  @staticmethod
  def _residue_ca_distances(iseq, adjacency, atoms):
    """Return ``(ca_iseq, {iseq: hops-to-CA})`` over the covalent graph
    restricted to the residue group of *iseq*.

    Returns ``(None, {})`` if the residue has no CA atom.

    Parameters
    ----------
    iseq : int
    adjacency : collections.defaultdict of set
    atoms : flex array of iotbx.pdb.hierarchy.atom

    Returns
    -------
    (int or None, dict)
    """
    rg = atoms[iseq].parent().parent()
    rg_iseqs = {a.i_seq for ag in rg.atom_groups() for a in ag.atoms()}
    ca_iseq = None
    for ag in rg.atom_groups():
      for a in ag.atoms():
        if a.name.strip().upper() == 'CA':
          ca_iseq = a.i_seq
          break
      if ca_iseq is not None:
        break
    if ca_iseq is None:
      return None, {}
    dist = {ca_iseq: 0}
    queue = deque([ca_iseq])
    while queue:
      cur = queue.popleft()
      for nb in _neighbour_iseqs(adjacency, cur):
        if nb in dist or nb not in rg_iseqs:
          continue
        dist[nb] = dist[cur] + 1
        queue.append(nb)
    return ca_iseq, dist

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
    looked up in *visited* as a node ``(i_seq, sym_op)``, so two
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
