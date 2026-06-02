"""Breadth-first QM-region growth over the covalent graph, with cap-candidate
detection."""

from __future__ import absolute_import, division, print_function

from collections import deque

from cctbx import sgtbx

from mmtbx.geometry_restraints.endoexo.util import _canon_op


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
