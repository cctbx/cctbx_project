"""QMRegionBuilder -- the per-seed QM-region extraction pipeline.

Orchestrates seed discovery, covalent-graph construction, BFS region growth,
hydrogen capping, charge estimation, and (optionally) writing a PDB, mmCIF,
and sidecar PHIL file per seed.  Decoupled from the ProgramTemplate / data
manager so it can be driven directly in memory: construct with a params
object, then call :meth:`run` with a model.
"""

from __future__ import absolute_import, division, print_function

import os
from collections import defaultdict

import libtbx.phil
from libtbx import Auto
from libtbx.utils import Sorry
from cctbx import sgtbx
from cctbx.array_family import flex
from iotbx.pdb import common_residue_names_get_class
from scitbx import matrix

from scipy.spatial import ConvexHull

from mmtbx.geometry_restraints.endoexo.util import _canon_op
from mmtbx.geometry_restraints.endoexo.seeds import SeedFinder
from mmtbx.geometry_restraints.endoexo.graph import AtomGraphBuilder
from mmtbx.geometry_restraints.endoexo.cutting import BondCutDetector
from mmtbx.geometry_restraints.endoexo.grow import QMRegionGrower
from mmtbx.geometry_restraints.endoexo.capping import HydrogenCapper
from mmtbx.geometry_restraints.endoexo.charge import ChargeEstimator


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


class QMRegionBuilder(object):
  """Build QM regions around seed sites for a single model.

  Parameters
  ----------
  params : group_args-like
      Resolved parameters (the ``master_phil_str`` scope extract from
      :class:`mmtbx.programs.endoexo.Program`, or any object exposing the
      same attributes).
  logger : file-like or None, optional
      Destination for diagnostic messages.
  """

  def __init__(self, params, logger=None):
    self.params = params
    self.logger = logger
    self._results = []
    self._include_nodes = set()
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

  def run(self, model, model_name, default_output_filename):
    """Locate seeds, grow QM regions, and (optionally) write outputs.

    Parameters
    ----------
    model : mmtbx.model.manager
        The input model.  Restraints are processed in place.
    model_name : str
        Source model filename; used to derive output filename stems.
    default_output_filename : callable
        ``f(prefix=..., suffix=...) -> stem`` returning an output filename
        stem (no extension); typically
        ``ProgramTemplate.get_default_output_filename``.

    Returns
    -------
    list of dict
        One per-seed result dict (see :meth:`_process_seed`).
    """
    self.model_name = model_name
    self.default_output_filename = default_output_filename

    if self.params.preferred_cut_fallback and not self.params.use_preferred_cuts:
      print(
        'Note: preferred_cut_fallback=True has no effect because '
        'use_preferred_cuts=False (the geometric heuristic already applies '
        'to every bond).',
        file=self.logger,
      )

    model.process(
      pdb_interpretation_params=model.get_current_pdb_interpretation_params(),
      make_restraints=True,
    )

    model = self._apply_altloc_filter(model)

    self._include_nodes = self._resolve_residues_to_include(model)

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

    self._results = []
    for seed_index, (sel_str, seeds) in enumerate(seed_groups, start=1):
      result = self._process_seed(
        seed_index, seeds, model, adjacency, selection_str=sel_str
      )
      self._results.append(result)
    return self._results

  def get_results(self):
    """Return the per-seed result dicts produced by the last :meth:`run`."""
    return self._results

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
    visited_nodes, cap_nodes = self._region_grower.grow_region(
      qm_atoms, adjacency, model, max_depth=self.params.max_depth,
      preferred_cut_fallback=(
        self.params.use_preferred_cuts and self.params.preferred_cut_fallback),
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
      # In-memory hand-off for downstream consumers (e.g., mmtbx.geometry_restraints.qmi
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
    """Return the initial ``(iseq, sym_op)`` node set for the QM region.

    When *skip_radius_search* is False (default) all atoms within
    ``params.radius`` of every seed are included.  When *skip_radius_search*
    is True only the seed atoms themselves are added and BFS expansion
    (controlled by ``params.max_depth``) is relied upon to grow the region.

    The radius search is symmetry-aware: the KD-tree supplies ASU
    (identity-image) atoms within ``params.radius``, and
    :meth:`AtomGraphBuilder.seed_sym_nodes_within_radius` adds the
    symmetry-image atoms inside the same sphere.  Without the latter, a
    metal on a special position would seed only the identity copy of its
    coordinating residues, so the symmetry copies -- reached later by BFS
    -- would be truncated differently (e.g. cut at CA-CB where the ASU copy
    keeps CA/CB).

    Parameters
    ----------
    seeds : list of iotbx.pdb.hierarchy.atom
    model : mmtbx.model.manager

    Returns
    -------
    set of (int, sgtbx.rt_mx)
    """
    identity = _canon_op(sgtbx.rt_mx())
    qm_nodes = set()
    if self.params.skip_radius_search:
      for seed in seeds:
        qm_nodes.add((seed.i_seq, identity))
    else:
      for seed in seeds:
        mask = self._graph_builder.atoms_within_radius_best(
          seed, model, self.params.radius
        )
        for iseq in mask.iselection():
          qm_nodes.add((iseq, identity))
      qm_nodes |= self._graph_builder.seed_sym_nodes_within_radius(
        seeds, model, self.params.radius
      )
    qm_nodes |= self._include_nodes_for(seeds, model)
    return qm_nodes

  def _resolve_residues_to_include(self, model):
    """Identity-op nodes for every atom of every residue group touched by
    ``params.residues_to_include.selection``.

    The selection is expanded to whole residue groups, so a partial match
    (e.g. a single atom name) still pulls in the complete residue.  Returns
    an empty set when no selection is configured.

    Parameters
    ----------
    model : mmtbx.model.manager

    Returns
    -------
    set of (int, sgtbx.rt_mx)
    """
    selection = self.params.residues_to_include.selection
    if not selection:
      return set()
    identity = _canon_op(sgtbx.rt_mx())
    atoms = model.get_hierarchy().atoms()
    nodes = set()
    for iseq in model.selection(selection).iselection():
      residue_group = atoms[iseq].parent().parent()
      for atom_group in residue_group.atom_groups():
        for residue_atom in atom_group.atoms():
          nodes.add((residue_atom.i_seq, identity))
    return nodes

  def _include_nodes_for(self, seeds, model):
    """Return the ``residues_to_include`` nodes applicable to *seeds*.

    With ``scope=global`` every resolved include node applies to every seed
    region.  With ``scope=per_seed`` (default) an included residue is kept
    only when at least one of its atoms lies within ``proximity`` of a seed
    atom in this group; the whole residue is kept when any atom qualifies.

    Parameters
    ----------
    seeds : list of iotbx.pdb.hierarchy.atom
    model : mmtbx.model.manager

    Returns
    -------
    set of (int, sgtbx.rt_mx)
    """
    if not self._include_nodes:
      return set()
    scope = self.params.residues_to_include
    if scope.scope == 'global':
      return self._include_nodes

    # per_seed: union the proximity spheres of every seed in this group,
    # reusing the graph builder's cached KD-tree.
    near = flex.bool(model.get_number_of_atoms(), False)
    for seed in seeds:
      near = near | self._graph_builder.atoms_within_radius_best(
        seed, model, scope.proximity)
    near_iseqs = set(near.iselection())

    # Keep a residue whole if any of its atoms is inside a sphere.
    atoms = model.get_hierarchy().atoms()
    nodes_by_residue = defaultdict(list)
    for node in self._include_nodes:
      nodes_by_residue[id(atoms[node[0]].parent().parent())].append(node)
    kept = set()
    for residue_nodes in nodes_by_residue.values():
      if any(iseq in near_iseqs for (iseq, _op) in residue_nodes):
        kept.update(residue_nodes)
    return kept

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
    # is called per seed but the water set is the same every call.
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
      # waters.
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
    model_stem = os.path.splitext(os.path.basename(self.model_name))[0]

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

    return self.default_output_filename(
      prefix=model_stem,
      suffix=seed_suffix,
    )

