"""QMRegionBuilder: the per-seed QM-region extraction pipeline.

Orchestrates seed discovery, covalent-graph construction, BFS region growth
and hydrogen capping, and formats the sidecar PHIL that describes each region.
It is independent of the ProgramTemplate and the data manager, so it can be
driven directly in memory: construct it with a params object, then call
:meth:`run` with a model.
"""

import copy
import os
from collections import defaultdict

import libtbx.phil
from libtbx import Auto, group_args
from libtbx.utils import Sorry, null_out
from cctbx import sgtbx
from cctbx.array_family import flex
from cctbx.crystal import super_cell
from iotbx.pdb import common_residue_names_get_class
from scitbx import matrix

from mmtbx.geometry_restraints.endo_exo.util import neighbour_iseqs
from mmtbx.geometry_restraints.endo_exo.seeds import SeedFinder
from mmtbx.geometry_restraints.endo_exo.graph import AtomGraphBuilder
from mmtbx.geometry_restraints.endo_exo.cutting import BondCutDetector
from mmtbx.geometry_restraints.endo_exo.grow import QMRegionGrower
from mmtbx.geometry_restraints.endo_exo.capping import HydrogenCapper


# Sidecar PHIL emitted alongside each QM-region model file.  Indices are
# 0-based positional indices into the output model's atom list.
master_sidecar_phil_str = """
endo_exo_region {
  cap_iseqs = None
    .type = ints
    .help = "Indices of the boundary atoms in the QM-region output file; \
None when the region has none."
  cap_original_elements = None
    .type = strings
    .help = "Element each boundary atom carried before capping, parallel to \
cap_iseqs."
  completed_iseqs = None
    .type = ints
    .help = "Indices of atoms this run ADDED to repair a valence the input \
model left open at a chain break; they are not caps and carry no original \
element.  Their positions are constructed, not measured, so a consumer may \
wish to leave them free to relax."
  seed_iseqs = None
    .type = ints
    .help = "Indices of seed atoms in the QM-region output file (the metal \
atom in metal-scan mode, or all atoms matched by the selection string)."
  selection_string = None
    .type = str
    .help = "Original CCTBX selection string that seeded the region; \
None for metal-scan regions."
}
"""


def format_sidecar(result):
  """Return the sidecar PHIL text that describes one QM region.

  The sidecar records which atoms are hydrogen caps, the element each carried
  beforehand, which atoms were added to repair an open valence, which atoms
  seeded the region, and the selection string.  Every index it holds is
  positional into the region's model.

  Parameters
  ----------
  result : dict
      One region as returned by :meth:`QMRegionBuilder.run`.

  Returns
  -------
  str
  """
  # An empty list formats as a bare ``key =``, which does not parse back, so
  # empty lists are written as None (the declared default) instead.
  def _or_none(values):
    values = list(values)
    return values or None

  sidecar_phil = libtbx.phil.parse(master_sidecar_phil_str)
  sidecar_phil_extract = sidecar_phil.extract()
  region = sidecar_phil_extract.endo_exo_region
  region.cap_iseqs = _or_none(result['cap_iseqs'])
  region.cap_original_elements = _or_none(result['cap_original_elements'])
  region.completed_iseqs = _or_none(result['completed_iseqs'])
  region.seed_iseqs = _or_none(result['seed_iseqs'])
  region.selection_string = result['selection_string']
  return sidecar_phil.format(sidecar_phil_extract).as_str()


class QMRegionBuilder:
  """Build QM regions around seed sites for a single model.

  Parameters
  ----------
  params : group_args-like
      Resolved parameters (the ``master_phil_str`` scope extract from
      :class:`mmtbx.programs.endo_exo.Program`, or any object exposing the
      same attributes).
  logger : file-like or None, optional
      Destination for diagnostic messages; None discards them.
  """

  def __init__(self, params, logger=None):
    self.params = params
    if logger is None:
      logger = null_out()
    self.logger = logger
    self.model_name = None
    self.default_output_filename = None
    self._results = []
    self._include_nodes = set()
    self._hbond_cache = None
    self._waters_cache = None
    self._waters_cache_rep_frac = None
    self._graph_builder = AtomGraphBuilder()
    self._capper = HydrogenCapper(log=self.logger)
    self._bond_cut_detector = BondCutDetector(
      use_preferred_cuts=self.params.capping.preferred_cuts,
      log=self.logger,
    )
    self._region_grower = QMRegionGrower(
      self._bond_cut_detector,
      log=self.logger,
    )

  def run(self, model, model_name, default_output_filename):
    """Locate seeds and grow a QM region around each.

    Parameters
    ----------
    model : mmtbx.model.manager
        The input model.  Restraints are processed in place.
    model_name : str
        Source model filename; used to derive output filename stems.
    default_output_filename : callable
        ``f(prefix=..., suffix=...) -> stem`` returning an output filename
        stem (no extension), called with the model name as *prefix* and
        the region's own suffix.

    Returns
    -------
    list of dict
        One per-seed result dict (see :meth:`_process_seed`).

    Raises
    ------
    Sorry
        When no seed atoms are found, a selection matches no atoms, the model
        has no hydrogens where hydrogen bonds are needed, or the symmetry
        images need more chain IDs than are free.
    """
    self.model_name = model_name
    self.default_output_filename = default_output_filename

    # Only the bond graph is needed, so the CDL and symmetric-flip stages are
    # switched off: CDL ideals never reach the adjacency, and flipping would
    # move the input coordinates.  Metal auto-linking and the metal
    # coordination library are off as well, so the graph holds covalent bonds
    # only, with no metal-donor coordination bonds and no idealized Zn/Fe-S
    # proxies.
    interpretation_params = copy.deepcopy(
      model.get_current_pdb_interpretation_params())
    pdb_interpretation = interpretation_params.pdb_interpretation
    pdb_interpretation.restraints_library.cdl = False
    pdb_interpretation.restraints_library.mcl = False
    pdb_interpretation.flip_symmetric_amino_acids = False
    pdb_interpretation.automatic_linking.link_metals = False
    model.process(
      pdb_interpretation_params=interpretation_params,
      make_restraints=True,
    )

    model = self._apply_altloc_filter(model)

    self._include_nodes = self._resolve_residues_to_include(model)
    # Hydrogen bonds and waters are found once per run and reused across seed
    # groups.  They index this model's atoms, so they cannot outlive it.
    self._hbond_cache = None
    self._waters_cache = None

    seed_finder = SeedFinder()
    selection_strings = [s for s in (self.params.selection or []) if s]
    element_filter = [
      e for e in (self.params.element_filter or []) if e and e.strip()
    ]
    if element_filter and selection_strings:
      print(
        f"Note: ignoring element_filter={element_filter} because an "
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
      seeds_all = []
      for _, grp in seed_groups:
        seeds_all.extend(grp)
      self._print_seeds(seeds_all, label='seed atoms')

    # seeds_flat = [atom for _, grp in seed_groups for atom in grp]  # needed for contact_cutoff

    grm = model.get_restraints_manager().geometry
    sites_cart = model.get_sites_cart()
    bond_proxies_simple, bond_proxies_asu = grm.get_all_bond_proxies(
      sites_cart=sites_cart
    )
    asu_mappings = grm.pair_proxies(
      sites_cart=sites_cart).bond_proxies.asu_mappings()

    adjacency = self._graph_builder.build_adjacency(
      bond_proxies_simple, bond_proxies_asu, asu_mappings)

    # Seed-contact edges disabled (see contact_cutoff note in endo_exo.py).
    # cutoff = self.params.contact_cutoff
    # if not self.params.buffer.skip_search:
    #   added_edges = self._graph_builder.add_seed_contact_edges(
    #     seeds_flat, model, adjacency, cutoff=cutoff)
    #   print(f'Added {added_edges} seed-contact edges '
    #         f'(cutoff={cutoff:.2f} A)', file=self.logger)
    print(
      f'Always-included seed-centered radius: '
      f'{self.params.buffer.radius:.2f} A',
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

    * ``'all'``: no filtering; *model* is returned unchanged.
    * ``'auto'``: for each residue group containing non-blank altlocs,
      retain the altloc with the highest mean atom occupancy and drop
      the others.
    * any specific letter (e.g. ``'A'``): retain that letter where it
      is present; if a residue has non-blank altlocs but not the
      requested letter, fall back to the highest-occupancy altloc and
      emit a warning.

    Atoms with empty altloc (the shared backbone) are always kept.  The
    conformer kept is the one :meth:`iotbx.pdb.hierarchy.root.remove_alt_confs`
    keeps, here with its altloc label left in place.

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
        Keys: ``file_name``, ``n_atoms``, ``model``, ``seed_iseqs``,
        ``cap_iseqs``, ``cap_original_elements``, ``cap_anchor_iseqs``,
        ``completed_iseqs``, ``sym_image_provenance``, ``selection_string``.
    """
    qm_atoms = self._seed_qm_region(seeds, model)
    visited_nodes, cap_nodes = self._region_grower.grow_region(
      qm_atoms, adjacency, model, max_depth=self.params.max_search_depth)

    visited_nodes = self._add_hull_waters(model, visited_nodes)

    region = self._materialize_qm_region(
      model, visited_nodes, cap_nodes, seeds, adjacency)

    # Reported after materialisation and told which valences were repaired,
    # so that the report names only the atoms that are still a bond short.
    self._report_open_valences(
      model, visited_nodes, cap_nodes, adjacency, region.completed_nodes)

    file_name = self._make_output_filename(
      seed_index, seeds, selection_str=selection_str
    )

    return {
      'file_name': file_name,
      'n_atoms': region.model.get_number_of_atoms(),
      # In-memory hand-off: the truncated sub-model (no restraints manager)
      # and the positional indices of seed and cap atoms inside it, plus the
      # heavy-atom anchor each cap is bonded to (kept in memory only, not the
      # on-disk sidecar).
      'model': region.model,
      'seed_iseqs': region.seed_iseqs,
      'cap_iseqs': region.cap_iseqs,
      'cap_original_elements': region.cap_original_elements,
      'cap_anchor_iseqs': region.cap_anchor_iseqs,
      # Atoms added to repair an open valence.  They are not caps, and their
      # positions are constructed rather than measured, so a consumer that
      # freezes caps may still want to leave these free to relax.
      'completed_iseqs': region.completed_iseqs,
      # {sub-model i_seq: ((chain, resseq, resname, name, altloc),
      # symmetry_operation_xyz)} for symmetry-image atoms, so a metal-ligand
      # bond to one can be restrained against its ASU parent (kept in memory
      # only, not the on-disk sidecar).
      'sym_image_provenance': region.sym_image_provenance,
      'selection_string': selection_str,
    }

  # ------------------------------------------------------------------
  # Pipeline steps
  # ------------------------------------------------------------------

  def _seed_qm_region(self, seeds, model):
    """Return the initial ``(iseq, sym_op)`` node set for the QM region.

    When *buffer.skip_search* is False (default) all atoms within
    ``params.buffer.radius`` of every seed are included.  When
    *buffer.skip_search* is True only the seed atoms themselves are added and
    BFS expansion (controlled by ``params.max_search_depth``) is relied upon to
    grow the region.

    The radius search is symmetry-aware: the KD-tree supplies ASU
    (identity-image) atoms within ``params.buffer.radius``, and
    :meth:`AtomGraphBuilder.seed_sym_nodes_within_radius` adds the
    symmetry-image atoms inside the same sphere.  Both halves are needed for
    a seed on a special position, so that symmetry copies of its
    surroundings are seeded, and therefore truncated, exactly like the ASU
    copy.

    Parameters
    ----------
    seeds : list of iotbx.pdb.hierarchy.atom
    model : mmtbx.model.manager

    Returns
    -------
    set of (int, sgtbx.rt_mx)
    """
    identity = sgtbx.rt_mx()
    qm_nodes = set()
    if self.params.buffer.skip_search:
      for seed in seeds:
        qm_nodes.add((seed.i_seq, identity))
    else:
      for seed in seeds:
        mask = self._graph_builder.atoms_within_radius(
          seed, model, self.params.buffer.radius
        )
        for iseq in mask.iselection():
          qm_nodes.add((iseq, identity))
      qm_nodes |= self._graph_builder.seed_sym_nodes_within_radius(
        seeds, model, self.params.buffer.radius
      )
    qm_nodes |= self._include_nodes_for(seeds, model)
    qm_nodes |= self._hbond_partner_nodes(model, qm_nodes)
    return qm_nodes

  def _report_open_valences(self, model, visited_nodes, cap_nodes, adjacency,
                            completed_nodes=frozenset()):
    """Name region atoms whose backbone partner is missing from the model.

    Capping covers only bonds this code severs.  A residue whose neighbour
    is absent from the input, in a disordered loop or at the edge of an
    extract, has no bond to sever and none to cap, so its nitrogen or
    carbonyl carbon reaches the region a bond short and is a radical centre
    in whatever is run on it afterwards.  This method cannot place the
    missing atom, so it reports the defect rather than repairing it.

    Parameters
    ----------
    model : mmtbx.model.manager
    visited_nodes : set of (int, sgtbx.rt_mx)
    cap_nodes : dict
    adjacency : collections.defaultdict of set
    completed_nodes : frozenset of (int, sgtbx.rt_mx), optional
        Nodes whose valence was already repaired; excluded from the report.
    """
    atoms = model.get_hierarchy().atoms()
    # Keyed on the atom, not the node: symmetry images of one atom are one
    # problem, reported once.
    open_valences = set()
    for (iseq, _op) in visited_nodes:
      if (iseq, _op) in cap_nodes:
        continue
      atom = atoms[iseq]
      name = atom.name.strip().upper()
      if (name not in ('N', 'C')
          or atom.element.strip().upper() not in ('N', 'C')):
        continue
      # An atom whose valence materialisation completed is not short of
      # anything, so it is not reported.
      if (iseq, _op) in completed_nodes:
        continue
      residue_group = atom.parent().parent()
      partner = 'C' if name == 'N' else 'N'
      hydrogens = 0
      oxygens = 0
      on_backbone = False
      chain_continues = False
      for neighbour_iseq in neighbour_iseqs(adjacency, iseq):
        other = atoms[neighbour_iseq]
        if other.element_is_hydrogen():
          hydrogens += 1
          continue
        if other.parent().parent() != residue_group:
          if other.name.strip().upper() == partner:
            chain_continues = True    # nothing missing
            break
          continue                    # some other link, not the backbone one
        other_name = other.name.strip().upper()
        if other_name == 'CA' and other.element.strip().upper() == 'C':
          on_backbone = True
        if other.element.strip().upper() == 'O':
          oxygens += 1
      if chain_continues:
        continue
      # Only a peptide backbone has a neighbouring residue that can be
      # missing.  A ligand that happens to name an atom N or C is bonded
      # entirely within itself, and is short of nothing.
      if not on_backbone:
        continue
      # A real terminus is not short of anything: an amine carries its own
      # hydrogens, a carboxylate its second oxygen.
      if name == 'C' and oxygens >= 2:
        continue
      if name == 'N' and hydrogens >= 2:
        continue
      open_valences.add(
        f'{atom.parent().resname.strip()} '
        f'{residue_group.resseq.strip()} {atom.name.strip()}')

    if open_valences:
      print(
        f'Note: {len(open_valences)} atom(s) in this region are a bond '
        f'short because the neighbouring residue is absent from the model: '
        f'{", ".join(sorted(open_valences))}.',
        file=self.logger,
      )

  def _hbond_partner_nodes(self, model, qm_nodes):
    """Return nodes for atoms hydrogen-bonded to one already in *qm_nodes*.

    A region cut by distance runs through hydrogen bonds, leaving a donor
    without its acceptor.  The far atom joins as a seed rather than as a
    finished fragment, so the cut rules shape it like any other part of the
    region and a partner sidechain arrives truncated.

    Only partners of the nodes as they stand are taken, and those partners
    are not themselves followed, so a region cannot walk across a surface
    one hydrogen bond at a time.  Reading the seeds rather than the grown
    region keeps the result independent of BFS order and makes a wider
    sphere yield a superset of a narrower one.

    The result is monotone in the radius but not in the option itself:
    enabling it can make a region smaller.  Seeds are pre-loaded into the
    visited set, and the backbone CA-N rule withholds its cut until the
    next residue's amide is visited, so a partner seeded anywhere can
    satisfy that guard and enable a cut that would not otherwise happen.
    The cut prunes backwards, taking with it any residue reached only
    through it.

    No-op when ``params.include_hbond_partners`` is ``False``.

    Parameters
    ----------
    model : mmtbx.model.manager
    qm_nodes : set of (int, sgtbx.rt_mx)

    Returns
    -------
    set of (int, sgtbx.rt_mx)
        Partner nodes not already in *qm_nodes*.
    """
    if not self.params.include_hbond_partners:
      return set()

    ops_by_iseq = defaultdict(list)
    for iseq, op in qm_nodes:
      ops_by_iseq[iseq].append(op)

    partners = set()
    for record in self._hydrogen_bonds(model):
      donor = record.atom_D.index
      acceptor = record.atom_A.index
      # The operation relates the two halves of the bond, and which half it
      # moves depends on which atom of the pair is the hydrogen.
      if record.atom_H.index == record.i:
        near, far = donor, acceptor
      else:
        near, far = acceptor, donor
      for this, other, forward in ((near, far, True), (far, near, False)):
        if this not in ops_by_iseq:
          continue
        step = record.symop if forward else record.symop.inverse()
        for op in ops_by_iseq[this]:
          node = (other, op.multiply(step).new_denominators(1, 12))
          if node not in qm_nodes:
            partners.add(node)

    if partners:
      print(
        f'Adding {len(partners)} atom nodes hydrogen-bonded to the QM '
        f'region.',
        file=self.logger,
      )
    return partners

  def _hydrogen_bonds(self, model):
    """Return the model's hydrogen bonds, found once and cached.

    The search covers the whole structure and does not depend on which seed
    is being grown, so its result is shared by every seed group.  The
    records index this model's atoms and are invalid for any other model.

    Parameters
    ----------
    model : mmtbx.model.manager

    Returns
    -------
    list of libtbx.introspection.group_args
    """
    from mmtbx.nci import hbond

    if self._hbond_cache is None:
      if not model.has_hd():
        raise Sorry('Model has no hydrogens. Add them and run again.')
      self._hbond_cache = hbond.find(model=model).result
    return self._hbond_cache

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
    identity = sgtbx.rt_mx()
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
      near = near | self._graph_builder.atoms_within_radius(
        seed, model, scope.proximity)
    near_iseqs = set(near.iselection())

    # Keep a residue whole if any of its atoms is inside a sphere.
    atoms = model.get_hierarchy().atoms()
    nodes_by_residue = defaultdict(list)
    for node in self._include_nodes:
      nodes_by_residue[atoms[node[0]].parent().parent()].append(node)
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
    the visited nodes, not from ASU coordinates, so symmetry
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
    from scipy.spatial import ConvexHull, QhullError
    try:
      hull = ConvexHull(hull_cart)
    except QhullError:
      print('No hull waters added: the QM region encloses no volume.',
            file=self.logger)
      return visited_nodes

    centroid_frac = uc.fractionalize(hull_cart.mean())

    visited_iseqs_by_op = defaultdict(set)
    for (iseq, op) in visited_nodes:
      visited_iseqs_by_op[op].add(iseq)

    # Collect water residues once per run: ``_add_hull_waters`` is called per
    # seed but the water set is the same every call.
    if self._waters_cache is None:
      cached = []
      for chain in model.get_hierarchy().chains():
        for rg in chain.residue_groups():
          if common_residue_names_get_class(
              rg.atom_groups()[0].resname) != 'common_water':
            continue
          rep_iseq = rg.atom_groups()[0].atoms()[0].i_seq
          iseqs = []
          for ag in rg.atom_groups():
            for a in ag.atoms():
              iseqs.append(a.i_seq)
          cached.append((rep_iseq, iseqs))
      self._waters_cache = cached
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
        final_op = shift_op.multiply(sg_op).new_denominators(1, 12)
        already = visited_iseqs_by_op.get(final_op, set())
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
                              seeds, adjacency):
    """Build a materialized QM sub-model from the BFS nodes.

    Nodes carrying the same ``i_seq`` but different ``sym_op`` are
    realised as separate atoms in the output: the parent atom's
    position is transformed by the ``sym_op``'s fractional rotation
    plus translation, and the result lives in its own chain block.
    A symmetry image that positionally duplicates a kept atom group
    (special-position case: the fixing op maps the group onto itself)
    is dropped whole, keyed on its heavy atoms so an on-axis water is
    not left with its two disorder H orientations as extra protons.

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
    adjacency : collections.defaultdict of set
        Full covalent graph of the parent model, used to find the
        neighbours of an atom whose valence needs completing.

    Returns
    -------
    libtbx.group_args
        ``model`` (mmtbx.model.manager)
            Sub-model containing the materialized atoms.
        ``seed_iseqs`` (list of int)
            Positional indices of seed atoms in ``model``.
        ``cap_iseqs`` (list of int)
            Positional indices of cap atoms in ``model``.
        ``cap_original_elements`` (list of str)
            Element each cap atom carried before capping (parallel to
            ``cap_iseqs``); unchanged from the emitted element when
            ``capping.enable`` is False.
        ``cap_anchor_iseqs`` (list of int)
            Sorted, unique positional indices of the QM-region heavy atoms
            the caps are bonded to (their anchors).
        ``completed_nodes`` (frozenset of (int, sgtbx.rt_mx))
            Nodes whose open valence was repaired.
        ``completed_iseqs`` (list of int)
            Sorted positional indices of the atoms those repairs added.
        ``sym_image_provenance`` (dict)
            ``{sub-model i_seq: ((chain, resseq, resname, name, altloc),
            symmetry_operation_xyz)}`` for each materialized symmetry-image
            atom: the ASU parent's selection identity (not the atom object,
            whose parent chain dangles once this model is freed) plus the
            symmetry op.
    """
    import mmtbx.model as mmtbx_model

    parent_atoms = model.get_hierarchy().atoms()
    out_root, node_of_atom = self._copy_images(model, visited_nodes)
    atom_for_node = self._drop_coincident_images(
      out_root, node_of_atom, parent_atoms)
    out_root.atoms().reset_serial()
    orig_element_by_cap, anchor_atom_by_cap = self._cap_boundary(
      cap_nodes, atom_for_node)
    completed_nodes, completed_atoms = self._repair_open_valences(
      visited_nodes, cap_nodes, adjacency, parent_atoms, atom_for_node,
      set(orig_element_by_cap))
    if completed_nodes:
      # Sorted before anything is indexed.  The recorded indices are
      # positional into the emitted file, and a consumer re-processing the
      # region runs pdb_interpretation with sort_atoms on by default: an atom
      # appended at the end of its group would be moved to the front of the
      # hydrogen block, shifting every index recorded after it.  Sorting here
      # leaves that later sort with nothing to do.
      out_root.sort_atoms_in_place()
      # Renumbered afterwards, so serials follow the final order rather than
      # being blank on the appended atoms.
      out_root.atoms().reset_serial()
      print(
        f'Repaired {len(completed_nodes)} open valence(s) left by a gap in '
        f'the input model, adding {len(completed_atoms)} atom(s).',
        file=self.logger,
      )

    # Assemble the final mmtbx.model.manager
    model_sel = mmtbx_model.manager(
      model_input=None,
      pdb_hierarchy=out_root,
      crystal_symmetry=model.crystal_symmetry())

    # ``atom_for_node`` already holds the materialized atoms.  Every index
    # below is read only after the wrapping above, which is what gives the
    # copied atoms their final ``.i_seq``.
    identity = sgtbx.rt_mx()
    seed_indices = []
    for s in seeds:
      a = atom_for_node.get((s.i_seq, identity))
      if a is not None:
        seed_indices.append(a.i_seq)
    seed_indices.sort()
    cap_atoms_sorted = []
    for (iseq, op) in cap_nodes:
      a = atom_for_node.get((iseq, op))
      if a is not None:
        cap_atoms_sorted.append(a)
    cap_atoms_sorted.sort(key=lambda a: a.i_seq)
    # Indexed rather than fetched with ``.get``: both lists are built from
    # ``cap_nodes`` under the same not-None filter, so a missing entry is a
    # bug and should raise.
    cap_original_elements = [orig_element_by_cap[a]
                             for a in cap_atoms_sorted]
    cap_anchor_indices = sorted(set(
      anchor_atom_by_cap[a].i_seq for a in cap_atoms_sorted
      if a in anchor_atom_by_cap))

    return group_args(
      model=model_sel,
      seed_iseqs=seed_indices,
      cap_iseqs=[a.i_seq for a in cap_atoms_sorted],
      cap_original_elements=cap_original_elements,
      cap_anchor_iseqs=cap_anchor_indices,
      completed_nodes=frozenset(completed_nodes),
      completed_iseqs=sorted(a.i_seq for a in completed_atoms),
      sym_image_provenance=self._sym_image_provenance(
        atom_for_node, parent_atoms))

  def _copy_images(self, model, visited_nodes):
    """Copy the atoms of *visited_nodes* into a new hierarchy, one block of
    chains per symmetry operation, the identity first.

    Each image's sites and anisotropic U are transformed by its operation and
    its chains get the IDs :meth:`_image_chain_ids` assigns.

    Parameters
    ----------
    model : mmtbx.model.manager
    visited_nodes : set of (int, sgtbx.rt_mx)

    Returns
    -------
    out_root : iotbx.pdb.hierarchy.root
    node_of_atom : dict
        ``{output atom: (parent_iseq, sym_op)}``.
    """
    import iotbx.pdb

    unit_cell = model.crystal_symmetry().unit_cell()
    n_parent_atoms = model.get_number_of_atoms()

    # Group iseqs by symmetry operation, the identity first.
    by_op = defaultdict(set)
    for (iseq, op) in visited_nodes:
      by_op[op].add(iseq)
    ops = sorted(by_op, key=lambda op: (not op.is_unit_mx(), op.as_xyz()))
    chain_id_for = self._image_chain_ids(model.get_hierarchy(), by_op, ops)

    # Build an empty target hierarchy and append per-sym_op sub-models
    out_root = iotbx.pdb.hierarchy.root()
    out_hier_model = iotbx.pdb.hierarchy.model(id="1")
    out_root.append_model(out_hier_model)

    # As we copy each per-sym_op sub-hierarchy into the master,
    # remember which parent node each fresh atom came from.
    # ``model.select(mask)`` preserves atom order, so the k-th atom in
    # the copied sub-hierarchy corresponds to the k-th parent i_seq in
    # ascending order from ``by_op[op]``.
    node_of_atom = {}  # output atom -> (parent_iseq, op)
    for op in ops:
      iseqs_in_order = sorted(by_op[op])
      mask = flex.bool(n_parent_atoms, False)
      for i in iseqs_in_order:
        mask[i] = True
      sub_model = model.select(mask)
      sub_hier = sub_model.get_hierarchy()

      # Transform coordinates and anisotropic U from ASU positions to the
      # image.
      if not op.is_unit_mx():
        sub_model.set_sites_cart(super_cell.sym_equiv_sites_cart(
          sites_cart=sub_model.get_sites_cart(), unit_cell=unit_cell,
          rt_mx=op))
        sub_hier.atoms().rotate_uij_and_clear_sigmas(
          r_cart=unit_cell.matrix_cart(op.r()))
        # Rename chains so this image is distinguishable
        for ch in sub_hier.chains():
          ch.id = chain_id_for[(ch.id.strip(), op)]

      # Detached-copy each chain into the master and record the parent
      # node for every fresh atom object.
      k = 0
      for ch in sub_hier.chains():
        new_ch = ch.detached_copy()
        out_hier_model.append_chain(new_ch)
        for new_atom in new_ch.atoms():
          node_of_atom[new_atom] = (iseqs_in_order[k], op)
          k += 1
    return out_root, node_of_atom

  @staticmethod
  def _image_chain_ids(parent_hier, by_op, ops):
    """Return the chain ID each chain of the region takes under each
    symmetry operation.

    The identity keeps the original IDs; every image takes a free
    single-character ID per chain it contains.  PDB's chain field is one
    character wide; multi-char IDs overflow into the residue-name field and
    break the bond graph that ``make_restraints`` builds.

    Parameters
    ----------
    parent_hier : iotbx.pdb.hierarchy.root
    by_op : dict
        ``{sym_op: set of parent i_seq}``.
    ops : list of sgtbx.rt_mx
        The keys of *by_op*, the identity first.

    Returns
    -------
    dict
        ``{(original_chain_id, sym_op): chain_id}``.

    Raises
    ------
    Sorry
        When the images need more chain IDs than the model leaves free.
    """
    parent_atoms = parent_hier.atoms()
    orig_chain_ids = sorted({ch.id.strip() for ch in parent_hier.chains()})
    available_chain_ids = iter(c for c in
      "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
      if c not in orig_chain_ids)
    chain_id_for = {}
    for op in ops:
      region_chain_ids = sorted({
        parent_atoms[i].parent().parent().parent().id.strip()
        for i in by_op[op]})
      for orig_id in region_chain_ids:
        if op.is_unit_mx():
          chain_id_for[(orig_id, op)] = orig_id
        else:
          try:
            chain_id_for[(orig_id, op)] = next(available_chain_ids)
          except StopIteration:
            raise Sorry(
              'The symmetry images in the QM region need more '
              'single-character chain IDs than the model leaves free.')
    return chain_id_for

  @staticmethod
  def _drop_coincident_images(out_root, node_of_atom, parent_atoms):
    """Drop the symmetry images that duplicate a kept atom group, and index
    the surviving atoms by node.

    This is the special-position case, where the fixing op maps the group
    onto itself.  Emptied residue groups and chains are dropped too.  The
    test is keyed per group on its heavy atoms, because an on-axis water
    shares its oxygen with the kept image but places its two disorder
    hydrogens apart; a group whose heavy atoms coincide with kept ones is
    dropped whole.  The tolerance absorbs coordinate drift.

    Parameters
    ----------
    out_root : iotbx.pdb.hierarchy.root
        Changed in place.
    node_of_atom : dict
        ``{output atom: (parent_iseq, sym_op)}``.
    parent_atoms : iotbx.pdb.hierarchy.af_shared_atom

    Returns
    -------
    dict
        ``{(parent_iseq, sym_op): surviving atom}``.
    """
    POS_TOL = 0.2
    # An image is compared only against other images of the same parent atom
    # group, so two altlocs of one residue, which share backbone positions,
    # never cancel each other.  A single coincidence drops the whole group,
    # since a group straddling the axis has some atoms on it and some off,
    # and keeping it would put two nuclei in one place.
    out_hier_model = out_root.models()[0]
    seen_images = defaultdict(list)  # parent group -> [(op, heavy atom), ]
    atom_for_node = {}  # (parent_iseq, op) -> surviving atom
    for ch in list(out_hier_model.chains()):
      for rg in list(ch.residue_groups()):
        for ag in list(rg.atom_groups()):
          heavy = [a for a in ag.atoms() if not a.element_is_hydrogen()]
          is_duplicate = False
          if heavy:
            parent_iseq, image_op = node_of_atom[heavy[0]]
            group_key = parent_atoms[parent_iseq].parent().memory_id()
            other_images = [kept for (kept_op, kept) in seen_images[group_key]
                            if kept_op != image_op]
            for a in heavy:
              for prev in other_images:
                if a.distance(prev) < POS_TOL:
                  is_duplicate = True
                  break
              if is_duplicate:
                break
          if is_duplicate:
            rg.remove_atom_group(ag)
            continue
          for atom in ag.atoms():
            if not atom.element_is_hydrogen():
              iseq_of, op_of = node_of_atom[atom]
              seen_images[
                parent_atoms[iseq_of].parent().memory_id()].append(
                  (op_of, atom))
            atom_for_node[node_of_atom[atom]] = atom
        if len(list(rg.atom_groups())) == 0:
          ch.remove_residue_group(rg)
      if len(list(ch.residue_groups())) == 0:
        out_hier_model.remove_chain(ch)
    return atom_for_node

  def _cap_boundary(self, cap_nodes, atom_for_node):
    """Record each boundary atom's element and anchor, then cap it.

    The record describes the region boundary, which exists whether or not a
    hydrogen is placed there, so only the placement is gated on
    ``capping.enable``.

    Parameters
    ----------
    cap_nodes : dict
        ``{(cap_iseq, cap_op): (anchor_iseq, anchor_op)}``.
    atom_for_node : dict
        ``{(parent_iseq, sym_op): materialized atom}``.

    Returns
    -------
    orig_element_by_cap : dict
        ``{cap atom: element before capping}``.
    anchor_atom_by_cap : dict
        ``{cap atom: anchor atom}``.
    """
    orig_element_by_cap = {}
    anchor_atom_by_cap = {}
    for cap_node, anchor_node in cap_nodes.items():
      cap_atom = atom_for_node.get(cap_node)
      anchor_atom = atom_for_node.get(anchor_node)
      if cap_atom is None:
        continue
      orig_element_by_cap[cap_atom] = cap_atom.element.strip()
      if anchor_atom is not None:
        # Hold the atom, not its i_seq: the copies do not carry a final
        # i_seq until the hierarchy is wrapped in a model, so reading it here
        # yields the pre-copy value for every symmetry image.
        anchor_atom_by_cap[cap_atom] = anchor_atom
      if self.params.capping.enable:
        self._capper.cap_atom(anchor_atom, cap_atom)
    return orig_element_by_cap, anchor_atom_by_cap

  def _repair_open_valences(self, visited_nodes, cap_nodes, adjacency,
                            parent_atoms, atom_for_node, cap_atoms):
    """Complete the backbone valences a gap in the input model leaves open.

    Gated on ``capping.enable``, because capping and these repairs settle
    the same question, what to do with a severed valence.

    Parameters
    ----------
    visited_nodes : set of (int, sgtbx.rt_mx)
    cap_nodes : dict
    adjacency : collections.defaultdict of set
    parent_atoms : iotbx.pdb.hierarchy.af_shared_atom
    atom_for_node : dict
        ``{(parent_iseq, sym_op): materialized atom}``.
    cap_atoms : set of iotbx.pdb.hierarchy.atom
        The materialized cap atoms, which are never repaired.

    Returns
    -------
    completed_nodes : list of (int, sgtbx.rt_mx)
        The nodes repaired.
    completed_atoms : list of iotbx.pdb.hierarchy.atom
        The atoms the repairs added.
    """
    completed_nodes = []
    completed_atoms = []
    if not self.params.capping.enable:
      return completed_nodes, completed_atoms

    # Complete any two-coordinate backbone nitrogen the input left short.
    # The cut rules never reach these: a nitrogen inside the buffer sphere is
    # a seed, so its CA-N bond is never cut-tested and it would otherwise
    # arrive as an aminyl centre.
    for node, nitrogen, neighbours in self._open_valence_sites(
        QMRegionGrower._is_chain_break_amine, visited_nodes, cap_nodes,
        adjacency, parent_atoms, atom_for_node, cap_atoms):
      added = self._capper.complete_amine(nitrogen, neighbours)
      if added:
        completed_nodes.append(node)
        completed_atoms.extend(added)

    # The other half of the same severed peptide bond.  This repair guesses
    # nothing: the carbon holds its CA and O, so the missing nitrogen is
    # fixed by sp2 planarity.  Cutting the carbon away instead is not an
    # option, because that would delete a carbonyl oxygen that may
    # coordinate the metal.
    for node, carbon, neighbours in self._open_valence_sites(
        QMRegionGrower._is_chain_break_carbonyl, visited_nodes, cap_nodes,
        adjacency, parent_atoms, atom_for_node, cap_atoms):
      # Which repair depends on whether the chain continues past this
      # residue.  Later polymer means an interior gap, so the missing
      # partner is the next residue's nitrogen; nothing later means the
      # chain ends here and what is missing is the second carboxylate
      # oxygen, which depositors often leave unmodelled.
      if QMRegionGrower._has_later_polymer_residue(parent_atoms[node[0]]):
        added = self._capper.complete_carbonyl(carbon, neighbours)
      else:
        added = self._capper.complete_carboxylate(carbon, neighbours)
      if added:
        completed_nodes.append(node)
        completed_atoms.extend(added)
    return completed_nodes, completed_atoms

  @staticmethod
  def _open_valence_sites(is_open, visited_nodes, cap_nodes, adjacency,
                          parent_atoms, atom_for_node, cap_atoms):
    """Return the materialized atoms *is_open* finds a bond short, with
    their materialized neighbours.

    Parameters
    ----------
    is_open : callable
        ``is_open(iseq, adjacency, parent_atoms) -> bool``.
    visited_nodes : set of (int, sgtbx.rt_mx)
    cap_nodes : dict
        Nodes that are caps, never repaired.
    adjacency : collections.defaultdict of set
    parent_atoms : iotbx.pdb.hierarchy.af_shared_atom
    atom_for_node : dict
        ``{(parent_iseq, sym_op): materialized atom}``.
    cap_atoms : set of iotbx.pdb.hierarchy.atom

    Returns
    -------
    list of tuple
        ``(node, atom, neighbours)`` per site, in node order.
    """
    sites = []
    for (iseq, op) in sorted(visited_nodes):
      if (iseq, op) in cap_nodes:
        continue
      if not is_open(iseq, adjacency, parent_atoms):
        continue
      atom = atom_for_node.get((iseq, op))
      if atom is None or atom in cap_atoms:
        continue
      neighbours = []
      for (j, edge_op) in adjacency[iseq]:
        neighbour = atom_for_node.get(
          (j, op.multiply(edge_op).new_denominators(1, 12)))
        if neighbour is not None:
          neighbours.append(neighbour)
      sites.append(((iseq, op), atom, neighbours))
    return sites

  @staticmethod
  def _sym_image_provenance(atom_for_node, parent_atoms):
    """Return the ASU parent of each materialized symmetry-image atom.

    The parent's selection identity is recorded rather than the atom object,
    whose parent chain dangles once this filtered model is freed after the
    region is returned.  Angles cannot cross symmetry, so only bond
    restraints can use this.  Call once the atoms carry their final i_seq.

    Parameters
    ----------
    atom_for_node : dict
        ``{(parent_iseq, sym_op): materialized atom}``.
    parent_atoms : iotbx.pdb.hierarchy.af_shared_atom

    Returns
    -------
    dict
        ``{sub-model i_seq: ((chain, resseq, resname, name, altloc),
        symmetry_operation_xyz)}``.
    """
    sym_image_provenance = {}
    for (parent_iseq, op), atom in atom_for_node.items():
      if op.is_unit_mx():
        continue
      pa = parent_atoms[parent_iseq]
      ag = pa.parent()
      rg = ag.parent()
      ch = rg.parent()
      ident = (ch.id.strip(), rg.resseq.strip(), ag.resname.strip(),
               pa.name.strip(), ag.altloc.strip())
      sym_image_provenance[atom.i_seq] = (ident, op.as_xyz())
    return sym_image_provenance

  def _make_output_filename(self, seed_index, seeds, selection_str=None):
    """Build the output filename stem from the model name and seed identity.

    When *selection_str* is given the stem is::

        {model_stem}_sel{seed_index:03d}_within{radius}A_depth{depth}

    When seeds are individual metal atoms the stem is::

        {model_stem}_{name}_chain{chain}_res{resseq}{icode}_within{radius}A_depth{depth}

    with ``_alt{altloc}`` after the residue when the seed has an altloc.

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
      f'_within{self.params.buffer.radius:.2f}A'
      f'_depth{self.params.max_search_depth}'
    )

    if selection_str:
      seed_suffix = f'_sel{seed_index:03d}' + common_suffix
    else:
      seed = seeds[0]
      atom_group = seed.parent()
      residue_group = atom_group.parent()
      chain = residue_group.parent()
      seed_suffix = (
        f'_{seed.name.strip()}_chain{chain.id.strip()}'
        f'_res{residue_group.resseq.strip()}{residue_group.icode.strip()}')
      if atom_group.altloc.strip():
        seed_suffix += f'_alt{atom_group.altloc.strip()}'
      seed_suffix += common_suffix

    return self.default_output_filename(
      prefix=model_stem,
      suffix=seed_suffix,
    )
