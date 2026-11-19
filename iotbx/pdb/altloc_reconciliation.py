"""
Make altloc labels of interacting residues consistent, so that the same
label means the same state. Implementation of
iotbx.pdb.hierarchy.root.reconcile_altlocs.

Each residue_group with alternative conformations is a site. Conformer
labels are those with at least one non-hydrogen atom; hydrogen-only labels
(e.g. neutron H/D exchange sites) are not conformers and are never
relabeled or scored as such. Two sites
interact if any of their altloc atoms are within distance_cutoff (crystal
symmetry included). For every pair of interacting sites, atoms of labels
paired after relabeling are scored: van der Waals overlaps beyond
clash_tolerance are penalized and covalent links (e.g. peptide bond between
residues with backbone altlocs) are rewarded; pairs separated by two to
four bonds are ignored. Altloc atoms of different residues are considered
bonded through polymer links (C-N, O3'-P) or disulfides: overlapping
alternate conformers often put other atoms at bonding distance, and current
labels cannot be trusted to tell which pairs are real bonds. Covalent
links to ligands (e.g. glycosylation) are also bonds, if their length is
within bond_tolerance of the sum of covalent radii and the same two atoms
are linked in every pairing of labels (overlapping conformers make such
contacts in one pairing only). Links count once per pair of labels of two
residues. Metal coordination is too ambiguous to pair labels: contacts of
metals are clashes only if much shorter than a bond. Connected
groups of sites are relabeled to the lowest-scoring assignment, changing as
few sites as possible; a site is relabeled only if that lowers the score by
more than significant_clash (small differences between nearly identical
conformers are noise). If the best assignment leaves a pair of sites with
a clash (overlaps beyond the tolerance, between two paired labels) worse
than the pair alone could reach, by more than unresolved_clash, the pair
cannot be made consistent: it is
reported and left unchanged together with the sites it interacts with, and
the rest of the group is relabeled consistently with them. Clashes that no
relabeling can remove are bad geometry, not inconsistent labels, and are
not reported.
Only labels change: coordinates, occupancies and B-factors do not.
"""
from __future__ import absolute_import, division, print_function
import itertools
from cctbx import crystal, uctbx
from cctbx.array_family import flex
from cctbx.eltbx import covalent_radii, van_der_waals_radii
from iotbx.pdb import common_residue_names_get_class
from libtbx import group_args
from libtbx.utils import null_out

bond_tolerance = 0.4
bond_search_cutoff = 3.0
polymer_links = [("C", "N"), ("O3'", "P")]
significant_clash = 0.1 # smaller score changes are noise
max_block_combinations = 256 # joint relabelings of a site and neighbors
unresolved_clash = 0.3 # smaller clashes that cannot be removed are not reported

def _element(atom):
  e = atom.element.strip()
  if (len(e) == 0): e = atom.name.strip()[:1]
  e = e.capitalize()
  if (e == "D"): e = "H"
  return e

def _vdw_radius(e):
  return van_der_waals_radii.vdw.table.get(e, 1.7)

non_metals = set(["H", "B", "C", "N", "O", "F", "SI", "P", "S", "CL", "AS",
  "SE", "BR", "TE", "I", "HE", "NE", "AR", "KR", "XE"])

def _kind(resname, element):
  """polymer, water, metal, ion (non-metal, e.g. halide) or other
  (ligands, saccharides...)"""
  c = common_residue_names_get_class(resname)
  if (c in ["common_amino_acid", "modified_amino_acid", "common_rna_dna",
            "modified_rna_dna"]): return "polymer"
  if (c == "common_water"): return "water"
  if (c == "common_element"):
    if (element.upper() in non_metals): return "ion"
    return "metal"
  return "other"

def _covalent_radius(e):
  try: return covalent_radii.table(e).radius()
  except RuntimeError: return 0.77

class _atoms(object):
  """Atoms used for scoring, with their site, label and radii"""
  def __init__(self, model, include_hydrogens):
    self.sites = [] # residue_groups with altlocs
    self.labels = [] # sorted conformer altlocs, per site
    self.residue = [] # per atom: index of residue_group in model
    self.chain = []
    self.kind = []
    self.name = []
    self.site = [] # per atom: index into self.sites or -1
    self.altloc = []
    self.vdw = []
    self.cov = []
    sites_cart = flex.vec3_double()
    residue_groups = [(i_chain, rg)
      for i_chain, chain in enumerate(model.chains())
      for rg in chain.residue_groups()]
    for i_rg, (i_chain, rg) in enumerate(residue_groups):
      labels = sorted(set([ag.altloc for ag in rg.atom_groups()
        if ag.altloc.strip() != "" and
           not all([a.element_is_hydrogen() for a in ag.atoms()])]))
      i_site = -1
      if (len(labels) > 0):
        i_site = len(self.sites)
        self.sites.append(rg)
        self.labels.append(labels)
      for ag in rg.atom_groups():
        for atom in ag.atoms():
          if (not include_hydrogens and atom.element_is_hydrogen()): continue
          e = _element(atom)
          self.residue.append(i_rg)
          self.chain.append(i_chain)
          self.kind.append(_kind(ag.resname, e))
          self.name.append(atom.name.strip())
          self.site.append(i_site if ag.altloc in labels else -1)
          self.altloc.append(ag.altloc)
          self.vdw.append(_vdw_radius(e))
          self.cov.append(_covalent_radius(e))
          sites_cart.append(atom.xyz)
    self.sites_cart = sites_cart

  def is_altloc(self, i):
    return self.site[i] >= 0 and self.altloc[i].strip() != ""

def _symmetry(crystal_symmetry, sites_cart, distance_cutoff):
  if (crystal_symmetry is not None and not crystal_symmetry.is_nonsense()):
    return crystal_symmetry
  # Isolated molecule (or placeholder/incomplete symmetry, e.g. CRYST1 1 1 1
  # in cryo-EM models): box large enough to exclude contacts between images
  return crystal.symmetry(
    unit_cell=uctbx.non_crystallographic_unit_cell(
      sites_cart=sites_cart, buffer_layer=distance_cutoff+1),
    space_group_symbol="P 1")

def _pairs(symmetry, sites_cart, selection, distance_cutoff):
  """Yield (i, j, distance) for selected atoms within distance_cutoff,
  once per symmetry operation, i and j indexing sites_cart"""
  i_seqs = selection.iselection()
  if (i_seqs.size() == 0): return
  sites = sites_cart.select(i_seqs)
  uc = symmetry.unit_cell()
  sites_frac = uc.fractionalize(sites)
  pst = symmetry.special_position_settings().pair_asu_table(
    distance_cutoff=distance_cutoff,
    sites_frac=sites_frac).extract_pair_sym_table()
  ii, jj = flex.size_t(), flex.size_t()
  for i, pair_sym_dict in enumerate(pst):
    for j, ops in pair_sym_dict.items():
      for k in range(len(ops)): # faster than iterating ops
        op = ops[k]
        if (op.is_unit_mx()):
          ii.append(i)
          jj.append(j)
        else:
          d = uc.distance(sites_frac[i], op*sites_frac[j])
          yield i_seqs[i], i_seqs[j], d
  distances = (sites.select(ii) - sites.select(jj)).norms()
  for i, j, d in zip(i_seqs.select(ii), i_seqs.select(jj), distances):
    yield i, j, d

def _linked(atoms, i, j):
  """i and j may be covalently bonded"""
  if (atoms.residue[i] == atoms.residue[j]): return True
  if (atoms.altloc[i].strip() == "" or atoms.altloc[j].strip() == ""):
    return True
  if (atoms.name[i] == "SG" and atoms.name[j] == "SG"): return True
  if (atoms.chain[i] != atoms.chain[j]): return False
  for k, l in [(i, j), (j, i)]:
    if (atoms.residue[l] == atoms.residue[k] + 1 and
        (atoms.name[k], atoms.name[l]) in polymer_links):
      return True
  return False

def _ligand_link(atoms, i, j):
  """i and j may be a covalent bond with a ligand or between ligands (e.g.
  glycosylation)"""
  kinds = set([atoms.kind[i], atoms.kind[j]])
  if (len(kinds & set(["metal", "ion", "water"])) > 0): return False
  return kinds != set(["polymer"])

def _coordination(atoms, i, j, d):
  """Contact of a metal not much shorter than a bond: may be coordination,
  neither a bond nor a clash. Coordination is not used to pair labels:
  overlapping alternate conformers often put atoms at coordination
  distance."""
  return ((atoms.kind[i] == "metal" or atoms.kind[j] == "metal") and
          d > atoms.cov[i] + atoms.cov[j] - bond_tolerance)

def _bonded(atoms, i, j, d, links):
  """links: accepted ligand links (i, j), i < j, from _ligand_links"""
  if (d >= atoms.cov[i] + atoms.cov[j] + bond_tolerance): return False
  return _linked(atoms, i, j) or (min(i, j), max(i, j)) in links

def _ligand_links(atoms, pairs):
  """Links to ligands and ions among (i, j, d) pairs, as a set of (i, j),
  i < j. Distance must be within bond_tolerance of a bond. A real link
  joins the same two atoms in every pairing of labels of the residues;
  alternate conformers that happen to overlap at bonding distance do it in
  one pairing only."""
  candidates = {}
  for i, j, d in pairs:
    bond = atoms.cov[i] + atoms.cov[j]
    if (i == j or abs(d - bond) >= bond_tolerance or _linked(atoms, i, j) or
        not _ligand_link(atoms, i, j)):
      continue
    if (atoms.residue[i] > atoms.residue[j]): i, j = j, i
    key = (atoms.residue[i], atoms.residue[j], atoms.name[i], atoms.name[j])
    candidates.setdefault(key, set()).add((i, j))
  def n_labels(i):
    if (atoms.site[i] < 0): return 1
    return len(atoms.labels[atoms.site[i]])
  result = set()
  for key, ij in candidates.items():
    i, j = list(ij)[0]
    label_pairs = set([(atoms.altloc[i], atoms.altloc[j]) for i, j in ij])
    if (len(label_pairs) >= min(n_labels(i), n_labels(j))):
      result.update([(min(i, j), max(i, j)) for i, j in ij])
  return result

def _compatible(atoms, x, a, b):
  """Atom x may be in the same conformer as a and b (a, b paired)"""
  for y in (a, b):
    if (atoms.residue[x] == atoms.residue[y]):
      return atoms.altloc[x].strip() == "" or \
        atoms.altloc[x] == atoms.altloc[y]
  return True

def _excluded(atoms, neighbors, a, b):
  """a and b (paired) are separated by two, three or four bonds"""
  for x in neighbors[a]:
    if (not _compatible(atoms, x, a, b)): continue
    for y in neighbors[b]:
      if (not _compatible(atoms, y, a, b)): continue
      if (x == y or y in neighbors[x]): return True
      for z in neighbors[x] & neighbors[y]:
        if (_compatible(atoms, z, a, b)): return True
  return False

def _edge_terms(atoms, symmetry, distance_cutoff, clash_tolerance):
  """For each pair of interacting sites (s, t), s < t:
  {(label_s, label_t): [clash, n_bonds]} over atom pairs with these labels"""
  n = len(atoms.site)
  # covalent neighbors (any atom) for 1-3 to 1-5 exclusions
  neighbors = [set() for _ in range(n)]
  pairs = list(_pairs(symmetry, atoms.sites_cart, flex.bool(n, True),
    bond_search_cutoff))
  links = _ligand_links(atoms, pairs)
  for i, j, d in pairs:
    if (i != j and _bonded(atoms, i, j, d, links)):
      neighbors[i].add(j)
      neighbors[j].add(i)
  selection = flex.bool([atoms.is_altloc(i) for i in range(n)])
  # farther pairs can neither clash nor be bonded
  max_contact = max(2*max(atoms.vdw) - clash_tolerance,
                    2*max(atoms.cov) + bond_tolerance)
  edges = {}
  for a, b, d in _pairs(symmetry, atoms.sites_cart, selection,
      min(distance_cutoff, max_contact)):
    s, t = atoms.site[a], atoms.site[b]
    if (s == t): continue
    if (_bonded(atoms, a, b, d, links)):
      clash, n_bonds = 0, 1
    else:
      clash = atoms.vdw[a] + atoms.vdw[b] - clash_tolerance - d
      if (clash <= 0 or _coordination(atoms, a, b, d)): continue
      if (_excluded(atoms, neighbors, a, b)): continue
      n_bonds = 0
    if (s > t):
      s, t, a, b = t, s, b, a
    terms = edges.setdefault((s, t), {}).setdefault(
      (atoms.altloc[a], atoms.altloc[b]), [0, 0])
    terms[0] += clash
    terms[1] += n_bonds
  return edges

def _clusters(pairs):
  """Connected groups of sites of pairs"""
  sites = sorted(set([k for pair in pairs for k in pair]))
  position = dict([(k, i) for i, k in enumerate(sites)])
  return [[sites[i] for i in group] for group in _components(len(sites),
    [(position[k], position[l]) for k, l in pairs])]

def _components(n_sites, edges):
  parent = list(range(n_sites))
  def root(i):
    while (parent[i] != i):
      parent[i] = parent[parent[i]]
      i = parent[i]
    return i
  for s, t in edges:
    parent[root(s)] = root(t)
  groups = {}
  for s in sorted(set([s for e in edges for s in e])):
    groups.setdefault(root(s), []).append(s)
  return sorted(groups.values())

class _component(object):
  """Search over label permutations of a connected group of sites"""
  def __init__(self, sites, labels, edges, bond_weight):
    self.sites = sites
    self.options = [list(itertools.permutations(labels[s])) for s in sites]
    self.labels = [labels[s] for s in sites]
    position = dict([(s, k) for k, s in enumerate(sites)])
    self.edges = []
    for (s, t), terms in edges.items():
      if (s in position):
        # a link between two residues is one bond: more bonds between two
        # labels come from overlapping alternate conformers
        self.edges.append((position[s], position[t], [
          (ls, lt, clash - bond_weight*min(n_bonds, 1), clash)
          for (ls, lt), (clash, n_bonds) in terms.items()]))
    # per site: (other site, terms, index of this site's label in terms)
    self.incident = [[] for k in sites]
    for k, l, terms in self.edges:
      self.incident[k].append((l, terms, 0))
      self.incident[l].append((k, terms, 1))

  def mapping(self, k, option):
    return dict(zip(self.labels[k], self.options[k][option]))

  def score(self, assignment):
    """(cost, residual clash) for assignment: option index per site"""
    maps = [self.mapping(k, o) for k, o in enumerate(assignment)]
    cost, clash = 0, 0
    for k, l, terms in self.edges:
      for ls, lt, c, cl in terms:
        if (maps[k][ls] == maps[l][lt]):
          cost += c
          clash += cl
    return cost, clash

  def local_cost(self, k, option, assignment):
    """Cost of the edges of site k with option, other sites as assigned"""
    m = self.mapping(k, option)
    result = 0
    for l, terms, side in self.incident[k]:
      ml = self.mapping(l, assignment[l])
      for term in terms:
        if (m[term[side]] == ml[term[1-side]]):
          result += term[2]
    return result

  def frustrated(self, assignment, threshold):
    """Pairs of sites whose worst clash between paired labels exceeds, by
    more than threshold, the lowest the pair alone can reach (labels mapped
    to any label of the pair). Clashes that no relabeling can remove are bad
    geometry, not inconsistent labels."""
    maps = [self.mapping(k, o) for k, o in enumerate(assignment)]
    result = []
    for k, l, terms in self.edges:
      clash = max([0] + [cl for ls, lt, c, cl in terms
        if maps[k][ls] == maps[l][lt]])
      labels = sorted(set(self.labels[k]) | set(self.labels[l]))
      lowest = min([max([0] + [cl for ls, lt, c, cl in terms
        if ls == f[self.labels[l].index(lt)]])
        for f in itertools.permutations(labels, len(self.labels[l]))])
      if (clash - lowest > threshold):
        result.append((k, l))
    return result

  def key(self, assignment):
    changed = tuple([int(o != 0) for o in assignment])
    return (round(self.score(assignment)[0], 6), sum(changed), changed)

  def best(self, max_combinations, min_gain, fixed=()):
    """Lowest-cost assignment; sites in fixed, and sites whose relabeling
    lowers the cost by min_gain or less, are left unchanged"""
    choices = [[0] if k in fixed else range(len(o))
      for k, o in enumerate(self.options)]
    n_combinations = 1
    for c in choices: n_combinations *= len(c)
    if (n_combinations <= max_combinations):
      assignment = list(min(itertools.product(*choices), key=self.key))
    else:
      assignment = self.steepest_descent(fixed=fixed,
        max_combinations=max_combinations)
    # undo insignificant changes, least significant first
    while True:
      gains = [(self.local_cost(k, 0, assignment) -
                self.local_cost(k, o, assignment), k)
        for k, o in enumerate(assignment) if o != 0]
      if (len(gains) == 0): break
      gain, k = min(gains)
      if (gain > min_gain): break
      assignment[k] = 0
    return tuple(assignment)

  def steepest_descent(self, fixed, max_combinations):
    """Relabel one site at a time, the move that lowers the cost most first
    (later site on ties, as in key)"""
    n = len(self.sites)
    assignment = [0]*n
    def best_move(k):
      if (len(self.options[k]) == 1 or k in fixed): return (0, 0)
      current = self.local_cost(k, assignment[k], assignment)
      return min([
        (round(self.local_cost(k, o, assignment) - current, 6), o)
        for o in range(len(self.options[k])) if o != assignment[k]])
    moves = [best_move(k) for k in range(n)]
    while True:
      k = min(range(n), key=lambda k: (moves[k][0], -k))
      if (moves[k][0] >= 0): break
      assignment[k] = moves[k][1]
      for l in [k] + [l for l, terms, side in self.incident[k]]:
        moves[l] = best_move(l)
    # Single moves cannot fix neighbors that are wrong together (e.g. two
    # linked residues both mislabeled). Relabel jointly connected sites of
    # pairs not at their best pairing, or, if too many, each such site with
    # its neighbors.
    for cycle in range(50):
      changed = False
      strained = self.strained(assignment)
      in_strained = set([k for pair in strained for k in pair])
      for cluster in _clusters([(k, l) for k, l, terms in self.edges
          if k in in_strained and l in in_strained]):
        cluster = [k for k in cluster if k not in fixed]
        n_combinations = 1
        for k in cluster: n_combinations *= len(self.options[k])
        if (n_combinations <=
            max(max_combinations, max_block_combinations)):
          blocks = [cluster]
        else:
          blocks = []
          for k in cluster:
            block = [k]
            n_combinations = len(self.options[k])
            neighbors = sorted(set([l for l, t, s in self.incident[k]]),
              key=lambda l: l not in in_strained)
            for l in neighbors:
              if (l in fixed or len(self.options[l]) == 1 or
                  n_combinations * len(self.options[l]) >
                    max_block_combinations):
                continue
              block.append(l)
              n_combinations *= len(self.options[l])
            blocks.append(block)
        for block in blocks:
          if (self.optimize_block(block, assignment)): changed = True
      if (not changed): break
    return assignment

  def optimize_block(self, block, assignment):
    """Best options for sites in block, others fixed; True if changed"""
    current = [assignment[b] for b in block]
    def block_key(options):
      trial = assignment[:]
      for b, o in zip(block, options): trial[b] = o
      return (round(self.block_cost(block, trial), 6),
        sum([int(o != c) for o, c in zip(options, current)]))
    best = min(itertools.product(
      *[range(len(self.options[b])) for b in block]), key=block_key)
    if (block_key(best)[0] < block_key(current)[0]):
      for b, o in zip(block, best): assignment[b] = o
      return True
    return False

  def block_cost(self, block, assignment):
    """Cost of the edges with at least one site in block"""
    maps = dict([(k, self.mapping(k, assignment[k])) for k in block])
    result, seen = 0, set()
    for k in block:
      for l, terms, side in self.incident[k]:
        if (id(terms) in seen): continue
        seen.add(id(terms))
        if (l not in maps): maps[l] = self.mapping(l, assignment[l])
        for term in terms:
          if (maps[k][term[side]] == maps[l][term[1-side]]):
            result += term[2]
    return result

  def strained(self, assignment):
    """Pairs of sites whose current pairing is not their lowest-cost one"""
    result = []
    for k, l, terms in self.edges:
      def cost(ok, ol):
        mk, ml = self.mapping(k, ok), self.mapping(l, ol)
        return sum([c for ls, lt, c, cl in terms if mk[ls] == ml[lt]])
      lowest = min([cost(ok, ol) for ok in range(len(self.options[k]))
        for ol in range(len(self.options[l]))])
      if (cost(assignment[k], assignment[l]) > lowest + 1.e-6):
        result.append((k, l))
    return result

def run(hierarchy,
        crystal_symmetry,
        distance_cutoff=4.5,
        clash_tolerance=0.5,
        bond_weight=1.0,
        include_hydrogens=False,
        check_only=False,
        max_combinations=10000,
        log=None):
  if (log is None): log = null_out()
  n_sites, n_components, swaps, unresolved = 0, 0, [], []
  for model in hierarchy.models():
    atoms = _atoms(model=model, include_hydrogens=include_hydrogens)
    n_sites += len(atoms.sites)
    if (len(atoms.sites) < 2): continue
    symmetry = _symmetry(crystal_symmetry, atoms.sites_cart, distance_cutoff)
    edges = _edge_terms(atoms=atoms, symmetry=symmetry,
      distance_cutoff=distance_cutoff, clash_tolerance=clash_tolerance)
    for sites in _components(len(atoms.sites), edges):
      n_components += 1
      component = _component(sites=sites, labels=atoms.labels, edges=edges,
        bond_weight=bond_weight)
      best = component.best(max_combinations=max_combinations,
        min_gain=significant_clash)
      frustrated = component.frustrated(best, threshold=unresolved_clash)
      if (len(frustrated) > 0):
        # leave frustrated sites and their neighbors unchanged, relabel the
        # rest of the group consistently with them
        held = set([k for pair in frustrated for k in pair])
        held |= set([l for k in held for l, terms, side in
          component.incident[k]])
        best = component.best(max_combinations=max_combinations,
          min_gain=significant_clash, fixed=held)
        ids = [atoms.sites[sites[k]].id_str() for k in sorted(held)]
        unresolved.append(ids)
        print("  cannot make altlocs consistent: %s" % ", ".join(
          ['"%s"' % i for i in ids]), file=log)
      for k, option in enumerate(best):
        if (option == 0): continue
        rg = atoms.sites[sites[k]]
        mapping = component.mapping(k, option)
        swaps.append((rg.id_str(), mapping))
        print('  relabel "%s": %s' % (rg.id_str(), " ".join(
          ["%s->%s" % (l, mapping[l]) for l in sorted(mapping)])), file=log)
        if (not check_only):
          rg.relabel_altlocs(mapping)
  return group_args(
    n_sites      = n_sites,
    n_components = n_components,
    swaps        = swaps,
    unresolved   = unresolved)
