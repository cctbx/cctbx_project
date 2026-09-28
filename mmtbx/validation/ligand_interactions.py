"""
Ligand interaction profile from geometric tools (stage 1: H-bonds, clashes, vdW
contacts).

Input: an mmtbx.model.manager with electron-cloud H, as the ligand validation tool
prepares it (this module does not place H), processed with restraints, the
ligand's iselection and its selection string. The selection must be one molecule:
connected by the restraints' bonds (it may span residues, e.g. a glycan; one atom,
e.g. an ion, is connected; alternate conformers of an atom count as connected).
Sources:
  pnp     cctbx process_nonbonded_proxies on the ligand and the residues within
          3 A (ligand_overlaps; what validate_ligands reports): clashes and H-bonds
          with at least one ligand atom
  probe2  library call (mmtbx.programs.probe2, approach=both, ligand selection as
          source, the rest as target, raw output, no files written); per atom
          pair the dots per class (wc, cc, so, hb, bo), the pair class (the first
          class in pair_class_order with dots), the minimum atom-pair gap
Entries (type, subtype, atoms, operators, labels, residue, symop, geometry,
sources, cross_check, model_support), ligand-environment only:
  hbond   one per D-H...A: pnp's H-bonds and probe2's hb pairs
  clash   pnp's clashes and probe2's bo pairs
  vdw     probe2's wc, cc or so pairs (subtype), dot counts per class, minimum gap
cross_check (hbond, clash): "pnp and probe2", "pnp only", "probe2 only" (the last
two listed as disagreements), or "symmetry, probe2 not applicable" (probe2 does not
see symmetry-related copies). Ligand-internal pnp records are listed apart
(internal), not as entries. model_support is left unset (the measure is open).
The donor of an H in a probe2 H-bond is its bonded heavy atom from the restraints'
1-2 connectivity (same conformer); an H without exactly one bonded heavy atom is
listed in unresolved_donors and the donor left unset. probe2's hb geometry adds
d_HA and a_DHA from the model (probe2's hb class is overlap-based, without angles).
Symmetry convention: the selected ligand stays in place; an entry's operators (one
per atom, "x,y,z" for untransformed atoms) move its partner atoms (the
non-ligand side; for the ligand with its own symmetry copy, the acceptor of an
H-bond and the second atom of a clash). symop is the partner's operator. pnp's
operator belongs to its own atom pair (rt_mx_ji); the operator or its inverse is
taken, whichever applied to the partner reproduces pnp's distance. The partner
residue is labelled with the operator ("A 856 301 (x+1,y,z)"), also when it is
the ligand itself.
Contact patches: probe2's ligand -> environment dots, at their location on the
ligand's surface, grouped by spatial connectivity (dots closer than
patch_link_distance are linked); per patch the dot counts per class, the mean dot
location (center), atom pairs, ligand atoms and residues. No patch-level type and no rule for counting patches.
"""
from __future__ import absolute_import, division, print_function
import sys
from six.moves import cStringIO as StringIO
import iotbx.phil
import cctbx.geometry_restraints.process_nonbonded_proxies as pnp
from libtbx import group_args
from libtbx.utils import Sorry, null_out
from scitbx.array_family import flex
from cctbx import sgtbx

probe_classes = ("wc", "cc", "so", "hb", "bo")

master_phil_str = """
ligand_interactions
{
  pair_class_order = bo hb so cc wc
    .type = strings
    .help = "Pair class for an atom pair with dots of several probe2 classes: the \\
first class in this list that has dots."
  patch_link_distance = 0.5
    .type = float
    .help = "Contact patches: ligand -> environment dots closer than this (A) \\
are linked. At probe density 16 dots/A^2 the nearest-neighbour spacing is \\
0.256 A (median; 90th percentile 0.261); 0.5 A bridges the seams between \\
neighbouring atoms' dot sets. 5XH3 ligand 856: 11 patches at 0.4 A, 9 at \\
0.45-0.5 A, 6 at 0.55-0.7 A."
  use_neutron_distances = False
    .type = bool
    .help = "probe2's use_neutron_distances (electron-cloud H is the convention)."
  include scope mmtbx.probe.Helpers.probe_phil_parameters
}
"""

def master_params():
  return iotbx.phil.parse(master_phil_str, process_includes=True)

# ------------------------------------------------------------------------------
# atoms

def atom_key(atom):
  """(chain, resid+icode, resname, name, altloc) of a hierarchy atom, as probe2 names it."""
  ag = atom.parent()
  rg = ag.parent()
  return (rg.parent().id.strip(), str(rg.resseq_as_int()) + rg.icode.strip(),
    ag.resname.strip().upper(), atom.name.strip(), ag.altloc.strip())

def atom_label(atom):
  ag = atom.parent()
  rg = ag.parent()
  alt = ag.altloc.strip()
  return "%s %s %s %s%s" % (rg.parent().id.strip(), ag.resname.strip(),
    (rg.resseq + rg.icode).strip(), atom.name.strip(), (" alt " + alt) if alt else "")

def residue_label(atom):
  ag = atom.parent()
  rg = ag.parent()
  return "%s %s %s" % (rg.parent().id.strip(), ag.resname.strip(),
    (rg.resseq + rg.icode).strip())

# ------------------------------------------------------------------------------
# probe2

def probe_atom_field(atom):
  """The atom field of probe2's raw output for a hierarchy atom (probe2's format)."""
  ag = atom.parent()
  rg = ag.parent()
  icode = rg.icode if rg.icode != "" else " "
  return "{:>2s}{:>4s}{}{:>3s} {:<3s}{:1s}".format(rg.parent().id,
    str(rg.resseq_as_int()), icode, ag.resname.strip().upper(), atom.name, ag.altloc)

def probe_atom_table(atoms):
  """{probe2 atom field: i_seq}; Sorry if two atoms share a field."""
  result = {}
  for a in atoms:
    f = probe_atom_field(a)
    if f in result:
      raise Sorry("probe2 atom field %r is not unique." % f)
    result[f] = a.i_seq
  return result

def probe_key(field):
  """
  atom_key from a probe2 raw atom field, by probe2's columns
  "{:>2s}{:>4s}{}{:>3s} {:<3s}{:1s}": chain 2, resid 4, icode 1, resname (3 or
  more), a space, atom name, altloc 1. A resid wider than 4 characters shifts the
  columns; mapping to i_seqs uses probe_atom_table instead.
  """
  if len(field) < 15:
    raise Sorry("probe2 atom field %r not understood." % field)
  t = field[7:-1].lstrip().split(" ", 1)
  if len(t) != 2:
    raise Sorry("probe2 atom field %r not understood." % field)
  return (field[0:2].strip(), field[2:6].strip() + field[6].strip(), t[0],
    t[1].strip(), field[-1].strip())

def pair_class(dots, order=("bo", "hb", "so", "cc", "wc")):
  """The first class in order with dots (dots: {class: count})."""
  for c in order:
    if dots.get(c, 0):
      return c
  return None

def parse_probe_raw(text, order=("bo", "hb", "so", "cc", "wc")):
  """
  probe2 raw output (':direction:class:source:target:gap:dotGap:spike x:y:z:
  spikeLen:score:srcClass:targetClass:loc x:y:z:srcB:targetB'). Returns
  group_args(pairs, dots): pairs {frozenset((source field, target field)):
  dict(dots={class: n}, min_gap, pair_class)} over both directions (min_gap: the
  atom-pair gap); dots [dict(direction, cls, source, target, loc, spike, gap)]
  (loc: the dot on the source atom's surface; gap: the dot's gap).
  """
  pairs, dots = {}, []
  for line in text.splitlines():
    f = line.split(":")
    if len(f) < 17 or f[2] not in probe_classes:
      continue
    source, target = f[3], f[4]
    p = pairs.setdefault(frozenset([source, target]),
      dict(dots=dict([(c, 0) for c in probe_classes]), min_gap=None))
    p["dots"][f[2]] += 1
    gap = float(f[5])
    p["min_gap"] = gap if p["min_gap"] is None else min(p["min_gap"], gap)
    dots.append(dict(direction=f[1], cls=f[2], source=source, target=target,
      spike=(float(f[7]), float(f[8]), float(f[9])),
      loc=(float(f[14]), float(f[15]), float(f[16])), gap=float(f[6])))
  for p in pairs.values():
    p["pair_class"] = pair_class(p["dots"], order)
  return group_args(pairs=pairs, dots=dots)

def probe2_parameters(source_selection, target_selection, probe=None,
                      use_neutron_distances=False):
  """probe2's master PHIL and params: approach=both, raw, no files; probe scope copied."""
  import iotbx.cli_parser
  from mmtbx.programs import probe2
  parser = iotbx.cli_parser.CCTBXParser(program_class=probe2.Program, logger=null_out())
  parser.parse_args([
    "source_selection=%s" % source_selection,
    "target_selection=%s" % target_selection,
    "approach=both", "output.format=raw", "output.write_files=False",
    "output.filename=probe2_ligand_interactions.txt",
    "use_neutron_distances=%s" % use_neutron_distances])
  params = parser.working_phil.extract()
  if probe is not None:
    for name in [n for n in dir(probe) if not n.startswith("_")]:
      if hasattr(params.probe, name):
        setattr(params.probe, name, getattr(probe, name))
  return parser.master_phil, params

def run_probe2(model, source_selection, target_selection, probe=None,
               use_neutron_distances=False):
  """
  probe2 as a library call on a deep copy of model (run() can add phantom H to
  waters): returns the raw output string. The model must carry its H. Raw, not
  JSON: probe2's approach=both writes no JSON.
  """
  from iotbx.data_manager import DataManager
  from mmtbx.programs import probe2
  master_phil, params = probe2_parameters(source_selection, target_selection, probe,
    use_neutron_distances)
  dm = DataManager(["model"])
  dm.add_model("ligand_interactions_model", model)
  p2 = probe2.Program(dm, params, master_phil=master_phil, logger=null_out())
  p2.overrideModel(model.deep_copy(), processed=False)
  results, output = p2.run()
  return output

# ------------------------------------------------------------------------------
# process_nonbonded_proxies

def ligand_overlaps(model, sel_str, within_radius=3.0):
  """
  Clashes and H-bonds with at least one ligand atom (cctbx
  process_nonbonded_proxies on the ligand and the residues within within_radius).
  Returns group_args: n_clashes, clashscore, n_clashes_sym, clashes_str,
  n_hbonds (validate_ligands' report), clash_records [dict(i_seq, j_seq,
  distance, sum_vdw_radii, overlap, symop)], hbond_records [dict(d_seq, h_seq,
  a_seq, d_HA, d_DA, a_DHA, symop)] (full-model i_seqs), hbond_criteria and
  clash_criteria (the values applied). The model must have H.
  """
  sel_within_str = '%s or (residues_within (%s, %s))' \
    % (sel_str, within_radius, sel_str)
  sel_within = model.selection(sel_within_str)
  model_within = model.select(sel_within)
  isel_ligand_within = model_within.iselection(sel_str)

  processed_nbps = pnp.manager(model = model_within)
  clashes = processed_nbps.get_clashes()
  hbonds = processed_nbps.get_hbonds()

  clashes_dict = clashes._clashes_dict
  hbonds_dict = hbonds._hbonds_dict

  ligand_clashes_dict = {}
  for iseq_tuple, record in clashes_dict.items():
    if (iseq_tuple[0] in isel_ligand_within or
        iseq_tuple[1] in isel_ligand_within):
      ligand_clashes_dict[iseq_tuple] = record

  ligand_clashes = pnp.clashes(
                  clashes_dict = ligand_clashes_dict,
                  model        = model_within)

  ligand_hbonds_dict = {}
  # iseq_tuple is (donor, H, acceptor)
  for iseq_tuple, record in hbonds_dict.items():
    if any(i_seq in isel_ligand_within for i_seq in iseq_tuple):
      ligand_hbonds_dict[iseq_tuple] = record

  ligand_hbonds = pnp.hbonds(
                  hbonds_dict  = ligand_hbonds_dict,
                  model        = model_within)

  results_hbonds = ligand_hbonds.get_results()

  string_io = StringIO()
  ligand_clashes.show(log=string_io, show_clashscore=False,
    show_header=False)

  results = ligand_clashes.get_results()

  # i_seqs from model_within to the full model; symop from the rt_mx (r[4]; r[3]
  # is a display flag)
  isel_within = sel_within.iselection()
  def xyz(rt_mx):
    return "" if rt_mx is None else rt_mx.as_xyz()
  clash_records = [dict(
    i_seq = int(isel_within[i]), j_seq = int(isel_within[j]),
    distance = r[0], sum_vdw_radii = r[1], overlap = r[2], symop = xyz(r[4]))
    for (i, j), r in ligand_clashes_dict.items()]
  hbond_records = [dict(
    d_seq = int(isel_within[d]), h_seq = int(isel_within[h]),
    a_seq = int(isel_within[a]), d_HA = r[0], d_DA = r[1], a_DHA = r[2],
    symop = xyz(r[4]))
    for (d, h, a), r in ligand_hbonds_dict.items()]
  p = processed_nbps
  # a_YAH_cutoff is stored by pnp but not applied
  hbond_criteria = dict(Hs=list(p.Hs), As=list(p.As), Ds=list(p.Ds),
    d_HA_cutoff=list(p.d_HA_cutoff), d_DA_cutoff=list(p.d_DA_cutoff),
    a_DHA_cutoff=p.a_DHA_cutoff, min_bonds_H_A=p.min_bonds_H_A)
  # hard-coded in pnp.manager: model_distance - vdw_sum < -0.40
  clash_criteria = dict(min_overlap=0.4, within_radius=within_radius)
  hbond_criteria["within_radius"] = within_radius

  return group_args(
    n_clashes      = results.n_clashes,
    clashscore     = results.clashscore,
    n_clashes_sym  = results.n_clashes_sym,
    clashes_str    = string_io.getvalue(),
    n_hbonds       = results_hbonds.n_hbonds,
    clash_records  = clash_records,
    hbond_records  = hbond_records,
    hbond_criteria = hbond_criteria,
    clash_criteria = clash_criteria)

def _identity(symop):
  return symop in (None, "", "x,y,z")

# ------------------------------------------------------------------------------
# patches

def dot_patches(xyz, link_distance):
  """Connected components of points (flex.vec3_double) linked below link_distance: list of index lists."""
  n = xyz.size()
  parent = list(range(n))
  def find(i):
    while parent[i] != i:
      parent[i] = parent[parent[i]]
      i = parent[i]
    return i
  cells = {}
  for i, p in enumerate(xyz):
    cells.setdefault(tuple([int(c // link_distance) for c in p]), []).append(i)
  d2 = link_distance * link_distance
  for (cx, cy, cz), members in cells.items():
    for dx in (-1, 0, 1):
      for dy in (-1, 0, 1):
        for dz in (-1, 0, 1):
          other = cells.get((cx + dx, cy + dy, cz + dz))
          if other is None:
            continue
          for i in members:
            for j in other:
              if j <= i:
                continue
              p, q = xyz[i], xyz[j]
              if (p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 + (p[2] - q[2]) ** 2 < d2:
                ri, rj = find(i), find(j)
                if ri != rj:
                  parent[ri] = rj
  groups = {}
  for i in range(n):
    groups.setdefault(find(i), []).append(i)
  return sorted(groups.values(), key=lambda g: (-len(g), g[0]))

# ------------------------------------------------------------------------------

class manager(object):
  """
  The ligand interaction profile (stage 1). model: mmtbx.model.manager with H,
  processed with restraints; ligand_isel: flex.size_t; sel_str: the ligand's
  selection string (selects exactly ligand_isel); params: extract of
  master_phil_str (the ligand_interactions scope), default master values.
  """
  def __init__(self, model, ligand_isel, sel_str, params=None, log=None):
    if params is None:
      params = master_params().extract().ligand_interactions
    self.model = model
    self.ligand_isel = flex.size_t(list(ligand_isel))
    self.sel_str = sel_str
    self.params = params
    self.log = log if log is not None else null_out()
    if not model.has_hd():
      raise Sorry("ligand_interactions needs a model with H (electron-cloud positions).")
    bad = [c for c in params.pair_class_order if c not in probe_classes]
    if bad or sorted(params.pair_class_order) != sorted(probe_classes):
      raise Sorry("pair_class_order must list each of %s once, got %s." % (
        " ".join(probe_classes), " ".join(params.pair_class_order)))
    self.entries = []
    self.internal = []
    self.unresolved_donors = []
    self.patches = []
    self.disagreements = []

  def run(self):
    atoms = self.model.get_hierarchy().atoms()
    self._atoms = atoms
    lig = set(self.ligand_isel)
    self._lig = lig
    self._fsc0 = self.model.get_restraints_manager().geometry.shell_sym_tables[0] \
      .full_simple_connectivity()
    if set(self.model.selection(self.sel_str).iselection()) != lig:
      raise Sorry("sel_str %r does not select the ligand_isel atoms." % self.sel_str)
    self._check_one_molecule()
    table = probe_atom_table(atoms)
    # probe2; by selection, not resname: other copies of the ligand are environment
    self.probe_output = run_probe2(self.model, "(%s)" % self.sel_str, "not (%s)" % self.sel_str,
      probe=self.params.probe, use_neutron_distances=self.params.use_neutron_distances)
    parsed = parse_probe_raw(self.probe_output, self.params.pair_class_order)
    # fields that are not model atoms (e.g. probe2's phantom water H) are skipped
    self.probe_pairs = {}
    self.probe_unmapped = set()
    for k, v in parsed.pairs.items():
      seqs = [table.get(x) for x in k]
      if None in seqs:
        self.probe_unmapped.update([x for x in k if x not in table])
        continue
      a, b = seqs if len(seqs) == 2 else seqs * 2
      i, j = (a, b) if a in lig else (b, a)
      if (i in lig) != (j in lig):
        self.probe_pairs[(i, j)] = v
    self.probe_dots = [dict(d, source=table[d["source"]], target=table[d["target"]])
      for d in parsed.dots if d["source"] in table and d["target"] in table]
    self.overlaps = ligand_overlaps(self.model, self.sel_str)
    self._build_entries()
    self._build_patches()
    return self

  def _check_one_molecule(self):
    """Sorry unless the ligand is one molecule by the restraints' bonds (altlocs of an atom joined)."""
    atoms = self._atoms
    parent = dict([(i, i) for i in self._lig])
    def find(i):
      while parent[i] != i:
        parent[i] = parent[parent[i]]
        i = parent[i]
      return i
    def union(i, j):
      parent[find(i)] = find(j)
    same_atom = {}
    for i in self._lig:
      for j in self._fsc0[i]:
        if j in parent:
          union(i, j)
      rg = atoms[i].parent().parent()
      k = (rg.parent().id, rg.resseq, rg.icode, atoms[i].name)
      if k in same_atom:
        union(i, same_atom[k])
      same_atom[k] = i
    parts = {}
    for i in sorted(self._lig):
      parts.setdefault(find(i), []).append(i)
    if len(parts) > 1:
      raise Sorry("The ligand selection %r is not one molecule: %d fragments not "
        "connected by the restraints' bonds (%s). Select one molecule." % (self.sel_str,
        len(parts), "; ".join([residue_label(atoms[g[0]]) + " " + atoms[g[0]].name.strip()
        for g in sorted(parts.values())])))

  # -- entries -----------------------------------------------------------------

  def _donor_of(self, h):
    """The heavy atom bonded to H (restraints' 1-2 connectivity); None, listed, if not exactly one."""
    heavy = [k for k in self._fsc0[h] if not self._is_h(k)]
    if len(heavy) == 1:
      return heavy[0]
    self.unresolved_donors.append(dict(h=atom_label(self._atoms[h]),
      bonded_heavy=[atom_label(self._atoms[k]) for k in heavy]))
    return None

  def _moved(self, i, rt_mx):
    uc = self.model.crystal_symmetry().unit_cell()
    return uc.orthogonalize(rt_mx * uc.fractionalize(self._atoms[i].xyz))

  def _partner_operator(self, fixed, partner, symop, distance):
    """pnp's operator or its inverse: the one that, applied to partner, gives distance from fixed."""
    rt = sgtbx.rt_mx(symop)
    best = None
    for c in (rt, rt.inverse()):
      x = self._moved(partner, c)
      err = abs(self._atoms[fixed].distance(x) - distance)
      if best is None or err < best[0]:
        best = (err, c.as_xyz())
    if best[0] > 1.e-3:
      raise Sorry("symmetry operator %s does not reproduce the distance %.3f A of %s ... %s."
        % (symop, distance, atom_label(self._atoms[fixed]), atom_label(self._atoms[partner])))
    return best[1]

  def _entry(self, type_, subtype, atom_seqs, geometry, sources, operators=None,
             cross_check=None):
    """operators: one per atom (None: all "x,y,z"); non-identity ones mark the partner."""
    atoms = self._atoms
    lig = self._lig
    if operators is None:
      operators = ["x,y,z"] * len(atom_seqs)
    moved = [k for k, i in enumerate(atom_seqs) if i is not None and
      not _identity(operators[k])]
    partner = moved or [k for k, i in enumerate(atom_seqs) if i is not None and
      i not in lig]
    ligand_atoms = [i for k, i in enumerate(atom_seqs) if i is not None and i in lig
      and k not in moved]
    def op(k):
      return "" if k not in moved else " (%s)" % operators[k]
    residue = None
    if partner:
      residue = residue_label(atoms[atom_seqs[partner[0]]]) + op(partner[0])
    return dict(type=type_, subtype=subtype,
      atoms=[None if i is None else int(i) for i in atom_seqs],
      operators=list(operators),
      labels=[None if i is None else atom_label(atoms[i]) + op(k)
        for k, i in enumerate(atom_seqs)],
      ligand_atoms=[atom_label(atoms[i]) for i in ligand_atoms],
      residue=residue, symop=(operators[moved[0]] if moved else None),
      geometry=geometry, sources=sorted(sources),
      cross_check=cross_check, model_support=None)

  def _is_h(self, i):
    return self._atoms[i].element.strip().upper() in ("H", "D")

  def _add_checked(self, type_, atom_seqs, e, operators, i, j):
    """hbond/clash entry with its pnp-probe2 cross-check; disagreements listed."""
    if operators is not None:
      check = "symmetry, probe2 not applicable"
    elif e["sources"] == set(["pnp", "probe2"]):
      check = "pnp and probe2"
    else:
      check = "%s only" % list(e["sources"])[0]
    entry = self._entry(type_, None, atom_seqs, e["geometry"], e["sources"],
      operators, check)
    self.entries.append(entry)
    if check.endswith(" only"):
      self.disagreements.append(dict(type=type_, labels=entry["labels"],
        sources=entry["sources"], missing=sorted(set(["pnp", "probe2"]) - e["sources"]),
        probe_class=self._probe_class_of(i, j)))

  def _build_entries(self):
    order = self.params.pair_class_order
    lig = self._lig
    atoms = self._atoms
    # H-bonds: key (H, A, operators of D, H, A)
    hb = {}
    for r in self.overlaps.hbond_records:
      d, h, a = r["d_seq"], r["h_seq"], r["a_seq"]
      g = dict(d_HA=r["d_HA"], d_DA=r["d_DA"], a_DHA=r["a_DHA"])
      if _identity(r["symop"]):
        ops = ""
        if d in lig and h in lig and a in lig:
          self.internal.append(dict(type="hbond", labels=[atom_label(atoms[k])
            for k in (d, h, a)], geometry=dict(pnp=g)))
          continue
      elif h in lig:
        op = self._partner_operator(h, a, r["symop"], r["d_HA"])
        ops = ("x,y,z", "x,y,z", op)
      else:
        op = self._partner_operator(a, h, r["symop"], r["d_HA"])
        ops = (op, op, "x,y,z")
      e = hb.setdefault((h, a, ops), dict(d=d, geometry={}, sources=set()))
      e["sources"].add("pnp")
      e["geometry"]["pnp"] = g
    for (i, j), p in self.probe_pairs.items():
      if p["pair_class"] != "hb":
        continue
      g = dict(dots=dict(p["dots"]), min_gap=p["min_gap"])
      if self._is_h(i) != self._is_h(j):
        h, a = (i, j) if self._is_h(i) else (j, i)
        d = self._donor_of(h)
        # from the model: probe2's hb class does not use angles
        g["d_HA"] = atoms[h].distance(atoms[a])
        g["a_DHA"] = None if d is None else atoms[h].angle(atoms[a], atoms[d], deg=True)
      else:
        h, a, d = None, j, i
      e = hb.setdefault((h, a, ""), dict(d=d, geometry={}, sources=set()))
      e["sources"].add("probe2")
      e["geometry"]["probe2"] = g
    for (h, a, ops), e in sorted(hb.items(), key=lambda x: (str(x[0][0]), x[0][1], str(x[0][2]))):
      self._add_checked("hbond", [e["d"], h, a], e, list(ops) if ops else None, h, a)
    # clashes: key (ligand atom, partner atom, partner operator); for the ligand
    # with its own copy the lower i_seq stays
    cl = {}
    for r in self.overlaps.clash_records:
      i, j = sorted([r["i_seq"], r["j_seq"]])
      if j in lig and i not in lig:
        i, j = j, i
      g = dict(distance=r["distance"], sum_vdw_radii=r["sum_vdw_radii"],
        overlap=r["overlap"])
      if _identity(r["symop"]):
        op = ""
        if i in lig and j in lig:
          self.internal.append(dict(type="clash", labels=[atom_label(atoms[i]),
            atom_label(atoms[j])], geometry=dict(pnp=g)))
          continue
      else:
        op = self._partner_operator(i, j, r["symop"], r["distance"])
      e = cl.setdefault((i, j, op), dict(geometry={}, sources=set()))
      e["sources"].add("pnp")
      e["geometry"]["pnp"] = g
    for (i, j), p in self.probe_pairs.items():
      if p["pair_class"] != "bo":
        continue
      e = cl.setdefault((i, j, ""), dict(geometry={}, sources=set()))
      e["sources"].add("probe2")
      e["geometry"]["probe2"] = dict(dots=dict(p["dots"]), min_gap=p["min_gap"])
    for (i, j, op), e in sorted(cl.items()):
      self._add_checked("clash", [i, j], e, ["x,y,z", op] if op else None, i, j)
    # vdW contacts: probe2's wc, cc, so pairs
    for (i, j), p in sorted(self.probe_pairs.items()):
      if p["pair_class"] in ("wc", "cc", "so"):
        self.entries.append(self._entry("vdw", p["pair_class"], [i, j],
          dict(probe2=dict(dots=dict(p["dots"]), min_gap=p["min_gap"])), ["probe2"]))
    self.pair_class_order = list(order)

  def _probe_class_of(self, i, j):
    if i is None or j is None:
      return None
    k = (i, j) if i in self._lig else (j, i)
    p = self.probe_pairs.get(k)
    return p["pair_class"] if p else None

  # -- patches -----------------------------------------------------------------

  def _build_patches(self):
    atoms = self._atoms
    lig = self._lig
    dots = [d for d in self.probe_dots if d["source"] in lig and d["target"] not in lig]
    xyz = flex.vec3_double([d["loc"] for d in dots])
    for group in dot_patches(xyz, self.params.patch_link_distance):
      counts = dict([(c, 0) for c in probe_classes])
      pairs, lig_atoms, residues = {}, set(), set()
      for k in group:
        d = dots[k]
        counts[d["cls"]] += 1
        i, j = d["source"], d["target"]
        pairs[(i, j)] = pairs.get((i, j), 0) + 1
        lig_atoms.add(atom_label(atoms[i]))
        residues.add(residue_label(atoms[j]))
      center = xyz.select(flex.size_t(group)).mean()
      self.patches.append(dict(n_dots=len(group), dots=counts, center=center,
        pairs=[dict(ligand=atom_label(atoms[i]), environment=atom_label(atoms[j]), dots=n)
          for (i, j), n in sorted(pairs.items(), key=lambda x: -x[1])],
        ligand_atoms=sorted(lig_atoms), residues=sorted(residues)))

  # -- summaries ---------------------------------------------------------------

  def counts(self):
    """Counts per type (vdW per subtype), per residue and per ligand atom."""
    per_type, per_residue, per_atom = {}, {}, {}
    for e in self.entries:
      t = e["type"] if e["subtype"] is None else "%s:%s" % (e["type"], e["subtype"])
      per_type[t] = per_type.get(t, 0) + 1
      r = e["residue"]
      per_residue.setdefault(r, {})
      per_residue[r][t] = per_residue[r].get(t, 0) + 1
      for a in e["ligand_atoms"]:
        per_atom.setdefault(a, {})
        per_atom[a][t] = per_atom[a].get(t, 0) + 1
    return group_args(per_type=per_type, per_residue=per_residue, per_ligand_atom=per_atom)

  def probe_parameters(self):
    p = self.params.probe
    return dict([(n, getattr(p, n)) for n in dir(p) if not n.startswith("_")])

  def as_dict(self):
    c = self.counts()
    return dict(ligand=self.sel_str, pair_class_order=list(self.params.pair_class_order),
      patch_link_distance=self.params.patch_link_distance,
      use_neutron_distances=self.params.use_neutron_distances,
      probe=self.probe_parameters(), hbond_criteria=self.overlaps.hbond_criteria,
      clash_criteria=self.overlaps.clash_criteria, entries=self.entries,
      counts=dict(per_type=c.per_type, per_residue=c.per_residue,
        per_ligand_atom=c.per_ligand_atom),
      disagreements=self.disagreements, internal=self.internal, patches=self.patches,
      unresolved_donors=self.unresolved_donors,
      probe_unmapped=sorted(self.probe_unmapped))

  def show(self, log=None):
    if log is None:
      log = sys.stdout
    c = self.counts()
    print("Ligand interactions: %s" % self.sel_str, file=log)
    print("  pair class order: %s; patch link distance %.2f A" % (
      " ".join(self.params.pair_class_order), self.params.patch_link_distance), file=log)
    print("  pnp H-bond criteria: %s" % ", ".join(["%s=%s" % (k, v) for k, v in
      sorted(self.overlaps.hbond_criteria.items())]), file=log)
    print("  counts: %s" % ", ".join(["%s %d" % (k, v) for k, v in sorted(c.per_type.items())]),
      file=log)
    for e in self.entries:
      g = e["geometry"]
      detail = []
      for s in sorted(g):
        detail.append("%s(%s)" % (s, ", ".join(["%s=%s" % (k, ("%.2f" % v)
          if isinstance(v, float) else v) for k, v in sorted(g[s].items()) if k != "dots"])))
      print("  %-6s %-3s %s  [%s] %s" % (e["type"], e["subtype"] or "",
        " ... ".join([l for l in e["labels"] if l]),
        e["cross_check"] or ", ".join(e["sources"]), " ".join(detail)), file=log)
    if self.disagreements:
      print("  disagreements:", file=log)
      for d in self.disagreements:
        print("    %s %s: reported by %s, not by %s (probe class %s)" % (d["type"],
          " ... ".join([l for l in d["labels"] if l]), ", ".join(d["sources"]),
          ", ".join(d["missing"]), d["probe_class"]), file=log)
    if self.unresolved_donors:
      print("  H without exactly one bonded heavy atom (donor unset):", file=log)
      for d in self.unresolved_donors:
        print("    %s: %s" % (d["h"], ", ".join(d["bonded_heavy"]) or "none"), file=log)
    if self.internal:
      print("  ligand-internal (not entries):", file=log)
      for d in self.internal:
        print("    %s %s" % (d["type"], " ... ".join(d["labels"])), file=log)
    print("  contact patches (ligand -> environment dots):", file=log)
    for k, p in enumerate(self.patches):
      print("    %d: %d dots (%s); ligand %s; residues %s" % (k + 1, p["n_dots"],
        " ".join(["%s %d" % (c, p["dots"][c]) for c in probe_classes if p["dots"][c]]),
        " ".join([a.split()[3] for a in p["ligand_atoms"]]), ", ".join(p["residues"])),
        file=log)
