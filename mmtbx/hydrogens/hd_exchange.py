"""Exchangeable H/D: atom and site identity, H parent table.

Pure hierarchy functions except h_parent_table(), which needs a model manager
with restraints. Used by reduce2 (neutron_option): input check and
normalization, output model and output check.
"""
from __future__ import absolute_import, division, print_function
from libtbx import group_args
from scitbx import matrix
from libtbx.utils import Sorry
import iotbx.pdb
from cctbx import sgtbx
from cctbx.eltbx import covalent_radii
from cctbx.geometry_restraints.auto_linking_types import iterate_covalent
from cctbx.geometry_restraints.linking_class import linking_class

get_class = iotbx.pdb.common_residue_names_get_class

# ------------------------------------------------------------------------------
# Identity

def canonical_resname(resname):
  # normalization renames DOD -> HOH and the output renames it back
  if get_class(name=resname) == "common_water":
    return "HOH"
  return resname

def normalized_name(atom):
  """D atoms only: the first D becomes H (DD21 -> HD21, 1DB -> 1HB)"""
  if atom.element.strip().upper() == "D":
    return atom.name.replace("D", "H", 1)
  return atom.name

def _index_by_memory_id(objects, target):
  for i, o in enumerate(objects):
    if o.memory_id() == target.memory_id():
      return i
  raise RuntimeError("internal error: object not found in its parent")

def residue_part(ag):
  """(model index, chain id, chain occurrence, resseq, icode, residue_group
  occurrence, canonical resname): the residue part of atom_key()/site_key()"""
  rg = ag.parent()
  chain = rg.parent()
  model = chain.parent()
  root = model.parent()
  model_index = _index_by_memory_id(root.models(), model)
  chain_occurrence = 0
  for c in model.chains():
    if c.memory_id() == chain.memory_id(): break
    if c.id == chain.id: chain_occurrence += 1
  rg_occurrence = 0
  for r in chain.residue_groups():
    if r.memory_id() == rg.memory_id(): break
    if r.resseq == rg.resseq and r.icode == rg.icode: rg_occurrence += 1
  return (model_index, chain.id, chain_occurrence, rg.resseq, rg.icode,
          rg_occurrence, canonical_resname(ag.resname))

def residue_parts(hierarchy):
  """residue part of every atom_group in one pass: {memory_id: tuple}"""
  result = {}
  for model_index, model in enumerate(hierarchy.models()):
    chain_count = {}
    for chain in model.chains():
      chain_occurrence = chain_count.get(chain.id, 0)
      chain_count[chain.id] = chain_occurrence + 1
      rg_count = {}
      for rg in chain.residue_groups():
        rid = (rg.resseq, rg.icode)
        rg_occurrence = rg_count.get(rid, 0)
        rg_count[rid] = rg_occurrence + 1
        for ag in rg.atom_groups():
          result[ag.memory_id()] = (model_index, chain.id, chain_occurrence,
            rg.resseq, rg.icode, rg_occurrence, canonical_resname(ag.resname))
  return result

# the keyword residue_part of atom_key()/site_key() shadows the function there
_residue_part, _residue_parts = residue_part, residue_parts

def atom_key(atom, residue_part=None):
  """residue part + own altloc + name as written + element; residue_part
  (from residue_parts()) saves the scan of preceding residues"""
  ag = atom.parent()
  if residue_part is None: residue_part = _residue_part(ag)
  return residue_part + (ag.altloc, atom.name, atom.element.strip().upper())

def site_key(atom, parent_altloc, residue_part=None):
  """residue part + parent's altloc + normalized name; for H or D of a
  validated pair (both give the same key) or a single"""
  if residue_part is None: residue_part = _residue_part(atom.parent())
  return residue_part + (parent_altloc, normalized_name(atom))

def rekey_records(records, renames):
  """records (site_key -> values) after placement renamed input H.
  renames: (residue part, parent altloc, old name, new name) in the order
  made (place_hydrogens.renamed_existing_h). Each exchanges the two names:
  both recorded -> the records swap; one recorded -> it moves to the other
  name. Returns a new dict."""
  result = dict(records)
  for part, parent_altloc, old, new in renames:
    key_old, key_new = part + (parent_altloc, old), part + (parent_altloc, new)
    value_old, value_new = result.pop(key_old, None), result.pop(key_new, None)
    if value_old is not None: result[key_new] = value_old
    if value_new is not None: result[key_old] = value_new
  return result

def atom_keys(hierarchy):
  """atom_key() of every atom in one pass, in hierarchy.atoms() order; a
  duplicate key is an internal error naming both atoms"""
  parts = _residue_parts(hierarchy)
  result = []
  seen = {}
  for atom in hierarchy.atoms():
    key = atom_key(atom, parts[atom.parent().memory_id()])
    if key in seen:
      raise Sorry("Internal error: duplicate atom identity:\n  %s\n  %s"
        % (seen[key], atom.id_str()))
    seen[key] = atom.id_str()
    result.append(key)
  return result

# ------------------------------------------------------------------------------
# Covalent partners and parent table

_EDITS_SLACK = 0.4

def _is_hd(atom):
  return atom.element.strip().upper() in ("H", "D")

def _covalent_radius(element):
  e = element.strip().upper()
  if e == "D": e = "H"
  try:
    return covalent_radii.table(e).radius()
  except (RuntimeError, ValueError):
    return None

def _edits_is_covalent(atoms, i, j, distance_ideal):
  """edits bonds count only up to r_cov(X) + r_cov(H) + 0.4 A"""
  r_i = _covalent_radius(atoms[i].element)
  r_j = _covalent_radius(atoms[j].element)
  if r_i is None or r_j is None: return False
  return distance_ideal <= r_i + r_j + _EDITS_SLACK

def covalent_h_partners(model):
  """{i_seq of H/D: {partner i_seq: [group_args(origin, distance_ideal,
  sym_op)]}} from the covalent bond proxies, simple and asu.
  sym_op is the operator that, applied to the partner (parent), gives the
  partner image bonded to this H/D; x,y,z for simple proxies. rt_mx_ji maps
  atom j into atom i's frame, so it is stored as is when the H/D is i and
  inverted when the H/D is j. Every asu proxy is included: j_sym indexes the
  j atom's mappings, so j_sym == 0 can still be a bond to a symmetry mate."""
  grm = model.get_restraints_manager()
  if grm is None:
    raise Sorry("Internal error: covalent_h_partners() needs restraints.")
  geometry = grm.geometry
  atoms = model.get_hierarchy().atoms()
  assert geometry.bond_params_table.size() == atoms.size()
  origin_ids = linking_class()
  metal = origin_ids.get_origin_id("metal coordination")
  edits = origin_ids.get_origin_id("edits")
  covalent = set(iterate_covalent()) - set([metal])
  simple, asu = geometry.get_all_bond_proxies(
    sites_cart=model.get_sites_cart())
  asu_mappings = None
  if asu.size() > 0:
    asu_mappings = geometry.pair_proxies().bond_proxies.asu_mappings()
  result = {}
  def add(i, j, proxy, rt_mx_ji):
    if proxy.origin_id not in covalent: return
    if proxy.origin_id == edits and not _edits_is_covalent(
        atoms, i, j, proxy.distance_ideal):
      return
    for h, partner, rt_mx in ((i, j, rt_mx_ji), (j, i, rt_mx_ji.inverse())):
      if not _is_hd(atoms[h]): continue
      info = group_args(
        origin         = origin_ids.get_origin_key(proxy.origin_id),
        distance_ideal = proxy.distance_ideal,
        sym_op         = str(rt_mx))
      result.setdefault(h, {}).setdefault(partner, []).append(info)
  for p in simple:
    i, j = p.i_seqs
    add(i, j, p, sgtbx.rt_mx())
  for p in asu:
    add(p.i_seq, p.j_seq, p, asu_mappings.get_rt_mx_ji(pair=p))
  return result

def h_parent_table(model):
  """{atom_key of each H/D: group_args(h_id, parent_id, parent_element,
  parent_occ, parent_altloc, exchangeable, status, source, parent_key,
  parent_sym_ops)}
  status: ok (one covalent partner), floating (none), ambiguous (several);
  source: the covalent restraints found, as (partner id, origin, ideal,
  sym_op); parent_key and parent_sym_ops (sorted operators of the parent
  images bonded to this H/D) identify the parent for pairing"""
  hierarchy = model.get_hierarchy()
  atoms = hierarchy.atoms()
  keys = atom_keys(hierarchy)
  partners = covalent_h_partners(model)
  result = {}
  for i, atom in enumerate(atoms):
    if not _is_hd(atom): continue
    found = partners.get(i, {})
    source = []
    for j in sorted(found):
      for info in found[j]:
        source.append((atoms[j].id_str(), info.origin, info.distance_ideal,
                       info.sym_op))
    entry = group_args(
      h_id           = atom.id_str(),
      parent_id      = None,
      parent_element = None,
      parent_occ     = None,
      parent_altloc  = None,
      exchangeable   = None,
      status         = None,
      source         = source,
      parent_key     = None,
      parent_sym_ops = None)
    if len(found) == 0:
      entry.status = "floating"
    elif len(found) > 1:
      entry.status = "ambiguous"
    else:
      j = list(found)[0]
      parent = atoms[j]
      entry.status         = "ok"
      entry.parent_key     = keys[j]
      entry.parent_sym_ops = tuple(sorted(set(
        info.sym_op for info in found[j])))
      entry.parent_id      = parent.id_str()
      entry.parent_element = parent.element.strip().upper()
      entry.parent_occ     = parent.occ
      entry.parent_altloc  = parent.parent().altloc
      entry.exchangeable   = entry.parent_element != "C"
    result[keys[i]] = entry
  return result

def _water_o(ag):
  """the O of a water H/D: in its own atom_group, else in the blank one"""
  groups = [ag]
  for other in ag.parent().atom_groups():
    if other.altloc.strip() == "" and _is_water(other):
      groups.append(other)
  for group in groups:
    atoms = group.atoms()
    for k in range(atoms.size()):
      if atoms[k].element.strip().upper() == "O":
        return atoms[k]
  return None

def add_water_h_parents(parent_table, hierarchy, oh_length):
  """Add the water H/D missing from a parent table, each bonded to its O.

  For water H placed after the table was made from the restraints, as reduce2
  does: the entries are those h_parent_table() gives a water H/D bonded to its
  O. The source records the bond as placed, (O id, "water H placement",
  oh_length, "x,y,z").

  Parameters
  ----------
  parent_table : dict
      From h_parent_table(); changed in place. Entries already in it stay.
  hierarchy : iotbx.pdb.hierarchy.root
      The hierarchy the table belongs to, with the water H/D.
  oh_length : float
      O-H length the water H were placed at (A).

  Returns
  -------
  int
      Number of entries added.
  """
  parts = _residue_parts(hierarchy)
  atoms = hierarchy.atoms()
  n_added = 0
  for k in range(atoms.size()):
    atom = atoms[k]
    ag = atom.parent()
    if not _is_hd(atom) or not _is_water(ag): continue
    key = atom_key(atom, parts[ag.memory_id()])
    if key in parent_table: continue
    o = _water_o(ag)
    if o is None: continue
    parent_table[key] = group_args(
      h_id           = atom.id_str(),
      parent_id      = o.id_str(),
      parent_element = "O",
      parent_occ     = o.occ,
      parent_altloc  = o.parent().altloc,
      exchangeable   = True,
      status         = "ok",
      source         = [(o.id_str(), "water H placement", oh_length,
                         "x,y,z")],
      parent_key     = atom_key(o, parts[o.parent().memory_id()]),
      parent_sym_ops = ("x,y,z",))
    n_added += 1
  return n_added

# ------------------------------------------------------------------------------
# Candidates and collapse check

_NAME_MISMATCH_MARGIN = 0.1
# distinct H sites on one parent are >= ~1.34 A apart
_SAME_SITE_RADIUS = 0.5

def _distance(a, b):
  return (matrix.col(a.xyz) - matrix.col(b.xyz)).length()

def _hd_only(ag):
  # by index: af_shared_atom has no __iter__ (each loop would end in a C++
  # IndexError)
  atoms = ag.atoms()
  return all(_is_hd(atoms[k]) for k in range(atoms.size()))

def find_hd_candidates(hierarchy, parent_table, provenance=None, model=None):
  """Every H/D as a candidate pair, a single or unresolved, plus same-name H/D
  on different parents and H/D pairs in one conformer (one of the two is
  floating there: pdb_interpretation restrains only one). No distance drops a
  name-matched pair; collapse_check() decides.
  Parent lookup: parent_table[atom_key] or, with provenance (output check),
  parent_table[provenance[atom_key]]. A pair is one H and one D in the same
  residue_group and resname with the same parent atom and sym_op, matched
  one-to-one by (D-normalized) name, then checked geometrically (another
  opposite isotope within 0.5 A and > 0.1 A closer than the name partner).
  Leftovers within 0.5 A of each other are matched by position; farther apart
  they are singles (one isotope per site on a multi-H parent). Both kinds of
  geometric match are flagged name_mismatch.
  Atoms are given by label (h_id, d_id, geometric_partner, atom_id) and by
  index into hierarchy.atoms() (i_h, i_d, i_geometric_partner, i_seq):
  geometric_partner is the near D, the near H, or for a leftover matched by
  position the pair's own D (both None without a name mismatch)."""
  atoms = hierarchy.atoms()
  keys = atom_keys(hierarchy)
  rotatable = None
  if model is not None:
    # fails when H and D of one site share a conformer (duplicate names in
    # rotamer_utils.tardy_model); such input is not collapsible anyway
    try:
      sel = model.rotatable_hd_selection(iselection=False)
    except AssertionError:
      sel = None
    if sel is not None and sel.size() == atoms.size(): rotatable = sel
  # entries per H/D
  groups = {}       # (residue part, parent key, sym ops) -> ([H], [D])
  by_name = {}      # (residue part, normalized name) -> [(atom, entry)]
  unresolved = []
  info = {}
  for i, atom in enumerate(atoms):
    if not _is_hd(atom): continue
    key = keys[i]
    table_key = key
    if provenance is not None:
      if key not in provenance:
        raise Sorry("Internal error: no provenance for %s" % atom.id_str())
      table_key = provenance[key]
    if table_key not in parent_table:
      raise Sorry("Internal error: %s not in the H parent table" %
        atom.id_str())
    e = parent_table[table_key]
    info[i] = (key, e)
    if e.status != "ok":
      continue
    residue = key[:7]
    g = groups.setdefault((residue, e.parent_key, e.parent_sym_ops), ([], []))
    g[atom.element.strip().upper() == "D"].append(i)
    by_name.setdefault((residue, normalized_name(atom)), []).append(i)
  pairs = []
  singles = []
  # per atom_group, once: H/D only, and its atom names with their counts
  hd_only = {}
  def ag_hd_only(ag):
    k = ag.memory_id()
    if k not in hd_only: hd_only[k] = _hd_only(ag)
    return hd_only[k]
  ag_names = {}
  def ag_name_count(ag, name):
    k = ag.memory_id()
    if k not in ag_names:
      counts = {}
      ag_atoms = ag.atoms()
      for j in range(ag_atoms.size()):
        n = ag_atoms[j].name
        counts[n] = counts.get(n, 0) + 1
      ag_names[k] = counts
    return ag_names[k].get(name, 0)
  def candidate(ih, id_, name_mismatch, i_geometric, n_competing):
    h, d = atoms[ih], atoms[id_]
    e = info[ih][1]
    occ_sum = h.occ + d.occ
    uij_diff = None
    if h.uij_is_defined() and d.uij_is_defined():
      uij_diff = max(abs(x - y) for x, y in zip(h.uij, d.uij))
    return group_args(
      h_key             = info[ih][0],
      d_key             = info[id_][0],
      h_id              = h.id_str(),
      d_id              = d.id_str(),
      i_h               = ih,
      i_d               = id_,
      parent_id         = e.parent_id,
      parent_occ        = e.parent_occ,
      parent_altloc     = e.parent_altloc,
      distance          = _distance(h, d),
      name_mismatch     = name_mismatch,
      geometric_partner = (None if i_geometric is None else
                           atoms[i_geometric].id_str()),
      i_geometric_partner = i_geometric,
      n_competing       = n_competing,
      layout_ok         = (e.parent_altloc == "" and
                           h.parent().memory_id() != d.parent().memory_id() and
                           h.parent().altloc != "" and
                           ag_hd_only(h.parent()) and
                           d.parent().altloc != "" and
                           ag_hd_only(d.parent())),
      occ_sum           = occ_sum,
      b_diff            = h.b - d.b,
      uij_diff          = uij_diff,
      rotatable         = None if rotatable is None else
                          bool(rotatable[ih] or rotatable[id_]))
  for gkey in sorted(groups):
    hs, ds = groups[gkey]
    used = set()
    # one-to-one by name
    names = {}
    for i in hs + ds:
      names.setdefault(normalized_name(atoms[i]), ([], []))[
        atoms[i].element.strip().upper() == "D"].append(i)
    for name in sorted(names):
      nh, nd = names[name]
      if not nh or not nd: continue
      n_competing = len(nh) + len(nd) - 2
      for ih in nh:
        free = [x for x in nd if x not in used]
        if not free: break
        id_ = min(free, key=lambda x: _distance(atoms[ih], atoms[x]))
        used.update([ih, id_])
        # nearest opposite isotope on the same parent, from either side
        d_name = _distance(atoms[ih], atoms[id_])
        geometric = None
        # (another opposite isotope on the same site, closer than the partner)
        near_d = min(ds, key=lambda x: _distance(atoms[ih], atoms[x]))
        near_h = min(hs, key=lambda x: _distance(atoms[x], atoms[id_]))
        d_near_d = _distance(atoms[ih], atoms[near_d])
        d_near_h = _distance(atoms[near_h], atoms[id_])
        if near_d != id_ and d_near_d <= _SAME_SITE_RADIUS and \
           d_near_d < d_name - _NAME_MISMATCH_MARGIN:
          geometric = near_d
        elif near_h != ih and d_near_h <= _SAME_SITE_RADIUS and \
           d_near_h < d_name - _NAME_MISMATCH_MARGIN:
          geometric = near_h
        pairs.append(candidate(ih, id_, geometric is not None, geometric,
                               n_competing))
    # leftovers by position, same site only: one isotope per site on a
    # multi-H parent (water O/H1/D2, Asn HD21 + DD22) is two singles
    lh = [i for i in hs if i not in used]
    ld = [i for i in ds if i not in used]
    while lh and ld:
      ih, id_ = min(((a, b) for a in lh for b in ld),
                    key=lambda p: _distance(atoms[p[0]], atoms[p[1]]))
      if _distance(atoms[ih], atoms[id_]) > _SAME_SITE_RADIUS: break
      pairs.append(candidate(ih, id_, True, id_, 0))
      lh.remove(ih)
      ld.remove(id_)
    for i in lh + ld:
      atom = atoms[i]
      collision = False
      if atom.element.strip().upper() == "D":
        nn = normalized_name(atom)
        # another atom of its atom_group has that name
        collision = ag_name_count(atom.parent(), nn) - \
          (1 if atom.name == nn else 0) > 0
      singles.append(group_args(atom_key=info[i][0], atom_id=atom.id_str(),
        element=atom.element.strip().upper(), i_seq=i,
        parent_id=info[i][1].parent_id, name_collision=collision))
  # same-name H/D on different parents
  different_parents = []
  for (residue, name), members in sorted(by_name.items()):
    hs = [i for i in members if atoms[i].element.strip().upper() == "H"]
    ds = [i for i in members if atoms[i].element.strip().upper() == "D"]
    for ih in hs:
      for id_ in ds:
        eh, ed = info[ih][1], info[id_][1]
        if (eh.parent_key, eh.parent_sym_ops) == \
           (ed.parent_key, ed.parent_sym_ops):
          continue
        different_parents.append(group_args(
          h_id=atoms[ih].id_str(), d_id=atoms[id_].id_str(),
          h_key=keys[ih], d_key=keys[id_], i_h=ih, i_d=id_,
          h_parent=(eh.parent_id, eh.parent_sym_ops),
          d_parent=(ed.parent_id, ed.parent_sym_ops),
          distance=_distance(atoms[ih], atoms[id_])))
  # floating or ambiguous. A floating H/D with a bonded same-name opposite
  # isotope in a compatible atom_group (same altloc, or either blank; same
  # residue_group and resname) is an H/D pair in one conformer:
  # pdb_interpretation restrains only one of the two
  same_conformer = []
  index = None
  for i in sorted(info):
    key, e = info[i]
    if e.status == "ok": continue
    atom = atoms[i]
    partner = None
    if e.status == "floating":
      if index is None:
        index = dict((a.memory_id(), k) for k, a in enumerate(atoms))
      element = atom.element.strip().upper()
      ag = atom.parent()
      for other_ag in ag.parent().atom_groups():
        if other_ag.resname != ag.resname: continue
        if not (other_ag.altloc == ag.altloc or other_ag.altloc == "" or
                ag.altloc == ""):
          continue
        other_atoms = other_ag.atoms()
        for k in range(other_atoms.size()):
          other = other_atoms[k]
          j = index[other.memory_id()]
          if (j != i and j in info and info[j][1].status == "ok" and
              other.element.strip().upper() != element and
              normalized_name(other) == normalized_name(atom)):
            partner = other
            break
        if partner is not None: break
    if partner is not None:
      same_conformer.append(group_args(atom_key=key, atom_id=atom.id_str(),
        i_seq=i, partner_id=partner.id_str()))
    else:
      unresolved.append(group_args(atom_key=key, atom_id=atom.id_str(),
        i_seq=i, status=e.status, source=e.source))
  return group_args(pairs=pairs, singles=singles, unresolved=unresolved,
    different_parents=different_parents, same_conformer=same_conformer)

_LAYOUT_TEXT = "H and D not in H/D-only altloc groups"

def _apart_text(distance):
  return "H and D %.2f A apart" % distance

def collapse_check(candidates, xyz_tolerance=0.2, accepted=None):
  """reduce2's policy -> group_args(pairs, problems, unresolved, notes).
  A pair collapses only if it is within xyz_tolerance, has no competing
  partner and no name mismatch, its parent is in the blank altloc and H and D
  are in different H/D-only altloc groups. Problems (stop): any other pair,
  same-name H/D on different parents, a single D whose H name exists in its
  atom_group, an H/D pair in one conformer (one of them floating).
  Unresolved: floating (remove) and ambiguous (stop, restraints named). Notes (never a stop): occupancy sum vs parent occupancy, B.
  accepted: (h_key, d_key) accepted by repair_plan(): pairs that collapse
  beyond xyz_tolerance when the distance is their only problem (case e), and
  same-name H and D on altloc copies of one parent kept as two singles (case
  c, no problem listed); None for none."""
  pairs = []
  problems = []
  notes = []
  for c in candidates.pairs:
    reasons = []
    if c.distance > xyz_tolerance and not (accepted is not None and
        (c.h_key, c.d_key) in accepted):
      reasons.append(_apart_text(c.distance))
    if c.n_competing > 0:
      reasons.append("%d competing partner(s)" % c.n_competing)
    if c.name_mismatch:
      if c.geometric_partner is not None and c.geometric_partner != c.d_id:
        reasons.append("names swapped (coincides with %s)" %
          c.geometric_partner)
      else:
        reasons.append("names do not match")
    if c.parent_altloc != "":
      reasons.append("parent in conformer group '%s'" % c.parent_altloc)
    elif not c.layout_ok:
      reasons.append(_LAYOUT_TEXT)
    if accepted is not None and (c.h_key, c.d_key) in accepted and \
       c.distance > xyz_tolerance and reasons:
      # only the distance was accepted
      reasons.insert(0, _apart_text(c.distance))
    if reasons:
      problems.append(group_args(kind="pair", h_id=c.h_id, d_id=c.d_id,
        reasons=reasons, candidate=c))
    else:
      pairs.append(c)
    if abs(c.occ_sum - c.parent_occ) > 0.01:
      notes.append(group_args(kind="occupancy", h_id=c.h_id, d_id=c.d_id,
        text="H + D occupancy %.2f, parent %.2f%s" % (c.occ_sum,
          c.parent_occ, " (sum 1)" if abs(c.occ_sum - 1) <= 0.01 else ""),
        sum_is_one=abs(c.occ_sum - 1) <= 0.01))
    if abs(c.b_diff) > 0.01:
      notes.append(group_args(kind="b", h_id=c.h_id, d_id=c.d_id,
        text="H and D B differ by %.2f" % c.b_diff))
  for dp in candidates.different_parents:
    if accepted is not None and (dp.h_key, dp.d_key) in accepted: continue
    problems.append(group_args(kind="different parents", h_id=dp.h_id,
      d_id=dp.d_id, reasons=["same name on different parents: %s %s / %s %s"
        % (dp.h_parent + dp.d_parent)], candidate=dp))
  for x in candidates.same_conformer:
    problems.append(group_args(kind="same conformer", h_id=None,
      d_id=None, reasons=["H/D pair in one conformer not supported: %s / %s"
        % (x.atom_id, x.partner_id)], candidate=x))
  # the bonded partner of such a pair is reported there, not as a collision
  explained = set(x.partner_id for x in candidates.same_conformer)
  for s in candidates.singles:
    if s.name_collision and s.atom_id not in explained:
      problems.append(group_args(kind="name collision", h_id=None,
        d_id=s.atom_id, reasons=["renamed D collides with an existing H"],
        candidate=s))
  unresolved = []
  for u in candidates.unresolved:
    action = "remove" if u.status == "floating" else "stop"
    unresolved.append(group_args(atom_id=u.atom_id, status=u.status,
      action=action, source=u.source))
  return group_args(pairs=pairs, problems=problems, unresolved=unresolved,
    notes=notes)

# ------------------------------------------------------------------------------
# Repairs of deposited H/D layouts (B16)

# (a): the parent copies' ideal H positions within this of each other; the
# supplied atom within _REPAIR_ANGLE of each copy's ideal direction, its
# length in [L - below, L + above] at either dictionary length L
_REPAIR_IDEAL_TOLERANCE = 0.10
_REPAIR_ANGLE = 15.0
_REPAIR_BOND_BELOW = 0.20
_REPAIR_BOND_ABOVE = 0.05
# (b), (d): H and D of one site
_REPAIR_PAIR_TOLERANCE = 0.20
# (e): an empirical repair limit, the largest H-D distance on a one-H parent
# in the deposited models examined (7T5E Arg B 66 HE/DE 0.59 A), not a
# geometric bound
REPAIR_E_LIMIT = 0.6
_OCC_TOLERANCE = 0.01

def _label(atom):
  return atom.id_str().replace("pdb=", "").replace('"', "").strip()

class _chemistry(object):
  """X-H lengths and H counts from the model's effective chemistry: its
  loaded monomer server (user CIFs included) and its bond restraints"""
  def __init__(self, model):
    self.mon_lib_srv = model.get_mon_lib_srv()
    geometry = model.get_restraints_manager().geometry
    self.bonds = geometry.bond_params_table
    self.fsc0 = geometry.shell_sym_tables[0].full_simple_connectivity()
    self.atoms = model.get_hierarchy().atoms()
    self._comps = {}

  def comp(self, resname):
    resname = resname.strip()
    if resname not in self._comps:
      try:
        c = self.mon_lib_srv.get_comp_comp_id_direct(resname)
      except Exception:
        c = None
      self._comps[resname] = c
    return self._comps[resname]

  def lengths(self, h, parent):
    """X-H lengths the H/D may have: the dictionary's X-ray and neutron
    values and the restraint target"""
    result = set()
    comp = self.comp(h.parent().resname)
    if comp is not None:
      names = set((parent.name.strip(), normalized_name(h).strip()))
      for b in comp.bond_list:
        if set((b.atom_id_1, b.atom_id_2)) != names: continue
        for v in (b.value_dist, b.value_dist_neutron):
          if v is not None and v > 0: result.add(v)
    bond = self.bonds.lookup(h.i_seq, parent.i_seq)
    if bond is not None: result.add(bond.distance_ideal)
    return sorted(result)

  def n_h_in_dictionary(self, parent):
    """number of H bonded to the parent's name in its dictionary, or None"""
    comp = self.comp(parent.parent().resname)
    if comp is None: return None
    hs = set(a.atom_id for a in comp.atom_list
             if a.type_symbol.strip().upper() in ("H", "D"))
    n = parent.name.strip()
    return len([b for b in comp.bond_list if
      (b.atom_id_1 == n and b.atom_id_2 in hs) or
      (b.atom_id_2 == n and b.atom_id_1 in hs)])

  def heavy_neighbours(self, atom):
    """the restrained bond partners of atom that are not H/D"""
    return [self.atoms[j] for j in self.fsc0[atom.i_seq]
            if not _is_hd(self.atoms[j])]

  def ideal_direction(self, parent):
    """unit vector of a one-H parent's ideal X-H, from its own heavy
    neighbours (opposite their sum: three neighbours, or two around an sp2
    parent), with the reason when there is none"""
    nb = self.heavy_neighbours(parent)
    if len(nb) == 2:
      if self.n_h_in_dictionary(parent) != 1:
        return None, "%s with 2 heavy neighbours carries %s H in the " \
          "dictionary (no single ideal direction)" % (_label(parent),
          self.n_h_in_dictionary(parent))
    elif len(nb) != 3:
      return None, "%s has %d heavy neighbours (no ideal direction)" % (
        _label(parent), len(nb))
    s = matrix.col((0, 0, 0))
    p = matrix.col(parent.xyz)
    for a in nb: s += (matrix.col(a.xyz) - p).normalize()
    if s.length() < 1.e-6:
      return None, "no ideal direction on %s" % _label(parent)
    return -s.normalize(), None

  def out_of_plane(self, parent, atom):
    """degrees by which atom is out of the plane of parent, its one heavy
    neighbour X and a neighbour of X in one dictionary plane with parent and
    atom; None where this is not measurable"""
    nb = self.heavy_neighbours(parent)
    comp = self.comp(parent.parent().resname)
    if len(nb) != 1 or comp is None: return None
    x = nb[0]
    name = normalized_name(atom).strip()
    planes = {}
    for p in comp.get_planes():
      planes[p.plane_id] = set(a.atom_id for a in p.plane_atoms)
    for r in self.heavy_neighbours(x):
      if r.i_seq == parent.i_seq: continue
      for members in planes.values():
        if not set((parent.name.strip(), x.name.strip(), r.name.strip(),
                    name)) <= members:
          continue
        t = matrix.dihedral_angle(sites=[r.xyz, x.xyz, parent.xyz, atom.xyz],
          deg=True)
        if t is None: return None
        return min(abs(t), 180 - abs(t)), "%s-%s-%s" % (r.name.strip(),
          x.name.strip(), parent.name.strip())
    return None

def _in_window(length, lengths):
  return any(L - _REPAIR_BOND_BELOW <= length <= L + _REPAIR_BOND_ABOVE
             for L in lengths)

def _lengths_text(lengths):
  return "/".join("%.2f" % L for L in lengths)

def repair_plan(model, candidates, match):
  """Repairs of deposited H/D layouts that the input check refuses, each
  proven from the model alone; everything else is refused with a reason.

  Parameters
  ----------
  model : mmtbx.model.manager
      The interpreted model the check ran on, with restraints; lengths and
      H counts come from its own monomer server and bond restraints (user
      CIFs included).
  candidates : group_args
      find_hd_candidates() on model.
  match : group_args
      collapse_check() of candidates (the check's own first match).

  Returns
  -------
  group_args
      split, rename, relabel, accept: lists of group_args, one per repair,
      each with the atom keys (atom_keys() of model) it changes, its atom
      label (atom) and what was done and why (reason):

      - split (a): a blank H/D whose covalent partners are the altloc copies
        of one parent, one copy per real conformer, at a valid position
        from every copy (each copy's ideal H from its own neighbours, at
        both dictionary lengths, within 0.10 A of the other copies' ideal
        positions; direction within 15 deg and length in [L - 0.20,
        L + 0.05] at either length): one copy per conformer, occupancy
        shared by the parent copies' occupancy, B/ANISOU copied.
      - rename (b): names permuted between the H set and the D set of one
        parent: each H with a mutual unique opposite isotope within 0.20 A
        takes that D's name; a leftover H whose name is then taken gets the
        one name left free when no D carries it.
      - relabel (d): H and D in heavy-atom conformer groups on a blank
        parent, a mutual unique pair within 0.20 A, H + D = parent
        occupancy (0.01), not both occupancies equal to their conformers'
        (0.01): moved into two free H/D-only altloc letters, values
        unchanged. That the two form an exchange pair is an inference from
        the occupancies, not proof of how refinement constrained them.
      - accept (e): H and D too far apart on a parent with one H site in
        the dictionary, carrying only these two, each at a valid length,
        at most REPAIR_E_LIMIT apart (an empirical repair limit). The pair
        collapses with collapse_check(accepted=...), the H's position kept.
      - conformer (c): same-name H and D on the altloc copies of one parent
        atom in two real conformers of one residue_group, each in its own
        parent's atom_group, neither in a pair, each occupancy equal to its
        conformer's (0.01): kept as two singles, no atom changed (with
        collapse_check(accepted=...)). That each conformer carries one
        isotope is an inference from the occupancies. Each entry carries
        what conformer_site_isotope() needs. The same layout with H + D =
        the copies' occupancy otherwise is consistent with an exchange pair
        over the two copies, but the occupancies do not determine how H and
        D are distributed over the conformers: refused.

      refused: {atom label: reason} for every refused atom (labels as in
      the check's lines, without "pdb=" and quotes).
  """
  hierarchy = model.get_hierarchy()
  atoms = hierarchy.atoms()
  keys = atom_keys(hierarchy)
  chem = _chemistry(model)
  by_id = dict((a.id_str(), i) for i, a in enumerate(atoms))
  plan = group_args(split=[], rename=[], relabel=[], accept=[], conformer=[],
    refused={})
  def refuse(atom_ids, reason):
    for a in atom_ids:
      if a is not None: plan.refused[_label_of_id(a)] = reason
  _plan_split(plan, candidates, atoms, keys, by_id, chem, refuse)
  done = set()
  _plan_rename(plan, match, model, atoms, keys, refuse, done)
  table = None
  singles = set(s.i_seq for s in candidates.singles)
  for x in match.problems:
    if x.kind != "pair":
      if x.kind == "different parents":
        _plan_conformer(plan, x.candidate, atoms, keys, by_id, singles,
          refuse)
      else:
        refuse([x.h_id, x.d_id] + [getattr(x.candidate, "atom_id", None),
          getattr(x.candidate, "partner_id", None)],
          "no repair for this pattern")
      continue
    c = x.candidate
    if (c.i_h, c.i_d) in done: continue
    if x.reasons == [_LAYOUT_TEXT]:
      _plan_relabel(plan, c, atoms, keys, refuse)
    elif len(x.reasons) == 1 and x.reasons[0] == _apart_text(c.distance):
      if table is None: table = _hd_on_parents(model, atoms, keys)
      _plan_accept(plan, c, atoms, keys, by_id, chem, table, refuse)
    else:
      refuse([x.h_id, x.d_id], "no repair for: %s" % "; ".join(x.reasons))
  return plan

def _label_of_id(id_str):
  return id_str.replace("pdb=", "").replace('"', "").strip()

def _hd_on_parents(model, atoms, keys):
  """{parent id_str: [i_seq]} of the H/D with one covalent partner, in the
  asymmetric unit"""
  table = h_parent_table(model)
  result = {}
  for i, a in enumerate(atoms):
    if not _is_hd(a): continue
    e = table[keys[i]]
    if e.status == "ok" and e.parent_sym_ops == ("x,y,z",):
      result.setdefault(e.parent_id, []).append(i)
  return result

def _plan_split(plan, candidates, atoms, keys, by_id, chem, refuse):
  # (a)
  for u in candidates.unresolved:
    if u.status != "ambiguous": continue
    a = atoms[u.i_seq]
    why = None
    parents = sorted(set(by_id.get(s[0]) for s in u.source))
    ps = [atoms[i] for i in parents if i is not None]
    rg = a.parent().parent()
    if a.parent().altloc != "":
      why = "H/D not in the blank altloc"
    elif None in parents or len(ps) < 2:
      why = "partners not found"
    elif any(s[3] != "x,y,z" for s in u.source):
      why = "a partner is a symmetry mate"
    elif len(set(p.name for p in ps)) != 1:
      why = "partners are different atoms"
    elif any(p.parent().parent().memory_id() != rg.memory_id() for p in ps):
      why = "a partner is in another residue"
    else:
      alts = [p.parent().altloc for p in ps]
      # real conformers carry heavy atoms; H/D-only groups are layout letters
      conformers = sorted(ag.altloc for ag in rg.atom_groups()
        if ag.altloc != "" and ag.resname == a.parent().resname and
        not _hd_only(ag))
      if "" in alts or len(set(alts)) != len(alts):
        why = "partners are not altloc copies"
      elif sorted(alts) != conformers:
        why = "not one partner copy per conformer (%s; conformers %s)" % (
          " ".join(sorted(alts)), " ".join(conformers))
      elif sum(p.occ for p in ps) <= 0:
        why = "parent copies with total occupancy 0 (no basis for sharing it)"
    per_copy, spread = [], {}
    if why is None:
      lengths = chem.lengths(a, ps[0])
      if not lengths:
        why = "no X-H length for %s" % _label(a)
      ideal_at = dict((L, []) for L in lengths)
      for p in ps:
        if why is not None: break
        d, no_direction = chem.ideal_direction(p)
        if d is None:
          why = no_direction
          break
        v = matrix.col(a.xyz) - matrix.col(p.xyz)
        angle = v.angle(d, deg=True)
        if angle is None:
          # on the parent's position: no direction
          why = "bond %.3f A to %s (dictionary %s A)" % (v.length(),
            _label(p), _lengths_text(lengths))
          break
        per_copy.append((p.parent().altloc, v.length(), angle))
        if angle > _REPAIR_ANGLE:
          why = "%.1f deg from the ideal direction on %s (limit %.0f)" % (
            angle, _label(p), _REPAIR_ANGLE)
        elif not _in_window(v.length(), lengths):
          why = "bond %.3f A to %s (dictionary %s A)" % (v.length(),
            _label(p), _lengths_text(lengths))
        for L in lengths:
          ideal_at[L].append(matrix.col(p.xyz) + d * L)
      if why is None:
        for L, pos in ideal_at.items():
          spread[L] = max((x - y).length() for x in pos for y in pos)
        worst = max(spread.values())
        if worst > _REPAIR_IDEAL_TOLERANCE:
          why = "ideal positions on the copies %.3f A apart (limit %.2f)" % (
            worst, _REPAIR_IDEAL_TOLERANCE)
    if why is None:
      n = normalized_name(a)
      for p in ps:
        for other in p.parent().atoms():
          if _is_hd(other) and normalized_name(other) == n:
            why = "name %s already in conformer %s" % (n.strip(),
              p.parent().altloc)
    if why is not None:
      refuse([u.atom_id], why)
      continue
    total = sum(p.occ for p in ps)
    plan.split.append(group_args(key=keys[u.i_seq],
      parent_keys=[keys[p.i_seq] for p in ps],
      shares=[p.occ / total for p in ps],
      atom=_label(a),
      reason="split over the copies of %s (conformers %s): ideal positions "
        "%s A apart; per copy %s" % (ps[0].name.strip(),
        " ".join(p.parent().altloc for p in ps),
        "/".join("%.3f" % spread[L] for L in sorted(spread)),
        ", ".join("%s %.3f A %.1f deg" % c for c in per_copy))))

def _plan_rename(plan, match, model, atoms, keys, refuse, done):
  # (b): name-mismatch pair problems grouped by parent
  groups = {}
  for x in match.problems:
    if x.kind != "pair" or not x.candidate.name_mismatch: continue
    c = x.candidate
    g = groups.setdefault(c.parent_id, group_args(bad=[], ids=set(),
      pairs=set()))
    g.ids.update([x.h_id, x.d_id])
    g.pairs.add((c.i_h, c.i_d))
    other = [r for r in x.reasons if not (r.startswith("names ") or
             r == _apart_text(c.distance))]
    if other or c.n_competing or c.parent_altloc != "" or not c.layout_ok:
      g.bad.append("; ".join(x.reasons))
  if not groups: return
  on_parent = _hd_on_parents(model, atoms, keys)
  for parent_id in sorted(groups):
    g = groups[parent_id]
    done.update(g.pairs)
    members = on_parent.get(parent_id, [])
    hs = [i for i in members if atoms[i].element.strip().upper() == "H"]
    ds = [i for i in members if atoms[i].element.strip().upper() == "D"]
    why = None
    if g.bad: why = "not only a name permutation: " + " | ".join(g.bad)
    renames = []
    def near(i, pool):
      return [j for j in pool if _distance(atoms[i], atoms[j]) <=
              _REPAIR_PAIR_TOLERANCE]
    if why is None:
      for i in hs:
        nd = near(i, ds)
        if len(nd) > 1:
          why = "%s: %d D within %.2f A" % (_label(atoms[i]), len(nd),
            _REPAIR_PAIR_TOLERANCE)
          break
        if not nd: continue
        j = nd[0]
        if near(j, hs) != [i]:
          why = "%s: %d H within %.2f A" % (_label(atoms[j]),
            len(near(j, hs)), _REPAIR_PAIR_TOLERANCE)
          break
        new = normalized_name(atoms[j])
        if new != atoms[i].name:
          renames.append((i, new, "the name of the D at its site, %s, "
            "%.3f A" % (atoms[j].name.strip(), _distance(atoms[i], atoms[j]))))
      if why is None and not renames:
        why = "no partner within %.2f A" % _REPAIR_PAIR_TOLERANCE
    if why is None:
      # a leftover H whose name a renamed H now takes gets the one name left
      # free, when no D carries it (that D sits elsewhere: refused)
      final = dict((i, atoms[i].name) for i in hs)
      for i, new, _ in renames: final[i] = new
      renamed = set(i for i, _, _ in renames)
      clash = [i for i in hs if i not in renamed and
               any(final[j] == final[i] for j in renamed)]
      d_names = set(normalized_name(atoms[j]) for j in ds)
      pool = set(atoms[i].name for i in hs) | d_names
      free = sorted(pool - set(final.values()) - d_names)
      if len(clash) == 1 and len(free) == 1:
        renames.append((clash[0], free[0], "the one name left free; no D "
          "carries it"))
        final[clash[0]] = free[0]
      elif clash:
        why = "%d leftover H, %d free names" % (len(clash), len(free))
    if why is None:
      seen = set()
      for i in hs:
        k = (atoms[i].parent().memory_id(), final[i])
        if k in seen:
          why = "name %s twice after renaming" % final[i].strip()
          break
        seen.add(k)
    if why is not None:
      refuse(g.ids, why)
      continue
    for i, new, reason in renames:
      plan.rename.append(group_args(key=keys[i], new=new,
        atom=_label(atoms[i]),
        reason="renamed %s (%s)" % (new.strip(), reason)))

def _plan_relabel(plan, c, atoms, keys, refuse):
  # (d)
  h, d = atoms[c.i_h], atoms[c.i_d]
  ids = [c.h_id, c.d_id]
  ha, da = h.parent().altloc, d.parent().altloc
  co = (_conformer_occupancy(h.parent()), _conformer_occupancy(d.parent()))
  why = None
  if c.n_competing or c.name_mismatch or c.parent_altloc != "":
    why = "not a plain pair"
  elif ha == "" or da == "" or ha == da:
    why = "H and D not in two altloc groups"
  elif c.distance > _REPAIR_PAIR_TOLERANCE:
    why = "%s (limit %.2f)" % (_apart_text(c.distance),
      _REPAIR_PAIR_TOLERANCE)
  elif abs(c.occ_sum - c.parent_occ) > _OCC_TOLERANCE:
    why = "H + D occupancy %.2f, parent %.2f (no exchange pair)" % (
      c.occ_sum, c.parent_occ)
  elif not all(isinstance(o, float) for o in co):
    why = "conformer occupancy not one value (%s)" % (co,)
  elif abs(h.occ - co[0]) <= _OCC_TOLERANCE and \
       abs(d.occ - co[1]) <= _OCC_TOLERANCE:
    why = ("H and D occupancies equal their conformers' (%.2f/%.2f): one "
           "isotope per conformer or an exchange pair" % co)
  if why is not None:
    refuse(ids, why)
    return
  plan.relabel.append(group_args(h_key=keys[c.i_h], d_key=keys[c.i_d],
    atom="%s / %s" % (_label(h), _label(d)),
    reason="moved into H/D-only altloc groups: H %.2f in %s (conformer "
      "%.2f), D %.2f in %s (conformer %.2f), %.3f A apart; read as an exchange "
      "pair (H + D = parent occupancy), an inference" % (h.occ, ha, co[0],
      d.occ, da, co[1], c.distance)))

def _plan_accept(plan, c, atoms, keys, by_id, chem, on_parent, refuse):
  # (e)
  h, d = atoms[c.i_h], atoms[c.i_d]
  ids = [c.h_id, c.d_id]
  i_parent = by_id.get(c.parent_id)
  why = None
  if c.n_competing or c.name_mismatch or c.parent_altloc != "" or \
     not c.layout_ok or i_parent is None:
    why = "not a plain pair"
  if why is None:
    p = atoms[i_parent]
    n = chem.n_h_in_dictionary(p)
    if n != 1:
      why = "%s carries %s H in the dictionary (one H site needed)" % (
        p.name.strip(), n)
      measured = [chem.out_of_plane(p, a) for a in (h, d)]
      if all(m is not None for m in measured):
        k = 0 if measured[0][0] >= measured[1][0] else 1
        why += "; %s %.0f deg out of the %s plane" % (
          (h, d)[k].element.strip(), measured[k][0], measured[k][1])
    elif c.distance > REPAIR_E_LIMIT:
      why = "%s (limit %.2f A, empirical)" % (_apart_text(c.distance),
        REPAIR_E_LIMIT)
  if why is None:
    others = [i for i in on_parent.get(c.parent_id, [])
              if i not in (c.i_h, c.i_d)]
    if others:
      why = "other H/D on %s: %s" % (p.name.strip(),
        ", ".join(_label(atoms[i]) for i in others))
  words = []
  if why is None:
    lengths = chem.lengths(h, p)
    for a in (h, d):
      v = matrix.col(a.xyz) - matrix.col(p.xyz)
      if not lengths or not _in_window(v.length(), lengths):
        why = "bond %.3f A to %s (dictionary %s A)" % (v.length(),
          _label(p), _lengths_text(lengths))
        break
      words.append("%s %.2f at %.3f A" % (a.element.strip(), a.occ,
        v.length()))
  if why is not None:
    refuse(ids, why)
    return
  plan.accept.append(group_args(h_key=keys[c.i_h], d_key=keys[c.i_d],
    atom="%s / %s" % (_label(h), _label(d)),
    reason="paired %.3f A apart on a one-H parent (limit %.2f A, an "
      "empirical repair limit); %s; the H position is kept" % (c.distance,
      REPAIR_E_LIMIT, ", ".join(words))))

def _plan_conformer(plan, dp, atoms, keys, by_id, singles, refuse):
  # (c)
  h, d = atoms[dp.i_h], atoms[dp.i_d]
  ids = [dp.h_id, dp.d_id]
  i_ph, i_pd = by_id.get(dp.h_parent[0]), by_id.get(dp.d_parent[0])
  why = None
  if i_ph is None or i_pd is None:
    why = "parents not found"
  elif dp.h_parent[1] != ("x,y,z",) or dp.d_parent[1] != ("x,y,z",):
    why = "a parent is a symmetry mate"
  else:
    ph, pd = atoms[i_ph], atoms[i_pd]
    gh, gd = ph.parent(), pd.parent()
    if ph.name != pd.name or gh.resname != gd.resname or \
       gh.parent().memory_id() != gd.parent().memory_id():
      why = "parents are different atoms"
    elif gh.altloc == "" or gd.altloc == "" or gh.altloc == gd.altloc:
      why = "parents are not altloc copies in two conformers"
    elif h.parent().memory_id() != gh.memory_id() or \
         d.parent().memory_id() != gd.memory_id():
      why = "H or D not in its parent's conformer group"
    elif dp.i_h not in singles or dp.i_d not in singles:
      why = "H or D also in a pair"
  if why is not None:
    refuse(ids, "no repair for H and D on different parents: " + why)
    return
  co = (_conformer_occupancy(gh), _conformer_occupancy(gd))
  if not all(isinstance(o, float) for o in co):
    refuse(ids, "conformer occupancy not one value (%s)" % (co,))
    return
  words = "%s %.2f on %s (conformer %.2f), %s %.2f on %s (conformer %.2f)" % (
    h.name.strip(), h.occ, _label(ph), co[0], d.name.strip(), d.occ,
    _label(pd), co[1])
  if abs(h.occ - co[0]) <= _OCC_TOLERANCE and \
     abs(d.occ - co[1]) <= _OCC_TOLERANCE:
    plan.conformer.append(group_args(h_key=keys[dp.i_h], d_key=keys[dp.i_d],
      atom="%s / %s" % (_label(h), _label(d)),
      reason="kept as two singles: %s; read as one isotope per conformer "
        "(each at its conformer's occupancy), an inference" % words,
      water=_is_water(gh), exchangeable=ph.element.strip().upper() != "C",
      name=normalized_name(d)))
  elif abs(h.occ + d.occ - (ph.occ + pd.occ)) <= _OCC_TOLERANCE:
    refuse(ids, "%s, H + D = %.2f: consistent with an exchange pair over the "
      "two copies, but the occupancies do not determine how H and D are "
      "distributed over the conformers (%.2f/%.2f), so the site is not "
      "repaired; remove the H and D at this site and rerun with "
      "approach=add keep_existing_H=True, which places them anew "
      "(approach=optimize adds no H)" % (words, h.occ + d.occ, co[0], co[1]))
  else:
    refuse(ids, "%s: neither one isotope per conformer nor an exchange "
      "pair" % words)

def repairs_made(plan):
  """True if apply_repairs() changes anything (case e changes no atom)"""
  return bool(plan.split or plan.rename or plan.relabel)

def apply_repairs(hierarchy, plan, key_map=None):
  """Apply repair_plan()'s split, rename and relabel to hierarchy.

  Parameters
  ----------
  hierarchy : iotbx.pdb.hierarchy.root
      The hierarchy to change in place: the planned model's or one that
      contains its atoms.
  plan : group_args
      repair_plan() result.
  key_map : dict or None
      {plan atom key: atom key in hierarchy} when the two differ (the plan
      was made on a selection of hierarchy).

  Returns
  -------
  dict
      {atom key after the repair: atom key before} for every atom whose
      identity changed (renamed, relabelled, split copies); keys in
      hierarchy's terms.
  """
  def k(key):
    return key if key_map is None else key_map[key]
  keys = atom_keys(hierarchy)
  index = dict((key, i) for i, key in enumerate(keys))
  atoms = hierarchy.atoms()
  changed = []     # (atom, key before)
  for r in plan.rename:
    a = atoms[index[k(r.key)]]
    a.set_name(r.new)
    changed.append((a, k(r.key)))
  # (d): every relabelled pair of a residue_group into the same two free
  # letters, H in the first, D in the second
  by_rg = {}
  for x in plan.relabel:
    h, d = atoms[index[k(x.h_key)]], atoms[index[k(x.d_key)]]
    by_rg.setdefault(h.parent().parent().memory_id(), []).append(
      (h, d, k(x.h_key), k(x.d_key)))
  for items in by_rg.values():
    rg = items[0][0].parent().parent()
    used = set(ag.altloc for ag in rg.atom_groups())
    free = [c for c in _ALTLOC_LETTERS if c not in used][:2]
    if len(free) < 2:
      raise Sorry("Internal error: no free altloc letters for %s" %
        items[0][0].id_str())
    made = {}
    for h, d, h_key, d_key in items:
      for a, alt, key in ((h, free[0], h_key), (d, free[1], d_key)):
        group = (alt, a.parent().resname)
        if group not in made:
          made[group] = iotbx.pdb.hierarchy.atom_group(altloc=alt,
            resname=a.parent().resname)
          rg.append_atom_group(made[group])
        moved = a.detached_copy()
        made[group].append_atom(moved)
        a.parent().remove_atom(a)
        changed.append((moved, key))
  # (a)
  splits = []
  for s in plan.split:
    a = atoms[index[k(s.key)]]
    targets = [atoms[index[k(p)]].parent() for p in s.parent_keys]
    splits.append((a, targets, s.shares, k(s.key)))
  for a, targets, shares, key in splits:
    for ag, share in zip(targets, shares):
      c = a.detached_copy()
      c.set_occ(a.occ * share)
      ag.append_atom(c)
      changed.append((c, key))
    ag = a.parent()
    rg = ag.parent()
    ag.remove_atom(a)
    # the H/D may have been the only blank atom (6TAE Gln A 13)
    _remove_empty_atom_groups(rg)
  result = {}
  for a, before in changed:
    result[atom_key(a)] = before
  _finish(hierarchy)
  return result

# ------------------------------------------------------------------------------
# Normalization and output

NEUTRON_OPTIONS = ("none", "all_h", "all_d", "hd_and_h", "hd_and_d", "all_hd")

_ALTLOC_LETTERS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"

def _values(atom):
  uij = None
  if atom.uij_is_defined(): uij = tuple(atom.uij)
  return (atom.occ, atom.b, uij)

def _set_values(atom, occ, b, uij):
  atom.set_occ(occ)
  atom.set_b(b)
  if uij is None:
    atom.uij_erase()
  else:
    atom.set_uij(uij)

def _finish(hierarchy):
  hierarchy.sort_atoms_in_place()
  hierarchy.atoms().reset_serial()
  hierarchy.atoms().reset_i_seq()

def _remove_empty_atom_groups(rg):
  for ag in list(rg.atom_groups()):
    if ag.atoms_size() == 0:
      rg.remove_atom_group(ag)

def _is_water(ag):
  return get_class(name=ag.resname) == "common_water"

def normalize_to_h(hierarchy, match):
  """H-only hierarchy from validated H/D input -> group_args(records, counts).
  match: collapse_check() result without problems or ambiguous atoms, made on
  this hierarchy (floating atoms already removed). For each pair: record
  {site_key: {"H": (occ, b, uij), "D": (occ, b, uij)}} (uij None when
  isotropic), keep the H in the blank-altloc atom_group of its residue at the
  parent occupancy, delete the D. Every other D (a single) is renamed to H in
  place, altloc kept. DOD -> HOH."""
  if len(match.problems) > 0 or \
     any(u.action == "stop" for u in match.unresolved):
    raise Sorry("Internal error: normalize_to_h() on input with problems.")
  atoms = hierarchy.atoms()
  # site keys of all pairs before any atom moves
  parts = _residue_parts(hierarchy)
  pairs = []
  for c in match.pairs:
    h, d = atoms[c.i_h], atoms[c.i_d]
    key = site_key(h, c.parent_altloc, parts[h.parent().memory_id()])
    assert key == site_key(d, c.parent_altloc, parts[d.parent().memory_id()])
    pairs.append((h, d, c, key))
  records = {}
  n_pairs = 0
  for h, d, c, key in pairs:
    if key in records:
      raise Sorry("Internal error: duplicate H/D site %s" % h.id_str())
    records[key] = {"H": _values(h), "D": _values(d)}
    h_ag = h.parent()
    rg = h_ag.parent()
    blank = [ag for ag in rg.atom_groups()
             if ag.altloc == "" and ag.resname == h_ag.resname]
    if len(blank) != 1:
      raise Sorry("Internal error: no blank-altloc atom_group for %s" %
        h.id_str())
    moved = h.detached_copy()
    moved.set_occ(c.parent_occ)
    blank[0].append_atom(moved)
    h_ag.remove_atom(h)
    d.parent().remove_atom(d)
    _remove_empty_atom_groups(rg)
    n_pairs += 1
  n_renamed = 0
  for atom in hierarchy.atoms():
    if atom.element.strip().upper() != "D": continue
    atom.set_name(normalized_name(atom))
    atom.set_element(" H")
    n_renamed += 1
  n_water = 0
  for ag in hierarchy.atom_groups():
    if ag.resname == "DOD":
      ag.resname = "HOH"
      n_water += 1
  _finish(hierarchy)
  return group_args(records=records,
    counts=group_args(n_pairs=n_pairs, n_d_renamed=n_renamed,
                      n_dod=n_water))

def conformer_site_isotope(option, entry):
  """'H' or 'D': the isotope apply_neutron_option() writes in both conformers
  of a repair_plan() conformer (c) entry, per neutron option"""
  kind = _output_kind(option, entry.water, entry.exchangeable, True)[0]
  if not entry.name.strip().startswith("H"): kind = "H"   # NOT_PAIRABLE
  return kind

def _output_kind(option, water, exchangeable, conformer):
  """'H', 'D' or 'split' per "What each mode writes"; plus a reason when an
  exchangeable or all_hd site is written as one isotope because of its
  conformer"""
  if water:
    return ("H" if option in ("none", "all_h") else "D"), None
  if option in ("none", "all_h"): return "H", None
  if option == "all_d": return "D", None
  if exchangeable:
    if conformer: return "D", "conformer"
    return "split", None
  if option == "hd_and_h": return "H", None
  if option == "hd_and_d": return "D", None
  # all_hd, non-exchangeable
  if conformer: return "H", "conformer"
  return "split", None

_MERGE_TOLERANCE = 0.05

def _merge_conformer_copies(hierarchy, parent_table, unchanged_keys=()):
  """Placement copies every H of a residue with alternates into each
  conformer group, also H whose parent is in the blank group. In a
  residue_group with at least two conformers, copies with exactly one per
  conformer, the same name and all within 0.05 A of the first are
  merged into one blank-group H (parent occupancy, B/uij of the first copy,
  the lowest altloc). Returns {memory_id of merged H: (atom_key, entry) of
  the first copy} and the merged sites."""
  keys = atom_keys(hierarchy)
  key_of = dict((a.memory_id(), k) for a, k in zip(hierarchy.atoms(), keys))
  merged = {}
  sites = []
  site_atoms = []
  for rg in hierarchy.residue_groups():
    conformers = sorted(set(ag.altloc for ag in rg.atom_groups()) - set([""]))
    # a lone alternate H (one conformer) is not a copy
    if len(conformers) < 2: continue
    copies = {}
    for ag in rg.atom_groups():
      if ag.altloc == "": continue
      for atom in ag.atoms():
        if atom.element.strip().upper() != "H": continue
        if key_of[atom.memory_id()] in unchanged_keys: continue
        e = parent_table.get(key_of[atom.memory_id()])
        if e is None or e.status != "ok" or e.parent_altloc != "": continue
        copies.setdefault((ag.resname, atom.name), []).append(atom)
    for (resname, name), atoms in sorted(copies.items()):
      if sorted(a.parent().altloc for a in atoms) != conformers: continue
      atoms.sort(key=lambda a: a.parent().altloc)
      first = atoms[0]
      if any(_distance(first, a) > _MERGE_TOLERANCE for a in atoms[1:]):
        continue
      blank = [ag for ag in rg.atom_groups()
               if ag.altloc == "" and ag.resname == resname]
      if len(blank) != 1: continue
      first_key = key_of[first.memory_id()]
      e = parent_table[first_key]
      h = first.detached_copy()
      h.set_occ(e.parent_occ)
      blank[0].append_atom(h)
      merged[h.memory_id()] = (first_key, e)
      sites.append(group_args(site=h.id_str(),
        reason="merged conformer copies",
        copies=[a.id_str() for a in atoms]))
      site_atoms.append(h)
      for a in atoms:
        a.parent().remove_atom(a)
    _remove_empty_atom_groups(rg)
  return merged, sites, site_atoms

NOT_PAIRABLE = "name not pairable by phenix.refine, kept as H"

def _conformer_occupancy(ag):
  """occupancy of the heavy atoms with ag's altloc in its residue_group if
  they agree on one value (0.01), None if they do not, "no heavy atoms" if
  there are none"""
  occ = [a.occ for g in ag.parent().atom_groups() if g.altloc == ag.altloc
         for a in g.atoms() if not _is_hd(a)]
  if not occ:
    return "no heavy atoms"
  if max(occ) - min(occ) > 0.01:
    return None
  return occ[0]

def output_label(atom):
  """atom label as in the reduce2 reports, e.g. ' DG ASER A  37 '"""
  return atom.id_str().replace("pdb=", "").replace('"', "")

def _to_d(atom):
  atom.set_name(atom.name.replace("H", "D", 1))
  atom.set_element(" D")

def apply_neutron_option(hierarchy, option, parent_table, records=None,
                         log=None, unchanged_keys=None):
  """H-only hierarchy -> H/D output for option (see NEUTRON_OPTIONS).
  Exchangeability, parent occupancy and altloc come from parent_table (keyed
  by atom_key). Split sites of a residue_group go into the first two unused
  altloc letters (H group, D group). Coincident conformer copies of an H on
  a blank-group parent are first merged into one blank-group H. Not split,
  one isotope, listed with a reason: conformer sites and H in an altloc
  group ("conformer"); residues with fewer than two free letters ("no free
  altloc letters"); names not starting with H, never split or D in any
  mode (NOT_PAIRABLE). option "none" returns at once, hierarchy unchanged.
  Values follow the restoration policy: recorded pair as pair
  -> both restored; recorded pair as one isotope -> parent occupancy, that
  isotope's B/uij; new pair -> 0.5 x parent occupancy each, B/uij of the H;
  single without record -> unchanged (occupancy note if not the parent's).
  Water names follow the content: D -> DOD, H -> HOH where it was DOD.
  Returns group_args(n_exchangeable (exchangeable output sites after the
  merge, not water), n_split, n_split_new (split without a record, at
  0.5 x parent occupancy), n_to_d, n_water (residue_groups renamed to
  DOD), n_water_to_hoh (DOD renamed to HOH), n_restored,
  partial_parent_sites, conformer_sites, unclassified_sites, merged_sites,
  notes, provenance={output atom_key: source atom_key}, unused_records
  (sorted record keys whose site is not in the hierarchy, e.g. an H deleted
  by the Optimizer)); merged H map to the first copy's atom_key. Every listed entry has site (input label),
  reason and output (label of the atom as written; the H of a split site).
  unchanged_keys: atom keys written as they are (H or D, values kept), not
  required in parent_table, mapped to themselves in the provenance."""
  if option not in NEUTRON_OPTIONS:
    raise Sorry("Internal error: unknown neutron option %s" % option)
  result = group_args(n_exchangeable=0, n_split=0, n_split_new=0, n_to_d=0,
    n_water=0,
    n_water_to_hoh=0, n_restored=0,
    partial_parent_sites=[], conformer_sites=[], unclassified_sites=[],
    merged_sites=[], notes=[], provenance={}, unused_records=[])
  if option == "none": return result
  if records is None: records = {}
  if unchanged_keys is None: unchanged_keys = set()
  used_records = set()
  merged, result.merged_sites, merged_atoms = _merge_conformer_copies(
    hierarchy, parent_table, unchanged_keys)
  # listed entries and the atom each describes; their output labels are
  # filled in once the atoms are written
  listed = list(zip(result.merged_sites, merged_atoms))
  written_as = {}   # memory_id of an input H -> the atom written for it
  # identity after the merge, one pass each
  atoms = hierarchy.atoms()
  keys = atom_keys(hierarchy)
  parts = _residue_parts(hierarchy)
  # decide per H; group split sites by residue_group
  written = []   # (atom, source key)
  plan = []
  split_by_rg = {}
  for i, atom in enumerate(atoms):
    element = atom.element.strip().upper()
    if keys[i] in unchanged_keys:
      if element in ("H", "D"): written.append((atom, keys[i]))
      continue
    if element == "D":
      raise Sorry("Internal error: apply_neutron_option() input has D: %s" %
        atom.id_str())
    if element != "H": continue
    source_key = keys[i]
    if atom.memory_id() in merged:
      source_key, e = merged[atom.memory_id()]
    elif source_key not in parent_table:
      raise Sorry("Internal error: %s not in the H parent table" %
        atom.id_str())
    else:
      e = parent_table[source_key]
    if e.status != "ok":
      raise Sorry("Internal error: %s has parent status %s" %
        (atom.id_str(), e.status))
    ag = atom.parent()
    water = _is_water(ag)
    conformer = e.parent_altloc != "" or ag.altloc != ""
    kind, reason = _output_kind(option, water, e.exchangeable, conformer)
    # an H whose name does not start with H is never split and never D:
    # phenix.refine cannot pair it, and pdb_interpretation has no D alias for
    # it (cif_types.hydrogen_deuterium_aliases(): H names starting with H)
    if kind in ("split", "D") and not atom.name.strip().startswith("H"):
      kind, reason = "H", NOT_PAIRABLE
    # one per output site, after the merge (a split pair counts once)
    if e.exchangeable and not water:
      result.n_exchangeable += 1
    record_key = site_key(atom, e.parent_altloc, parts[ag.memory_id()])
    if record_key in records: used_records.add(record_key)
    item = group_args(atom=atom, source_key=source_key, entry=e, kind=kind,
      reason=reason, water=water, record=records.get(record_key))
    plan.append(item)
    if kind == "split":
      split_by_rg.setdefault(ag.parent().memory_id(), []).append(item)
  # altloc letters, once per residue_group
  groups = {}
  for rg_id, items in split_by_rg.items():
    rg = items[0].atom.parent().parent()
    used = set(ag.altloc for ag in rg.atom_groups())
    free = [c for c in _ALTLOC_LETTERS if c not in used]
    if len(free) < 2:
      for item in items:
        item.kind = _output_kind(option, item.water, item.entry.exchangeable,
                                 True)[0]
        item.reason = "no free altloc letters"
      continue
    groups[rg_id] = (rg, free[0], free[1], {})
  def new_group(rg_id, altloc, resname):
    rg, h_alt, d_alt, made = groups[rg_id]
    key = (altloc, resname)
    if key not in made:
      ag = iotbx.pdb.hierarchy.atom_group(altloc=altloc, resname=resname)
      rg.append_atom_group(ag)
      made[key] = ag
    return made[key]
  for item in plan:
    atom, e, record = item.atom, item.entry, item.record
    site = atom.id_str()
    if item.reason is not None:
      target = result.unclassified_sites if \
        item.reason == NOT_PAIRABLE else \
        result.conformer_sites
      target.append(group_args(site=site, reason=item.reason,
        written=item.kind))
      listed.append((target[-1], atom))
    if item.kind == "split":
      rg_id = atom.parent().parent().memory_id()
      rg, h_alt, d_alt, made = groups[rg_id]
      resname = atom.parent().resname
      h = atom.detached_copy()
      d = atom.detached_copy()
      _to_d(d)
      if record is not None:
        _set_values(h, *record["H"])
        _set_values(d, *record["D"])
        result.n_restored += 1
      else:
        occ = 0.5 * e.parent_occ
        _set_values(h, occ, atom.b, _values(atom)[2])
        _set_values(d, occ, atom.b, _values(atom)[2])
      atom.parent().remove_atom(atom)
      new_group(rg_id, h_alt, resname).append_atom(h)
      new_group(rg_id, d_alt, resname).append_atom(d)
      written.append((h, item.source_key))
      written.append((d, item.source_key))
      written_as[atom.memory_id()] = h
      result.n_split += 1
      if record is None: result.n_split_new += 1
      if abs(e.parent_occ - 1) > 0.01:
        result.partial_parent_sites.append(group_args(site=site,
          reason="partial parent occupancy", parent_occ=e.parent_occ))
        listed.append((result.partial_parent_sites[-1], atom))
      continue
    if item.kind == "D":
      _to_d(atom)
      result.n_to_d += 1
    if record is not None:
      occ, b, uij = record[item.kind]
      _set_values(atom, e.parent_occ, b, uij)
      result.n_restored += 1
    else:
      # a conformer copy (H in an altloc group, parent in the blank group) is
      # compared with its conformer's occupancy, other H with the parent's
      ref, what = e.parent_occ, "parent"
      value = atom.occ
      if atom.parent().altloc != "" and e.parent_altloc == "":
        occ = _conformer_occupancy(atom.parent())
        # no heavy atom in that altloc: a lone alternate H, keeps the parent;
        # copies of one H over such groups share the parent's population
        if occ != "no heavy atoms":
          ref, what = occ, "conformer"
        else:
          name = normalized_name(atom)
          value = sum(x.occ for g in atom.parent().parent().atom_groups()
            if g.altloc != "" for x in g.atoms()
            if normalized_name(x) == name)
      if ref is not None and abs(value - ref) > 0.01:
        result.notes.append(group_args(site=atom.id_str(),
          reason="occupancy differs from %s" % what,
          text="occupancy %.2f, %s %.2f" % (value, what, ref)))
        listed.append((result.notes[-1], atom))
    written.append((atom, item.source_key))
    written_as[atom.memory_id()] = atom
  # water names follow the content, per residue_group and for every water
  # atom_group of it (including a blank-altloc O, so that the residue stays
  # one residue): D -> DOD; H -> HOH where it was DOD (other names such as
  # WAT are kept); no H/D -> unchanged
  for rg in hierarchy.residue_groups():
    water = [ag for ag in rg.atom_groups() if _is_water(ag)]
    elements = set(a.element.strip().upper() for ag in water
                   for a in ag.atoms())
    if "D" in elements and "H" in elements:
      raise Sorry("Internal error: water with H and D: %s" %
        water[0].atoms()[0].id_str())
    if "D" in elements:
      for ag in water:
        ag.resname = "DOD"
      result.n_water += 1
    elif "H" in elements and any(ag.resname == "DOD" for ag in water):
      for ag in water:
        if ag.resname == "DOD": ag.resname = "HOH"
      result.n_water_to_hoh += 1
  for rg_id in groups:
    _remove_empty_atom_groups(groups[rg_id][0])
  # output label of each listed site: the atom as written (for a split site
  # its H); site keeps the input label
  for entry, atom in listed:
    entry.output = output_label(written_as.get(atom.memory_id(), atom))
  # provenance, keyed by the written atoms' final identity
  parts = _residue_parts(hierarchy)
  for atom, source_key in written:
    key = atom_key(atom, parts[atom.parent().memory_id()])
    if key in result.provenance:
      raise Sorry("Internal error: duplicate output atom %s" % atom.id_str())
    result.provenance[key] = source_key
  result.unused_records = sorted(set(records) - used_records)
  _finish(hierarchy)
  if log is not None:
    print("Neutron option %s: %d split, %d to D, %d water, %d restored" % (
      option, result.n_split, result.n_to_d, result.n_water,
      result.n_restored), file=log)
  return result
