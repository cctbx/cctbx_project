"""Exchangeable H/D: atom and site identity, H parent table.

Pure hierarchy functions except h_parent_table(), which needs a model manager
with restraints. Used by reduce2 (neutron_option); not wired in yet.
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

def _residue_part(ag):
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

def _residue_parts(hierarchy):
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

def atom_key(atom, residue_part=None):
  """residue part + own altloc + name as written + element; residue_part
  (from _residue_parts()) saves the scan of preceding residues"""
  ag = atom.parent()
  if residue_part is None: residue_part = _residue_part(ag)
  return residue_part + (ag.altloc, atom.name, atom.element.strip().upper())

def site_key(atom, parent_altloc, residue_part=None):
  """residue part + parent's altloc + normalized name; for H or D of a
  validated pair (both give the same key) or a single"""
  if residue_part is None: residue_part = _residue_part(atom.parent())
  return residue_part + (parent_altloc, normalized_name(atom))

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

# ------------------------------------------------------------------------------
# Candidates and collapse check

_NAME_MISMATCH_MARGIN = 0.1
# distinct H sites on one parent are >= ~1.34 A apart
_SAME_SITE_RADIUS = 0.5

def _distance(a, b):
  return (matrix.col(a.xyz) - matrix.col(b.xyz)).length()

def _hd_only(ag):
  return all(_is_hd(a) for a in ag.atoms())

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
  geometric match are flagged name_mismatch."""
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
  def candidate(ih, id_, name_mismatch, geometric_partner, n_competing):
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
      geometric_partner = geometric_partner,
      n_competing       = n_competing,
      layout_ok         = (e.parent_altloc == "" and
                           h.parent().memory_id() != d.parent().memory_id() and
                           h.parent().altloc != "" and _hd_only(h.parent()) and
                           d.parent().altloc != "" and _hd_only(d.parent())),
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
          geometric = atoms[near_d].id_str()
        elif near_h != ih and d_near_h <= _SAME_SITE_RADIUS and \
           d_near_h < d_name - _NAME_MISMATCH_MARGIN:
          geometric = atoms[near_h].id_str()
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
      pairs.append(candidate(ih, id_, True, atoms[id_].id_str(), 0))
      lh.remove(ih)
      ld.remove(id_)
    for i in lh + ld:
      atom = atoms[i]
      collision = False
      if atom.element.strip().upper() == "D":
        nn = normalized_name(atom)
        collision = any(a.name == nn for a in atom.parent().atoms()
                        if a.memory_id() != atom.memory_id())
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
        for other in other_ag.atoms():
          j = index[other.memory_id()]
          if (j != i and j in info and info[j][1].status == "ok" and
              other.element.strip().upper() != element and
              normalized_name(other) == normalized_name(atom)):
            partner = other
            break
        if partner is not None: break
    if partner is not None:
      same_conformer.append(group_args(atom_key=key, atom_id=atom.id_str(),
        partner_id=partner.id_str()))
    else:
      unresolved.append(group_args(atom_key=key, atom_id=atom.id_str(),
        status=e.status, source=e.source))
  return group_args(pairs=pairs, singles=singles, unresolved=unresolved,
    different_parents=different_parents, same_conformer=same_conformer)

def collapse_check(candidates, xyz_tolerance=0.2):
  """reduce2's policy -> group_args(pairs, problems, unresolved, notes).
  A pair collapses only if it is within xyz_tolerance, has no competing
  partner and no name mismatch, its parent is in the blank altloc and H and D
  are in different H/D-only altloc groups. Problems (stop): any other pair,
  same-name H/D on different parents, a single D whose H name exists in its
  atom_group, an H/D pair in one conformer (one of them floating).
  Unresolved: floating (remove) and ambiguous (stop, restraints named). Notes (never a stop): occupancy sum vs parent occupancy, B."""
  pairs = []
  problems = []
  notes = []
  for c in candidates.pairs:
    reasons = []
    if c.distance > xyz_tolerance:
      reasons.append("H and D %.2f A apart" % c.distance)
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
      reasons.append("H and D not in H/D-only altloc groups")
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

def _merge_conformer_copies(hierarchy, parent_table):
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
  for rg in hierarchy.residue_groups():
    conformers = sorted(set(ag.altloc for ag in rg.atom_groups()) - set([""]))
    # a lone alternate H (one conformer) is not a copy
    if len(conformers) < 2: continue
    copies = {}
    for ag in rg.atom_groups():
      if ag.altloc == "": continue
      for atom in ag.atoms():
        if atom.element.strip().upper() != "H": continue
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
      for a in atoms:
        a.parent().remove_atom(a)
    _remove_empty_atom_groups(rg)
  return merged, sites

def _to_d(atom):
  atom.set_name(atom.name.replace("H", "D", 1))
  atom.set_element(" D")

def apply_neutron_option(hierarchy, option, parent_table, records=None,
                         log=None):
  """H-only hierarchy -> H/D output for option (see NEUTRON_OPTIONS).
  Exchangeability, parent occupancy and altloc come from parent_table (keyed
  by atom_key). Split sites of a residue_group go into the first two unused
  altloc letters (H group, D group). Coincident conformer copies of an H on
  a blank-group parent are first merged into one blank-group H. Not split,
  one isotope, listed with a reason: conformer sites and H in an altloc
  group ("conformer"); residues with fewer than two free letters ("no free
  altloc letters"); names not starting with H ("name not pairable by
  phenix.refine"). option "none" returns at once, hierarchy unchanged. Values follow the restoration policy: recorded pair as pair
  -> both restored; recorded pair as one isotope -> parent occupancy, that
  isotope's B/uij; new pair -> 0.5 x parent occupancy each, B/uij of the H;
  single without record -> unchanged (occupancy note if not the parent's).
  Water names follow the content: D -> DOD, H -> HOH where it was DOD.
  Returns group_args(n_split, n_to_d, n_water (residue_groups renamed to
  DOD), n_water_to_hoh (DOD renamed to HOH), n_restored,
  partial_parent_sites, conformer_sites, unclassified_sites, merged_sites,
  notes, provenance={output atom_key: source atom_key}); merged H map to
  the first copy's atom_key."""
  if option not in NEUTRON_OPTIONS:
    raise Sorry("Internal error: unknown neutron option %s" % option)
  result = group_args(n_split=0, n_to_d=0, n_water=0, n_water_to_hoh=0,
    n_restored=0,
    partial_parent_sites=[], conformer_sites=[], unclassified_sites=[],
    merged_sites=[], notes=[], provenance={})
  if option == "none": return result
  if records is None: records = {}
  merged, result.merged_sites = _merge_conformer_copies(hierarchy,
    parent_table)
  # identity after the merge, one pass each
  atoms = hierarchy.atoms()
  keys = atom_keys(hierarchy)
  parts = _residue_parts(hierarchy)
  # decide per H; group split sites by residue_group
  plan = []
  split_by_rg = {}
  for i, atom in enumerate(atoms):
    element = atom.element.strip().upper()
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
    # sites that cannot be split are written as conformer sites are
    if kind == "split" and not atom.name.strip().startswith("H"):
      kind = _output_kind(option, water, e.exchangeable, True)[0]
      reason = "name not pairable by phenix.refine"
    item = group_args(atom=atom, source_key=source_key, entry=e, kind=kind,
      reason=reason, water=water,
      record=records.get(site_key(atom, e.parent_altloc,
                                  parts[ag.memory_id()])))
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
  written = []   # (atom, source key)
  for item in plan:
    atom, e, record = item.atom, item.entry, item.record
    site = atom.id_str()
    if item.reason is not None:
      target = result.unclassified_sites if \
        item.reason == "name not pairable by phenix.refine" else \
        result.conformer_sites
      target.append(group_args(site=site, reason=item.reason,
        written=item.kind))
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
      result.n_split += 1
      if abs(e.parent_occ - 1) > 0.01:
        result.partial_parent_sites.append(group_args(site=site,
          reason="partial parent occupancy", parent_occ=e.parent_occ))
      continue
    if item.kind == "D":
      _to_d(atom)
      result.n_to_d += 1
    if record is not None:
      occ, b, uij = record[item.kind]
      _set_values(atom, e.parent_occ, b, uij)
      result.n_restored += 1
    elif abs(atom.occ - e.parent_occ) > 0.01:
      result.notes.append(group_args(site=atom.id_str(),
        reason="occupancy differs from parent",
        text="occupancy %.2f, parent %.2f" % (atom.occ, e.parent_occ)))
    written.append((atom, item.source_key))
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
  # provenance, keyed by the written atoms' final identity
  parts = _residue_parts(hierarchy)
  for atom, source_key in written:
    key = atom_key(atom, parts[atom.parent().memory_id()])
    if key in result.provenance:
      raise Sorry("Internal error: duplicate output atom %s" % atom.id_str())
    result.provenance[key] = source_key
  _finish(hierarchy)
  if log is not None:
    print("Neutron option %s: %d split, %d to D, %d water, %d restored" % (
      option, result.n_split, result.n_to_d, result.n_water,
      result.n_restored), file=log)
  return result
