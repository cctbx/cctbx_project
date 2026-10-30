from __future__ import absolute_import, division, print_function
import os, shutil, sys, tempfile
from libtbx import group_args
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex
from six.moves import cStringIO as StringIO
from iotbx.map_model_manager import map_model_manager

# Supplied TLS in local anisotropic sharpening (model_sharpen with
#  tlso_group_info and replace_aniso_with_tls_equiv) and averaging of
#  scaling information over locations (_average_scale_factor_info_over_xyz).
# Each case can be run alone:
#  libtbx.python tst_map_model_manager_supplied_tls.py options_passed_on

# ---------------------------------------------------------------------------
# model_sharpen with two supplied TLS groups (residue ranges split at
#  SPLIT_A or SPLIT_B residues).  Runs are cached and shared by the cases.
#  The generated model has residues 219-228 (86 atoms); with d_min 3 and
#  n_boxes 8 there is one original grid entry, at about (10.8, 17.7, 11.4).
#   Split A (219:222 | 223:228; 36 and 50 atoms): the grid entry is inside
#    both masks and nearest to an atom of 223:228 (overlap_nearest); 50
#    atoms give a different number of model points under the rule using
#    the original entry count and under the old rule using the growing
#    entry count when 219:222 is processed first (sampling_count).
#   Split B (219:220 | 221:228): the grid entry is not inside the mask of
#    219:220 (empty_mask_group).

def get_mmm(d_min = 3):
  mmm = map_model_manager()
  mmm.generate_map(d_min = d_min, wrapping = False)
  mmm.set_log(StringIO())
  mmm.set_resolution(d_min)
  return mmm

def run_in_temp_dir(f):
  '''Run f() in a new temporary directory (library code writes files such
     as new_mask.ccp4 to the current directory), then remove it'''
  cwd = os.getcwd()
  tmp_dir = tempfile.mkdtemp()
  try:
    os.chdir(tmp_dir)
    return f()
  finally:
    os.chdir(cwd)
    shutil.rmtree(tmp_dir, ignore_errors = True)

T1 = (0.30, 0.10, 0.15, 0.02, 0.01, 0.0)
L1 = (0.002, 0.001, 0.0015, 0.0, 0.0, 0.0)
T2 = (0.05, 0.40, 0.08, 0.0, -0.02, 0.01)
L2 = (0.001, 0.003, 0.001, 0.0, 0.0, 0.0)
U_EPS = 1.e-6       # agreement of U values (A**2) with the expected ones
U_DIFFERENT = 0.01  # T1 and T2 differ by up to 0.30 A**2 on the diagonal
SPLIT_A = 4         # first group: first 4 residues
SPLIT_B = 2         # first group: first 2 residues
MIN_MODEL_POINTS = 25 # minimum number_to_add in the code under test

def get_resseqs(mmm):
  ph = mmm.model().get_hierarchy()
  return sorted(set([rg.resseq_as_int() for rg in ph.residue_groups()]))

def get_groups(mmm, split):
  '''Two TLS groups on adjacent, disjoint residue ranges (the first split
     residues and the rest), each with its own T, L and origin at the group
     centre'''
  resseqs = get_resseqs(mmm)
  mid = resseqs[split]
  return get_groups_from_ranges(mmm,
    ((resseqs[0], mid - 1), (mid, resseqs[-1])))

def get_groups_from_ranges(mmm, ranges):
  '''Two TLS groups on the two residue ranges (first, last), with T1, L1
     and T2, L2 and origin at the group centre'''
  from mmtbx_tls_ext import tlso
  groups = []
  for (first, last), T, L in zip(ranges, (T1, T2), (L1, L2)):
    sel = "resseq %s:%s" %(first, last)
    center = mmm.model().apply_selection_string(sel).get_sites_cart().mean()
    groups.append(group_args(selection = sel,
      tlso = tlso(t = T, l = L, s = (0.,) * 9, origin = tuple(center))))
  return groups

_tls_runs = {}

def run_two_groups(split = SPLIT_A, reverse = False):
  '''model_sharpen (n_boxes 8) with the two groups of get_groups(split);
     see run_tls_groups'''
  return run_tls_groups(('split', split), lambda: model_sharpen_setup(
    get_groups_f = lambda mmm: get_groups(mmm, split)), reverse = reverse)

def model_sharpen_setup(get_groups_f, n_boxes = 8):
  '''(mmm, groups, sharpen) for model_sharpen of get_mmm() with the
     groups returned by get_groups_f(mmm)'''
  mmm = get_mmm()
  def sharpen(tlso_group_info):
    mmm.model_sharpen(local_sharpen = True, anisotropic_sharpen = True,
      tlso_group_info = tlso_group_info,
      replace_aniso_with_tls_equiv = True,
      n_boxes = n_boxes, nproc = 1)
  return mmm, get_groups_f(mmm), sharpen

def run_tls_groups(key, setup, reverse = False):
  '''Run the sharpening given by setup() (returns mmm, the two TLS groups
     and a function sharpen(tlso_group_info)) with the two groups (in
     reversed order if reverse) and record what the TLS step did.  Class
     methods are wrapped and always restored.  Recorded:
      - options seen by _run_group_of_anisotropic_sharpen
      - for each entry passed to _update_scale_factor_info_from_aniso,
        whether it has overall_scale and whether its target_scale_factors
        were changed by the update
      - in _analyze_aniso_replace_with_supplied: its inputs, the number of
        original entries, the mask returned by
        _create_mask_from_selection_as_string for each selection, the
        model sites of each selection, and the xyz and
        uu_b_cart_as_u_cart of every entry when it returns (all of these
        from the last call), and for every call the selections, T and L
        lists it received (tl_calls).
     Nothing here depends on the order in which that function makes masks
     or processes groups; group membership is worked out afterwards from
     these records (see get_ownership).'''
  key = (key, bool(reverse))
  if key in _tls_runs:
    return _tls_runs[key]
  M = map_model_manager
  rec = group_args(run_group_kw_list = [], update_calls = [],
     replace_calls = 0, active = False, kw = None, masks = {},
     n_original = None, entries = None, group_definitions = None,
     selections = None, sites = None, has_model = None, tl_calls = [])
  originals = dict(
    _update_scale_factor_info_from_aniso =
      M._update_scale_factor_info_from_aniso,
    _run_group_of_anisotropic_sharpen = M._run_group_of_anisotropic_sharpen,
    _analyze_aniso_replace_with_supplied =
      M._analyze_aniso_replace_with_supplied,
    _create_mask_from_selection_as_string =
      M._create_mask_from_selection_as_string)

  def wrapped_update(self, scale_factor_info, **kw):
    before = []
    for sgi in scale_factor_info.value_list:
      before.append((bool(sgi.overall_scale), [list(si.target_scale_factors)
        for si in (sgi.scaling_info_list or [])]))
    result = originals['_update_scale_factor_info_from_aniso'](
       self, scale_factor_info, **kw)
    entries = []
    for (has_overall_scale, b), sgi in zip(before,
        scale_factor_info.value_list):
      after = [list(si.target_scale_factors)
        for si in (sgi.scaling_info_list or [])]
      entries.append(group_args(has_overall_scale = has_overall_scale,
        changed = (after != b)))
    rec.update_calls.append(group_args(kw = dict(kw), entries = entries))
    return result

  def wrapped_run(self, **kw):
    rec.run_group_kw_list.append(dict(
      get_scale_as_aniso_u = kw.get('get_scale_as_aniso_u'),
      get_tls_from_u = kw.get('get_tls_from_u'),
      replace_aniso_with_tls_equiv = kw.get('replace_aniso_with_tls_equiv'),
      tlso_group_info = bool(kw.get('tlso_group_info'))))
    return originals['_run_group_of_anisotropic_sharpen'](self, **kw)

  def wrapped_mask(self, *args, **kw):
    result = originals['_create_mask_from_selection_as_string'](
       self, *args, **kw)
    if rec.active:
      rec.masks[kw.get('selection_string')] = result
    return result

  def wrapped_replace(self, scale_factor_info, **kw):
    rec.replace_calls += 1
    rec.active = True
    rec.kw = dict(kw)
    tgi = kw['tlso_group_info']
    rec.tl_calls.append(group_args(
      selections = list(getattr(tgi, 'tlso_selection_list', None) or []),
      T_list = [tuple(x) for x in (getattr(tgi, 'T_list', None) or [])],
      L_list = [tuple(x) for x in (getattr(tgi, 'L_list', None) or [])]))
    rec.n_original = len(scale_factor_info.value_list)
    try:
      result = originals['_analyze_aniso_replace_with_supplied'](
         self, scale_factor_info, **kw)
    finally:
      rec.active = False
    rec.selections = list(kw['tlso_group_info'].tlso_selection_list)
    model = self.get_model_by_id(kw.get('model_id') or 'model')
    rec.has_model = model is not None
    rec.sites = dict([(s, model.apply_selection_string(s).get_sites_cart()
       if model is not None else None) for s in rec.selections])
    rec.entries = [group_args(index = i, xyz = tuple(xyz),
       uu = v.get('uu_b_cart_as_u_cart'))
       for i, (v, xyz) in enumerate(zip(scale_factor_info.value_list,
         scale_factor_info.xyz_list))]
    return result

  for name, f in (('_update_scale_factor_info_from_aniso', wrapped_update),
      ('_run_group_of_anisotropic_sharpen', wrapped_run),
      ('_analyze_aniso_replace_with_supplied', wrapped_replace),
      ('_create_mask_from_selection_as_string', wrapped_mask)):
    setattr(M, name, f)

  def sharpen():
    mmm, groups, sharpen_f = setup()
    rec.group_definitions = dict([(g.selection, g) for g in groups])
    if reverse:
      groups = list(reversed(groups))
    tlso_group_info = group_args(
      tlso_list = [g.tlso for g in groups],
      selection_as_text_list = [g.selection for g in groups])
    sharpen_f(tlso_group_info)
  try:
    run_in_temp_dir(sharpen)
  finally:
    for name, f in originals.items():
      setattr(M, name, f)
  _tls_runs[key] = rec
  return rec

def tst_options_passed_on():
  '''With supplied TLS and replace_aniso_with_tls_equiv (get_scale_as_aniso_u
     and get_tls_from_u not given), the sharpening driver receives
     get_scale_as_aniso_u = True and get_tls_from_u = True and updates
     scale factors from the aniso U values'''
  record = run_two_groups()
  print("run_group options:", record.run_group_kw_list)
  print("update calls:", len(record.update_calls))
  assert record.run_group_kw_list, "no call to _run_group_of_anisotropic_sharpen"
  for kw in record.run_group_kw_list:
    assert kw['get_scale_as_aniso_u'] is True, kw
    assert kw['get_tls_from_u'] is True, kw
  assert len(record.update_calls) > 0, \
    "_update_scale_factor_info_from_aniso not called"

def tst_all_tls_entries_updated():
  '''Every entry passed to _update_scale_factor_info_from_aniso has
     overall_scale and has its target_scale_factors changed'''
  record = run_two_groups()
  assert len(record.update_calls) > 0, \
    "_update_scale_factor_info_from_aniso not called"
  for call in record.update_calls:
    n = len(call.entries)
    n_with = len([e for e in call.entries if e.has_overall_scale])
    n_changed = len([e for e in call.entries if e.changed])
    print("update call: entries %s with overall_scale %s changed %s" %(
      n, n_with, n_changed))
    assert n > 0
    assert n_with == n, (n_with, n)
    assert n_changed == n, (n_changed, n)

def inside_mask(mask_map_manager, xyz):
  from iotbx.map_model_manager import is_inside_mask
  site_frac = mask_map_manager.crystal_symmetry().unit_cell().fractionalize(
    xyz)
  return is_inside_mask(mask_map_manager, site_frac = site_frac,
     inside = True)

def nearest_distance(sites, xyz):
  return (sites - xyz).norms().min_max_mean().min

def get_ownership(rec, model_points = True):
  '''Group owning each entry written by the TLS step, worked out in the
     test from the records of the run (model-point entries are left out
     unless model_points):
       model-point entry (index >= number of original entries): the group
         whose selection created it.  Model points are placed exactly at
         atom sites of their selection and the selections are disjoint, so
         the group is the one whose selected sites include the entry xyz.
       original (grid) entry inside one or more groups' masks: the group
         with the nearest selected atom; on an exact tie the group listed
         first.  Without a model (no selected atoms) the distance is taken
         as 0 for every group, so the group listed first.  Entries inside
         no mask are not owned (not written).
     Returns {index: group_args(selection, kind, in_masks)}'''
  owner = {}
  for e in rec.entries:
    if e.index >= rec.n_original:
      if not model_points:
        continue
      found = [s for s in rec.selections if
        nearest_distance(rec.sites[s], e.xyz) == 0.]
      assert len(found) == 1, (e.xyz, found)
      owner[e.index] = group_args(selection = found[0], kind = 'model',
        in_masks = None)
    else:
      in_masks = [s for s in rec.selections if inside_mask(rec.masks[s],
         e.xyz)]
      if not in_masks:
        continue
      best = None
      for s in in_masks: # list order, so an exact tie keeps the first
        d = nearest_distance(rec.sites[s], e.xyz) if rec.has_model else 0.
        if best is None or d < best[1]:
          best = (s, d)
      owner[e.index] = group_args(selection = best[0], kind = 'grid',
        in_masks = in_masks)
  return owner

def check_tls_run(rec, model_points = True, n_calls = 1):
  '''Common checks of a two-group run (n_calls calls of the TLS step; the
     records are those of the last call); returns the overall U (removed
     from each TLS U) and the ownership of entries (get_ownership)'''
  from cctbx import adptbx
  assert rec.replace_calls == n_calls, (rec.replace_calls, n_calls)
  # get_tls_from_u None is taken as True by _analyze_aniso_replace_with_supplied
  assert rec.kw.get('get_tls_from_u') in (None, True), \
    rec.kw.get('get_tls_from_u')
  tgi = rec.kw['tlso_group_info']
  # No coordinate shift between TLS and map (sites used as they are)
  assert not [x for x in tgi.tlso_shift_cart_list if x], \
    tgi.tlso_shift_cart_list
  assert len(rec.selections) == 2
  assert sorted(rec.masks.keys()) == sorted(rec.selections), (
    sorted(rec.masks.keys()))
  aniso_b_cart = rec.kw.get('aniso_b_cart')
  b_iso = rec.kw.get('b_iso')
  if aniso_b_cart:
    overall_u = flex.double(adptbx.b_as_u(tuple(flex.double(aniso_b_cart) -
      flex.double((b_iso, b_iso, b_iso, 0, 0, 0)))))
  else:
    overall_u = flex.double(6, 0)
  return overall_u, get_ownership(rec, model_points = model_points)

def expected_u(rec, selection, xyz, overall_u):
  '''U from the TLS of the group with this selection at xyz, minus the
     overall U (computed here from the TLS definitions of the test)'''
  from mmtbx_tls_ext import uaniso_from_tls_one_group
  t = rec.group_definitions[selection].tlso
  return flex.double(uaniso_from_tls_one_group(tlso = t,
     sites_cart = flex.vec3_double([xyz]), zeroize_trace = False)[0]
     ) - overall_u

def check_entry_u(rec, e, selection, kind, overall_u):
  other = [s for s in rec.group_definitions.keys() if s != selection][0]
  own_u = expected_u(rec, selection, e.xyz, overall_u)
  other_u = expected_u(rec, other, e.xyz, overall_u)
  assert e.uu is not None, (selection, kind, e.xyz)
  assert approx_equal(list(e.uu), list(own_u), eps = U_EPS, out = None), (
    "%s entry of '%s' at %s: U %s expected %s (other group: %s)" %(
    kind, selection, str(e.xyz), list(e.uu), list(own_u), list(other_u)))
  assert flex.max(flex.abs(flex.double(e.uu) - other_u)) > U_DIFFERENT, (
    selection, kind, e.xyz, list(e.uu), list(other_u))

def check_all_owned(rec, overall_u, owner):
  counts = {}
  for e in rec.entries:
    if e.index in owner:
      o = owner[e.index]
      check_entry_u(rec, e, o.selection, o.kind, overall_u)
      k = (o.selection, o.kind)
      counts[k] = counts.get(k, 0) + 1
  return counts

def tst_two_tls_groups():
  '''With two TLS groups and get_tls_from_u, uu_b_cart_as_u_cart of every
     owned entry equals the U from its own group's TLS at the entry's site
     minus the overall U (and differs from the other group's).  Model-point
     entries keep their own group's TLS also when they lie inside the
     mask of the other group.'''
  rec = run_two_groups()
  overall_u, owner = check_tls_run(rec)
  counts = check_all_owned(rec, overall_u, owner)
  n_in_other = {}
  for e in rec.entries:
    o = owner.get(e.index)
    if o and o.kind == 'model':
      other = [s for s in rec.selections if s != o.selection][0]
      if inside_mask(rec.masks[other], e.xyz):
        n_in_other[o.selection] = n_in_other.get(o.selection, 0) + 1
  print("Entries checked (selection, kind): %s" %(
     ", ".join(["'%s' %s: %s" %(k[0], k[1], counts[k])
       for k in sorted(counts.keys())])))
  print("Model-point entries inside the other group's mask: %s" %(
     ", ".join(["'%s': %s" %(s, n_in_other.get(s, 0))
        for s in rec.selections])))
  for s in rec.selections:
    assert counts.get((s, 'model')), (s, counts)
  assert [1 for k in counts.keys() if k[1] == 'grid'], counts
  # Required for the check above to cover model points in the other mask
  assert sum(n_in_other.values()) > 0, n_in_other

def tst_overlap_nearest():
  '''An original grid entry inside both masks takes the TLS of the group
     with the nearest selected atom, in both group orders; in at least one
     order that is not the group processed last (so "later group wins"
     would give a different result)'''
  n_differs_from_later = 0
  for reverse in (False, True):
    rec = run_two_groups(reverse = reverse)
    overall_u, owner = check_tls_run(rec)
    both = [e for e in rec.entries if e.index in owner and
       owner[e.index].kind == 'grid' and len(owner[e.index].in_masks) == 2]
    print("Order %s: original entries inside both masks: %s" %(
      rec.selections, len(both)))
    assert len(both) > 0, "test setup: no original entry inside both masks"
    for e in both:
      o = owner[e.index]
      check_entry_u(rec, e, o.selection, 'grid (both masks)', overall_u)
      later = o.in_masks[-1]
      print("  entry at (%.2f, %.2f, %.2f): nearest atom distances %s; "
        "nearest group '%s', group processed last '%s'" %(
        e.xyz + (", ".join(["'%s' %.2f" %(s, nearest_distance(rec.sites[s],
          e.xyz)) for s in o.in_masks]), o.selection, later)))
      if later != o.selection:
        n_differs_from_later += 1
  assert n_differs_from_later > 0, \
    "test setup: nearest group is always the group processed last"

def tst_reversed_order():
  '''Reversing the order of the two groups gives the same U for every
     entry (original grid entries and model points) and the same set of
     entries.  Matching between runs: original entries are matched by
     index (the original xyz lists must be identical in the two runs);
     model-point entries are placed exactly at atom sites of their own
     selection and are matched by (selection, xyz).'''
  compare_orders(run_two_groups(reverse = False),
    run_two_groups(reverse = True))

def compare_orders(fwd, rev):
  '''Same U for every entry and same set of entries in two runs with the
     groups in opposite orders (see tst_reversed_order)'''
  u_fwd, own_fwd = check_tls_run(fwd)
  u_rev, own_rev = check_tls_run(rev)
  assert rev.selections == list(reversed(fwd.selections))
  assert approx_equal(list(u_fwd), list(u_rev), out = None)
  assert fwd.n_original == rev.n_original
  def keyed(rec, owner):
    result = {}
    for e in rec.entries:
      if e.index < rec.n_original:
        result[('original', e.index)] = e
      else:
        result[(owner[e.index].selection, e.xyz)] = e
    return result
  k_fwd = keyed(fwd, own_fwd)
  k_rev = keyed(rev, own_rev)
  n_compared = {'original': 0, 'model': 0}
  for k in sorted(set(k_fwd.keys()) & set(k_rev.keys()), key = str):
    e_fwd, e_rev = k_fwd[k], k_rev[k]
    if k[0] == 'original':
      assert approx_equal(e_fwd.xyz, e_rev.xyz, out = None), (e_fwd.xyz,
         e_rev.xyz)
    assert (e_fwd.uu is None) == (e_rev.uu is None), (k, e_fwd.uu, e_rev.uu)
    if e_fwd.uu is not None:
      assert approx_equal(list(e_fwd.uu), list(e_rev.uu), eps = U_EPS,
        out = None), (k, list(e_fwd.uu), list(e_rev.uu))
    n_compared['original' if k[0] == 'original' else 'model'] += 1
  print("Entries compared: original %s, model points %s" %(
    n_compared['original'], n_compared['model']))
  for s in fwd.selections:
    n_f = len([k for k in k_fwd.keys() if k[0] == s])
    n_r = len([k for k in k_rev.keys() if k[0] == s])
    print("'%s': model-point entries forward %s reversed %s" %(s, n_f, n_r))
  assert set(k_fwd.keys()) == set(k_rev.keys()), (
    len(k_fwd), len(k_rev), len(set(k_fwd.keys()) & set(k_rev.keys())))
  assert n_compared['original'] > 0 and n_compared['model'] > 0
  return n_compared

def tst_sampling_count():
  '''The model points added for each group are its selected atom sites
     taken every i_ratio-th, with i_ratio = max(1, int((1 + n_atoms) /
     max(25, n_original))) from the number of original entries, in both
     group orders.  Under the old rule (growing entry count, including
     model points added for groups processed earlier) at least one group
     in one order would get a different number (so the test
     discriminates).'''
  n_old_differs = 0
  counts = {}
  for reverse in (False, True):
    rec = run_two_groups(reverse = reverse)
    overall_u, owner = check_tls_run(rec)
    assert rec.n_original < MIN_MODEL_POINTS, rec.n_original
    n_entries_old_rule = rec.n_original
    for s in rec.selections:
      sites = rec.sites[s]
      n_atoms = sites.size()
      i_ratio = max(1, int((1 + n_atoms) / max(MIN_MODEL_POINTS,
         rec.n_original)))
      expected_sites = set([tuple(sites[i]) for i in range(0, n_atoms,
         i_ratio)])
      old_ratio = max(1, int((1 + n_atoms) / max(MIN_MODEL_POINTS,
         n_entries_old_rule)))
      n_old = len(range(0, n_atoms, old_ratio))
      got_sites = set([e.xyz for e in rec.entries if e.index in owner and
         owner[e.index].kind == 'model' and owner[e.index].selection == s])
      print("Order %s: '%s' atoms %s model points %s expected %s "
        "(old rule: %s)" %(rec.selections, s, n_atoms, len(got_sites),
        len(expected_sites), n_old))
      assert got_sites == expected_sites, (s, len(got_sites),
         len(expected_sites))
      counts.setdefault(s, []).append(len(got_sites))
      if n_old != len(expected_sites):
        n_old_differs += 1
      n_entries_old_rule += n_old
  for s, c in counts.items():
    assert len(set(c)) == 1, (s, c)
  assert n_old_differs > 0, \
    "test setup: old and new sampling rules give the same numbers"

def tst_empty_mask_group():
  '''A group whose mask contains none of the original grid entries gets
     model-point entries carrying its own TLS, not that of the group
     processed before it.  Real configuration, no wrapper: split B in
     reversed order, so the group of the first two residues is processed
     second and the only original grid entry is not inside its mask.'''
  rec = run_two_groups(split = SPLIT_B, reverse = True)
  overall_u, owner = check_tls_run(rec)
  s = rec.selections[1]
  n_original_inside = len([e for e in rec.entries if
    e.index < rec.n_original and inside_mask(rec.masks[s], e.xyz)])
  model_entries = [e for e in rec.entries if e.index in owner and
    owner[e.index].kind == 'model' and owner[e.index].selection == s]
  print("Group processed second: '%s'; original entries inside its mask: "
     "%s; model-point entries: %s" %(s, n_original_inside,
     len(model_entries)))
  assert n_original_inside == 0, (
    "test setup: expected no original entries inside mask of '%s'" %(s),
     n_original_inside)
  assert model_entries
  for e in model_entries:
    check_entry_u(rec, e, s, 'model', overall_u)

# ---------------------------------------------------------------------------
# More rules for original entries inside more than one mask: an exact tie
#  and no model (the group listed first in both), and ownership with many
#  original grid entries.

def tie_ranges():
  '''Overlapping residue ranges (first:r) and (r:last), where r is the
     residue of the atom nearest to the first original grid entry of the
     split A run (original grid entries depend on the map and n_boxes, not
     on the groups; checked in tst_tie_first_listed).  That atom is selected
     by both groups, so for that entry the nearest selected atom of each
     group is the same atom and the two distances are exactly equal.
     Returns the ranges and the original xyz list of the split A run.'''
  probe = run_two_groups()
  probe_xyz = [e.xyz for e in probe.entries[:probe.n_original]]
  mmm = get_mmm()
  atoms = mmm.model().get_hierarchy().atoms()
  dist = (atoms.extract_xyz() - probe_xyz[0]).norms()
  i_nearest = list(dist).index(flex.min(dist))
  r = atoms[i_nearest].parent().parent().resseq_as_int()
  resseqs = get_resseqs(mmm)
  return ((resseqs[0], r), (r, resseqs[-1])), probe_xyz

def tst_tie_first_listed():
  '''Two groups on overlapping residue ranges sharing the residue of the
     atom nearest to an original grid entry (tie_ranges): that entry is
     inside both masks at exactly equal nearest-atom distance and takes the
     TLS of the group listed first, so the other group's TLS when the order
     is reversed.  Other original entries follow the nearest-atom rule.
     Model-point entries are not checked here (the selections overlap, so
     a model point can be at an atom of both groups).'''
  ranges, probe_xyz = tie_ranges()
  def setup():
    return model_sharpen_setup(
      get_groups_f = lambda mmm: get_groups_from_ranges(mmm, ranges))
  n_ties = []
  for reverse in (False, True):
    rec = run_tls_groups(('tie', ranges), setup, reverse = reverse)
    overall_u, owner = check_tls_run(rec, model_points = False)
    # Same original grid entries as in the run used to choose the ranges
    assert rec.n_original == len(probe_xyz), (rec.n_original,
       len(probe_xyz))
    for e in rec.entries[:rec.n_original]:
      assert approx_equal(e.xyz, probe_xyz[e.index], out = None), (e.xyz,
        probe_xyz[e.index])
    ties = []
    for e in rec.entries[:rec.n_original]:
      in_masks = [s for s in rec.selections if inside_mask(rec.masks[s],
         e.xyz)]
      d = [nearest_distance(rec.sites[s], e.xyz) for s in rec.selections]
      print("Order %s: original entry at (%.2f, %.2f, %.2f) inside masks of "
        "%s; nearest atom distances %s" %((rec.selections,) + e.xyz +
        (in_masks, ", ".join(["%.6f" %(x) for x in d]))))
      if len(in_masks) == 2 and d[0] == d[1]:
        ties.append(e)
    assert ties, ("test setup: no original entry inside both masks at "
      "exactly equal nearest-atom distance")
    for e in ties:
      # The group listed first in this run
      check_entry_u(rec, e, rec.selections[0], 'grid (tie)', overall_u)
    counts = check_all_owned(rec, overall_u, owner)
    print("Order %s: entries at an exact tie: %s; original entries checked "
       "%s" %(rec.selections, len(ties), sum(counts.values())))
    n_ties.append(len(ties))
  assert n_ties[0] == n_ties[1], n_ties

NO_MODEL_NOISE = 0.3 # uniform noise added to each half map (x rms of map)

def no_model_setup(n_boxes = 8):
  '''half_map_sharpen in a map_model_manager with two half maps and no
     model.  Half maps: the map of get_mmm() plus independent uniform
     noise (own generator with a fixed seed, so the global random state
     used by other code is not changed).  The TLS groups (get_groups
     split A) are defined from the model of get_mmm(), which is then not
     used.'''
  mmm = get_mmm()
  groups = get_groups(mmm, SPLIT_A)
  mm = mmm.map_manager()
  rms = mm.map_data().as_1d().standard_deviation_of_the_sample()
  import random
  generator = random.Random(7)
  half_maps = []
  for i in range(2):
    m = mm.deep_copy()
    data = m.map_data()
    noise = (flex.double([generator.random() for j in range(data.size())])
      - 0.5) * (NO_MODEL_NOISE * rms)
    noise.reshape(data.accessor())
    m.set_map_data(data + noise)
    half_maps.append(m)
  new_mmm = map_model_manager(map_manager_1 = half_maps[0],
     map_manager_2 = half_maps[1])
  new_mmm.set_log(StringIO())
  new_mmm.set_resolution(3)
  assert new_mmm.get_model_by_id('model') is None
  def sharpen(tlso_group_info):
    new_mmm.half_map_sharpen(local_sharpen = True, anisotropic_sharpen = True,
      tlso_group_info = tlso_group_info, n_boxes = n_boxes, nproc = 1)
  return new_mmm, groups, sharpen

def tst_no_model_first_listed():
  '''Supplied TLS without a model (half_map_sharpen): each group's mask
     is the mask around all density, so an original entry inside a mask is
     inside every group's mask and takes the TLS of the group listed first
     (the other group's when the order is reversed); no model-point
     entries are added.'''
  n_inside = []
  for reverse in (False, True):
    rec = run_tls_groups(('no_model',), no_model_setup, reverse = reverse)
    assert rec.has_model is False, rec.has_model
    overall_u, owner = check_tls_run(rec)
    assert len(rec.entries) == rec.n_original, (len(rec.entries),
       rec.n_original)
    inside = []
    for e in rec.entries:
      in_masks = [s for s in rec.selections if inside_mask(rec.masks[s],
         e.xyz)]
      assert len(in_masks) in (0, len(rec.selections)), in_masks
      if in_masks:
        inside.append(e)
    print("Order %s: original entries %s, inside the masks %s, entries "
       "after the TLS step %s" %(rec.selections, rec.n_original,
       len(inside), len(rec.entries)))
    assert inside, "test setup: no original entry inside the masks"
    for e in inside:
      check_entry_u(rec, e, rec.selections[0], 'grid (no model)', overall_u)
    n_inside.append(len(inside))
  assert n_inside[0] == n_inside[1], n_inside

MULTI_N_BOXES = 64
MULTI_SPLIT = 5

def run_multi(reverse):
  return run_tls_groups(('multi', MULTI_SPLIT, MULTI_N_BOXES),
    lambda: model_sharpen_setup(
      get_groups_f = lambda mmm: get_groups(mmm, MULTI_SPLIT),
      n_boxes = MULTI_N_BOXES), reverse = reverse)

def tst_multi_grid_ownership():
  '''model_sharpen with n_boxes 64 (several original grid entries, some
     inside both masks): every owned entry has the U of its group by the
     nearest-atom rule (get_ownership), and every entry has the same U with
     the group order reversed'''
  runs = []
  n_differs_from_later = 0
  for reverse in (False, True):
    rec = run_multi(reverse)
    overall_u, owner = check_tls_run(rec)
    counts = check_all_owned(rec, overall_u, owner)
    grid = [owner[i] for i in sorted(owner.keys()) if owner[i].kind == 'grid']
    both = [o for o in grid if len(o.in_masks) > 1]
    n_differs_from_later += len([o for o in both
       if o.selection != o.in_masks[-1]])
    print("Order %s: original entries %s, owned %s, inside both masks %s "
       "(owner %s); model-point entries %s" %(rec.selections, rec.n_original,
       len(grid), len(both), ", ".join(["'%s'" %(o.selection) for o in both]),
       len(rec.entries) - rec.n_original))
    print("  Entries checked (selection, kind): %s" %(
       ", ".join(["'%s' %s: %s" %(k[0], k[1], counts[k])
         for k in sorted(counts.keys())])))
    assert rec.n_original >= 5, rec.n_original
    assert len(both) >= 2, len(both)
    runs.append(rec)
  compare_orders(runs[0], runs[1])
  print("Entries inside both masks whose nearest group is not the group "
     "processed last: %s" %(n_differs_from_later))
  assert n_differs_from_later > 0, \
    "test setup: nearest group is always the group processed last"

# ---------------------------------------------------------------------------
# Supplied TLS through external_sharpen and half_map_sharpen on a manager
#  with the model 'model' (group masks around each group's atoms, nearest-
#  atom rule and model points, as in model_sharpen), and the same
#  tlso_group_info object converted more than once (used twice, and
#  model_sharpen global-local-global).

def check_with_model(rec, n_calls = 1):
  '''The TLS step (last call) sees the model; every owned entry has its own
     group's U (nearest-atom rule, get_ownership); each group has
     model-point entries and at least one original grid entry is owned.
     Returns the counts of check_all_owned.'''
  assert rec.has_model is True, rec.has_model
  overall_u, owner = check_tls_run(rec, n_calls = n_calls)
  counts = check_all_owned(rec, overall_u, owner)
  print("Order %s: original entries %s, entries after the TLS step %s; "
     "entries checked (selection, kind): %s" %(rec.selections,
     rec.n_original, len(rec.entries),
     ", ".join(["'%s' %s: %s" %(k[0], k[1], counts[k])
       for k in sorted(counts.keys())])))
  for s in rec.selections:
    assert counts.get((s, 'model')), (s, counts)
  assert [1 for k in counts.keys() if k[1] == 'grid'], counts
  return counts

def check_tl_calls(rec, n_calls = None):
  '''Every call of the TLS step received both groups (in the order of the
     run) with their own T and L'''
  if n_calls is not None:
    assert len(rec.tl_calls) == n_calls, (len(rec.tl_calls), n_calls)
  assert rec.tl_calls
  for i, c in enumerate(rec.tl_calls):
    print("TLS step call %s: selections %s" %(i + 1, c.selections))
    assert c.selections == rec.selections, (i, c.selections, rec.selections)
    assert len(c.T_list) == 2 and len(c.L_list) == 2, (i, c.T_list,
       c.L_list)
    for s, T, L in zip(c.selections, c.T_list, c.L_list):
      t = rec.group_definitions[s].tlso
      assert approx_equal(T, t.t, out = None), (i, s, T, t.t)
      assert approx_equal(L, t.l, out = None), (i, s, L, t.l)

def external_setup():
  '''external_sharpen of get_mmm() (with its model 'model') against an
     'external_map' that is an identical deep copy of the map.  A differing
     (for example noisy) external map hits a separate, older crash in
     cctbx/maptbx/refine_sharpening.py: for external-map sharpening the
     overall cc is fixed at 1, so the overall ssqr is zero and bb_values are
     infinite (deepcopy in _average_scale_factor_info_over_xyz then fails).
     So this case tests the TLS conversion and group rules only, not
     realistic external sharpening.  TLS groups: get_groups split A.'''
  mmm = get_mmm()
  groups = get_groups(mmm, SPLIT_A)
  mmm.add_map_manager_by_id(map_manager = mmm.map_manager().deep_copy(),
     map_id = 'external_map')
  assert mmm.get_model_by_id('model') is not None
  def sharpen(tlso_group_info):
    mmm.external_sharpen(local_sharpen = True, anisotropic_sharpen = True,
      tlso_group_info = tlso_group_info, n_boxes = 8, nproc = 1)
  return mmm, groups, sharpen

def tst_external_tls():
  '''external_sharpen with supplied TLS (group_args with tlso_list and
     selection_as_text_list) and an 'external_map' (identical copy of the
     map, see external_setup) runs; the manager has a
     model, so every entry owned by the TLS step has the U of its group by
     the nearest-atom rule, and each group gets model-point entries (the
     checks of the model_sharpen cases)'''
  rec = run_tls_groups(('external',), external_setup)
  check_tl_calls(rec, n_calls = 1)
  check_with_model(rec)

HALF_MAP_N_BOXES = 8
HALF_MAP_SPLIT = SPLIT_A

def half_map_with_model_setup(n_boxes = HALF_MAP_N_BOXES,
    split = HALF_MAP_SPLIT):
  '''no_model_setup (two half maps; n_boxes) with the model of get_mmm()
     added as 'model' (same map grid and origin); TLS groups get_groups
     (split)'''
  new_mmm, groups, sharpen = no_model_setup(n_boxes = n_boxes)
  model = get_mmm().model()
  new_mmm.add_model_by_id(model = model, model_id = 'model')
  assert new_mmm.get_model_by_id('model') is not None
  if split != SPLIT_A:
    groups = get_groups(new_mmm, split)
  return new_mmm, groups, sharpen

def run_half_map_with_model(reverse, n_boxes = HALF_MAP_N_BOXES,
    split = HALF_MAP_SPLIT):
  return run_tls_groups(('half_map_model', n_boxes, split),
    lambda: half_map_with_model_setup(n_boxes = n_boxes, split = split),
    reverse = reverse)

def tst_half_map_with_model():
  '''half_map_sharpen with supplied TLS on a manager with two half maps and
     the model 'model': the TLS step sees the model, the two groups' masks
     differ, every owned entry has its own group's U by the nearest-atom
     rule, each group gets model-point entries, and the result is the same
     with the group order reversed.  At least one original entry is inside
     both masks, and for at least one such entry (over both orders) the
     nearest group is not the group processed last (so "later group wins"
     would give a different result).'''
  runs = []
  n_both = []
  n_differs_from_later = 0
  for reverse in (False, True):
    rec = run_half_map_with_model(reverse)
    check_tl_calls(rec, n_calls = 1)
    check_with_model(rec)
    m0, m1 = [(rec.masks[s].map_data().as_1d().as_double() > 0.5).as_int()
       for s in rec.selections]
    n_inside = [m0.count(1), m1.count(1)]
    n_diff = (m0 != m1).count(True)
    print("Order %s: grid points inside the masks %s; inside one mask only "
       "%s of %s" %(rec.selections, n_inside, n_diff, m0.size()))
    assert min(n_inside) > 0, n_inside
    assert n_diff > 0, "masks of the two groups are identical"
    overall_u, owner = check_tls_run(rec)
    both = [o for i, o in sorted(owner.items()) if o.kind == 'grid' and
       len(o.in_masks) > 1]
    n_not_later = len([o for o in both if o.selection != o.in_masks[-1]])
    print("Order %s: original entries %s, inside both masks %s (owner %s); "
       "nearest group not the group processed last: %s" %(rec.selections,
       rec.n_original, len(both),
       ", ".join(["'%s'" %(o.selection) for o in both]), n_not_later))
    n_both.append(len(both))
    n_differs_from_later += n_not_later
    runs.append(rec)
  assert min(n_both) > 0, ("test setup: no original entry inside both masks",
     n_both)
  assert n_differs_from_later > 0, \
    "test setup: nearest group is always the group processed last"
  compare_orders(runs[0], runs[1])

def reused_setup():
  '''model_sharpen (as model_sharpen_setup, split A) run twice with the same
     tlso_group_info object'''
  mmm, groups, sharpen = model_sharpen_setup(
    get_groups_f = lambda mmm: get_groups(mmm, SPLIT_A))
  def sharpen_twice(tlso_group_info):
    sharpen(tlso_group_info)
    sharpen(tlso_group_info)
  return mmm, groups, sharpen_twice

def tst_tls_reused():
  '''The same tlso_group_info object passed to model_sharpen twice: both
     runs complete, and the TLS step of each (in particular the second)
     receives both groups' T and L; the second run's entries follow the
     nearest-atom rule'''
  rec = run_tls_groups(('reused',), reused_setup)
  check_tl_calls(rec, n_calls = 2)
  check_with_model(rec, n_calls = 2)

def global_local_global_setup():
  '''model_sharpen as model_sharpen_setup (split A) with
     overall_sharpen_before_and_after_local'''
  mmm = get_mmm()
  def sharpen(tlso_group_info):
    mmm.model_sharpen(local_sharpen = True, anisotropic_sharpen = True,
      overall_sharpen_before_and_after_local = True,
      tlso_group_info = tlso_group_info,
      replace_aniso_with_tls_equiv = True,
      n_boxes = 8, nproc = 1)
  return mmm, get_groups(mmm, SPLIT_A), sharpen

def tst_global_local_global_model():
  '''model_sharpen global-local-global with supplied TLS (the object is
     copied and converted again by the inner runs) completes; every call
     of the TLS step receives both groups' T and L, and the entries of the
     last call follow the nearest-atom rule'''
  rec = run_tls_groups(('global_local_global',), global_local_global_setup)
  print("TLS step calls: %s" %(rec.replace_calls))
  check_tl_calls(rec)
  check_with_model(rec, n_calls = rec.replace_calls)

def glg_setup(setup, method_name):
  '''(mmm, groups, sharpen) from setup(), with sharpen replaced by the
     method method_name of mmm (for example half_map_sharpen) run
     global-local-global (local_sharpen, anisotropic_sharpen and
     overall_sharpen_before_and_after_local) with supplied TLS'''
  mmm, groups, unused_sharpen = setup()
  def sharpen(tlso_group_info):
    getattr(mmm, method_name)(local_sharpen = True,
      anisotropic_sharpen = True,
      overall_sharpen_before_and_after_local = True,
      tlso_group_info = tlso_group_info, n_boxes = 8, nproc = 1)
  return mmm, groups, sharpen

def tst_global_local_global_half_map():
  '''half_map_sharpen global-local-global with supplied TLS
     (half_map_with_model_setup) completes; every call of the TLS step
     receives both groups' T and L, and the entries of the last call follow
     the nearest-atom rule'''
  rec = run_tls_groups(('global_local_global_half_map',),
    lambda: glg_setup(half_map_with_model_setup, 'half_map_sharpen'))
  print("TLS step calls: %s" %(rec.replace_calls))
  check_tl_calls(rec)
  check_with_model(rec, n_calls = rec.replace_calls)

# ---------------------------------------------------------------------------
# Synthetic averaging over locations

TSF = (1., 0.9, 0.8, 0.7)    # target_scale_factors
CC = (0.9, 0.8, 0.7, 0.6)    # cc_list
RMS_FO = (100., 80., 50., 30.)
OVERALL_SCALE = (1.2, 1.0, 0.7, 0.5)
U_CART = (1., 2., 3., 0.1, 0.2, 0.3)

def scaled(values, scale):
  return flex.double([scale * x for x in values])

def make_si(scale = 1.):
  sthol2 = flex.double([0.25 / (d * d) for d in (8., 6., 4., 3.)])
  return group_args(
    group_args_type = 'synthetic scaling_info',
    target_sthol2 = sthol2,
    cc_list = scaled(CC, scale),
    rms_fo_list = scaled(RMS_FO, scale),
    rms_fc_list = flex.double([110., 85., 55., 32.]),
    target_scale_factors = scaled(TSF, scale),
   )

def make_location(scale = 1., overall_si = True):
  '''Synthetic complete location; all values are proportional to scale'''
  return group_args(
    group_args_type = 'synthetic scaling_group_info',
    direction_vectors = [(1., 0., 0.), (0., 1., 0.)],
    scaling_info_list = [make_si(scale = scale) for i in range(2)],
    overall_si = make_si(scale = scale) if overall_si else None,
    overall_scale = scaled(OVERALL_SCALE, scale),
    aa_b_cart_as_u_cart = tuple(scaled(U_CART, 1. * scale)),
    bb_b_cart_as_u_cart = tuple(scaled(U_CART, 2. * scale)),
    ss_b_cart_as_u_cart = tuple(scaled(U_CART, 3. * scale)),
    uu_b_cart_as_u_cart = tuple(scaled(U_CART, 4. * scale)),
    fo_b_cart_as_u_cart = tuple(scaled(U_CART, 5. * scale)),
   )

XYZ = [(10., 11., 12.), (20., 21., 22.), (30., 31., 32.)]

def make_scale_factor_info(location_list):
  return group_args(
    group_args_type = 'synthetic scale_factor_info',
    value_list = location_list,
    xyz_list = flex.vec3_double(XYZ[:len(location_list)]),
   )

def average(location_list):
  mmm = map_model_manager()
  mmm.set_log(StringIO())
  sfi = make_scale_factor_info(location_list)
  return mmm._average_scale_factor_info_over_xyz(sfi).value_list[0]

def check_values(values, expected, label):
  assert values is not None, (label, values)
  assert approx_equal(list(values), list(expected)), (label, list(values),
    list(expected))

def tst_average_overall_scale_stored():
  '''Two complete locations (scales 1 and 3): averaged overall_scale is
     their mean (scale 2)'''
  avg = average([make_location(scale = 1.), make_location(scale = 3.)])
  print("averaged overall_scale:", avg.overall_scale if avg.overall_scale
    is None else list(avg.overall_scale))
  check_values(avg.overall_scale, scaled(OVERALL_SCALE, 2.), 'overall_scale')

def tst_mixed_overall_si():
  '''Middle of three locations has overall_si None: it is skipped and the
     averages are over the other two (scales 1 and 4 -> 2.5).  If all
     locations have overall_si None, averaging completes.'''
  avg = average([make_location(scale = 1.),
     make_location(scale = 2., overall_si = False),
     make_location(scale = 4.)])
  for si in avg.scaling_info_list:
    check_values(si.target_scale_factors, scaled(TSF, 2.5),
       'target_scale_factors')
    check_values(si.cc_list, scaled(CC, 2.5), 'cc_list')
    check_values(si.rms_fo_list, scaled(RMS_FO, 2.5), 'rms_fo_list')
  check_values(avg.overall_si.target_scale_factors, scaled(TSF, 2.5),
    'overall_si.target_scale_factors')
  check_values(avg.overall_scale, scaled(OVERALL_SCALE, 2.5),
    'overall_scale')
  check_values(avg.ss_b_cart_as_u_cart, scaled(U_CART, 3. * 2.5),
    'ss_b_cart_as_u_cart')

  # All locations without overall_si
  avg = average([make_location(scale = 1., overall_si = False),
     make_location(scale = 2., overall_si = False)])
  print("All overall_si None: scaling_info_list=%s overall_si=%s "
     "overall_scale=%s direction_vectors=%s" %(avg.scaling_info_list,
      avg.overall_si, avg.overall_scale, avg.direction_vectors))
  assert avg.scaling_info_list is None, avg.scaling_info_list
  assert avg.overall_si is None, avg.overall_si
  assert avg.overall_scale is None, avg.overall_scale

CC_A = (0.9, 0.8, 0.6, 0.4)
CC_B = (0.7, 0.5, 0.3, 0.2)

def tst_averaged_ssqr():
  '''Two locations with different cc_list (CC*) and ssqr_values derived
     from it as the producer does (ssqr = 2 * (1/cc - 1)): the averaged
     ssqr_values follow the averaged CC*.  Then the get_tls_from_u update
     (with zero U, so that only the error factor 1/(1 + ssqr) acts; this
     update does not apply overall_scale) gives target_scale_factors of
     1/(1 + ssqr) = cc / (2 - cc) at the mean cc'''
  def location(cc):
    loc = make_location(scale = 1.)
    for si in loc.scaling_info_list + [loc.overall_si]:
      si.cc_list = flex.double(cc)
      si.add(key = 'ssqr_values', value = flex.double(
        [2. * (1. / x - 1.) for x in cc]))
    return loc
  avg = average([location(CC_A), location(CC_B)])
  mean_cc = [0.5 * (a + b) for a, b in zip(CC_A, CC_B)]
  expected_ssqr = [2. * (1. / c - 1.) for c in mean_cc]
  first_ssqr = [2. * (1. / c - 1.) for c in CC_A]
  for si in avg.scaling_info_list:
    print("averaged ssqr_values:", list(si.ssqr_values),
       "expected:", expected_ssqr)
    assert approx_equal(list(si.cc_list), mean_cc, out = None), list(
       si.cc_list)
    assert approx_equal(list(si.ssqr_values), expected_ssqr, out = None), (
      list(si.ssqr_values), expected_ssqr, "first location:", first_ssqr)

  # Resolution dependence of the get_tls_from_u update for this entry
  mmm = get_mmm()
  avg.uu_b_cart_as_u_cart = (0., 0., 0., 0., 0., 0.)
  assert avg.overall_scale is not None
  sfi = group_args(group_args_type = 'synthetic scale_factor_info',
    value_list = [avg], xyz_list = flex.vec3_double([XYZ[0]]),
    d_min = 3.)
  mmm._update_scale_factor_info_from_aniso(sfi, get_tls_from_u = True)
  expected_factors = [c / (2. - c) for c in mean_cc]  # 1/(1 + ssqr)
  for si in avg.scaling_info_list:
    print("updated target_scale_factors:", list(si.target_scale_factors),
       "expected:", expected_factors)
    assert approx_equal(list(si.target_scale_factors), expected_factors,
       eps = 1.e-4, out = None), (list(si.target_scale_factors),
        expected_factors)

# ---------------------------------------------------------------------------
# No TLS: local anisotropic sharpening is unchanged.  Expected fingerprint
#  recorded identically with the code before and after supplied TLS was
#  made to take effect.  Tolerance 1e-4 absolute: values are of order
#  0.01-1.5 and were recorded to 6 decimals; this allows for platform /
#  compiler floating-point differences while still detecting any real
#  change in the sharpening.

EXPECTED_SIZE = 38400
EXPECTED_RMS = 0.151635
EXPECTED_MIN = -0.268431
EXPECTED_MAX = 1.467591
EXPECTED_SAMPLE = (-0.016822, -0.055808, -0.043456, -0.036190, -0.027208,
   -0.029304)
FINGERPRINT_EPS = 1.e-4

def tst_no_tls_unchanged():
  '''model_sharpen local anisotropic without TLS gives the recorded
     sharpened-map fingerprint'''
  def sharpen():
    mmm = get_mmm()
    mmm.model_sharpen(local_sharpen = True, anisotropic_sharpen = True,
        n_boxes = 8, nproc = 1)
    return mmm
  mmm = run_in_temp_dir(sharpen)
  data = mmm.map_manager().map_data().as_1d()
  mmm_stats = data.min_max_mean()
  rms = data.standard_deviation_of_the_sample()
  n = data.size()
  sample = [data[i] for i in (0, n // 7, n // 3, n // 2, (2 * n) // 3,
     n - 1)]
  print("FINGERPRINT size %d mean %.6f rms %.6f min %.6f max %.6f "
    "values %s" %(n, mmm_stats.mean, rms, mmm_stats.min, mmm_stats.max,
     " ".join(["%.6f" %(x) for x in sample])))
  assert n == EXPECTED_SIZE, n
  for label, value, expected in (('rms', rms, EXPECTED_RMS),
     ('min', mmm_stats.min, EXPECTED_MIN),
     ('max', mmm_stats.max, EXPECTED_MAX)):
    assert approx_equal(value, expected, eps = FINGERPRINT_EPS), (
      label, value, expected)
  assert approx_equal(sample, EXPECTED_SAMPLE, eps = FINGERPRINT_EPS), (
    sample, EXPECTED_SAMPLE)

def run(args = None):
  # Each case is listed separately in iotbx/run_tests.py; add new ones there
  all_tests = [
     ('options_passed_on', tst_options_passed_on),
     ('all_tls_entries_updated', tst_all_tls_entries_updated),
     ('two_tls_groups', tst_two_tls_groups),
     ('overlap_nearest', tst_overlap_nearest),
     ('reversed_order', tst_reversed_order),
     ('sampling_count', tst_sampling_count),
     ('empty_mask_group', tst_empty_mask_group),
     ('tie_first_listed', tst_tie_first_listed),
     ('no_model_first_listed', tst_no_model_first_listed),
     ('multi_grid_ownership', tst_multi_grid_ownership),
     ('half_map_with_model', tst_half_map_with_model),
     ('tls_reused', tst_tls_reused),
     ('global_local_global_model', tst_global_local_global_model),
     ('global_local_global_half_map', tst_global_local_global_half_map),
     ('external_tls', tst_external_tls),
     ('average_overall_scale_stored', tst_average_overall_scale_stored),
     ('mixed_overall_si', tst_mixed_overall_si),
     ('averaged_ssqr', tst_averaged_ssqr),
     ('no_tls_unchanged', tst_no_tls_unchanged),
    ]
  if args:
    unknown = set(args) - set([name for name, f in all_tests])
    assert not unknown, "Unknown test name(s): %s" %(
      " ".join(sorted(unknown)))
    tests = [(name, f) for name, f in all_tests if name in args]
  else:
    tests = all_tests
  for name, f in tests:
    print("Running %s" %(name))
    f()
  print("OK")

if __name__ == "__main__":
  run(sys.argv[1:])
