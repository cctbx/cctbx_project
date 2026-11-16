from __future__ import absolute_import, division, print_function
import time
import iotbx.pdb
from libtbx import group_args
from scitbx.array_family import flex
from mmtbx.pair_interaction import pair_interaction
from mmtbx.pair_interaction.tst_04 import pdb_str

"""
Interaction search of single-conformer views on a shared grid: the flags of
points_and_pairs (near/far points of a region of differing atoms), the fixed
grid of prepare/run_prepared, and run_prepared_views against independent
per-view runs (exact by construction); the pairs of a whole hierarchy
(interaction_pairs, views_pairs) and the fragment rule applied to them
(fragment_from_pairs). Nothing is written to disk.
"""

ALTLOC_RESIDUES = (4, 12, 21)   # side chains duplicated as A and B
B_ONLY_RESIDUE  = 26            # every atom altloc B: absent from view A
BACKBONE = ("N", "CA", "C", "O", "CB", "H", "HA")

def altloc_fixture(shift=0.7, b_only=True):
  """
  group_args(ph, views): the 2lvr hierarchy of tst_04 with the side chains
  of TYR 4, GLN 12 and ARG 21 as altloc A (original) and B (shifted by
  shift A along x; 0: duplicate conformers) and, with b_only, HIS 26 as
  altloc B only; views: for labels A and B, group_args(ph, iselection,
  residues) of parent.select(blank or label), residues the 1-based parent
  residue group index of each view residue. Built from the PDB text so atom
  order is the parent order.
  """
  lines = []
  for l in pdb_str.splitlines():
    if(not l.startswith("ATOM") and not l.startswith("HETATM")):
      lines.append(l)
      continue
    resseq = int(l[22:26])
    name = l[12:16].strip()
    if(resseq in ALTLOC_RESIDUES and name not in BACKBONE):
      lines.append(l[:16] + "A" + l[17:])
      x = float(l[30:38]) + shift
      lines.append(l[:16] + "B" + l[17:30] + "%8.3f" % x + l[38:])
    elif(resseq == B_ONLY_RESIDUE and b_only):
      lines.append(l[:16] + "B" + l[17:])
    else:
      lines.append(l)
  ph = iotbx.pdb.input(source_info=None,
    lines=flex.std_string(lines)).construct_hierarchy()
  views = []
  for label in ("A", "B"):
    sel = flex.bool([a.parent().altloc.strip() in ("", label)
                     for a in ph.atoms()])
    residues = []
    offset = 0
    for i, rg in enumerate(ph.residue_groups()):
      n = rg.atoms_size()
      if(sel[offset:offset+n].count(True) > 0): residues.append(i+1)
      offset += n
    views.append(group_args(ph=ph.select(sel), iselection=sel.iselection(),
      residues=residues))
  assert ph.atoms_size() > 492
  assert (26 not in views[0].residues) == b_only
  assert 26 in views[1].residues
  return group_args(ph=ph, views=views)

def exercise_flags():
  """points_and_pairs flags through _grid_search: near/far points."""
  f = altloc_fixture()
  p = pair_interaction.prepare(f.views[0].ph)        # view A, no grid
  air = list(p.atom_in_residue)
  wfc = p.element_wfc_dict
  step = 0.3*pair_interaction.A2B
  full = pair_interaction._grid_search(p.ph, air, step, 'sedd', None, wfc)
  assert len(full) > 0 and set(full.values()) == set([0])
  ref = pair_interaction.get_interactions(p.ph, air, silva_type='sedd',
    element_wfc_dict=wfc)
  assert sorted(full) == sorted(ref)                 # bit identity
  # region: the side chain of TYR 4 (parent residue 4) in view A, Bohr
  region = p.xyz.select(p.residue_atoms[3])
  flagged = pair_interaction._grid_search(p.ph, air, step, 'sedd', None,
    wfc, region_sites=region)
  assert sorted(flagged) == sorted(full)             # same pairs
  only_near = set([k for k, v in flagged.items() if v == 1])
  assert len(only_near) > 0
  near = pair_interaction._grid_search(p.ph, air, step, 'sedd', None,
    wfc, region_sites=region, near_only=True)
  assert set(near.values()) == set([1])
  assert only_near <= set(near) < set(full)          # strict subset
  # with a core, as run_prepared uses it
  core = [3, 4, 5]
  c_full = pair_interaction._grid_search(p.ph, air, step, 'sedd', core, wfc)
  c_flag = pair_interaction._grid_search(p.ph, air, step, 'sedd', core, wfc,
    region_sites=region)
  assert sorted(c_full) == sorted(c_flag)
  # Review Focus 5: near_only with nothing near searches nothing
  empty = pair_interaction._grid_search(p.ph, air, step, 'sedd', None, wfc,
    region_sites=flex.vec3_double(), near_only=True)
  assert empty == {}

def exercise_grid():
  """run_prepared on a given grid; atoms outside the grid are an error."""
  f = altloc_fixture()
  p = pair_interaction.prepare(f.views[0].ph)
  air = list(p.atom_in_residue)
  step = 0.3*pair_interaction.A2B
  # a grid equal to the atoms' own box samples the same points
  own = group_args(xyz_min=p.xyz.min(), xyz_max=p.xyz.max())
  d0 = pair_interaction._grid_search(p.ph, air, step, 'sedd', None,
    p.element_wfc_dict)
  d1 = pair_interaction._grid_search(p.ph, air, step, 'sedd', None,
    p.element_wfc_dict, grid=own)
  assert d0 == d1
  # run_prepared: the grid of the core's own sub-selection changes nothing,
  # a wider grid runs and returns the three lists
  core = [3, 4, 5]
  cs = pair_interaction._core_selection(p, core)
  sub = p.xyz.select(cs.selection)
  own = group_args(xyz_min=sub.min(), xyz_max=sub.max())
  r0 = pair_interaction.run_prepared(p, core)
  r1 = pair_interaction.run_prepared(p, core, grid=own)
  assert [sorted(x) for x in r0] == [sorted(x) for x in r1]
  c, a, m = pair_interaction.run_prepared(p, core, grid=own)
  assert len(c) > 0 and set(c) <= set(a) and set(core) <= set(m)
  wide = group_args(xyz_min=tuple([v-7.0 for v in own.xyz_min]),
    xyz_max=tuple([v+7.0 for v in own.xyz_max]))
  c, a, m = pair_interaction.run_prepared(p, core, grid=wide)
  assert len(c) > 0 and set(c) <= set(a) and set(core) <= set(m)
  # Review Focus 1: an atom outside the grid is an error
  small = group_args(xyz_min=own.xyz_min,
    xyz_max=tuple([v-5.0 for v in own.xyz_max]))
  try:
    pair_interaction.run_prepared(p, core, grid=small)
  except AssertionError as e:
    assert str(e).find("grid") > -1, str(e)
  else:
    raise AssertionError("atom outside the grid was accepted")

def mapped_run_prepared(view, core, grid=None):
  """run_prepared of one view, mapped to the parent (the reference)."""
  in_view = dict([(r, k+1) for k, r in enumerate(view.residues)])
  core_v = [in_view[r] for r in core if r in in_view]
  if(len(core_v) == 0): return None
  c, a, m = pair_interaction.run_prepared(view.prepared, core_v, grid=grid)
  isel = view.iselection
  return (sorted([isel[i-1]+1 for i in c]),
          sorted([isel[i-1]+1 for i in a]),
          sorted([view.residues[i-1] for i in m]))

def prepared_views(f):
  return [group_args(prepared=pair_interaction.prepare(v.ph),
    iselection=v.iselection, residues=v.residues) for v in f.views]

CORES = [[1, 2], [3, 4, 5], [13, 14, 15], [24, 25, 26, 27]]

def exercise_views():
  """run_prepared_views equals independent per-view runs on the per-core
  grid (the union box of the views' sub-selections), mapped."""
  views = prepared_views(altloc_fixture())
  for core in CORES:
    for vs in [views, views[::-1]]:          # either view as the reference
      g = pair_interaction.views_grid(vs, core)
      got = pair_interaction.run_prepared_views(vs, core)
      expected = [mapped_run_prepared(v, core, grid=g) for v in vs]
      assert got == expected, (core, got, expected)
  # the altloc cores do distinguish the views (non-vacuous)
  for core in [[3, 4, 5], [24, 25, 26, 27]]:
    a, b = pair_interaction.run_prepared_views(views, core)
    assert a[1] != b[1], core
  # Review Focus 4: one view equals run_prepared on its own box
  assert pair_interaction.run_prepared_views(views[:1], [3, 4, 5]) == \
    [mapped_run_prepared(views[0], [3, 4, 5])]
  # Review Focus 3: a core absent from a view, or from all views
  got = pair_interaction.run_prepared_views(views, [26])   # HIS 26: B only
  assert got[0] is None and got[1] == mapped_run_prepared(views[1], [26])
  assert pair_interaction.run_prepared_views(views, [999]) == [None, None]

def exercise_duplicates():
  """Duplicate conformers (A and B at the same coordinates): every view
  equals the single-conformer run on its own box, bit for bit."""
  views = prepared_views(altloc_fixture(shift=0.0, b_only=False))
  for core in CORES:
    got = pair_interaction.run_prepared_views(views, core)
    expected = [mapped_run_prepared(v, core) for v in views]  # own box
    assert got == expected, (core, got, expected)
    assert got[0][2] == got[1][2], core

def mapped_pairs(view, silva_type, grid=None):
  """interaction_pairs of one view in parent residue indices."""
  res = view.residues
  result = set()
  for a, b in pair_interaction.interaction_pairs(view.prepared, silva_type,
                                                 grid=grid):
    pa, pb = res[a-1], res[b-1]
    result.add((min(pa, pb), max(pa, pb)))
  return result

def exercise_pairs():
  """The pairs of the whole hierarchy: the core-restricted search returns
  exactly the pairs of the unrestricted search that touch the core, and
  fragment_from_pairs applies the two-pass rule of run_prepared to them."""
  ph = iotbx.pdb.input(source_info=None,
    lines=flex.std_string(pdb_str.splitlines())).construct_hierarchy()
  p = pair_interaction.prepare(ph)
  air = list(p.atom_in_residue)
  pairs = {}
  for st in ("sedd", "dori"):
    pairs[st] = pair_interaction.interaction_pairs(p, st)
    assert len(pairs[st]) > 0
    for core in CORES:
      restricted = set(pair_interaction.get_interactions(p.ph, air,
        silva_type=st, core=core, element_wfc_dict=p.element_wfc_dict))
      expected = set([q for q in pairs[st] if q[0] in core or q[1] in core])
      assert restricted == expected, (st, core, restricted, expected)
  for core in CORES:
    new_core = set(core)
    for a, b in pairs["sedd"]:
      if(a in core or b in core): new_core.update((a, b))
    mols = set(new_core)
    for a, b in pairs["dori"]:
      if(a in core or b in core): mols.update((a, b))
    assert len(mols) > len(core), core          # non-vacuous
    got = pair_interaction.fragment_from_pairs(p, core, pairs["sedd"],
                                               pairs["dori"])
    atoms = lambda rs: sorted([i+1 for r in rs for i in p.residue_atoms[r-1]])
    assert got == (atoms(core), atoms(mols), sorted(mols)), (core, got)
  # against the two passes of run_prepared on the whole hierarchy (no core
  # of this peptide has a 15 Bohr sub-selection keeping every atom, so
  # run_prepared itself searches other points): its second pass searches
  # around the new core and keeps the pairs touching the original core
  for core in CORES:
    sedd = pair_interaction.get_interactions(p.ph, air, silva_type="sedd",
      core=core, element_wfc_dict=p.element_wfc_dict)
    new_core = set(core)
    for q in sedd:
      if(set(q) & set(core)): new_core.update(q)
    dori = pair_interaction.get_interactions(p.ph, air, silva_type="dori",
      core=sorted(new_core), element_wfc_dict=p.element_wfc_dict)
    mols = set(new_core)
    for q in dori:
      if(set(q) & set(core)): mols.update(q)
    got = pair_interaction.fragment_from_pairs(p, core, pairs["sedd"],
                                               pairs["dori"])
    assert got[2] == sorted(mols), (core, got[2], sorted(mols))
  # a core residue beyond the hierarchy contributes nothing (run_prepared)
  got = pair_interaction.fragment_from_pairs(p, [1, 999], pairs["sedd"],
                                             pairs["dori"])
  assert got[0] == sorted([i+1 for i in p.residue_atoms[0]])
  assert 999 in got[2] and 1 in got[2]

def exercise_views_pairs():
  """views_pairs equals the union over the views of independent searches
  on the shared grid of all view atoms (views_grid without a core), with
  either view as the reference; duplicate conformers give the pairs of
  the single-conformer hierarchy on its own box, bit for bit."""
  views = prepared_views(altloc_fixture())
  g = pair_interaction.views_grid(views)
  for st in ("sedd", "dori"):
    per_view = [mapped_pairs(v, st, grid=g) for v in views]
    expected = per_view[0] | per_view[1]
    for vs in [views, views[::-1]]:
      got = pair_interaction.views_pairs(vs, st)
      assert got == expected, (st, got ^ expected)
    assert per_view[0] != per_view[1], st             # non-vacuous
    assert len([q for q in expected if 26 in q]) > 0  # HIS 26: B only
  # one view: the search on its own box
  assert pair_interaction.views_pairs(views[:1], "dori") == \
    mapped_pairs(views[0], "dori")
  # a view without atoms contributes nothing
  empty = group_args(prepared=pair_interaction.prepare(
    iotbx.pdb.hierarchy.root()), iselection=flex.size_t(), residues=[])
  assert pair_interaction.views_pairs(views + [empty], "dori") == \
    pair_interaction.views_pairs(views, "dori")
  views = prepared_views(altloc_fixture(shift=0.0, b_only=False))
  for st in ("sedd", "dori"):
    single = mapped_pairs(views[0], st)
    assert mapped_pairs(views[1], st) == single
    assert pair_interaction.views_pairs(views, st) == single, st

def run():
  exercise_flags()
  exercise_grid()
  exercise_views()
  exercise_duplicates()
  exercise_pairs()
  exercise_views_pairs()

if(__name__ == "__main__"):
  t0 = time.time()
  run()
  print("Time: %6.2f" % (time.time()-t0))
  print("OK")
