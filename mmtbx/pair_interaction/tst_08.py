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
per-view runs (exact by construction). Nothing is written to disk.
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

def run():
  exercise_flags()
  exercise_grid()
  exercise_views()
  exercise_duplicates()

if(__name__ == "__main__"):
  t0 = time.time()
  run()
  print("Time: %6.2f" % (time.time()-t0))
  print("OK")
