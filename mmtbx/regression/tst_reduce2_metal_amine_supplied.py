"""reduce2 metal_bound_amines_as_nh2 (A5): supplied H (keep_existing_H and
approach=optimize), reruns and conformers.
Placement and geometry: tst_reduce2_metal_amine.py.
"""
from __future__ import absolute_import, division, print_function
from mmtbx.regression.reduce2_metal_amine_utils import *
from mmtbx.regression.reduce2_metal_amine_utils import shared, shared_text

def test_chain_start_rerun():
  # add -> optimize -> add keep_existing_H=True, and plain add on the output:
  # NH2 every time, at the same positions
  h0, d, l = shared('x1_on')
  previous = 'tst_reduce2_metal_amine_x1_on_out.pdb'
  for i, args in enumerate((['approach=optimize'], ['keep_existing_H=True'],
                            [])):
    with open(previous) as f:
      text = f.read()
    prefix = 'tst_reduce2_metal_amine_gr%d' % (i + 1)
    h, d, l = run(prefix, text, [ON] + args)
    assert h_on(h, 288, 'N') == [('H2', ''), ('H3', '')], (args,
      h_on(h, 288, 'N'))
    for name in ('H2', 'H3'):
      dist = atom(h, 288, name).distance(atom(h0, 288, name))
      assert dist < 1.e-3, (args, name, dist)
    previous = prefix + '_out.pdb'

def test_chain_start_supplied_nh2_all_chain_starts():
  # converted Gly 288 rerun with keep_existing_H=True under all_chain_starts
  # (the setting builds NH3): the supplied NH2 stays, one line, no fallback
  h0, d, l = shared('x1_on')
  text = shared_text('x1_on')
  h, d, l = run('tst_reduce2_metal_amine_gacs1', text,
    [ON, 'keep_existing_H=True', 'n_terminal_charge=all_chain_starts'])
  assert h_on(h, 288, 'N') == [('H2', ''), ('H3', '')], h_on(h, 288, 'N')
  for name in ('H2', 'H3'):
    assert atom(h, 288, name).distance(atom(h0, 288, name)) < 1.e-3, name
  assert section(d, CHAIN_START) == ['N GLY A 288: metal-bound chain start, '
    'unchanged: supplied NH2 kept (H2, H3)'], d
  assert 'not judged' not in d

def test_chain_start_supplied_h():
  # supplied H at a metal-bound chain start: the option builds no NH3 and
  # deletes nothing; every supplied atom keeps its name and coordinates
  one_h = shared_text('x1_off')
  h2 = [l for l in shared_text('x1_on').splitlines() if l[12:16] == ' H2 '][0]
  h1 = [l for l in one_h.splitlines() if l[12:26] == ' H   GLY A 288'][0]
  two_h = one_h.replace(h1, h1 + '\n' + h2)
  for tag, text, names in (('a', one_h, ['H']), ('b', two_h, ['H', 'H2'])):
    h_in = iotbx.pdb.input(source_info=None, lines=text).construct_hierarchy()
    h, d, l = run('tst_reduce2_metal_amine_gs_' + tag, text,
      [ON, 'keep_existing_H=True'])
    off, d_off, l_off = run('tst_reduce2_metal_amine_gs_off', text,
      [OFF, 'keep_existing_H=True'])
    assert [x for x, alt in h_on(h, 288, 'N')] == names, (tag,
      h_on(h, 288, 'N'))
    for name in names:
      dist = atom(h, 288, name).distance(atom(h_in, 288, name))
      assert dist < 1.e-3, (tag, name, dist)
    assert section(d, CHAIN_START) == ['N GLY A 288: metal-bound chain start, '
      'unchanged: supplied H kept (%s)' % ', '.join(names)], (tag, d)
    # every other H as without the option
    a, b = h_xyz(h), h_xyz(off)
    for k in a:
      if ' GLY A 288' in k: continue
      assert k in b, (tag, k)
      dist = (matrix.col(a[k]) - matrix.col(b[k])).length()
      assert dist < 1.e-3, (tag, k, dist)

def test_chain_start_conformer_supplied():
  # a blank N, CA/C/O in A and B, a blank Ni, a supplied H in A only: A keeps
  # it, B gets what the option-off run places
  run('tst_reduce2_metal_amine_gcs0', gly288_conformers(None, False), [OFF])
  with open('tst_reduce2_metal_amine_gcs0_out.pdb') as f:
    text = '\n'.join(l for l in f.read().splitlines()
      if l[12:26] != ' H  BGLY A 288') + '\n'
  h_in = iotbx.pdb.input(source_info=None, lines=text).construct_hierarchy()
  off, d_off, l_off = run('tst_reduce2_metal_amine_gcs_off', text,
    [OFF, 'keep_existing_H=True'])
  h, d, l = run('tst_reduce2_metal_amine_gcs', text,
    [ON, 'keep_existing_H=True'])
  assert atom(h, 288, 'H', 'A').distance(atom(h_in, 288, 'H', 'A')) < 1.e-3
  assert atom(h, 288, 'H', 'B').distance(atom(off, 288, 'H', 'B')) < 1.e-3
  assert section(d, CHAIN_START) == ['N GLY A 288: metal-bound chain start, '
    'unchanged: supplied H in conformer A kept (H); conformer B as '
    'n_terminal_charge places it'], d

def test_supplied():
  # supplied NH3 (reduce2's own output without the option): kept, reported
  nh3 = shared_text('lys_off')
  h, description, log = run('tst_reduce2_metal_amine_nh3_keep', nh3,
    [ON, 'keep_existing_H=True'])
  assert [n for n, alt in h_on(h, 359, 'NZ')] == ['HZ1', 'HZ2', 'HZ3']
  assert section(description, SUPPLIED)[0].startswith(
    'NZ  LYS A 359 (HZ1, HZ2, HZ3)  (NZ  LYS A 359 - CD'), description
  assert CONVERTED not in description

def test_rerun_stable():
  # add -> optimize -> add keep_existing_H=True: same H, same positions; a
  # supplied NH2 is not completed (placement's new HZ1 is dropped)
  h0, description, log = shared('lys_on')
  previous = 'tst_reduce2_metal_amine_lys_on_out.pdb'
  for i, args in enumerate((['approach=optimize'], ['keep_existing_H=True'])):
    prefix = 'tst_reduce2_metal_amine_rr%d' % (i + 1)
    with open(previous) as f:
      text = f.read()
    previous = prefix + '_out.pdb'
    h, description, log = run(prefix, text, [ON] + args)
    off, d_off, l_off = run(prefix + '_off', text, [OFF] + args)
    check_nh2(h, 359, 'NZ', 'CE', 396, 'CD', ['HZ2', 'HZ3'])
    # the NH2 where it was; every other H as without the option
    for name in ('HZ2', 'HZ3'):
      d = atom(h, 359, name).distance(atom(h0, 359, name))
      assert d < 1.e-3, (args, name, d)
    a, b = h_xyz(h), h_xyz(off)
    others = [k for k in a if ' HZ' not in k]
    assert sorted(others) == sorted(k for k in b if ' HZ' not in k), args
    for k in others:
      d = (matrix.col(a[k]) - matrix.col(b[k])).length()
      assert d < 1.e-3, (args, k, d)
    assert CONVERTED not in description, args
    if args == ['keep_existing_H=True']:
      assert section(description, SUPPLIED_NH2) == ['NZ  LYS A 359 (HZ2, '
        'HZ3)  (NZ  LYS A 359 - CD    CD A 396 2.21 A (geometric))'
        ], description
  # HZ2-NZ-HZ3 105 deg: riding H would place the reloaded NH2 (alg1a) and
  # move it; a fresh approach=optimize keeps both H
  lys_cif('tst_reduce2_metal_amine_lys105.cif', ((
    ' LYS  HZ2  NZ  HZ3  109.000  3.000', ' LYS  HZ2  NZ  HZ3  105.000  3.000'),))
  cif = 'tst_reduce2_metal_amine_lys105.cif'
  h0, description, log = run('tst_reduce2_metal_amine_rr105', pdb_lys_cd,
    [ON, cif])
  assert section(description, CONVERTED)[0].startswith('NZ  LYS A 359: HZ1')
  with open('tst_reduce2_metal_amine_rr105_out.pdb') as f:
    text = f.read()
  h, description, log = run('tst_reduce2_metal_amine_rr105_opt', text,
    [ON, cif, 'approach=optimize'])
  for name in ('HZ2', 'HZ3'):
    d = atom(h, 359, name).distance(atom(h0, 359, name))
    assert d < 1.e-3, (name, d)
  assert section(description, SUPPLIED_NH2) == ['NZ  LYS A 359 (HZ2, HZ3)  '
    '(NZ  LYS A 359 - CD    CD A 396 2.21 A (geometric))'], description

def test_supplied_nh2_any_metal():
  # a supplied NH2 is not completed to NH3 whatever the metal site: Cd one
  # cell away (contact through symmetry), Cd at occupancy 0.5
  nh2 = shared_text('lys_on')
  h0 = iotbx.pdb.input(source_info=None, lines=nh2).construct_hierarchy()
  cd = [l for l in nh2.splitlines() if l.startswith('HETATM') and
        l[17:20] == ' CD'][0]
  x = float(cd[30:38])
  for tag, new in (('sym', cd[:30] + '%8.3f' % (x + 35.557) + cd[38:]),
                   ('occ', cd[:54] + '  0.50' + cd[60:])):
    h, description, log = run('tst_reduce2_metal_amine_s_' + tag,
      nh2.replace(cd, new), [ON, 'keep_existing_H=True'])
    assert h_on(h, 359, 'NZ') == [('HZ2', ''), ('HZ3', '')], (tag,
      h_on(h, 359, 'NZ'))
    for name in ('HZ2', 'HZ3'):
      assert atom(h, 359, name).distance(atom(h0, 359, name)) < 1.e-3, name
    assert section(description, SUPPLIED_NH2)[0].startswith(
      'NZ  LYS A 359 (HZ2, HZ3)'), (tag, description)
    assert CONVERTED not in description, tag

def test_one_supplied_h():
  # one supplied H stays; the second is placed around it, or (supplied H
  # towards the metal) the group is completed as without the option
  nh2 = shared_text('lys_on')
  one = '\n'.join(l for l in nh2.splitlines() if l[12:16] != ' HZ3') + '\n'
  h0 = iotbx.pdb.input(source_info=None, lines=one).construct_hierarchy()
  h, description, log = run('tst_reduce2_metal_amine_one', one,
    [ON, 'keep_existing_H=True'])
  assert section(description, CONVERTED) == ['NZ  LYS A 359: HZ1 not placed; '
    'HZ2 kept as supplied, HZ3 placed  (NZ  LYS A 359 - CD    CD A 396 2.21 A '
    '(geometric))'], description
  assert atom(h, 359, 'HZ2').distance(atom(h0, 359, 'HZ2')) < 1.e-3
  check_nh2(h, 359, 'NZ', 'CE', 396, 'CD', ['HZ2', 'HZ3'])
  # the supplied HZ2 points at the Cd
  n, cd = [matrix.col(atom(h0, r, x).xyz) for r, x in ((359, 'NZ'),
    (396, 'CD'))]
  bad = n + (cd - n).normalize() * 0.89
  hz2 = [l for l in one.splitlines() if l[12:16] == ' HZ2'][0]
  bad_one = one.replace(hz2, hz2[:30] + '%8.3f%8.3f%8.3f' % bad.elems +
    hz2[54:])
  off, d_off, l_off = run('tst_reduce2_metal_amine_bad_off', bad_one,
    [OFF, 'keep_existing_H=True'])
  h, description, log = run('tst_reduce2_metal_amine_bad', bad_one,
    [ON, 'keep_existing_H=True'])
  lines = section(description, UNRESOLVED)
  assert len(lines) == 1 and lines[0].startswith('NZ  LYS A 359: geometry not '
    'placeable') and lines[0].endswith('HZ2 retained and re-idealized '
    '(completed as with the option off)'), lines
  assert CONVERTED not in description
  assert h_on(h, 359, 'NZ') == h_on(off, 359, 'NZ') == [('HZ1', ''),
    ('HZ2', ''), ('HZ3', '')]
  # as without the option, the NH3 rebuild moves the supplied HZ2 too
  for name in ('HZ1', 'HZ2', 'HZ3'):
    assert atom(h, 359, name).distance(atom(off, 359, name)) < 1.e-3, name

def test_names():
  # new H: the first kept name (HZ2) on the side of (CE-NZ) x (Cd-NZ);
  # supplied H keep their names with their coordinates
  h, description, log = shared('lys_on')
  n, ce, cd = [matrix.col(atom(h, r, x).xyz)
    for r, x in ((359, 'NZ'), (359, 'CE'), (396, 'CD'))]
  side = (ce - n).cross(cd - n)
  assert (matrix.col(atom(h, 359, 'HZ2').xyz) - n).dot(side) > 0
  assert (matrix.col(atom(h, 359, 'HZ3').xyz) - n).dot(side) < 0
  lines = shared_text('lys_on').splitlines()
  i2 = [i for i, l in enumerate(lines) if l[12:16] == ' HZ2'][0]
  i3 = [i for i, l in enumerate(lines) if l[12:16] == ' HZ3'][0]
  l2, l3 = lines[i2], lines[i3]
  lines[i2] = l2[:30] + l3[30:54] + l2[54:]
  lines[i3] = l3[:30] + l2[30:54] + l3[54:]
  swapped = '\n'.join(lines) + '\n'
  for args in (['keep_existing_H=True'], ['approach=optimize']):
    hs, description, log = run('tst_reduce2_metal_amine_nm_s', swapped,
      [ON] + args)
    for a, b in (('HZ2', 'HZ3'), ('HZ3', 'HZ2')):
      d = atom(hs, 359, a).distance(atom(h, 359, b))
      assert d < 1.e-3, (args, a, d)

def test_nh3_across_conformers():
  # HZ1 blank, HZ2/HZ3 in A and in B (B moved 0.1 A): three H in each
  # conformer, an NH3: handled as without the option
  nh3 = shared_text('lys_off')
  lines = []
  for l in nh3.splitlines():
    if l[12:16] in (' HZ2', ' HZ3'):
      lines.append(l[:16] + 'A' + l[17:54] + '  0.50' + l[60:])
      x = float(l[30:38]) + 0.1
      lines.append(l[:16] + 'B' + l[17:30] + '%8.3f' % x + l[38:54] + '  0.50'
        + l[60:])
    else:
      lines.append(l)
  text = '\n'.join(lines) + '\n'
  for args in (['approach=optimize'], ['keep_existing_H=True']):
    off, d_off, l_off = run('tst_reduce2_metal_amine_cf_off', text,
      [OFF] + args)
    on, d_on, l_on = run('tst_reduce2_metal_amine_cf_on', text, [ON] + args)
    assert SUPPLIED_NH2 not in d_on and CONVERTED not in d_on, (args, d_on)
    a = dict((x.id_str(), x.xyz) for x in off.atoms())
    b = dict((x.id_str(), x.xyz) for x in on.atoms())
    assert sorted(a) == sorted(b), args
    for k in a:
      d = (matrix.col(a[k]) - matrix.col(b[k])).length()
      assert d < 1.e-3, (args, k, d)
  # documented limits
  def split(alts):
    # alts: {name: altlocs}; '' keeps the blank H; B moved 0.1 A
    out = []
    for l in nh3.splitlines():
      name = l[12:16].strip()
      if l[17:20] == 'LYS' and name in alts:
        for alt in alts[name]:
          if alt == '':
            out.append(l)
            continue
          x = float(l[30:38]) + (0.1 if alt == 'B' else 0)
          out.append(l[:16] + alt + l[17:30] + '%8.3f' % x + l[38:54] + '  0.50'
            + l[60:])
      elif not (l[17:20] == 'LYS' and name == 'HZ3' and 'HZ3' not in alts):
        out.append(l)
    return '\n'.join(out) + '\n'
  def same(h1, h2, names, tag):
    a = dict((x.id_str(), x.xyz) for x in h1.atoms())
    b = dict((x.id_str(), x.xyz) for x in h2.atoms())
    for k in a:
      if k[5:9].strip() not in names: continue
      d = (matrix.col(a[k]) - matrix.col(b[k])).length()
      assert d < 1.e-3, (tag, k, d)
  mixed = ('NZ  LYS A 359: H in blank and alternate conformers: completed and '
    're-idealized as with the option off')
  # 1. blank HZ1, HZ2 in A and B: an NH2 in each conformer
  text = split({'HZ1': [''], 'HZ2': ['A', 'B']})
  h_in = iotbx.pdb.input(source_info=None, lines=text).construct_hierarchy()
  h, d, l = run('tst_reduce2_metal_amine_cf1_opt', text,
    [ON, 'approach=optimize'])
  same(h, h_in, ('HZ1', 'HZ2'), 'cf1 optimize')
  assert section(d, SUPPLIED_NH2) == [
    'NZ  LYS A 359 (HZ1, HZ2 altloc A)  (NZ  LYS A 359 - CD    CD A 396 2.21 A '
    '(geometric))', 'NZ  LYS A 359 (HZ1, HZ2 altloc B)  (NZ  LYS A 359 - CD    '
    'CD A 396 2.21 A (geometric))'], d
  off, d_off, l_off = run('tst_reduce2_metal_amine_cf1_off', text,
    [OFF, 'keep_existing_H=True'])
  h, d, l = run('tst_reduce2_metal_amine_cf1_keep', text,
    [ON, 'keep_existing_H=True'])
  assert sorted((x.name.strip(), x.parent().altloc) for x in h.atoms()
    if x.name.strip().startswith('HZ')) == [('HZ1', ''), ('HZ2', 'A'),
    ('HZ2', 'B'), ('HZ3', 'A'), ('HZ3', 'B')]
  same(h, off, ('HZ1', 'HZ2', 'HZ3'), 'cf1 keep')
  assert section(d, UNRESOLVED) == [mixed, mixed], d
  assert SUPPLIED_NH2 not in d
  # 2. blank HZ1, HZ2/HZ3 in A, HZ2 only in B: B an NH2, A an NH3. The
  # shared HZ1 belongs to A's NH3 and may be moved by its rotator (a limit);
  # B's own HZ2 stays
  text = split({'HZ1': [''], 'HZ2': ['A', 'B'], 'HZ3': ['A']})
  h_in = iotbx.pdb.input(source_info=None, lines=text).construct_hierarchy()
  h, d, l = run('tst_reduce2_metal_amine_cf2_opt', text,
    [ON, 'approach=optimize'])
  hz = lambda hh: sorted((x.name.strip(), x.parent().altloc)
    for x in hh.atoms() if x.element.strip() == 'H')
  assert hz(h) == hz(h_in)
  get = lambda hh, name, alt: [x for x in hh.atoms()
    if x.name.strip() == name and x.parent().altloc == alt][0]
  assert get(h, 'HZ2', 'B').distance(get(h_in, 'HZ2', 'B')) < 1.e-3
  assert 'Added MoverNH3Rotator 1 to chain A LYS 359' in d, d
  assert get(h, 'HZ1', '').distance(get(h_in, 'HZ1', '')) > 0.5
  assert SUPPLIED_NH2 + ': kept as supplied' not in d
  assert section(d, SUPPLIED_NH2 + ', partly preserved') == [
    "NZ  LYS A 359 conformer B: supplied NH2 partly preserved (HZ2 kept; "
    "shared HZ1 belongs to conformer A's NH3 and may move)  (NZ  LYS A 359 - "
    "CD    CD A 396 2.21 A (geometric))"], d

if __name__ == '__main__':
  test_chain_start_rerun()
  test_chain_start_supplied_nh2_all_chain_starts()
  test_chain_start_supplied_h()
  test_chain_start_conformer_supplied()
  test_supplied()
  test_rerun_stable()
  test_supplied_nh2_any_metal()
  test_one_supplied_h()
  test_names()
  test_nh3_across_conformers()
  print('OK')
