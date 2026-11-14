"""reduce2: metal-bound neutral NH2 (A17 step B2). A neutral NH2 that its
dictionary puts in riding H's pyramidal scope and whose N binds a metal keeps
both H, placed with the lone pair at the metal (metal_bound_amines_as_nh2):
EDN (ethylenediamine, CCD ideal heavy atoms) with Cu 2.00 A from NAA, C-N-Cu
109.5 deg, at Cu-NAA-CAB-CAC 60, 180 and -60 deg. Before (riding's dictionary
torsion, no metal): an H 1.11-1.13 A from Cu at 180 and -60 deg. The H keep
the dictionary angles, so a fresh riding idealization of the output and
approach=optimize leave them where they are. With the option off they are
placed as without a metal and listed. Two Cu positions (alternate
conformers) are unresolved: the NH2 is turned away from both. A Cu beyond
the coordination window (2.45 A) orients the NH2 too, and is avoided when
another Cu is in the window. One supplied H that the construction cannot
keep is placed with the other, away from Cu.
"""
from __future__ import absolute_import, division, print_function
import math
import os
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out
from scitbx import matrix
from scitbx.array_family import flex
from six.moves import cStringIO as StringIO
from iotbx.cli_parser import run_program
from mmtbx.programs import reduce2
from mmtbx.hydrogens import reduce_hydrogen

EDN = [('NAA', (16.861, 14.833, 15.000)), ('CAB', (15.562, 15.519, 15.000)),
       ('CAC', (14.438, 14.481, 15.000)), ('NAD', (13.139, 15.167, 15.000))]

def place_atom(a, b, c, bond, angle, torsion):
  """the site d with |cd| = bond, angle bcd, torsion abcd"""
  a, b, c = [matrix.col(v) for v in (a, b, c)]
  bc = (c - b).normalize(); n = (b - a).cross(bc).normalize(); m = n.cross(bc)
  th, ph = math.radians(angle), math.radians(torsion)
  d2 = matrix.col((-bond * math.cos(th), bond * math.sin(th) * math.cos(ph),
    bond * math.sin(th) * math.sin(ph)))
  return c + matrix.sqr(tuple(bc) + tuple(m) + tuple(n)).transpose() * d2

def edn_cu(torsions, occ=None, distance=2.00):
  lines = ["CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1"]
  for i, (n, x) in enumerate(EDN):
    lines.append("HETATM%5d  %-3s EDN A   1    %8.3f%8.3f%8.3f  1.00 20.00"
      "           %s" % (i + 1, n, x[0], x[1], x[2], n[0]))
  for k, t in enumerate(torsions):
    cu = place_atom(EDN[2][1], EDN[1][1], EDN[0][1], distance, 109.5, t)
    alt = ' ' if len(torsions) == 1 else 'AB'[k]
    o = 1.0 if len(torsions) == 1 else 0.5
    if occ is not None: o = occ
    lines.append("HETATM%5d CU  %s CU A   2    %8.3f%8.3f%8.3f%6.2f 20.00"
      "          CU" % (5 + k, alt, cu[0], cu[1], cu[2], o))
  return '\n'.join(lines) + '\nEND\n'

CRYST1 = "CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1\n"

# HOA (hydroxylamine, X = O) with Cu 2.00 A from N, O-N-Cu 109.5 deg
def hoa_cu():
  n, o = matrix.col((20.740, 20.031, 19.924)), matrix.col((19.284, 19.933,
    20.030))
  ref = matrix.col((20.740, 21.5, 19.924))
  cu = place_atom(ref, o, n, 2.00, 109.5, 180)
  return '\n'.join([CRYST1.strip(),
    "HETATM    1  N   HOA A   1    %8.3f%8.3f%8.3f  1.00 20.00           N"
    % tuple(n),
    "HETATM    2  O   HOA A   1    %8.3f%8.3f%8.3f  1.00 20.00           O"
    % tuple(o),
    "HETATM    3 CU    CU A   2    %8.3f%8.3f%8.3f  1.00 20.00          CU"
    % tuple(cu)]) + '\nEND\n'

def run(prefix, pdb_str, args=()):
  fn, out = prefix + '_in.pdb', prefix + '_out.pdb'
  with open(fn, 'w') as f:
    f.write(pdb_str)
  if os.path.exists(out): os.remove(out)
  run_program(program_class=reduce2.Program, logger=StringIO(), args=[fn,
    'output.filename=' + out, 'output.description_file_name=%s.txt' % prefix,
    'output.overwrite=True'] + list(args))
  with open(prefix + '.txt') as f:
    description = f.read()
  return out, description

def sites(h):
  return dict((a.name.strip(), matrix.col(a.xyz)) for a in h.atoms()
    if a.name.strip() != 'CU')

def cu_sites(h):
  return [matrix.col(a.xyz) for a in h.atoms() if a.name.strip() == 'CU']

def angle(a, b, c):
  return math.degrees((a - b).angle(c - b))

def check(s, cus):
  """both H on NAA, dictionary H-N-H and C-N-H; H...Cu and Cu-N-H"""
  h1, h2 = s['HAA1'], s['HAA2']
  assert abs(angle(h1, s['NAA'], h2) - 105.5) < 0.2, angle(h1, s['NAA'], h2)
  for h in (h1, h2):
    assert abs(angle(s['CAB'], s['NAA'], h) - 109.4) < 0.3, angle(s['CAB'],
      s['NAA'], h)
  return min((h - cu).length() for h in (h1, h2) for cu in cus), [
    angle(cu, s['NAA'], h) for h in (h1, h2) for cu in cus]

def test_three_torsions():
  # the lone pair at Cu: H...Cu 2.4-2.6 A, Cu-N-H 105-115 deg, whatever the
  # torsion; the other N (NAD, no metal) unchanged pyramidal
  for t in (60, 180, -60):
    out, description = run('tst_reduce2_metal_neutral_nh2_%d' % (t % 360),
      edn_cu([t]))
    h = iotbx.pdb.input(out).construct_hierarchy()
    d, cnh = check(sites(h), cu_sites(h))
    assert 2.4 < d < 2.6, (t, d)
    for v in cnh:
      assert 105 < v < 115, (t, cnh)
    assert reduce_hydrogen.METAL_NEUTRAL_NH2_ORIENTED in description, t
    assert 'NAA EDN A   1 (HAA1, HAA2)' in description, description

def test_reload_and_optimize():
  # a fresh riding idealization (current torsion) and approach=optimize on
  # the output leave the H where placement put them
  out, _ = run('tst_reduce2_metal_neutral_nh2_reload', edn_cu([180]))
  before = sites(iotbx.pdb.input(out).construct_hierarchy())
  m = mmtbx.model.manager(model_input=iotbx.pdb.input(out), log=null_out())
  reduce_hydrogen.add_ccd_restraint_objects(m)
  m.process(make_restraints=True)
  m.setup_riding_h_manager(idealize=True)
  after = sites(m.get_hierarchy())
  for name in ('HAA1', 'HAA2', 'HAD1', 'HAD2'):
    assert (after[name] - before[name]).length() < 0.005, (name,
      (after[name] - before[name]).length())
  with open(out) as f:
    placed = f.read()
  out2, description = run('tst_reduce2_metal_neutral_nh2_optimize', placed,
    ['approach=optimize'])
  again = sites(iotbx.pdb.input(out2).construct_hierarchy())
  for name in ('HAA1', 'HAA2'):
    assert (again[name] - before[name]).length() < 0.005, (name,
      (again[name] - before[name]).length())
  # supplied H off the dictionary angles (HAA2 turned 10 deg about N-Cu):
  # approach=optimize keeps both where they are (not re-idealized by riding)
  h = iotbx.pdb.input(out).construct_hierarchy()
  s = sites(h)
  cu = cu_sites(h)[0]
  r = (cu - s['NAA']).axis_and_angle_as_r3_rotation_matrix(10, deg=True)
  for a in h.atoms():
    if a.name.strip() == 'HAA2':
      a.xyz = tuple(r * (matrix.col(a.xyz) - s['NAA']) + s['NAA'])
  moved = sites(h)['HAA2']
  out3, description = run('tst_reduce2_metal_neutral_nh2_supplied',
    CRYST1 + h.as_pdb_string(), ['approach=optimize'])
  kept = sites(iotbx.pdb.input(out3).construct_hierarchy())['HAA2']
  assert (kept - moved).length() < 0.005, (kept - moved).length()

def test_option_off():
  # metal_bound_amines_as_nh2=False: placed as without a metal (an H 1.11 A
  # from Cu at 180 deg), and listed
  out, description = run('tst_reduce2_metal_neutral_nh2_off', edn_cu([180]),
    ['metal_bound_amines_as_nh2=False'])
  h = iotbx.pdb.input(out).construct_hierarchy()
  d, cnh = check(sites(h), cu_sites(h))
  assert d < 1.3, d
  assert ('Metal-bound neutral NH2, not oriented '
    '(metal_bound_amines_as_nh2=False)') in description, description

def test_unresolved_two_positions():
  # Cu in A (60 deg) and B (-60 deg), each 0.5: no unique metal position;
  # listed, turned to the largest smallest H...Cu over both positions. They
  # are 120 deg apart about C-N, so no NH2 avoids both: 1.71 A, the most a
  # pyramidal pair gets (master: an H 1.11 A from Cu B)
  out, description = run('tst_reduce2_metal_neutral_nh2_two',
    edn_cu([60, -60]))
  h = iotbx.pdb.input(out).construct_hierarchy()
  d, cnh = check(sites(h), cu_sites(h))
  assert d > 1.65, d
  assert 'Metal-bound neutral NH2, unresolved' in description, description
  assert 'turned away from the metal positions' in description, description

def test_near_metal():
  # Cu 2.45 A from NAA, beyond the coordination window (radical-SAM SAM/SAH
  # bind Fe at 2.2-2.6 A): oriented the same way, listed as outside the
  # window (before: an H about 1.4 A from Cu at 180 deg)
  out, description = run('tst_reduce2_metal_neutral_nh2_near',
    edn_cu([180], distance=2.45))
  h = iotbx.pdb.input(out).construct_hierarchy()
  d, cnh = check(sites(h), cu_sites(h))
  assert 2.8 < d < 3.1, d
  assert 'metal outside the coordination window' in description, description
  # a Cu at occupancy 0 is not there (as in the window): not oriented
  out, description = run('tst_reduce2_metal_neutral_nh2_near0',
    edn_cu([180], occ=0.0, distance=2.45))
  assert 'Metal-bound neutral NH2' not in description, description

def test_hydroxylamine():
  # X = O (HOA, no anchor for riding H): oriented with the lone pair at Cu
  out, description = run('tst_reduce2_metal_neutral_nh2_hoa', hoa_cu())
  h = iotbx.pdb.input(out).construct_hierarchy()
  s = sites(h)
  cu = cu_sites(h)[0]
  d = min((s[x] - cu).length() for x in ('HN1', 'HN2'))
  assert 2.3 < d < 2.7, d
  assert reduce_hydrogen.METAL_NEUTRAL_NH2_ORIENTED in description, \
    description

def test_one_supplied_h():
  # keep_existing_H with HAA1 supplied anti to CAC about C-N (on the Cu at
  # 180 deg) and HAA2 missing: the construction around the supplied H is
  # refused, so both H are turned away from Cu and listed (before: the
  # supplied H at 1.11 A from Cu stayed, as on master)
  from mmtbx.hydrogens import reduce_hydrogen as rh
  def place(pdb_str, keep):
    m = mmtbx.model.manager(model_input=iotbx.pdb.input(source_info=None,
      lines=pdb_str.split('\n')), log=null_out())
    o = rh.place_hydrogens(model=m, stop_for_unknowns=False,
      metal_bound_amines_as_nh2=True, keep_existing_H=keep)
    o.run()
    return o
  h = place(edn_cu([180]), False).get_model().get_hierarchy()
  s = sites(h)
  hx = place_atom(s['CAC'], s['CAB'], s['NAA'], 0.89, 109.5, 180)
  for a in list(h.atoms()):
    if a.name.strip() == 'HAA2': a.parent().remove_atom(a)
    elif a.name.strip() == 'HAA1': a.xyz = tuple(hx)
  o = place(CRYST1 + h.as_pdb_string(), True)
  out = o.get_model().get_hierarchy()
  s, cu = sites(out), cu_sites(out)[0]
  d = min((s[x] - cu).length() for x in ('HAA1', 'HAA2'))
  assert d > 2.3, d
  assert [t for t, c in o.metal_amine_neutral_unresolved
          if 'supplied HAA1 re-placed' in t], o.metal_amine_neutral_unresolved

def test_window_and_near_metal():
  # two Cu present at once: 2.00 A from NAA at +60 deg (in the window) and
  # 2.45 A at -60 deg (outside it). Built around the first, an H lies 1.56 A
  # from the second; turned away from both the closest H...Cu is 1.92 A. So
  # the site is unresolved and turned away from both, not oriented around
  # one metal
  def two_cu(distance, angle, torsion):
    lines = [l for l in edn_cu([60]).splitlines() if not l.startswith('END')]
    cu = place_atom(EDN[2][1], EDN[1][1], EDN[0][1], distance, angle, torsion)
    lines.append("HETATM    6 CU    CU A   3    %8.3f%8.3f%8.3f  1.00 20.00"
      "          CU" % tuple(cu))
    m = mmtbx.model.manager(model_input=iotbx.pdb.input(source_info=None,
      lines=lines), log=null_out())
    o = reduce_hydrogen.place_hydrogens(model=m, stop_for_unknowns=False,
      metal_bound_amines_as_nh2=True)
    o.run()
    return o
  o = two_cu(2.45, 109.5, -60)
  h = o.get_model().get_hierarchy()
  d, cnh = check(sites(h), cu_sites(h))
  assert len(cu_sites(h)) == 2
  assert d > 1.90, d
  assert not [t for t, c in o.metal_amine_neutral if 'NAA' in t], \
    o.metal_amine_neutral
  hits = [(t, c) for t, c in o.metal_amine_neutral_unresolved if 'NAA' in t]
  assert len(hits) == 1 and 'turned away' in hits[0][0], \
    o.metal_amine_neutral_unresolved
  # both Cu listed with the site
  assert 'CU A   2 2.00 A' in hits[0][1] and 'CU A   3 2.45 A' in hits[0][1], \
    hits
  # the second Cu 2.50 A from NAA between CAB and the first Cu (C-N-Cu 55
  # deg), clear of the H built around the first: accepted, both Cu listed
  o = two_cu(2.50, 55, 60)
  h = o.get_model().get_hierarchy()
  d, cnh = check(sites(h), cu_sites(h)[:1])
  assert 2.4 < d < 2.6, d
  hits = [(t, c) for t, c in o.metal_amine_neutral if 'NAA' in t]
  assert len(hits) == 1 and 'CU A   3 2.50 A' in hits[0][1], \
    (o.metal_amine_neutral, o.metal_amine_neutral_unresolved)

if __name__ == '__main__':
  test_three_torsions()
  test_reload_and_optimize()
  test_option_off()
  test_unresolved_two_positions()
  test_near_metal()
  test_hydroxylamine()
  test_one_supplied_h()
  test_window_and_near_metal()
  print('OK')
