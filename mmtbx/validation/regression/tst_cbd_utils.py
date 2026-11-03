from __future__ import absolute_import, division, print_function
from math import pi, cos, sin, sqrt, atan2, radians, degrees
from libtbx.test_utils import approx_equal
from mmtbx.validation import cbd_utils
from mmtbx.validation.cbd_utils import radial_deviation

def as_xy(v):
  return (v.r*cos(v.t), v.r*sin(v.t))

def exercise_angle_wrapping():
  # 0.3 at 175 deg minus 0.1 at 90 deg lands at ~193.9 deg if not wrapped
  rc = radial_deviation(0.3, 175) - radial_deviation(0.1, 90)
  assert -pi < rc.t <= pi, rc.t/pi*180
  assert approx_equal(rc.t/pi*180, 193.88072251867072-360, eps=1e-6)
  # wrapping must not change the vector itself
  a = radial_deviation(0.3, 175)
  b = radial_deviation(0.1, 90)
  assert approx_equal(as_xy(rc), (as_xy(a)[0]-as_xy(b)[0],
                                  as_xy(a)[1]-as_xy(b)[1]), eps=1e-9)
  # sweep: results always in [-180, 180] (+/-180 equivalent, fp noise)
  for t1 in range(-180, 181, 15):
    for t2 in range(-180, 181, 15):
      for op in ['add', 'sub']:
        x = radial_deviation(0.3, t1)
        y = radial_deviation(0.2, t2)
        z = x+y if op=='add' else x-y
        assert -pi-1e-12 <= z.t <= pi+1e-12, (t1, t2, op, z.t/pi*180)

def exercise_input_in_cartesian():
  # (r, t) are read as (x, y); all four quadrants and the axes
  for x, y in [(3, 4), (-3, 4), (-3, -4), (3, -4), (0, 2), (-2, 0), (2, 0)]:
    v = radial_deviation(x, y, input_in_cartesian=True)
    assert approx_equal(v.r, (x*x+y*y)**0.5, eps=1e-9), (x, y, v.r)
    assert approx_equal(as_xy(v), (x, y), eps=1e-9), (x, y, as_xy(v))

def exercise_tables():
  assert sorted(cbd_utils.tables) == ['PRO', 'VAL_THR_ILE', 'VAL_THR_ILE_pre_PRO',
                                      'others', 'others_pre_PRO']
  for name, table in cbd_utils.tables.items():
    assert table, name
    for (phi, psi), (r, t) in table.items():
      assert 180 not in (phi, psi), (name, phi, psi)
      assert -pi-1e-6 <= t <= pi+1e-6, (name, phi, psi, t)
  # the pre-PRO helix offset is much larger than for other residues
  assert approx_equal(cbd_utils.tables['others_pre_PRO'][(-60, -50)][0], 0.110,
                      eps=0.003)
  assert approx_equal(cbd_utils.tables['VAL_THR_ILE_pre_PRO'][(-60, -50)][0],
                      0.155, eps=0.003)
  assert cbd_utils.tables['others'][(-60, -50)][0] < 0.01

def exercise_table_selection():
  T = cbd_utils.tables
  for resname, next_resname, name in [
      ('PRO', 'PRO', 'PRO'),
      ('PRO', 'ALA', 'PRO'),
      ('VAL', 'PRO', 'VAL_THR_ILE_pre_PRO'),
      ('THR', None, 'VAL_THR_ILE'),
      ('ALA', 'PRO', 'others_pre_PRO'),
      ('ALA', 'GLY', 'others'),
      ('ALA', 'HYP', 'others'),
      ]:
    assert cbd_utils.get_correction_table(resname, next_resname) is T[name], \
      (resname, next_resname)

def exercise_wrap_lookup():
  # psi 176 rounds to 180, which is the -180 key
  table = cbd_utils.tables['others']
  key = [k for k in sorted(table) if k[1] == -180][0]
  v = cbd_utils.get_correction_vector('ALA', 'GLY', (key[0], 176.))
  r, t = table[key]
  assert approx_equal(v, (r*cos(t), r*sin(t)))
  # D-amino acids are not corrected
  assert cbd_utils.get_correction_vector('DAL', 'GLY', (key[0], 176.)) is None

def synthetic(n, scale):
  # v = scale*c for every residue, c at varied angles
  devs, dihedrals, corrections = [], [], []
  for i in range(n):
    t = radians(i*37 % 360)
    c = (0.05*cos(t), 0.05*sin(t))
    v = (scale*c[0], scale*c[1])
    devs.append(sqrt(v[0]**2 + v[1]**2))
    dihedrals.append(degrees(atan2(v[1], v[0])))
    corrections.append(c)
  return devs, dihedrals, corrections

def exercise_apply_correction():
  d, t = cbd_utils.apply_correction(0.3, 45., (0.1, 0.), 0.5)
  x, y = 0.3*cos(radians(45)) - 0.05, 0.3*sin(radians(45))
  assert approx_equal((d, t), (sqrt(x*x + y*y), degrees(atan2(y, x))))
  assert approx_equal(cbd_utils.apply_correction(0.3, 45., (0.1, 0.), 0.),
                      (0.3, 45.))
  assert -180 < cbd_utils.apply_correction(0.1, 179., (0., 0.05), 1.)[1] <= 180

def exercise_choose_strength():
  s = cbd_utils.choose_strength(*synthetic(3000, 0.6))
  assert approx_equal(s.fitted, 0.6) and abs(s.applied - 0.6) < 0.005
  assert s.n == 3000 and not s.user_set and not s.check_failed
  # a small model is pulled toward 0.5
  s = cbd_utils.choose_strength(*synthetic(30, 0.6))
  assert approx_equal(s.pulled, (18 + 50)/130.) and 0.5 < s.applied < 0.6
  s = cbd_utils.choose_strength(*synthetic(3000, 1.4))
  assert s.fitted > 1 and s.applied == 1.
  s = cbd_utils.choose_strength(*synthetic(3000, -0.4))
  assert approx_equal(s.fitted, -0.4) and s.applied == 0.
  s = cbd_utils.choose_strength([], [], [])
  assert s.n == 0 and s.fitted is None and s.applied == 0.

def exercise_never_worse():
  # one large residue drives the fit; 200 residues with v = 0 get worse
  devs, dihedrals, corrections = [1.5], [0.], [(1.5, 0.)]
  devs += [0.]*200
  dihedrals += [0.]*200
  corrections += [(0.1, 0.)]*200
  s = cbd_utils.choose_strength(devs, dihedrals, corrections)
  assert s.pulled > 0.5 and s.check_failed and s.applied == 0.
  assert s.mean_corrected > s.mean_uncorrected
  s = cbd_utils.choose_strength(devs, dihedrals, corrections, user_strength=1.)
  assert s.user_set and not s.check_failed and s.applied == 1.
  assert approx_equal(s.fitted, 2.25/4.25)

def exercise_gross_deviation_not_fitted():
  # one residue 2.5 A off, against its table vector, must not drive the fit
  devs, dihedrals, corrections = synthetic(3000, 0.6)
  devs.append(2.5)
  dihedrals.append(180.)
  corrections.append((0.05, 0.))
  s = cbd_utils.choose_strength(devs, dihedrals, corrections)
  assert abs(s.fitted - 0.6) < 0.01, s.fitted
  assert s.n == 3000, s.n

def run():
  exercise_angle_wrapping()
  exercise_input_in_cartesian()
  exercise_tables()
  exercise_table_selection()
  exercise_wrap_lookup()
  exercise_apply_correction()
  exercise_choose_strength()
  exercise_never_worse()
  exercise_gross_deviation_not_fitted()
  print('OK')

if __name__ == '__main__':
  run()
