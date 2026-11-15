"""Riding H torsion targets from dihedral proxies the registry stored
reordered (A62). The proxy registry stores sort_i_seqs() order (i1 < i4,
i2 < i3) and negates the ideal per swap, so the atom next to the H may be
its parent or the parent's neighbour. Riding's convention is
dihedral(third, a1, a0, H). Before: a Var target was taken in the proxy's
order (wrong sign when the atom next to the H is not its parent) and a CONST
ideal was negated whenever the H was the last atom (wrong sign when that
atom is its parent). Fixture: ethanol with torsions H11-C1-C2-O1 (methyl,
propeller) and HO1-O1-C2-C1 (hydroxyl, alg1b); four atom orders give the
four stored forms."""
from __future__ import absolute_import, division, print_function
import math
import iotbx.pdb
import iotbx.cif
import mmtbx.model
from cctbx import geometry_restraints
from scitbx import matrix
from scitbx.math import dihedral_angle
from libtbx.utils import null_out

HEAVY = ['C1', 'C2', 'O1']
HYD = ['H11', 'H12', 'H13', 'H21', 'H22', 'HO1']
BONDS = [('C1', 'C2', 1.52), ('C2', 'O1', 1.43), ('C1', 'H11', 0.97),
  ('C1', 'H12', 0.97), ('C1', 'H13', 0.97), ('C2', 'H21', 0.97),
  ('C2', 'H22', 0.97), ('O1', 'HO1', 0.85)]
ANGLES = [('C1', 'C2', 'O1'), ('C2', 'C1', 'H11'), ('C2', 'C1', 'H12'),
  ('C2', 'C1', 'H13'), ('H11', 'C1', 'H12'), ('H11', 'C1', 'H13'),
  ('H12', 'C1', 'H13'), ('C1', 'C2', 'H21'), ('C1', 'C2', 'H22'),
  ('O1', 'C2', 'H21'), ('O1', 'C2', 'H22'), ('H21', 'C2', 'H22'),
  ('C2', 'O1', 'HO1')]
# both in riding's order read backwards: dihedral(H, a0, a1, third)
TORSIONS = {'H11': ('H11', 'C1', 'C2', 'O1'), 'HO1': ('HO1', 'O1', 'C2', 'C1')}
# stored forms (sort_i_seqs) of the two torsions in each atom order
ORDERS = {
  'A':  HEAVY + HYD,                  # (O1,C1,C2,H11)  (C1,C2,O1,HO1)
  'B':  HEAVY[::-1] + HYD,            # (O1,C2,C1,H11)  (C1,O1,C2,HO1)
  'HA': HYD + HEAVY,                  # (H11,C1,C2,O1)  (HO1,C2,O1,C1)
  'HB': HYD + HEAVY[::-1]}            # (H11,C2,C1,O1)  (HO1,O1,C2,C1)

def nerf(a, b, c, d, ang, tor):
  a, b, c = [matrix.col(v) for v in (a, b, c)]
  bc = (c - b).normalize(); n = (b - a).cross(bc).normalize(); m = n.cross(bc)
  A, T = math.radians(ang), math.radians(tor)
  return tuple(c + (-d * math.cos(A)) * bc + (d * math.sin(A) * math.cos(T)) * m
    + (d * math.sin(A) * math.sin(T)) * n)

def sites(seed):
  """seed: torsion of H11 and HO1, as in TORSIONS"""
  s = {'C1': (10., 10., 10.), 'C2': (11.52, 10., 10.)}
  s['O1'] = nerf((10., 11., 10.), s['C1'], s['C2'], 1.43, 109.5, 180.)
  for n, t in (('H11', seed), ('H12', seed + 120), ('H13', seed - 120)):
    s[n] = nerf(s['O1'], s['C2'], s['C1'], 0.97, 109.5, t)
  s['HO1'] = nerf(s['C1'], s['C2'], s['O1'], 0.85, 109.5, seed)
  for n, t in (('H21', 120.), ('H22', -120.)):
    s[n] = nerf(s['O1'], s['C1'], s['C2'], 0.97, 109.5, t)
  return s

def pdb(order, seed):
  s = sites(seed)
  lines = ['CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1']
  for k, n in enumerate(ORDERS[order]):
    e = 'H' if n.startswith('H') else n[0]
    lines.append('HETATM%5d %-4s LG8 A   1    %8.3f%8.3f%8.3f  1.00 20.00'
      '          %2s' % ((k + 1, ' ' + n) + s[n] + (e,)))
  return '\n'.join(lines + ['END']) + '\n'

def cif(v, period, const):
  L = ['data_comp_list', 'loop_', '_chem_comp.id',
    '_chem_comp.three_letter_code', '_chem_comp.name', '_chem_comp.group',
    '_chem_comp.number_atoms_all', '_chem_comp.number_atoms_nh',
    '_chem_comp.desc_level', "LG8 LG8 'ethanol' ligand 9 3 .", '',
    'data_comp_LG8', 'loop_', '_chem_comp_atom.comp_id',
    '_chem_comp_atom.atom_id', '_chem_comp_atom.type_symbol',
    '_chem_comp_atom.type_energy']
  te = dict(C1='CH3', C2='CH2', O1='OH1')
  L += ['LG8 %s %s %s' % (n, 'H' if n.startswith('H') else n[0],
    te.get(n, 'H')) for n in HEAVY + HYD]
  L += ['loop_', '_chem_comp_bond.comp_id', '_chem_comp_bond.atom_id_1',
    '_chem_comp_bond.atom_id_2', '_chem_comp_bond.type',
    '_chem_comp_bond.value_dist', '_chem_comp_bond.value_dist_esd']
  L += ['LG8 %s %s single %.3f 0.020' % b for b in BONDS]
  L += ['loop_', '_chem_comp_angle.comp_id', '_chem_comp_angle.atom_id_1',
    '_chem_comp_angle.atom_id_2', '_chem_comp_angle.atom_id_3',
    '_chem_comp_angle.value_angle', '_chem_comp_angle.value_angle_esd']
  L += ['LG8 %s %s %s 109.5 3.0' % a for a in ANGLES]
  L += ['loop_', '_chem_comp_tor.comp_id', '_chem_comp_tor.id',
    '_chem_comp_tor.atom_id_1', '_chem_comp_tor.atom_id_2',
    '_chem_comp_tor.atom_id_3', '_chem_comp_tor.atom_id_4',
    '_chem_comp_tor.value_angle', '_chem_comp_tor.value_angle_esd',
    '_chem_comp_tor.period']
  for k, h in enumerate(('H11', 'HO1')):
    L.append('LG8 %s_%d %s %s %s %s %.2f %s %d' % (('CONST' if const else
      'var', k + 1) + TORSIONS[h] + (v, '0.0' if const else '10.0', period)))
  return '\n'.join(L) + '\n'

def model(order, v, period=3, const=False, seed=40., riding=True, **kw):
  """unsorted, so that the file order sets the i_seqs (sorting puts the H
  last and the heavy atoms by name)"""
  inp = iotbx.pdb.input(source_info=None, lines=pdb(order, seed))
  m = mmtbx.model.manager(model_input=None,
    pdb_hierarchy=inp.construct_hierarchy(sort_atoms=False),
    crystal_symmetry=inp.crystal_symmetry(), log=null_out())
  m.set_restraint_objects([('lg8.cif',
    iotbx.cif.reader(input_string=cif(v, period, const)).model())])
  p = m.get_default_pdb_interpretation_params()
  p.pdb_interpretation.sort_atoms = False
  m.process(make_restraints=True, pdb_interpretation_params=p)
  if riding:
    m.setup_riding_h_manager(**kw)
  return m

def atoms(m):
  return dict((a.name.strip(), a) for a in m.get_hierarchy().atoms())

def torsion(m, h):
  at = atoms(m)
  return math.degrees(dihedral_angle(sites=[at[n].xyz for n in TORSIONS[h]]))

def stored(m):
  g = m.get_restraints_manager().geometry
  names = [a.name.strip() for a in m.get_hierarchy().atoms()]
  return sorted(tuple(names[i] for i in p.i_seqs) for p in
    list(g.dihedral_proxies) + list(g.const_dihedral_proxies or []))

def check_placed(order, v, period, const, tol=0.01):
  m = model(order, v, period, const, use_ideal_dihedral=True)
  par = m.riding_h_manager.h_parameterization
  at = atoms(m)
  assert par[at['H11'].i_seq].htype == 'prop'
  assert par[at['HO1'].i_seq].htype == 'alg1b'
  for h in ('H11', 'HO1'):
    t = torsion(m, h)
    r = geometry_restraints.angle_delta_deg(angle_1=t, angle_2=v,
      periodicity=1 if const else period)
    assert abs(r) < tol, (order, v, period, const, h, t, r)
  return m

def exercise_stored_forms():
  # the four orders give the four stored forms of each torsion
  forms = set()
  for order in ORDERS:
    s = stored(model(order, 64., riding=False))
    forms.update(s)
  assert len(forms) == 8, forms

def placed_twice(m):
  # a second placement (new manager) moves nothing
  before = dict((n, matrix.col(a.xyz)) for n, a in atoms(m).items())
  m.unset_riding_h_manager()
  m.setup_riding_h_manager(use_ideal_dihedral=True)
  for n, a in atoms(m).items():
    d = (matrix.col(a.xyz) - before[n]).length()
    assert d < 1.e-4, (n, d)

def exercise_var_64():
  # 64 deg, period 3: before, H11 in A and HB, HO1 in B and HA at 56 deg
  # (residual 8 = 2 x (64 - 60)), and the next placement moved them again
  for order in ORDERS:
    placed_twice(check_placed(order, 64., 3, False))

def exercise_const_64():
  # CONST 64: before, HO1 in A and HA, H11 in B and HB at -64 deg
  for order in ORDERS:
    placed_twice(check_placed(order, 64., 0, True))

def exercise_sign_insensitive():
  # ideals whose negative has the same minima within 1 deg (Var 60/3, 60.3/3,
  # 170/72 as GeoStd's backbone C-CA-N-H, CONST 180 and 0): unchanged, to
  # the bit. Seed 40: in the reordered forms the H stays at the mirror (-60,
  # -40) of the minimum nearest the seed, as before; that choice flips at the
  # next placement (not idempotent), kept so these H do not move
  mirror = {'H11': ('A', 'HB'), 'HO1': ('B', 'HA')}
  for v, period, const, nearest in ((60., 3, False, 60.),
      (60.3, 3, False, 60.3), (170., 72, False, 40.), (180., 0, True, 180.),
      (0., 0, True, 0.)):
    for order in ORDERS:
      m = check_placed(order, v, period, const, tol=0.7)
      for h in ('H11', 'HO1'):
        expected = nearest
        if not const and order in mirror[h]: expected = -nearest
        assert abs(torsion(m, h) - expected) < 1.e-6, (v, order, h,
          torsion(m, h))

def master_const_target(m, h):
  """master's CONST target of h: the stored ideal, negated if h is i4"""
  i = atoms(m)[h].i_seq
  for p in m.get_restraints_manager().geometry.const_dihedral_proxies:
    if i == p.i_seqs[0]: return p.angle_ideal
    if i == p.i_seqs[3]: return -p.angle_ideal

def exercise_near_const():
  # CONST within 0.5 deg of 0 or 180 (sign_sensitive's 1 deg tolerance, a
  # compatibility choice): master's target to the bit in every stored form,
  # 0.6 deg from the dictionary's in half of them
  for v in (0.3, -0.3, 179.7, -179.7):
    for order in ORDERS:
      m = check_placed(order, v, 0, True, tol=0.7)
      for h in ('H11', 'HO1'):
        d = geometry_restraints.angle_delta_deg(angle_1=torsion(m, h),
          angle_2=master_const_target(m, h), periodicity=1)
        assert abs(d) < 1.e-6, (v, order, h, torsion(m, h))

def exercise_var_refinement():
  # without use_ideal_dihedral (refinement) the H keep their current torsion
  for order in ORDERS:
    for v in (64., 60.):
      m = model(order, v, seed=40., use_ideal_dihedral=False)
      for h in ('H11', 'HO1'):
        assert abs(torsion(m, h) - 40.) < 0.1, (order, v, h, torsion(m, h))

def exercise_helper():
  from mmtbx.hydrogens.connectivity import riding_torsion
  # (H, a0, a1, t) and (t, a1, a0, H): riding order; the swapped forms: -1
  assert riding_torsion((5, 1, 2, 7), 5, 1) == (7, 1)
  assert riding_torsion((0, 2, 1, 5), 5, 1) == (0, 1)
  assert riding_torsion((0, 1, 2, 5), 5, 1) == (0, -1)
  assert riding_torsion((5, 2, 1, 7), 5, 1) == (7, -1)
  assert riding_torsion((0, 5, 1, 7), 5, 1) is None   # H in the middle
  assert riding_torsion((5, 2, 3, 7), 5, 1) is None   # parent not in it

# methanesulfonamide: sulfonamide NH2 (outside A17's pyramidal scope: planar
# alg1a pair), torsion HN1-N1-S1-C1 -67.16 deg, period 3 (AG4's HN9). Sorted
# as any input: N1 before S1, so the registry stores (C1, N1, S1, HN1)
def _msa():
  s = {'S1': (10., 10., 10.), 'C1': (11.77, 10., 10.)}
  s['N1'] = nerf((11.77, 11., 10.), s['C1'], s['S1'], 1.62, 106.0, 180.)
  s['O1'] = nerf(s['N1'], s['C1'], s['S1'], 1.43, 108.0, 120.)
  s['O2'] = nerf(s['N1'], s['C1'], s['S1'], 1.43, 108.0, -120.)
  for n, t in (('H1', 180.), ('H2', 60.), ('H3', -60.)):
    s[n] = nerf(s['N1'], s['S1'], s['C1'], 0.97, 109.5, t)
  s['HN1'] = nerf(s['C1'], s['S1'], s['N1'], 0.86, 112.0, -40.)
  s['HN2'] = nerf(s['C1'], s['S1'], s['N1'], 0.86, 112.0, 100.)
  lines = ['CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1']
  for k, n in enumerate(['S1', 'C1', 'N1', 'O1', 'O2', 'H1', 'H2', 'H3', 'HN1',
      'HN2']):
    e = 'H' if n.startswith('H') else n[0]
    lines.append('HETATM%5d %-4s LG9 A   1    %8.3f%8.3f%8.3f  1.00 20.00'
      '          %2s' % ((k + 1, ' ' + n) + s[n] + (e,)))
  pdb_str = '\n'.join(lines + ['END']) + '\n'
  atoms = [('S1', 'S'), ('C1', 'C'), ('N1', 'N'), ('O1', 'O'), ('O2', 'O'),
    ('H1', 'H'), ('H2', 'H'), ('H3', 'H'), ('HN1', 'H'), ('HN2', 'H')]
  bonds = [('S1', 'C1', 'single', 1.77), ('S1', 'N1', 'single', 1.62),
    ('S1', 'O1', 'double', 1.43), ('S1', 'O2', 'double', 1.43),
    ('C1', 'H1', 'single', 0.97), ('C1', 'H2', 'single', 0.97),
    ('C1', 'H3', 'single', 0.97), ('N1', 'HN1', 'single', 0.86),
    ('N1', 'HN2', 'single', 0.86)]
  angles = [('C1', 'S1', 'N1', 106.0), ('C1', 'S1', 'O1', 108.0),
    ('C1', 'S1', 'O2', 108.0), ('N1', 'S1', 'O1', 107.0),
    ('N1', 'S1', 'O2', 107.0), ('O1', 'S1', 'O2', 119.0),
    ('S1', 'N1', 'HN1', 112.0), ('S1', 'N1', 'HN2', 112.0),
    ('HN1', 'N1', 'HN2', 115.0)] + [(a, 'C1', b, 109.5) for a, b in [
    ('S1', 'H1'), ('S1', 'H2'), ('S1', 'H3'), ('H1', 'H2'), ('H1', 'H3'),
    ('H2', 'H3')]]
  L = ['data_comp_list', 'loop_', '_chem_comp.id',
    '_chem_comp.three_letter_code', '_chem_comp.name', '_chem_comp.group',
    '_chem_comp.number_atoms_all', '_chem_comp.number_atoms_nh',
    '_chem_comp.desc_level', "LG9 LG9 'methanesulfonamide' ligand 10 5 .",
    '', 'data_comp_LG9', 'loop_', '_chem_comp_atom.comp_id',
    '_chem_comp_atom.atom_id', '_chem_comp_atom.type_symbol',
    '_chem_comp_atom.type_energy']
  L += ['LG9 %s %s %s' % (n, e, e) for n, e in atoms]
  L += ['loop_', '_chem_comp_bond.comp_id', '_chem_comp_bond.atom_id_1',
    '_chem_comp_bond.atom_id_2', '_chem_comp_bond.type',
    '_chem_comp_bond.value_dist', '_chem_comp_bond.value_dist_esd']
  L += ['LG9 %s %s %s %.3f 0.020' % b for b in bonds]
  L += ['loop_', '_chem_comp_angle.comp_id', '_chem_comp_angle.atom_id_1',
    '_chem_comp_angle.atom_id_2', '_chem_comp_angle.atom_id_3',
    '_chem_comp_angle.value_angle', '_chem_comp_angle.value_angle_esd']
  L += ['LG9 %s %s %s %.1f 3.0' % a for a in angles]
  L += ['loop_', '_chem_comp_tor.comp_id', '_chem_comp_tor.id',
    '_chem_comp_tor.atom_id_1', '_chem_comp_tor.atom_id_2',
    '_chem_comp_tor.atom_id_3', '_chem_comp_tor.atom_id_4',
    '_chem_comp_tor.value_angle', '_chem_comp_tor.value_angle_esd',
    '_chem_comp_tor.period', 'LG9 var_1 HN1 N1 S1 C1 -67.16 30.0 3']
  return pdb_str, '\n'.join(L) + '\n'

def exercise_planar_nh2_refinement():
  # the planar alg1a path takes its target also in refinement (use_ideal_
  # dihedral False). Before: the target of the reordered proxy was the mirror
  # of the minimum nearest the current torsion, so the NH2 turned by ~120 deg
  # at every new manager (AG4, A17 step B1 item N2)
  pdb_str, cif_str = _msa()
  m = mmtbx.model.manager(model_input=iotbx.pdb.input(source_info=None,
    lines=pdb_str), log=null_out())
  m.set_restraint_objects([('lg9.cif',
    iotbx.cif.reader(input_string=cif_str).model())])
  m.process(make_restraints=True)
  names = [a.name.strip() for a in m.get_hierarchy().atoms()]
  assert [p for p in stored(m) if 'HN1' in p] == [
    ('C1', 'N1', 'S1', 'HN1')], stored(m)
  par = None
  for k in range(3):
    before = m.get_sites_cart().deep_copy()
    m.unset_riding_h_manager()
    m.setup_riding_h_manager()
    par = m.riding_h_manager.h_parameterization
    assert par[names.index('HN1')].htype == 'alg1a'
    after = m.get_sites_cart()
    t = math.degrees(dihedral_angle(sites=[after[names.index(n)]
      for n in ('HN1', 'N1', 'S1', 'C1')]))
    r = geometry_restraints.angle_delta_deg(angle_1=t, angle_2=-67.16,
      periodicity=3)
    assert abs(r) < 0.01, (k, t, r)
    if k > 0:
      d = max((matrix.col(after[i]) - matrix.col(before[i])).length()
        for i in range(len(names)))
      assert d < 1.e-4, (k, d)

if __name__ == '__main__':
  exercise_stored_forms()
  exercise_helper()
  exercise_var_64()
  exercise_const_64()
  exercise_sign_insensitive()
  exercise_near_const()
  exercise_var_refinement()
  exercise_planar_nh2_refinement()
  print('OK')
