"""Riding frame of an H whose torsion restraints include one with an H third
atom (A62c). In placement every H is still at its seed, so a target measured
from another H is arbitrary. A torsion with a heavy third atom sets the
target over one with an H third atom, whatever their order in the
dictionary; with H-third torsions only, an alg1b H or a propeller is
framed on a heavy atom bonded to the parent's neighbour, with the staggered
default (as for the amine NH2, A17); with no such heavy atom the H third
stays. A planar NH2 keeps its H third: its target is read in refinement too.
Fixtures: ethanol (HO1 with C1-C2-O1-HO1 and H21-C2-O1-HO1, as TRE's HO6),
methanol (no heavy atom on C1 but O1) and formamide (planar NH2, torsion
HN1-N1-C1-H1)."""
from __future__ import absolute_import, division, print_function
import math
import iotbx.pdb
import iotbx.cif
import mmtbx.model
from cctbx import geometry_restraints
from scitbx import matrix
from scitbx.math import dihedral_angle
from libtbx.utils import null_out

def nerf(a, b, c, d, ang, tor):
  a, b, c = [matrix.col(v) for v in (a, b, c)]
  bc = (c - b).normalize(); n = (b - a).cross(bc).normalize(); m = n.cross(bc)
  A, T = math.radians(ang), math.radians(tor)
  return tuple(c + (-d * math.cos(A)) * bc + (d * math.sin(A) * math.cos(T)) * m
    + (d * math.sin(A) * math.sin(T)) * n)

def ethanol(seed):
  """C1 methyl, C2 with H21/H22, hydroxyl O1-HO1; HO1 at seed from C1"""
  s = {'C1': (10., 10., 10.), 'C2': (11.52, 10., 10.)}
  s['O1'] = nerf((10., 11., 10.), s['C1'], s['C2'], 1.43, 109.5, 180.)
  for n, t in (('H11', 60.), ('H12', 180.), ('H13', -60.)):
    s[n] = nerf(s['O1'], s['C2'], s['C1'], 0.97, 109.5, t)
  for n, t in (('H21', 120.), ('H22', -120.)):
    s[n] = nerf(s['O1'], s['C1'], s['C2'], 0.97, 109.5, t)
  s['HO1'] = nerf(s['C1'], s['C2'], s['O1'], 0.85, 109.5, seed)
  bonds = [('C1', 'C2', 1.52), ('C2', 'O1', 1.43), ('C1', 'H11', 0.97),
    ('C1', 'H12', 0.97), ('C1', 'H13', 0.97), ('C2', 'H21', 0.97),
    ('C2', 'H22', 0.97), ('O1', 'HO1', 0.85)]
  return s, bonds

def methanol(seed):
  """C1 methyl on O1-HO1; HO1 at seed from H11"""
  s = {'C1': (10., 10., 10.), 'O1': (11.43, 10., 10.)}
  for n, t in (('H11', 60.), ('H12', 180.), ('H13', -60.)):
    s[n] = nerf((10., 11., 10.), s['O1'], s['C1'], 0.97, 109.5, t)
  s['HO1'] = nerf(s['H11'], s['C1'], s['O1'], 0.85, 109.5, seed)
  bonds = [('C1', 'O1', 1.43), ('C1', 'H11', 0.97), ('C1', 'H12', 0.97),
    ('C1', 'H13', 0.97), ('O1', 'HO1', 0.85)]
  return s, bonds

def formamide(t_hn1):
  """H1-C1(=O1)-N1(HN1)(HN2), planar, all angles 120; HN1 at O1-C1-N1-HN1
  t_hn1, HN2 opposite"""
  s = {'C1': (10., 10., 10.), 'N1': (11.33, 10., 10.)}
  s['O1'] = nerf((10., 11., 10.), s['N1'], s['C1'], 1.23, 120., 0.)
  s['H1'] = nerf((10., 11., 10.), s['N1'], s['C1'], 0.97, 120., 180.)
  s['HN1'] = nerf(s['O1'], s['C1'], s['N1'], 0.86, 120., t_hn1)
  s['HN2'] = nerf(s['O1'], s['C1'], s['N1'], 0.86, 120., t_hn1 + 180.)
  bonds = [('C1', 'N1', 1.33), ('C1', 'O1', 1.23), ('C1', 'H1', 0.97),
    ('N1', 'HN1', 0.86), ('N1', 'HN2', 0.86)]
  return s, bonds

FORMAMIDE_TE = dict(C1='C', O1='O', N1='NH2')

def angles(bonds):
  nb = {}
  for a, b, d in bonds:
    nb.setdefault(a, []).append(b); nb.setdefault(b, []).append(a)
  result = []
  for c, ns in sorted(nb.items()):
    for i in range(len(ns)):
      for j in range(i + 1, len(ns)):
        result.append((ns[i], c, ns[j]))
  return result

def pdb(s):
  lines = ['CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1']
  for k, n in enumerate(sorted(s)):
    e = 'H' if n.startswith('H') else n[0]
    lines.append('HETATM%5d %-4s LG9 A   1    %8.3f%8.3f%8.3f  1.00 20.00'
      '          %2s' % ((k + 1, ' ' + n) + s[n] + (e,)))
  return '\n'.join(lines + ['END']) + '\n'

def cif(s, bonds, torsions, angle=109.5, te=None):
  """torsions: (atoms, ideal, period), period 0 for CONST, in this order"""
  te = te or dict(C1='CH3', C2='CH2', O1='OH1')
  L = ['data_comp_list', 'loop_', '_chem_comp.id',
    '_chem_comp.three_letter_code', '_chem_comp.name', '_chem_comp.group',
    '_chem_comp.number_atoms_all', '_chem_comp.number_atoms_nh',
    '_chem_comp.desc_level', "LG9 LG9 'test' ligand %d 3 ." % len(s), '',
    'data_comp_LG9', 'loop_', '_chem_comp_atom.comp_id',
    '_chem_comp_atom.atom_id', '_chem_comp_atom.type_symbol',
    '_chem_comp_atom.type_energy']
  L += ['LG9 %s %s %s' % (n, 'H' if n.startswith('H') else n[0],
    te.get(n, 'H')) for n in sorted(s)]
  L += ['loop_', '_chem_comp_bond.comp_id', '_chem_comp_bond.atom_id_1',
    '_chem_comp_bond.atom_id_2', '_chem_comp_bond.type',
    '_chem_comp_bond.value_dist', '_chem_comp_bond.value_dist_esd']
  L += ['LG9 %s %s single %.3f 0.020' % b for b in bonds]
  L += ['loop_', '_chem_comp_angle.comp_id', '_chem_comp_angle.atom_id_1',
    '_chem_comp_angle.atom_id_2', '_chem_comp_angle.atom_id_3',
    '_chem_comp_angle.value_angle', '_chem_comp_angle.value_angle_esd']
  L += ['LG9 %s %s %s %.1f 3.0' % (a + (angle,)) for a in angles(bonds)]
  L += ['loop_', '_chem_comp_tor.comp_id', '_chem_comp_tor.id',
    '_chem_comp_tor.atom_id_1', '_chem_comp_tor.atom_id_2',
    '_chem_comp_tor.atom_id_3', '_chem_comp_tor.atom_id_4',
    '_chem_comp_tor.value_angle', '_chem_comp_tor.value_angle_esd',
    '_chem_comp_tor.period']
  for k, (names, v, p) in enumerate(torsions):
    L.append('LG9 %s_%d %s %s %s %s %.2f %s %d' % (('CONST' if p == 0 else
      'var', k + 1) + tuple(names) + (v, '0.0' if p == 0 else '10.0', p)))
  return '\n'.join(L) + '\n'

def model(s, bonds, torsions, angle=109.5, te=None, **kw):
  inp = iotbx.pdb.input(source_info=None, lines=pdb(s))
  m = mmtbx.model.manager(model_input=inp, log=null_out())
  m.set_restraint_objects([('lg9.cif', iotbx.cif.reader(
    input_string=cif(s, bonds, torsions, angle, te)).model())])
  m.process(make_restraints=True)
  m.setup_riding_h_manager(**kw)
  return m

def atoms(m):
  return dict((a.name.strip(), a) for a in m.get_hierarchy().atoms())

def torsion(m, names):
  at = atoms(m)
  return math.degrees(dihedral_angle(sites=[at[n].xyz for n in names]))

def residual(m, names, v, period):
  return geometry_restraints.angle_delta_deg(angle_1=torsion(m, names),
    angle_2=v, periodicity=max(1, period))

def frame(m, h):
  """the riding frame atom (a2) of h"""
  at = atoms(m)
  names = dict((a.i_seq, a.name.strip()) for a in at.values())
  return names[m.riding_h_manager.h_parameterization[at[h].i_seq].a2]

HEAVY_T = (('C1', 'C2', 'O1', 'HO1'), 70., 1)
H_T = (('H21', 'C2', 'O1', 'HO1'), -20., 1)

def exercise_heavy_third_wins():
  # TRE HO6: a heavy-third and an H-third torsion for the same H, in either
  # order in the dictionary and with either seed: the heavy one is the target
  for torsions in ((HEAVY_T, H_T), (H_T, HEAVY_T)):
    for seed in (-100., 40., 170.):
      s, bonds = ethanol(seed)
      m = model(s, bonds, torsions, use_ideal_dihedral=True)
      assert frame(m, 'HO1') == 'C1', (torsions, seed, frame(m, 'HO1'))
      r = residual(m, *HEAVY_T)
      assert abs(r) < 0.01, (torsions, seed, r)

def exercise_heavy_third_wins_const():
  # the same with CONST torsions, heavy first and H-third last
  for torsions in (((HEAVY_T[0], 70., 0), (H_T[0], -20., 0)),
                   ((H_T[0], -20., 0), (HEAVY_T[0], 70., 0))):
    s, bonds = ethanol(40.)
    m = model(s, bonds, torsions, use_ideal_dihedral=True)
    assert frame(m, 'HO1') == 'C1', (torsions, frame(m, 'HO1'))
    assert abs(residual(m, HEAVY_T[0], 70., 0)) < 0.01

def exercise_methyl_heavy_third_wins():
  # a propeller: one H with a heavy-third torsion, a sibling with an H-third
  # one after it; the group's target is the heavy one
  torsions = [(('H11', 'C1', 'C2', 'O1'), 75., 1),
              (('H12', 'C1', 'C2', 'H21'), 10., 1)]
  for order in (torsions, torsions[::-1]):
    s, bonds = ethanol(40.)
    m = model(s, bonds, order, use_ideal_dihedral=True)
    assert abs(residual(m, ('H11', 'C1', 'C2', 'O1'), 75., 1)) < 0.01, order

def exercise_h_third_only():
  # H-third torsions only: frame C1 (bonded to C2), staggered default 180,
  # whatever the seed
  for seed in (-100., 40., 170.):
    s, bonds = ethanol(seed)
    m = model(s, bonds, [H_T], use_ideal_dihedral=True)
    assert frame(m, 'HO1') == 'C1', (seed, frame(m, 'HO1'))
    assert abs(residual(m, ('C1', 'C2', 'O1', 'HO1'), 180., 1)) < 0.01

def exercise_no_heavy_alternative():
  # methanol: C1 has no heavy neighbour but O1, so the H third stays
  t = (('H11', 'C1', 'O1', 'HO1'), 180., 1)
  s, bonds = methanol(40.)
  m = model(s, bonds, [t], use_ideal_dihedral=True)
  assert frame(m, 'HO1') == 'H11', frame(m, 'HO1')

def exercise_refinement_keeps_sites():
  # without use_ideal_dihedral (refinement) the H keep their torsion,
  # whichever frame atom, and a second idealization moves nothing
  for torsions in ((HEAVY_T, H_T), (H_T, HEAVY_T), (H_T,)):
    s, bonds = ethanol(40.)
    m = model(s, bonds, list(torsions))
    t = torsion(m, ('C1', 'C2', 'O1', 'HO1'))
    assert abs(t - 40.) < 0.1, (torsions, t)
    before = m.get_sites_cart().deep_copy()
    m.idealize_h_riding()
    d = (m.get_sites_cart() - before).norms()
    assert d.all_lt(1.e-4), (torsions, max(d))

def exercise_planar_nh2_h_third():
  # a planar NH2 (alg1a) whose only torsion has an H third atom: the H third
  # stays its frame. The planar path reads b1's target in every
  # idealization, refinement included, so a constant (the heavy frame at
  # 180) would swap HN1 and HN2 when HN1 is cis to O1, at every new manager.
  # Input at the restraint: nothing moves
  for seed, v in ((0., 180.), (180., 0.)):
    t = (('HN1', 'N1', 'C1', 'H1'), v, 1)
    s, bonds = formamide(seed)
    m = model(s, bonds, [t], angle=120., te=FORMAMIDE_TE)
    at = atoms(m)
    assert m.riding_h_manager.h_parameterization[at['HN1'].i_seq].htype == \
      'alg1a'
    assert frame(m, 'HN1') == 'H1', (seed, frame(m, 'HN1'))
    for n in ('HN1', 'HN2'):
      d = (matrix.col(at[n].xyz) - matrix.col(s[n])).length()
      assert d < 0.01, (seed, n, d)
    before = m.get_sites_cart().deep_copy()
    m.unset_riding_h_manager()
    m.setup_riding_h_manager()
    d = (m.get_sites_cart() - before).norms()
    assert d.all_lt(1.e-4), (seed, max(d))

if __name__ == '__main__':
  exercise_heavy_third_wins()
  exercise_heavy_third_wins_const()
  exercise_methyl_heavy_third_wins()
  exercise_h_third_only()
  exercise_no_heavy_alternative()
  exercise_refinement_keeps_sites()
  exercise_planar_nh2_h_third()
  print('OK')
