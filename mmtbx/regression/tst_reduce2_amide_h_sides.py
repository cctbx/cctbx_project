"""reduce2: the Asn/Gln flip Mover keeps each supplied amide H on its own
side of the NH2 (A48). Before, MoverAmideFlip put the first H of the bonded
list at +120 degrees and the second at -120 degrees from C-N, so supplied
names opposite reduce2's convention were exchanged, and with them the
occupancy and B of each site.
Fragment: 7T5E Ile B 128 - Gly B 130 as deposited (H in A, D in B); Asn B 129
HD21/DD21 trans to CB, HD22/DD22 cis.
"""
from __future__ import absolute_import, division, print_function
import math
import os
import iotbx.pdb
from scitbx import matrix
from six.moves import cStringIO as StringIO
from iotbx.cli_parser import run_program
from mmtbx.programs import reduce2
from mmtbx.reduce import Movers

pdb_7t5e = """\
CRYST1   68.296   42.269   70.413  90.00  98.47  90.00 P 1 21 1
ATOM      1  N   ILE B 128     -11.939 -10.707   6.992  1.00 14.18           N
ATOM      2  CA  ILE B 128     -11.852 -11.885   6.133  1.00 18.40           C
ATOM      3  C   ILE B 128     -11.674 -13.157   6.960  1.00 22.51           C
ATOM      4  O   ILE B 128     -12.300 -14.198   6.698  1.00 16.15           O
ATOM      5  CB  ILE B 128     -10.696 -11.707   5.127  1.00 17.74           C
ATOM      6  CG1 ILE B 128     -11.076 -10.703   4.034  1.00 13.07           C
ATOM      7  CG2 ILE B 128     -10.357 -13.013   4.481  1.00 18.96           C
ATOM      8  CD1 ILE B 128      -9.833 -10.231   3.207  1.00 16.33           C
ATOM      9  HA  ILE B 128     -12.769 -12.006   5.573  1.00 16.68           H
ATOM     10  HB  ILE B 128      -9.828 -11.326   5.668  1.00 14.50           H
ATOM     11 HG12 ILE B 128     -11.773 -11.179   3.354  1.00 13.03           H
ATOM     12 HG13 ILE B 128     -11.536  -9.806   4.469  1.00 15.84           H
ATOM     13 HG21 ILE B 128      -9.457 -12.847   3.883  1.00 21.15           H
ATOM     14 HG22 ILE B 128     -10.142 -13.792   5.212  1.00 20.43           H
ATOM     15 HG23 ILE B 128     -11.192 -13.307   3.844  1.00 18.46           H
ATOM     16 HD11 ILE B 128      -9.769 -10.727   2.239  1.00 21.58           H
ATOM     17 HD12 ILE B 128      -9.924  -9.164   3.040  1.00 14.73           H
ATOM     18 HD13 ILE B 128      -8.920 -10.415   3.771  1.00 22.14           H
ATOM     19  H  AILE B 128     -11.234  -9.998   6.763  0.76 17.57           H
ATOM     20  D  BILE B 128     -11.242  -9.992   6.801  0.24 17.42           D
ATOM     21  N   ASN B 129     -10.830 -13.095   7.975  1.00 20.98           N
ATOM     22  CA  ASN B 129     -10.524 -14.260   8.767  1.00 22.96           C
ATOM     23  C   ASN B 129     -11.381 -14.327   9.999  1.00 21.46           C
ATOM     24  O   ASN B 129     -11.168 -15.212  10.819  1.00 25.93           O
ATOM     25  CB  ASN B 129      -9.046 -14.275   9.144  1.00 21.29           C
ATOM     26  CG  ASN B 129      -8.185 -14.465   7.935  1.00 31.63           C
ATOM     27  OD1 ASN B 129      -8.225 -15.500   7.277  1.00 33.36           O
ATOM     28  ND2 ASN B 129      -7.429 -13.464   7.616  1.00 27.73           N
ATOM     29  HA  ASN B 129     -10.705 -15.159   8.185  1.00 20.54           H
ATOM     30  HB2 ASN B 129      -8.775 -13.330   9.623  1.00 20.76           H
ATOM     31  HB3 ASN B 129      -8.837 -15.085   9.840  1.00 32.92           H
ATOM     32  H  AASN B 129     -10.341 -12.249   8.269  0.04 21.45           H
ATOM     33 HD21AASN B 129      -6.821 -13.533   6.799  0.20 30.81           H
ATOM     34 HD22AASN B 129      -7.430 -12.614   8.182  0.06 24.66           H
ATOM     35  D  BASN B 129     -10.363 -12.236   8.288  0.96 21.65           D
ATOM     36 DD21BASN B 129      -6.816 -13.542   6.805  0.80 31.44           D
ATOM     37 DD22BASN B 129      -7.441 -12.630   8.197  0.94 24.42           D
ATOM     38  N   GLY B 130     -12.321 -13.400  10.147  1.00 18.70           N
ATOM     39  CA  GLY B 130     -13.229 -13.375  11.280  1.00 18.93           C
ATOM     40  C   GLY B 130     -14.660 -13.609  10.885  1.00 15.95           C
ATOM     41  O   GLY B 130     -15.588 -13.021  11.456  1.00 19.10           O
ATOM     42  HA2 GLY B 130     -12.938 -14.157  11.986  1.00 26.22           H
ATOM     43  HA3 GLY B 130     -13.170 -12.416  11.804  1.00 19.20           H
ATOM     44  D  BGLY B 130     -12.458 -12.622   9.509  1.00 20.06           D
END
"""

def xray_pdb():
  '''pdb_7t5e without D; the H in blank atom groups at occupancy 1; Asn B 129
  HD21 and HD22 ideal on their input sides: in the CG-ND2-OD1 plane, 120
  degrees from ND2-CG, 0.86 A (deposited: 1.4 degrees off, 0.02 A).'''
  lines = []
  for l in pdb_7t5e.splitlines():
    if l.startswith('ATOM'):
      if l[76:78].strip() == 'D': continue
      l = l[:16] + ' ' + l[17:54] + '  1.00' + l[60:]
    lines.append(l)
  h = iotbx.pdb.input(source_info=None, lines=lines).construct_hierarchy()
  at = asn_129(h)
  n = matrix.col(at[('ND2', '')].xyz)
  c = matrix.col(at[('CG', '')].xyz)
  o = matrix.col(at[('OD1', '')].xyz)
  normal = (n - c).cross(o - c).normalize()
  to_c = (c - n).normalize() * 0.86
  for name in ('HD21', 'HD22'):
    a = at[(name, '')]
    ideal = [n + to_c.rotate_around_origin(normal, s * math.pi / 180)
      for s in (120, -120)]
    a.xyz = tuple(min(ideal, key=lambda p: dist(p, a.xyz)))
  return h.as_pdb_string(crystal_symmetry=iotbx.pdb.input(source_info=None,
    lines=pdb_7t5e).crystal_symmetry())

NEUTRON = ['neutron_option=hd_and_h', 'use_neutron_distances=True']

def run(prefix, pdb_str, args=()):
  '''reduce2 approach=optimize add_flip_movers=True on pdb_str.
  :return: output hierarchy.'''
  fn, out = prefix + '_in.pdb', prefix + '_out.pdb'
  with open(fn, 'w') as f:
    f.write(pdb_str)
  if os.path.exists(out): os.remove(out)
  run_program(program_class=reduce2.Program, logger=StringIO(), args=[fn,
    'approach=optimize', 'add_flip_movers=True', 'output.filename=' + out,
    'output.description_file_name=' + prefix + '.txt',
    'output.overwrite=True'] + list(args))
  return iotbx.pdb.input(out).construct_hierarchy()

def asn_129(hierarchy):
  '''{(name, altloc): atom} of Asn B 129.'''
  return dict(((a.name.strip(), a.parent().altloc), a)
    for rg in hierarchy.residue_groups() if rg.resseq.strip() == '129'
    for a in rg.atoms())

def dist(a, b):
  return (matrix.col(a) - matrix.col(b)).length()

def torsion(ref, name):
  '''CB-CG-ND2-<name> in degrees, -180..180.'''
  return matrix.dihedral_angle(sites=[ref[('CB', '')].xyz,
    ref[('CG', '')].xyz, ref[('ND2', '')].xyz, ref[name].xyz], deg=True)

def side(degrees):
  '''"trans" or "cis" to CB.'''
  if abs(degrees) > 90: return 'trans'
  return 'cis'

def test_xray_names_keep_their_sites():
  # Unflipped: each name stays on its input position.
  inp = asn_129(iotbx.pdb.input(source_info=None,
    lines=xray_pdb()).construct_hierarchy())
  out = asn_129(run('tst_r2_amide_xray', xray_pdb()))
  for heavy in ('ND2', 'OD1'):
    assert dist(out[(heavy, '')].xyz, inp[(heavy, '')].xyz) < 0.001, heavy
  for name in ('HD21', 'HD22'):
    d = dist(out[(name, '')].xyz, inp[(name, '')].xyz)
    assert d < 0.01, (name, d)

def test_hd_and_h_values_keep_their_sites():
  # Each written H/D is on the site that carries its name in the input, with
  # that site's occupancy and B.
  inp = asn_129(iotbx.pdb.input(source_info=None,
    lines=pdb_7t5e).construct_hierarchy())
  out = asn_129(run('tst_r2_amide_hd', pdb_7t5e, NEUTRON))
  sites = [k for k in inp if k[0] in ('HD21', 'HD22', 'DD21', 'DD22')]
  assert len(sites) == 4
  for k in sites:
    a = out[k]
    nearest = min(sites, key=lambda s: dist(a.xyz, inp[s].xyz))
    assert nearest[0][2:] == k[0][2:], (k, nearest)
    assert dist(a.xyz, inp[k].xyz) < 0.1, (k, dist(a.xyz, inp[k].xyz))
    assert (a.occ, a.b) == (inp[k].occ, inp[k].b), k

def test_flipped_names_keep_their_sides():
  # Flipped: ND2 on the input OD1 position, each H on the same side of CB as
  # in the input.
  inp = asn_129(iotbx.pdb.input(source_info=None,
    lines=xray_pdb()).construct_hierarchy())
  out = asn_129(run('tst_r2_amide_flip', xray_pdb(),
    ['set_flip_states=1 . B ASN 129 Flipped AnglesAdjusted']))
  assert dist(out[('ND2', '')].xyz, inp[('OD1', '')].xyz) < 0.1
  assert side(torsion(inp, ('HD21', ''))) == 'trans'
  assert side(torsion(inp, ('HD22', ''))) == 'cis'
  for name in ('HD21', 'HD22'):
    t_in, t_out = torsion(inp, (name, '')), torsion(out, (name, ''))
    assert side(t_out) == side(t_in), (name, t_in, t_out)

def test_mover_unequal_lengths():
  # Constructor level: HD21 at 0.93 A and HD22 at 1.07 A along their input
  # directions, listed in i_seq order. In both coarse states each H keeps its
  # own N-H length and its side of CB.
  h = iotbx.pdb.input(source_info=None, lines=xray_pdb()).construct_hierarchy()
  at = asn_129(h)
  for a in at.values():
    a.element = a.element.strip()
  n = matrix.col(at[('ND2', '')].xyz)
  length = {'HD21': 0.93, 'HD22': 1.07}
  for name in length:
    a = at[(name, '')]
    a.xyz = tuple(n + (matrix.col(a.xyz) - n).normalize() * length[name])
  inp = dict((k, a.xyz) for k, a in at.items())
  g = lambda name: at[(name, '')]
  bonded = {
    g('ND2'): [g('HD21'), g('HD22'), g('CG')],
    g('HD21'): [g('ND2')], g('HD22'): [g('ND2')],
    g('CG'): [g('CB'), g('OD1'), g('ND2')], g('OD1'): [g('CG')],
    g('CB'): [g('CA'), g('CG'), g('HB2'), g('HB3')],
    g('HB2'): [g('CB')], g('HB3'): [g('CB')], g('CA'): [g('CB')]}
  mover = Movers.MoverAmideFlip(g('ND2'), 'CA', bonded, 0.5)
  coarse = mover.CoarsePositions()
  index = dict((a.name.strip(), i) for i, a in enumerate(coarse.atoms))
  def side_of(pos, name):
    return side(matrix.dihedral_angle(sites=[inp[('CB', '')],
      inp[('CG', '')], pos[index['ND2']], pos[index[name]]], deg=True))
  ref = {'HD21': 'trans', 'HD22': 'cis'}
  for state, n_name in ((0, 'ND2'), (1, 'OD1')):
    pos = coarse.positions[state]
    for name in length:
      d = dist(pos[index[name]], pos[index['ND2']])
      assert abs(d - length[name]) < 1e-4, (state, name, d)
      assert side_of(pos, name) == ref[name], (state, name)
    # the flip moves ND2 onto the input OD1 site
    assert dist(pos[index['ND2']], inp[(n_name, '')]) < 1e-6, state

if __name__ == '__main__':
  test_xray_names_keep_their_sites()
  test_hd_and_h_values_keep_their_sites()
  test_flipped_names_keep_their_sides()
  test_mover_unequal_lengths()
  print('OK')
