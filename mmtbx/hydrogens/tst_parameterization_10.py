"""Riding H: an amine NH2 (cif: N with 2 H and one sp3 C, N or O neighbour,
no plane restraint on the H) is two alg1a H, each at its own X-N-H and N-H,
with the dictionary H-N-H (A17). Before: planar alg1a pair (H-N-H 360 - 2
C-N-H, ~141 deg), or not parameterized when the ideal H-N-H is 107-111 deg.
Placement (use_ideal_dihedral) takes the torsion restraint minimum nearest
180 and does not depend on the seed sites or the frame; otherwise the NH2
keeps its current torsion. Not riding with ignore_h_with_dof (neutron
refinement). Unchanged: planar NH2 (Asn ND2), an aniline without a plane
restraint, and any XH2 without cif information (Lys with HZ3 missing).
Also master's path: an NH2 H in any plane restraint and X-N-H/H-N-H targets
no pyramid can meet. With a heavy atom bonded to X but only H-third torsions,
that atom is the frame (first H anti to it); with none (methylamine with
H-C-N-H torsions only) the frame is an H: pyramidal, orientation arbitrary."""
from __future__ import absolute_import, division, print_function
import math
import iotbx.pdb
import iotbx.cif
import mmtbx.model
from mmtbx.hydrogens import riding, parameterization
from scitbx.array_family import flex
from scitbx import matrix
from libtbx.utils import null_out

# putrescine, CCD ideal coordinates; GeoStd: H-N-H 105.45, C-N-H 109.0/109.5,
# one torsion per N
pdb_put = """\
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1
HETATM    1  N2  PUT A   1       9.750   9.980  13.089  1.00 20.00           N
HETATM    2  C4  PUT A   1      10.535   9.980  11.849  1.00 20.00           C
HETATM    3  C3  PUT A   1       9.591   9.980  10.646  1.00 20.00           C
HETATM    4  C2  PUT A   1      10.409   9.980   9.354  1.00 20.00           C
HETATM    5  C1  PUT A   1       9.465   9.980   8.151  1.00 20.00           C
HETATM    6  N1  PUT A   1      10.250   9.980   6.911  1.00 20.00           N
HETATM    7 HN21 PUT A   1       9.266  10.864  13.126  1.00 20.00           H
HETATM    8 HN22 PUT A   1      10.410   9.980  13.852  1.00 20.00           H
HETATM    9  H41 PUT A   1      11.164  10.869  11.817  1.00 20.00           H
HETATM   10  H42 PUT A   1      11.164   9.090  11.817  1.00 20.00           H
HETATM   11  H31 PUT A   1       8.962   9.090  10.678  1.00 20.00           H
HETATM   12  H32 PUT A   1       8.962  10.869  10.678  1.00 20.00           H
HETATM   13  H21 PUT A   1      11.038  10.869   9.322  1.00 20.00           H
HETATM   14  H22 PUT A   1      11.038   9.090   9.322  1.00 20.00           H
HETATM   15  H11 PUT A   1       8.836   9.090   8.183  1.00 20.00           H
HETATM   16  H12 PUT A   1       8.836  10.869   8.183  1.00 20.00           H
HETATM   17 HN11 PUT A   1       9.590   9.980   6.148  1.00 20.00           H
HETATM   18 HN12 PUT A   1      10.735  10.864   6.874  1.00 20.00           H
END
"""

# ethanolamine (CCD ETA coordinates) with a user CIF: angles 109.5 at N, no
# torsions (an eLBOW-like dictionary)
pdb_lg1 = """\
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1
HETATM    1  CA  LG1 A   1       9.289  10.612   9.692  1.00 20.00           C
HETATM    2  N   LG1 A   1       8.532   9.426  10.114  1.00 20.00           N
HETATM    3  C   LG1 A   1      10.681  10.581  10.327  1.00 20.00           C
HETATM    4  O   LG1 A   1      11.412   9.464   9.816  1.00 20.00           O
HETATM    5  HA1 LG1 A   1       9.385  10.615   8.606  1.00 20.00           H
HETATM    6  HA2 LG1 A   1       8.764  11.511  10.013  1.00 20.00           H
HETATM    7  H   LG1 A   1       8.485   9.366  11.120  1.00 20.00           H
HETATM    8  H2  LG1 A   1       7.608   9.423   9.707  1.00 20.00           H
HETATM    9  HB1 LG1 A   1      11.211  11.502  10.086  1.00 20.00           H
HETATM   10  HB2 LG1 A   1      10.585  10.488  11.408  1.00 20.00           H
HETATM   11  HO  LG1 A   1      12.306   9.381  10.176  1.00 20.00           H
END
"""

cif_lg1 = """\
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
LG1 LG1 'ethanolamine' ligand 11 4 .
data_comp_LG1
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge
LG1 N   N NT2 0
LG1 CA  C CH2 0
LG1 C   C CH2 0
LG1 O   O OH1 0
LG1 H   H HNT2 0
LG1 H2  H HNT2 0
LG1 HA1 H HCH2 0
LG1 HA2 H HCH2 0
LG1 HB1 H HCH2 0
LG1 HB2 H HCH2 0
LG1 HO  H HOH1 0
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
LG1 N   CA  single 1.470 0.020
LG1 CA  C   single 1.520 0.020
LG1 C   O   single 1.420 0.020
LG1 N   H   single 0.890 0.020
LG1 N   H2  single 0.890 0.020
LG1 CA  HA1 single 0.970 0.020
LG1 CA  HA2 single 0.970 0.020
LG1 C   HB1 single 0.970 0.020
LG1 C   HB2 single 0.970 0.020
LG1 O   HO  single 0.840 0.020
loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
LG1 CA  N   H   109.5 3.0
LG1 CA  N   H2  109.5 3.0
LG1 H   N   H2  109.5 3.0
LG1 N   CA  C   109.5 3.0
LG1 N   CA  HA1 109.5 3.0
LG1 N   CA  HA2 109.5 3.0
LG1 C   CA  HA1 109.5 3.0
LG1 C   CA  HA2 109.5 3.0
LG1 HA1 CA  HA2 109.5 3.0
LG1 CA  C   O   109.5 3.0
LG1 CA  C   HB1 109.5 3.0
LG1 CA  C   HB2 109.5 3.0
LG1 O   C   HB1 109.5 3.0
LG1 O   C   HB2 109.5 3.0
LG1 HB1 C   HB2 109.5 3.0
LG1 C   O   HO  109.5 3.0
"""

# Asn, planar side-chain amide (control)
pdb_asn = """\
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1
ATOM      1  N   ASN A   1       8.991   9.346  10.623  1.00 10.00           N
ATOM      2  CA  ASN A   1      10.137  10.220  10.389  1.00 10.00           C
ATOM      3  C   ASN A   1      11.317   9.432   9.836  1.00 10.00           C
ATOM      4  O   ASN A   1      11.210   8.240   9.551  1.00 10.00           O
ATOM      5  CB  ASN A   1       9.759  11.333   9.417  1.00 10.00           C
ATOM      6  CG  ASN A   1      10.857  12.360   9.230  1.00 10.00           C
ATOM      7  OD1 ASN A   1      11.896  12.307   9.885  1.00 10.00           O
ATOM      8  ND2 ASN A   1      10.632  13.300   8.322  1.00 10.00           N
ATOM      9  H   ASN A   1       8.180   9.616  10.111  1.00 10.00           H
ATOM     10  HA  ASN A   1      10.423  10.629  11.231  1.00 10.00           H
ATOM     11  HB2 ASN A   1       8.973  11.778   9.736  1.00 10.00           H
ATOM     12  HB3 ASN A   1       9.537  10.952   8.564  1.00 10.00           H
ATOM     13 HD21 ASN A   1      11.235  14.010   8.165  1.00 10.00           H
ATOM     14 HD22 ASN A   1       9.865  13.270   7.873  1.00 10.00           H
END
"""

def model(pdb_str, cif_str=None, riding_h=True, **kw):
  m = mmtbx.model.manager(model_input=iotbx.pdb.input(source_info=None,
    lines=pdb_str), log=null_out())
  if cif_str is not None:
    m.set_restraint_objects([('lg.cif',
      iotbx.cif.reader(input_string=cif_str).model())])
  m.process(make_restraints=True)
  if riding_h:
    m.setup_riding_h_manager(**kw)
  return m

def ang(a, b, c):
  return math.degrees((matrix.col(a) - matrix.col(b)).angle(
    matrix.col(c) - matrix.col(b)))

def nh2(m, n_name, x_name, h_names, altloc=None):
  atoms = dict((a.name.strip(), a) for a in m.get_hierarchy().atoms()
    if altloc is None or a.parent().altloc in ('', altloc))
  par = m.riding_h_manager.h_parameterization
  n, x = atoms[n_name], atoms[x_name]
  h1, h2 = [atoms[h] for h in h_names]
  hnh = ang(h1.xyz, n.xyz, h2.xyz)
  xnh = [ang(x.xyz, n.xyz, h.xyz) for h in (h1, h2)]
  p = [matrix.col(v.xyz) for v in (x, h1, h2)]
  nrm = (p[1] - p[0]).cross(p[2] - p[0]).normalize()
  height = abs((matrix.col(n.xyz) - p[0]).dot(nrm))
  types = [par[h.i_seq].htype if par[h.i_seq] is not None else None
    for h in (h1, h2)]
  return hnh, xnh, height, types

def exercise_geostd_put():
  # alg1a placed both H planar (H-N-H ~141)
  m = model(pdb_put)
  for n, x, hs in (('N1', 'C1', ('HN11', 'HN12')),
                   ('N2', 'C4', ('HN21', 'HN22'))):
    hnh, xnh, height, types = nh2(m, n, x, hs)
    assert types == ['alg1a', 'alg1a'], types
    assert abs(hnh - 105.45) < 0.1, (n, hnh)
    for v in xnh:
      assert 108.9 < v < 109.7, (n, xnh)
    assert height > 0.3, (n, height)

def exercise_user_cif_109():
  # H-N-H 109.5 was refused (unk_list): H not parameterized
  m = model(pdb_lg1, cif_lg1)
  hnh, xnh, height, types = nh2(m, 'N', 'CA', ('H', 'H2'))
  assert types == ['alg1a', 'alg1a'], types
  assert abs(hnh - 109.5) < 0.1, hnh
  for v in xnh:
    assert abs(v - 109.5) < 0.1, xnh
  assert height > 0.3, height

def exercise_asn_planar():
  m = model(pdb_asn)
  hnh, xnh, height, types = nh2(m, 'ND2', 'CG', ('HD21', 'HD22'))
  assert types == ['alg1a', 'alg1a'], types
  assert abs(hnh + sum(xnh) - 360) < 0.5, (hnh, xnh)
  assert height < 0.02, height

def cif(code, atoms, bonds, angles, planes=(), torsions=()):
  """monomer cif: atoms (name, element), bonds (a, b, type, x-ray, neutron),
  angles (a, b, c, value), planes [atom names], torsions (a, b, c, d, value,
  period)"""
  lines = ["data_comp_list", "loop_", "_chem_comp.id",
    "_chem_comp.three_letter_code", "_chem_comp.name", "_chem_comp.group",
    "_chem_comp.number_atoms_all", "_chem_comp.number_atoms_nh",
    "_chem_comp.desc_level", "%s %s 'test' ligand %d %d ." % (code, code,
    len(atoms), len([a for a in atoms if a[1] != 'H'])), "",
    "data_comp_" + code, "loop_", "_chem_comp_atom.comp_id",
    "_chem_comp_atom.atom_id", "_chem_comp_atom.type_symbol",
    "_chem_comp_atom.type_energy"]
  lines += ["%s %s %s %s" % (code, n, e, e) for n, e in atoms]
  lines += ["loop_", "_chem_comp_bond.comp_id", "_chem_comp_bond.atom_id_1",
    "_chem_comp_bond.atom_id_2", "_chem_comp_bond.type",
    "_chem_comp_bond.value_dist", "_chem_comp_bond.value_dist_esd",
    "_chem_comp_bond.value_dist_neutron"]
  lines += ["%s %s %s %s %.3f 0.020 %.3f" % ((code,) + b) for b in bonds]
  lines += ["loop_", "_chem_comp_angle.comp_id", "_chem_comp_angle.atom_id_1",
    "_chem_comp_angle.atom_id_2", "_chem_comp_angle.atom_id_3",
    "_chem_comp_angle.value_angle", "_chem_comp_angle.value_angle_esd"]
  lines += ["%s %s %s %s %.1f 3.0" % ((code,) + a) for a in angles]
  if planes:
    lines += ["loop_", "_chem_comp_plane_atom.comp_id",
      "_chem_comp_plane_atom.plane_id", "_chem_comp_plane_atom.atom_id",
      "_chem_comp_plane_atom.dist_esd"]
    for k, plane in enumerate(planes):
      lines += ["%s plan-%d %s 0.020" % (code, k + 1, n) for n in plane]
  if torsions:
    lines += ["loop_", "_chem_comp_tor.comp_id", "_chem_comp_tor.id",
      "_chem_comp_tor.atom_id_1", "_chem_comp_tor.atom_id_2",
      "_chem_comp_tor.atom_id_3", "_chem_comp_tor.atom_id_4",
      "_chem_comp_tor.value_angle", "_chem_comp_tor.value_angle_esd",
      "_chem_comp_tor.period"]
    lines += ["%s var_%d %s %s %s %s %.1f 10.0 %d" % ((code, k + 1) + t)
      for k, t in enumerate(torsions)]
  return "\n".join(lines) + "\n"

def tetra(center, a, b, h1, h2):
  return [(a, center, h1), (a, center, h2), (b, center, h1), (b, center, h2),
          (h1, center, h2)]

# ethanolamine (pdb_lg1 coordinates) with unequal restraints on the two H:
# C-N-H 106.0 / 113.0, N-H 0.89 / 0.95, H-N-H 107.0
pdb_lg2 = pdb_lg1.replace('LG1', 'LG2')
cif_lg2 = cif('LG2',
  [('N', 'N'), ('CA', 'C'), ('C', 'C'), ('O', 'O'), ('H', 'H'), ('H2', 'H'),
   ('HA1', 'H'), ('HA2', 'H'), ('HB1', 'H'), ('HB2', 'H'), ('HO', 'H')],
  [('N', 'CA', 'single', 1.47, 1.47), ('CA', 'C', 'single', 1.52, 1.52),
   ('C', 'O', 'single', 1.42, 1.42), ('N', 'H', 'single', 0.89, 1.01),
   ('N', 'H2', 'single', 0.95, 1.04), ('CA', 'HA1', 'single', 0.97, 1.09),
   ('CA', 'HA2', 'single', 0.97, 1.09), ('C', 'HB1', 'single', 0.97, 1.09),
   ('C', 'HB2', 'single', 0.97, 1.09), ('O', 'HO', 'single', 0.84, 0.97)],
  [('CA', 'N', 'H', 106.0), ('CA', 'N', 'H2', 113.0), ('H', 'N', 'H2', 107.0)]
  + [(a, b, c, 109.5) for a, b, c in [('N', 'CA', 'C'), ('CA', 'C', 'O'),
     ('C', 'O', 'HO')] + tetra('CA', 'N', 'C', 'HA1', 'HA2') +
     tetra('C', 'CA', 'O', 'HB1', 'HB2')])

# ethanolamine with the NH2 H named HN1/HN2 and one plane holding HA1 first,
# then both NH2 H (connectivity.plane_h is keyed by the first H only)
_lg1_atoms = [('N', 'N'), ('CA', 'C'), ('C', 'C'), ('O', 'O'), ('HN1', 'H'),
  ('HN2', 'H'), ('HA1', 'H'), ('HA2', 'H'), ('HB1', 'H'), ('HB2', 'H'),
  ('HO', 'H')]
_lg1_bonds = [('N', 'CA', 'single', 1.47, 1.47),
  ('CA', 'C', 'single', 1.52, 1.52), ('C', 'O', 'single', 1.42, 1.42),
  ('N', 'HN1', 'single', 0.89, 1.01), ('N', 'HN2', 'single', 0.89, 1.01),
  ('CA', 'HA1', 'single', 0.97, 1.09), ('CA', 'HA2', 'single', 0.97, 1.09),
  ('C', 'HB1', 'single', 0.97, 1.09), ('C', 'HB2', 'single', 0.97, 1.09),
  ('O', 'HO', 'single', 0.84, 0.97)]
_lg1_heavy_angles = [(a, b, c, 109.5) for a, b, c in [('N', 'CA', 'C'),
  ('CA', 'C', 'O'), ('C', 'O', 'HO')] + tetra('CA', 'N', 'C', 'HA1', 'HA2') +
  tetra('C', 'CA', 'O', 'HB1', 'HB2')]
def _lg1_named(code):
  return (pdb_lg1.replace('LG1', code).replace(' H   %s' % code,
    ' HN1 %s' % code).replace(' H2  %s' % code, ' HN2 %s' % code))
pdb_lg5 = _lg1_named('LG5')
cif_lg5 = cif('LG5', _lg1_atoms, _lg1_bonds,
  [('CA', 'N', 'HN1', 109.5), ('CA', 'N', 'HN2', 109.5),
   ('HN1', 'N', 'HN2', 109.5)] + _lg1_heavy_angles,
  planes=[['N', 'CA', 'HA1', 'HN1', 'HN2']])

# X-N-H 125/125 and H-N-H 120: angle sum 370, no torsion gives it
pdb_lg7 = _lg1_named('LG7')
cif_lg7 = cif('LG7', _lg1_atoms, _lg1_bonds,
  [('CA', 'N', 'HN1', 125.0), ('CA', 'N', 'HN2', 125.0),
   ('HN1', 'N', 'HN2', 120.0)] + _lg1_heavy_angles)

# ethanolamine whose only NH2 torsion has an H third atom (HA1-CA-N-HN1, as
# CCD-fallback torsions H-C-N-H); C is a heavy atom bonded to CA
pdb_lg6 = _lg1_named('LG6')
cif_lg6 = cif('LG6', _lg1_atoms, _lg1_bonds,
  [('CA', 'N', 'HN1', 109.5), ('CA', 'N', 'HN2', 109.5),
   ('HN1', 'N', 'HN2', 106.0)] + _lg1_heavy_angles,
  torsions=[('HA1', 'CA', 'N', 'HN1', 175.0, 3)])

# methylamine, the same torsion: CA has no other heavy neighbour
pdb_lg4 = """\
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1
HETATM    1  N   LG4 A   1       9.551   9.823   9.872  1.00 20.00           N
HETATM    2  C   LG4 A   1      10.668   9.688  10.794  1.00 20.00           C
HETATM    3  HN1 LG4 A   1       8.693   9.927  10.398  1.00 20.00           H
HETATM    4  HN2 LG4 A   1       9.652  10.682   9.348  1.00 20.00           H
HETATM    5  H1  LG4 A   1      11.597   9.574  10.234  1.00 20.00           H
HETATM    6  H2  LG4 A   1      10.805  10.526  11.488  1.00 20.00           H
HETATM    7  H3  LG4 A   1      10.544   8.785  11.392  1.00 20.00           H
END
"""
cif_lg4 = cif('LG4', [('N', 'N'), ('C', 'C'), ('HN1', 'H'), ('HN2', 'H'),
  ('H1', 'H'), ('H2', 'H'), ('H3', 'H')],
  [('N', 'C', 'single', 1.455, 1.455), ('N', 'HN1', 'single', 0.89, 1.01),
   ('N', 'HN2', 'single', 0.89, 1.01), ('C', 'H1', 'single', 0.97, 1.09),
   ('C', 'H2', 'single', 0.97, 1.09), ('C', 'H3', 'single', 0.97, 1.09)],
  [('C', 'N', 'HN1', 109.5), ('C', 'N', 'HN2', 109.5),
   ('HN1', 'N', 'HN2', 106.0)] + [(a, 'C', b, 109.5) for a, b in (
   ('N', 'H1'), ('N', 'H2'), ('N', 'H3'), ('H1', 'H2'), ('H1', 'H3'),
   ('H2', 'H3'))],
  torsions=[('H1', 'C', 'N', 'HN1', 175.0, 3)])

# aniline whose cif makes the NH2 pyramidal (C-N-H 113, H-N-H 112, sum 338;
# 13 GeoStd anilines are like this) without a plane restraint on the H
def _aniline():
  atoms, lines = [], ["CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1"]
  for k in range(6):
    t = math.radians(60 * k)
    atoms.append(('C%d' % (k + 1), 'C', 1.39 * math.cos(t), 1.39 * math.sin(t)))
    if k > 0:
      atoms.append(('H%d' % (k + 1), 'H', 2.47 * math.cos(t), 2.47 * math.sin(t)))
  atoms += [('N', 'N', 2.79, 0.0), ('HN1', 'H', 3.235, 0.771),
            ('HN2', 'H', 3.235, -0.771)]
  for i, (n, e, x, y) in enumerate(atoms):
    lines.append("HETATM%5d %-4s LG3 A   1    %8.3f%8.3f%8.3f  1.00 20.00"
      "          %2s" % (i + 1, ' ' + n, x + 10, y + 10, 10.0, e))
  return "\n".join(lines) + "\nEND\n"
pdb_anl = _aniline()
_ring = ['C%d' % k for k in range(1, 7)]
cif_anl = cif('LG3', [(n, 'C') for n in _ring] + [('N', 'N')] +
  [('H%d' % k, 'H') for k in range(2, 7)] + [('HN1', 'H'), ('HN2', 'H')],
  [(_ring[k], _ring[(k + 1) % 6], 'aromatic', 1.39, 1.39) for k in range(6)]
  + [('C1', 'N', 'single', 1.40, 1.40), ('N', 'HN1', 'single', 0.89, 1.01),
     ('N', 'HN2', 'single', 0.89, 1.01)]
  + [('C%d' % k, 'H%d' % k, 'single', 0.93, 1.08) for k in range(2, 7)],
  [(_ring[k - 1], _ring[k], _ring[(k + 1) % 6], 120.0) for k in range(6)]
  + [('C2', 'C1', 'N', 120.0), ('C6', 'C1', 'N', 120.0)]
  + [(a, 'C%d' % k, 'H%d' % k, 120.0) for k in range(2, 7)
     for a in (_ring[k - 2], _ring[k % 6])]
  + [('C1', 'N', 'HN1', 113.0), ('C1', 'N', 'HN2', 113.0),
     ('HN1', 'N', 'HN2', 112.0)],
  planes=[_ring + ['N'] + ['H%d' % k for k in range(2, 7)]])

# Lys (tst_add_hydrogen_1) with HZ3 missing: an incomplete NH3, not an amine
pdb_lys_no_hz3 = """\
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1
ATOM     28  N   LYS A   3     -50.756 -15.600 -14.375  1.00105.98           N
ATOM     29  CA  LYS A   3     -50.135 -16.230 -15.528  1.00115.21           C
ATOM     30  C   LYS A   3     -50.110 -17.740 -15.357  1.00105.26           C
ATOM     31  O   LYS A   3     -50.331 -18.255 -14.264  1.00118.65           O
ATOM     32  CB  LYS A   3     -48.714 -15.704 -15.737  1.00110.03           C
ATOM     33  CG  LYS A   3     -47.784 -15.954 -14.563  1.00100.14           C
ATOM     34  CD  LYS A   3     -46.379 -15.452 -14.853  1.00101.77           C
ATOM     35  CE  LYS A   3     -45.456 -15.678 -13.663  1.00112.21           C
ATOM     36  NZ  LYS A   3     -44.063 -15.229 -13.941  1.00110.68           N
ATOM     37  H   LYS A   3     -50.849 -16.122 -13.698  1.00105.98           H
ATOM     38  HA  LYS A   3     -50.648 -16.015 -16.323  1.00115.21           H
ATOM     39  HB2 LYS A   3     -48.332 -16.141 -16.514  1.00110.03           H
ATOM     40  HB3 LYS A   3     -48.755 -14.746 -15.882  1.00110.03           H
ATOM     41  HG2 LYS A   3     -48.120 -15.488 -13.782  1.00100.14           H
ATOM     42  HG3 LYS A   3     -47.738 -16.907 -14.386  1.00100.14           H
ATOM     43  HD2 LYS A   3     -46.017 -15.929 -15.616  1.00101.77           H
ATOM     44  HD3 LYS A   3     -46.409 -14.501 -15.041  1.00101.77           H
ATOM     45  HE2 LYS A   3     -45.789 -15.178 -12.902  1.00112.21           H
ATOM     46  HE3 LYS A   3     -45.432 -16.625 -13.453  1.00112.21           H
ATOM     47  HZ1 LYS A   3     -43.731 -15.679 -14.633  1.00110.68           H
ATOM     48  HZ2 LYS A   3     -44.056 -14.360 -14.131  1.00110.68           H
END
"""

PUT_NH2 = (('N1', 'C1', ('HN11', 'HN12')), ('N2', 'C4', ('HN21', 'HN22')))

def check_put(m, altloc=None):
  for n, x, hs in PUT_NH2:
    hnh, xnh, height, types = nh2(m, n, x, hs, altloc=altloc)
    assert types == ['alg1a', 'alg1a'], (n, altloc, types)
    assert abs(hnh - 105.45) < 0.1, (n, altloc, hnh)
    assert height > 0.3, (n, altloc, height)

def exercise_unequal_restraints():
  # each H at its own X-N-H and N-H; H-N-H from the dictionary
  m = model(pdb_lg2, cif_lg2)
  hnh, xnh, height, types = nh2(m, 'N', 'CA', ('H', 'H2'))
  assert types == ['alg1a', 'alg1a'], types
  assert abs(hnh - 107.0) < 0.05, hnh
  assert abs(xnh[0] - 106.0) < 0.05 and abs(xnh[1] - 113.0) < 0.05, xnh
  atoms = dict((a.name.strip(), a) for a in m.get_hierarchy().atoms())
  for name, d in (('H', 0.89), ('H2', 0.95)):
    assert abs(atoms['N'].distance(atoms[name]) - d) < 1.e-3, name

def exercise_alt_conformers():
  # each conformer of a split residue on its own
  lines = []
  for alt, shift in (('A', 0.0), ('B', 0.4)):
    for l in pdb_put.splitlines():
      if not l.startswith('HETATM'): continue
      x, y = float(l[30:38]) + shift, float(l[38:46]) + 0.75 * shift
      lines.append(l[:16] + alt + l[17:30] + '%8.3f%8.3f' % (x, y) + l[46:54]
        + '  0.50' + l[60:])
  m = model(pdb_put.splitlines()[0] + '\n' + '\n'.join(lines) + '\nEND\n')
  for alt in ('A', 'B'):
    check_put(m, altloc=alt)

def _put_sites(m):
  return dict((a.name.strip(), matrix.col(a.xyz))
    for a in m.get_hierarchy().atoms())

def exercise_repeated_idealization():
  # outside placement the NH2 keeps its torsion; a second idealization with
  # a new manager (next macro-cycle) moves nothing
  from scitbx.math import dihedral_angle
  m = model(pdb_put, use_ideal_dihedral=True)
  check_put(m)
  h = m.get_hierarchy()
  s = _put_sites(m)
  r = (s['N1'] - s['C1']).axis_and_angle_as_r3_rotation_matrix(40, deg=True)
  for a in h.atoms():
    if a.name.strip() in ('HN11', 'HN12'):
      a.xyz = tuple(r * (matrix.col(a.xyz) - s['N1']) + s['N1'])
  def torsion(m):
    s = _put_sites(m)
    return dihedral_angle(sites=[s['C2'], s['C1'], s['N1'], s['HN11']],
      deg=True)
  m2 = model(h.as_pdb_string(crystal_symmetry=m.crystal_symmetry()),
    riding_h=False)
  t0 = torsion(m2)
  m2.setup_riding_h_manager(idealize=True)
  check_put(m2)
  assert abs(torsion(m2) - t0) < 0.05, (t0, torsion(m2))
  sc = m2.get_sites_cart().deep_copy()
  m2.unset_riding_h_manager()
  m2.setup_riding_h_manager(idealize=True)
  d = flex.sqrt((m2.get_sites_cart() - sc).dot())
  assert flex.max(d) < 1.e-4, flex.max(d)

def exercise_seed_frame_invariance():
  # placement seeds every new H at one site (the residue's mean + 0.5 A per
  # axis in the model frame): the NH2 H, by name, may not depend on the seed
  # or on a rigid-body move of the input
  rots = [matrix.sqr((1, 0, 0, 0, 1, 0, 0, 0, 1))] + [
    matrix.col(axis).normalize().axis_and_angle_as_r3_rotation_matrix(
      angle, deg=True) for axis, angle in (((1, 2, 3), 70), ((-2, 1, 0.5), 155),
      ((0, 0, 1), 250))]
  seeds = [(0.5, 0.5, 0.5), (-0.4, 0.9, 0.2)]
  h0 = iotbx.pdb.input(source_info=None, lines=pdb_put).construct_hierarchy()
  center = matrix.col(h0.atoms().extract_xyz().mean())
  nh2_h = ('HN11', 'HN12', 'HN21', 'HN22')
  ref = None
  for rot in rots:
    for seed in seeds:
      h = h0.deep_copy()
      for a in h.atoms():
        a.xyz = tuple(rot * (matrix.col(a.xyz) - center) + center)
      heavy = flex.vec3_double([a.xyz for a in h.atoms()
        if a.element.strip() != 'H'])
      for a in h.atoms():
        if a.name.strip() in nh2_h:
          a.xyz = tuple(matrix.col(heavy.mean()) + matrix.col(seed))
      m = model(h.as_pdb_string(crystal_symmetry=iotbx.pdb.input(
        source_info=None, lines=pdb_put).crystal_symmetry()),
        use_ideal_dihedral=True)
      check_put(m)
      s = _put_sites(m)
      back = dict((n, rot.transpose() * (s[n] - center) + center)
        for n in nh2_h)
      if ref is None:
        ref = back
        continue
      for n in nh2_h:
        # PDB coordinates: 0.001 A rounding of the rotated input
        assert (back[n] - ref[n]).length() < 0.005, (n, rot, seed,
          (back[n] - ref[n]).length())

def exercise_ignore_h_with_dof():
  # neutron refinement: rotatable H are refined individually, not ridden
  for pdb_str, cif_str, hs, ch2 in ((pdb_put, None,
      ('HN11', 'HN12', 'HN21', 'HN22'), 'H11'),
      (pdb_lg1, cif_lg1, ('H', 'H2'), 'HA1')):
    m = model(pdb_str, cif_str, idealize=False, ignore_h_with_dof=True)
    par = m.riding_h_manager.h_parameterization
    atoms = dict((a.name.strip(), a) for a in m.get_hierarchy().atoms())
    assert [par[atoms[h].i_seq] for h in hs] == [None] * len(hs), \
      [par[atoms[h].i_seq].htype for h in hs if par[atoms[h].i_seq]]
    assert par[atoms[ch2].i_seq] is not None, ch2

def exercise_no_cif():
  # without cif information an XH2 keeps the old path: an incomplete NH3
  # (Lys with HZ3 missing) is not parameterized, PUT stays planar
  m = model(pdb_lys_no_hz3, riding_h=False)
  atoms = dict((a.name.strip(), a) for a in m.get_hierarchy().atoms())
  for srv in (None, m.get_mon_lib_srv()):
    rm = riding.manager(pdb_hierarchy=m.get_hierarchy(),
      geometry_restraints=m.get_restraints_manager().geometry,
      mon_lib_srv=srv)
    for h in ('HZ1', 'HZ2'):
      assert rm.h_parameterization[atoms[h].i_seq] is None, (srv, h)
  m = model(pdb_put, riding_h=False)
  rm = riding.manager(pdb_hierarchy=m.get_hierarchy(),
    geometry_restraints=m.get_restraints_manager().geometry, mon_lib_srv=None)
  sc = m.get_sites_cart().deep_copy()
  rm.idealize_riding_h_positions(sites_cart=sc)
  m.set_sites_cart(sc)
  m.riding_h_manager = rm
  hnh, xnh, height, types = nh2(m, 'N1', 'C1', ('HN11', 'HN12'))
  assert types == ['alg1a', 'alg1a'], types
  assert height < 0.02, height

def exercise_aniline_no_plane():
  # a conjugated X (aromatic C) keeps the planar pair whatever the angles
  m = model(pdb_anl, cif_anl)
  hnh, xnh, height, types = nh2(m, 'N', 'C1', ('HN1', 'HN2'))
  assert types == ['alg1a', 'alg1a'], types
  assert height < 0.02, height

def exercise_plane_any_h():
  # an NH2 H in a plane restraint keeps master's path (here: H-N-H 109.5,
  # not parameterized), also when another H comes first in the plane
  m = model(pdb_lg5, cif_lg5)
  hnh, xnh, height, types = nh2(m, 'N', 'CA', ('HN1', 'HN2'))
  assert types == [None, None], types

def exercise_infeasible_angles():
  # angle sum above 360: no pyramidal NH2 gives it; master's path (planar
  # pair, rides also in neutron refinement)
  r = math.radians
  assert parameterization.manager.nh2_dphi(r(125), r(125), r(120)) \
    is None
  assert parameterization.manager.nh2_dphi(r(100), r(130), r(25)) \
    is None
  assert abs(parameterization.manager.nh2_dphi(r(120), r(120),
    r(120)) - math.pi) < 1.e-6
  for kw in ({}, {'ignore_h_with_dof': True}):
    m = model(pdb_lg7, cif_lg7, **kw)
    hnh, xnh, height, types = nh2(m, 'N', 'CA', ('HN1', 'HN2'))
    assert types == ['alg1a', 'alg1a'], (kw, types)

def _seeded_runs(pdb_str, cif_str, seeded):
  """placement with the seeded H at one site (the heavy atoms' mean + seed),
  input rotated rigidly: yields (model, rotation, center)"""
  rots = [matrix.sqr((1, 0, 0, 0, 1, 0, 0, 0, 1))] + [
    matrix.col(axis).normalize().axis_and_angle_as_r3_rotation_matrix(
      angle, deg=True) for axis, angle in (((1, 2, 3), 70), ((-2, 1, 0.5), 155))]
  inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  h0 = inp.construct_hierarchy()
  center = matrix.col(h0.atoms().extract_xyz().mean())
  for rot in rots:
    for seed in ((0.5, 0.5, 0.5), (-0.4, 0.9, 0.2)):
      h = h0.deep_copy()
      for a in h.atoms():
        a.xyz = tuple(rot * (matrix.col(a.xyz) - center) + center)
      heavy = flex.vec3_double([a.xyz for a in h.atoms()
        if a.element.strip() != 'H'])
      for a in h.atoms():
        if a.name.strip() in seeded:
          a.xyz = tuple(matrix.col(heavy.mean()) + matrix.col(seed))
      yield (model(h.as_pdb_string(crystal_symmetry=inp.crystal_symmetry()),
        cif_str, use_ideal_dihedral=True), rot, center)

def exercise_h_third_torsion_frame():
  # placement seeds every H at one site: a frame from an H third atom (still
  # at its seed) made the NH2 depend on the seed. With a heavy atom bonded
  # to X the frame is that atom: identical sites by name
  all_h = ('HN1', 'HN2', 'HA1', 'HA2', 'HB1', 'HB2', 'HO')
  ref = None
  for m, rot, center in _seeded_runs(pdb_lg6, cif_lg6, all_h):
    hnh, xnh, height, types = nh2(m, 'N', 'CA', ('HN1', 'HN2'))
    assert types == ['alg1a', 'alg1a'], types
    assert abs(hnh - 106.0) < 0.1, hnh
    s = _put_sites(m)
    # the first H anti to the frame atom C (no heavy-third torsion)
    from scitbx.math import dihedral_angle
    anti = [abs(abs(dihedral_angle(sites=[s['C'], s['CA'], s['N'], s[n]],
      deg=True)) - 180) for n in ('HN1', 'HN2')]
    assert min(anti) < 0.5, anti
    back = dict((n, rot.transpose() * (s[n] - center) + center)
      for n in ('HN1', 'HN2'))
    if ref is None:
      ref = back
      continue
    for n in back:
      assert (back[n] - ref[n]).length() < 0.005, (n,
        (back[n] - ref[n]).length())
  # methylamine: no heavy atom bonded to X, the frame is an H (at its seed
  # in placement): pyramidal at the dictionary angles, first H anti to the
  # frame H; the orientation is arbitrary (not compared by name). Before:
  # master's planar pair (H-N-H 141)
  for m, rot, center in _seeded_runs(pdb_lg4, cif_lg4,
      ('HN1', 'HN2', 'H1', 'H2', 'H3')):
    hnh, xnh, height, types = nh2(m, 'N', 'C', ('HN1', 'HN2'))
    assert types == ['alg1a', 'alg1a'], types
    assert abs(hnh - 106.0) < 0.1 and height > 0.3, (height, hnh)
    for v in xnh:
      assert abs(v - 109.5) < 0.1, xnh
  m = model(pdb_lg4, cif_lg4)
  hnh, xnh, height, types = nh2(m, 'N', 'C', ('HN1', 'HN2'))
  assert types == ['alg1a', 'alg1a'] and abs(hnh - 106.0) < 0.1, (types, hnh)

if __name__ == '__main__':
  exercise_geostd_put()
  exercise_user_cif_109()
  exercise_asn_planar()
  exercise_unequal_restraints()
  exercise_alt_conformers()
  exercise_repeated_idealization()
  exercise_seed_frame_invariance()
  exercise_ignore_h_with_dof()
  exercise_no_cif()
  exercise_aniline_no_plane()
  exercise_plane_any_h()
  exercise_infeasible_angles()
  exercise_h_third_torsion_frame()
  print('OK')
