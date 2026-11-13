"""reduce2: an H on a blank parent whose covalent partner in another residue is
split goes into one atom group per letter of the partner, so riding builds each
copy from that conformer (A26). Occupancy from the partner's conformer, never
above the parent's.
Fragments: 1us0 (Leu A 15, split, Gly A 16), 7T5E (Ser A 102 C alt A/B 1.5 A
apart, Gly A 103), 6vw1 (Asn A 90 and its N-glycan, ND2 split by hand), 1us0
Lys A 100 / Leu A 101 with other letters on Leu 101 (B7/B8), 1xmk (Ni site).
"""
from __future__ import absolute_import, division, print_function
import os
import iotbx.pdb
import mmtbx.model
import mmtbx.hydrogens.riding
from libtbx.utils import null_out
from six.moves import cStringIO as StringIO
from iotbx.cli_parser import run_program
from mmtbx.programs import reduce2
from mmtbx.hydrogens import reduce_hydrogen
from mmtbx_hydrogens_ext import apply_new_H_positions
from scitbx.matrix import col

pdb_1us0 = """\
CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1
ATOM    114  N  ALEU A  15      14.664   4.195  33.948  0.65  3.72           N
ATOM    115  CA ALEU A  15      15.410   2.941  33.949  0.65  3.38           C
ATOM    116  C  ALEU A  15      14.371   1.839  33.819  0.65  3.80           C
ATOM    117  O  ALEU A  15      13.415   1.754  34.585  0.65  6.14           O
ATOM    118  CB ALEU A  15      16.182   2.809  35.268  0.65  4.04           C
ATOM    119  CG ALEU A  15      16.961   1.505  35.381  0.65  5.05           C
ATOM    120  CD1ALEU A  15      17.953   1.319  34.251  0.65  5.38           C
ATOM    121  CD2ALEU A  15      17.733   1.535  36.697  0.65 10.96           C
ATOM    122  N  BLEU A  15      14.681   4.161  34.207  0.35  3.71           N
ATOM    123  CA BLEU A  15      15.248   2.820  34.194  0.35  4.46           C
ATOM    124  C  BLEU A  15      14.217   1.788  33.773  0.35  4.08           C
ATOM    125  O  BLEU A  15      13.163   1.756  34.386  0.35  4.24           O
ATOM    126  CB BLEU A  15      15.814   2.293  35.555  0.35  5.69           C
ATOM    127  CG BLEU A  15      17.326   2.315  35.630  0.35  8.88           C
ATOM    128  CD1BLEU A  15      17.933   1.215  34.800  0.35  9.92           C
ATOM    129  CD2BLEU A  15      17.824   3.713  35.303  0.35 10.13           C
ATOM    130  N   GLY A  16      14.579   0.964  32.822  1.00  3.49           N
ATOM    131  CA  GLY A  16      13.750  -0.200  32.616  1.00  3.60           C
ATOM    132  C   GLY A  16      14.629  -1.435  32.425  1.00  3.24           C
ATOM    133  O   GLY A  16      15.820  -1.362  32.229  1.00  4.16           O
ATOM    134  N   LEU A  17      13.937  -2.592  32.470  1.00  3.21           N
ATOM    135  CA  LEU A  17      14.588  -3.873  32.248  1.00  2.99           C
ATOM    136  C   LEU A  17      14.376  -4.319  30.802  1.00  2.76           C
ATOM    137  O   LEU A  17      13.234  -4.478  30.352  1.00  3.50           O
ATOM    138  CB  LEU A  17      13.979  -4.930  33.171  1.00  3.46           C
ATOM    139  CG  LEU A  17      14.593  -6.328  33.031  1.00  3.74           C
ATOM    140  CD1 LEU A  17      15.997  -6.376  33.577  1.00  5.85           C
ATOM    141  CD2 LEU A  17      13.702  -7.359  33.696  1.00  4.98           C
END
"""

pdb_7t5e = """\
CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1
ATOM    750  N   ALA A 101     -22.484 -23.151 -27.063  1.00 21.54           N
ATOM    751  CA  ALA A 101     -21.115 -23.421 -26.638  1.00 27.76           C
ATOM    752  C   ALA A 101     -20.637 -22.484 -25.513  1.00 23.05           C
ATOM    753  O   ALA A 101     -20.935 -21.283 -25.496  1.00 19.90           O
ATOM    754  CB  ALA A 101     -20.202 -23.319 -27.853  1.00 18.19           C
ATOM    755  N  ASER A 102     -19.771 -22.934 -24.587  0.61 25.21           N
ATOM    756  CA ASER A 102     -19.503 -22.013 -23.445  0.61 25.64           C
ATOM    757  C  ASER A 102     -18.406 -20.956 -23.680  0.61 15.74           C
ATOM    758  O  ASER A 102     -18.427 -19.911 -23.041  0.61 21.72           O
ATOM    759  CB ASER A 102     -19.139 -22.762 -22.161  0.61 21.34           C
ATOM    760  OG ASER A 102     -17.928 -23.487 -22.297  0.61 25.13           O
ATOM    761  N  BSER A 102     -20.064 -23.173 -24.539  0.39 24.66           N
ATOM    762  CA BSER A 102     -19.033 -22.685 -23.659  0.39 23.83           C
ATOM    763  C  BSER A 102     -17.776 -22.423 -24.468  0.39  9.65           C
ATOM    764  O  BSER A 102     -17.128 -23.330 -25.003  0.39 22.38           O
ATOM    765  CB BSER A 102     -18.761 -23.701 -22.560  0.39 28.21           C
ATOM    766  OG BSER A 102     -18.282 -23.032 -21.420  0.39 23.61           O
ATOM    767  N   GLY A 103     -17.521 -21.161 -24.626  1.00 20.91           N
ATOM    768  CA  GLY A 103     -16.243 -20.493 -24.729  1.00 20.77           C
ATOM    769  C   GLY A 103     -15.085 -20.969 -23.887  1.00 18.14           C
ATOM    770  O   GLY A 103     -13.935 -20.719 -24.253  1.00 17.90           O
ATOM    771  N   VAL A 104     -15.357 -21.633 -22.765  1.00 16.19           N
ATOM    772  CA  VAL A 104     -14.306 -22.058 -21.850  1.00 20.76           C
ATOM    773  C   VAL A 104     -13.507 -23.193 -22.471  1.00 21.56           C
ATOM    774  O   VAL A 104     -14.073 -24.196 -22.915  1.00 24.16           O
ATOM    775  CB  VAL A 104     -14.907 -22.496 -20.509  1.00 21.71           C
ATOM    776  CG1 VAL A 104     -13.813 -23.029 -19.561  1.00 26.01           C
ATOM    777  CG2 VAL A 104     -15.605 -21.320 -19.887  1.00 21.63           C
END
"""

pdb_glycan = """\
CRYST1   80.000   80.000   80.000  90.00  90.00  90.00 P 1
ATOM      1  N   ASN A  90      76.017 -20.785 155.731  1.00142.97           N
ATOM      1  CA  ASN A  90      75.942 -19.415 156.227  1.00148.71           C
ATOM      1  C   ASN A  90      77.350 -18.831 156.235  1.00143.80           C
ATOM      1  O   ASN A  90      78.239 -19.359 156.912  1.00146.03           O
ATOM      1  CB  ASN A  90      75.320 -19.390 157.626  1.00158.76           C
ATOM      1  CG AASN A  90      75.068 -17.982 158.153  0.60168.44           C
ATOM      1  CG BASN A  90      74.890 -17.998 158.074  0.40168.44           C
ATOM      1  OD1AASN A  90      75.922 -17.100 158.061  0.60161.90           O
ATOM      1  OD1BASN A  90      75.625 -17.021 157.928  0.40161.90           O
ATOM      1  ND2AASN A  90      73.886 -17.776 158.730  0.60187.09           N
ATOM      1  ND2BASN A  90      73.691 -17.911 158.646  0.40187.09           N
HETATM    1  C1  NAG A   1      73.523 -16.504 159.309  1.00103.52           C
HETATM    1  C2  NAG A   1      72.427 -15.749 158.564  1.00114.27           C
HETATM    1  C3  NAG A   1      72.112 -14.434 159.280  1.00117.99           C
HETATM    1  C4  NAG A   1      71.829 -14.658 160.763  1.00124.42           C
HETATM    1  C5  NAG A   1      72.930 -15.509 161.403  1.00121.57           C
HETATM    1  C6  NAG A   1      72.619 -15.937 162.820  1.00122.81           C
HETATM    1  C7  NAG A   1      71.992 -15.651 156.147  1.00125.03           C
HETATM    1  C8  NAG A   1      72.570 -15.344 154.798  1.00127.15           C
HETATM    1  N2  NAG A   1      72.819 -15.497 157.186  1.00121.51           N
HETATM    1  O3  NAG A   1      70.990 -13.814 158.662  1.00116.89           O
HETATM    1  O4  NAG A   1      71.791 -13.386 161.404  1.00133.35           O
HETATM    1  O5  NAG A   1      73.123 -16.717 160.654  1.00114.93           O
HETATM    1  O6  NAG A   1      73.745 -15.787 163.673  1.00121.65           O
HETATM    1  O7  NAG A   1      70.831 -16.020 156.289  1.00125.47           O
END
"""

pdb_letters = """\
CRYST1   49.280   66.590   47.260  90.00  92.40  90.00 P 1
ATOM      1  CB  LEU A  99       2.042  -8.443  28.374  1.00  4.39           C
ATOM      1  CG  LEU A  99       3.288  -8.366  27.496  1.00  4.12           C
ATOM      1  CD1 LEU A  99       3.490  -6.934  27.006  1.00  5.01           C
ATOM      1  CD2 LEU A  99       4.542  -8.832  28.221  1.00  4.87           C
ATOM      1  N  ALEU A  99       1.157 -10.590  27.621  0.66  4.43           N
ATOM      1  CA ALEU A  99       1.636  -9.865  28.775  0.66  4.35           C
ATOM      1  C  ALEU A  99       0.647  -9.804  29.954  0.66  5.10           C
ATOM      1  O  ALEU A  99       1.001  -9.301  31.013  0.66  5.75           O
ATOM      1  N  BLEU A  99       1.128 -10.561  27.518  0.34  4.21           N
ATOM      1  CA BLEU A  99       1.589  -9.867  28.703  0.34  4.68           C
ATOM      1  C  BLEU A  99       0.511  -9.920  29.786  0.34  5.07           C
ATOM      1  O  BLEU A  99       0.826  -9.609  30.921  0.34  6.87           O
ATOM      1  N  ALYS A 100      -0.576 -10.322  29.747  0.66  5.43           N
ATOM      1  CA ALYS A 100      -1.619 -10.320  30.778  0.66  5.43           C
ATOM      1  C  ALYS A 100      -2.069  -8.897  31.115  0.66  5.93           C
ATOM      1  O  ALYS A 100      -2.394  -8.582  32.262  0.66  8.82           O
ATOM      1  CB ALYS A 100      -1.190 -11.097  32.066  0.66  6.55           C
ATOM      1  CG ALYS A 100      -1.608 -12.556  32.067  0.66  6.01           C
ATOM      1  CD ALYS A 100      -0.951 -13.432  31.019  0.66  4.90           C
ATOM      1  CE ALYS A 100       0.545 -13.638  31.353  0.66  4.24           C
ATOM      1  NZ ALYS A 100       1.166 -14.539  30.345  0.66  3.46           N
ATOM      1  N  BLYS A 100      -0.704 -10.316  29.364  0.34  5.05           N
ATOM      1  CA BLYS A 100      -1.886 -10.422  30.230  0.34  4.48           C
ATOM      1  C  BLYS A 100      -2.280  -9.075  30.819  0.34  4.68           C
ATOM      1  O  BLYS A 100      -2.690  -8.931  31.957  0.34  7.12           O
ATOM      1  CB BLYS A 100      -1.678 -11.489  31.328  0.34  4.37           C
ATOM      1  CG BLYS A 100      -1.577 -12.883  30.725  0.34  3.94           C
ATOM      1  CD BLYS A 100      -0.951 -13.867  31.731  0.34  3.67           C
ATOM      1  CE BLYS A 100       0.563 -13.923  31.648  0.34  4.64           C
ATOM      1  NZ BLYS A 100       1.046 -14.494  30.363  0.34  6.20           N
ATOM      1  N   LEU A 101      -2.148  -8.042  30.050  1.00  5.54           N
ATOM      1  CA  LEU A 101      -2.483  -6.661  30.271  1.00  5.20           C
ATOM      1  C   LEU A 101      -3.724  -6.281  29.497  1.00  4.90           C
ATOM      1  O   LEU A 101      -4.048  -6.860  28.467  1.00  5.82           O
ATOM      1  CB CLEU A 101      -1.327  -5.767  29.810  0.60  7.02           C
ATOM      1  CB DLEU A 101      -1.177  -5.767  29.810  0.40  7.02           C
ATOM      1  CG CLEU A 101      -0.031  -6.041  30.578  0.60 11.07           C
ATOM      1  CG DLEU A 101       0.119  -6.041  30.578  0.40 11.07           C
ATOM      1  CD1CLEU A 101       1.079  -5.212  29.980  0.60 16.98           C
ATOM      1  CD1DLEU A 101       1.229  -5.212  29.980  0.40 16.98           C
ATOM      1  CD2CLEU A 101      -0.164  -5.762  32.061  0.60 15.86           C
ATOM      1  CD2DLEU A 101      -0.014  -5.762  32.061  0.40 15.86           C
END
"""

# 1xmk: Gly 288 blank, Ser 289 N blank with the rest in A/B, His 290 N on the Ni
pdb_ni = """\
CRYST1   35.557   43.526   45.471  90.00  90.00  90.00 P 21 21 21
ATOM      1  N   GLY A 288      12.334  23.080  30.892  1.00 20.00           N
ATOM      2  CA  GLY A 288      11.259  22.629  29.968  1.00 20.00           C
ATOM      3  C   GLY A 288      11.734  21.310  29.346  1.00 20.00           C
ATOM      4  O   GLY A 288      11.046  20.536  28.686  1.00 20.00           O
ATOM      5  N   SER A 289      12.992  21.082  29.673  1.00 20.00           N
ATOM      6  CA ASER A 289      13.766  20.002  29.117  0.65 20.00           C
ATOM      7  C  ASER A 289      14.970  19.728  30.022  0.65 20.00           C
ATOM      8  O  ASER A 289      15.660  18.700  29.826  0.65 20.00           O
ATOM      9  CB ASER A 289      14.288  20.329  27.695  0.65 20.00           C
ATOM     10  OG ASER A 289      14.983  21.589  27.736  0.65 20.00           O
ATOM     11  CA BSER A 289      13.766  20.002  29.117  0.35 20.00           C
ATOM     12  C  BSER A 289      14.900  19.628  30.022  0.35 20.00           C
ATOM     13  O  BSER A 289      15.660  18.700  29.826  0.35 20.00           O
ATOM     14  CB BSER A 289      14.388  20.129  27.695  0.35 20.00           C
ATOM     15  OG BSER A 289      13.583  21.089  27.036  0.35 20.00           O
ATOM     16  N   HIS A 290      15.185  20.681  30.963  1.00 20.00           N
ATOM     17  CA  HIS A 290      16.475  20.629  31.642  1.00 20.00           C
ATOM     18  C   HIS A 290      16.332  19.771  32.904  1.00 20.00           C
ATOM     19  O   HIS A 290      16.691  20.326  33.953  1.00 20.00           O
ATOM     20  CB  HIS A 290      17.030  22.019  31.962  1.00 20.00           C
ATOM     21  CG  HIS A 290      16.054  22.940  32.675  1.00 20.00           C
ATOM     22  ND1 HIS A 290      14.740  23.013  32.337  1.00 20.00           N
ATOM     23  CD2 HIS A 290      16.315  23.854  33.689  1.00 20.00           C
ATOM     24  CE1 HIS A 290      14.147  23.928  33.126  1.00 20.00           C
ATOM     25  NE2 HIS A 290      15.081  24.445  33.946  1.00 20.00           N
TER
HETATM   26 NI    NI A 398      13.849  21.950  30.990  NIOC 20.00          NI
END
"""

CONFORMER_SITES = 'Conformer sites written as one isotope'
NOT_PLACED = ('N-H not placed: backbone amide N or thymine/uracil N3 bound to '
  'a metal, taken as deprotonated')
UNRESOLVED = 'N-H placed at metal sites with unresolved protonation'
LETTERS = ("H placed once, not split: the residue's altloc letters differ "
  "from its covalent partner's")

def run(prefix, pdb_str, args=()):
  '''reduce2 on pdb_str. :return: (output hierarchy, description, log).'''
  fn, out, txt = prefix + '_in.pdb', prefix + '_out.pdb', prefix + '.txt'
  with open(fn, 'w') as f:
    f.write(pdb_str)
  for n in (out, txt):
    if os.path.exists(n): os.remove(n)
  log = StringIO()
  run_program(program_class=reduce2.Program, logger=log, args=[fn,
    'output.filename=' + out, 'output.description_file_name=' + txt]
    + list(args))
  with open(txt) as f:
    description = f.read()
  return iotbx.pdb.input(out).construct_hierarchy(), description, log.getvalue()

def rerun(prefix, h_in, args):
  '''reduce2 on a written hierarchy (its file), for reruns.'''
  return run(prefix, h_in.as_pdb_string(), args)

def atoms_of(hierarchy, resseq, names):
  '''{(name, altloc): atom} in residue resseq (a string).'''
  return dict(((a.name.strip(), a.parent().altloc), a)
    for rg in hierarchy.residue_groups() if rg.resseq.strip() == resseq
    for a in rg.atoms() if a.name.strip() in names)

def conformer_lines(pdb_str, letter):
  '''pdb_str restricted to blank + letter.'''
  return ''.join(l + '\n' for l in pdb_str.splitlines()
    if not (l[:6] in ('ATOM  ', 'HETATM') and l[16] not in (' ', letter)))

def blank_letters(pdb_str, letter):
  '''pdb_str restricted to blank + letter, that letter written blank, full
  occupancy.'''
  result = ''
  for l in conformer_lines(pdb_str, letter).splitlines():
    if l[:6] in ('ATOM  ', 'HETATM') and l[16] == letter:
      l = l[:16] + ' ' + l[17:54] + '  1.00' + l[60:]
    result += l + '\n'
  return result

def fresh_riding_distances(hierarchy, names, neutron=False):
  '''Written H/D with names -> distance to the position a fresh riding manager
  gives after reading the file, as refinement does.'''
  m = mmtbx.model.manager(
    model_input=iotbx.pdb.input(lines=hierarchy.as_pdb_string().splitlines(),
      source_info=None),
    stop_for_unknowns=False, log=null_out())
  reduce_hydrogen.add_ccd_restraint_objects(m)
  p = reduce_hydrogen.get_reduce_pdb_interpretation_params(neutron)
  p.pdb_interpretation.flip_symmetric_amino_acids = False
  m.process(make_restraints=True, pdb_interpretation_params=p,
    plain_pairs_radius=None)
  h = m.get_hierarchy()
  x0 = m.get_sites_cart()
  rm = mmtbx.hydrogens.riding.manager(pdb_hierarchy=h,
    geometry_restraints=m.get_restraints_manager().geometry,
    use_ideal_dihedral=False, mon_lib_srv=m.get_mon_lib_srv())
  x1 = x0.deep_copy()
  apply_new_H_positions(sites_cart=x1, parameterization=rm.parameterization_cpp)
  atoms = h.atoms()
  result = {}
  for rc in rm.parameterization_cpp:
    a = atoms[rc.ih]
    if a.name.strip() in names:
      result[(a.parent().parent().resseq.strip(), a.name.strip(),
        a.parent().altloc)] = (col(x0[rc.ih]) - col(x1[rc.ih])).length()
  return result

def test_amide_after_split_c():
  h, d, l = run('tst_a26_amide', pdb_1us0)
  a = atoms_of(h, '16', ('H',))
  assert sorted(a) == [('H', 'A'), ('H', 'B')], sorted(a)
  assert abs(a[('H', 'A')].occ - 0.65) < 1.e-3 and \
         abs(a[('H', 'B')].occ - 0.35) < 1.e-3
  # Leu 17 follows a blank C: one blank H
  assert sorted(atoms_of(h, '17', ('H',))) == [('H', '')]
  # each copy where the model of its conformer alone puts it
  for x in 'AB':
    w = atoms_of(run('tst_a26_amide_' + x, conformer_lines(pdb_1us0, x))[0],
      '16', ('H',))
    assert a[('H', x)].distance(w[('H', '')]) < 0.01, x
  # the two conformers really differ here
  assert a[('H', 'A')].distance(a[('H', 'B')]) > 0.05
  for prefix, args in (('tst_a26_amide_opt', ['approach=optimize']),
      ('tst_a26_amide_keep', ['approach=add', 'keep_existing_H=True'])):
    b = atoms_of(rerun(prefix, h, args)[0], '16', ('H',))
    assert sorted(b) == sorted(a), (prefix, sorted(b))
    for k in a:
      assert a[k].distance(b[k]) < 0.01, (prefix, k)
      assert abs(a[k].occ - b[k].occ) < 1.e-3, (prefix, k, b[k].occ)

def test_7t5e_gly_103_hd_and_h():
  h, d, l = run('tst_a26_7t5e', pdb_7t5e,
    ['neutron_option=hd_and_h', 'use_neutron_distances=True'])
  a = atoms_of(h, '103', ('H', 'D'))
  assert sorted(a) == [('D', 'A'), ('D', 'B')], sorted(a)
  assert abs(a[('D', 'A')].occ - 0.61) < 1.e-3 and \
         abs(a[('D', 'B')].occ - 0.39) < 1.e-3
  assert a[('D', 'A')].distance(a[('D', 'B')]) > 1.0
  for k, dist in fresh_riding_distances(h, ('D',), neutron=True).items():
    if k[0] == '103': assert dist < 0.01, (k, dist)
  # listed as sites without an H/D occupancy pair, no occupancy notes for them
  assert CONFORMER_SITES in d, d
  assert 'no H/D occupancy pair' in d, d
  notes = d[d.index('Occupancy notes'):] if 'Occupancy notes' in d else ''
  assert 'GLY A 103' not in notes, notes

def test_glycan_h1():
  # A20 limit: H1 has no angle to ND2 (2neigbs on O5, C2), so a fresh riding
  # manager moves it 0.05 A also without any split: positions are compared
  # with the unsplit control, not with fresh riding.
  h, d, l = run('tst_a26_glycan', pdb_glycan)
  a = atoms_of(h, '1', ('H1', 'C1'))
  nd2 = atoms_of(h, '90', ('ND2',))
  assert sorted(k[1] for k in a if k[0] == 'H1') == ['A', 'B'], sorted(a)
  assert abs(a[('H1', 'A')].occ - 0.6) < 1.e-3
  assert abs(a[('H1', 'B')].occ - 0.4) < 1.e-3
  c1 = col(a[('C1', '')].xyz)
  for x in 'AB':
    angle = (col(a[('H1', x)].xyz) - c1).angle(col(nd2[('ND2', x)].xyz) - c1,
      deg=True)
    assert 95 < angle < 120, (x, angle)
  # unsplit control: one H1, on the site of the copy built from that conformer
  for x in 'AB':
    c = atoms_of(run('tst_a26_glycan_ctl' + x, blank_letters(pdb_glycan, x))[0],
      '1', ('H1',))
    assert sorted(c) == [('H1', '')]
    assert c[('H1', '')].distance(a[('H1', x)]) < 0.01, x

def test_glycan_contact_at_other_atom():
  # NAG C2 in contact with ND2 A/B, 1.80 A: the NAG-ASN link is C1-ND2, so H2
  # is one blank H at the parent's occupancy; H1 still split
  lines = []
  for l in pdb_glycan.splitlines():
    if l[12:16] == ' C2 ' and l[17:20] == 'NAG':
      l = l[:30] + '%8.3f%8.3f%8.3f' % (72.836775, -16.318304, 158.610623) + \
        l[54:]
    lines.append(l)
  h = run('tst_a26_glycan_c2', '\n'.join(lines) + '\n')[0]
  a = atoms_of(h, '1', ('H1', 'H2'))
  assert sorted(k for k in a if k[0] == 'H2') == [('H2', '')], sorted(a)
  assert abs(a[('H2', '')].occ - 1) < 1.e-3
  assert sorted(k[1] for k in a if k[0] == 'H1') == ['A', 'B'], sorted(a)
  assert abs(a[('H1', 'A')].occ - 0.6) < 1.e-3
  assert abs(a[('H1', 'B')].occ - 0.4) < 1.e-3

def test_validate_ligands_one_blank_glycan():
  # the H1 copies make no conformer of the NAG (fails without the change in
  # validate_ligands.py)
  from mmtbx.programs import validate_ligands
  with open('tst_a26_val.pdb', 'w') as f:
    f.write(pdb_glycan)
  result = run_program(program_class=validate_ligands.Program,
    args=['tst_a26_val.pdb'], logger=null_out())
  found = sorted((lr.resname, lr.altloc) for lr in result.ligand_manager)
  assert found == [('NAG', '')], found

def water_lines(blank=False):
  # O 1.4 and 1.5 A from Gly 16 N, along N-H; no bond to anything
  return ('HETATM  301  O  AHOH A 201      15.614   1.083  31.880  0.10  5.00'
    '           O\n'
    'HETATM  302  O  BHOH A 201      15.688   1.091  31.813  0.20  5.00'
    '           O\n')

def test_water_is_no_partner():
  # with a split C: letters and occupancies from the C (fails on master: blank H)
  h = run('tst_a26_water', pdb_1us0.replace('END\n', water_lines() + 'END\n'))[0]
  a = atoms_of(h, '16', ('H',))
  assert sorted(a) == [('H', 'A'), ('H', 'B')], sorted(a)
  assert abs(a[('H', 'A')].occ - 0.65) < 1.e-3
  assert abs(a[('H', 'B')].occ - 0.35) < 1.e-3

def test_water_no_split_of_its_own():
  # all-blank peptide: one blank H at the parent's occupancy (passes on master)
  blank = blank_letters(pdb_1us0, 'A')
  h = run('tst_a26_water_blank', blank.replace('END\n',
    water_lines() + 'END\n'))[0]
  a = atoms_of(h, '16', ('H',))
  assert sorted(a) == [('H', '')], sorted(a)
  assert abs(a[('H', '')].occ - 1) < 1.e-3

def test_letters_differ_unchanged():
  # Leu 101 side chain in C/D, Lys 100 C in A/B: as on master (B7/B8)
  h, d, l = run('tst_a26_letters', pdb_letters)
  a = atoms_of(h, '101', ('H',))
  assert sorted(a) == [('H', 'C'), ('H', 'D')], sorted(a)
  assert abs(a[('H', 'C')].occ - 0.6) < 1.e-3
  assert abs(a[('H', 'D')].occ - 0.4) < 1.e-3
  return d, l

def test_letters_differ_logged():
  # the full line in the description, a summary in the log
  d, l = test_letters_differ_unchanged()
  assert [x for x in d.splitlines() if x.startswith(LETTERS)], d
  line = [x for x in d.splitlines() if 'LEU A 101' in x and 'A, B' in x]
  assert len(line) == 1, line
  assert 'placed, not split per conformer: 1 H, 1 residue: LEU A 101' \
    in l, l

pdb_1pq7 = """\
CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1
ATOM    328  N   TRP A  40       6.387  10.256   1.511  1.00  5.49           N
ATOM    329  CA  TRP A  40       7.193   9.080   1.247  1.00  5.04           C
ATOM    330  C   TRP A  40       6.321   7.884   0.831  1.00  4.88           C
ATOM    331  O   TRP A  40       6.480   7.354  -0.277  1.00  5.86           O
ATOM    332  CB  TRP A  40       8.095   8.763   2.458  1.00  5.03           C
ATOM    333  CG  TRP A  40       9.027   7.641   2.206  1.00  5.03           C
ATOM    334  CD1 TRP A  40       8.993   6.404   2.791  1.00  5.31           C
ATOM    335  CD2 TRP A  40      10.125   7.632   1.305  1.00  5.34           C
ATOM    336  NE1 TRP A  40      10.003   5.631   2.317  1.00  6.10           N
ATOM    337  CE2 TRP A  40      10.722   6.354   1.386  1.00  5.95           C
ATOM    338  CE3 TRP A  40      10.668   8.596   0.434  1.00  7.16           C
ATOM    339  CZ2 TRP A  40      11.829   6.034   0.622  1.00  7.79           C
ATOM    340  CZ3 TRP A  40      11.776   8.226  -0.326  1.00  9.41           C
ATOM    341  CH2 TRP A  40      12.325   6.985  -0.228  1.00  9.46           C
ATOM    352  N  ACYS A  41       5.508   7.360   1.752  0.59  3.97           N
ATOM    353  N  BCYS A  41       5.374   7.503   1.660  0.40  4.85           N
ATOM    354  CA ACYS A  41       4.721   6.164   1.517  0.60  3.70           C
ATOM    355  CA BCYS A  41       4.595   6.308   1.311  0.39  5.11           C
ATOM    356  C  ACYS A  41       3.333   6.308   2.146  0.60  3.80           C
ATOM    357  C  BCYS A  41       3.340   6.305   2.209  0.40  4.25           C
ATOM    358  O  ACYS A  41       2.979   7.330   2.757  0.60  3.74           O
ATOM    359  O  BCYS A  41       3.124   7.213   3.049  0.40  3.20           O
ATOM    360  CB ACYS A  41       5.431   4.898   2.010  0.61  3.50           C
ATOM    361  CB BCYS A  41       5.483   5.140   1.450  0.40  7.12           C
ATOM    362  SG ACYS A  41       6.731   4.335   0.827  0.60  3.82           S
ATOM    363  SG BCYS A  41       5.836   4.613   3.064  0.39 10.10           S
ATOM    372  N   GLY A  42       2.514   5.271   1.971  1.00  4.10           N
ATOM    373  CA  GLY A  42       1.322   5.042   2.727  1.00  4.02           C
ATOM    374  C   GLY A  42       1.417   3.729   3.505  1.00  3.74           C
ATOM    375  O   GLY A  42       2.392   2.983   3.382  1.00  4.06           O
ATOM    570  N   HIS A  56       9.728  -2.142   0.726  1.00  4.02           N
ATOM    571  CA  HIS A  56      10.687  -1.249   0.118  1.00  4.42           C
ATOM    572  C   HIS A  56      10.576   0.171   0.697  1.00  4.38           C
ATOM    573  O   HIS A  56      11.534   0.956   0.573  1.00  5.25           O
ATOM    574  CB  HIS A  56      10.592  -1.249  -1.407  1.00  4.99           C
ATOM    575  CG  HIS A  56       9.480  -0.462  -1.999  1.00  5.02           C
ATOM    576  ND1 HIS A  56       8.294  -0.988  -2.415  1.00  5.14           N
ATOM    577  CD2 HIS A  56       9.435   0.875  -2.263  1.00  5.76           C
ATOM    578  CE1 HIS A  56       7.560  -0.010  -2.923  1.00  5.72           C
ATOM    579  NE2 HIS A  56       8.216   1.121  -2.838  1.00  5.92           N
ATOM    587  N   CYS A  57       9.445   0.516   1.315  1.00  4.44           N
ATOM    588  CA  CYS A  57       9.313   1.815   1.946  1.00  4.72           C
ATOM    589  C   CYS A  57      10.241   1.988   3.148  1.00  4.55           C
ATOM    590  O   CYS A  57      10.553   3.132   3.509  1.00  5.19           O
ATOM    591  CB  CYS A  57       7.875   2.020   2.439  1.00  5.19           C
ATOM    592  SG  CYS A  57       6.674   2.294   1.132  1.00  6.27           S
ATOM    597  N   VAL A  58      10.599   0.888   3.802  1.00  4.24           N
ATOM    598  CA  VAL A  58      11.264   0.950   5.096  1.00  4.21           C
ATOM    599  C   VAL A  58      12.556   0.168   5.196  1.00  4.23           C
ATOM    600  O   VAL A  58      13.370   0.456   6.092  1.00  5.18           O
ATOM    601  CB  VAL A  58      10.295   0.532   6.237  1.00  4.70           C
ATOM    602  CG1 VAL A  58       9.187   1.577   6.396  1.00  5.80           C
ATOM    603  CG2 VAL A  58       9.702  -0.860   6.006  1.00  5.12           C
END
"""

def test_disulfide_one_conformer():
  """1pq7: CYS A 41 SG A bonds the blank CYS A 57 SG (2.06 A), SG B is 3.13 A
  away. HG 57 exists in conformer B only, at the occupancy of CYS 41 B; the
  disulfide takes HG A 41 and leaves HG B 41. Deposited H kept: the same."""
  for prefix, args in (('tst_a26_ss', ()),
      ('tst_a26_ss_keep', ('keep_existing_H=True',))):
    h, d, l = run(prefix, pdb_1pq7, args)
    a = atoms_of(h, '57', ('HG',))
    assert sorted(a) == [('HG', 'B')], (prefix, sorted(a))
    assert abs(a[('HG', 'B')].occ - 0.40) < 0.02, a[('HG', 'B')].occ
    assert sorted(atoms_of(h, '41', ('HG',))) == [('HG', 'B')], prefix
    sg = atoms_of(h, '57', ('SG',))[('SG', '')]
    assert abs(a[('HG', 'B')].distance(sg) - 1.2) < 0.2
    # the two conformers are neutral: no Sorry from the charge check
    from mmtbx.ligands import electrons
    m = mmtbx.model.manager(model_input=iotbx.pdb.input(prefix + '_out.pdb'),
      log=null_out())
    m.process(make_restraints=True)
    charges = electrons.electron_distribution(m.get_hierarchy(),
      m.get_restraints_manager().geometry).get_conformer_charges()
    assert charges == {'A': 0, 'B': 0}, charges

def test_metal_per_conformer():
  # Ni at 0.50: covers the conformer at 0.35, not the one at 0.65
  h, d, l = run('tst_a26_ni', pdb_ni.replace('NIOC', '0.50'))
  a = atoms_of(h, '290', ('H',))
  assert sorted(a) == [('H', 'A')], sorted(a)
  assert abs(a[('H', 'A')].occ - 0.65) < 1.e-3
  assert 'H  BHIS A 290' in d and NOT_PLACED in d, d
  assert 'H  AHIS A 290' in d[d.index(UNRESOLVED):], d
  assert 'coverage below population' in d[d.index(UNRESOLVED):] or \
    'population' in d[d.index(UNRESOLVED):], d
  # fully occupied Ni: neither copy
  h, d, l = run('tst_a26_ni_full', pdb_ni.replace('NIOC', '1.00'))
  assert atoms_of(h, '290', ('H',)) == {}
  assert 'H  AHIS A 290' in d and 'H  BHIS A 290' in d, d

if __name__ == '__main__':
  test_amide_after_split_c()
  test_7t5e_gly_103_hd_and_h()
  test_glycan_h1()
  test_glycan_contact_at_other_atom()
  test_validate_ligands_one_blank_glycan()
  test_water_is_no_partner()
  test_water_no_split_of_its_own()
  test_letters_differ_unchanged()
  test_letters_differ_logged()
  test_metal_per_conformer()
  test_disulfide_one_conformer()
  print('OK')
