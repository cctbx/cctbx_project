"""reduce2 metal_bound_amines_as_nh2=True (A5 C1b): a primary amine that the
dictionary makes NH3+ (Lys NZ, an N-terminal NH3 that placement builds) and
that coordinates a metal of NH_METAL_ELEMENTS is placed as NH2, the two H
oriented with the metal in place of the H not placed. Neutral NH2 (Cu(en))
and planar NH2 (Asn) keep both H; supplied H are kept; sites without a unique
metal position (symmetry mate, alternate metal positions) are unresolved.
"""
from __future__ import absolute_import, division, print_function
import os, math
import iotbx.pdb
from scitbx import matrix
from six.moves import cStringIO as StringIO
from iotbx.cli_parser import run_program
from mmtbx.programs import reduce2
from mmtbx.regression.tst_reduce2_metal_nh import pdb_1xmk

# 1xmk Lys 359 NZ on Cd 396 (2.21 A)
pdb_lys_cd = """\
CRYST1   35.557   43.526   45.471  90.00  90.00  90.00 P 21 21 21
ATOM      1  N   LYS A 359       6.646   3.174  41.907  1.00  8.43           N
ATOM      2  CA  LYS A 359       7.297   4.207  42.714  1.00  7.92           C
ATOM      3  C   LYS A 359       7.023   4.045  44.227  1.00  7.33           C
ATOM      4  O   LYS A 359       7.955   4.152  45.067  1.00  7.65           O
ATOM      5  CB  LYS A 359       6.866   5.583  42.191  1.00  8.93           C
ATOM      6  CG  LYS A 359       7.454   6.699  43.061  1.00 11.26           C
ATOM      7  CD  LYS A 359       7.702   7.990  42.363  1.00 16.15           C
ATOM      8  CE  LYS A 359       8.185   9.046  43.369  1.00 18.69           C
ATOM      9  NZ  LYS A 359       9.511   8.759  44.042  1.00 18.51           N
TER
HETATM   10 CD    CD A 396      10.293  10.651  44.885  1.00 10.07          CD
END
"""

# 7T5E His 1 (chain A): N-terminal amine and ND1 on Cu
pdb_his_cu = """\
CRYST1   68.296   42.269   70.413  90.00  98.47  90.00 P 1 21 1
ATOM      1  N   HIS A   1     -14.833   0.543 -16.983  1.00 13.01           N
ATOM      2  CA  HIS A   1     -16.179  -0.029 -16.928  1.00 13.82           C
ATOM      3  C   HIS A   1     -16.832   0.202 -18.278  1.00 16.50           C
ATOM      4  O   HIS A   1     -16.979   1.340 -18.717  1.00 13.17           O
ATOM      5  CB  HIS A   1     -16.979   0.595 -15.795  1.00 17.49           C
ATOM      6  CG  HIS A   1     -16.603   0.040 -14.458  1.00 16.41           C
ATOM      7  ND1 HIS A   1     -15.299   0.015 -14.007  1.00 16.27           N
ATOM      8  CD2 HIS A   1     -17.329  -0.630 -13.547  1.00 14.74           C
ATOM      9  CE1 HIS A   1     -15.261  -0.585 -12.829  1.00 15.40           C
ATOM     10  NE2 HIS A   1     -16.488  -0.965 -12.516  1.00 20.05           N
TER
HETATM   11 CU    CU A 301     -13.829   0.757 -15.122  1.00 19.10          CU
END
"""

# Cu(en): ethylenediamine is a neutral diamine in its dictionary (NH2)
pdb_edn_cu = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
HETATM    1  NAA EDN A   1      14.136  16.021  15.475  1.00 20.00           N
HETATM    2  CAB EDN A   1      15.606  16.021  15.475  1.00 20.00           C
HETATM    3  CAC EDN A   1      16.113  14.588  15.475  1.00 20.00           C
HETATM    4  NAD EDN A   1      15.528  13.860  14.339  1.00 20.00           N
HETATM    5 CU    CU A   2      13.618  14.511  14.237  1.00 20.00          CU
END
"""

# Asn ND2 (1r27 Asn A 52) with a Zn 2.05 A away on the CG-ND2 axis
pdb_asn_zn = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
ATOM      1  N   ASN A   1      14.609  16.377  14.719  1.00 20.00           N
ATOM      2  CA  ASN A   1      15.000  15.000  15.000  1.00 20.00           C
ATOM      3  C   ASN A   1      16.090  14.846  16.057  1.00 20.00           C
ATOM      4  O   ASN A   1      17.253  14.590  15.738  1.00 20.00           O
ATOM      5  CB  ASN A   1      15.424  14.321  13.699  1.00 20.00           C
ATOM      6  CG  ASN A   1      15.765  12.865  13.893  1.00 20.00           C
ATOM      7  OD1 ASN A   1      15.202  12.197  14.757  1.00 20.00           O
ATOM      8  ND2 ASN A   1      16.680  12.357  13.079  1.00 20.00           N
TER
HETATM    9 ZN    ZN B 101      18.095  11.572  11.820  1.00 20.00          ZN
END
"""

# 7T5E His A 1 on Cu A 301, with Asn A 60 and its glycan (NAG C 1, NAG C 2,
# MAN C 3)
pdb_his_cu_glycan = """\
CRYST1   68.296   42.269   70.413  90.00  98.47  90.00 P 1 21 1
ATOM      1  N   HIS A   1     -14.833   0.543 -16.983  1.00 13.01           N
ATOM      2  CA  HIS A   1     -16.179  -0.029 -16.928  1.00 13.82           C
ATOM      3  C   HIS A   1     -16.832   0.202 -18.278  1.00 16.50           C
ATOM      4  O   HIS A   1     -16.979   1.340 -18.717  1.00 13.17           O
ATOM      5  CB  HIS A   1     -16.979   0.595 -15.795  1.00 17.49           C
ATOM      6  CG  HIS A   1     -16.603   0.040 -14.458  1.00 16.41           C
ATOM      7  ND1 HIS A   1     -15.299   0.015 -14.007  1.00 16.27           N
ATOM      8  CD2 HIS A   1     -17.329  -0.630 -13.547  1.00 14.74           C
ATOM      9  CE1 HIS A   1     -15.261  -0.585 -12.829  1.00 15.40           C
ATOM     10  NE2 HIS A   1     -16.488  -0.965 -12.516  1.00 20.05           N
ATOM    927  N   ASN A  60     -15.178  -7.233 -42.292  1.00 18.01           N
ATOM    928  CA  ASN A  60     -15.029  -6.740 -40.929  1.00 16.90           C
ATOM    929  C   ASN A  60     -16.363  -6.747 -40.221  1.00 18.98           C
ATOM    930  O   ASN A  60     -17.406  -6.445 -40.810  1.00 23.95           O
ATOM    931  CB  ASN A  60     -14.456  -5.320 -40.895  1.00 20.16           C
ATOM    932  CG  ASN A  60     -12.957  -5.276 -41.123  1.00 29.51           C
ATOM    933  OD1 ASN A  60     -12.423  -5.994 -41.962  1.00 25.82           O
ATOM    934  ND2 ASN A  60     -12.263  -4.422 -40.366  1.00 20.64           N
TER
HETATM 7062  C1  NAG C   1     -10.860  -4.178 -40.436  1.00 27.06           C
HETATM 7063  C2  NAG C   1     -10.681  -2.675 -40.638  1.00 25.95           C
HETATM 7064  C3  NAG C   1      -9.195  -2.328 -40.662  1.00 32.04           C
HETATM 7065  C4  NAG C   1      -8.473  -2.921 -39.456  1.00 24.81           C
HETATM 7066  C5  NAG C   1      -8.814  -4.402 -39.257  1.00 26.90           C
HETATM 7067  C6  NAG C   1      -8.349  -4.931 -37.923  1.00 24.05           C
HETATM 7068  C7  NAG C   1     -12.458  -1.556 -41.962  1.00 32.11           C
HETATM 7069  C8  NAG C   1     -12.975  -1.282 -43.346  1.00 38.01           C
HETATM 7070  N2  NAG C   1     -11.331  -2.280 -41.886  1.00 29.56           N
HETATM 7071  O3  NAG C   1      -9.051  -0.912 -40.647  1.00 37.51           O
HETATM 7072  O4  NAG C   1      -7.066  -2.854 -39.649  1.00 32.97           O
HETATM 7073  O5  NAG C   1     -10.230  -4.577 -39.274  1.00 25.62           O
HETATM 7074  O6  NAG C   1      -8.560  -6.336 -37.806  1.00 34.61           O
HETATM 7075  O7  NAG C   1     -13.016  -1.117 -40.977  1.00 39.12           O
HETATM 7090  C1  NAG C   2      -6.436  -1.792 -39.032  1.00 39.57           C
HETATM 7091  C2  NAG C   2      -4.930  -2.007 -39.190  1.00 39.58           C
HETATM 7092  C3  NAG C   2      -4.176  -1.020 -38.317  1.00 46.39           C
HETATM 7093  C4  NAG C   2      -4.564   0.390 -38.715  1.00 39.21           C
HETATM 7094  C5  NAG C   2      -6.077   0.572 -38.599  1.00 43.91           C
HETATM 7095  C6  NAG C   2      -6.529   1.927 -39.101  1.00 54.83           C
HETATM 7096  C7  NAG C   2      -4.374  -4.309 -39.775  1.00 47.57           C
HETATM 7097  C8  NAG C   2      -3.968  -5.662 -39.261  1.00 35.21           C
HETATM 7098  N2  NAG C   2      -4.542  -3.365 -38.856  1.00 33.73           N
HETATM 7099  O3  NAG C   2      -2.774  -1.242 -38.441  1.00 49.87           O
HETATM 7100  O4  NAG C   2      -3.912   1.331 -37.868  1.00 72.64           O
HETATM 7101  O5  NAG C   2      -6.772  -0.429 -39.367  1.00 39.07           O
HETATM 7102  O6  NAG C   2      -7.521   1.844 -40.116  1.00 67.90           O
HETATM 7103  O7  NAG C   2      -4.563  -4.085 -40.967  1.00 41.41           O
HETATM 7117  C1  MAN C   3      -2.909   2.196 -38.339  1.00 72.39           C
HETATM 7118  C2  MAN C   3      -3.477   3.339 -37.453  1.00 78.40           C
HETATM 7119  C3  MAN C   3      -2.863   4.712 -37.881  1.00 85.45           C
HETATM 7120  C4  MAN C   3      -1.409   4.553 -38.486  1.00 84.13           C
HETATM 7121  C5  MAN C   3      -1.350   3.521 -39.675  1.00 83.85           C
HETATM 7122  C6  MAN C   3      -0.026   2.694 -39.722  1.00 66.63           C
HETATM 7127  O2  MAN C   3      -3.134   3.083 -36.062  1.00 53.62           O
HETATM 7128  O3  MAN C   3      -2.865   5.706 -36.798  1.00 62.87           O
HETATM 7129  O4  MAN C   3      -0.935   5.818 -38.957  1.00 80.75           O
HETATM 7130  O5  MAN C   3      -2.482   2.608 -39.621  1.00 89.03           O
HETATM 7131  O6  MAN C   3      -0.307   1.350 -40.109  1.00 68.09           O
HETATM 7217 CU    CU A 301     -13.829   0.757 -15.122  1.00 19.10          CU
END
"""

# YI1: an O-NH3+ (aminooxy) with a Zn 2.05 A from NAV (geostd coordinates)
pdb_yi1_zn = """\
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
HETATM    1  FAG YI1 A   1      14.072  20.211  18.737  1.00 20.00           F
HETATM    2  CAE YI1 A   1      15.125  20.820  19.390  1.00 20.00           C
HETATM    3  FAF YI1 A   1      15.502  21.850  18.546  1.00 20.00           F
HETATM    4  PAB YI1 A   1      14.542  21.511  21.048  1.00 20.00           P
HETATM    5  OAC YI1 A   1      13.360  22.400  20.669  1.00 20.00           O
HETATM    6  OAA YI1 A   1      14.179  20.249  21.821  1.00 20.00           O
HETATM    7  OAD YI1 A   1      15.769  22.260  21.557  1.00 20.00           O
HETATM    8  CAH YI1 A   1      16.261  19.851  19.482  1.00 20.00           C
HETATM    9  CAJ YI1 A   1      16.052  18.478  19.424  1.00 20.00           C
HETATM   10  CAK YI1 A   1      17.128  17.607  19.537  1.00 20.00           C
HETATM   11  CAL YI1 A   1      18.412  18.093  19.717  1.00 20.00           C
HETATM   12  CAM YI1 A   1      18.638  19.469  19.794  1.00 20.00           C
HETATM   13  CAI YI1 A   1      17.551  20.333  19.680  1.00 20.00           C
HETATM   14  CAN YI1 A   1      20.000  20.000  20.000  1.00 20.00           C
HETATM   15  CAQ YI1 A   1      20.920  19.323  20.803  1.00 20.00           C
HETATM   16  CAR YI1 A   1      22.195  19.825  21.001  1.00 20.00           C
HETATM   17  CAS YI1 A   1      22.586  21.021  20.407  1.00 20.00           C
HETATM   18  CAP YI1 A   1      21.676  21.700  19.604  1.00 20.00           C
HETATM   19  CAO YI1 A   1      20.402  21.197  19.403  1.00 20.00           C
HETATM   20  CAT YI1 A   1      23.969  21.547  20.602  1.00 20.00           C
HETATM   21  OAU YI1 A   1      24.759  20.996  19.527  1.00 20.00           O
HETATM   22  NAV YI1 A   1      26.069  21.421  19.619  1.00 20.00           N
HETATM   23 ZN    ZN A   2      27.085  20.582  18.049  1.00 20.00          ZN
END
"""

# the Lys-Cd fixture turned so that the Cd lies on the 2-fold of P 1 2 1
# (x = 0, z = 0) with NZ-Cd perpendicular to the axis
pdb_lys_cd_2fold = """\
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1 2 1
ATOM      1  N   LYS A 359       8.812  19.431   0.327  1.00  8.43           N
ATOM      2  CA  LYS A 359       7.392  19.240   0.029  1.00  7.92           C
ATOM      3  C   LYS A 359       7.051  17.779  -0.346  1.00  7.33           C
ATOM      4  O   LYS A 359       6.310  17.521  -1.331  1.00  7.65           O
ATOM      5  CB  LYS A 359       6.567  19.740   1.222  1.00  8.93           C
ATOM      6  CG  LYS A 359       5.074  19.480   0.992  1.00 11.26           C
ATOM      7  CD  LYS A 359       4.149  20.438   1.657  1.00 16.15           C
ATOM      8  CE  LYS A 359       2.693  20.000   1.437  1.00 18.69           C
ATOM      9  NZ  LYS A 359       2.214  20.000   0.000  1.00 18.51           N
TER
HETATM   10 CD    CD A 396       0.000  20.000   0.000  1.00 10.07          CD
END
"""

CONVERTED = 'Coordinating NH3+ placed as NH2'
UNRESOLVED = 'Coordinating NH3+ left as NH3+ (unresolved)'
SUPPLIED = 'Supplied NH3+ at metal sites'
SUPPLIED_NH2 = 'Supplied NH2 at metal sites'
CHAIN_START = 'Metal-bound chain starts'
GLY288 = ('N GLY A 288: metal-bound chain start, placed as NH2 '
  '(n_terminal_charge=%s would place 1 H, wrong for a terminus and for a '
  'chain break)')
REPORTED = 'N-H next to a metal, not changed'
ON = 'metal_bound_amines_as_nh2=True'

def run(prefix, pdb_str, args=(), result=False):
  '''reduce2 on pdb_str. :return: (output hierarchy, description, log
  [, program result]).'''
  fn, out, txt = prefix + '_in.pdb', prefix + '_out.pdb', prefix + '.txt'
  with open(fn, 'w') as f:
    f.write(pdb_str)
  for n in (out, txt):
    if os.path.exists(n): os.remove(n)
  log = StringIO()
  r = run_program(program_class=reduce2.Program, logger=log, args=[fn,
    'output.filename=' + out, 'output.description_file_name=' + txt]
    + list(args))
  with open(txt) as f:
    description = f.read()
  h = iotbx.pdb.input(out).construct_hierarchy()
  if result: return h, description, log.getvalue(), r
  return h, description, log.getvalue()

def atom(hierarchy, resseq, name, altloc=''):
  found = [a for a in hierarchy.atoms()
           if a.parent().parent().resseq_as_int() == resseq
           and a.name.strip() == name and a.parent().altloc == altloc]
  assert len(found) == 1, (resseq, name, altloc)
  return found[0]

def h_on(hierarchy, resseq, name):
  '''(name, altloc) of the H/D within 1.2 A of atom name in residue resseq.'''
  parent = atom(hierarchy, resseq, name)
  return sorted((a.name.strip(), a.parent().altloc)
    for a in hierarchy.atoms()
    if a.element.strip() in ('H', 'D') and a.distance(parent) < 1.2)

def h_xyz(hierarchy):
  return dict((a.id_str(), a.xyz) for a in hierarchy.atoms()
    if a.element.strip() in ('H', 'D'))

def angle(a, b, c):
  return (matrix.col(a.xyz) - matrix.col(b.xyz)).angle(
    matrix.col(c.xyz) - matrix.col(b.xyz), deg=True)

def section(text, title):
  '''The lines of a report section whose title starts with title.'''
  lines = text.splitlines()
  starts = [i for i, l in enumerate(lines) if l.startswith(title)]
  assert len(starts) == 1, (title, text)
  result = []
  for l in lines[starts[0] + 1:]:
    if not l.startswith('  ') or l.strip().startswith('Time to'): break
    result.append(l.strip())
  return result

def check_nh2(h, resseq, n_name, c_name, metal_resseq, metal_name, names):
  '''Two H, tetrahedral with the metal as the third substituent.'''
  assert [n for n, alt in h_on(h, resseq, n_name)] == names, \
    h_on(h, resseq, n_name)
  n, c = atom(h, resseq, n_name), atom(h, resseq, c_name)
  m = atom(h, metal_resseq, metal_name)
  h1, h2 = [atom(h, resseq, x) for x in names]
  assert 104 < angle(h1, n, h2) < 115, angle(h1, n, h2)
  for x in (h1, h2):
    assert 0.85 < x.distance(n) < 1.1, x.distance(n)
    assert 104 < angle(c, n, x) < 115, angle(c, n, x)
    assert 100 < angle(m, n, x) < 115, angle(m, n, x)
    assert x.distance(m) > 2.3, x.distance(m)

def test_lys_cd():
  h, description, log = run('tst_reduce2_metal_amine_lys', pdb_lys_cd, [ON])
  check_nh2(h, 359, 'NZ', 'CE', 396, 'CD', ['HZ2', 'HZ3'])
  for text in (description, log):
    assert section(text, CONVERTED) == ['NZ  LYS A 359: HZ1 not placed  '
      '(NZ  LYS A 359 - CD    CD A 396 2.21 A (geometric))'], text
    assert REPORTED not in text

def test_option_off():
  h, description, log = run('tst_reduce2_metal_amine_off', pdb_lys_cd)
  assert [n for n, alt in h_on(h, 359, 'NZ')] == ['HZ1', 'HZ2', 'HZ3']
  assert CONVERTED not in description and CONVERTED not in log
  assert section(description, REPORTED)[-1].startswith(
    'NZ  LYS A 359 (HZ1, HZ2, HZ3)'), description

def test_terminus_cu():
  # resseq 1: placement builds NH3 (n_terminal_charge=residue_one)
  h, description, log = run('tst_reduce2_metal_amine_his', pdb_his_cu, [ON])
  check_nh2(h, 1, 'N', 'CA', 301, 'CU', ['H2', 'H3'])
  assert section(description, CONVERTED) == ['N   HIS A   1: H1 not placed  '
    '(N   HIS A   1 - CU    CU A 301 2.13 A (geometric))'], description

def test_neutral_nh2_kept():
  h, description, log = run('tst_reduce2_metal_amine_edn', pdb_edn_cu, [ON])
  assert [n for n, alt in h_on(h, 1, 'NAA')] == ['HAA1', 'HAA2']
  assert [n for n, alt in h_on(h, 1, 'NAD')] == ['HAD1', 'HAD2']
  assert CONVERTED not in description

def test_chain_start_settings():
  # Gly 288 starts the chain in the model: NH2 whatever n_terminal_charge
  # says (no sequence records are read)
  for setting in ('residue_one', 'no_charge', 'all_chain_starts'):
    h, description, log = run('tst_reduce2_metal_amine_gly_' + setting,
      pdb_1xmk, [ON, 'n_terminal_charge=' + setting])
    assert h_on(h, 288, 'N') == [('H2', ''), ('H3', '')], (setting,
      h_on(h, 288, 'N'))
    check_nh2(h, 288, 'N', 'CA', 398, 'NI', ['H2', 'H3'])
    expected = (GLY288 % setting if setting != 'all_chain_starts' else
      'N GLY A 288: metal-bound chain start, placed as NH2')
    for text in (description, log):
      assert section(text, CHAIN_START) == [expected], (setting, text)

def test_chain_start_preceding_residue():
  # a residue before Gly 288 with its C bonded to N: no chain start, an amide
  # on the Ni (step B: no H). Passes without C1c.
  pdb = pdb_1xmk.replace('ATOM      1  N   GLY A 288', """\
ATOM     41  N   GLY A 287      11.272  25.005  34.136  1.00 20.00           N
ATOM     42  CA  GLY A 287      11.618  24.171  32.989  1.00 20.00           C
ATOM     43  C   GLY A 287      12.237  24.201  31.601  1.00 20.00           C
ATOM     44  O   GLY A 287      12.648  25.262  31.134  1.00 20.00           O
ATOM      1  N   GLY A 288""")
  h, description, log = run('tst_reduce2_metal_amine_gly_prev', pdb, [ON])
  assert h_on(h, 288, 'N') == [], h_on(h, 288, 'N')
  assert 'metal-bound chain start' not in description
  assert 'H   GLY A 288  (N   GLY A 288 - NI' in '\n'.join(
    section(description, 'N-H not placed')), description

def test_chain_start_refused():
  # two Ni positions (A, B): no unique frame, NH3 kept with the reason
  pdb = pdb_1xmk.replace(
    'HETATM   21 NI    NI A 398      13.849  21.950  30.990  1.00',
    'HETATM   21 NI  A NI A 398      13.849  21.950  30.990  0.50').replace(
    'END\n', 'HETATM   22 NI  B NI A 398      13.849  22.450  30.990  0.50 20.00'
    '          NI\nEND\n')
  h, description, log = run('tst_reduce2_metal_amine_gly_alt', pdb, [ON])
  assert [x for x, alt in h_on(h, 288, 'N')] == ['H1', 'H2', 'H3'], \
    h_on(h, 288, 'N')
  lines = section(description, CHAIN_START)
  assert len(lines) == 1 and lines[0].startswith('N GLY A 288: metal-bound '
    'chain start, NH3 kept: no unique metal position') and lines[0].endswith(
    '(n_terminal_charge=residue_one would place 1 H)'), lines

def test_chain_start_without_metal():
  # no metal: as n_terminal_charge places it (passes without C1c)
  pdb = '\n'.join(l for l in pdb_1xmk.splitlines() if ' NI ' not in l) + '\n'
  for setting in ('residue_one', 'no_charge'):
    h, description, log = run('tst_reduce2_metal_amine_gly_nometal', pdb,
      [ON, 'n_terminal_charge=' + setting])
    assert h_on(h, 288, 'N') == [('H', '')], (setting, h_on(h, 288, 'N'))
    assert CHAIN_START not in description

def test_chain_start_rerun():
  # add -> optimize -> add keep_existing_H=True, and plain add on the output:
  # NH2 every time, at the same positions
  h0, d, l = run('tst_reduce2_metal_amine_gr0', pdb_1xmk, [ON])
  previous = 'tst_reduce2_metal_amine_gr0_out.pdb'
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
  h0, d, l = run('tst_reduce2_metal_amine_gacs0', pdb_1xmk, [ON])
  with open('tst_reduce2_metal_amine_gacs0_out.pdb') as f:
    text = f.read()
  h, d, l = run('tst_reduce2_metal_amine_gacs1', text,
    [ON, 'keep_existing_H=True', 'n_terminal_charge=all_chain_starts'])
  assert h_on(h, 288, 'N') == [('H2', ''), ('H3', '')], h_on(h, 288, 'N')
  for name in ('H2', 'H3'):
    assert atom(h, 288, name).distance(atom(h0, 288, name)) < 1.e-3, name
  assert section(d, CHAIN_START) == ['N GLY A 288: metal-bound chain start, '
    'unchanged: supplied NH2 kept (H2, H3)'], d
  assert 'not judged' not in d

def test_chain_start_option_off():
  # without the option: one H on Gly 288 (residue_one), as on master
  h, description, log = run('tst_reduce2_metal_amine_gly_off', pdb_1xmk)
  assert h_on(h, 288, 'N') == [('H', '')]
  assert CHAIN_START not in description and CONVERTED not in description

def gly288_conformers(metal_alt, split_n=True):
  '''pdb_1xmk with Gly 288 in A and B (B 0.1 A away; with split_n False the
  N stays blank) and the Ni only in altloc metal_alt, all at 0.5 (the Ni
  blank and full if metal_alt is None).'''
  lines = []
  for l in pdb_1xmk.splitlines():
    if l.startswith('ATOM') and l[17:26] == 'GLY A 288':
      if l[12:16] == ' N  ' and not split_n:
        lines.append(l)
        continue
      lines.append(l[:16] + 'A' + l[17:54] + '  0.50' + l[60:])
      x = float(l[30:38]) + 0.1
      lines.append(l[:16] + 'B' + l[17:30] + '%8.3f' % x + l[38:54] + '  0.50'
        + l[60:])
    elif ' NI ' in l and metal_alt is not None:
      lines.append(l[:16] + metal_alt + l[17:54] + '  0.50' + l[60:])
    else:
      lines.append(l)
  return '\n'.join(lines) + '\n'

def test_chain_start_supplied_h():
  # supplied H at a metal-bound chain start: the option builds no NH3 and
  # deletes nothing; every supplied atom keeps its name and coordinates
  run('tst_reduce2_metal_amine_gs0', pdb_1xmk)
  with open('tst_reduce2_metal_amine_gs0_out.pdb') as f:
    one_h = f.read()
  run('tst_reduce2_metal_amine_gs1', pdb_1xmk, [ON])
  with open('tst_reduce2_metal_amine_gs1_out.pdb') as f:
    h2 = [l for l in f.read().splitlines() if l[12:16] == ' H2 '][0]
  h1 = [l for l in one_h.splitlines() if l[12:26] == ' H   GLY A 288'][0]
  two_h = one_h.replace(h1, h1 + '\n' + h2)
  for tag, text, names in (('a', one_h, ['H']), ('b', two_h, ['H', 'H2'])):
    h_in = iotbx.pdb.input(source_info=None, lines=text).construct_hierarchy()
    h, d, l = run('tst_reduce2_metal_amine_gs_' + tag, text,
      [ON, 'keep_existing_H=True'])
    off, d_off, l_off = run('tst_reduce2_metal_amine_gs_off', text,
      ['keep_existing_H=True'])
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

def test_chain_start_conformers():
  # Gly 288 N in A and B, the Ni in one of them: that conformer NH2, the
  # other one H (residue_one), the line names the conformer
  for metal_alt, other in (('B', 'A'), ('A', 'B')):
    h, d, l = run('tst_reduce2_metal_amine_gc' + metal_alt,
      gly288_conformers(metal_alt), [ON])
    n = dict((alt, atom(h, 288, 'N', alt)) for alt in ('A', 'B'))
    on = lambda alt: sorted(x.name.strip() for x in h.atoms()
      if x.element.strip() == 'H' and x.parent().altloc == alt
      and x.distance(n[alt]) < 1.2)
    assert on(metal_alt) == ['H2', 'H3'], (metal_alt, on(metal_alt))
    assert on(other) == ['H'], (metal_alt, on(other))
    assert section(d, CHAIN_START) == [GLY288.replace('N GLY A 288:',
      'N GLY A 288 conformer %s:' % metal_alt) % 'residue_one'], d
  # a blank N with the Ni only in B: no common frame, unchanged
  h, d, l = run('tst_reduce2_metal_amine_gcn', gly288_conformers('B', False),
    [ON])
  assert [x for x, alt in h_on(h, 288, 'N')] == ['H', 'H'], h_on(h, 288, 'N')
  assert section(d, CHAIN_START) == ['N GLY A 288: metal-bound chain start, '
    'unchanged: metal only in conformer B (no common frame)'], d
  # a blank N, CA/C/O in A and B, a blank Ni, a supplied H in A only: A keeps
  # it, B gets what the option-off run places
  run('tst_reduce2_metal_amine_gcs0', gly288_conformers(None, False))
  with open('tst_reduce2_metal_amine_gcs0_out.pdb') as f:
    text = '\n'.join(l for l in f.read().splitlines()
      if l[12:26] != ' H  BGLY A 288') + '\n'
  h_in = iotbx.pdb.input(source_info=None, lines=text).construct_hierarchy()
  off, d_off, l_off = run('tst_reduce2_metal_amine_gcs_off', text,
    ['keep_existing_H=True'])
  h, d, l = run('tst_reduce2_metal_amine_gcs', text,
    [ON, 'keep_existing_H=True'])
  assert atom(h, 288, 'H', 'A').distance(atom(h_in, 288, 'H', 'A')) < 1.e-3
  assert atom(h, 288, 'H', 'B').distance(atom(off, 288, 'H', 'B')) < 1.e-3
  assert section(d, CHAIN_START) == ['N GLY A 288: metal-bound chain start, '
    'unchanged: supplied H in conformer A kept (H); conformer B as '
    'n_terminal_charge places it'], d

def test_planar_nh2_guard():
  # Asn ND2 next to a Zn: an amide N, not an amine (passes without the option)
  h, description, log = run('tst_reduce2_metal_amine_asn', pdb_asn_zn, [ON])
  assert [n for n, alt in h_on(h, 1, 'ND2')] == ['HD21', 'HD22']
  assert CONVERTED not in description
  assert section(description, REPORTED)[0].startswith(
    'ND2 ASN A   1 (HD21, HD22)'), description

def test_unresolved():
  # two alternate Cd positions: full coverage, but no single frame
  pdb = pdb_lys_cd.replace(
    'HETATM   10 CD    CD A 396      10.293  10.651  44.885  1.00',
    'HETATM   10 CD  A CD A 396      10.293  10.651  44.885  0.50').replace(
    'END\n', 'HETATM   11 CD  B CD A 396      11.011  10.159  44.542  0.50 10.07'
    '          CD\nEND\n')
  h, description, log = run('tst_reduce2_metal_amine_alt', pdb, [ON])
  assert [n for n, alt in h_on(h, 359, 'NZ')] == ['HZ1', 'HZ2', 'HZ3']
  lines = section(description, UNRESOLVED)
  assert len(lines) == 1 and lines[0].startswith(
    'NZ  LYS A 359: no unique metal position'), lines

def test_supplied():
  # supplied NH3 (reduce2's own output without the option): kept, reported
  run('tst_reduce2_metal_amine_nh3', pdb_lys_cd)
  with open('tst_reduce2_metal_amine_nh3_out.pdb') as f:
    nh3 = f.read()
  h, description, log = run('tst_reduce2_metal_amine_nh3_keep', nh3,
    [ON, 'keep_existing_H=True'])
  assert [n for n, alt in h_on(h, 359, 'NZ')] == ['HZ1', 'HZ2', 'HZ3']
  assert section(description, SUPPLIED)[0].startswith(
    'NZ  LYS A 359 (HZ1, HZ2, HZ3)  (NZ  LYS A 359 - CD'), description
  assert CONVERTED not in description

def test_rerun_stable():
  # add -> optimize -> add keep_existing_H=True: same H, same positions; a
  # supplied NH2 is not completed (placement's new HZ1 is dropped)
  h0, description, log = run('tst_reduce2_metal_amine_rr0', pdb_lys_cd, [ON])
  for i, args in enumerate((['approach=optimize'], ['keep_existing_H=True'])):
    prefix = 'tst_reduce2_metal_amine_rr%d' % (i + 1)
    with open('tst_reduce2_metal_amine_rr%d_out.pdb' % i) as f:
      text = f.read()
    h, description, log = run(prefix, text, [ON] + args)
    off, d_off, l_off = run(prefix + '_off', text, args)
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

def test_qmready_charge():
  # one converted site: total charge one lower
  pdb = pdb_lys_cd.replace('          CD\n', '          CD2+\n')
  charges = []
  for args in ([], [ON]):
    h, description, log, r = run('tst_reduce2_metal_amine_qm', pdb,
      ['qmready=True'] + args, result=True)
    charges.append(r.total_charge)
  assert charges[1] == charges[0] - 1, charges

def test_hd_and_h():
  h, description, log = run('tst_reduce2_metal_amine_hd', pdb_lys_cd,
    [ON, 'neutron_option=hd_and_h'])
  assert h_on(h, 359, 'NZ') == [('DZ2', 'B'), ('DZ3', 'B'), ('HZ2', 'A'),
    ('HZ3', 'A')], h_on(h, 359, 'NZ')

def lys_cd_at(text, angle):
  '''text with Cd 396 moved to C-N-Cd angle (deg) at the same N-Cd distance,
  in the CE-NZ-Cd plane.'''
  h = iotbx.pdb.input(source_info=None, lines=text).construct_hierarchy()
  n, ce, cd = [matrix.col(atom(h, r, x).xyz)
    for r, x in ((359, 'NZ'), (359, 'CE'), (396, 'CD'))]
  u = (ce - n).normalize()
  w = (cd - n) - u * (cd - n).dot(u)
  w = w.normalize()
  a = math.radians(angle)
  new = n + (u * math.cos(a) + w * math.sin(a)) * (cd - n).length()
  lines = []
  for l in text.splitlines():
    if l.startswith('HETATM') and l[17:20] == ' CD':
      l = l[:30] + '%8.3f%8.3f%8.3f' % new.elems + l[54:]
    lines.append(l)
  return '\n'.join(lines) + '\n'

def test_geometry_not_placeable():
  # C-N-Cd 140 deg: riding would put both H on one site; 145: none. No H is
  # removed: new H as without the option, supplied H as they are
  run('tst_reduce2_metal_amine_g0', pdb_lys_cd, [ON])
  with open('tst_reduce2_metal_amine_g0_out.pdb') as f:
    nh2 = f.read()
  for angle in (140, 145):
    pdb = lys_cd_at(pdb_lys_cd, angle)
    off, d_off, l_off = run('tst_reduce2_metal_amine_g%d_off' % angle, pdb)
    h, description, log = run('tst_reduce2_metal_amine_g%d' % angle, pdb, [ON])
    assert CONVERTED not in description
    lines = section(description, UNRESOLVED)
    assert len(lines) == 1 and lines[0].startswith('NZ  LYS A 359: geometry '
      'not placeable, C-N-M %d deg' % angle), lines
    assert h_on(h, 359, 'NZ') == h_on(off, 359, 'NZ') == [('HZ1', ''),
      ('HZ2', ''), ('HZ3', '')]
    for name in ('HZ1', 'HZ2', 'HZ3'):
      assert atom(h, 359, name).distance(atom(off, 359, name)) < 1.e-3, name
    # supplied NH2
    supplied = lys_cd_at(nh2, angle)
    h0 = iotbx.pdb.input(source_info=None, lines=supplied).construct_hierarchy()
    h, description, log = run('tst_reduce2_metal_amine_g%d_keep' % angle,
      supplied, [ON, 'keep_existing_H=True'])
    assert h_on(h, 359, 'NZ') == [('HZ2', ''), ('HZ3', '')], h_on(h, 359, 'NZ')
    for name in ('HZ2', 'HZ3'):
      assert atom(h, 359, name).distance(atom(h0, 359, name)) < 1.e-3, name
    # kept as supplied: not placed, not checked
    assert UNRESOLVED not in description and CONVERTED not in description
    assert section(description, SUPPLIED_NH2)[0].startswith(
      'NZ  LYS A 359 (HZ2, HZ3)'), description

def test_supplied_nh2_any_metal():
  # a supplied NH2 is not completed to NH3 whatever the metal site: Cd one
  # cell away (contact through symmetry), Cd at occupancy 0.5
  run('tst_reduce2_metal_amine_s0', pdb_lys_cd, [ON])
  with open('tst_reduce2_metal_amine_s0_out.pdb') as f:
    nh2 = f.read()
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

def test_only_reference_angles_removed():
  # an angle already on C-N-M (here 107 deg, weight 123, origin edits)
  # survives placement unchanged
  import mmtbx.model
  from libtbx.utils import null_out
  from cctbx import geometry_restraints
  from mmtbx.hydrogens import reduce_hydrogen
  class placement(reduce_hydrogen.place_hydrogens):
    def add_link_h_restraints(self):
      reduce_hydrogen.place_hydrogens.add_link_h_restraints(self)
      atoms = self.model.get_hierarchy().atoms()
      i = dict((a.name.strip(), a.i_seq) for a in atoms)
      self.model.get_restraints_manager().geometry.add_angles_in_place([
        geometry_restraints.angle_proxy(i_seqs = [i['CE'], i['NZ'], i['CD']],
          angle_ideal = 107, weight = 123, origin_id = 4)])
  m = mmtbx.model.manager(model_input = iotbx.pdb.input(source_info=None,
    lines=pdb_lys_cd), log = null_out())
  p = placement(m, metal_bound_amines_as_nh2 = True)
  p.run()
  assert p.metal_amine_converted_text, p.metal_amine_unresolved
  m = p.get_model()
  i = dict((a.name.strip(), a.i_seq) for a in m.get_hierarchy().atoms())
  found = [(round(x.angle_ideal, 3), x.weight, x.origin_id)
    for x in m.get_restraints_manager().geometry.angle_proxies
    if set(x.i_seqs) == set([i['CE'], i['NZ'], i['CD']])]
  assert found == [(107, 123, 4)], found
  # and the injected H-N-M angles are gone
  assert len([x for x in m.get_restraints_manager().geometry.angle_proxies
    if i['CD'] in x.i_seqs]) == 1

def test_n_o_not_converted():
  # an O-NH3+ (YI1 NAV) next to a Zn: not a primary amine, unchanged
  h, description, log = run('tst_reduce2_metal_amine_yi1', pdb_yi1_zn, [ON])
  assert [n for n, alt in h_on(h, 1, 'NAV')] == ['HAV1', 'HAV2', 'HAV3']
  assert CONVERTED not in description and UNRESOLVED not in description

def test_lists_after_gc():
  # the lists are built at the end of placement from what step B saw; the
  # atoms of that earlier hierarchy may be collected by then
  import gc
  from mmtbx.hydrogens import reduce_hydrogen
  finish = reduce_hydrogen.place_hydrogens._finish_metal_lists
  def collected_first(self):
    gc.collect()
    return finish(self)
  reduce_hydrogen.place_hydrogens._finish_metal_lists = collected_first
  try:
    h, description, log = run('tst_reduce2_metal_amine_gc', pdb_1xmk)
  finally:
    reduce_hydrogen.place_hydrogens._finish_metal_lists = finish
  assert section(description, REPORTED) == [
    'N   GLY A 288 (H)  (N   GLY A 288 - NI    NI A 398 1.89 A (geometric))'
    ], description

def test_other_h_unchanged():
  # with a converted amine, every H away from it as without the option
  # (placement's riding H manager is kept; a new one moved the glycan H)
  off, d, l = run('tst_reduce2_metal_amine_gly_off', pdb_his_cu_glycan)
  on, d, l = run('tst_reduce2_metal_amine_gly_on', pdb_his_cu_glycan, [ON])
  assert section(d, CONVERTED)[0].startswith('N   HIS A   1: H1 not placed')
  n = atom(on, 1, 'N')
  xyz = lambda h: dict((a.id_str(), a.xyz) for a in h.atoms()
    if a.element.strip() == 'H' and a.distance(n) > 1.2)
  a, b = xyz(off), xyz(on)
  assert sorted(a) == sorted(b)
  moved = [k for k in a if (matrix.col(a[k]) - matrix.col(b[k])).length()
           > 1.e-3]
  assert not moved, moved

def lys_cif(fn, changes):
  '''geostd LYS with the rows changed (old, new), written to fn.'''
  import libtbx.load_env
  geostd = libtbx.env.find_in_repositories(
    relative_path='chem_data/geostd/l/data_LYS.cif', test=os.path.isfile)
  assert geostd is not None
  with open(geostd) as f:
    cif = f.read()
  for old, new in changes:
    assert cif.count(old) == 1, old
    cif = cif.replace(old, new)
  with open(fn, 'w') as f:
    f.write(cif)

def test_unequal_targets():
  # a LYS dictionary with C-N-HZ2 105, C-N-HZ3 115 deg and unequal N-H: the
  # coordinates place_metal_nh2() accepts are the final ones
  lys_cif('tst_reduce2_metal_amine_lys.cif', (
      (' LYS  NZ  HZ3  coval  0.890  0.020  1.040',
       ' LYS  NZ  HZ3  coval  0.950  0.020  1.100'),
      (' LYS  CE   NZ  HZ2  110.000  3.000', ' LYS  CE   NZ  HZ2  105.000  3.000'),
      (' LYS  CE   NZ  HZ3  110.000  3.000', ' LYS  CE   NZ  HZ3  115.000  3.000')))
  h, description, log = run('tst_reduce2_metal_amine_unequal', pdb_lys_cd,
    [ON, 'tst_reduce2_metal_amine_lys.cif'])
  assert section(description, CONVERTED) == ['NZ  LYS A 359: HZ1 not placed  '
    '(NZ  LYS A 359 - CD    CD A 396 2.21 A (geometric))'], description
  n, ce, cd = [matrix.col(atom(h, r, x).xyz)
    for r, x in ((359, 'NZ'), (359, 'CE'), (396, 'CD'))]
  from libtbx import group_args
  from mmtbx.hydrogens.reduce_hydrogen import place_metal_nh2
  placed, reason = place_metal_nh2(n, ce, cd, group_args(d = [0.89, 0.95],
    cnh = [105., 115.], hnh = 109.))
  assert reason is None, reason
  for name, xyz in zip(('HZ2', 'HZ3'), placed):
    d = (matrix.col(atom(h, 359, name).xyz) - xyz).length()
    assert d < 2.e-3, (name, d)
  assert abs(atom(h, 359, 'HZ3').distance(atom(h, 359, 'NZ')) - 0.95) < 2.e-3
  assert angle(atom(h, 359, 'CE'), atom(h, 359, 'NZ'), atom(h, 359, 'HZ2')) < \
    angle(atom(h, 359, 'CE'), atom(h, 359, 'NZ'), atom(h, 359, 'HZ3'))

def test_symmetry_contact():
  # the same Cd through a translation and through a rotation
  # (-x+1/2,-y,z+1/2): the same H as with the Cd itself
  from cctbx import crystal, sgtbx
  h0, d0, l0 = run('tst_reduce2_metal_amine_sym0', pdb_lys_cd, [ON])
  cs = crystal.symmetry(unit_cell=(35.557, 43.526, 45.471, 90, 90, 90),
    space_group_symbol='P 21 21 21')
  uc = cs.unit_cell()
  cd = matrix.col(atom(h0, 396, 'CD').xyz)
  for tag, op in (('t', 'x+1,y,z'), ('r', '-x+1/2,-y,z+1/2')):
    rt = sgtbx.rt_mx(op)
    xyz = uc.orthogonalize(rt.inverse() * uc.fractionalize(cd.elems))
    pdb = pdb_lys_cd.replace('CD A 396      10.293  10.651  44.885',
      'CD A 396    %8.3f%8.3f%8.3f' % tuple(xyz))
    h, description, log = run('tst_reduce2_metal_amine_sym' + tag, pdb, [ON])
    lines = section(description, CONVERTED)
    assert len(lines) == 1 and ' sym ' in lines[0], (tag, lines)
    for name in ('HZ2', 'HZ3'):
      d = atom(h, 359, name).distance(atom(h0, 359, name))
      assert d < 2.e-3, (tag, name, d)

def test_shared_frame_restraints():
  # blank N, C and Cd, H groups A and B (O split): both converted, and the
  # angle restraints of the rest, with one already on C-N-M, unchanged in
  # count and metadata
  import mmtbx.model
  from libtbx.utils import null_out
  from cctbx import geometry_restraints
  from mmtbx.hydrogens import reduce_hydrogen
  class placement(reduce_hydrogen.place_hydrogens):
    def add_link_h_restraints(self):
      reduce_hydrogen.place_hydrogens.add_link_h_restraints(self)
      i = dict((a.name.strip(), a.i_seq)
        for a in self.model.get_hierarchy().atoms())
      self.model.get_restraints_manager().geometry.add_angles_in_place([
        geometry_restraints.angle_proxy(i_seqs = [i['CE'], i['NZ'], i['CD']],
          angle_ideal = 107, weight = 123, origin_id = 4)])
  pdb = pdb_lys_cd.replace(
    'ATOM      4  O   LYS A 359       7.955   4.152  45.067  1.00',
    'ATOM      4  O  ALYS A 359       7.955   4.152  45.067  0.50').replace(
    'ATOM      5  CB', 'ATOM      4  O  BLYS A 359       7.855   4.252  45.067'
    '  0.50  7.65           O\nATOM      5  CB')
  result = []
  for on in (False, True):
    m = mmtbx.model.manager(model_input = iotbx.pdb.input(source_info=None,
      lines=pdb), log = null_out())
    p = placement(m, metal_bound_amines_as_nh2 = on)
    p.run()
    m = p.get_model()
    atoms = m.get_hierarchy().atoms()
    nz_h = set(a.i_seq for a in atoms if a.name.strip().startswith('HZ'))
    label = lambda i: atoms[i].id_str()
    result.append((sorted((tuple(label(i) for i in x.i_seqs), x.angle_ideal,
      x.weight, x.slack, x.origin_id, x.sym_ops is None)
      for x in m.get_restraints_manager().geometry.angle_proxies
      if not nz_h & set(x.i_seqs)), p, nz_h, atoms))
  assert result[0][0] == result[1][0]
  p, nz_h, atoms = result[1][1:]
  assert sorted((atoms[i].name.strip(), atoms[i].parent().altloc)
    for i in nz_h) == [('HZ2', 'A'), ('HZ2', 'B'), ('HZ3', 'A'), ('HZ3', 'B')]
  assert len(p.metal_amine_converted_text) == 2, p.metal_amine_converted_text

def test_one_supplied_h():
  # one supplied H stays; the second is placed around it, or (supplied H
  # towards the metal) the group is completed as without the option
  run('tst_reduce2_metal_amine_one0', pdb_lys_cd, [ON])
  with open('tst_reduce2_metal_amine_one0_out.pdb') as f:
    nh2 = f.read()
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
    ['keep_existing_H=True'])
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
  h, description, log = run('tst_reduce2_metal_amine_nm', pdb_lys_cd, [ON])
  n, ce, cd = [matrix.col(atom(h, r, x).xyz)
    for r, x in ((359, 'NZ'), (359, 'CE'), (396, 'CD'))]
  side = (ce - n).cross(cd - n)
  assert (matrix.col(atom(h, 359, 'HZ2').xyz) - n).dot(side) > 0
  assert (matrix.col(atom(h, 359, 'HZ3').xyz) - n).dot(side) < 0
  with open('tst_reduce2_metal_amine_nm_out.pdb') as f:
    text = f.read()
  lines = text.splitlines()
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
  run('tst_reduce2_metal_amine_cf0', pdb_lys_cd)
  with open('tst_reduce2_metal_amine_cf0_out.pdb') as f:
    nh3 = f.read()
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
    off, d_off, l_off = run('tst_reduce2_metal_amine_cf_off', text, args)
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
    ['keep_existing_H=True'])
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

def lys_cd_alternates(names):
  '''pdb_lys_cd with the atoms names and the Cd in A and B at 0.5, B moved
  0.3 A.'''
  out = []
  for l in pdb_lys_cd.splitlines():
    if ((l.startswith('ATOM') and l[12:16].strip() in names) or
        (l.startswith('HETATM') and ' CD ' in l)):
      out.append(l[:16] + 'A' + l[17:54] + '  0.50' + l[60:])
      out.append(l[:16] + 'B' + l[17:30] + '%8.3f' % (float(l[30:38]) + 0.3) +
        l[38:54] + '  0.50' + l[60:])
    else:
      out.append(l)
  return '\n'.join(out) + '\n'

def test_alternate_metals_per_conformer():
  # NZ (or only CE) and the Cd in A and B: each conformer NH2, oriented to
  # its own Cd (a guard: passes without C1d)
  for tag, names in (('nz', ('CE', 'NZ')), ('ce', ('CE',))):
    h, d, l = run('tst_reduce2_metal_amine_alt_' + tag, lys_cd_alternates(
      names), [ON])
    assert len(section(d, CONVERTED)) == 2, (tag, d)
    a = dict(((x.name.strip(), x.parent().altloc), x) for x in h.atoms())
    for alt in 'AB':
      n = a.get(('NZ', alt)) or a[('NZ', '')]
      ce, cd = a[('CE', alt)], a[('CD', alt)]
      hs = [a[(x, alt)] for x in ('HZ2', 'HZ3')]
      assert not ('HZ1', alt) in a, (tag, alt)
      for x in hs:
        assert 104 < angle(cd, n, x) < 115, (tag, alt, angle(cd, n, x))
        assert 104 < angle(ce, n, x) < 115, (tag, alt, angle(ce, n, x))
        assert x.distance(cd) > 2.3, (tag, alt, x.distance(cd))
      assert 104 < angle(hs[0], n, hs[1]) < 115, (tag, alt)

def test_special_position():
  # a Cd on a 2-fold (and 0.04 A off it): one metal position, one contact
  # (a guard: cctbx merges symmetry-equivalent images closer than 0.5 A,
  # min_distance_sym_equiv; for a 2-fold that is 0.25 A from the axis)
  for tag, pdb in (('on', pdb_lys_cd_2fold), ('off', pdb_lys_cd_2fold.replace(
      'CD A 396       0.000', 'CD A 396       0.040'))):
    h, d, l = run('tst_reduce2_metal_amine_2fold_' + tag, pdb, [ON])
    lines = section(d, CONVERTED)
    assert len(lines) == 1 and lines[0].startswith('NZ  LYS A 359: HZ1 not '
      'placed  (NZ  LYS A 359 - CD    CD A 396 2.') and ';' not in lines[0], \
      (tag, lines)
    check_nh2(h, 359, 'NZ', 'CE', 396, 'CD', ['HZ2', 'HZ3'])

if __name__ == '__main__':
  test_lys_cd()
  test_option_off()
  test_terminus_cu()
  test_neutral_nh2_kept()
  test_chain_start_settings()
  test_chain_start_preceding_residue()
  test_chain_start_refused()
  test_chain_start_without_metal()
  test_chain_start_rerun()
  test_chain_start_option_off()
  test_chain_start_supplied_nh2_all_chain_starts()
  test_chain_start_supplied_h()
  test_chain_start_conformers()
  test_planar_nh2_guard()
  test_unresolved()
  test_supplied()
  test_rerun_stable()
  test_qmready_charge()
  test_hd_and_h()
  test_geometry_not_placeable()
  test_supplied_nh2_any_metal()
  test_only_reference_angles_removed()
  test_n_o_not_converted()
  test_lists_after_gc()
  test_other_h_unchanged()
  test_unequal_targets()
  test_symmetry_contact()
  test_shared_frame_restraints()
  test_one_supplied_h()
  test_names()
  test_nh3_across_conformers()
  test_alternate_metals_per_conformer()
  test_special_position()
  print('OK')
