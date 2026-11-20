"""Fixtures, run helpers and a per-process cache of shared reduce2 runs for
tst_reduce2_metal_amine.py and tst_reduce2_metal_amine_supplied.py
(metal_bound_amines_as_nh2, A5). Not a test file."""
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

# 1MNZ Tris 392 N and O2 on Mg 389 (2.26 and 2.31 A)
pdb_trs_mg = """\
CRYST1   93.880   99.680  102.900  90.00  90.00  90.00 I 2 2 2
HETATM    1  C   TRS A 392      33.735  32.634  58.095  1.00 12.66           C
HETATM    2  C1  TRS A 392      34.731  31.531  58.435  1.00 18.88           C
HETATM    3  C2  TRS A 392      34.369  33.522  57.041  1.00 14.10           C
HETATM    4  C3  TRS A 392      32.500  31.852  57.652  1.00 12.06           C
HETATM    5  N   TRS A 392      33.373  33.550  59.189  1.00 10.61           N
HETATM    6  O1  TRS A 392      36.136  31.939  58.199  1.00 33.72           O
HETATM    7  O2  TRS A 392      35.424  34.387  57.627  1.00 12.13           O
HETATM    8  O3  TRS A 392      32.743  31.018  56.515  1.00 11.60           O
HETATM    9 MG    MG A 389      34.719  35.310  59.621  1.00  8.72          MG
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
OFF = 'metal_bound_amines_as_nh2=False'

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

# one reduce2 run per fixture and options for the whole file; the tests only
# read the results. x1_off also runs with the lists built after a forced
# gc.collect() (test_lists_after_gc)
SHARED_RUNS = {
  'lys_on':  ('pdb_lys_cd', [ON]),
  'lys_off': ('pdb_lys_cd', [OFF]),
  'x1_on':   ('pdb_1xmk', [ON]),
  'x1_off':  ('pdb_1xmk', [OFF]),
}
_shared = {}

def shared(key):
  '''(hierarchy, description, log) of SHARED_RUNS[key], run once.'''
  if key not in _shared:
    pdb, args = SHARED_RUNS[key]
    finish = None
    if key == 'x1_off':
      import gc
      from mmtbx.hydrogens import reduce_hydrogen
      finish = reduce_hydrogen.place_hydrogens._finish_metal_lists
      def collected_first(self):
        gc.collect()
        return finish(self)
      reduce_hydrogen.place_hydrogens._finish_metal_lists = collected_first
    try:
      _shared[key] = run('tst_reduce2_metal_amine_' + key, globals()[pdb], args)
    finally:
      if finish is not None:
        reduce_hydrogen.place_hydrogens._finish_metal_lists = finish
  return _shared[key]

def shared_text(key):
  '''The written model of shared(key).'''
  shared(key)
  with open('tst_reduce2_metal_amine_%s_out.pdb' % key) as f:
    return f.read()

def processed(pdb_str):
  '''A processed model of pdb_str, with restraints (no reduce2 run).'''
  import mmtbx.model
  from libtbx.utils import null_out
  m = mmtbx.model.manager(model_input = iotbx.pdb.input(source_info=None,
    lines=pdb_str), log = null_out())
  m.process(make_restraints = True)
  return m

def lys_targets(**kw):
  '''The geostd LYS targets of HZ2/HZ3 (X-ray N-H), for place_metal_nh2().'''
  from libtbx import group_args
  t = dict(d = [0.89, 0.89], cnh = [110., 110.], hnh = 109.)
  t.update(kw)
  return group_args(**t)

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

def classifier_model():
  '''EDN + Cu, Asn + Zn, YI1 + Zn and Lys + Cd in one processed model, 40 A
  apart (chains A-D), with the metal_amine_sites() of each N.'''
  out = ['CRYST1  200.000  200.000  200.000  90.00  90.00  90.00 P 1']
  for k, (chain, pdb) in enumerate((('A', pdb_edn_cu), ('B', pdb_asn_zn),
                                    ('C', pdb_yi1_zn), ('D', pdb_lys_cd))):
    for l in pdb.splitlines():
      if not l.startswith(('ATOM', 'HETATM')): continue
      x = float(l[30:38]) + 40 * k
      out.append(l[:21] + chain + l[22:30] + '%8.3f' % x + l[38:])
  return processed('\n'.join(out) + '\nEND\n')

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
