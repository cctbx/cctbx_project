"""reduce2: the names of a propeller (CH3, NH3) whose torsion restraints tell
its H apart (A62b). The CCD fallback gives every H its own torsion of period
1 from the CCD's coordinates. Placement names such a propeller by those
torsions, and after the Optimizer, which turns a methyl as a symmetric group,
it is named again (sites unchanged). Period-3 propellers, all standard
residues among them, keep the handedness rule. Before: only the handedness,
by a swap of the second and third names, so with the riding sign fixed
(A62) 1NA's methyl names came out 120 deg off; master's default output had
1NA HM1-3 120 deg off too. Fixtures: 1NA heavy atoms from the CCD's ideal
coordinates; ALA, ILE, LEU, MET, THR, VAL from 1us0 and DT from 5gsk.

Supplied H/D pairs (keep_existing_H=True neutron_option=all_hd) keep their
occupancies and B factors at their sites when their names change. A
three-name rename done as two swaps, whose records read the names after
both, would send the values round the methyl's sites.
"""
from __future__ import absolute_import, division, print_function
import os
from itertools import permutations
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out
from six.moves import cStringIO as StringIO
from cctbx import geometry_restraints
from scitbx.math import dihedral_angle
from scitbx import matrix
from iotbx.cli_parser import run_program
from mmtbx.programs import reduce2
from mmtbx.hydrogens import reduce_hydrogen

def run(prefix, pdb_str, args=()):
  fn, out = prefix + '_in.pdb', prefix + '_out.pdb'
  with open(fn, 'w') as f:
    f.write(pdb_str)
  if os.path.exists(out): os.remove(out)
  run_program(program_class=reduce2.Program, logger=StringIO(), args=[fn,
    'output.filename=' + out, 'output.description_file_name=%s.txt' % prefix,
    'output.overwrite=True'] + list(args))
  return out

def torsions(fn):
  """(sites by name, {H name: [(atom names, ideal, period)]}) of a fresh
  model of fn with the restraints placement uses (CCD fallback included)"""
  m = mmtbx.model.manager(model_input=iotbx.pdb.input(fn), log=null_out())
  reduce_hydrogen.add_ccd_restraint_objects(m)
  m.process(make_restraints=True)
  g = m.get_restraints_manager().geometry
  atoms = m.get_hierarchy().atoms()
  names = [a.name.strip() for a in atoms]
  result = {}
  for kind, proxies in (('var', g.dihedral_proxies),
                        ('const', g.const_dihedral_proxies or [])):
    for p in proxies:
      for ih in (p.i_seqs[0], p.i_seqs[3]):
        if atoms[ih].element.strip().upper() not in ('H', 'D'): continue
        result.setdefault(names[ih], []).append((tuple(names[i]
          for i in p.i_seqs), p.angle_ideal,
          p.periodicity if kind == 'var' else 1))
  return dict((n, a.xyz) for n, a in zip(names, atoms)), result

def residual(sites, t, h, site):
  names, ideal, period = t
  d = dihedral_angle(sites=[site if n == h else sites[n] for n in names],
    deg=True)
  return abs(geometry_restraints.angle_delta_deg(angle_1=d, angle_2=ideal,
    periodicity=max(1, period)))

def check_methyl(sites, tors, hs, tol=None):
  """the names hs are the permutation over their sites that meets their
  torsions best; with tol, each H within tol of its torsions"""
  def cost(perm):
    return sum(residual(sites, t, h, sites[s]) for h, s in zip(hs, perm)
      for t in tors[h])
  best = min(cost(p) for p in permutations(hs))
  assert cost(hs) <= best + 1.e-6, (hs, cost(hs), best)
  if tol is not None:
    for h in hs:
      r = max(residual(sites, t, h, sites[h]) for t in tors[h])
      assert r < tol, (h, r)

# ------------------------------------------------------------------------------

def test_place():
  # approach=place: both methyls of 1NA meet their per-H torsions by name
  sites, tors = torsions(run('tst_reduce2_propeller_names_place', pdb_1na,
    ['approach=place']))
  check_methyl(sites, tors, ('HM1', 'HM2', 'HM3'), tol=1.)
  check_methyl(sites, tors, ('H81', 'H82', 'H83'), tol=1.)

def test_add():
  # default (Optimizer): CM is staggered at its torsions, HM1-3 meet them by
  # name (master: 120 deg off); C8 (an acetyl methyl) is turned by the
  # Optimizer, and its names are the best over its sites
  sites, tors = torsions(run('tst_reduce2_propeller_names_add', pdb_1na))
  check_methyl(sites, tors, ('HM1', 'HM2', 'HM3'), tol=1.)
  check_methyl(sites, tors, ('H81', 'H82', 'H83'))

def test_standard():
  # period-3 methyls of standard residues: every H line as master wrote it
  out = run('tst_reduce2_propeller_names_std', pdb_std)
  lines = [l for l in open(out).read().splitlines()
    if l.startswith(('ATOM', 'HETATM')) and l[76:78].strip() == 'H']
  answer = answer_std.strip().splitlines()
  assert len(lines) == len(answer), (len(lines), len(answer))
  diff = [(a, b) for a, b in zip(answer, lines) if a != b]
  assert not diff, '\n'.join(a + '\n' + b for a, b in diff[:10])

def supplied_pairs(prefix, args, sites, names):
  """1NA with an H/D pair (altlocs A, B) at each of sites, the pair at site
  k named names[k] (H, D) with D occupancy 0.1 (k+1), H 1 - that, B factor
  11 (k+1). Returns, per site, the output's [(name, element, occ, b)] within
  0.2 A of it (neutron distances move the H along their bonds)."""
  lines = pdb_1na.strip().splitlines()[:-1]
  for k, (xyz, (h, d)) in enumerate(zip(sites, names)):
    occ_d, b = 0.1 * (k + 1), 11. * (k + 1)
    for name, alt, occ, e in ((h, 'A', 1 - occ_d, 'H'), (d, 'B', occ_d, 'D')):
      lines.append('HETATM%5d %-4s%s1NA A   1    %8.3f%8.3f%8.3f%6.2f%6.2f'
        '          %2s' % (30 + len(lines), ' ' + name, alt, xyz[0], xyz[1],
        xyz[2], occ, b, e))
  out = run(prefix, '\n'.join(lines) + '\nEND\n', ['keep_existing_H=True',
    'neutron_option=all_hd'] + list(args))
  atoms = iotbx.pdb.input(out).construct_hierarchy().atoms()
  result = []
  for xyz in sites:
    result.append(sorted((a.name.strip(), a.element.strip(), round(a.occ, 2),
      round(a.b, 2)) for a in atoms
      if (matrix.col(a.xyz) - matrix.col(xyz)).length() < 0.2))
  return result

def check_values(result, pairs=(0, 1, 2)):
  # site k holds the values of input pair pairs[k], whatever their names now
  for k, at_site in enumerate(result):
    occ_d, b = round(0.1 * (pairs[k] + 1), 2), round(11. * (pairs[k] + 1), 2)
    assert [(e, o, bb) for n, e, o, bb in at_site] == [
      ('D', occ_d, b), ('H', round(1 - occ_d, 2), b)], (k, at_site)

# the sites of 1NA's HM1, HM2, HM3: C1-O1-CM-H 180, -60, 60 deg
sites_cm = [(18.938, 14.860, 15.869), (18.142, 13.668, 15.194),
            (18.246, 15.067, 14.458)]

def test_supplied_cycle_place():
  # names cyclically permuted: placement names the methyl back by its
  # torsions (a three-name cycle), the values stay at their sites
  result = supplied_pairs('tst_reduce2_propeller_names_cycle', [
    'approach=place'], sites_cm, [('HM2', 'DM2'), ('HM3', 'DM3'),
    ('HM1', 'DM1')])
  check_values(result)
  assert [[n for n, e, o, b in r] for r in result] == [['DM1', 'HM1'],
    ['DM2', 'HM2'], ['DM3', 'HM3']], result

def test_supplied_add():
  # correctly named pairs, Optimizer: it staggers the methyl by a 120 deg
  # turn of the supplied atoms, which take their values along (as on master:
  # site k holds pair k+1). The names given after it by the torsions move no
  # values (as two swaps: one site further round).
  result = supplied_pairs('tst_reduce2_propeller_names_pairs_add', [],
    sites_cm, [('HM1', 'DM1'), ('HM2', 'DM2'), ('HM3', 'DM3')])
  check_values(result, pairs=(1, 2, 0))
  assert [[n for n, e, o, b in r] for r in result] == [['DM1', 'HM1'],
    ['DM2', 'HM2'], ['DM3', 'HM3']], result

def test_supplied_swap():
  # a CH2 with its names swapped: the handedness rule swaps them back (two
  # names, as on master), the values stay at their sites
  result = supplied_pairs('tst_reduce2_propeller_names_swap', [
    'approach=place'], [(14.607, 14.961, 10.968), (15.237, 16.289, 11.557)],
    [('H62', 'D62'), ('H61', 'D61')])
  check_values(result)
  assert [[n for n, e, o, b in r] for r in result] == [['D61', 'H61'],
    ['D62', 'H62']], result

# ------------------------------------------------------------------------------

pdb_1na = """
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
HETATM    1  C1  1NA A   1      15.887  14.694  15.119  1.00 20.00           C
HETATM    2  C2  1NA A   1      14.569  15.126  15.764  1.00 20.00           C
HETATM    3  C3  1NA A   1      13.407  14.697  14.862  1.00 20.00           C
HETATM    4  C4  1NA A   1      13.654  15.245  13.452  1.00 20.00           C
HETATM    5  C5  1NA A   1      15.042  14.805  12.982  1.00 20.00           C
HETATM    6  C6  1NA A   1      15.284  15.320  11.562  1.00 20.00           C
HETATM    7  C7  1NA A   1      14.902  15.107  18.178  1.00 20.00           C
HETATM    8  C8  1NA A   1      14.761  14.452  19.527  1.00 20.00           C
HETATM    9  CM  1NA A   1      18.172  14.629  15.320  1.00 20.00           C
HETATM   10  N2  1NA A   1      14.432  14.490  17.076  1.00 20.00           N
HETATM   11  O1  1NA A   1      16.975  15.060  15.970  1.00 20.00           O
HETATM   12  O3  1NA A   1      12.180  15.222  15.371  1.00 20.00           O
HETATM   13  O4  1NA A   1      12.663  14.737  12.558  1.00 20.00           O
HETATM   14  O5  1NA A   1      16.038  15.333  13.855  1.00 20.00           O
HETATM   15  O6  1NA A   1      16.576  14.904  11.118  1.00 20.00           O
HETATM   16  O7  1NA A   1      15.437  16.191  18.083  1.00 20.00           O
END
"""

pdb_std = """
CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1
ATOM     18  N   ALA A   1      15.555  12.152  39.893  1.00 20.00           N
ATOM     19  CA  ALA A   1      14.153  11.844  39.611  1.00 20.00           C
ATOM     20  C   ALA A   1      14.014  10.631  38.728  1.00 20.00           C
ATOM     21  O   ALA A   1      14.644  10.406  37.752  1.00 20.00           O
ATOM     22  CB  ALA A   1      13.343  12.950  38.894  1.00 20.00           C
ATOM     63  N   ILE A   4      10.519  10.546  33.584  1.00 20.00           N
ATOM     64  CA  ILE A   4      10.571  11.633  32.619  1.00 20.00           C
ATOM     65  C   ILE A   4       9.172  11.830  32.066  1.00 20.00           C
ATOM     66  O   ILE A   4       8.408  10.887  31.787  1.00 20.00           O
ATOM     67  CB  ILE A   4      11.629  11.351  31.519  1.00 20.00           C
ATOM     68  CG1 ILE A   4      11.902  12.558  30.625  1.00 20.00           C
ATOM     69  CG2 ILE A   4      11.309  10.131  30.694  1.00 20.00           C
ATOM     70  CD1 ILE A   4      13.175  12.395  29.776  1.00 20.00           C
ATOM     82  N   LEU A   5       8.793  13.088  31.859  1.00 20.00           N
ATOM     83  CA  LEU A   5       7.495  13.415  31.309  1.00 20.00           C
ATOM     84  C   LEU A   5       7.475  13.258  29.812  1.00 20.00           C
ATOM     85  O   LEU A   5       8.283  13.868  29.108  1.00 20.00           O
ATOM     86  CB  LEU A   5       7.178  14.879  31.691  1.00 20.00           C
ATOM     87  CG  LEU A   5       5.778  15.362  31.308  1.00 20.00           C
ATOM     88  CD1 LEU A   5       4.608  14.558  31.862  1.00 20.00           C
ATOM     89  CD2 LEU A   5       5.657  16.811  31.757  1.00 20.00           C
ATOM      1  N   MET A   0      18.709  11.104  41.491  1.00 20.00           N
ATOM      2  CA  MET A   0      17.937  11.685  40.368  1.00 20.00           C
ATOM      3  C   MET A   0      16.479  11.264  40.309  1.00 20.00           C
ATOM      4  O   MET A   0      16.161  10.141  40.705  1.00 20.00           O
ATOM      5  CB  MET A   0      18.672  11.229  39.111  1.00 20.00           C
ATOM      6  CG  MET A   0      18.004  11.534  37.807  1.00 20.00           C
ATOM      7  SD  MET A   0      18.852  11.085  36.316  1.00 20.00           S
ATOM      8  CE  MET A   0      20.086  12.267  35.949  1.00 20.00           C
ATOM    321  N   THR A  19      14.553  -7.200  27.782  1.00 20.00           N
ATOM    322  CA  THR A  19      14.176  -8.603  27.740  1.00 20.00           C
ATOM    323  C   THR A  19      14.578  -9.324  26.462  1.00 20.00           C
ATOM    324  O   THR A  19      14.284 -10.524  26.332  1.00 20.00           O
ATOM    325  CB  THR A  19      12.675  -8.802  27.996  1.00 20.00           C
ATOM    326  OG1 THR A  19      11.927  -8.210  26.942  1.00 20.00           O
ATOM    327  CG2 THR A  19      12.241  -8.267  29.356  1.00 20.00           C
ATOM    452  N   VAL A  27      13.214 -15.629  32.266  1.00 20.00           N
ATOM    453  CA  VAL A  27      12.767 -14.249  32.254  1.00 20.00           C
ATOM    454  C   VAL A  27      11.641 -14.003  33.257  1.00 20.00           C
ATOM    455  O   VAL A  27      11.487 -12.880  33.732  1.00 20.00           O
ATOM    456  CB  VAL A  27      12.357 -13.813  30.831  1.00 20.00           C
ATOM    457  CG1 VAL A  27      11.059 -14.438  30.367  1.00 20.00           C
ATOM    458  CG2 VAL A  27      12.301 -12.293  30.740  1.00 20.00           C
ATOM     65  P    DT B   3      10.011  19.045   4.108  1.00 20.00           P
ATOM     66  OP1  DT B   3      10.235  19.681   5.425  1.00 20.00           O
ATOM     67  OP2  DT B   3      10.906  17.959   3.643  1.00 20.00           O
ATOM     68  O5'  DT B   3       8.532  18.470   4.070  1.00 20.00           O
ATOM     69  C5'  DT B   3       7.452  19.242   4.612  1.00 20.00           C
ATOM     70  C4'  DT B   3       6.153  18.502   4.416  1.00 20.00           C
ATOM     71  O4'  DT B   3       5.899  18.337   3.006  1.00 20.00           O
ATOM     72  C3'  DT B   3       6.130  17.086   4.976  1.00 20.00           C
ATOM     73  O3'  DT B   3       5.798  17.096   6.359  1.00 20.00           O
ATOM     74  C2'  DT B   3       5.039  16.439   4.154  1.00 20.00           C
ATOM     75  C1'  DT B   3       5.228  17.096   2.794  1.00 20.00           C
ATOM     76  N1   DT B   3       6.037  16.289   1.874  1.00 20.00           N
ATOM     77  C2   DT B   3       5.382  15.302   1.178  1.00 20.00           C
ATOM     78  O2   DT B   3       4.185  15.095   1.285  1.00 20.00           O
ATOM     79  N3   DT B   3       6.183  14.567   0.341  1.00 20.00           N
ATOM     80  C4   DT B   3       7.540  14.721   0.135  1.00 20.00           C
ATOM     81  O4   DT B   3       8.121  13.998  -0.668  1.00 20.00           O
ATOM     82  C5   DT B   3       8.169  15.755   0.922  1.00 20.00           C
ATOM     83  C7   DT B   3       9.642  15.978   0.784  1.00 20.00           C
ATOM     84  C6   DT B   3       7.398  16.472   1.749  1.00 20.00           C
END
"""

answer_std = """
ATOM      6  H   ALA A   1      15.810  12.966  39.784  1.00 20.00           H
ATOM      7  HA  ALA A   1      13.782  11.722  40.499  1.00 20.00           H
ATOM      8  HB1 ALA A   1      13.758  13.141  38.038  1.00 20.00           H
ATOM      9  HB2 ALA A   1      12.434  12.639  38.759  1.00 20.00           H
ATOM     10  HB3 ALA A   1      13.341  13.748  39.446  1.00 20.00           H
ATOM     19  H   ILE A   4      11.110   9.935  33.455  1.00 20.00           H
ATOM     20  HA  ILE A   4      10.851  12.462  33.038  1.00 20.00           H
ATOM     21  HB  ILE A   4      12.446  11.166  32.008  1.00 20.00           H
ATOM     22 HG12 ILE A   4      11.152  12.680  30.022  1.00 20.00           H
ATOM     23 HG13 ILE A   4      12.011  13.345  31.182  1.00 20.00           H
ATOM     24 HG21 ILE A   4      10.452  10.260  30.258  1.00 20.00           H
ATOM     25 HG22 ILE A   4      12.004  10.009  30.028  1.00 20.00           H
ATOM     26 HG23 ILE A   4      11.270   9.357  31.277  1.00 20.00           H
ATOM     27 HD11 ILE A   4      13.077  11.616  29.207  1.00 20.00           H
ATOM     28 HD12 ILE A   4      13.295  13.189  29.232  1.00 20.00           H
ATOM     29 HD13 ILE A   4      13.936  12.280  30.367  1.00 20.00           H
ATOM     38  H   LEU A   5       9.280  13.775  32.033  1.00 20.00           H
ATOM     39  HA  LEU A   5       6.816  12.821  31.664  1.00 20.00           H
ATOM     40  HB2 LEU A   5       7.264  14.969  32.653  1.00 20.00           H
ATOM     41  HB3 LEU A   5       7.817  15.457  31.245  1.00 20.00           H
ATOM     42  HG  LEU A   5       5.697  15.251  30.348  1.00 20.00           H
ATOM     43 HD11 LEU A   5       4.644  14.575  32.831  1.00 20.00           H
ATOM     44 HD12 LEU A   5       3.778  14.955  31.555  1.00 20.00           H
ATOM     45 HD13 LEU A   5       4.675  13.644  31.545  1.00 20.00           H
ATOM     46 HD21 LEU A   5       6.336  17.338  31.308  1.00 20.00           H
ATOM     47 HD22 LEU A   5       4.775  17.141  31.525  1.00 20.00           H
ATOM     48 HD23 LEU A   5       5.784  16.856  32.718  1.00 20.00           H
ATOM     57  H   MET A   0      18.329  11.190  42.258  1.00 20.00           H
ATOM     58  HA  MET A   0      17.900  12.650  40.455  1.00 20.00           H
ATOM     59  HB2 MET A   0      19.540  11.661  39.095  1.00 20.00           H
ATOM     60  HB3 MET A   0      18.781  10.266  39.157  1.00 20.00           H
ATOM     61  HG2 MET A   0      17.153  11.068  37.794  1.00 20.00           H
ATOM     62  HG3 MET A   0      17.858  12.492  37.766  1.00 20.00           H
ATOM     63  HE1 MET A   0      20.720  12.300  36.683  1.00 20.00           H
ATOM     64  HE2 MET A   0      20.539  12.007  35.132  1.00 20.00           H
ATOM     65  HE3 MET A   0      19.668  13.135  35.833  1.00 20.00           H
ATOM     73  H   THR A  19      13.962  -6.669  27.454  1.00 20.00           H
ATOM     74  HA  THR A  19      14.672  -9.044  28.448  1.00 20.00           H
ATOM     75  HB  THR A  19      12.483  -9.753  28.017  1.00 20.00           H
ATOM     76  HG1 THR A  19      11.973  -7.373  26.992  1.00 20.00           H
ATOM     77 HG21 THR A  19      12.426  -7.317  29.413  1.00 20.00           H
ATOM     78 HG22 THR A  19      11.290  -8.410  29.481  1.00 20.00           H
ATOM     79 HG23 THR A  19      12.724  -8.725  30.062  1.00 20.00           H
ATOM     87  H   VAL A  27      12.708 -16.164  31.822  1.00 20.00           H
ATOM     88  HA  VAL A  27      13.511 -13.694  32.536  1.00 20.00           H
ATOM     89  HB  VAL A  27      13.037 -14.139  30.221  1.00 20.00           H
ATOM     90 HG11 VAL A  27      10.348 -14.178  30.974  1.00 20.00           H
ATOM     91 HG12 VAL A  27      10.860 -14.124  29.471  1.00 20.00           H
ATOM     92 HG13 VAL A  27      11.156 -15.403  30.366  1.00 20.00           H
ATOM     93 HG21 VAL A  27      13.177 -11.932  30.946  1.00 20.00           H
ATOM     94 HG22 VAL A  27      12.042 -12.040  29.840  1.00 20.00           H
ATOM     95 HG23 VAL A  27      11.649 -11.962  31.377  1.00 20.00           H
ATOM    116  H5'  DT B   3       7.602  19.390   5.559  1.00 20.00           H
ATOM    117 H5''  DT B   3       7.408  20.100   4.162  1.00 20.00           H
ATOM    118  H4'  DT B   3       5.447  19.025   4.826  1.00 20.00           H
ATOM    119  H3'  DT B   3       6.978  16.632   4.852  1.00 20.00           H
ATOM    120 HO3'  DT B   3       6.103  17.697   6.860  1.00 20.00           H
ATOM    121  H2'  DT B   3       5.158  15.478   4.102  1.00 20.00           H
ATOM    122 H2''  DT B   3       4.161  16.629   4.520  1.00 20.00           H
ATOM    123  H1'  DT B   3       4.360  17.291   2.407  1.00 20.00           H
ATOM    124  H3   DT B   3       5.797  13.941  -0.105  1.00 20.00           H
ATOM    125 H5M1  DT B   3      10.069  15.138   0.555  1.00 20.00           H
ATOM    126 H5M2  DT B   3       9.992  16.307   1.627  1.00 20.00           H
ATOM    127 H5M3  DT B   3       9.798  16.630   0.083  1.00 20.00           H
ATOM    128  H6   DT B   3       7.803  17.130   2.268  1.00 20.00           H
"""

if __name__ == '__main__':
  test_place()
  test_add()
  test_standard()
  test_supplied_cycle_place()
  test_supplied_add()
  test_supplied_swap()
  print('OK')
