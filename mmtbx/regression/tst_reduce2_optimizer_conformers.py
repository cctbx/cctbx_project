"""reduce2 Optimizer: each altloc pass uses the bonded neighbours of its own
conformer. Placement copies the H of a residue with altlocs into each altloc
group, so a blank parent carries H in A and in B; with neighbour lists built
from all atoms the parent counted both copies (and atoms that exist only in
another conformer) and got no rotator, NH3 rotator or methyl stagger.
Models: the Zn-Cys model of tst_add_hydrogen_1 with O split A/B and the Zn in
altloc A (A4 B0 case c2); 1US0 Tyr A 189, Lys A 85 and Ile A 14 (backbone
atoms split A/B, side chains blank); 3CP5 Thr A 46 (blank CB, OG1 and CG2
split A/B: CB counted both conformers' substituents).
A residue present only in another conformer (no blank atoms, one other altloc)
in a chain without the pass's altloc is left out of that pass:
the c2 Zn after TER; 4Y7E BMA A 414 (altloc A) and its superposed copy BMA C 1
(altloc B), each in its own chain; three altlocs; alt_id given.
"""
from __future__ import absolute_import, division, print_function
import os, re
from libtbx.utils import null_out
from iotbx.cli_parser import run_program
import iotbx.pdb
import scitbx.math
from mmtbx.programs import reduce2
from mmtbx.reduce import Optimizers

m_c2 = """\
CRYST1  108.910  108.910  108.910  90.00  90.00  90.00 I 4 3 2
ATOM      1  N   CYS A  48     -40.561 -16.057   6.830  1.00 99.72           N
ATOM      2  CA  CYS A  48     -41.351 -14.915   6.380  1.00106.82           C
ATOM      3  C   CYS A  48     -41.146 -13.719   7.305  1.00104.36           C
ATOM      4  CB  CYS A  48     -42.835 -15.275   6.303  1.00105.79           C
ATOM      5  SG  CYS A  48     -43.707 -15.186   7.879  1.00101.00           S
ATOM      6  H   CYS A  48     -39.952 -15.857   7.403  1.00 99.72           H
ATOM      7  HA  CYS A  48     -41.063 -14.664   5.488  1.00106.82           H
ATOM      8  HB2 CYS A  48     -43.271 -14.661   5.692  1.00105.79           H
ATOM      9  HB3 CYS A  48     -42.916 -16.184   5.974  1.00105.79           H
ATOM     10  O  ACYS A  48     -40.407 -13.794   8.283  0.50104.56           O
ATOM     11  O  BCYS A  48     -40.107 -13.794   8.283  0.50104.56           O
HETATM   12 ZN  A ZN A 102     -44.322 -17.446   8.351  0.50 98.67          Zn
"""

m_c2_ter = """\
CRYST1  108.910  108.910  108.910  90.00  90.00  90.00 I 4 3 2
ATOM      1  N   CYS A  48     -40.561 -16.057   6.830  1.00 99.72           N
ATOM      2  CA  CYS A  48     -41.351 -14.915   6.380  1.00106.82           C
ATOM      3  C   CYS A  48     -41.146 -13.719   7.305  1.00104.36           C
ATOM      4  CB  CYS A  48     -42.835 -15.275   6.303  1.00105.79           C
ATOM      5  SG  CYS A  48     -43.707 -15.186   7.879  1.00101.00           S
ATOM      6  H   CYS A  48     -39.952 -15.857   7.403  1.00 99.72           H
ATOM      7  HA  CYS A  48     -41.063 -14.664   5.488  1.00106.82           H
ATOM      8  HB2 CYS A  48     -43.271 -14.661   5.692  1.00105.79           H
ATOM      9  HB3 CYS A  48     -42.916 -16.184   5.974  1.00105.79           H
ATOM     10  O  ACYS A  48     -40.407 -13.794   8.283  0.50104.56           O
ATOM     11  O  BCYS A  48     -40.107 -13.794   8.283  0.50104.56           O
TER
HETATM   12 ZN  A ZN A 102     -44.322 -17.446   8.351  0.50 98.67          Zn
"""

m_tyr = """\
CRYST1   49.280   66.590   47.260  90.00  92.40  90.00 P 1 21 1
ATOM      1  CB  TYR A 189      31.161   2.318  20.873  1.00  4.91           C
ATOM      2  CG  TYR A 189      32.212   1.832  21.814  1.00  4.64           C
ATOM      3  CD1 TYR A 189      33.153   2.656  22.427  1.00  5.01           C
ATOM      4  CD2 TYR A 189      32.298   0.465  22.082  1.00  5.33           C
ATOM      5  CE1 TYR A 189      34.117   2.143  23.284  1.00  5.53           C
ATOM      6  CE2 TYR A 189      33.257  -0.053  22.913  1.00  5.78           C
ATOM      7  CZ  TYR A 189      34.192   0.781  23.527  1.00  5.71           C
ATOM      8  OH  TYR A 189      35.102   0.214  24.322  1.00  7.60           O
ATOM      9  N  ATYR A 189      30.179   3.957  22.382  0.70  3.81           N
ATOM     10  CA ATYR A 189      30.751   3.760  21.088  0.70  4.26           C
ATOM     11  C  ATYR A 189      29.859   4.180  19.894  0.70  4.04           C
ATOM     12  O  ATYR A 189      30.312   4.189  18.749  0.70  4.60           O
ATOM     13  N  BTYR A 189      30.072   4.022  22.367  0.30  4.47           N
ATOM     14  CA BTYR A 189      30.638   3.770  21.034  0.30  4.24           C
ATOM     15  C  BTYR A 189      29.628   4.269  19.967  0.30  4.19           C
ATOM     16  O  BTYR A 189      30.000   4.201  18.798  0.30  5.01           O
"""

m_lys = """\
CRYST1   49.280   66.590   47.260  90.00  92.40  90.00 P 1 21 1
ATOM      1  CA  LYS A  85      -3.766  -1.696   8.490  1.00  9.38           C
ATOM      2  CB  LYS A  85      -4.715  -0.995   7.464  1.00 11.25           C
ATOM      3  CG  LYS A  85      -5.987  -0.493   8.142  1.00 13.75           C
ATOM      4  CD  LYS A  85      -6.775   0.362   7.150  1.00 17.82           C
ATOM      5  CE  LYS A  85      -8.114   0.815   7.747  1.00 27.26           C
ATOM      6  NZ  LYS A  85      -8.971   1.246   6.577  1.00 49.79           N
ATOM      7  N  ALYS A  85      -2.477  -2.019   7.903  0.56  9.64           N
ATOM      8  C  ALYS A  85      -4.452  -2.813   9.223  0.56  8.85           C
ATOM      9  O  ALYS A  85      -4.798  -2.658  10.425  0.56  9.36           O
ATOM     10  N  BLYS A  85      -2.647  -2.343   7.781  0.44  7.16           N
ATOM     11  C  BLYS A  85      -4.579  -2.735   9.230  0.44  7.40           C
ATOM     12  O  BLYS A  85      -4.987  -2.510  10.383  0.44  8.51           O
"""

m_ile = """\
CRYST1   49.280   66.590   47.260  90.00  92.40  90.00 P 1 21 1
ATOM      1  N   ILE A  14      12.838   6.233  33.912  1.00  4.10           N
ATOM      2  CA  ILE A  14      14.047   6.362  33.112  1.00  4.43           C
ATOM      3  CB  ILE A  14      14.849   7.614  33.452  1.00  5.25           C
ATOM      4  CG1 ILE A  14      15.485   7.510  34.839  1.00  5.78           C
ATOM      5  CG2 ILE A  14      13.998   8.853  33.256  1.00  6.22           C
ATOM      6  CD1 ILE A  14      16.350   8.738  35.140  1.00  6.90           C
ATOM      7  C  AILE A  14      14.936   5.140  33.040  0.65  4.01           C
ATOM      8  O  AILE A  14      15.829   5.092  32.213  0.65  4.16           O
ATOM      9  C  BILE A  14      14.843   5.050  33.236  0.35  3.74           C
ATOM     10  O  BILE A  14      15.658   4.832  32.323  0.35  5.01           O
"""

m_thr = """\
CRYST1   39.592   32.644   42.064  90.00  94.56  90.00 P 1 21 1
ATOM    283  N   CYS A  45      -0.599  20.876  36.277  1.00  7.63           N
ATOM    284  CA  CYS A  45       0.525  20.659  35.381  1.00  7.94           C
ATOM    285  C   CYS A  45       0.890  21.874  34.534  1.00  7.46           C
ATOM    286  O   CYS A  45       2.101  22.043  34.206  1.00  8.28           O
ATOM    287  CB  CYS A  45       0.229  19.487  34.422  1.00  7.87           C
ATOM    288  SG  CYS A  45      -0.557  18.027  35.141  1.00  8.15           S
ATOM    289  N   THR A  46      -0.035  22.735  34.171  1.00  8.39           N
ATOM    290  CA  THR A  46       0.151  23.886  33.332  1.00  9.63           C
ATOM    291  C   THR A  46       1.093  24.917  33.930  1.00 11.11           C
ATOM    292  O   THR A  46       1.675  25.762  33.289  1.00 13.68           O
ATOM    293  CB  THR A  46      -1.217  24.509  32.944  1.00 14.06           C
ATOM    294  OG1ATHR A  46      -1.783  23.647  31.940  0.71 25.52           O
ATOM    295  OG1BTHR A  46      -2.100  24.715  34.036  0.29 16.93           O
ATOM    296  CG2ATHR A  46      -1.071  25.899  32.389  0.71 18.52           C
ATOM    297  CG2BTHR A  46      -1.909  23.534  31.986  0.29  6.58           C
ATOM    298  N   ALA A  47       1.305  24.909  35.232  1.00 10.58           N
ATOM    299  CA  ALA A  47       2.252  25.818  35.845  1.00 12.13           C
ATOM    300  C   ALA A  47       3.657  25.590  35.289  1.00 11.98           C
ATOM    301  O   ALA A  47       4.468  26.498  35.305  1.00 21.55           O
ATOM    302  CB  ALA A  47       2.251  25.606  37.352  1.00 14.76           C
"""

m_bma = """\
CRYST1   65.709  100.706  104.737  90.00  90.00  90.00 P 21 21 21
HETATM 4814  C1 ABMA A 414      58.277  85.591 -21.661  0.50 34.39           C
HETATM 4815  C2 ABMA A 414      59.360  84.534 -21.848  0.50 33.10           C
HETATM 4816  C3 ABMA A 414      58.795  83.150 -21.558  0.50 31.09           C
HETATM 4817  C4 ABMA A 414      58.059  83.120 -20.225  0.50 29.83           C
HETATM 4818  C5 ABMA A 414      57.075  84.278 -20.119  0.50 31.35           C
HETATM 4819  C6 ABMA A 414      56.337  84.300 -18.780  0.50 31.27           C
HETATM 4820  O1 ABMA A 414      58.807  86.901 -21.893  0.50 37.20           O
HETATM 4821  O2 ABMA A 414      60.438  84.791 -20.942  0.50 33.70           O
HETATM 4822  O3 ABMA A 414      59.851  82.191 -21.524  0.50 31.12           O
HETATM 4823  O4 ABMA A 414      57.299  81.913 -20.193  0.50 26.68           O
HETATM 4824  O5 ABMA A 414      57.760  85.510 -20.329  0.50 32.83           O
HETATM 4825  O6 ABMA A 414      57.234  84.116 -17.675  0.50 31.10           O
TER
HETATM 4645  C1 BBMA C   1      58.274  85.577 -21.658  0.50 34.39           C
HETATM 4646  C2 BBMA C   1      59.360  84.522 -21.844  0.50 33.10           C
HETATM 4647  C3 BBMA C   1      58.802  83.137 -21.549  0.50 31.08           C
HETATM 4648  C4 BBMA C   1      58.077  83.105 -20.210  0.50 29.83           C
HETATM 4649  C5 BBMA C   1      57.082  84.255 -20.113  0.50 31.35           C
HETATM 4650  C6 BBMA C   1      56.338  84.278 -18.776  0.50 31.27           C
HETATM 4651  O1 BBMA C   1      58.800  86.888 -21.887  0.50 37.21           O
HETATM 4652  O2 BBMA C   1      60.439  84.787 -20.942  0.50 33.70           O
HETATM 4653  O3 BBMA C   1      59.859  82.179 -21.524  0.50 31.16           O
HETATM 4654  O4 BBMA C   1      57.334  81.888 -20.166  0.50 26.65           O
HETATM 4655  O5 BBMA C   1      57.754  85.493 -20.327  0.50 32.82           O
HETATM 4656  O6 BBMA C   1      57.234  84.114 -17.668  0.50 31.08           O
"""

m_three = """\
ATOM      1  N  AALA A   1      10.000  10.000  10.000  0.40 10.00           N
ATOM      2  N  BALA A   1      10.300  10.000  10.000  0.30 10.00           N
ATOM      3  N  CALA A   1      10.600  10.000  10.000  0.30 10.00           N
TER
ATOM      4  N  AALA B  10      20.000  10.000  10.000  0.50 10.00           N
ATOM      5  N  BALA B  10      20.300  10.000  10.000  0.50 10.00           N
ATOM      6  N  AALA B  11      23.000  10.000  10.000  0.50 10.00           N
ATOM      7  N   ALA B  12      26.000  10.000  10.000  1.00 10.00           N
ATOM      8  CA AALA B  12      27.000  10.000  10.000  0.50 10.00           C
ATOM      9  CA BALA B  12      27.000  10.300  10.000  0.50 10.00           C
"""

def run(model_text, prefix, *extra):
  fn_in = prefix + ".pdb"
  fn_out = prefix + "_out.pdb"
  fn_desc = prefix + "_out.txt"
  with open(fn_in, "w") as f:
    f.write(model_text)
  for fn in (fn_out, fn_desc):
    if os.path.exists(fn): os.remove(fn)
  run_program(program_class=reduce2.Program, logger=null_out(),
    args=[fn_in, "approach=add", "output.filename=" + fn_out,
      "output.description_file_name=" + fn_desc, "overwrite=True"] + list(extra))
  with open(fn_desc) as f:
    return f.read()

def movers(desc):
  '''{(pass altloc, Mover line)}: Movers added, and methyls staggered'''
  out, alt = set(), None
  for l in desc.split('\n'):
    m = re.search(r"Running Reduce optimization on model index \d+, "
      r"alternate '(.*)'", l)
    if m:
      alt = m.group(1)
      continue
    m = re.match(r"\s*Added (Mover\w+) \d+ to (chain .*?)( with .*)?$", l)
    if m:
      out.add((alt, m.group(1) + ' ' + m.group(2)))
    m = re.match(r"\s*Used (Mover\w+) to stagger (chain .*)$", l)
    if m:
      out.add((alt, 'stagger ' + m.group(2)))
  return out

def swap_ab(model_text):
  return '\n'.join(l[:16] + {'A': 'B', 'B': 'A'}.get(l[16], l[16]) + l[17:]
    if l.startswith(('ATOM', 'HETATM')) else l for l in model_text.split('\n'))

def test_c2():
  '''blank SG, Zn only in A: HG/A is not placed (thiolate), HG/B is; in the B
  pass SG has CB and HG/B only, and HG/B gets a rotator. Labels swapped: HG/A
  in the A pass'''
  m = movers(run(m_c2, 'tst_reduce2_optimizer_conformers_c2'))
  assert ('B', 'MoverSingleHydrogenRotator chain A BCYS 48 HG') in m, m
  m = movers(run(swap_ab(m_c2), 'tst_reduce2_optimizer_conformers_c2_swap'))
  assert ('A', 'MoverSingleHydrogenRotator chain A ACYS 48 HG') in m, m

def test_c2_after_ter():
  '''The same with the Zn after TER, a chain of its own with only altloc A:
  the B pass leaves that residue out (only in another conformer), so the
  result is test_c2's: HG/B gets a rotator; labels swapped, HG/A'''
  m = movers(run(m_c2_ter, 'tst_reduce2_optimizer_conformers_c2_ter'))
  assert ('B', 'MoverSingleHydrogenRotator chain A BCYS 48 HG') in m, m
  m = movers(run(swap_ab(m_c2_ter), 'tst_reduce2_optimizer_conformers_c2_ter_swap'))
  assert ('A', 'MoverSingleHydrogenRotator chain A ACYS 48 HG') in m, m

def test_tyr_oh():
  '''1US0 Tyr A 189: backbone in A/B, ring blank; HH/A and HH/B each get a
  rotator in their pass, also with the labels swapped'''
  for model, name in ((m_tyr, 'tyr'), (swap_ab(m_tyr), 'tyr_swap')):
    m = movers(run(model, 'tst_reduce2_optimizer_conformers_' + name))
    for alt in ('A', 'B'):
      assert (alt, 'MoverSingleHydrogenRotator chain A %sTYR 189 HH' % alt) \
        in m, (name, alt, m)

def test_lys_nh3():
  '''1US0 Lys A 85: NZ blank, its three H in A and in B; NH3 rotator in both
  passes'''
  m = movers(run(m_lys, 'tst_reduce2_optimizer_conformers_lys'))
  for alt in ('A', 'B'):
    assert (alt, 'MoverNH3Rotator chain A LYS 85') in m, (alt, m)

def test_ile_methyl():
  '''1US0 Ile A 14: CD1 and CG2 blank, their H in A and in B; staggered in
  both passes'''
  m = movers(run(m_ile, 'tst_reduce2_optimizer_conformers_ile'))
  for alt in ('A', 'B'):
    for methyl in ('CD1', 'CG2'):
      assert any(x[0] == alt and x[1].startswith('stagger chain A') and
                 x[1].endswith('ILE 14 ' + methyl) for x in m), (alt, methyl, m)

def test_thr_split_branch():
  '''3CP5 Thr A 46: CB blank, OG1 and CG2 in A and B; in each pass HG1 gets a
  rotator and CG2 a stagger Mover, also with the labels swapped. Riding already
  staggers CG2, so the HG2x torsions are only a secondary check.'''
  for model, name in ((m_thr, 'thr'), (swap_ab(m_thr), 'thr_swap')):
    prefix = 'tst_reduce2_optimizer_conformers_' + name
    desc = run(model, prefix)
    refused = [l for l in desc.split('\n')
      if l.startswith('Could not add') and 'THR 46' in l]
    assert not refused, (name, refused)
    m = movers(desc)
    h = iotbx.pdb.input(prefix + '_out.pdb').construct_hierarchy()
    for alt in ('A', 'B'):
      assert (alt, 'MoverSingleHydrogenRotator chain A %sTHR 46 HG1' % alt) \
        in m, (name, alt, m)
      assert any(x[0] == alt and x[1].startswith('stagger chain A') and
                 x[1].endswith('THR 46 CG2') for x in m), (name, alt, m)
      xyz = {}
      for rg in h.residue_groups():
        if rg.resseq_as_int() != 46: continue
        for ag in rg.atom_groups():
          if ag.altloc in ('', alt):
            for a in ag.atoms():
              xyz[a.name.strip()] = a.xyz
      for hn in ('HG21', 'HG22', 'HG23'):
        t = scitbx.math.dihedral_angle(
          sites=[xyz['CA'], xyz['CB'], xyz['CG2'], xyz[hn]], deg=True)
        assert abs((t % 120.0) - 60.0) < 5.0, (name, alt, hn, t)

def test_glycan_copies():
  '''Two superposed copies of a sugar, each in its own chain with one altloc:
  each pass holds the Movers of its own copy only, and no Clique spans the
  copies (before: both copies in both passes, one Clique of 10 per pass)'''
  cliques = []
  orig = Optimizers.OptimizerC
  class Recorder(object):
    def __init__(self, *a, **k):
      self._optC = orig(*a, **k)
    def __getattr__(self, name):
      return getattr(self._optC, name)
    def OptimizeCliqueCoarse(self, movers, edges, seconds):
      cliques.append(set(a.parent().parent().parent().id
        for m in movers for a in m.CoarsePositions().atoms))
      return self._optC.OptimizeCliqueCoarse(movers, edges, seconds)
  Optimizers.OptimizerC = Recorder
  try:
    m = movers(run(m_bma, 'tst_reduce2_optimizer_conformers_bma'))
  finally:
    Optimizers.OptimizerC = orig
  hos = ('HO1', 'HO2', 'HO3', 'HO4', 'HO6')
  assert m == set([('A', 'MoverSingleHydrogenRotator chain A ABMA 414 ' + h) for h in hos] +
    [('B', 'MoverSingleHydrogenRotator chain C BBMA 1 ' + h) for h in hos]), m
  assert cliques and all(len(c) == 1 for c in cliques), cliques

def test_three_altlocs():
  '''Pass C, chain B without C: Ala 10 (A and B) and Ala 12 (blank and A, B)
  give their first conformer (A), Ala 11 (A only) is left out; pass A keeps
  it.'''
  model = iotbx.pdb.input(source_info=None, lines=m_three).construct_hierarchy().models()[0]
  def ids(alt):
    return sorted((a.parent().parent().parent().id, a.parent().parent().resseq_as_int(),
      a.name.strip(), a.parent().altloc) for a in Optimizers.GetAtomsForConformer(model, alt))
  assert ids('C') == [('A', 1, 'N', 'C'), ('B', 10, 'N', 'A'), ('B', 12, 'CA', 'A'),
    ('B', 12, 'N', '')], ids('C')
  assert ('B', 11, 'N', 'A') in ids('A'), ids('A')

def test_full_occupancy():
  '''Ala 11 (A only) at occupancy 0.99 is present in every state: pass C keeps
  it; at 0.98 it is left out.'''
  for occ, kept in (('0.99', True), ('0.98', False)):
    text = m_three.replace('23.000  10.000  10.000  0.50', '23.000  10.000  10.000  ' + occ)
    model = iotbx.pdb.input(source_info=None, lines=text).construct_hierarchy().models()[0]
    got = [a for a in Optimizers.GetAtomsForConformer(model, 'C')
      if a.parent().parent().resseq_as_int() == 11]
    assert len(got) == (1 if kept else 0), (occ, [a.id_str() for a in got])

def test_alt_id_single_pass():
  '''alt_id=B (one pass): the Zn after TER, only in A, is left out of it as in
  the full run; HG/B gets a rotator'''
  m = movers(run(m_c2_ter, 'tst_reduce2_optimizer_conformers_c2_ter_b', 'alt_id=B'))
  assert set(x[0] for x in m) == set(['B']), m
  assert ('B', 'MoverSingleHydrogenRotator chain A BCYS 48 HG') in m, m

if __name__ == "__main__":
  test_c2()
  test_c2_after_ter()
  test_tyr_oh()
  test_lys_nh3()
  test_ile_methyl()
  test_thr_split_branch()
  test_glycan_copies()
  test_three_altlocs()
  test_full_occupancy()
  test_alt_id_single_pass()
  print("OK")
