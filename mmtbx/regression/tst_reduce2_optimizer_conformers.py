"""reduce2 Optimizer: each altloc pass uses the bonded neighbours of its own
conformer. Placement copies the H of a residue with altlocs into each altloc
group, so a blank parent carries H in A and in B; with neighbour lists built
from all atoms the parent counted both copies (and atoms that exist only in
another conformer) and got no rotator, NH3 rotator or methyl stagger.
Models: the Zn-Cys model of tst_add_hydrogen_1 with O split A/B and the Zn in
altloc A (A4 B0 case c2); 1US0 Tyr A 189, Lys A 85 and Ile A 14 (backbone
atoms split A/B, side chains blank).
"""
from __future__ import absolute_import, division, print_function
import os, re
from libtbx.utils import null_out
from iotbx.cli_parser import run_program
from mmtbx.programs import reduce2

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

def run(model_text, prefix):
  fn_in = prefix + ".pdb"
  fn_out = prefix + "_out.pdb"
  fn_desc = prefix + "_out.txt"
  with open(fn_in, "w") as f:
    f.write(model_text)
  for fn in (fn_out, fn_desc):
    if os.path.exists(fn): os.remove(fn)
  run_program(program_class=reduce2.Program, logger=null_out(),
    args=[fn_in, "approach=add", "output.filename=" + fn_out,
      "output.description_file_name=" + fn_desc, "overwrite=True"])
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
  the B pass takes that chain's first conformer, so ZN/A is in the B pass and
  SG has three neighbours there; no rotator for HG/B. Pins today's result,
  see plan item B8 (passes before the per-pass lists too).'''
  m = movers(run(m_c2_ter, 'tst_reduce2_optimizer_conformers_c2_ter'))
  assert not [x for x in m if x[1].endswith('CYS 48 HG')], m

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

if __name__ == "__main__":
  test_c2()
  test_c2_after_ter()
  test_tyr_oh()
  test_lys_nh3()
  test_ile_methyl()
  print("OK")
