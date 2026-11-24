"""reduce2: a peptide bond that an altloc pass holds across two letters (A70).
A pass gives every residue one conformer (Optimizers.ResidueAtomsForConformer()),
so it can hold C of residue i in one letter next to N of residue i+1 in another
(Ala A/B, then Ser A/B/C: pass C holds C of Ala A). Restraints are made per
conformer: no proxy bonds that pair, riding treats the amide H as rotatable and
the Optimizer gave it a single-H rotator (master wrote it 0.52 A from C of Ala A).
Optimizers.CrossLetterPeptideBonds() decides, for placement, the pass's neighbor
lists and so its Movers: same names bonded in another letter combination, within
1.25 times that proxy's ideal distance. The control is the pass's heavy atoms as
one unlettered conformer.
"""
import os, re, math
from libtbx.utils import null_out
from iotbx.cli_parser import run_program
import iotbx.pdb
import mmtbx.model
from mmtbx.programs import reduce2
from mmtbx.reduce import Optimizers
from mmtbx.hydrogens import reduce_hydrogen

def run(model_text, prefix, *extra):
  """reduce2 (approach=add unless given); returns (description, output hierarchy)"""
  fn_in, fn_out, fn_desc = prefix + ".pdb", prefix + "_out.pdb", prefix + "_out.txt"
  with open(fn_in, "w") as f:
    f.write(model_text)
  for fn in (fn_out, fn_desc):
    if os.path.exists(fn): os.remove(fn)
  approach = [] if any(e.startswith('approach=') for e in extra) else ["approach=add"]
  run_program(program_class=reduce2.Program, logger=null_out(),
    args=[fn_in, "output.filename=" + fn_out,
      "output.description_file_name=" + fn_desc, "overwrite=True"] + approach + list(extra))
  with open(fn_desc) as f:
    desc = f.read()
  with open(fn_out) as f:
    return desc, f.read()

def sites(text):
  """{(resseq, name, altloc): xyz}"""
  h = iotbx.pdb.input(source_info=None, lines=text).construct_hierarchy()
  return dict(((int(a.parent().parent().resseq), a.name.strip(),
    a.parent().altloc.strip()), a.xyz) for a in h.atoms())

def passes(desc):
  """{pass altloc: its part of the description}"""
  out, alt = {}, None
  for l in desc.split('\n'):
    m = re.search(r"Running Reduce optimization on model index \d+, alternate '(.*)'", l)
    if m:
      alt = m.group(1)
      out[alt] = ''
    elif alt is not None:
      out[alt] += l + '\n'
  return out

def control(text, held):
  """the heavy atoms of residue r in letter held[r] (or blank), unlettered"""
  out = []
  for l in text.splitlines():
    if l.startswith(('ATOM', 'HETATM')):
      if l[16] not in (' ', held[int(l[22:26])]):
        continue
      l = l[:16] + ' ' + l[17:54] + '  1.00' + l[60:]
    out.append(l)
  return '\n'.join(out) + '\n'

def stretched(text, distance):
  """Ser C moved along C(Ala A)->N(Ser C) to that distance"""
  s = sites(text)
  c, n = s[1, 'C', 'A'], s[2, 'N', 'C']
  d = [(y - x) * (distance / math.dist(c, n) - 1) for x, y in zip(c, n)]
  out = []
  for l in text.splitlines():
    if l.startswith('ATOM') and l[16] == 'C':
      xyz = [float(l[k:k+8]) + e for k, e in zip((30, 38, 46), d)]
      l = l[:30] + ''.join('%8.3f' % v for v in xyz) + l[54:]
    out.append(l)
  return '\n'.join(out) + '\n'

def model_of(text):
  """model with restraints, as reduce2 placement interprets it"""
  m = mmtbx.model.manager(model_input=iotbx.pdb.input(source_info=None, lines=text),
    log=null_out())
  p = reduce_hydrogen.get_reduce_pdb_interpretation_params(False)
  m.process(make_restraints=True, pdb_interpretation_params=p)
  return m

def bonds(text):
  """{pass: [(C id, N id, accepted)]} of CrossLetterPeptideBonds()"""
  m = model_of(text)
  table = m.get_restraints_manager().geometry.bond_params_table
  out = {}
  for x in ('A', 'B', 'C'):
    out[x] = [(c.id_str(), n.id_str(), ok) for c, n, cr, nr, d, ok in
      Optimizers.CrossLetterPeptideBonds(m.get_hierarchy().models()[0], x, table)]
  return out

PEP = """CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00 10.00           N
ATOM      2  CA AALA A   1      11.458  10.000  10.000  0.50 10.00           C
ATOM      3  CA BALA A   1      11.458  10.000  10.200  0.50 10.00           C
ATOM      4  C  AALA A   1      12.009  11.420  10.000  0.50 10.00           C
ATOM      5  C  BALA A   1      12.009  11.420  10.200  0.50 10.00           C
ATOM      6  O  AALA A   1      11.249  12.390  10.000  0.50 10.00           O
ATOM      7  O  BALA A   1      11.249  12.390  10.200  0.50 10.00           O
ATOM      8  CB AALA A   1      11.988   9.236  11.200  0.50 10.00           C
ATOM      9  CB BALA A   1      11.988   9.236  11.400  0.50 10.00           C
ATOM     10  N  ASER A   2      13.330  11.540  10.000  0.33 10.00           N
ATOM     11  N  BSER A   2      13.330  11.540  10.200  0.33 10.00           N
ATOM     12  N  CSER A   2      13.330  11.540   9.750  0.33 10.00           N
ATOM     13  CA ASER A   2      13.980  12.840  10.000  0.33 10.00           C
ATOM     14  CA BSER A   2      13.980  12.840  10.200  0.33 10.00           C
ATOM     15  CA CSER A   2      13.980  12.840   9.750  0.33 10.00           C
ATOM     16  C  ASER A   2      15.490  12.700  10.000  0.33 10.00           C
ATOM     17  C  BSER A   2      15.490  12.700  10.200  0.33 10.00           C
ATOM     18  C  CSER A   2      15.490  12.700   9.750  0.33 10.00           C
ATOM     19  O  ASER A   2      16.050  11.600  10.000  0.33 10.00           O
ATOM     20  O  BSER A   2      16.050  11.600  10.200  0.33 10.00           O
ATOM     21  O  CSER A   2      16.050  11.600   9.750  0.33 10.00           O
ATOM     22  CB ASER A   2      13.530  13.640  11.220  0.33 10.00           C
ATOM     23  CB BSER A   2      13.530  13.640  11.420  0.33 10.00           C
ATOM     24  CB CSER A   2      13.530  13.640  10.970  0.33 10.00           C
ATOM     25  OG ASER A   2      14.020  14.970  11.180  0.33 10.00           O
ATOM     26  OG BSER A   2      14.020  14.970  11.380  0.33 10.00           O
ATOM     27  OG CSER A   2      14.020  14.970  10.930  0.33 10.00           O
ATOM     28  N   GLY A   3      16.150  13.850  10.000  1.00 10.00           N
ATOM     29  CA  GLY A   3      17.600  13.900  10.000  1.00 10.00           C
ATOM     30  C   GLY A   3      18.100  15.330  10.000  1.00 10.00           C
ATOM     31  O   GLY A   3      17.350  16.300  10.000  1.00 10.00           O
END
"""

def test_peptide():
  """H of Ser C as in the control, no rotator in pass C, clear of C(Ala A)"""
  desc, out = run(PEP, 'cl_pep')
  ctl = sites(run(control(PEP, {1: 'A', 2: 'C', 3: ' '}), 'cl_pep_ctl')[1])
  s = sites(out)
  assert math.dist(s[2, 'H', 'C'], ctl[2, 'H', '']) < 0.005, (s[2, 'H', 'C'], ctl[2, 'H', ''])
  assert math.dist(s[2, 'H', 'C'], s[1, 'C', 'A']) > 1.8
  part = passes(desc)['C']
  assert 'Peptide bond across letters: " C  AALA A   1 " - " N  CSER A   2 " 1.35 A, bonded' \
    in part, part
  assert 'CSER 2 H ' not in part, part
  assert 'Amide H set' not in desc
  # its own passes: unchanged, C(Ala A/B)-N(Ser A/B) have proxies
  assert 'across letters' not in passes(desc)['A'] + passes(desc)['B']

def test_place():
  """approach=place puts it there too (the call in placement)"""
  desc, out = run(PEP, 'cl_pep_place', 'approach=place')
  ctl = sites(run(control(PEP, {1: 'A', 2: 'C', 3: ' '}), 'cl_pep_ctl_place',
    'approach=place')[1])
  assert math.dist(sites(out)[2, 'H', 'C'], ctl[2, 'H', '']) < 0.005

def reverse(text):
  """Ala A/B/C (C a copy of A), Ser A/B: pass C holds C of Ala C, N of Ser A"""
  out = []
  for l in text.splitlines():
    if l.startswith('ATOM') and l[16] == 'C':
      continue
    out.append(l)
    if l.startswith('ATOM') and l[16] == 'A' and int(l[22:26]) == 1:
      z = float(l[46:54]) - 0.25
      out.append(l[:16] + 'C' + l[17:46] + '%8.3f' % z + l[54:])
  return '\n'.join(out) + '\n'

def test_reverse():
  """the pass's C is its own letter, N borrowed: bonded, H(A) left to pass A"""
  text = reverse(PEP)
  b = bonds(text)
  assert b['C'] == [('pdb=" C  CALA A   1 "', 'pdb=" N  ASER A   2 "', True)], b
  desc, out = run(text, 'cl_rev')
  assert 'Peptide bond across letters: " C  CALA A   1 " - " N  ASER A   2 "' \
    in passes(desc)['C']
  ctl = sites(run(control(text, {1: 'A', 2: 'A', 3: ' '}), 'cl_rev_ctl')[1])
  assert math.dist(sites(out)[2, 'H', 'A'], ctl[2, 'H', '']) < 0.005

def test_riding_entry():
  """the H gets a planar riding entry, so idealizing again keeps it; without
  the H of the bonded combination, the exterior bisector"""
  desc, out = run(PEP, 'cl_entry')
  for drop in ((), ('A', 'B')):
    lines = [l for l in out.splitlines() if not (l.startswith('ATOM') and
      l[12:16] == ' H  ' and l[16] in drop and l[17:26] == 'SER A   2')]
    lines = [l[:30] + ''.join('%8.3f' % v for v in MASTER_H) + l[54:]
      if l.startswith('ATOM') and l[12:17] == ' H  C' and l[17:26] == 'SER A   2'
      else l for l in lines]
    m = model_of('\n'.join(lines) + '\n')
    moved = reduce_hydrogen.place_amide_h_across_letters(m)
    assert len(moved) == 1 and moved[0][1] > 1.5, moved
    h = moved[0][0]
    p = m.get_riding_h_manager().h_parameterization[h.i_seq]
    assert p.htype == 'flat_2neigbs'
    s = sites(m.model_as_pdb())
    x = s[2, 'H', 'C']
    if not drop:
      assert math.dist(x, sites(out)[2, 'H', 'C']) < 0.005
    else:
      c, n, ca = s[1, 'C', 'A'], s[2, 'N', 'C'], s[2, 'CA', 'C']
      def angle(a, b):
        u = [p - q for p, q in zip(a, n)]; v = [p - q for p, q in zip(b, n)]
        return math.degrees(math.acos(sum(i * j for i, j in zip(u, v)) /
          (math.dist(a, n) * math.dist(b, n))))
      # 0.05 deg: coordinates written to 0.001 A
      assert abs(angle(c, x) - angle(ca, x)) < 0.05, (angle(c, x), angle(ca, x))
      assert abs(angle(c, x) + angle(ca, x) + angle(c, ca) - 360) < 0.05
    m.idealize_h_riding()
    assert math.dist(sites(m.model_as_pdb())[2, 'H', 'C'], x) < 0.001

def test_optimizer_alone():
  """the Optimizer on a model whose H did not go through placement (master's
  output, rotatable H): the pass bonds the pair, so no rotator"""
  from mmtbx.hydrogens import default_probe_phil
  desc, out = run(PEP, 'cl_alone')
  lines = [l[:30] + ''.join('%8.3f' % v for v in MASTER_H) + l[54:]
    if l.startswith('ATOM') and l[12:17] == ' H  C' and l[17:26] == 'SER A   2'
    else l for l in out.splitlines()]
  m = model_of('\n'.join(lines) + '\n')
  opt = Optimizers.Optimizer(default_probe_phil(), True, m, fillAtomDump=False)
  part = passes(opt.getInfo())['C']
  assert '" N  CSER A   2 " 1.35 A, bonded' in part, part
  assert 'CSER 2 H ' not in part, part

def test_stretched():
  """C-N 1.90 A: not bonded by placement nor by the Optimizer"""
  text = stretched(PEP, 1.90)
  b = bonds(text)
  assert b['C'] == [('pdb=" C  AALA A   1 "', 'pdb=" N  CSER A   2 "', False)], b
  assert b['A'] == [] and b['B'] == []
  # placement leaves it where master puts it (riding, rotatable)
  desc, out = run(text, 'cl_far_place', 'approach=place')
  assert math.dist(sites(out)[2, 'H', 'C'], (14.115, 11.075, 9.004)) < 0.002
  desc, out = run(text, 'cl_far')
  part = passes(desc)['C']
  assert 'not bonded (over 1.25 x 1.33 A)' in part, part
  assert 'Added MoverSingleHydrogenRotator 2 to chain A CSER 2 H ' in part, part

BLANK = """CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00 10.00           N
ATOM      2  CA  ALA A   1      11.458  10.000  10.000  1.00 10.00           C
ATOM      3  C   ALA A   1      12.009  11.420  10.000  1.00 10.00           C
ATOM      4  O   ALA A   1      11.249  12.390  10.000  1.00 10.00           O
ATOM      5  CB AALA A   1      11.988   9.236  11.200  0.50 10.00           C
ATOM      6  CB BALA A   1      11.988   9.236  11.400  0.50 10.00           C
ATOM      7  N  ASER A   2      13.330  11.540  10.000  0.33 10.00           N
ATOM      8  N  BSER A   2      13.330  11.540  10.200  0.33 10.00           N
ATOM      9  N  CSER A   2      13.330  11.540   9.750  0.33 10.00           N
ATOM     10  CA ASER A   2      13.980  12.840  10.000  0.33 10.00           C
ATOM     11  CA BSER A   2      13.980  12.840  10.200  0.33 10.00           C
ATOM     12  CA CSER A   2      13.980  12.840   9.750  0.33 10.00           C
ATOM     13  C  ASER A   2      15.490  12.700  10.000  0.33 10.00           C
ATOM     14  C  BSER A   2      15.490  12.700  10.200  0.33 10.00           C
ATOM     15  C  CSER A   2      15.490  12.700   9.750  0.33 10.00           C
ATOM     16  O  ASER A   2      16.050  11.600  10.000  0.33 10.00           O
ATOM     17  O  BSER A   2      16.050  11.600  10.200  0.33 10.00           O
ATOM     18  O  CSER A   2      16.050  11.600   9.750  0.33 10.00           O
ATOM     19  CB ASER A   2      13.530  13.640  11.220  0.33 10.00           C
ATOM     20  CB BSER A   2      13.530  13.640  11.420  0.33 10.00           C
ATOM     21  CB CSER A   2      13.530  13.640  10.970  0.33 10.00           C
ATOM     22  OG ASER A   2      14.020  14.970  11.180  0.33 10.00           O
ATOM     23  OG BSER A   2      14.020  14.970  11.380  0.33 10.00           O
ATOM     24  OG CSER A   2      14.020  14.970  10.930  0.33 10.00           O
ATOM     25  N   GLY A   3      16.150  13.850  10.000  1.00 10.00           N
ATOM     26  CA  GLY A   3      17.600  13.900  10.000  1.00 10.00           C
ATOM     27  C   GLY A   3      18.100  15.330  10.000  1.00 10.00           C
ATOM     28  O   GLY A   3      17.350  16.300  10.000  1.00 10.00           O
END
"""

NA = """CRYST1   80.000   80.000   80.000  90.00  90.00  90.00 P 1
ATOM      1  P     A A   2      10.342  -5.036   5.198  1.00  0.00           P
ATOM      2  OP1   A A   2      11.732  -5.326   5.616  1.00  0.00           O
ATOM      3  OP2   A A   2       9.245  -5.144   6.185  1.00  0.00           O
ATOM      4  O5'   A A   2      10.301  -3.545   4.582  1.00  0.00           O
ATOM      5  C5'   A A   2      11.387  -2.992   3.859  1.00  0.00           C
ATOM      6  C4'   A A   2      11.046  -1.560   3.437  1.00  0.00           C
ATOM      7  O4'   A A   2       9.807  -1.514   2.753  1.00  0.00           O
ATOM      8  C3'A  A A   2      10.900  -0.617   4.629  0.50  0.00           C
ATOM      9  C3'B  A A   2      10.900  -0.617   4.779  0.50  0.00           C
ATOM     10  O3'A  A A   2      12.167  -0.168   5.083  0.50  0.00           O
ATOM     11  O3'B  A A   2      12.167  -0.168   5.233  0.50  0.00           O
ATOM     12  C2'   A A   2      10.048   0.493   4.017  1.00  0.00           C
ATOM     13  O2'   A A   2      10.862   1.442   3.353  1.00  0.00           O
ATOM     14  C1'   A A   2       9.217  -0.242   2.962  1.00  0.00           C
ATOM     15  N9    A A   2       7.794  -0.369   3.344  1.00  0.00           N
ATOM     16  C8    A A   2       7.080  -1.486   3.707  1.00  0.00           C
ATOM     17  N7    A A   2       5.798  -1.279   3.832  1.00  0.00           N
ATOM     18  C5    A A   2       5.646   0.067   3.520  1.00  0.00           C
ATOM     19  C6    A A   2       4.521   0.913   3.401  1.00  0.00           C
ATOM     20  N6    A A   2       3.267   0.490   3.594  1.00  0.00           N
ATOM     21  N1    A A   2       4.723   2.204   3.071  1.00  0.00           N
ATOM     22  C2    A A   2       5.965   2.631   2.857  1.00  0.00           C
ATOM     23  N3    A A   2       7.096   1.935   2.908  1.00  0.00           N
ATOM     24  C4    A A   2       6.863   0.637   3.246  1.00  0.00           C
ATOM     25  P  A  A A   3      12.341   0.647   6.465  0.33  0.00           P
ATOM     26  P  B  A A   3      12.341   0.647   6.615  0.33  0.00           P
ATOM     27  P  C  A A   3      12.341   0.647   6.265  0.33  0.00           P
ATOM     28  OP1A  A A   3      13.769   1.006   6.611  0.33  0.00           O
ATOM     29  OP1B  A A   3      13.769   1.006   6.761  0.33  0.00           O
ATOM     30  OP1C  A A   3      13.769   1.006   6.411  0.33  0.00           O
ATOM     31  OP2A  A A   3      11.662  -0.117   7.534  0.33  0.00           O
ATOM     32  OP2B  A A   3      11.662  -0.117   7.684  0.33  0.00           O
ATOM     33  OP2C  A A   3      11.662  -0.117   7.334  0.33  0.00           O
ATOM     34  O5'A  A A   3      11.514   2.005   6.194  0.33  0.00           O
ATOM     35  O5'B  A A   3      11.514   2.005   6.344  0.33  0.00           O
ATOM     36  O5'C  A A   3      11.514   2.005   5.994  0.33  0.00           O
ATOM     37  C5'A  A A   3      10.977   2.765   7.257  0.33  0.00           C
ATOM     38  C5'B  A A   3      10.977   2.765   7.407  0.33  0.00           C
ATOM     39  C5'C  A A   3      10.977   2.765   7.057  0.33  0.00           C
ATOM     40  C4'A  A A   3      10.049   3.851   6.701  0.33  0.00           C
ATOM     41  C4'B  A A   3      10.049   3.851   6.851  0.33  0.00           C
ATOM     42  C4'C  A A   3      10.049   3.851   6.501  0.33  0.00           C
ATOM     43  O4'A  A A   3       9.026   3.266   5.913  0.33  0.00           O
ATOM     44  O4'B  A A   3       9.026   3.266   6.063  0.33  0.00           O
ATOM     45  O4'C  A A   3       9.026   3.266   5.713  0.33  0.00           O
ATOM     46  C3'A  A A   3       9.295   4.550   7.828  0.33  0.00           C
ATOM     47  C3'B  A A   3       9.295   4.550   7.978  0.33  0.00           C
ATOM     48  C3'C  A A   3       9.295   4.550   7.628  0.33  0.00           C
ATOM     49  O3'A  A A   3      10.040   5.574   8.463  0.33  0.00           O
ATOM     50  O3'B  A A   3      10.040   5.574   8.613  0.33  0.00           O
ATOM     51  O3'C  A A   3      10.040   5.574   8.263  0.33  0.00           O
ATOM     52  C2'A  A A   3       8.080   5.100   7.082  0.33  0.00           C
ATOM     53  C2'B  A A   3       8.080   5.100   7.232  0.33  0.00           C
ATOM     54  C2'C  A A   3       8.080   5.100   6.882  0.33  0.00           C
ATOM     55  O2'A  A A   3       8.371   6.366   6.525  0.33  0.00           O
ATOM     56  O2'B  A A   3       8.371   6.366   6.675  0.33  0.00           O
ATOM     57  O2'C  A A   3       8.371   6.366   6.325  0.33  0.00           O
ATOM     58  C1'A  A A   3       7.886   4.109   5.931  0.33  0.00           C
ATOM     59  C1'B  A A   3       7.886   4.109   6.081  0.33  0.00           C
ATOM     60  C1'C  A A   3       7.886   4.109   5.731  0.33  0.00           C
ATOM     61  N9 A  A A   3       6.637   3.330   6.084  0.33  0.00           N
ATOM     62  N9 B  A A   3       6.637   3.330   6.234  0.33  0.00           N
ATOM     63  N9 C  A A   3       6.637   3.330   5.884  0.33  0.00           N
ATOM     64  C8 A  A A   3       6.473   1.983   6.295  0.33  0.00           C
ATOM     65  C8 B  A A   3       6.473   1.983   6.445  0.33  0.00           C
ATOM     66  C8 C  A A   3       6.473   1.983   6.095  0.33  0.00           C
ATOM     67  N7 A  A A   3       5.227   1.602   6.366  0.33  0.00           N
ATOM     68  N7 B  A A   3       5.227   1.602   6.516  0.33  0.00           N
ATOM     69  N7 C  A A   3       5.227   1.602   6.166  0.33  0.00           N
ATOM     70  C5 A  A A   3       4.509   2.773   6.143  0.33  0.00           C
ATOM     71  C5 B  A A   3       4.509   2.773   6.293  0.33  0.00           C
ATOM     72  C5 C  A A   3       4.509   2.773   5.943  0.33  0.00           C
ATOM     73  C6 A  A A   3       3.130   3.053   6.015  0.33  0.00           C
ATOM     74  C6 B  A A   3       3.130   3.053   6.165  0.33  0.00           C
ATOM     75  C6 C  A A   3       3.130   3.053   5.815  0.33  0.00           C
ATOM     76  N6 A  A A   3       2.183   2.119   6.151  0.33  0.00           N
ATOM     77  N6 B  A A   3       2.183   2.119   6.301  0.33  0.00           N
ATOM     78  N6 C  A A   3       2.183   2.119   5.951  0.33  0.00           N
ATOM     79  N1 A  A A   3       2.754   4.315   5.725  0.33  0.00           N
ATOM     80  N1 B  A A   3       2.754   4.315   5.875  0.33  0.00           N
ATOM     81  N1 C  A A   3       2.754   4.315   5.525  0.33  0.00           N
ATOM     82  C2 A  A A   3       3.691   5.248   5.573  0.33  0.00           C
ATOM     83  C2 B  A A   3       3.691   5.248   5.723  0.33  0.00           C
ATOM     84  C2 C  A A   3       3.691   5.248   5.373  0.33  0.00           C
ATOM     85  N3 A  A A   3       5.012   5.114   5.663  0.33  0.00           N
ATOM     86  N3 B  A A   3       5.012   5.114   5.813  0.33  0.00           N
ATOM     87  N3 C  A A   3       5.012   5.114   5.463  0.33  0.00           N
ATOM     88  C4 A  A A   3       5.362   3.831   5.950  0.33  0.00           C
ATOM     89  C4 B  A A   3       5.362   3.831   6.100  0.33  0.00           C
ATOM     90  C4 C  A A   3       5.362   3.831   5.750  0.33  0.00           C
END
"""

def test_controls():
  """blank backbone between the letters; nucleic acid O3'(A)-P(C): no pair"""
  for text, prefix in ((BLANK, 'cl_blank'), (NA, 'cl_na')):
    assert bonds(text) == {'A': [], 'B': [], 'C': []}
    desc, out = run(text, prefix)
    assert 'across letters' not in desc

PQ7 = """CRYST1   32.869   37.022   39.780 102.89 104.59 102.37 P 1
ATOM    165  N   TRP A  40       6.387  10.256   1.511  1.00  5.49           N
ATOM    166  CA  TRP A  40       7.193   9.080   1.247  1.00  5.05           C
ATOM    167  C   TRP A  40       6.321   7.884   0.831  1.00  4.88           C
ATOM    168  O   TRP A  40       6.480   7.354  -0.277  1.00  5.86           O
ATOM    169  CB  TRP A  40       8.095   8.763   2.458  1.00  5.02           C
ATOM    170  CG  TRP A  40       9.027   7.641   2.206  1.00  5.02           C
ATOM    171  CD1 TRP A  40       8.993   6.404   2.791  1.00  5.31           C
ATOM    172  CD2 TRP A  40      10.125   7.632   1.305  1.00  5.33           C
ATOM    173  NE1 TRP A  40      10.003   5.631   2.317  1.00  6.10           N
ATOM    174  CE2 TRP A  40      10.722   6.354   1.386  1.00  5.95           C
ATOM    175  CE3 TRP A  40      10.668   8.596   0.434  1.00  7.16           C
ATOM    176  CZ2 TRP A  40      11.829   6.034   0.622  1.00  7.80           C
ATOM    177  CZ3 TRP A  40      11.776   8.226  -0.326  1.00  9.41           C
ATOM    178  CH2 TRP A  40      12.325   6.985  -0.228  1.00  9.46           C
ATOM    179  N  ACYS A  41       5.508   7.360   1.752  0.59  3.97           N
ATOM    180  CA ACYS A  41       4.721   6.164   1.517  0.60  3.70           C
ATOM    181  C  ACYS A  41       3.333   6.308   2.146  0.60  3.80           C
ATOM    182  O  ACYS A  41       2.979   7.330   2.757  0.60  3.75           O
ATOM    183  CB ACYS A  41       5.431   4.898   2.010  0.61  3.50           C
ATOM    184  SG ACYS A  41       6.731   4.335   0.827  0.60  3.83           S
ATOM    185  N  BCYS A  41       5.374   7.503   1.660  0.40  4.85           N
ATOM    186  CA BCYS A  41       4.595   6.308   1.311  0.39  5.11           C
ATOM    187  C  BCYS A  41       3.340   6.305   2.209  0.40  4.25           C
ATOM    188  O  BCYS A  41       3.124   7.213   3.049  0.40  3.20           O
ATOM    189  CB BCYS A  41       5.483   5.140   1.450  0.40  7.12           C
ATOM    190  SG BCYS A  41       5.836   4.613   3.064  0.39 10.10           S
ATOM    191  N   GLY A  42       2.514   5.271   1.971  1.00  4.11           N
ATOM    192  CA  GLY A  42       1.322   5.042   2.727  1.00  4.02           C
ATOM    193  C   GLY A  42       1.417   3.729   3.505  1.00  3.74           C
ATOM    194  O   GLY A  42       2.392   2.983   3.382  1.00  4.06           O
ATOM    299  N   CYS A  57       9.445   0.516   1.315  1.00  4.44           N
ATOM    300  CA  CYS A  57       9.313   1.815   1.946  1.00  4.72           C
ATOM    301  C   CYS A  57      10.241   1.988   3.148  1.00  4.55           C
ATOM    302  O   CYS A  57      10.553   3.132   3.509  1.00  5.19           O
ATOM    303  CB  CYS A  57       7.875   2.020   2.439  1.00  5.19           C
ATOM    304  SG  CYS A  57       6.674   2.294   1.132  1.00  6.27           S
END
"""

ZN = """CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
ATOM      1  N   GLY A   1       8.546  11.010  10.000  1.00 10.00           N
ATOM      2  CA  GLY A   1       8.500  12.468  10.000  1.00 10.00           C
ATOM      3  C   GLY A   1       9.900  13.050  10.000  1.00 10.00           C
ATOM      4  O   GLY A   1      10.880  12.320  10.000  1.00 10.00           O
ATOM      5  N  ACYS A   2       9.990  14.380  10.000  0.50 10.00           N
ATOM      6  N  BCYS A   2       9.990  14.380  10.150  0.50 10.00           N
ATOM      7  CA ACYS A   2      11.280  15.040  10.000  0.50 10.00           C
ATOM      8  CA BCYS A   2      11.280  15.040  10.150  0.50 10.00           C
ATOM      9  C  ACYS A   2      11.100  16.550  10.000  0.50 10.00           C
ATOM     10  C  BCYS A   2      11.100  16.550  10.150  0.50 10.00           C
ATOM     11  O  ACYS A   2       9.980  17.050  10.000  0.50 10.00           O
ATOM     12  O  BCYS A   2       9.980  17.050  10.150  0.50 10.00           O
ATOM     13  CB ACYS A   2      12.130  14.600  11.190  0.50 10.00           C
ATOM     14  CB BCYS A   2      12.130  14.600  11.340  0.50 10.00           C
ATOM     15  SG ACYS A   2      13.780  15.330  11.170  0.50 10.00           S
ATOM     16  SG BCYS A   2      11.600  15.200  12.930  0.50 10.00           S
HETATM   17 ZN  A ZN A 101      15.300  14.000  11.900  0.50 10.00          ZN
END
"""

def test_links_not_peptide():
  """1pq7 Cys 41 SG B 3.1 A from SG of Cys 57 (bonded in A); a Zn in A on
  SG A: no pair, H as on master"""
  desc, out = run(PQ7, 'cl_1pq7')
  assert 'across letters' not in desc
  s = sites(out)
  assert (41, 'HG', 'A') not in s and (57, 'HG', '') not in s
  assert math.dist(s[41, 'HG', 'B'], (6.022, 3.427, 3.056)) < 0.002
  assert math.dist(s[57, 'HG', 'B'], (7.054, 3.193, 0.434)) < 0.002
  desc, out = run(ZN, 'cl_zn')
  assert 'across letters' not in desc
  s = sites(out)
  assert (2, 'HG', 'A') not in s
  assert math.dist(s[2, 'HG', 'B'], (10.815, 16.094, 12.769)) < 0.002

def test_neutron():
  """neutron_option=all_d: N-D 1.02 A, D as in the control"""
  desc, out = run(PEP, 'cl_pep_d', 'neutron_option=all_d')
  ctl = sites(run(control(PEP, {1: 'A', 2: 'C', 3: ' '}), 'cl_pep_d_ctl',
    'neutron_option=all_d')[1])
  s = sites(out)
  assert abs(math.dist(s[2, 'D', 'C'], s[2, 'N', 'C']) - 1.02) < 0.001
  assert math.dist(s[2, 'D', 'C'], ctl[2, 'D', '']) < 0.005

# master's H of Ser C (single-H rotator), 0.52 A from C of Ala A
MASTER_H = (12.475, 11.576, 9.837)

def test_optimize():
  """a right input stays; master's clashing input is repaired; a rerun moves
  nothing"""
  desc, out = run(PEP, 'cl_opt_in')
  right = sites(out)
  desc, again = run(out, 'cl_opt_right', 'approach=optimize')
  again = sites(again)
  # the amide H; the methyl and Ser HG B are re-optimized on master too
  for x in 'ABC':
    assert math.dist(right[2, 'H', x], again[2, 'H', x]) < 0.005, x
  assert 'Amide H set' not in desc
  lines = []
  for l in out.splitlines():
    if l.startswith('ATOM') and l[12:17] == ' H  C' and l[17:26] == 'SER A   2':
      l = l[:30] + ''.join('%8.3f' % v for v in MASTER_H) + l[54:]
    lines.append(l)
  desc, fixed = run('\n'.join(lines) + '\n', 'cl_opt_master', 'approach=optimize')
  assert 'Amide H set on peptide bonds across altloc letters' in desc
  assert '" H  CSER A   2 " moved 1.57 A' in desc, desc
  fixed = sites(fixed)
  assert math.dist(fixed[2, 'H', 'C'], right[2, 'H', 'C']) < 0.005
  desc, rerun = run(open('cl_opt_master_out.pdb').read(), 'cl_opt_rerun',
    'approach=optimize')
  rerun = sites(rerun)
  assert 'Amide H set' not in desc
  # the amide H; Ser HG B has two poses of equal score (-130 and 180 deg) and
  # may still move, on master too
  for x in 'ABC':
    assert math.dist(fixed[2, 'H', x], rerun[2, 'H', x]) < 0.005, x

def omit_amide_h(text, letters):
  """text without the amide H of Ser in those letters"""
  return '\n'.join(l for l in text.splitlines() if not (l.startswith('ATOM') and
    l[12:16] == ' H  ' and l[16] in letters and l[17:26] == 'SER A   2')) + '\n'

def test_reference_h_missing():
  """optimize without the amide H of Ser A, of B, of both: the first two take
  the planar entry of the reference that has one, as in the control; without
  any, the exterior bisector"""
  desc, out = run(PEP, 'cl_ref_in')
  ctl = sites(run(control(PEP, {1: 'A', 2: 'C', 3: ' '}), 'cl_ref_ctl')[1])
  for letters in ('A', 'B'):
    desc, opt = run(omit_amide_h(out, letters), 'cl_ref_omit_' + letters,
      'approach=optimize')
    h = sites(opt)[2, 'H', 'C']
    assert math.dist(h, ctl[2, 'H', '']) < 0.005, (letters, h, ctl[2, 'H', ''])
    assert 'Amide H set' not in desc, (letters, desc)
  desc, opt = run(omit_amide_h(out, 'AB'), 'cl_ref_omit_AB', 'approach=optimize')
  s = sites(opt)
  h, n = s[2, 'H', 'C'], s[2, 'N', 'C']
  def angle(x):
    u = [a - b for a, b in zip(h, n)]
    v = [a - b for a, b in zip(x, n)]
    return math.degrees(math.acos(sum(a * b for a, b in zip(u, v)) /
      math.dist(h, n) / math.dist(x, n)))
  a_c, a_ca = angle(s[1, 'C', 'A']), angle(s[2, 'CA', 'C'])
  assert abs(a_c - a_ca) < 0.05, (a_c, a_ca)
  # in the plane of C(Ala A), N and CA
  c, ca = s[1, 'C', 'A'], s[2, 'CA', 'C']
  u = [a - b for a, b in zip(c, n)]
  v = [a - b for a, b in zip(ca, n)]
  c_n_ca = math.degrees(math.acos(sum(a * b for a, b in zip(u, v)) /
    math.dist(c, n) / math.dist(ca, n)))
  assert abs(a_c + a_ca + c_n_ca - 360) < 0.1, (a_c, a_ca, c_n_ca)

# DNE G 2 - HT7 G 3 of tst_add_hydrogen_23 (9eor), split into DNE A/B and
# HT7 A/B/C: pass C holds C of DNE A next to N of HT7 C. HT7 is a beta-amino
# acid (N bonded to CB, two H on N in conformer C): not an ordinary amide
HT7 = """CRYST1   80.000   80.000   80.000  90.00  90.00  90.00 P 1
HETATM14142  N  ADNE G   2     -19.536 -15.401  62.306  0.50 19.63           N
HETATM14143  CA ADNE G   2     -20.187 -14.689  63.387  0.50 20.76           C
HETATM14144  C  ADNE G   2     -21.461 -15.483  63.644  0.50 19.22           C
HETATM14145  O  ADNE G   2     -21.379 -16.625  64.070  0.50 20.98           O
HETATM14146  CB ADNE G   2     -19.292 -14.498  64.611  0.50 21.84           C
HETATM14147  CG ADNE G   2     -17.813 -14.327  64.273  0.50 23.05           C
HETATM14148  CD ADNE G   2     -17.239 -12.964  64.370  0.50 23.31           C
HETATM14149  CE ADNE G   2     -16.122 -12.745  63.372  0.50 25.60           C
HETATM14142  N  BDNE G   2     -19.536 -15.401  62.506  0.50 19.63           N
HETATM14143  CA BDNE G   2     -20.187 -14.689  63.587  0.50 20.76           C
HETATM14144  C  BDNE G   2     -21.461 -15.483  63.844  0.50 19.22           C
HETATM14145  O  BDNE G   2     -21.379 -16.625  64.270  0.50 20.98           O
HETATM14146  CB BDNE G   2     -19.292 -14.498  64.811  0.50 21.84           C
HETATM14147  CG BDNE G   2     -17.813 -14.327  64.473  0.50 23.05           C
HETATM14148  CD BDNE G   2     -17.239 -12.964  64.570  0.50 23.31           C
HETATM14149  CE BDNE G   2     -16.122 -12.745  63.572  0.50 25.60           C
HETATM14150  N  AHT7 G   3     -22.623 -14.925  63.307  0.33 17.64           N
HETATM14151  CA AHT7 G   3     -24.843 -15.694  62.547  0.33 21.80           C
HETATM14152  C  AHT7 G   3     -24.240 -16.569  61.471  0.33 24.00           C
HETATM14153  O  AHT7 G   3     -23.164 -17.136  61.638  0.33 27.40           O
HETATM14154  CB AHT7 G   3     -23.893 -15.480  63.744  0.33 18.87           C
HETATM14155  CG AHT7 G   3     -24.479 -14.557  64.804  0.33 16.66           C
HETATM14156  CD AHT7 G   3     -23.598 -14.420  66.002  0.33 18.30           C
HETATM14157  CE1AHT7 G   3     -23.619 -15.222  67.111  0.33 19.71           C
HETATM14158  CE2AHT7 G   3     -22.550 -13.456  66.245  0.33 18.01           C
HETATM14159  NZ1AHT7 G   3     -22.684 -14.812  68.030  0.33 19.29           N
HETATM14160  CZ2AHT7 G   3     -22.029 -13.716  67.535  0.33 18.92           C
HETATM14161  CZ3AHT7 G   3     -22.008 -12.396  65.526  0.33 18.29           C
HETATM14162  CH2AHT7 G   3     -20.986 -12.984  68.092  0.33 18.64           C
HETATM14163  CH3AHT7 G   3     -20.988 -11.657  66.089  0.33 19.09           C
HETATM14164  CT2AHT7 G   3     -20.469 -11.960  67.343  0.33 18.11           C
HETATM14150  N  BHT7 G   3     -22.623 -14.925  63.507  0.33 17.64           N
HETATM14151  CA BHT7 G   3     -24.843 -15.694  62.747  0.33 21.80           C
HETATM14152  C  BHT7 G   3     -24.240 -16.569  61.671  0.33 24.00           C
HETATM14153  O  BHT7 G   3     -23.164 -17.136  61.838  0.33 27.40           O
HETATM14154  CB BHT7 G   3     -23.893 -15.480  63.944  0.33 18.87           C
HETATM14155  CG BHT7 G   3     -24.479 -14.557  65.004  0.33 16.66           C
HETATM14156  CD BHT7 G   3     -23.598 -14.420  66.202  0.33 18.30           C
HETATM14157  CE1BHT7 G   3     -23.619 -15.222  67.311  0.33 19.71           C
HETATM14158  CE2BHT7 G   3     -22.550 -13.456  66.445  0.33 18.01           C
HETATM14159  NZ1BHT7 G   3     -22.684 -14.812  68.230  0.33 19.29           N
HETATM14160  CZ2BHT7 G   3     -22.029 -13.716  67.735  0.33 18.92           C
HETATM14161  CZ3BHT7 G   3     -22.008 -12.396  65.726  0.33 18.29           C
HETATM14162  CH2BHT7 G   3     -20.986 -12.984  68.292  0.33 18.64           C
HETATM14163  CH3BHT7 G   3     -20.988 -11.657  66.289  0.33 19.09           C
HETATM14164  CT2BHT7 G   3     -20.469 -11.960  67.543  0.33 18.11           C
HETATM14150  N  CHT7 G   3     -22.623 -14.925  63.207  0.33 17.64           N
HETATM14151  CA CHT7 G   3     -24.843 -15.694  62.447  0.33 21.80           C
HETATM14152  C  CHT7 G   3     -24.240 -16.569  61.371  0.33 24.00           C
HETATM14153  O  CHT7 G   3     -23.164 -17.136  61.538  0.33 27.40           O
HETATM14154  CB CHT7 G   3     -23.893 -15.480  63.644  0.33 18.87           C
HETATM14155  CG CHT7 G   3     -24.479 -14.557  64.704  0.33 16.66           C
HETATM14156  CD CHT7 G   3     -23.598 -14.420  65.902  0.33 18.30           C
HETATM14157  CE1CHT7 G   3     -23.619 -15.222  67.011  0.33 19.71           C
HETATM14158  CE2CHT7 G   3     -22.550 -13.456  66.145  0.33 18.01           C
HETATM14159  NZ1CHT7 G   3     -22.684 -14.812  67.930  0.33 19.29           N
HETATM14160  CZ2CHT7 G   3     -22.029 -13.716  67.435  0.33 18.92           C
HETATM14161  CZ3CHT7 G   3     -22.008 -12.396  65.426  0.33 18.29           C
HETATM14162  CH2CHT7 G   3     -20.986 -12.984  67.992  0.33 18.64           C
HETATM14163  CH3CHT7 G   3     -20.988 -11.657  65.989  0.33 19.09           C
HETATM14164  CT2CHT7 G   3     -20.469 -11.960  67.243  0.33 18.11           C
END"""

# master's H of HT7 C (reduce2 3.24.1, approach=add)
MASTER_HT7_C = {
  'HA': (-23.778, -16.368, 64.017),
  'HA1': (-24.638, -13.686, 64.308),
  'HA2': (-22.217, -15.498, 62.661),
  'HB1': (-25.670, -16.092, 62.760),
  'HB2': (-25.348, -14.896, 64.969),
  'HE1': (-24.158, -15.958, 67.192),
  'HZ1': (-22.539, -15.187, 68.690),
  'HZ3': (-22.327, -12.186, 64.578),
  'HH2': (-20.653, -13.180, 68.838),
  'HH3': (-20.635, -10.935, 65.521),
  'HC2': (-25.081, -14.834, 62.067),
  'HN1': (-22.083, -14.833, 63.908),
  'HT2': (-19.760, -11.460, 67.578)}

def test_not_ordinary_amide():
  """HT7: the bond is decided, but conformer C is written as on master"""
  desc, out = run(HT7, 'cl_ht7')
  assert 'Peptide bond across letters: " C  ADNE G   2 " - " N  CHT7 G   3 "' \
    in passes(desc)['C'], passes(desc)['C']
  s = sites(out)
  hs = dict((k[1], v) for k, v in s.items() if k[0] == 3 and k[2] == 'C' and
    k[1].startswith('H'))
  assert sorted(hs) == sorted(MASTER_HT7_C), sorted(hs)
  for k in hs:
    assert math.dist(hs[k], MASTER_HT7_C[k]) < 0.002, (k, hs[k], MASTER_HT7_C[k])

if __name__ == "__main__":
  test_peptide()
  test_place()
  test_reverse()
  test_riding_entry()
  test_optimizer_alone()
  test_stretched()
  test_controls()
  test_links_not_peptide()
  test_neutron()
  test_optimize()
  test_reference_h_missing()
  test_not_ordinary_amide()
  print("OK")
