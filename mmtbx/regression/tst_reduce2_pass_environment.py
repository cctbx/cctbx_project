"""reduce2 Optimizer: what an altloc pass holds and what it writes (A33).
In the pass for X every residue gives ResidueAtomsForConformer(): its blank
atoms; its X atoms; without X, one other letter only if all its atoms have
occupancy >= 0.99, several other letters their alphabetically first; the same
in every chain. A residue group whose site (chain object, resseq, icode) has
another group with X gives its blank atoms only. A residue with altlocs but without X is the pass's environment
(EnvironmentResidues()): it starts from its initial state, gets Movers, and is
put back after the pass; the pass records nothing for it (deletions, flips,
His records, report lines). So every atom is written by its own pass, a blank
atom of a residue with letters by the last pass that has one of them.
Fixtures: the selection table (chain P: lead Gly 1 sets the chain's letters,
probe Ala 2; chain Z carries the pass letter B); 1us0 His A 306 (blank ring,
N/CA B/C) with a Zn in A only; 1us0 Asn A 7 with a
blank side chain and N/CA B/C; a Cys with SG B/C and a water in A only; the
1ejg residues within 6 A of Ser A 22 OG (Pro 22 A and Ser 22 B/C are separate
residue groups); a Ser in B only at occupancy 1.00 (6lk1 Thr A 75/81, Ser A
82); 1exr His A 107 with a C conformer elsewhere (His charge reruns).
"""
import os, re, math
from libtbx.utils import null_out
from iotbx.cli_parser import run_program
import iotbx.pdb
from mmtbx.programs import reduce2
from mmtbx.reduce import Optimizers

def run(model_text, prefix, *extra):
  """reduce2 approach=add; returns (description, output hierarchy)"""
  fn_in, fn_out, fn_desc = prefix + ".pdb", prefix + "_out.pdb", prefix + "_out.txt"
  with open(fn_in, "w") as f:
    f.write(model_text)
  for fn in (fn_out, fn_desc):
    if os.path.exists(fn): os.remove(fn)
  run_program(program_class=reduce2.Program, logger=null_out(),
    args=[fn_in, "approach=add", "output.filename=" + fn_out,
      "output.description_file_name=" + fn_desc, "overwrite=True"] + list(extra))
  with open(fn_desc) as f:
    desc = f.read()
  return desc, iotbx.pdb.input(fn_out).construct_hierarchy()

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

def report(part):
  """the result lines of a pass's report"""
  return [l for l in part.split('BEGIN REPORT')[1].split('END REPORT')[0].split('\n')
    if 'final score' in l]

def h_sites(h):
  """{(chain, resseq+icode, altloc, name, resname): xyz} of the H"""
  return dict(((a.parent().parent().parent().id, a.parent().parent().resseq.strip() +
    a.parent().parent().icode.strip(), a.parent().altloc.strip(), a.name.strip(),
    a.parent().resname.strip()), a.xyz)
    for a in h.models()[0].atoms() if a.element_is_hydrogen())

def dist(h, resseq, alt, n1, n2):
  xyz = {}
  for rg in h.models()[0].residue_groups():
    if rg.resseq.strip() != str(resseq): continue
    for ag in rg.atom_groups():
      if ag.altloc.strip() in ('', alt):
        for a in ag.atoms():
          xyz[a.name.strip()] = a.xyz
  return math.dist(xyz[n1], xyz[n2])

def written_by_own_pass(model_text, prefix, *extra):
  """Every H as the single pass of its writer (alt_id=) has it: its own letter;
  a blank H of a residue with letters, the residue's first letter (the last of
  its passes); a residue without letters, the last pass. Returns the full
  run's description."""
  desc, h = run(model_text, prefix, *extra)
  full = h_sites(h)
  model = iotbx.pdb.input(source_info=None, lines=model_text).construct_hierarchy().models()[0]
  last = sorted(x for x in Optimizers.AlternatesInModel(model) if x.strip())[0]
  first = {}
  for rg in model.residue_groups():
    for ag in rg.atom_groups():
      k = (rg.parent().id, rg.resseq.strip() + rg.icode.strip(), ag.resname.strip())
      first[k] = min([first.get(k, 'zz')] + Optimizers.ResidueLetters(rg))
  single = {}
  for k in full:
    w = k[2] or first.get((k[0], k[1], k[4]), 'zz')
    w = last if w == 'zz' else w
    if w not in single:
      single[w] = h_sites(run(model_text, prefix + '_' + w, 'alt_id=' + w, *extra)[1])
    assert k in single[w], (k, w)
    assert math.dist(full[k], single[w][k]) < 0.0015, (k, w, full[k], single[w][k])
  return desc

# --- the selection table --------------------------------------------------------

CLASSES = {
  'blank': [('', 1.0)],
  'blank+1': [('', 1.0), ('Y', 0.5)],
  'blank+1 full': [('', 1.0), ('Y', 1.0)],
  'blank+2': [('', 1.0), ('Y', 0.5), ('W', 0.5)],
  '1': [('Y', 0.5)],
  '1 full': [('Y', 1.0)],
  '2': [('Y', 0.5), ('W', 0.5)],
}
# chain P's letters (lead residue) and the probe's letters
COLUMNS = {
  'has X': (['A', 'B'], {'Y': 'A', 'W': 'C'}),
  'lacks X, L in residue': (['A'], {'Y': 'A', 'W': 'C'}),
  'lacks X, L not in residue': (['A'], {'Y': 'C', 'W': 'D'}),
}

def table_fixture(groups, leadAlts, mapping, ter=True, reverse=False):
  s = []
  n = [0]
  def add(name, alt, resname, chain, resseq, occ):
    n[0] += 1
    s.append("ATOM  %5d  %-3s%1s%3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f           %s\n" % (
      n[0], name, alt, resname, chain, resseq, 1.5 * n[0], 0.0, 0.0, occ, 20.0, name[0]))
  for nm in ('N', 'C', 'O'):
    add(nm, '', 'GLY', 'P', 1, 1.0)
  for a in leadAlts:
    add('CA', a, 'GLY', 'P', 1, 1.0 / len(leadAlts))
  blank = any(g[0] == '' for g in groups)
  for alt, occ in (reversed(groups) if reverse else groups):
    if alt == '':
      for nm in ('N', 'CA', 'C', 'O'):
        add(nm, '', 'ALA', 'P', 2, occ)
    else:
      for nm in (('CB',) if blank else ('N', 'CA', 'C', 'O', 'CB')):
        add(nm, mapping[alt], 'ALA', 'P', 2, occ)
  if ter:
    s.append("TER\n")
  for a in ('A', 'B'):
    add('CA', a, 'GLY', 'Z', 1, 0.5)
  return ''.join(s)

def probe(text, conf='B'):
  m = iotbx.pdb.input(source_info=None, lines=text).construct_hierarchy().models()[0]
  return ' '.join(a.name.strip() + ('/' + a.parent().altloc if a.parent().altloc.strip() else '')
    for a in Optimizers.GetAtomsForConformer(m, conf)
    if a.parent().parent().parent().id == 'P' and a.parent().parent().resseq.strip() == '2')

def test_table():
  """Pass B, probe residue Ala 2. '*': not what master gave (master: a chain
  with B gave the residue's blank atoms only, a chain without B its first
  conformer, minus single-letter residues below 0.99)."""
  bb = 'N CA C O'
  full = lambda x: 'N/%s CA/%s C/%s O/%s CB/%s' % ((x,) * 5)
  expected = {
    ('blank', 'has X'): bb,
    ('blank', 'lacks X, L in residue'): bb,
    ('blank', 'lacks X, L not in residue'): bb,
    ('blank+1', 'has X'): bb,
    ('blank+1', 'lacks X, L in residue'): bb,                  # * blank + A
    ('blank+1', 'lacks X, L not in residue'): bb,
    ('blank+1 full', 'has X'): bb + ' CB/A',                   # * blank
    ('blank+1 full', 'lacks X, L in residue'): bb + ' CB/A',
    ('blank+1 full', 'lacks X, L not in residue'): bb + ' CB/C',  # * blank
    ('blank+2', 'has X'): bb + ' CB/A',                        # * blank
    ('blank+2', 'lacks X, L in residue'): bb + ' CB/A',
    ('blank+2', 'lacks X, L not in residue'): bb + ' CB/C',    # * blank
    ('1', 'has X'): '',
    ('1', 'lacks X, L in residue'): '',
    ('1', 'lacks X, L not in residue'): '',
    ('1 full', 'has X'): full('A'),                            # * nothing
    ('1 full', 'lacks X, L in residue'): full('A'),
    ('1 full', 'lacks X, L not in residue'): full('C'),        # * nothing
    ('2', 'has X'): full('A'),                                 # * nothing
    ('2', 'lacks X, L in residue'): full('A'),
    ('2', 'lacks X, L not in residue'): full('C'),             # * nothing
  }
  for (cls, col), want in sorted(expected.items()):
    lead, mapping = COLUMNS[col]
    for ter in (True, False):
      got = probe(table_fixture(CLASSES[cls], lead, mapping, ter=ter))
      assert got == want, (cls, col, ter, got, want)

def test_first_letter_alphabetical():
  """Letters given B before A in the file: pass C holds A, as the passes run
  C, B, A."""
  text = table_fixture([('W', 0.5), ('Y', 0.5)], ['A', 'B'], {'Y': 'A', 'W': 'B'})
  text = text.replace(' CA BGLY Z   1', ' CA CGLY Z   1')    # pass C
  assert probe(text, 'C') == 'N/A CA/A C/A O/A CB/A', probe(text, 'C')

# --- written by the own pass --------------------------------------------------

m_his306_zn = """\
CRYST1   49.280   66.590   47.260  90.00  92.40  90.00 P 1 21 1
ATOM    119  C   HIS A 306      14.677  -0.747   1.918  1.00  5.64           C
ATOM    120  O   HIS A 306      15.585  -0.087   2.444  1.00  6.43           O
ATOM    121  CB  HIS A 306      14.146  -2.665   3.467  1.00  5.27           C
ATOM    122  CG  HIS A 306      12.697  -2.414   3.685  1.00  4.67           C
ATOM    123  ND1 HIS A 306      12.135  -1.146   3.697  1.00  4.78           N
ATOM    124  CD2 HIS A 306      11.710  -3.317   3.918  1.00  4.82           C
ATOM    125  CE1 HIS A 306      10.859  -1.311   3.914  1.00  4.70           C
ATOM    126  NE2 HIS A 306      10.554  -2.615   4.056  1.00  4.79           N
ATOM    127  N  BHIS A 306      16.032  -2.645   1.879  0.62  5.79           N
ATOM    128  CA BHIS A 306      14.642  -2.250   2.067  0.62  5.62           C
ATOM    129  N  CHIS A 306      15.963  -2.917   1.965  0.38  5.83           N
ATOM    130  CA CHIS A 306      14.666  -2.254   2.070  0.38  5.27           C
HETATM  999 ZN  A ZN A 999       9.329   0.402   4.038  0.70 10.00          ZN
"""

def test_shared_ring_ion_lock():
  """Pass A holds His 306 as B (environment) and its Zn locks the blank ring
  flipped there; the ring is put back, so the ring H of B and C stay bonded
  (with only the borrowed B atoms put back: 3.0, 3.2, 2.0 A)."""
  desc, h = run(m_his306_zn, 'tst_reduce2_pass_environment_his306', 'add_flip_movers=True')
  for alt in ('B', 'C'):
    for hn, n, ideal in (('HD1', 'ND1', 0.86), ('HD2', 'CD2', 0.93), ('HE1', 'CE1', 0.93)):
      d = dist(h, 306, alt, hn, n)
      assert abs(d - ideal) < 0.01, (alt, hn, n, d)
  p = passes(desc)
  assert 'Set MoverHisFlip on chain A HIS 306 to state' in p['A'], p['A']
  assert not report(p['A']), report(p['A'])
  assert ' Environment Movers (residues without alternate \'A\'; written by their own passes): 1' in p['A']

m_asn7 = """\
CRYST1   49.280   66.590   47.260  90.00  92.40  90.00 P 1 21 1
ATOM     55  N  BASN A   7       5.327  13.486  26.043  0.62  5.73           N
ATOM     55  N  CASN A   7       5.377  13.486  26.043  0.38  5.73           N
ATOM     56  CA BASN A   7       4.427  14.403  25.339  0.62  5.82           C
ATOM     56  CA CASN A   7       4.477  14.403  25.339  0.38  5.82           C
ATOM     57  C   ASN A   7       3.007  13.907  25.217  1.00  5.86           C
ATOM     58  O   ASN A   7       2.154  14.594  24.639  1.00  7.05           O
ATOM     59  CB  ASN A   7       5.011  14.800  23.994  1.00  5.69           C
ATOM     60  CG  ASN A   7       5.113  13.637  23.035  1.00  5.20           C
ATOM     61  OD1 ASN A   7       4.574  12.559  23.328  1.00  5.45           O
ATOM     62  ND2 ASN A   7       5.775  13.847  21.904  1.00  5.97           N
HETATM  900  O  AHOH A 401      25.000  25.000  25.000  0.50 10.00           O
HETATM  901  O  BHOH A 401      25.300  25.000  25.000  0.50 10.00           O
"""

def test_shared_amide_flip():
  """Asn 7: side chain blank, N/CA B/C, the flip set for every pass. Before,
  pass A (chain with A, the residue without it) held the blank side chain
  without H, reset it unflipped and left the flipped H of B and C 2.4-3.0 A
  from ND2."""
  desc, h = run(m_asn7, 'tst_reduce2_pass_environment_asn7', 'add_flip_movers=True',
    'set_flip_states=1 b A ASN 7 Flipped AnglesNotAdjusted')
  for alt in ('B', 'C'):
    for hn in ('HD21', 'HD22'):
      d = dist(h, 7, alt, hn, 'ND2')
      assert d < 1.1, (alt, hn, d)

m_cys_w = """\
CRYST1  108.910  108.910  108.910  90.00  90.00  90.00 I 4 3 2
ATOM      1  N   CYS A  48     -40.561 -16.057   6.830  1.00 99.72           N
ATOM      2  CA  CYS A  48     -41.351 -14.915   6.380  1.00106.82           C
ATOM      3  C   CYS A  48     -41.146 -13.719   7.305  1.00104.36           C
ATOM      4  O   CYS A  48     -40.407 -13.794   8.283  1.00104.56           O
ATOM      5  CB  CYS A  48     -42.835 -15.275   6.303  1.00105.79           C
ATOM      6  SG BCYS A  48     -43.707 -15.186   7.879  0.50101.00           S
ATOM      7  SG CCYS A  48     -43.407 -14.886   7.979  0.50101.00           S
TER
HETATM    8  O  AHOH W   1     -46.501 -16.214   7.508  0.70 20.00           O
"""

def test_environment_not_written():
  """Cys 48 (SG B/C) is environment in pass A, whose water (A only, 0.70) is
  an acceptor of HG/B there: HG/B as pass B placed it, and pass A reports
  nothing for the Cys (before: pass A's pose, 0.002 A off)."""
  desc = written_by_own_pass(m_cys_w, 'tst_reduce2_pass_environment_cys')
  p = passes(desc)
  assert [l for l in report(p['B']) if 'CYS 48 HG' in l], p['B']
  assert not [l for l in report(p['A']) if 'CYS 48' in l], p['A']
  assert re.search(r"Environment Movers \(residues without alternate 'A'; written by their own passes\): \d+", p['A']), p['A']

def test_environment_starts_initial():
  """In pass A the environment (Cys 48 B) starts from its initial state, as in
  its own pass B, not from where pass B left it."""
  seen = {}
  orig = Optimizers.Optimizer._PlaceMovers
  def place(self, atoms, *a, **k):
    alt = a[-1]
    for x in atoms:
      if x.name.strip() == 'HG' and x.parent().altloc == 'B':
        seen[alt] = x.xyz
    return orig(self, atoms, *a, **k)
  Optimizers.Optimizer._PlaceMovers = place
  try:
    desc, h = run(m_cys_w, 'tst_reduce2_pass_environment_cys_start')
  finally:
    Optimizers.Optimizer._PlaceMovers = orig
  assert sorted(seen) == ['A', 'B'], seen
  assert seen['A'] == seen['B'], seen
  final = [a.xyz for a in h.atoms() if a.name.strip() == 'HG' and a.parent().altloc == 'B'][0]
  assert math.dist(final, seen['B']) > 0.01, (final, seen)   # pass B moved it

m_1ejg_22 = """\
CRYST1   40.824   18.498   22.371  90.00  90.47  90.00 P 1 21 1
ATOM    315  N   ARG A  17       4.980   6.166  -0.135  1.00  1.46           N
ATOM    316  CA  ARG A  17       5.850   6.474  -1.244  1.00  1.51           C
ATOM    317  C   ARG A  17       5.641   5.554  -2.455  1.00  1.61           C
ATOM    318  O   ARG A  17       6.033   5.932  -3.563  1.00  2.20           O
ATOM    319  CB  ARG A  17       7.317   6.413  -0.802  1.00  1.71           C
ATOM    320  CG  ARG A  17       7.692   7.571   0.124  1.00  1.91           C
ATOM    321  CD  ARG A  17       7.903   8.893  -0.615  1.00  2.03           C
ATOM    322  NE  ARG A  17       9.156   8.833  -1.364  1.00  2.13           N
ATOM    323  CZ  ARG A  17       9.472   9.600  -2.396  1.00  1.95           C
ATOM    324  NH1 ARG A  17       8.585  10.416  -2.966  1.00  2.59           N
ATOM    325  NH2 ARG A  17      10.713   9.561  -2.889  1.00  2.55           N
ATOM    386  N   THR A  21       4.056   9.404  -4.216  1.00  1.67           N
ATOM    387  CA  THR A  21       3.664  10.576  -3.461  1.00  1.43           C
ATOM    388  C   THR A  21       4.915  11.347  -3.001  1.00  1.58           C
ATOM    389  O   THR A  21       5.844  10.729  -2.478  1.00  2.25           O
ATOM    390  CB  THR A  21       2.932  10.113  -2.178  1.00  1.52           C
ATOM    391  OG1 THR A  21       1.964   9.124  -2.552  1.00  1.89           O
ATOM    392  CG2 THR A  21       2.280  11.274  -1.459  1.00  2.12           C
ATOM    400  N   PRO A  22       4.915  12.683  -3.102  1.00  1.83           N
ATOM    401  CA APRO A  22       6.042  13.429  -2.601  0.57  1.82           C
ATOM    402  C  APRO A  22       6.387  13.122  -1.160  0.57  1.66           C
ATOM    403  O  APRO A  22       5.480  13.006  -0.345  0.57  2.09           O
ATOM    404  CB APRO A  22       5.655  14.896  -2.744  0.57  2.86           C
ATOM    405  CG APRO A  22       4.661  14.854  -4.058  0.57  2.66           C
ATOM    406  CD APRO A  22       3.957  13.505  -3.910  0.57  2.27           C
ATOM    414  CA BSER A  22       6.034  13.399  -2.687  0.21  1.55           C
ATOM    415  C  BSER A  22       6.367  13.062  -1.223  0.21  2.92           C
ATOM    416  O  BSER A  22       5.412  13.050  -0.345  0.21  1.87           O
ATOM    417  CB BSER A  22       5.409  14.835  -2.876  0.21  2.60           C
ATOM    418  OG BSER A  22       4.760  15.243  -1.635  0.21  2.11           O
ATOM    420  CA CSER A  22       6.112  13.653  -2.656  0.22  2.31           C
ATOM    421  C  CSER A  22       6.354  13.275  -1.187  0.22  1.92           C
ATOM    422  O  CSER A  22       5.636  12.705  -0.270  0.22  1.77           O
ATOM    423  CB CSER A  22       5.605  15.097  -2.687  0.22  3.56           C
ATOM    424  OG CSER A  22       6.750  15.771  -2.280  0.22  6.38           O
ATOM    426  N   GLU A  23       7.651  13.190  -0.860  1.00  1.86           N
ATOM    427  CA  GLU A  23       8.108  13.002   0.513  1.00  1.81           C
ATOM    428  C   GLU A  23       7.394  13.968   1.459  1.00  1.66           C
ATOM    429  O   GLU A  23       7.064  13.585   2.589  1.00  1.97           O
ATOM    430  CB  GLU A  23       9.620  13.168   0.681  1.00  2.12           C
ATOM    435  CG AGLU A  23      10.459  12.148  -0.048  0.58  2.73           C
ATOM    436  CD AGLU A  23      11.942  12.208   0.242  0.58  4.60           C
ATOM    437  OE1AGLU A  23      12.418  13.330   0.541  0.58  7.57           O
ATOM    438  OE2AGLU A  23      12.629  11.211   0.128  0.58  7.42           O
ATOM    441  CG BGLU A  23      10.414  11.969   0.088  0.42  1.79           C
ATOM    442  CD BGLU A  23      11.911  12.048   0.452  0.42  2.33           C
ATOM    443  OE1BGLU A  23      12.430  13.113   0.795  0.42  2.70           O
ATOM    444  OE2BGLU A  23      12.514  10.945   0.316  0.42  4.43           O
ATOM    447  N   ALA A  24       7.188  15.227   1.041  1.00  1.63           N
ATOM    448  CA  ALA A  24       6.619  16.195   1.970  1.00  1.65           C
ATOM    449  C   ALA A  24       5.192  15.819   2.387  1.00  1.51           C
ATOM    450  O   ALA A  24       4.795  16.037   3.532  1.00  1.85           O
ATOM    451  CB  ALA A  24       6.663  17.597   1.363  1.00  2.43           C
ATOM    457  N  ALEU A  25       4.398  15.289   1.441  0.57  1.61           N
ATOM    458  CA ALEU A  25       3.048  14.822   1.781  0.57  1.48           C
ATOM    459  C  ALEU A  25       3.113  13.568   2.640  0.57  1.50           C
ATOM    460  O  ALEU A  25       2.332  13.396   3.580  0.57  2.01           O
ATOM    461  CB ALEU A  25       2.209  14.621   0.512  0.57  1.70           C
ATOM    462  CG ALEU A  25       1.717  15.929  -0.108  0.57  1.80           C
ATOM    463  CD1ALEU A  25       1.207  15.625  -1.514  0.57  3.58           C
ATOM    464  CD2ALEU A  25       0.595  16.583   0.794  0.57  2.15           C
ATOM    476  N  BILE A  25       4.398  15.289   1.441  0.21  1.61           N
ATOM    477  CA BILE A  25       3.048  14.822   1.781  0.21  1.48           C
ATOM    478  C  BILE A  25       3.113  13.568   2.640  0.21  1.50           C
ATOM    479  O  BILE A  25       2.332  13.396   3.580  0.21  2.01           O
ATOM    480  CB BILE A  25       2.225  14.590   0.511  0.21  1.67           C
ATOM    481  CG1BILE A  25       1.931  15.856  -0.269  0.21  2.53           C
ATOM    482  CG2BILE A  25       0.955  13.756   0.691  0.21  2.16           C
ATOM    483  CD1BILE A  25       0.819  16.720   0.588  0.21  1.98           C
ATOM    495  N  CILE A  25       4.398  15.289   1.441  0.22  1.61           N
ATOM    496  CA CILE A  25       3.048  14.822   1.781  0.22  1.48           C
ATOM    497  C  CILE A  25       3.113  13.568   2.640  0.22  1.50           C
ATOM    498  O  CILE A  25       2.332  13.396   3.580  0.22  2.01           O
ATOM    499  CB CILE A  25       2.194  14.590   0.526  0.22  2.08           C
ATOM    500  CG1CILE A  25       1.693  15.919  -0.061  0.22  2.48           C
ATOM    501  CG2CILE A  25       0.939  13.672   0.877  0.22  2.83           C
ATOM    502  CD1CILE A  25       1.188  15.934  -1.404  0.22  3.60           C
ATOM    511  N   CYS A  26       4.028  12.634   2.294  1.00  1.57           N
ATOM    512  CA  CYS A  26       4.137  11.461   3.154  1.00  1.52           C
ATOM    513  C   CYS A  26       4.607  11.823   4.564  1.00  1.58           C
ATOM    514  O   CYS A  26       4.231  11.137   5.527  1.00  2.03           O
ATOM    515  CB  CYS A  26       5.038  10.395   2.541  1.00  1.82           C
ATOM    516  SG  CYS A  26       4.363   9.697   1.004  1.00  2.01           S
"""

def test_microheterogeneity():
  """1ejg Pro A 22 / Ser B,C 22: the Ser is environment in pass A, which runs
  last (the rule without ownership: HG/B 0.12 A as pass A placed it); every H
  as its own pass places it."""
  desc = written_by_own_pass(m_1ejg_22, 'tst_reduce2_pass_environment_1ejg')
  assert not [l for l in report(passes(desc)['A']) if 'SER 22' in l], passes(desc)['A']

def atoms22(m, conf, resseq='22'):
  """(resname, altloc, name) of residue resseq's atoms in the pass for conf"""
  return [(a.parent().resname.strip(), a.parent().altloc.strip(), a.name.strip())
    for a in Optimizers.GetAtomsForConformer(m, conf)
    if a.parent().parent().resseq.strip() + a.parent().parent().icode.strip() == resseq]

def test_same_site_1ejg():
  """1ejg residue 22: Pro (blank N + A) and Ser (B, C) are separate residue
  groups at the same site. A residue group without the pass's letter gives only
  its blank atoms when another group at its site has the letter: no Ser in pass
  A; the shared N (stored in the Pro group) in every pass; Pro's A atoms only
  in pass A."""
  m = iotbx.pdb.input(source_info=None, lines=m_1ejg_22).construct_hierarchy().models()[0]
  assert len([rg for rg in m.residue_groups() if rg.resseq.strip() == '22']) == 2
  pro = [('PRO', 'A', n) for n in ('CA', 'C', 'O', 'CB', 'CG', 'CD')]
  ser = lambda x: [('SER', x, n) for n in ('CA', 'C', 'O', 'CB', 'OG')]
  assert atoms22(m, 'A') == [('PRO', '', 'N')] + pro, atoms22(m, 'A')
  assert atoms22(m, 'B') == [('PRO', '', 'N')] + ser('B'), atoms22(m, 'B')
  assert atoms22(m, 'C') == [('PRO', '', 'N')] + ser('C'), atoms22(m, 'C')
  # the Ser group is still environment in pass A (letters, none of them A)
  env = Optimizers.EnvironmentResidues(m, 'A', ['C', 'B', 'A'])
  assert ('22', 'SER') in [(rg.resseq.strip(), rg.atom_groups()[0].resname) for rg in env]

m_same_site = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
ATOM      1  N   PRO A   7      10.000  10.000  10.000  1.00 10.00           N
ATOM      2  CA APRO A   7      11.400  10.000  10.000  0.60 10.00           C
ATOM      3  CA BSER A   7A     11.400  10.100  10.000  0.20 10.00           C
ATOM      4  CA CSER A   7A     11.400  10.200  10.000  0.20 10.00           C
ATOM      5  N   PRO A   8      14.000  10.000  10.000  1.00 10.00           N
ATOM      6  CA APRO A   8      15.400  10.000  10.000  0.60 10.00           C
ATOM      7  CA BSER A   8      15.400  10.100  10.000  0.20 10.00           C
ATOM      8  CA CSER A   8      15.400  10.200  10.000  0.20 10.00           C
ATOM      9  CA APRO A   9      18.000  10.000  10.000  0.60 10.00           C
TER
ATOM     10  CA  GLY B   1      20.000  20.000  20.000  1.00 10.00           C
TER
HETATM   11  O  BHOH A   9      18.000  10.100  10.000  0.30 10.00           O
HETATM   12  O  CHOH A   9      18.000  10.200  10.000  0.30 10.00           O
"""

def test_same_site_icode_and_chains():
  """The same-site exception needs the same chain object, resseq and icode:
  Ser B,C 7A next to Pro A 7 (other icode) and HOH B,C 9 in a second chain
  object A next to Pro A 9 give their first letter in pass A; Ser B,C 8 at
  Pro A 8's site gives nothing."""
  m = iotbx.pdb.input(source_info=None, lines=m_same_site).construct_hierarchy().models()[0]
  assert [c.id for c in m.chains()] == ['A', 'B', 'A']
  assert len([rg for rg in m.chains()[0].residue_groups() if rg.resseq.strip() == '8']) == 2
  assert atoms22(m, 'A', '7A') == [('SER', 'B', 'CA')], atoms22(m, 'A', '7A')
  assert atoms22(m, 'A', '8') == [('PRO', '', 'N'), ('PRO', 'A', 'CA')], atoms22(m, 'A', '8')
  assert atoms22(m, 'A', '9') == [('PRO', 'A', 'CA'), ('HOH', 'B', 'O')], atoms22(m, 'A', '9')
  assert atoms22(m, 'B', '8') == [('PRO', '', 'N'), ('SER', 'B', 'CA')], atoms22(m, 'B', '8')
  assert atoms22(m, 'C', '7') == [('PRO', '', 'N')], atoms22(m, 'C', '7')

class _OwnMover:
  def CoarsePositions(self):
    from types import SimpleNamespace
    return SimpleNamespace(atoms=[SimpleNamespace(i_seq=2)])

def test_his_cliques_whole():
  """His charge reruns: the pass's Cliques are recorded whole, environment
  Movers included, while the scores are the pass's own. An owned neighbour's
  loss next to a pinned environment His (same Mover id as in its own pass) is
  then the His neighbourhood's cost, not unattributed."""
  from types import SimpleNamespace
  from mmtbx.reduce import Movers, HisCharges
  his = Movers.MoverHisFlip.__new__(Movers.MoverHisFlip)
  his.CoarsePositions = lambda: SimpleNamespace(atoms=[SimpleNamespace(i_seq=1)])
  own = _OwnMover()
  opt = SimpleNamespace(_movers=[his, own], _envMovers=set([his]), _envISeqs=set([1]),
    _finalScore={his: 1., own: 1.}, _environmentAtoms=[], _atoms=[])
  rec = {'passes': [], 'cliques': []}
  Optimizers.Optimizer._HisRecordPass(opt, 0, 'C', {}, None, rec, [[his, own]])
  assert rec['passes'] == [('C', {(2,): 1.})], rec
  assert rec['cliques'] == [set([(1,), (2,)])], rec
  key = (0, 'A', '107', '')
  before = [('C', {(2,): 1.}), ('A', {(1,): 1.})]
  after = [('C', {(2,): .9}), ('A', {(1,): 1.})]
  t = HisCharges.trial(before, after, {key: set([(1,)])}, rec['cliques'], [key], {})
  assert not t['unattributed'], t
  assert abs(t['components'][0]['cost']['total'] - 0.1) < 1e-6, t

m_full_b = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00 10.00           N
ATOM      2  CA  ALA A   1      11.458  10.000  10.000  1.00 10.00           C
ATOM      3  C   ALA A   1      12.009  11.420  10.000  1.00 10.00           C
ATOM      4  O   ALA A   1      11.249  12.390  10.000  1.00 10.00           O
ATOM      5  CB AALA A   1      11.988   9.236  11.200  0.50 10.00           C
ATOM      6  CB BALA A   1      11.988   9.236   8.800  0.50 10.00           C
ATOM      7  N  BSER A   2      13.330  11.540  10.000  1.00 10.00           N
ATOM      8  CA BSER A   2      13.980  12.840  10.000  1.00 10.00           C
ATOM      9  C  BSER A   2      15.490  12.700  10.000  1.00 10.00           C
ATOM     10  O  BSER A   2      16.050  11.600  10.000  1.00 10.00           O
ATOM     11  CB BSER A   2      13.530  13.640  11.220  1.00 10.00           C
ATOM     12  OG BSER A   2      14.020  14.970  11.180  1.00 10.00           O
"""

def test_full_occupancy_in_every_pass():
  """6lk1's case: Ser 2 labelled B only at occupancy 1.00, in a chain with A,
  is in pass A (master: left out) as environment; every H as its own pass
  places it."""
  m = iotbx.pdb.input(source_info=None, lines=m_full_b).construct_hierarchy().models()[0]
  got = [a.name.strip() for a in Optimizers.GetAtomsForConformer(m, 'A')
    if a.parent().parent().resseq.strip() == '2']
  assert got == ['N', 'CA', 'C', 'O', 'CB', 'OG'], got
  env = Optimizers.EnvironmentResidues(m, 'A', ['B', 'A'])
  assert [rg.resseq.strip() for rg in env] == ['2'], env
  assert Optimizers.EnvironmentResidues(m, 'B', ['B', 'A']) == []
  assert Optimizers.EnvironmentResidues(m, 'A', ['A']) == []    # alt_id=A: all its own
  desc = written_by_own_pass(m_full_b, 'tst_reduce2_pass_environment_full')
  assert [l for l in report(passes(desc)['B']) if 'SER 2 OG' in l or 'SER 2 HG' in l], desc
  assert not [l for l in report(passes(desc)['A']) if 'SER 2' in l], passes(desc)['A']

def test_his_rerun_environment():
  """His charge reruns (1exr His A 107, A/B) with a C conformer elsewhere: the
  His is environment in pass C, which records no His for it; the result is
  the one without the C conformer."""
  from mmtbx.regression.tst_reduce2_his_charges_optimizer import pdb_1exr
  lines = [l for l in pdb_1exr.split('\n') if l.startswith(('ATOM', 'HETATM'))]
  far = ("HETATM 9001  O  BHOH Z   1      60.000  60.000  60.000  0.50 20.00           O\n"
         "HETATM 9002  O  CHOH Z   1      60.300  60.000  60.000  0.50 20.00           O\n")
  results = []
  for name, text in (('ab', pdb_1exr), ('abc', pdb_1exr.rstrip('\n') + '\n' + far)):
    made = []
    orig = Optimizers.Optimizer.__init__
    def init(self, *a, **k):
      made.append(self)
      orig(self, *a, **k)
    Optimizers.Optimizer.__init__ = init
    try:
      desc, h = run(text, 'tst_reduce2_pass_environment_1exr_' + name, 'reconcile_his_charges=True')
    finally:
      Optimizers.Optimizer.__init__ = orig
    his = made[0].getHisCharges()
    assert len(his) == 1 and len(his[0]['cases']) == 1, his
    c = his[0]['cases'][0]
    assert c['free'] == {'A': 1, 'B': 0} and c['final'] == {'A': 0, 'B': 0}, c
    assert list(c['free']) == ['A', 'B'], c['free']
    results.append((his[0]['reruns'], c['final'],
      dict((k, v) for k, v in h_sites(h).items() if k[0] != 'Z')))
    if name == 'abc':
      p = passes(desc)
      assert not [l for l in report(p['C']) if 'HIS 107' in l], p['C']
      assert 'Environment Movers' in p['C'], p['C']
  assert results[0][:2] == results[1][:2], results
  assert sorted(results[0][2]) == sorted(results[1][2])
  for k in results[0][2]:
    assert math.dist(results[0][2][k], results[1][2][k]) < 0.0015, k

if __name__ == "__main__":
  test_table()
  test_first_letter_alphabetical()
  test_shared_ring_ion_lock()
  test_shared_amide_flip()
  test_environment_not_written()
  test_environment_starts_initial()
  test_microheterogeneity()
  test_same_site_1ejg()
  test_same_site_icode_and_chains()
  test_his_cliques_whole()
  test_full_occupancy_in_every_pass()
  test_his_rerun_environment()
  print("OK")
