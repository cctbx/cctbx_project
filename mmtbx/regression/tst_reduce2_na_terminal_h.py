"""reduce2 add and place put the hydroxyl H on free nucleic-acid ends: HO5' on
an O5' without phosphate, HO3' on an O3' that is not linked, on the first and
last nucleotide of a chain object only (a TER or an mmCIF chain change ends
one); a standard nucleotide all HETATM and unlinked is its own first and last.
A free end at a chain break inside a chain object gets no H and is reported;
qmready puts H on every free end. Fragments: 5gsk A 1-3 (DNA), 1rrr A 1-2
(RNA).
"""
from __future__ import absolute_import, division, print_function
import os
import iotbx.cif
import iotbx.pdb
from six.moves import cStringIO as StringIO
from iotbx.cli_parser import run_program
from mmtbx.programs import reduce2

CAVEAT = 'does not show whether a phosphate was present'
GAP_CAVEAT = 'do not establish experimental termini'
NEUTRON_NOTE = ('The listed names are H labels from placement, before any '
  'requested isotope conversion.')

dna = """\
ATOM      1  O5'  DG A   1      14.243  20.552  -5.680  1.00 21.12           O
ATOM      2  C5'  DG A   1      14.519  21.926  -5.348  1.00 20.86           C
ATOM      3  C4'  DG A   1      13.391  22.547  -4.553  1.00 20.27           C
ATOM      4  O4'  DG A   1      12.161  22.635  -5.311  1.00 18.57           O
ATOM      5  C3'  DG A   1      13.012  21.844  -3.255  1.00 19.81           C
ATOM      6  O3'  DG A   1      13.903  22.237  -2.215  1.00 21.62           O
ATOM      7  C2'  DG A   1      11.586  22.313  -3.035  1.00 18.93           C
ATOM      8  C1'  DG A   1      11.045  22.370  -4.456  1.00 18.30           C
ATOM      9  N9   DG A   1      10.411  21.138  -4.910  1.00 18.13           N
ATOM     10  C8   DG A   1      10.965  20.144  -5.682  1.00 18.17           C
ATOM     11  N7   DG A   1      10.132  19.173  -5.940  1.00 18.20           N
ATOM     12  C5   DG A   1       8.953  19.554  -5.313  1.00 17.34           C
ATOM     13  C6   DG A   1       7.690  18.909  -5.256  1.00 17.10           C
ATOM     14  O6   DG A   1       7.353  17.833  -5.757  1.00 18.28           O
ATOM     15  N1   DG A   1       6.771  19.646  -4.516  1.00 17.32           N
ATOM     16  C2   DG A   1       7.031  20.851  -3.908  1.00 18.67           C
ATOM     17  N2   DG A   1       6.014  21.408  -3.232  1.00 19.76           N
ATOM     18  N3   DG A   1       8.205  21.460  -3.951  1.00 18.17           N
ATOM     19  C4   DG A   1       9.108  20.767  -4.677  1.00 17.68           C
ATOM     20  P    DG A   2      14.088  21.308  -0.949  1.00 22.88           P
ATOM     21  OP1  DG A   2      15.234  21.838  -0.174  1.00 26.15           O
ATOM     22  OP2  DG A   2      14.113  19.892  -1.404  1.00 23.55           O
ATOM     23  O5'  DG A   2      12.760  21.605  -0.127  1.00 21.40           O
ATOM     24  C5'  DG A   2      12.130  20.581   0.642  1.00 21.02           C
ATOM     25  C4'  DG A   2      10.669  20.918   0.812  1.00 19.22           C
ATOM     26  O4'  DG A   2      10.024  20.958  -0.480  1.00 20.06           O
ATOM     27  C3'  DG A   2       9.881  19.895   1.618  1.00 18.59           C
ATOM     28  O3'  DG A   2       9.982  20.208   3.008  1.00 20.77           O
ATOM     29  C2'  DG A   2       8.480  20.049   1.058  1.00 19.09           C
ATOM     30  C1'  DG A   2       8.752  20.320  -0.417  1.00 19.00           C
ATOM     31  N9   DG A   2       8.800  19.119  -1.251  1.00 18.41           N
ATOM     32  C8   DG A   2       9.891  18.580  -1.893  1.00 18.28           C
ATOM     33  N7   DG A   2       9.606  17.497  -2.563  1.00 18.16           N
ATOM     34  C5   DG A   2       8.251  17.298  -2.336  1.00 17.93           C
ATOM     35  C6   DG A   2       7.376  16.281  -2.810  1.00 17.95           C
ATOM     36  O6   DG A   2       7.641  15.301  -3.517  1.00 17.71           O
ATOM     37  N1   DG A   2       6.078  16.469  -2.349  1.00 18.47           N
ATOM     38  C2   DG A   2       5.670  17.498  -1.535  1.00 19.03           C
ATOM     39  N2   DG A   2       4.374  17.506  -1.197  1.00 19.60           N
ATOM     40  N3   DG A   2       6.473  18.453  -1.092  1.00 19.35           N
ATOM     41  C4   DG A   2       7.737  18.296  -1.537  1.00 17.37           C
ATOM     42  P    DT A   3      10.011  19.045   4.108  1.00 21.81           P
ATOM     43  OP1  DT A   3      10.235  19.681   5.425  1.00 24.71           O
ATOM     44  OP2  DT A   3      10.906  17.959   3.643  1.00 23.07           O
ATOM     45  O5'  DT A   3       8.532  18.470   4.070  1.00 21.48           O
ATOM     46  C5'  DT A   3       7.452  19.242   4.612  1.00 21.79           C
ATOM     47  C4'  DT A   3       6.153  18.502   4.416  1.00 21.68           C
ATOM     48  O4'  DT A   3       5.899  18.337   3.006  1.00 22.26           O
ATOM     49  C3'  DT A   3       6.130  17.086   4.976  1.00 22.13           C
ATOM     50  O3'  DT A   3       5.798  17.096   6.359  1.00 24.02           O
ATOM     51  C2'  DT A   3       5.039  16.439   4.154  1.00 22.71           C
ATOM     52  C1'  DT A   3       5.228  17.096   2.794  1.00 22.01           C
ATOM     53  N1   DT A   3       6.037  16.289   1.874  1.00 20.23           N
ATOM     54  C2   DT A   3       5.382  15.302   1.178  1.00 21.86           C
ATOM     55  O2   DT A   3       4.185  15.095   1.285  1.00 25.73           O
ATOM     56  N3   DT A   3       6.183  14.567   0.341  1.00 19.17           N
ATOM     57  C4   DT A   3       7.540  14.721   0.135  1.00 18.53           C
ATOM     58  O4   DT A   3       8.121  13.998  -0.668  1.00 20.53           O
ATOM     59  C5   DT A   3       8.169  15.755   0.922  1.00 18.77           C
ATOM     60  C7   DT A   3       9.642  15.978   0.784  1.00 18.48           C
ATOM     61  C6   DT A   3       7.398  16.472   1.749  1.00 19.98           C
"""

rna = """\
ATOM      1  O5'   G A   1       8.651  -8.743   4.346  1.00  0.00           O
ATOM      2  C5'   G A   1       7.735  -8.253   3.392  1.00  0.00           C
ATOM      3  C4'   G A   1       8.396  -7.192   2.507  1.00  0.00           C
ATOM      4  O4'   G A   1       7.511  -6.793   1.472  1.00  0.00           O
ATOM      5  C3'   G A   1       8.758  -5.911   3.261  1.00  0.00           C
ATOM      6  O3'   G A   1       9.987  -6.010   3.960  1.00  0.00           O
ATOM      7  C2'   G A   1       8.855  -4.944   2.089  1.00  0.00           C
ATOM      8  O2'   G A   1      10.066  -5.143   1.384  1.00  0.00           O
ATOM      9  C1'   G A   1       7.701  -5.411   1.204  1.00  0.00           C
ATOM     10  N9    G A   1       6.454  -4.653   1.460  1.00  0.00           N
ATOM     11  C8    G A   1       5.260  -5.082   1.985  1.00  0.00           C
ATOM     12  N7    G A   1       4.308  -4.191   1.942  1.00  0.00           N
ATOM     13  C5    G A   1       4.912  -3.083   1.357  1.00  0.00           C
ATOM     14  C6    G A   1       4.377  -1.795   1.047  1.00  0.00           C
ATOM     15  O6    G A   1       3.226  -1.394   1.204  1.00  0.00           O
ATOM     16  N1    G A   1       5.336  -0.936   0.510  1.00  0.00           N
ATOM     17  C2    G A   1       6.651  -1.291   0.265  1.00  0.00           C
ATOM     18  N2    G A   1       7.464  -0.354  -0.231  1.00  0.00           N
ATOM     19  N3    G A   1       7.144  -2.512   0.523  1.00  0.00           N
ATOM     20  C4    G A   1       6.230  -3.353   1.075  1.00  0.00           C
ATOM     21  P     A A   2      10.342  -5.036   5.198  1.00  0.00           P
ATOM     22  OP1   A A   2      11.732  -5.326   5.616  1.00  0.00           O
ATOM     23  OP2   A A   2       9.245  -5.144   6.185  1.00  0.00           O
ATOM     24  O5'   A A   2      10.301  -3.545   4.582  1.00  0.00           O
ATOM     25  C5'   A A   2      11.387  -2.992   3.859  1.00  0.00           C
ATOM     26  C4'   A A   2      11.046  -1.560   3.437  1.00  0.00           C
ATOM     27  O4'   A A   2       9.807  -1.514   2.753  1.00  0.00           O
ATOM     28  C3'   A A   2      10.900  -0.617   4.629  1.00  0.00           C
ATOM     29  O3'   A A   2      12.167  -0.168   5.083  1.00  0.00           O
ATOM     30  C2'   A A   2      10.048   0.493   4.017  1.00  0.00           C
ATOM     31  O2'   A A   2      10.862   1.442   3.353  1.00  0.00           O
ATOM     32  C1'   A A   2       9.217  -0.242   2.962  1.00  0.00           C
ATOM     33  N9    A A   2       7.794  -0.369   3.344  1.00  0.00           N
ATOM     34  C8    A A   2       7.080  -1.486   3.707  1.00  0.00           C
ATOM     35  N7    A A   2       5.798  -1.279   3.832  1.00  0.00           N
ATOM     36  C5    A A   2       5.646   0.067   3.520  1.00  0.00           C
ATOM     37  C6    A A   2       4.521   0.913   3.401  1.00  0.00           C
ATOM     38  N6    A A   2       3.267   0.490   3.594  1.00  0.00           N
ATOM     39  N1    A A   2       4.723   2.204   3.071  1.00  0.00           N
ATOM     40  C2    A A   2       5.965   2.631   2.857  1.00  0.00           C
ATOM     41  N3    A A   2       7.096   1.935   2.908  1.00  0.00           N
ATOM     42  C4    A A   2       6.863   0.637   3.246  1.00  0.00           C
"""

def residues(pdb_str, *resseqs):
  keep = [l for l in pdb_str.splitlines() if int(l[22:26]) in resseqs]
  return "CRYST1   90.000   90.000   90.000  90.00  90.00  90.00 P 1\n" + \
    "\n".join(keep) + "\nEND\n"

def ligand(resname, resseq, dx=25.0, extra=()):
  '''1rrr G A 1 as a HETATM residue resname resseq of chain A, dx A away;
  extra: (name, element, x, y, z) before the shift.'''
  out = []
  for l in rna.splitlines():
    if int(l[22:26]) != 1: continue
    out.append(l)
  for name, element, x, y, z in extra:
    out.append("ATOM     99 %-4s   G A   1    %8.3f%8.3f%8.3f  1.00  0.00"
               "          %2s" % (name, x, y, z, element))
  return "\n".join("HETATM" + l[6:17] + "%3s A%4d" % (resname, resseq) +
    l[26:30] + "%8.3f" % (float(l[30:38]) + dx) + l[38:] for l in out)

def renumber(pdb_str, offset, dx=25.0):
  '''The ATOM records of pdb_str with resseq + offset, x + dx.'''
  return "\n".join(l[:22] + "%4d" % (int(l[22:26]) + offset) + l[26:30] +
    "%8.3f" % (float(l[30:38]) + dx) + l[38:]
    for l in pdb_str.splitlines() if l.startswith("ATOM"))

def hetatm(pdb_str, resseq, resname=None, extra=()):
  '''pdb_str with residue resseq as HETATM, renamed to resname; extra:
  (name, element, x, y, z) added to it.'''
  out = []
  for l in pdb_str.splitlines():
    if l.startswith('ATOM') and int(l[22:26]) == resseq:
      l = 'HETATM' + l[6:]
      if resname is not None: l = l[:17] + '%3s' % resname + l[20:]
      last = l
    out.append(l)
  i = out.index(last) + 1
  for name, element, x, y, z in extra:
    out.insert(i, last[:12] + ' %-3s' % name + last[16:30] +
      '%8.3f%8.3f%8.3f' % (x, y, z) + last[54:76] + '%2s' % element)
    i += 1
  return '\n'.join(out) + '\n'

def after_ter(pdb_str, other):
  '''pdb_str (from residues()), TER, the records of other.'''
  return pdb_str.replace("END\n", "TER\n" + other + "\nEND\n")

def as_mmcif(pdb_str):
  '''pdb_str as mmCIF: each chain object of the hierarchy gets its own
  label_asym_id, all keep author chain A.'''
  inp = iotbx.pdb.input(lines=pdb_str.splitlines(), source_info=None)
  text = inp.construct_hierarchy().as_mmcif_string(
    crystal_symmetry=inp.crystal_symmetry())
  block = list(iotbx.cif.reader(input_string=text).model().values())[0]
  assert set(block['_atom_site.auth_asym_id']) == set(['A'])
  assert len(set(block['_atom_site.label_asym_id'])) > 1
  return text

def split(pdb_str, resseq, names):
  '''Atoms names of residue resseq in conformers A and B, B shifted 0.05 A.'''
  out = []
  for l in pdb_str.splitlines():
    if not (l.startswith('ATOM') and int(l[22:26]) == resseq and
            l[12:16].strip() in names):
      out.append(l)
      continue
    out.append(l[:16] + 'A' + l[17:54] + '  0.60' + l[60:])
    x = '%8.3f' % (float(l[30:38]) + 0.05)
    out.append(l[:16] + 'B' + l[17:30] + x + l[38:54] + '  0.40' + l[60:])
  return "\n".join(out) + "\n"

def run(prefix, pdb_str, *args, **kw):
  ext = kw.get('ext', 'pdb')
  fn = prefix + '.' + ext
  with open(fn, 'w') as f:
    f.write(pdb_str)
  fn_out, fn_txt = prefix + '_out.' + ext, prefix + '.txt'
  for p in (fn_out, fn_txt):
    if os.path.exists(p): os.remove(p)
  log = StringIO()
  run_program(program_class=reduce2.Program, logger=log,
    args=[fn, 'output.filename=' + fn_out, 'output.description_file_name='
          + fn_txt, 'overwrite=True'] + list(args))
  with open(fn_out) as f:
    out = f.read()
  with open(fn_txt) as f:
    txt = f.read()
  return out, txt, log.getvalue()

def hierarchy(pdb_str):
  return iotbx.pdb.input(lines=pdb_str.splitlines(),
    source_info=None).construct_hierarchy()

def on_oxygen(pdb_str):
  '''{(resseq, O name): sorted (element, altloc)} of H/D within 1.2 A of the
  O5' and O3' of each residue.'''
  root = hierarchy(pdb_str)
  atoms = list(root.atoms())
  result = {}
  for o in atoms:
    if o.name.strip() not in ("O5'", "O3'"): continue
    key = (o.parent().parent().resseq.strip(), o.name.strip())
    for h in atoms:
      if h.element.strip() in ('H', 'D') and h.distance(o) < 1.2:
        result.setdefault(key, set()).add((h.element.strip(),
          h.parent().altloc.strip()))
  return dict((k, sorted(v)) for k, v in result.items())

def ends_only(pdb_str, *expected):
  '''The terminal H are exactly on the (resseq, O) of expected, once each.'''
  found = on_oxygen(pdb_str)
  assert sorted(found) == sorted(expected), found
  for k in expected:
    assert found[k] == [('H', '')], found

def reported(text, *h_labels):
  assert CAVEAT in text, text
  lines = [l for l in text.splitlines()
           if l.strip().startswith(("HO5'", "HO3'"))]
  assert len(lines) == len(h_labels), lines
  for name, resseq in h_labels:
    assert [l for l in lines if name in l and (' %s' % resseq) in l], lines

def not_reported(text):
  assert CAVEAT not in text, text

def gap_reported(text, *o_labels):
  '''The free ends at chain breaks, (O name, resseq), once each.'''
  assert GAP_CAVEAT in text, text
  i = text.index(GAP_CAVEAT)
  assert GAP_CAVEAT not in text[i+1:], text
  # the section's title, not the lines that quote it
  head = [l for l in text.splitlines()
          if 'inside a chain' in l and not l.startswith(' ')]
  assert len(head) == 1 and head[0].startswith('%d ' % len(o_labels)), head
  lines = [l for l in text.splitlines()
           if l.strip().startswith(("O5'", "O3'"))]
  assert len(lines) == len(o_labels), lines
  for name, resseq in o_labels:
    assert [l for l in lines if name in l and (' %s' % resseq) in l], lines

def test_free_ends():
  for approach in ('add', 'place'):
    for name, pdb_str in (('dna', residues(dna, 1, 2)),
                          ('rna', residues(rna, 1, 2))):
      out, txt, log = run('free_%s_%s' % (name, approach), pdb_str,
        'approach=' + approach)
      ends_only(out, ('1', "O5'"), ('2', "O3'"))
      for text in (txt, log):
        reported(text, ("HO5'", 1), ("HO3'", 2))

def test_phosphorylated_end():
  # 5' phosphate at residue 2, O3' of 2 linked to P of 3
  for approach in ('add', 'place'):
    out, txt, log = run('phosphate_' + approach, residues(dna, 2, 3),
      'approach=' + approach)
    ends_only(out, ('3', "O3'"))
    for text in (txt, log):
      reported(text, ("HO3'", 3))

def test_chain_break():
  # residues 1 and 3: O3' of 1 is free but inside the chain
  pdb_str = residues(dna, 1, 3)
  for approach in ('add', 'place'):
    out, txt, log = run('gap_' + approach, pdb_str, 'approach=' + approach)
    ends_only(out, ('1', "O5'"), ('3', "O3'"))
    assert "DG A   1 HO3'" not in log, log
    for text in (txt, log):
      reported(text, ("HO5'", 1), ("HO3'", 3))
      gap_reported(text, ("O3'", 1))
  # every free end, as before
  out, txt, log = run('gap_qmready', pdb_str, 'qmready=True')
  ends_only(out, ('1', "O5'"), ('1', "O3'"), ('3', "O3'"))
  assert GAP_CAVEAT not in txt and GAP_CAVEAT not in log

def test_ter_at_gap():
  # a TER at the gap ends the chain object: the O3' of 1 is a chain end
  pdb_str = residues(dna, 1, 3)
  lines = pdb_str.splitlines()
  i = [k for k, l in enumerate(lines)
       if l.startswith('ATOM') and int(l[22:26]) == 3][0]
  pdb_str = "\n".join(lines[:i] + ['TER'] + lines[i:]) + "\n"
  out, txt, log = run('gap_ter', pdb_str)
  ends_only(out, ('1', "O5'"), ('1', "O3'"), ('3', "O3'"))
  reported(txt, ("HO5'", 1), ("HO3'", 1), ("HO3'", 3))
  assert GAP_CAVEAT not in txt

def test_nucleotide_ligand_after_polymer():
  # G A 101 after the polymer, same author chain: its own chain object (PDB:
  # TER; mmCIF: another label_asym_id), so DG 2 keeps its HO3' and G 101 gets
  # both of its own
  pdb_str = after_ter(residues(dna, 1, 2), ligand('G', 101))
  for ext, text in (('pdb', pdb_str), ('cif', as_mmcif(pdb_str))):
    out, txt, log = run('g_ligand_' + ext, text, ext=ext)
    ends_only(out, ('1', "O5'"), ('2', "O3'"), ('101', "O5'"), ('101', "O3'"))
    reported(txt, ("HO5'", 1), ("HO3'", 2), ("HO5'", 101), ("HO3'", 101))
    assert GAP_CAVEAT not in txt, txt

def test_modified_nucleotide_ligand():
  # OMG (modified_rna_dna) A 101 after the polymer: DG 2 keeps its HO3'; OMG
  # gets the H of its own dictionary (HO3', no HO5': a polymer residue)
  cm2 = (('CM2', 'C', 10.472, -3.922, 0.785),)
  pdb_str = after_ter(residues(dna, 1, 2), ligand('OMG', 101, extra=cm2))
  for ext, text in (('pdb', pdb_str), ('cif', as_mmcif(pdb_str))):
    out, txt, log = run('omg_ligand_' + ext, text, ext=ext)
    ends_only(out, ('1', "O5'"), ('2', "O3'"), ('101', "O3'"))
    reported(txt, ("HO5'", 1), ("HO3'", 2))
    assert GAP_CAVEAT not in txt, txt

def test_two_strands_one_chain_id():
  # two strands with author chain A, separated by TER (PDB) or label_asym_id
  # (mmCIF): each has its own 5' and 3' end
  pdb_str = after_ter(residues(dna, 1, 2), renumber(residues(dna, 1, 2), 10))
  for ext, text in (('pdb', pdb_str), ('cif', as_mmcif(pdb_str))):
    out, txt, log = run('two_strands_' + ext, text, ext=ext)
    ends_only(out, ('1', "O5'"), ('2', "O3'"), ('11', "O5'"), ('12', "O3'"))
    assert GAP_CAVEAT not in txt, txt

def test_free_het_nucleotides():
  # two free G after one TER: one chain object, but each G is all HETATM and
  # unlinked, so each is its own first and last
  pdb_str = after_ter(residues(dna, 1, 2),
    ligand('G', 101) + '\n' + ligand('G', 102, dx=40.0))
  assert len(list(hierarchy(pdb_str).chains())) == 2
  out, txt, log = run('two_free_g', pdb_str)
  ends_only(out, ('1', "O5'"), ('2', "O3'"), ('101', "O5'"), ('101', "O3'"),
            ('102', "O5'"), ('102', "O3'"))
  assert GAP_CAVEAT not in txt and GAP_CAVEAT not in log, txt

def test_het_nucleotide_without_ter():
  # a HETATM G right after the polymer, no TER: same chain object; DG 2
  # keeps its HO3'
  pdb_str = residues(dna, 1, 2).replace('END\n', ligand('G', 101) + '\nEND\n')
  assert len(list(hierarchy(pdb_str).chains())) == 1
  out, txt, log = run('g_no_ter', pdb_str)
  ends_only(out, ('1', "O5'"), ('2', "O3'"), ('101', "O5'"), ('101', "O3'"))
  assert GAP_CAVEAT not in txt, txt

def test_het_residues_inside_chain():
  # HETATM residue 2, linked on both sides, standard (DG) or modified (8OG):
  # not an end, so the chain's ends stay 1 and 3
  o8 = (('O8', 'O', 11.026, 19.029, -1.870),)
  for name, pdb_str in (
      ('dg', hetatm(residues(dna, 1, 2, 3), 2)),
      ('8og', hetatm(residues(dna, 1, 2, 3), 2, '8OG', o8))):
    out, txt, log = run('het_inside_' + name, pdb_str)
    ends_only(out, ('1', "O5'"), ('3', "O3'"))
    assert GAP_CAVEAT not in txt, txt
  # HETATM DG 2 linked to 1, then a gap (3 moved away): its free O3' is a
  # break inside the chain, not the end of a free heterogen
  moved = residues(dna, 1, 2).replace('END\n',
    renumber(residues(dna, 3), 0) + '\nEND\n')
  out, txt, log = run('het_inside_gap', hetatm(moved, 2))
  ends_only(out, ('1', "O5'"), ('3', "O3'"))
  gap_reported(txt, ("O3'", 2))

def test_same_chain_id_after_polymer():
  # a water of chain A after TER (another chain object) is not the last
  # nucleotide; DG 2 is
  water = ("TER\nHETATM   62  O   HOH A 101      20.000  20.000  20.000"
           "  1.00 20.00           O\nEND\n")
  pdb_str = residues(dna, 1, 2).replace("END\n", water)
  root = hierarchy(pdb_str)
  assert [c.id for c in root.chains()] == ['A', 'A']
  out, txt, log = run('water_after', pdb_str)
  ends_only(out, ('1', "O5'"), ('2', "O3'"))
  assert GAP_CAVEAT not in txt

def test_altlocs():
  out, txt, log = run('alt_5', split(residues(dna, 1, 2), 1, ("O5'", "C5'")))
  found = on_oxygen(out)
  assert found[('1', "O5'")] == [('H', 'A'), ('H', 'B')], found
  assert found[('2', "O3'")] == [('H', '')], found
  out, txt, log = run('alt_3', split(residues(dna, 1, 2), 2, ("O3'", "C3'")))
  found = on_oxygen(out)
  assert found[('1', "O5'")] == [('H', '')], found
  assert found[('2', "O3'")] == [('H', 'A'), ('H', 'B')], found

def test_keep_existing_h():
  pdb_str = residues(dna, 1, 2)
  out, txt, log = run('keep_none', pdb_str, 'keep_existing_H=True')
  ends_only(out, ('1', "O5'"), ('2', "O3'"))
  reported(txt, ("HO5'", 1), ("HO3'", 2))
  # all H supplied: nothing is new
  out2, txt, log = run('keep_all', out, 'keep_existing_H=True')
  ends_only(out2, ('1', "O5'"), ('2', "O3'"))
  not_reported(txt)
  not_reported(log)
  # HO5' left out: only that one is new
  partial = "\n".join(l for l in out.splitlines()
    if not (l.startswith('ATOM') and l[12:16].strip() == "HO5'"))
  out3, txt, log = run('keep_partial', partial, 'keep_existing_H=True')
  ends_only(out3, ('1', "O5'"), ('2', "O3'"))
  reported(txt, ("HO5'", 1))

def test_neutron_options():
  pdb_str = residues(dna, 1, 2)
  for option, kinds in (('all_h',    ['H']),
                        ('all_d',    ['D']),
                        ('hd_and_h', ['D', 'H']),
                        ('hd_and_d', ['D', 'H']),
                        ('all_hd',   ['D', 'H'])):
    out, txt, log = run('neutron_' + option, pdb_str,
      'neutron_option=' + option)
    found = on_oxygen(out)
    assert sorted(found) == [('1', "O5'"), ('2', "O3'")], (option, found)
    for k, v in found.items():
      assert sorted(e for e, altloc in v) == kinds, (option, k, v)
    # the report lists HO5'/HO3' also where the output has D or H/D pairs
    if option in ('all_h', 'all_d', 'hd_and_h'):
      for text in (txt, log):
        reported(text, ("HO5'", 1), ("HO3'", 2))
        assert NEUTRON_NOTE in ' '.join(text.split()), (option, text)
  # neutron distances without a conversion: the same sentence
  for args in (('use_neutron_distances=True',), ('qmready=True',)):
    out, txt, log = run('neutron_note_' + args[0].split('=')[0], pdb_str, *args)
    for text in (txt, log):
      reported(text, ("HO5'", 1), ("HO3'", 2))
      assert NEUTRON_NOTE in ' '.join(text.split()), (args, text)

def test_rerun_and_optimize():
  pdb_str = residues(rna, 1, 2)
  out, txt, log = run('first', pdb_str)
  ends_only(out, ('1', "O5'"), ('2', "O3'"))
  out2, txt, log = run('rerun', out)
  ends_only(out2, ('1', "O5'"), ('2', "O3'"))
  assert hierarchy(out).atoms_size() == hierarchy(out2).atoms_size()
  out3, txt, log = run('optimized', out, 'approach=optimize')
  ends_only(out3, ('1', "O5'"), ('2', "O3'"))
  assert hierarchy(out).atoms_size() == hierarchy(out3).atoms_size()

if (__name__ == "__main__"):
  test_free_ends()
  test_phosphorylated_end()
  test_chain_break()
  test_same_chain_id_after_polymer()
  test_free_het_nucleotides()
  test_het_nucleotide_without_ter()
  test_het_residues_inside_chain()
  test_ter_at_gap()
  test_nucleotide_ligand_after_polymer()
  test_modified_nucleotide_ligand()
  test_two_strands_one_chain_id()
  test_altlocs()
  test_keep_existing_h()
  test_neutron_options()
  test_rerun_and_optimize()
  print("OK")
