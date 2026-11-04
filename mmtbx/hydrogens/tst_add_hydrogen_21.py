"""
Negative controls for tst_add_hydrogen_20.py. exclude_H_on_esterified_O() counts
the heavy neighbours of an O by residue group and atom name: the same atom in
two conformers is one neighbour, another atom is a second one whatever its
altloc letters. So an O bonded to a second atom is an ester or an ether and keeps
no H, also if that atom has other letters than the H copies. The mirror of
tst_add_hydrogen_20 on RNA: a blank 3' O on C3' A/B keeps its H.
Fragments: 1rrr A 2-3 with A 3 in A/B (A 2 with its base in A/B), 4d3w NGA B 1
with C6, O6 in A/B and GAL B 2 in A/B. A 2 is called AMP where O3' is bonded to
P: the polymer dictionary and the 3' end rule never offer an HO3' that bonds to
P, the AMP dictionary does, and only exclude_H_on_esterified_O() takes it away.
"""
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out
from mmtbx.hydrogens import reduce_hydrogen

pdb_rna = """\
CRYST1   80.000   80.000   80.000  90.00  90.00  90.00 P 1
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
ATOM     33  N9 A  A A   2       7.794  -0.369   3.344  0.60  0.00           N
ATOM     33  N9 B  A A   2       7.894  -0.369   3.344  0.40  0.00           N
ATOM     34  C8 A  A A   2       7.080  -1.486   3.707  0.60  0.00           C
ATOM     34  C8 B  A A   2       7.180  -1.486   3.707  0.40  0.00           C
ATOM     35  N7 A  A A   2       5.798  -1.279   3.832  0.60  0.00           N
ATOM     35  N7 B  A A   2       5.898  -1.279   3.832  0.40  0.00           N
ATOM     36  C5 A  A A   2       5.646   0.067   3.520  0.60  0.00           C
ATOM     36  C5 B  A A   2       5.746   0.067   3.520  0.40  0.00           C
ATOM     37  C6 A  A A   2       4.521   0.913   3.401  0.60  0.00           C
ATOM     37  C6 B  A A   2       4.621   0.913   3.401  0.40  0.00           C
ATOM     38  N6 A  A A   2       3.267   0.490   3.594  0.60  0.00           N
ATOM     38  N6 B  A A   2       3.367   0.490   3.594  0.40  0.00           N
ATOM     39  N1 A  A A   2       4.723   2.204   3.071  0.60  0.00           N
ATOM     39  N1 B  A A   2       4.823   2.204   3.071  0.40  0.00           N
ATOM     40  C2 A  A A   2       5.965   2.631   2.857  0.60  0.00           C
ATOM     40  C2 B  A A   2       6.065   2.631   2.857  0.40  0.00           C
ATOM     41  N3 A  A A   2       7.096   1.935   2.908  0.60  0.00           N
ATOM     41  N3 B  A A   2       7.196   1.935   2.908  0.40  0.00           N
ATOM     42  C4 A  A A   2       6.863   0.637   3.246  0.60  0.00           C
ATOM     42  C4 B  A A   2       6.963   0.637   3.246  0.40  0.00           C
ATOM     43  P  A  A A   3      12.341   0.647   6.465  0.60  0.00           P
ATOM     43  P  B  A A   3      12.441   0.647   6.465  0.40  0.00           P
ATOM     44  OP1A  A A   3      13.769   1.006   6.611  0.60  0.00           O
ATOM     44  OP1B  A A   3      13.869   1.006   6.611  0.40  0.00           O
ATOM     45  OP2A  A A   3      11.662  -0.117   7.534  0.60  0.00           O
ATOM     45  OP2B  A A   3      11.762  -0.117   7.534  0.40  0.00           O
ATOM     46  O5'A  A A   3      11.514   2.005   6.194  0.60  0.00           O
ATOM     46  O5'B  A A   3      11.614   2.005   6.194  0.40  0.00           O
ATOM     47  C5'A  A A   3      10.977   2.765   7.257  0.60  0.00           C
ATOM     47  C5'B  A A   3      11.077   2.765   7.257  0.40  0.00           C
ATOM     48  C4'A  A A   3      10.049   3.851   6.701  0.60  0.00           C
ATOM     48  C4'B  A A   3      10.149   3.851   6.701  0.40  0.00           C
ATOM     49  O4'A  A A   3       9.026   3.266   5.913  0.60  0.00           O
ATOM     49  O4'B  A A   3       9.126   3.266   5.913  0.40  0.00           O
ATOM     50  C3'A  A A   3       9.295   4.550   7.828  0.60  0.00           C
ATOM     50  C3'B  A A   3       9.395   4.550   7.828  0.40  0.00           C
ATOM     51  O3'   A A   3      10.040   5.574   8.463  1.00  0.00           O
ATOM     52  C2'A  A A   3       8.080   5.100   7.082  0.60  0.00           C
ATOM     52  C2'B  A A   3       8.180   5.100   7.082  0.40  0.00           C
ATOM     53  O2'A  A A   3       8.371   6.366   6.525  0.60  0.00           O
ATOM     53  O2'B  A A   3       8.471   6.366   6.525  0.40  0.00           O
ATOM     54  C1'A  A A   3       7.886   4.109   5.931  0.60  0.00           C
ATOM     54  C1'B  A A   3       7.986   4.109   5.931  0.40  0.00           C
ATOM     55  N9 A  A A   3       6.637   3.330   6.084  0.60  0.00           N
ATOM     55  N9 B  A A   3       6.737   3.330   6.084  0.40  0.00           N
ATOM     56  C8 A  A A   3       6.473   1.983   6.295  0.60  0.00           C
ATOM     56  C8 B  A A   3       6.573   1.983   6.295  0.40  0.00           C
ATOM     57  N7 A  A A   3       5.227   1.602   6.366  0.60  0.00           N
ATOM     57  N7 B  A A   3       5.327   1.602   6.366  0.40  0.00           N
ATOM     58  C5 A  A A   3       4.509   2.773   6.143  0.60  0.00           C
ATOM     58  C5 B  A A   3       4.609   2.773   6.143  0.40  0.00           C
ATOM     59  C6 A  A A   3       3.130   3.053   6.015  0.60  0.00           C
ATOM     59  C6 B  A A   3       3.230   3.053   6.015  0.40  0.00           C
ATOM     60  N6 A  A A   3       2.183   2.119   6.151  0.60  0.00           N
ATOM     60  N6 B  A A   3       2.283   2.119   6.151  0.40  0.00           N
ATOM     61  N1 A  A A   3       2.754   4.315   5.725  0.60  0.00           N
ATOM     61  N1 B  A A   3       2.854   4.315   5.725  0.40  0.00           N
ATOM     62  C2 A  A A   3       3.691   5.248   5.573  0.60  0.00           C
ATOM     62  C2 B  A A   3       3.791   5.248   5.573  0.40  0.00           C
ATOM     63  N3 A  A A   3       5.012   5.114   5.663  0.60  0.00           N
ATOM     63  N3 B  A A   3       5.112   5.114   5.663  0.40  0.00           N
ATOM     64  C4 A  A A   3       5.362   3.831   5.950  0.60  0.00           C
ATOM     64  C4 B  A A   3       5.462   3.831   5.950  0.40  0.00           C
END
"""

pdb_glycan = """\
CRYST1   67.843   26.790   30.439  90.00  90.00  90.00 P 1
HETATM  110  C1  NGA B   1      52.184  16.633  15.213  1.00 10.82           C
HETATM  111  C2  NGA B   1      51.328  15.445  15.648  1.00 10.31           C
HETATM  112  C3  NGA B   1      50.232  15.929  16.563  1.00  9.54           C
HETATM  113  C4  NGA B   1      50.840  16.728  17.712  1.00 10.02           C
HETATM  114  C5  NGA B   1      51.754  17.850  17.194  1.00 10.27           C
HETATM  115  C6 ANGA B   1      52.419  18.596  18.349  0.60 10.78           C
HETATM  115  C6 BNGA B   1      52.519  18.596  18.349  0.40 10.78           C
HETATM  116  C7  NGA B   1      51.064  13.601  14.053  1.00 10.97           C
HETATM  117  C8  NGA B   1      50.412  13.261  12.752  1.00 11.18           C
HETATM  118  N2  NGA B   1      50.777  14.837  14.457  1.00 10.28           N
HETATM  119  O1  NGA B   1      53.223  16.116  14.399  1.00 13.50           O
HETATM  120  O3  NGA B   1      49.470  14.809  17.068  1.00  7.48           O
HETATM  121  O4  NGA B   1      51.582  15.870  18.572  1.00 10.21           O
HETATM  122  O5  NGA B   1      52.747  17.304  16.339  1.00 10.70           O
HETATM  123  O6 ANGA B   1      52.951  19.847  17.881  0.60 11.05           O
HETATM  123  O6 BNGA B   1      53.051  19.847  17.881  0.40 11.05           O
HETATM  124  O7  NGA B   1      51.809  12.830  14.679  1.00 12.42           O
HETATM  140  C1 AGAL B   2      48.317  14.425  16.335  0.60  7.03           C
HETATM  141  C2 AGAL B   2      47.709  13.209  17.044  0.60  6.63           C
HETATM  142  C3 AGAL B   2      46.407  12.813  16.374  0.60  6.30           C
HETATM  143  C4 AGAL B   2      45.528  14.023  16.145  0.60  6.50           C
HETATM  144  C5 AGAL B   2      46.290  15.197  15.510  0.60  6.90           C
HETATM  145  C6 AGAL B   2      45.452  16.426  15.354  0.60  7.16           C
HETATM  146  O2 AGAL B   2      48.594  12.104  17.056  0.60  6.77           O
HETATM  147  O3 AGAL B   2      45.706  11.865  17.180  0.60  6.72           O
HETATM  148  O4 AGAL B   2      44.917  14.359  17.391  0.60  6.48           O
HETATM  149  O5 AGAL B   2      47.392  15.490  16.361  0.60  6.90           O
HETATM  150  O6 AGAL B   2      46.180  17.431  14.607  0.60  8.08           O
HETATM  140  C1 BGAL B   2      48.467  14.425  16.335  0.40  7.03           C
HETATM  141  C2 BGAL B   2      47.859  13.209  17.044  0.40  6.63           C
HETATM  142  C3 BGAL B   2      46.557  12.813  16.374  0.40  6.30           C
HETATM  143  C4 BGAL B   2      45.678  14.023  16.145  0.40  6.50           C
HETATM  144  C5 BGAL B   2      46.440  15.197  15.510  0.40  6.90           C
HETATM  145  C6 BGAL B   2      45.602  16.426  15.354  0.40  7.16           C
HETATM  146  O2 BGAL B   2      48.744  12.104  17.056  0.40  6.77           O
HETATM  147  O3 BGAL B   2      45.856  11.865  17.180  0.40  6.72           O
HETATM  148  O4 BGAL B   2      45.067  14.359  17.391  0.40  6.48           O
HETATM  149  O5 BGAL B   2      47.542  15.490  16.361  0.40  6.90           O
HETATM  150  O6 BGAL B   2      46.330  17.431  14.607  0.40  8.08           O
END
"""

def blank_residue(pdb_str, resseq):
  '''Residue resseq in conformer A alone, written blank at full occupancy.'''
  result = []
  for l in pdb_str.splitlines():
    if l[:6] in ('ATOM  ', 'HETATM') and l[22:26].strip() == resseq:
      if l[16] == 'B': continue
      if l[16] == 'A': l = l[:16] + ' ' + l[17:54] + '  1.00' + l[60:]
    result.append(l)
  return '\n'.join(result) + '\n'

def as_amp(pdb_str):
  '''Residue 2 of pdb_rna as AMP.'''
  return ''.join((l[:17] + 'AMP' + l[20:]
    if l[:4] == 'ATOM' and l[22:26].strip() == '2' else l) + '\n'
    for l in pdb_str.splitlines())

def other_letters(pdb_str, resseq_o, name_c, resseq_p):
  '''
  Own atom name_c in A and B (B 0.1 A off), residue resseq_p from A/B to C/D,
  as if the partner of the blank O had been built apart from its own residue.
  '''
  result = []
  for l in pdb_str.splitlines():
    if l[:6] in ('ATOM  ', 'HETATM'):
      if l[22:26].strip() == resseq_p and l[16] in 'AB':
        l = l[:16] + {'A':'C', 'B':'D'}[l[16]] + l[17:]
      if l[22:26].strip() == resseq_o and l[12:16].strip() == name_c:
        result.append(l[:16] + 'A' + l[17:54] + '  0.60' + l[60:])
        x = float(l[30:38]) + 0.1
        l = l[:16] + 'B' + l[17:30] + '%8.3f' % x + l[38:54] + '  0.40' + l[60:]
    result.append(l)
  return '\n'.join(result) + '\n'

def placed(pdb_str, **kw):
  '''H as (resseq, name, altloc) after placement.'''
  model = mmtbx.model.manager(
    model_input = iotbx.pdb.input(source_info=None, lines=pdb_str.splitlines()),
    stop_for_unknowns = False, log = null_out())
  obj = reduce_hydrogen.place_hydrogens(model = model, **kw)
  obj.run()
  return [
    (a.parent().parent().resseq.strip(), a.name.strip(), a.parent().altloc)
    for a in obj.get_model().get_hierarchy().atoms()
    if a.element.strip() == 'H']

def test_nucleotide_o3_with_split_p():
  # O3' of AMP 2: C3' and P of A 3 (A and B), two atoms. HO2' keeps its H.
  for pdb_str in (pdb_rna, blank_residue(pdb_rna, '2')):
    hs = placed(as_amp(pdb_str))
    assert [h for h in hs if h[:2] == ('2', "HO3'")] == [], hs
    assert [h for h in hs if h[:2] == ('2', "HO2'")] != [], hs

def test_glycosidic_o_with_split_c1():
  # O3 of NGA 1: C3 and C1 of GAL 2 (A and B), two atoms, C6 split or not.
  # exclude_H_on_links() removes this H too, so this guards the result only.
  for pdb_str in (pdb_glycan, blank_residue(pdb_glycan, '1')):
    hs = placed(pdb_str)
    assert [h for h in hs if h[:2] == ('1', 'HO3')] == [], hs
    assert [h for h in hs if h[:2] == ('1', 'HO4')] != [], hs

def test_partner_with_other_letters():
  # O bonded to C A/B of its own residue and to the next residue in C/D; the H
  # copies are in A and B. Not per altloc: H A sees C A alone and H would stay.
  # Here exclude_H_on_links() leaves the H: a link H is judged per conformer.
  for pdb_str in (pdb_rna, blank_residue(pdb_rna, '2')):
    hs = placed(as_amp(other_letters(pdb_str, '2', "C3'", '3')))
    assert [h for h in hs if h[:2] == ('2', "HO3'")] == [], hs
    assert [h for h in hs if h[:2] == ('2', "HO2'")] != [], hs
  for pdb_str in (pdb_glycan, blank_residue(pdb_glycan, '1')):
    hs = placed(other_letters(pdb_str, '1', 'C3', '2'))
    assert [h for h in hs if h[:2] == ('1', 'HO3')] == [], hs
    assert [h for h in hs if h[:2] == ('1', 'HO4')] != [], hs

def test_rna_terminal_o3_on_split_c3():
  # C3' A and B are one neighbour, as CZ A and B in tst_add_hydrogen_20
  hs = placed(pdb_rna, na_terminal_h = True)
  assert sorted(h[2] for h in hs if h[:2] == ('3', "HO3'")) == ['A', 'B'], hs

if (__name__ == "__main__"):
  test_nucleotide_o3_with_split_p()
  test_glycosidic_o_with_split_c1()
  test_partner_with_other_letters()
  test_rna_terminal_o3_on_split_c3()
  print("OK")
