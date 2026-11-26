"""reduce2 placement: the free-form H on N of a peptide-like residue (A56).
All but one dictionary H on N are left out only where a bond restraint joins
that N to a C of another residue (a peptide bond), not by proximity or by
the residue being peptide-like. Before, every peptide-like residue lost them:
SAM and SAH (GeoStd neutral NH2, H-N-H 106.6 / 107.1) kept one N-H, and so did
DNE alone. Unchanged: a peptide-bonded N (HT7 after DNE keeps HN1 only), a
chain start (n_terminal_charge decides: H1, H2, H3), and the H on OXT, never
added. A metal-bound chain start handled in one conformer leaves the free-form
H in the others. Fixtures: CCD ideal heavy atoms of SAM and SAH as residue 406 (not a
chain start for n_terminal_charge=residue_one); DNE and HT7 of
tst_add_hydrogen_23.
"""
from __future__ import absolute_import, division, print_function
import math
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out
from scitbx import matrix
from mmtbx.hydrogens import reduce_hydrogen
from mmtbx.hydrogens import tst_add_hydrogen_23 as t23
from mmtbx.regression import tst_monomer_library_terminal_nh3 as tml

pdb_sam = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
HETATM    1  N   SAM A 406       7.950  14.496  15.920  1.00 20.00           N
HETATM    2  CA  SAM A 406       8.588  15.790  15.642  1.00 20.00           C
HETATM    3  C   SAM A 406       7.601  16.699  14.955  1.00 20.00           C
HETATM    4  O   SAM A 406       6.730  16.227  14.244  1.00 20.00           O
HETATM    5  OXT SAM A 406       7.674  17.906  15.111  1.00 20.00           O
HETATM    6  CB  SAM A 406       9.801  15.576  14.735  1.00 20.00           C
HETATM    7  CG  SAM A 406      10.857  14.756  15.480  1.00 20.00           C
HETATM    8  SD  SAM A 406      12.296  14.503  14.405  1.00 20.00           S
HETATM    9  CE  SAM A 406      11.662  13.181  13.337  1.00 20.00           C
HETATM   10  C5' SAM A 406      13.352  13.532  15.515  1.00 20.00           C
HETATM   11  C4' SAM A 406      14.662  13.191  14.804  1.00 20.00           C
HETATM   12  O4' SAM A 406      15.421  14.390  14.574  1.00 20.00           O
HETATM   13  C3' SAM A 406      15.517  12.268  15.694  1.00 20.00           C
HETATM   14  O3' SAM A 406      15.769  11.023  15.039  1.00 20.00           O
HETATM   15  C2' SAM A 406      16.840  13.055  15.891  1.00 20.00           C
HETATM   16  O2' SAM A 406      17.968  12.181  15.813  1.00 20.00           O
HETATM   17  C1' SAM A 406      16.815  14.036  14.691  1.00 20.00           C
HETATM   18  N9  SAM A 406      17.627  15.222  14.975  1.00 20.00           N
HETATM   19  C8  SAM A 406      17.191  16.380  15.548  1.00 20.00           C
HETATM   20  N7  SAM A 406      18.176  17.224  15.651  1.00 20.00           N
HETATM   21  C5  SAM A 406      19.304  16.663  15.153  1.00 20.00           C
HETATM   22  C6  SAM A 406      20.635  17.087  15.003  1.00 20.00           C
HETATM   23  N6  SAM A 406      21.028  18.347  15.419  1.00 20.00           N
HETATM   24  N1  SAM A 406      21.509  16.249  14.454  1.00 20.00           N
HETATM   25  C2  SAM A 406      21.143  15.044  14.056  1.00 20.00           C
HETATM   26  N3  SAM A 406      19.908  14.605  14.176  1.00 20.00           N
HETATM   27  C4  SAM A 406      18.966  15.370  14.718  1.00 20.00           C
END
"""

pdb_sah = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
HETATM    1  N   SAH A 406      15.726  15.720  22.094  1.00 20.00           N
HETATM    2  CA  SAH A 406      14.395  15.536  21.502  1.00 20.00           C
HETATM    3  CB  SAH A 406      14.510  14.662  20.251  1.00 20.00           C
HETATM    4  CG  SAH A 406      15.428  15.342  19.236  1.00 20.00           C
HETATM    5  SD  SAH A 406      15.564  14.309  17.752  1.00 20.00           S
HETATM    6  C   SAH A 406      13.490  14.866  22.502  1.00 20.00           C
HETATM    7  O   SAH A 406      13.951  14.105  23.321  1.00 20.00           O
HETATM    8  OXT SAH A 406      12.170  15.110  22.483  1.00 20.00           O
HETATM    9  C5' SAH A 406      16.679  15.332  16.754  1.00 20.00           C
HETATM   10  C4' SAH A 406      16.946  14.640  15.415  1.00 20.00           C
HETATM   11  O4' SAH A 406      15.715  14.492  14.690  1.00 20.00           O
HETATM   12  C3' SAH A 406      17.887  15.507  14.557  1.00 20.00           C
HETATM   13  O3' SAH A 406      19.089  14.799  14.250  1.00 20.00           O
HETATM   14  C2' SAH A 406      17.072  15.777  13.265  1.00 20.00           C
HETATM   15  O2' SAH A 406      17.907  15.688  12.108  1.00 20.00           O
HETATM   16  C1' SAH A 406      16.033  14.629  13.288  1.00 20.00           C
HETATM   17  N9  SAH A 406      14.840  14.992  12.521  1.00 20.00           N
HETATM   18  C8  SAH A 406      13.766  15.693  12.984  1.00 20.00           C
HETATM   19  N7  SAH A 406      12.887  15.837  12.036  1.00 20.00           N
HETATM   20  C5  SAH A 406      13.339  15.239  10.907  1.00 20.00           C
HETATM   21  C6  SAH A 406      12.835  15.075   9.606  1.00 20.00           C
HETATM   22  N6  SAH A 406      11.601  15.595   9.254  1.00 20.00           N
HETATM   23  N1  SAH A 406      13.570  14.412   8.720  1.00 20.00           N
HETATM   24  C2  SAH A 406      14.749  13.912   9.047  1.00 20.00           C
HETATM   25  N3  SAH A 406      15.259  14.041  10.253  1.00 20.00           N
HETATM   26  C4  SAH A 406      14.600  14.694  11.203  1.00 20.00           C
END
"""

def place(pdb_str, **kw):
  inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  model = mmtbx.model.manager(model_input=inp,
    restraint_objects=tml.old_dictionaries(pdb_str), log=null_out())
  o = reduce_hydrogen.place_hydrogens(model=model, stop_for_unknowns=False,
    **kw)
  o.run()
  return o

def n_h(o, resname, n_name='N'):
  """H/D names on the N of resname in the final model, and the N and its
  heavy neighbour's sites"""
  h = o.get_model().get_hierarchy()
  for ag in h.atom_groups():
    n = ag.get_atom(n_name)
    if ag.resname.strip() != resname or n is None: continue
    hs = [a for a in ag.atoms() if a.element.strip().upper() in ('H', 'D')
          and a.distance(n) < 1.3]
    return n, sorted(hs, key=lambda a: a.name)
  return None, []

def angle(a, b, c):
  u, v = matrix.col(a) - matrix.col(b), matrix.col(c) - matrix.col(b)
  return math.degrees(math.acos(max(-1, min(1, u.dot(v) / (u.length() *
    v.length())))))

def test_sam_sah():
  # no peptide bond to N: both dictionary H, pyramidal at the dictionary
  # H-N-H; nothing counted by the rule
  for pdb_str, resname, hnh in ((pdb_sam, 'SAM', 106.6),
                                (pdb_sah, 'SAH', 107.1)):
    o = place(pdb_str)
    n, hs = n_h(o, resname)
    assert [a.name.strip() for a in hs] == ['HN1', 'HN2'], (resname,
      [a.name for a in hs])
    value = angle(hs[0].xyz, n.xyz, hs[1].xyz)
    assert abs(value - hnh) < 0.3, (resname, value)
    assert o.not_placed_rules['free_form'] == 0, o.not_placed_rules

def test_dne_alone():
  # DNE G 2 without HT7: no peptide bond, H and H2
  dne = '\n'.join(l for l in t23.pdb_str.splitlines() if 'HT7' not in l)
  o = place(dne + '\n')
  n, hs = n_h(o, 'DNE')
  assert [a.name.strip() for a in hs] == ['H', 'H2'], [a.name for a in hs]
  assert o.not_placed_rules['free_form'] == 0, o.not_placed_rules

@t23.old_ht7_library
def test_peptide_bond_and_chain_start():
  # HT7 G 3 after DNE G 2: a bond restraint N-C to DNE, HA2 left out and
  # counted; DNE alone as residue 1 is a chain start without a next residue:
  # ambiguous, its dictionary H and H2 kept (A57)
  o = place(t23.pdb_str)
  n, hs = n_h(o, 'HT7')
  assert [a.name.strip() for a in hs] == ['HN1'], [a.name for a in hs]
  assert o.not_placed_rules['free_form'] == 1, o.not_placed_rules
  start = t23.pdb_str.replace('DNE G   2', 'DNE G   1')
  start = '\n'.join(l for l in start.splitlines() if 'HT7' not in l) + '\n'
  o = place(start)
  n, hs = n_h(o, 'DNE')
  assert [a.name.strip() for a in hs] == ['H', 'H2'], [a.name for a in hs]
  assert o.not_placed_rules['free_form'] == 0, o.not_placed_rules

def test_sam_chain_start():
  # SAM as residue 1 is a chain start without a peptide bond: standalone, so
  # chain-start placement leaves it and the free-form rule keeps both
  # dictionary H, HN1 and HN2 (A57)
  o = place(pdb_sam.replace('SAM A 406', 'SAM A   1'))
  n, hs = n_h(o, 'SAM')
  assert [a.name.strip() for a in hs] == ['HN1', 'HN2'], [a.name for a in hs]
  assert o.not_placed_rules['free_form'] == 0, o.not_placed_rules

def dne_split_cu(ht7=False):
  """DNE G 2 alone (with ht7: peptide-bonded to HT7 G 3, a terminus), every
  atom in conformers A and B (occupancy 0.5), and a Cu in conformer A only,
  2.00 A from N (CB-CA-N-Cu 180, CA-N-Cu 109.5)"""
  lines = [l for l in t23.pdb_str.splitlines() if ht7 or 'HT7' not in l]
  out, sites = [], {}
  for l in lines:
    if not l.startswith('HETATM'):
      if not l.startswith('END'): out.append(l)
      continue
    if l[17:20] == 'DNE':
      sites[l[12:16].strip()] = matrix.col((float(l[30:38]),
        float(l[38:46]), float(l[46:54])))
    for alt in 'AB':
      out.append(l[:16] + alt + l[17:54] + '  0.50' + l[60:])
  a, b, c = sites['CB'], sites['CA'], sites['N']
  bc = (c - b).normalize(); nn = (b - a).cross(bc).normalize(); mm = nn.cross(bc)
  th, ph = math.radians(109.5), math.radians(180)
  cu = c + 2.00 * (-math.cos(th) * bc + math.sin(th) * (math.cos(ph) * mm +
    math.sin(ph) * nn))
  out.append("HETATM14170 CU  ACU  G   4    %8.3f%8.3f%8.3f  0.50 20.00"
    "          CU" % tuple(cu))
  return '\n'.join(out) + '\nEND\n'

def split_ab(lines):
  """every HETATM line in conformers A and B at occupancy 0.5"""
  out = []
  for l in lines:
    if not l.startswith('HETATM'):
      if not l.startswith('END'): out.append(l)
      continue
    for alt in 'AB':
      out.append(l[:16] + alt + l[17:54] + '  0.50' + l[60:])
  return out

def cu_on_n(pdb_str, resname):
  """Cu 2.00 A from N of resname (CB-CA-N-Cu 180, CA-N-Cu 109.5)"""
  s = {}
  for l in pdb_str.splitlines():
    if l.startswith('HETATM') and l[17:20] == resname:
      s[l[12:16].strip()] = matrix.col((float(l[30:38]), float(l[38:46]),
        float(l[46:54])))
  a, b, c = s['CB'], s['CA'], s['N']
  bc = (c - b).normalize(); nn = (b - a).cross(bc).normalize()
  mm = nn.cross(bc)
  th, ph = math.radians(109.5), math.radians(180)
  return c + 2.00 * (-math.cos(th) * bc + math.sin(th) * (math.cos(ph) * mm +
    math.sin(ph) * nn))

def h_on_n(o, resname):
  """{altloc: H names on N} per atom group of resname, N looked up in the
  residue group"""
  found = {}
  for ag in o.get_model().get_hierarchy().atom_groups():
    if ag.resname.strip() != resname: continue
    n = [x.get_atom('N') for x in ag.parent().atom_groups()
         if x.get_atom('N') is not None and x.altloc in ('', ag.altloc)][0]
    found[ag.altloc] = sorted(a.name.strip() for a in ag.atoms()
      if a.element.strip().upper() in ('H', 'D') and a.distance(n) < 1.3)
  return found

def test_metal_chain_start_one_conformer():
  # DNE peptide-bonded to HT7 (a terminus, A57), Cu in conformer A only:
  # chain-start placement makes N/A an NH3, of which the metal-bound amine
  # conversion keeps H2, H3; N/B, not handled and without a peptide bond to
  # its N, keeps both dictionary H (H, H2). Before, the whole residue counted
  # as handled and N/B kept H only
  o = place(dne_split_cu(ht7=True), metal_bound_amines_as_nh2=True)
  h = o.get_model().get_hierarchy()
  found = {}
  for ag in h.atom_groups():
    if ag.resname.strip() != 'DNE': continue
    n = ag.get_atom('N')
    found[ag.altloc] = sorted(a.name.strip() for a in ag.atoms()
      if a.element.strip().upper() in ('H', 'D') and a.distance(n) < 1.3)
  assert found['B'] == ['H', 'H2'], found
  assert found['A'] == ['H2', 'H3'], found
  # the rule counts HT7's peptide-bonded N, in A and B, not DNE's
  assert o.not_placed_rules['free_form'] == 2, o.not_placed_rules
  # one conformer (blank) with the Cu: the whole residue is handled, H2, H3
  # as before (A57 not worse)
  blank = []
  for l in dne_split_cu(ht7=True).splitlines():
    if l.startswith('HETATM') and l[16] == 'B': continue
    if l.startswith('HETATM'): l = l[:16] + ' ' + l[17:54] + '  1.00' + l[60:]
    blank.append(l)
  o = place('\n'.join(blank) + '\n', metal_bound_amines_as_nh2=True)
  n, hs = n_h(o, 'DNE')
  assert [a.name.strip() for a in hs] == ['H2', 'H3'], [a.name for a in hs]
  assert o.not_placed_rules['free_form'] == 1, o.not_placed_rules # HT7
  # SAM 406 in A and B, Cu on N/A: standalone, so the metal option does not
  # convert it either (A57): both conformers keep HN1, HN2
  cu = cu_on_n(pdb_sam, 'SAM')
  sam = split_ab(pdb_sam.splitlines()) + ["HETATM   28 CU  ACU  A 407    "
    "%8.3f%8.3f%8.3f  0.50 20.00          CU" % tuple(cu)]
  o = place('\n'.join(sam) + '\nEND\n', metal_bound_amines_as_nh2=True)
  found = h_on_n(o, 'SAM')
  assert found == {'A': ['HN1', 'HN2'], 'B': ['HN1', 'HN2']}, found
  # supplied chain start in one conformer (keep_existing_H): DNE (bonded to
  # HT7) with N, CA, C, O blank and the side chain in A and B, Cu blank; the
  # supplied H2 on N in A only stays, conformer B gets H, H2
  lines = []
  for l in dne_split_cu(ht7=True).splitlines():
    if l.startswith('HETATM') and (l[12:16].strip() in ('N', 'CA', 'C', 'O')
        or ' CU ' in l):
      if l[16] == 'B': continue
      l = l[:16] + ' ' + l[17:54] + '  1.00' + l[60:]
    lines.append(l)
  o = place('\n'.join(lines) + '\n', metal_bound_amines_as_nh2=True)
  h = o.get_model().get_hierarchy()
  for rg in h.residue_groups():
    if rg.atom_groups()[0].resname.strip() != 'DNE': continue
    blank = [ag for ag in rg.atom_groups() if ag.altloc == ''][0]
    ag_a = [ag for ag in rg.atom_groups() if ag.altloc == 'A'][0]
    n = blank.get_atom('N')
    on_n = [(ag, a) for ag in rg.atom_groups() for a in ag.atoms()
            if a.element.strip() == 'H' and a.distance(n) < 1.3]
    for ag, a in on_n:
      ag.remove_atom(a)
      if ag.altloc == '' and a.name.strip() == 'H2': ag_a.append_atom(a)
  o = place(t23.pdb_str.splitlines()[0] + '\n' + h.as_pdb_string(),
    metal_bound_amines_as_nh2=True, keep_existing_H=True)
  found = h_on_n(o, 'DNE')
  assert found['A'] == ['H2'] and found['B'] == ['H', 'H2'], found
  # the same fixture without Cu: both conformers H, H2
  no_cu = '\n'.join(l for l in dne_split_cu().splitlines() if ' CU ' not in l)
  o = place(no_cu + '\n', metal_bound_amines_as_nh2=True)
  h = o.get_model().get_hierarchy()
  for ag in h.atom_groups():
    n = ag.get_atom('N')
    hs = sorted(a.name.strip() for a in ag.atoms()
      if a.element.strip().upper() in ('H', 'D') and a.distance(n) < 1.3)
    assert hs == ['H', 'H2'], (ag.altloc, hs)

if __name__ == '__main__':
  test_sam_sah()
  test_dne_alone()
  test_peptide_bond_and_chain_start()
  test_sam_chain_start()
  test_metal_chain_start_one_conformer()
  print('OK')
