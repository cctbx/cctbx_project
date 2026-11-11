"""
place_hydrogens records every H it adds and then removes (not_placed): model
id, chain, residue with insertion code, conformer, H, parent atom, cause and
detail; residues without a dictionary are not assessed; waters and free-form
H are rules with a count. reduce2 writes the full list, grouped by cause and
residue, into the description file and counts per cause into the log.
Fixtures from the tests of each cause.
"""
from __future__ import absolute_import, division, print_function
import collections
import math
import os
import iotbx.cif
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out
from scitbx import matrix
from six.moves import cStringIO as StringIO
from iotbx.cli_parser import run_program
from mmtbx.hydrogens import reduce_hydrogen
from mmtbx.programs import reduce2
from mmtbx.regression import reduce2_neutron_utils as u
from mmtbx.hydrogens import tst_add_hydrogen_11 as t11
from mmtbx.hydrogens import tst_add_hydrogen_16 as t16
from mmtbx.hydrogens import tst_add_hydrogen_17 as t17
from mmtbx.hydrogens import tst_add_hydrogen_19 as t19
from mmtbx.hydrogens import tst_add_hydrogen_23 as t23
from mmtbx.regression import tst_reduce2_metal_nh as tnh
from mmtbx.regression import tst_reduce2_metal_thiol as tth
from mmtbx.regression import tst_reduce2_blank_parent_split as tbp
from mmtbx.regression import tst_reduce2_na_terminal_h as tna

FIELDS = ('model_id', 'chain', 'resseq', 'icode', 'resname', 'altloc',
          'h_name', 'parent', 'cause', 'detail')

# ethanolamine with an sp3 NH2 (H-N-H 109.5): riding H refuses both H
eta_atoms = [
  ("N1", "N", "NH2", -0.700, 1.210, 0.0), ("C1", "C", "CH2", 0.0, 0.0, 0.0),
  ("C2", "C", "CH2", 1.520, 0.0, 0.0), ("O1", "O", "OH1", 2.220, 1.210, 0.0),
  ("HN1", "H", "HNH2", -0.200, 1.900, 0.500),
  ("HN2", "H", "HNH2", -1.200, 1.500, -0.800),
  ("H11", "H", "HCH2", -0.330, -0.550, 0.890),
  ("H12", "H", "HCH2", -0.330, -0.550, -0.890),
  ("H21", "H", "HCH2", 1.850, -0.550, 0.890),
  ("H22", "H", "HCH2", 1.850, -0.550, -0.890),
  ("HO1", "H", "HOH1", 3.120, 1.000, 0.0)]
eta_cif = u.make_cif("01", eta_atoms,
  [("C1", "N1", 1.47, 1.47), ("C1", "C2", 1.52, 1.52),
   ("C2", "O1", 1.43, 1.43), ("N1", "HN1", 0.89, 1.01),
   ("N1", "HN2", 0.89, 1.01), ("O1", "HO1", 0.85, 0.98),
   ("C1", "H11", 0.97, 1.09), ("C1", "H12", 0.97, 1.09),
   ("C2", "H21", 0.97, 1.09), ("C2", "H22", 0.97, 1.09)],
  [(x[0], x[1], x[2], 109.5) for x in
    [("C1", "N1", "HN1"), ("C1", "N1", "HN2"), ("HN1", "N1", "HN2"),
     ("N1", "C1", "C2"), ("C1", "C2", "O1"), ("C2", "O1", "HO1")] +
    u._tetra("C1", "N1", "C2") + u._tetra("C2", "C1", "O1")])

# 7T5E NAG C 1 - NAG C 2, linked BETA1-4 (C1 of NAG 2 to O4 of NAG 1); O1
# is in neither
pdb_nag_nag = """\
CRYST1   68.296   42.269   70.413  90.00  98.47  90.00 P 1 21 1
HETATM 3383  C1  NAG C   1     -10.860  -4.178 -40.436  1.00 27.06           C
HETATM 3384  C2  NAG C   1     -10.681  -2.675 -40.638  1.00 25.95           C
HETATM 3385  C3  NAG C   1      -9.195  -2.328 -40.662  1.00 32.04           C
HETATM 3386  C4  NAG C   1      -8.473  -2.921 -39.456  1.00 24.81           C
HETATM 3387  C5  NAG C   1      -8.814  -4.402 -39.257  1.00 26.90           C
HETATM 3388  C6  NAG C   1      -8.349  -4.931 -37.923  1.00 24.05           C
HETATM 3389  C7  NAG C   1     -12.458  -1.556 -41.962  1.00 32.11           C
HETATM 3390  C8  NAG C   1     -12.975  -1.282 -43.346  1.00 38.01           C
HETATM 3391  N2  NAG C   1     -11.331  -2.280 -41.886  1.00 29.56           N
HETATM 3392  O3  NAG C   1      -9.051  -0.912 -40.647  1.00 37.51           O
HETATM 3393  O4  NAG C   1      -7.066  -2.854 -39.649  1.00 32.97           O
HETATM 3394  O5  NAG C   1     -10.230  -4.577 -39.274  1.00 25.62           O
HETATM 3395  O6  NAG C   1      -8.560  -6.336 -37.806  1.00 34.61           O
HETATM 3396  O7  NAG C   1     -13.016  -1.117 -40.977  1.00 39.12           O
HETATM 3397  C1  NAG C   2      -6.436  -1.792 -39.032  1.00 39.57           C
HETATM 3398  C2  NAG C   2      -4.930  -2.007 -39.190  1.00 39.58           C
HETATM 3399  C3  NAG C   2      -4.176  -1.020 -38.317  1.00 46.39           C
HETATM 3400  C4  NAG C   2      -4.564   0.390 -38.715  1.00 39.21           C
HETATM 3401  C5  NAG C   2      -6.077   0.572 -38.599  1.00 43.91           C
HETATM 3402  C6  NAG C   2      -6.529   1.927 -39.101  1.00 54.83           C
HETATM 3403  C7  NAG C   2      -4.374  -4.309 -39.775  1.00 47.57           C
HETATM 3404  C8  NAG C   2      -3.968  -5.662 -39.261  1.00 35.21           C
HETATM 3405  N2  NAG C   2      -4.542  -3.365 -38.856  1.00 33.73           N
HETATM 3406  O3  NAG C   2      -2.774  -1.242 -38.441  1.00 49.87           O
HETATM 3407  O4  NAG C   2      -3.912   1.331 -37.868  1.00 72.64           O
HETATM 3408  O5  NAG C   2      -6.772  -0.429 -39.367  1.00 39.07           O
HETATM 3409  O6  NAG C   2      -7.521   1.844 -40.116  1.00 67.90           O
HETATM 3410  O7  NAG C   2      -4.563  -4.085 -40.967  1.00 41.41           O
END
"""

def dne(resseq='2'):
  """DNE G of tst_add_hydrogen_23 alone (no peptide bond to its N), as
  residue resseq."""
  return '\n'.join(l.replace('DNE G   2', 'DNE G %3s' % resseq)
    for l in t23.pdb_str.splitlines() if 'HT7' not in l) + '\n'

# geostd ideal coordinates of 0A0 and 004 as residue 1: their dictionaries
# name the H on N HN1, HN2 (0A0) and HN, HNA (004); chain-start placement adds
# H1, H2, H3, names these dictionaries do not know
pdb_0a0 = """\
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
HETATM    1  N   0A0 A   1      24.101   8.113 -14.105  1.00 20.00           N
HETATM    2  CA  0A0 A   1      23.417   9.374 -13.816  1.00 20.00           C
HETATM    3  C   0A0 A   1      24.034  10.606 -14.545  1.00 20.00           C
HETATM    4  O   0A0 A   1      23.726  11.736 -14.108  1.00 20.00           O
HETATM    5  CB  0A0 A   1      23.377   9.577 -12.304  1.00 20.00           C
HETATM    6  CG  0A0 A   1      24.721   9.728 -11.566  1.00 20.00           C
HETATM    7  OD1 0A0 A   1      25.780   9.685 -12.231  1.00 20.00           O
HETATM    8  OD2 0A0 A   1      24.641   9.879 -10.329  1.00 20.00           O
HETATM    9  OXT 0A0 A   1      24.725  10.376 -15.560  1.00 20.00           O
HETATM   10  CM  0A0 A   1      21.979   9.248 -14.335  1.00 20.00           C
END
"""

pdb_004 = """\
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
HETATM    1  C   004 A   1      30.780   2.418  52.694  1.00 20.00           C
HETATM    2  N   004 A   1      32.246   0.620  53.600  1.00 20.00           N
HETATM    3  O   004 A   1      31.476   3.170  53.401  1.00 20.00           O
HETATM    4  CA  004 A   1      31.205   0.915  52.634  1.00 20.00           C
HETATM    5  CB  004 A   1      29.996   0.035  52.845  1.00 20.00           C
HETATM    6  CE  004 A   1      27.695  -1.524  53.236  1.00 20.00           C
HETATM    7  CD1 004 A   1      28.399  -1.024  54.325  1.00 20.00           C
HETATM    8  CD2 004 A   1      28.141  -1.242  51.952  1.00 20.00           C
HETATM    9  CG1 004 A   1      29.538  -0.254  54.130  1.00 20.00           C
HETATM   10  CG2 004 A   1      29.281  -0.471  51.764  1.00 20.00           C
HETATM   11  OXT 004 A   1      29.801   2.728  51.993  1.00 20.00           O
END
"""

def with_atoms(pdb_str, resname, resseq, atoms):
  """pdb_str with atoms (name, element, x, y, z) added after the last line
  of resname resseq."""
  lines = pdb_str.rstrip('\n').splitlines()
  i = max(k for k, l in enumerate(lines)
          if l[17:20] == resname and l[22:26].strip() == resseq)
  new = ['%-6s%5d %4s %3s %1s%4s    %8.3f%8.3f%8.3f  1.00 20.00          %2s'
    % (lines[i][:6], 90 + j, (' %-3s' % n) if len(n) < 4 else n, resname,
       lines[i][21], resseq, x, y, z, e)
    for j, (n, e, x, y, z) in enumerate(atoms)]
  return '\n'.join(lines[:i+1] + new + lines[i+1:]) + '\n'

def gln_supplied(*names):
  """The Gln cut at CB (pdb_gln()) with supplied H/D of these names."""
  return with_atoms(pdb_gln(), 'GLN', '2', [(n, n[0], 10. + 2 * i, 10., 10.)
    for i, n in enumerate(names)])

def pdb_0a0_supplied_h1():
  """pdb_0a0 with a supplied H1 1 A from N."""
  n = [l for l in pdb_0a0.splitlines() if l[12:16] == ' N  '][0]
  return with_atoms(pdb_0a0, '0A0', '1', [('H1', 'H', float(n[30:38]),
    float(n[38:46]) + 1.0, float(n[46:54]))])

def split_ab(pdb_str, keep_blank=lambda l: False):
  """The HETATM lines of pdb_str in conformers A and B (B 0.1 A along x),
  except those keep_blank selects."""
  lines = []
  for l in pdb_str.splitlines():
    if l.startswith('HETATM') and not keep_blank(l):
      lines.append(l[:16] + 'A' + l[17:54] + '  0.50' + l[60:])
      lines.append(l[:16] + 'B' + l[17:30] + '%8.3f' % (float(l[30:38]) + 0.1)
        + l[38:54] + '  0.50' + l[60:])
    else:
      lines.append(l)
  return '\n'.join(lines) + '\n'

def nag_019(drop=(), asn=True):
  """pdb_str_019_control (Asn-NAG), without the NAG atoms drop (names), and
  without the Asn unless asn."""
  return '\n'.join(l for l in t16.pdb_str_019_control.splitlines()
    if not (l[17:20] == 'NAG' and l[12:16].strip() in drop)
    and (asn or l[17:20] != 'ASN')) + '\n'

def met_gln():
  '''pdb_str_005's MET A 1 - GLN A 2 (GLN cut at CB), with CRYST1.'''
  lines = t11.pdb_str_005.strip().splitlines()
  return [lines[0]] + [l for l in lines[1:]
    if l.startswith('ATOM') and int(l[22:26]) <= 2]

def gln_cg(altloc=' ', resname='GLN'):
  '''A CG for GLN A 2: 1.52 A from CB, 109.5 deg from CA, anti to N.'''
  xyz = {}
  for l in met_gln()[1:]:
    if int(l[22:26]) == 2:
      xyz[l[12:16].strip()] = matrix.col((float(l[30:38]), float(l[38:46]),
        float(l[46:54])))
  ca, cb, n = xyz['CA'], xyz['CB'], xyz['N']
  e = (cb - ca).normalize()
  p = ((n - ca) - e * (n - ca).dot(e)).normalize()
  a = math.radians(70.5)
  cg = cb + (e * math.cos(a) - p * math.sin(a)) * 1.52
  return ('ATOM     99  CG %s%s A   2    %8.3f%8.3f%8.3f  1.00115.53'
    '           C' % ((altloc, resname) + cg.elems))

def pdb_gln():
  return '\n'.join(met_gln()) + '\n'

def pdb_lys():
  # Lys cut at CG
  return '\n'.join([l.replace('GLN', 'LYS') for l in met_gln()] +
    [gln_cg(resname='LYS')]) + '\n'

def pdb_gln_cg_in_a():
  # CB in A and B, CG in A only: in B the parent CG is not in the model
  out = []
  for l in met_gln():
    if l.startswith('ATOM') and int(l[22:26]) == 2 and l[12:16] == ' CB ':
      out.append(l[:16] + 'A' + l[17:54] + '  0.50' + l[60:])
      out.append(l[:16] + 'B' + l[17:30] + '%8.3f' % (float(l[30:38]) + 0.3)
        + l[38:54] + '  0.50' + l[60:])
    else:
      out.append(l)
  cg = gln_cg('A')
  return '\n'.join(out + [cg[:54] + '  0.50' + cg[60:]]) + '\n'

def pdb_water_no_o():
  return t17.cryst1 + t17.water(1, 5.0, o=False) + t17.water(2, 15.0)

def place(pdb_str, restraint_objects=None, **kw):
  inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  model = mmtbx.model.manager(model_input=inp,
    restraint_objects=restraint_objects, log=null_out())
  o = reduce_hydrogen.place_hydrogens(model=model, stop_for_unknowns=False,
    **kw)
  o.run()
  return o

def h_keys(hierarchy):
  return set((a.parent().parent().parent().parent().id.strip(),
    a.parent().parent().parent().id.strip(), a.parent().parent().resseq.strip(),
    a.parent().parent().icode.strip(), a.parent().resname.strip(),
    a.parent().altloc.strip(), a.name.strip()) for a in hierarchy.atoms()
    if a.element.strip().upper() in ('H', 'D'))

def added_h(bogus, after_propeller=None):
  '''The H placement added (or kept from the input): those of the
  bogus-position pass and of the propeller, without the H that the
  propeller replaced (a chain start's H by H1, H2, H3).'''
  if after_propeller is None: return set(bogus)
  new = after_propeller - bogus
  replaced = set(k for k in bogus - after_propeller
                 if [x for x in new if x[:6] == k[:6]])
  return (bogus | after_propeller) - replaced

def place_counted(pdb_str, restraint_objects=None, **kw):
  '''place(), and added_h() of that run.'''
  added = []
  cls = reduce_hydrogen.place_hydrogens
  bogus, prop = cls.add_missing_H_atoms_at_bogus_position, \
    cls.place_n_terminal_propeller
  def bogus_counted(self, *a, **k):
    h = bogus(self, *a, **k)
    added.append(h_keys(h))
    return h
  def prop_counted(self, pdb_hierarchy):
    r = prop(self, pdb_hierarchy)
    added.append(h_keys(pdb_hierarchy))
    return r
  cls.add_missing_H_atoms_at_bogus_position = bogus_counted
  cls.place_n_terminal_propeller = prop_counted
  try:
    o = place(pdb_str, restraint_objects, **kw)
  finally:
    cls.add_missing_H_atoms_at_bogus_position = bogus
    cls.place_n_terminal_propeller = prop
  return o, added_h(*added)

def unaccounted(records, added, final):
  '''H added and removed without exactly one record, records of H that are
  in the final model, and records repeated (letters_differ aside).'''
  n = collections.Counter((r.model_id, r.chain, r.resseq, r.icode, r.resname,
    r.altloc, r.h_name) for r in records
    if r.h_name and r.cause != 'letters_differ')
  return sorted([('no record', k) for k in added - final if n[k] != 1] +
    [('placed', k) for k in final if n[k]] +
    [('repeated', k) for k in n if n[k] > 1])

def records(o, cause=None):
  return [tuple(getattr(r, f) for f in FIELDS) for r in o.not_placed
          if cause is None or r.cause == cause]

def check(o, cause, expected, lines):
  '''The records of cause are expected (all fields), and the description
  text has its header and each of lines.'''
  found = records(o, cause)
  assert found == expected, found
  text = o.not_placed_text()
  title = [c.title for c in reduce_hydrogen.NOT_PLACED_CAUSES
           if c.key == cause][0]
  assert [l for l in text.splitlines() if l.startswith(title + ':')], text
  for line in lines:
    assert line in text.splitlines(), (line, text)

def test_missing_neighbour_and_parent():
  # Gln cut at CB: HB2/HB3 on CB miss CG; the H of CG and NE2 have no parent
  o, _ = placed(pdb_gln())
  check(o, 'missing_neighbour',
    [('', 'A', '2', '', 'GLN', '', 'HB2', 'CB', 'missing_neighbour', ''),
     ('', 'A', '2', '', 'GLN', '', 'HB3', 'CB', 'missing_neighbour', '')],
    ['  GLN A 2: HB2 HB3 on CB'])
  check(o, 'parent_missing',
    [('', 'A', '2', '', 'GLN', '', h, p, 'parent_missing', '')
     for h, p in (('HG2', 'CG'), ('HG3', 'CG'), ('HE21', 'NE2'),
                  ('HE22', 'NE2'))],
    ['  GLN A 2: HG2 HG3 on CG; HE21 HE22 on NE2',
     '  What to do: Check whether the parent is missing from an incomplete '
     'residue or absent because of a link or modification not identified '
     'here.'])

def test_parent_missing_lys():
  o, _ = placed(pdb_lys())
  check(o, 'missing_neighbour',
    [('', 'A', '2', '', 'LYS', '', h, 'CG', 'missing_neighbour', '')
     for h in ('HG2', 'HG3')], ['  LYS A 2: HG2 HG3 on CG'])
  check(o, 'parent_missing',
    [('', 'A', '2', '', 'LYS', '', h, p, 'parent_missing', '')
     for h, p in (('HD2', 'CD'), ('HD3', 'CD'), ('HE2', 'CE'), ('HE3', 'CE'),
                  ('HZ1', 'NZ'), ('HZ2', 'NZ'), ('HZ3', 'NZ'))],
    ['  LYS A 2: HD2 HD3 on CD; HE2 HE3 on CE; HZ1 HZ2 HZ3 on NZ'])

def test_parent_in_one_conformer():
  # CG only in A: HG2/HG3 of B have no parent, those of A miss CD
  o, _ = placed(pdb_gln_cg_in_a())
  assert [r[5:9] for r in records(o, 'parent_missing')] == [
    ('A', 'HE21', 'NE2', 'parent_missing'),
    ('A', 'HE22', 'NE2', 'parent_missing'),
    ('B', 'HG2', 'CG', 'parent_missing'), ('B', 'HG3', 'CG', 'parent_missing'),
    ('B', 'HE21', 'NE2', 'parent_missing'),
    ('B', 'HE22', 'NE2', 'parent_missing')], records(o, 'parent_missing')
  assert [r[5:9] for r in records(o, 'missing_neighbour')] == [
    ('A', 'HG2', 'CG', 'missing_neighbour'),
    ('A', 'HG3', 'CG', 'missing_neighbour'),
    ('B', 'HB2', 'CB', 'missing_neighbour'),
    ('B', 'HB3', 'CB', 'missing_neighbour')], records(o, 'missing_neighbour')
  text = o.not_placed_text().splitlines()
  assert '  GLN A 2 conformer B: HG2 HG3 on CG; HE21 HE22 on NE2' in text

def test_water_without_o():
  # supplied H of a water without O: no parent; one record each, nothing
  # under missing_neighbour
  o, _ = placed(pdb_water_no_o(), exclude_water=False, keep_existing_H=True)
  assert records(o) == [
    ('', 'A', '1', '', 'HOH', '', h, 'O', 'parent_missing', 'water without O')
    for h in ('H1', 'H2')], records(o)
  assert '  HOH A 1: H1 H2 on O (water without O)' in \
    o.not_placed_text().splitlines()

_placed = {}

def placed(pdb_str, **kw):
  """place_counted() once per input and options, shared by the tests that
  only read the result."""
  key = (pdb_str, tuple(sorted(kw.items())))
  if key not in _placed:
    _placed[key] = place_counted(pdb_str, **kw)
  return _placed[key]

def test_accounting():
  # every H added and removed has exactly one record; H that are placed have
  # none; the check fails with a record removed or repeated
  # an empty parent (an H name the dictionary does not know) is never
  # parent_missing
  cases = [(pdb_gln(), {}), (pdb_lys(), {}), (pdb_gln_cg_in_a(), {}),
    (pdb_water_no_o(), {'exclude_water': False, 'keep_existing_H': True}),
    (t16.pdb_str_009, {}), (t16.pdb_str_015, {}), (t16.pdb_str_019_control, {}),
    (pdb_nag_nag, {}), (nag_019(drop=('O6',)), {}), (t16.pdb_str_021, {}),
    (tnh.pdb_1xmk, {'metal_bound_amines_as_nh2': True}), (pdb_0a0, {}),
    (pdb_004, {}), (gln_supplied('HE21', 'DE21'), {'keep_existing_H': True}),
    (gln_supplied('HXX'), {'keep_existing_H': True}),
    (pdb_0a0_supplied_h1(), {'keep_existing_H': True}),
    (split_ab(pdb_0a0), {}),
    (split_ab(pdb_0a0, lambda l: l[12:16] == ' N  '), {})]
  results = [placed(pdb_str, **kw) for pdb_str, kw in cases]
  # restraint objects are not a key of placed()
  results.append(place_counted(t19.cryst1 + t19.pdb_protein + t19.pdb_04,
    restraint_objects = ro_04()))
  for o, added in results:
    final = h_keys(o.get_model().get_hierarchy())
    assert added - final, 'nothing removed'
    assert unaccounted(o.not_placed, added, final) == [], \
      unaccounted(o.not_placed, added, final)
    rs = [r for r in o.not_placed if r.h_name]
    assert unaccounted(rs[1:], added, final)
    assert unaccounted(rs + rs[:1], added, final)
    assert [r for r in o.not_placed
            if r.cause == 'parent_missing' and not r.parent] == []
  # a gap end's HO3' is never added, so added minus final cannot see it: its
  # record is asserted explicitly
  o, added = placed(tna.residues(tna.dna, 1, 3), na_terminal_h='chain_ends')
  final = h_keys(o.get_model().get_hierarchy())
  assert unaccounted(o.not_placed, added, final) == []
  gap = [('', 'A', '1', '', 'DG', '', "HO3'", "O3'", 'na_gap_end', '')]
  assert records(o, 'na_gap_end') == gap, records(o)
  # DNE alone, before HT7, as a chain start: no record of an H in the final
  # model, every H removed recorded (their free-form H is never added)
  for pdb_str, kw in ((dne(), {}), (t23.pdb_str, {}),
      (dne(), {'n_terminal_charge': 'first_in_chain'}), (dne('1'), {})):
    o, added = placed(pdb_str, **kw)
    final = h_keys(o.get_model().get_hierarchy())
    assert unaccounted(o.not_placed, added, final) == [], \
      unaccounted(o.not_placed, added, final)

def test_supplied_names():
  # with keep_existing_H, the parent is looked up by the supplied name, then
  # as H; the record keeps the name: HE21 and DE21 on a Gln without NE2
  o, _ = placed(gln_supplied('HE21', 'DE21'), keep_existing_H=True)
  assert [r for r in records(o, 'parent_missing') if r[6] in
    ('HE21', 'DE21')] == [('', 'A', '2', '', 'GLN', '', h, 'NE2',
    'parent_missing', '') for h in ('DE21', 'HE21')], records(o)
  # an element-D atom with an H name, HD21 on an Asn without ND2: found by
  # its own name (as H it would be HH21)
  asn = '\n'.join(l for l in t16.pdb_str_019_control.splitlines()
    if l[17:20] == 'ASN' and l[12:16] != ' ND2') + '\n'
  o, _ = placed(with_atoms(asn, 'ASN', '114', [('HD21', 'D', 10., 10., 10.)]),
    keep_existing_H=True)
  assert [r for r in records(o) if r[6] == 'HD21'] == [('', 'A', '114', '',
    'ASN', '', 'HD21', 'ND2', 'parent_missing', '')], records(o)
  # a supplied name the dictionary does not know, also as H
  o, _ = placed(gln_supplied('HXX'), keep_existing_H=True)
  assert records(o, 'unknown_atoms') == [('', 'A', '2', '', 'GLN', '', 'HXX',
    '', 'unknown_atoms', '')], records(o)
  assert o.site_labels_unknown_atoms == [], o.site_labels_unknown_atoms
  # a supplied H1 on the chain-start N of 0A0: not from chain-start
  # placement, which adds H2 and H3
  o, _ = placed(pdb_0a0_supplied_h1(), keep_existing_H=True)
  detail = ('added by chain-start placement, name not in the 0A0 '
    'dictionary; no bond restraint in the processed model')
  assert records(o) == [('', 'A', '1', '', '0A0', '', 'H1', '',
    'unknown_atoms', '')] + [('', 'A', '1', '', '0A0', '', h, 'N',
    'not_parameterized', detail) for h in ('H2', 'H3')], records(o)

def test_unknown_atoms_supplied():
  # a lone supplied H with a name the dictionary does not know, in a residue
  # with an unknown heavy atom (ASP CGX): listed with the residue's unknown
  # atoms like HB2 and HB3, once each
  pdb_str = with_atoms(t16.pdb_str_021, 'ASP', '11',
    [('HXX', 'H', 10., 10., 10.)])
  o, _ = placed(pdb_str, keep_existing_H=True)
  labels = [l.split()[0] for l in o.site_labels_unknown_atoms]
  assert labels == ['HB2', 'HB3', 'HXX'], o.site_labels_unknown_atoms
  assert [r[6] for r in records(o, 'unknown_atoms')] == labels, records(o)
  verbose = o.unknown_atoms_text(verbose=True).splitlines()
  assert [l.split()[0] for l in verbose if l.startswith('   ')] == labels, \
    verbose
  count = '1 residue; 3 H atoms on them were not placed:'
  assert count in verbose[1], verbose
  description, log = run('unknown_supplied', pdb_str, 'keep_existing_H=True')
  lines = description.splitlines()
  assert [l for l in lines if count in l], description
  assert '  ASP A 11: HB2 HB3 on CB; HXX' in lines, description
  assert [l for l in log.splitlines() if count in l], log
  assert ('  residue with atom names not in its dictionary: 3 H, 1 residue: '
    'ASP A 11' in log.splitlines()), log

def test_chain_start_unknown_names():
  # 0A0 and 004 at residue 1: chain-start placement adds H1, H2, H3 on N,
  # names their dictionaries do not know; riding H has no bond for them and
  # they are removed: not_parameterized on N, not parent_missing
  for pdb_str, resname in ((pdb_0a0, '0A0'), (pdb_004, '004')):
    o, _ = placed(pdb_str)
    detail = ('added by chain-start placement, name not in the %s '
      'dictionary; no bond restraint in the processed model' % resname)
    assert records(o, 'not_parameterized') == [('', 'A', '1', '', resname,
      '', h, 'N', 'not_parameterized', detail) for h in ('H1', 'H2', 'H3')], \
      records(o)
    assert records(o, 'parent_missing') == [], records(o)
    assert ('  %s A 1: H1 H2 H3 on N (%s)' % (resname, detail) in
      o.not_placed_text().splitlines()), o.not_placed_text()
  # 0A0 in A and B, and with a blank N and the rest in A and B: H1, H2, H3
  # in each conformer, on N
  detail = ('added by chain-start placement, name not in the 0A0 '
    'dictionary; no bond restraint in the processed model')
  for pdb_str in (split_ab(pdb_0a0),
                  split_ab(pdb_0a0, lambda l: l[12:16] == ' N  ')):
    o, _ = placed(pdb_str)
    assert records(o) == [('', 'A', '1', '', '0A0', alt, h, 'N',
      'not_parameterized', detail) for alt in 'AB'
      for h in ('H1', 'H2', 'H3')], records(o)

def test_tertiary_amide():
  o = place(t11.pdb_str_004)
  check(o, 'tertiary_amide',
    [('', 'A', '10', '', 'HYP', '', 'H', 'N', 'tertiary_amide', '')],
    ['  HYP A 10: H on N'])

def test_link():
  o, _ = placed(t16.pdb_str_009)
  check(o, 'link',
    [('', 'D', '4', '', '9FZ', '', 'H15', 'C9', 'link', 'Misc. bond'),
     ('', 'D', '10', '', '9G2', '', 'H17', 'C10', 'link', 'Misc. bond')],
    ['  9FZ D 4: H15 on C9 (Misc. bond)',
     '  9G2 D 10: H17 on C10 (Misc. bond)'])

def test_esterified_o():
  o, _ = placed(t16.pdb_str_015)
  check(o, 'esterified_o',
    [('', 'A', '14', '', 'BGM', '', "HO3'", "O3'", 'esterified_o',
      "bonded to C3', P"),
     ('', 'A', '15', '', 'BGM', '', 'HOP2', 'OP2', 'esterified_o',
      'phosphodiester')],
    ["  BGM A 14: HO3' on O3' (bonded to C3', P)",
     '  BGM A 15: HOP2 on OP2 (phosphodiester)'])

def test_link_angle():
  o, _ = placed(t16.pdb_str_019_control)
  check(o, 'link_angle',
    [('', 'A', '114', '', 'ASN', '', 'HD22', 'ND2', 'link_angle',
      'link to C1 NAG A 401')],
    ['  ASN A 114: HD22 on ND2 (link to C1 NAG A 401)'])

def test_link_leaving():
  # Asn-NAG: O1 is the leaving atom of the NAG-ASN link; no completion hint
  o, _ = placed(t16.pdb_str_019_control)
  check(o, 'link_leaving',
    [('', 'A', '401', '', 'NAG', '', 'HO1', 'O1', 'link_leaving',
      'NAG-ASN link to ND2 ASN A 114')],
    ['  NAG A 401: HO1 on O1 (NAG-ASN link to ND2 ASN A 114)'])
  lines = o.not_placed_text().splitlines()
  i = [i for i, l in enumerate(lines)
       if l.startswith('The parent atom is a leaving atom of a link')][0]
  assert lines[i+1].startswith('  NAG A 401: '), lines
  assert records(o, 'parent_missing') == [], records(o, 'parent_missing')
  # NAG-NAG (BETA1-4): NAG 2's O1 leaves; NAG 1 has no link through C1
  o, _ = placed(pdb_nag_nag)
  assert records(o, 'link_leaving') == [('', 'C', '2', '', 'NAG', '', 'HO1',
    'O1', 'link_leaving', 'BETA1-4 link to O4 NAG C 1')], records(o)
  assert records(o, 'parent_missing') == [('', 'C', '1', '', 'NAG', '',
    'HO1', 'O1', 'parent_missing', '')], records(o)
  # an unlinked NAG without O1: parent_missing
  o = place(nag_019(asn=False))
  assert records(o, 'link_leaving') == [], records(o)
  assert records(o, 'parent_missing') == [('', 'A', '401', '', 'NAG', '',
    'HO1', 'O1', 'parent_missing', '')], records(o)
  # NAG in A and B, B 1.5 A away along -x: only A is linked to the Asn
  lines = []
  for l in t16.pdb_str_019_control.splitlines():
    if l[17:20] == 'NAG':
      lines.append(l[:16] + 'A' + l[17:54] + '  0.50' + l[60:])
      lines.append(l[:16] + 'B' + l[17:30] + '%8.3f' % (float(l[30:38]) - 1.5)
        + l[38:54] + '  0.50' + l[60:])
    else:
      lines.append(l)
  o = place('\n'.join(lines) + '\n')
  assert [r for r in records(o) if r[6] == 'HO1'] == [
    ('', 'A', '401', '', 'NAG', 'A', 'HO1', 'O1', 'link_leaving',
     'NAG-ASN link to ND2 ASN A 114'),
    ('', 'A', '401', '', 'NAG', 'B', 'HO1', 'O1', 'parent_missing', '')], \
    records(o)
  # linked NAG without O6: HO6 parent_missing, HO1 link_leaving
  o, _ = placed(nag_019(drop=('O6',)))
  assert records(o, 'parent_missing') == [('', 'A', '401', '', 'NAG', '',
    'HO6', 'O6', 'parent_missing', '')], records(o)
  assert [r[6] for r in records(o, 'link_leaving')] == ['HO1'], records(o)

def test_unknown_atoms():
  o, _ = placed(t16.pdb_str_021)
  check(o, 'unknown_atoms',
    [('', 'A', '11', '', 'ASP', '', 'HB2', 'CB', 'unknown_atoms', ''),
     ('', 'A', '11', '', 'ASP', '', 'HB3', 'CB', 'unknown_atoms', '')],
    ['  ASP A 11: HB2 HB3 on CB'])

def pointer(title):
  return "  Listed in full under '%s'." % title

def test_metal_nh_and_amine():
  # the pointers quote the titles of metal_nh_text() and metal_amine_text()
  o, _ = placed(tnh.pdb_1xmk, metal_bound_amines_as_nh2=True)
  check(o, 'metal_nh',
    [('', 'A', '289', '', 'SER', '', 'H', 'N', 'metal_nh',
      'N   SER A 289 - NI    NI A 398 1.80 A (geometric)'),
     ('', 'A', '290', '', 'HIS', '', 'H', 'N', 'metal_nh',
      'N   HIS A 290 - NI    NI A 398 1.84 A (geometric)')],
    [pointer(reduce_hydrogen.METAL_NH_NOT_PLACED)])
  check(o, 'metal_amine',
    [('', 'A', '288', '', 'GLY', '', 'H1', 'N', 'metal_amine',
      'N   GLY A 288 - NI    NI A 398 1.89 A (geometric)')],
    [pointer(reduce_hydrogen.METAL_AMINE_CONVERTED)])
  assert reduce_hydrogen.METAL_NH_NOT_PLACED in \
    o.metal_nh_text().splitlines()
  assert reduce_hydrogen.METAL_AMINE_CONVERTED in \
    o.metal_amine_text().splitlines()

def test_metal_thiolate():
  o = place(tth.m_zn_fixture)
  check(o, 'metal_thiolate',
    [('', 'A', '48', '', 'CYS', '', 'HG', 'SG', 'metal_thiolate',
      'SG  CYS A  48 - ZN    ZN A 102 2.39 A (metal coordination proxy)')],
    [pointer(reduce_hydrogen.METAL_THIOL_NOT_PLACED)])
  assert reduce_hydrogen.METAL_THIOL_NOT_PLACED in \
    o.metal_thiol_text().splitlines()

def test_na_gap_end():
  o, _ = placed(tna.residues(tna.dna, 1, 3), na_terminal_h='chain_ends')
  title = '1 free nucleic-acid end inside a chain (chain breaks) got no ' \
    'hydroxyl H'
  check(o, 'na_gap_end',
    [('', 'A', '1', '', 'DG', '', "HO3'", "O3'", 'na_gap_end', '')],
    [pointer(title)])
  assert o.na_terminal_gap_text().splitlines()[0] == title + ':'

def test_letters_differ():
  # placed, not split: own heading, no 'not placed' header
  o, _ = placed(tbp.pdb_letters)
  check(o, 'letters_differ',
    [('', 'A', '101', '', 'LEU', 'C', 'H', 'N', 'letters_differ',
      'residue letters C, D; partner C LYS A 100 letters A, B')],
    ['  LEU A 101 conformer C: H on N (residue letters C, D; partner C LYS '
     'A 100 letters A, B)'])
  text = o.not_placed_text()
  assert 'not placed' not in text, text
  assert text.splitlines()[0].startswith('H placed once, not split'), text

def test_letters_differ_summary():
  # the log: placed H under their own unindented heading, alone and next to
  # H not placed
  letters = placed(tbp.pdb_letters)[0].not_placed
  heading = 'H placed, not split per conformer: 1 H, 1 residue: LEU A 101'
  summary = reduce_hydrogen.not_placed_summary(letters, {})
  assert summary.splitlines() == [heading], summary
  gln = placed(pdb_gln())[0].not_placed
  summary = reduce_hydrogen.not_placed_summary(gln + letters, {})
  assert summary.splitlines() == [
    'H atoms not placed: 6 H, 1 residue',
    '  parent atom not in the model: 4 H, 1 residue: GLN A 2',
    '  incomplete side chain: 2 H, 1 residue: GLN A 2',
    heading], summary
  # the rules and the pointer are not under the placed heading
  rules = {'waters': 2, 'free_form': 0}
  summary = reduce_hydrogen.not_placed_summary(letters, rules, see='x',
    water_rule='add_h_to_water=False')
  assert summary.splitlines() == [
    'Rules, not listed: 2 waters get no H (add_h_to_water=False).',
    'Full list: see x.', heading], summary
  summary = reduce_hydrogen.not_placed_summary(gln + letters, rules,
    see='x', water_rule='add_h_to_water=False')
  assert summary.splitlines()[-3:] == [
    '  Rules, not listed: 2 waters get no H (add_h_to_water=False).',
    '  Full list: see x.', heading], summary

def test_no_dictionary():
  # not assessed: counted as a residue, not as H
  o = place(t19.pdb_04)
  check(o, 'no_dictionary',
    [('', 'B', '1', '', '04', '', '', '', 'no_dictionary',
      'no dictionary in the library, the CCD or the restraint files')],
    ['  04 B 1: all H (no dictionary in the library, the CCD or the '
     'restraint files)'])
  assert o.not_placed_text().splitlines()[0] == \
    'not assessed (no dictionary): 1 residue', o.not_placed_text()

def ro_04():
  return [("04.cif", iotbx.cif.reader(
    input_string=t19.cif_04_no_bonds).model())]

def test_no_bonds():
  # a dictionary without bonds: one record per H, by full identity, in each
  # model and next to a residue with an insertion code
  h_04 = ('H11', 'H12', 'H13', 'H21', 'H22', 'HO1')
  def no_bonds(model_id, icode):
    return [(model_id, 'B', '1', icode, '04', '', h, '', 'no_restraints',
             'no bonds in the processed model') for h in h_04]
  o = place(t19.cryst1 + t19.pdb_protein + t19.pdb_04, ro_04())
  check(o, 'no_restraints', no_bonds('', ''),
    ['  04 B 1: %s (no bonds in the processed model)' % ' '.join(h_04)])
  two = (t19.cryst1 + 'MODEL        1\n' + t19.pdb_protein + t19.pdb_04 +
    'ENDMDL\nMODEL        2\n' + t19.pdb_protein + t19.pdb_04 + 'ENDMDL\n')
  o = place(two, ro_04())
  assert records(o, 'no_restraints') == no_bonds('1', '') + no_bonds('2', '')
  icode = '\n'.join(l[:26] + 'A' + l[27:30] + '%8.3f' % (float(l[30:38]) + 6)
    + l[38:] for l in t19.pdb_04.splitlines())
  o = place(t19.cryst1 + t19.pdb_protein + t19.pdb_04 + icode + '\n', ro_04())
  assert records(o, 'no_restraints') == no_bonds('', '') + no_bonds('', 'A')

def test_riding_refused():
  o = place(u.heavy_pdb('01', eta_atoms),
    [('01.cif', iotbx.cif.reader(input_string=eta_cif).model())])
  check(o, 'riding_refused',
    [('', 'A', '1', '', '01', '', 'HN1', 'N1', 'riding_refused', ''),
     ('', 'A', '1', '', '01', '', 'HN2', 'N1', 'riding_refused', '')],
    ['  01 A 1: HN1 HN2 on N1'])

def test_not_parameterized():
  # two candidates for X at N (tst_add_hydrogen_23): no C-N-X, not placed
  o = place(t23.pdb_str.replace('HT7', 'X7X'), t23.ht7_cif(
    ' X7X  N    CG   single    1.456  0.020  1.456', 'X7X'))
  check(o, 'not_parameterized',
    [('', 'G', '3', '', 'X7X', '', 'HA2', 'N', 'not_parameterized', ''),
     ('', 'G', '3', '', 'X7X', '', 'HN1', 'N', 'not_parameterized', '')],
    ['  X7X G 3: HA2 HN1 on N'])

def gln_conformers():
  '''MET-GLN, GLN as A 2A in conformers A and B.'''
  out = []
  for l in met_gln()[1:]:
    if int(l[22:26]) == 1:
      out.append(l)
      continue
    l = l[:26] + 'A' + l[27:]
    out.append(l[:16] + 'A' + l[17:54] + '  0.60' + l[60:])
    x = '%8.3f' % (float(l[30:38]) + 0.1)
    out.append(l[:16] + 'B' + l[17:30] + x + l[38:54] + '  0.40' + l[60:])
  return '\n'.join(out)

def test_labels():
  # model id (two models), insertion code, conformer
  pdb_str = '\n'.join([met_gln()[0], 'MODEL        1', gln_conformers(),
    'ENDMDL', 'MODEL        2', gln_conformers(), 'ENDMDL', 'END'])
  o = place(pdb_str)
  found = records(o, 'missing_neighbour')
  assert found == [
    (m, 'A', '2', 'A', 'GLN', c, h, 'CB', 'missing_neighbour', '')
    for m in ('1', '2') for c in ('A', 'B') for h in ('HB2', 'HB3')], found
  text = o.not_placed_text().splitlines()
  for m in ('1', '2'):
    for c in ('A', 'B'):
      line = '  model %s GLN A 2A conformer %s: HB2 HB3 on CB' % (m, c)
      assert line in text, (line, text)
  assert text[0] == 'H atoms not placed: 24 H, 2 residues', text[0]
  # the log: residues without their conformers
  assert ('incomplete side chain: 8 H, 2 residues: model 1 GLN A 2A, '
    'model 2 GLN A 2A') in o.not_placed_summary(), o.not_placed_summary()

def run(prefix, pdb_str, *args):
  fn, out, txt = prefix + '.pdb', prefix + '_out.pdb', prefix + '.txt'
  with open(fn, 'w') as f:
    f.write(pdb_str)
  for p in (out, txt):
    if os.path.exists(p): os.remove(p)
  log = StringIO()
  run_program(program_class=reduce2.Program, logger=log, args=[fn,
    'output.filename=' + out, 'output.description_file_name=' + txt,
    'overwrite=True'] + list(args))
  with open(txt) as f:
    description = f.read()
  return description, log.getvalue()

def gln_chains(n, extra=''):
  '''MET-GLN n times, chains A, B, ..., 30 A apart; extra: more records.'''
  lines = [met_gln()[0]]
  ids = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  for i in range(n):
    for l in met_gln()[1:]:
      lines.append(l[:21] + ids[i] + l[22:30] +
        '%8.3f' % (float(l[30:38]) + 30. * i) + l[38:])
    lines.append('TER')
  return '\n'.join(lines + ([extra.rstrip('\n')] if extra else []) +
    ['END']) + '\n'

def test_reduce2_description_and_log():
  description, log = run('not_placed_gln', gln_chains(1))
  lines = description.splitlines()
  i = lines.index('A heavy neighbour of the parent atom is missing '
    '(incomplete side chain): 2 H, 1 residue')
  assert lines[i+1] == ('  What to do: If the data support the missing '
    'atoms, complete the residue to get these H.'), lines[i+1]
  assert '  GLN A 2: HB2 HB3 on CB' in lines
  assert 'incomplete side chain: 2 H, 1 residue: GLN A 2' in log, log
  assert 'parent atom not in the model: 4 H, 1 residue: GLN A 2' in log, log
  assert 'Full list: see the description file.' in log, log
  assert 'Full list' not in description

def test_no_truncation():
  # 25 residues: all in the description; the log shows the first 3
  description, log = run('not_placed_25', gln_chains(25))
  lines = [l for l in description.splitlines()
           if l.startswith('  GLN ') and l.endswith(': HB2 HB3 on CB')]
  assert len(lines) == 25, len(lines)
  assert 'more' not in description and '...' not in description
  assert ('incomplete side chain: 50 H, 25 residues: GLN A 2, GLN B 2, '
    'GLN C 2, ...') in log, log
  assert 'H atoms not placed: 150 H, 25 residues' in description

def waters(n, with_h=()):
  '''n waters W 1..n; with H1/H2 for the resseq in with_h.'''
  lines = []
  for i in range(n):
    lines.append('HETATM%5d  O   HOH W%4d    %8.3f%8.3f%8.3f  1.00 20.00'
      '           O' % (100 + i, i + 1, -30. + 3 * i, -10., 60.))
    if i + 1 in with_h:
      for name, dx in (('H1', 0.95), ('H2', -0.3)):
        lines.append('HETATM%5d  %-3s HOH W%4d    %8.3f%8.3f%8.3f  1.00 20.00'
          '           H' % (200 + i, name, i + 1, -30. + 3 * i + dx,
          -10. + (0.9 if name == 'H2' else 0.), 60.))
  return '\n'.join(lines) + '\n'

def test_rule_counts():
  # three waters, and the free-form H of DNE and HT7 on N: counted, not
  # listed; placement's own wording names no option
  o = place(t23.pdb_str.replace('END\n', waters(3) + 'END\n'))
  assert o.not_placed_rules == {'waters': 3, 'free_form': 2}, \
    o.not_placed_rules
  assert records(o) == [], records(o)
  rules = ('Rules, not listed: 3 waters get no H from placement; 2 dictionary '
    'H entries on N or OXT of peptide-like residues skipped by the free-form '
    'rule during initial placement (counted per model and conformer). '
    'Chain-start placement may add H on these N atoms later.')
  text = o.not_placed_text()
  assert text.splitlines()[0] == rules, text
  assert 'exclude_water' not in text
  # the log: the rules alone, on a line of their own
  assert o.not_placed_summary() == rules, o.not_placed_summary()
  one = place(t23.pdb_str.replace('END\n', waters(1) + 'END\n'))
  assert '1 water gets no H from placement' in one.not_placed_text()
  summary = reduce_hydrogen.not_placed_summary([], {'free_form': 1})
  assert summary.startswith('Rules, not listed: 1 dictionary H entry on N or '
    'OXT of peptide-like residues skipped by the free-form rule'), summary

def test_free_form():
  # the count is placement's policy, _terminal_h(): one dictionary H entry on
  # N of DNE, whether or not a peptide bond replaces it; DNE alone (its N
  # keeps one H) and DNE as a chain start where chain-start placement puts
  # H1, H2, H3 on N count it too
  for pdb_str, kw in ((dne(), {}), (dne(), {'n_terminal_charge':
      'first_in_chain'}), (dne('1'), {})):
    o, _ = placed(pdb_str, **kw)
    assert o.not_placed_rules['free_form'] == 1, o.not_placed_rules
    assert records(o) == [], records(o)
  o, _ = placed(dne('1'))
  names = [a.name.strip() for a in o.get_model().get_hierarchy().atoms()]
  assert ['H1', 'H2', 'H3'] == [n for n in names if n in ('H1', 'H2', 'H3')], \
    names

def test_free_form_per_conformer():
  # DNE and HT7 in A and B: four sites
  o = place(split_ab(t23.pdb_str))
  assert o.not_placed_rules['free_form'] == 4, o.not_placed_rules
  # DNE with a blank N and the rest in A and B: its H site in A and B; HT7
  # whole (its own blank N next to split HT7 atoms names the H differently,
  # left as it is): one
  o = place(split_ab(t23.pdb_str,
    lambda l: l[17:20] == 'HT7' or l[12:16] == ' N  '))
  assert o.not_placed_rules['free_form'] == 3, o.not_placed_rules

def test_water_rule_reduce2():
  # add_h_to_water=False: the waters that end without H, by the option;
  # add_h_to_water=True: no water line (their H come after the Optimizer)
  description, log = run('water_rule_off', gln_chains(1, waters(3)))
  assert 'Rules, not listed: 3 waters get no H (add_h_to_water=False).' in \
    description, description
  assert '3 waters get no H (add_h_to_water=False)' in log
  description, log = run('water_rule_on', gln_chains(1, waters(3)),
    'add_h_to_water=True')
  assert 'waters get no H' not in description + log
  # supplied H of water 2 remain with keep_existing_H: two waters without H
  description, log = run('water_rule_kept', gln_chains(1, waters(3, (2,))),
    'keep_existing_H=True')
  assert '2 waters get no H (add_h_to_water=False)' in description, \
    description

def test_run_twice():
  # each run() reports its own drops
  o = place(pdb_gln())
  first = records(o)
  o.run()
  assert records(o) == first, records(o)
  o = place(tna.residues(tna.dna, 1, 3), na_terminal_h='chain_ends')
  first, text = records(o), o.na_terminal_gap_text()
  assert records(o, 'na_gap_end') and text
  o.run()
  assert records(o) == first, records(o)
  assert o.na_terminal_gap_text() == text, o.na_terminal_gap_text()

if (__name__ == "__main__"):
  test_missing_neighbour_and_parent()
  test_parent_missing_lys()
  test_parent_in_one_conformer()
  test_water_without_o()
  test_accounting()
  test_supplied_names()
  test_unknown_atoms_supplied()
  test_chain_start_unknown_names()
  test_tertiary_amide()
  test_link()
  test_esterified_o()
  test_link_angle()
  test_link_leaving()
  test_unknown_atoms()
  test_metal_nh_and_amine()
  test_metal_thiolate()
  test_na_gap_end()
  test_letters_differ()
  test_letters_differ_summary()
  test_no_dictionary()
  test_no_bonds()
  test_riding_refused()
  test_not_parameterized()
  test_labels()
  test_reduce2_description_and_log()
  test_no_truncation()
  test_rule_counts()
  test_free_form()
  test_free_form_per_conformer()
  test_water_rule_reduce2()
  test_run_twice()
  print("OK")
