"""
An N terminus of a modified residue whose dictionary does not fit
mod_NH3 (no atom H: 0A0, 004, SAM; or H1-H3 on another atom: 2KZ H2 on CB,
5CA H2 on C2, NPH H3 on C3) keeps the dictionary's N-H names and gets the
missing H under names computed from the dictionary
(mmtbx.monomer_library.terminal_nh3); pdb_interpretation restrains it with
the generated modification NH3-<comp>, also when the output is read again.
Standard residues keep H1-H3 and mod_NH3.
"""
from __future__ import absolute_import, division, print_function
import math
import os
import iotbx.pdb
import mmtbx.model
import mmtbx.monomer_library.server
from libtbx.utils import null_out, Sorry
from scitbx import matrix
from six.moves import cStringIO as StringIO
from iotbx.cli_parser import run_program
from mmtbx.programs import reduce2
from mmtbx.monomer_library import pdb_interpretation, terminal_nh3
from mmtbx.regression import tst_reduce2_not_placed as tnp
from mmtbx.regression import tst_reduce2_chain_start_ligands as tcs

pdb_2kz = """\
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
HETATM    1  OD1 2KZ A   1      17.085  20.247  20.036  1.00 20.00           O
HETATM    2  CG  2KZ A   1      18.411  20.018  20.124  1.00 20.00           C
HETATM    3  OD2 2KZ A   1      18.822  18.892  20.330  1.00 20.00           O
HETATM    4  CB  2KZ A   1      19.356  21.125  19.970  1.00 20.00           C
HETATM    5  CA  2KZ A   1      20.688  20.877  19.924  1.00 20.00           C
HETATM    6  N   2KZ A   1      21.584  21.923  19.934  1.00 20.00           N
HETATM    7  C   2KZ A   1      21.181  19.479  19.862  1.00 20.00           C
HETATM    8  O   2KZ A   1      20.664  18.687  19.102  1.00 20.00           O
HETATM    9  OXT 2KZ A   1      22.202  19.090  20.650  1.00 20.00           O
END
"""

pdb_5ca = """\
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
HETATM    1  N   5CA A   1      17.247  23.037  19.482  1.00 20.00           N
HETATM    2  CA  5CA A   1      16.549  22.018  20.277  1.00 20.00           C
HETATM    3  CB  5CA A   1      15.076  22.407  20.424  1.00 20.00           C
HETATM    4  SG  5CA A   1      14.955  24.013  21.258  1.00 20.00           S
HETATM    5  C   5CA A   1      16.650  20.684  19.583  1.00 20.00           C
HETATM    6  O   5CA A   1      16.907  20.634  18.399  1.00 20.00           O
HETATM    7  N3S 5CA A   1      16.454  19.547  20.279  1.00 20.00           N
HETATM    8  S   5CA A   1      16.565  18.081  19.517  1.00 20.00           S
HETATM    9  O1S 5CA A   1      16.364  17.103  20.529  1.00 20.00           O
HETATM   10  O2S 5CA A   1      15.795  18.204  18.329  1.00 20.00           O
HETATM   11  O5' 5CA A   1      18.007  17.913  19.060  1.00 20.00           O
HETATM   12  C5' 5CA A   1      18.776  17.608  20.225  1.00 20.00           C
HETATM   13  C4' 5CA A   1      20.244  17.426  19.833  1.00 20.00           C
HETATM   14  O4' 5CA A   1      20.788  18.678  19.382  1.00 20.00           O
HETATM   15  C3' 5CA A   1      21.063  16.982  21.060  1.00 20.00           C
HETATM   16  O3' 5CA A   1      21.676  15.712  20.828  1.00 20.00           O
HETATM   17  C2' 5CA A   1      22.142  18.086  21.212  1.00 20.00           C
HETATM   18  O2' 5CA A   1      23.410  17.516  21.546  1.00 20.00           O
HETATM   19  C1' 5CA A   1      22.172  18.708  19.793  1.00 20.00           C
HETATM   20  N9  5CA A   1      22.667  20.086  19.840  1.00 20.00           N
HETATM   21  C8  5CA A   1      21.910  21.210  19.993  1.00 20.00           C
HETATM   22  N7  5CA A   1      22.674  22.264  19.992  1.00 20.00           N
HETATM   23  C5  5CA A   1      23.965  21.884  19.839  1.00 20.00           C
HETATM   24  C6  5CA A   1      25.193  22.564  19.766  1.00 20.00           C
HETATM   25  N6  5CA A   1      25.244  23.944  19.858  1.00 20.00           N
HETATM   26  N1  5CA A   1      26.303  21.851  19.607  1.00 20.00           N
HETATM   27  C2  5CA A   1      26.264  20.535  19.519  1.00 20.00           C
HETATM   28  N3  5CA A   1      25.137  19.856  19.582  1.00 20.00           N
HETATM   29  C4  5CA A   1      23.977  20.482  19.745  1.00 20.00           C
END
"""

pdb_sam = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1
HETATM    1  N   SAM A   1       7.950  14.496  15.920  1.00 20.00           N
HETATM    2  CA  SAM A   1       8.588  15.790  15.642  1.00 20.00           C
HETATM    3  C   SAM A   1       7.601  16.699  14.955  1.00 20.00           C
HETATM    4  O   SAM A   1       6.730  16.227  14.244  1.00 20.00           O
HETATM    5  OXT SAM A   1       7.674  17.906  15.111  1.00 20.00           O
HETATM    6  CB  SAM A   1       9.801  15.576  14.735  1.00 20.00           C
HETATM    7  CG  SAM A   1      10.857  14.756  15.480  1.00 20.00           C
HETATM    8  SD  SAM A   1      12.296  14.503  14.405  1.00 20.00           S
HETATM    9  CE  SAM A   1      11.662  13.181  13.337  1.00 20.00           C
HETATM   10  C5' SAM A   1      13.352  13.532  15.515  1.00 20.00           C
HETATM   11  C4' SAM A   1      14.662  13.191  14.804  1.00 20.00           C
HETATM   12  O4' SAM A   1      15.421  14.390  14.574  1.00 20.00           O
HETATM   13  C3' SAM A   1      15.517  12.268  15.694  1.00 20.00           C
HETATM   14  O3' SAM A   1      15.769  11.023  15.039  1.00 20.00           O
HETATM   15  C2' SAM A   1      16.840  13.055  15.891  1.00 20.00           C
HETATM   16  O2' SAM A   1      17.968  12.181  15.813  1.00 20.00           O
HETATM   17  C1' SAM A   1      16.815  14.036  14.691  1.00 20.00           C
HETATM   18  N9  SAM A   1      17.627  15.222  14.975  1.00 20.00           N
HETATM   19  C8  SAM A   1      17.191  16.380  15.548  1.00 20.00           C
HETATM   20  N7  SAM A   1      18.176  17.224  15.651  1.00 20.00           N
HETATM   21  C5  SAM A   1      19.304  16.663  15.153  1.00 20.00           C
HETATM   22  C6  SAM A   1      20.635  17.087  15.003  1.00 20.00           C
HETATM   23  N6  SAM A   1      21.028  18.347  15.419  1.00 20.00           N
HETATM   24  N1  SAM A   1      21.509  16.249  14.454  1.00 20.00           N
HETATM   25  C2  SAM A   1      21.143  15.044  14.056  1.00 20.00           C
HETATM   26  N3  SAM A   1      19.908  14.605  14.176  1.00 20.00           N
HETATM   27  C4  SAM A   1      18.966  15.370  14.718  1.00 20.00           C
END
"""

# NPH A 117 of a local model, cut to NPH-VAL-MET, renumbered from 1
pdb_nph = """\
CRYST1   80.020   97.150   49.850  90.00  90.00  90.00 C 2 2 21
HETATM 1819  N   NPH A   1      23.870  15.268 -50.490  1.00 25.06           N
HETATM 1820  CA  NPH A   1      23.515  14.664 -49.210  1.00 22.94           C
HETATM 1821  CB  NPH A   1      23.658  15.702 -48.113  1.00 21.44           C
HETATM 1822  SG  NPH A   1      25.281  16.410 -47.839  1.00 31.29           S
HETATM 1823  CD  NPH A   1      25.498  16.310 -46.059  0.50 33.70           C
HETATM 1824  CE  NPH A   1      26.971  16.492 -45.820  0.50 39.92           C
HETATM 1825  OZ  NPH A   1      27.815  15.348 -45.806  0.50 43.26           O
HETATM 1826  NZ  NPH A   1      27.411  17.709 -45.649  0.50 43.87           N
HETATM 1827  C6  NPH A   1      28.830  18.015 -45.436  0.50 51.14           C
HETATM 1828  C5  NPH A   1      29.800  17.009 -45.411  0.50 51.26           C
HETATM 1829  C6A NPH A   1      29.195  19.344 -45.262  0.50 52.54           C
HETATM 1830  C4A NPH A   1      31.137  17.332 -45.219  0.50 52.42           C
HETATM 1831  C10 NPH A   1      30.543  19.684 -45.064  0.50 54.23           C
HETATM 1832  C7  NPH A   1      28.201  20.331 -45.288  0.50 53.28           C
HETATM 1833  C4  NPH A   1      32.087  16.314 -45.202  0.50 53.40           C
HETATM 1834  C1A NPH A   1      31.522  18.669 -45.046  0.50 55.19           C
HETATM 1835  N10 NPH A   1      30.892  21.027 -44.889  0.50 53.87           N
HETATM 1836  C8  NPH A   1      28.549  21.670 -45.111  0.50 54.83           C
HETATM 1837  C3  NPH A   1      33.439  16.599 -45.020  0.50 53.96           C
HETATM 1838  N1  NPH A   1      32.879  18.963 -44.861  0.50 57.78           N
HETATM 1839  C9  NPH A   1      29.889  22.011 -44.911  0.50 55.29           C
HETATM 1840  C2  NPH A   1      33.832  17.924 -44.852  0.50 55.79           C
HETATM 1841  C   NPH A   1      22.061  14.206 -49.245  1.00 23.82           C
HETATM 1842  O   NPH A   1      21.158  15.023 -49.423  1.00 24.62           O
ATOM   1855  N   VAL A   2      21.831  12.913 -49.037  1.00 28.67           N
ATOM   1856  CA  VAL A   2      20.474  12.371 -49.062  1.00 33.96           C
ATOM   1857  C   VAL A   2      20.069  11.729 -47.737  1.00 36.17           C
ATOM   1858  O   VAL A   2      20.854  11.018 -47.109  1.00 40.12           O
ATOM   1859  CB  VAL A   2      20.297  11.336 -50.209  1.00 39.91           C
ATOM   1860  CG1 VAL A   2      18.831  10.915 -50.336  1.00 38.68           C
ATOM   1861  CG2 VAL A   2      20.791  11.918 -51.524  1.00 43.67           C
ATOM   1871  N   MET A   3      18.829  11.986 -47.328  1.00 39.11           N
ATOM   1872  CA  MET A   3      18.273  11.449 -46.092  1.00 39.84           C
ATOM   1873  C   MET A   3      16.761  11.361 -46.256  1.00 45.51           C
ATOM   1874  O   MET A   3      16.066  12.375 -46.173  1.00 45.08           O
ATOM   1875  CB  MET A   3      18.616  12.370 -44.916  1.00 46.77           C
ATOM   1876  CG  MET A   3      19.314  11.688 -43.740  1.00 45.80           C
ATOM   1877  SD  MET A   3      18.225  10.980 -42.492  1.00 41.41           S
ATOM   1878  CE  MET A   3      17.987   9.359 -43.118  1.00 44.09           C
END
"""

# residue 1, peptide-bonded to a GLY; dictionary N-H names; carbon H and
# its carbon; N-CA ideal (mod_NH3's 1.491 for an sp3 CA, the dictionary's
# for 2KZ's sp2 CA)
CONVERTED = [
  ('2KZ', tnp.with_next_gly(pdb_2kz), ['H', 'H4', 'H5'], ('H2', 'CB'), 1.402),
  ('5CA', tnp.with_next_gly(pdb_5ca), ['H', 'HN2', 'HN3'], ('H2', 'C2'),
    1.491),
  ('NPH', pdb_nph, ['H', 'HN3', 'HNT'], ('H3', 'C3'), 1.491),
  ('0A0', tnp.with_next_gly(tnp.pdb_0a0), ['HN1', 'HN2', 'HN3'], None, 1.491),
  ('004', tnp.with_next_gly(tnp.pdb_004), ['HN', 'HNA', 'HNB'], None, 1.491),
  ('SAM', tnp.with_next_gly(pdb_sam), ['HN1', 'HN2', 'HN3'], None, 1.491)]

_runs = {}

def run(prefix, pdb_str, *args):
  """reduce2 on pdb_str (once per input and arguments): the output
  hierarchy, the description, the output file name."""
  key = (pdb_str, args)
  if key not in _runs:
    prefix = 'tst_reduce2_terminal_nh3_' + prefix
    fn, out, txt = prefix + '.pdb', prefix + '_out.pdb', prefix + '.txt'
    with open(fn, 'w') as f:
      f.write(pdb_str)
    for p in (out, txt):
      if os.path.exists(p): os.remove(p)
    run_program(program_class=reduce2.Program, logger=StringIO(), args=[fn,
      'output.filename=' + out, 'output.description_file_name=' + txt,
      'overwrite=True'] + list(args))
    with open(txt) as f:
      description = f.read()
    _runs[key] = (iotbx.pdb.input(out).construct_hierarchy(), description,
      out)
  return _runs[key]

_mon_lib_srv = []
def reread(fn):
  """The output processed as phenix.refine reads it: the model, the
  residue names with their modifications."""
  m = mmtbx.model.manager(model_input=iotbx.pdb.input(fn), log=null_out())
  m.process(make_restraints=True)
  if not _mon_lib_srv:
    _mon_lib_srv.append(mmtbx.monomer_library.server.server())
    _mon_lib_srv.append(mmtbx.monomer_library.server.ener_lib())
  p = pdb_interpretation.process(mon_lib_srv=_mon_lib_srv[0],
    ener_lib=_mon_lib_srv[1], log=null_out(), pdb_inp=iotbx.pdb.input(fn),
    substitute_non_crystallographic_unit_cell_if_necessary=True)
  names = set(mm.residue_name for mm in
    p.all_chain_proxies.all_monomer_mappings)
  return m, names

def partners(m):
  """{i_seq: [bonded i_seq]} from the bond restraints"""
  result = {}
  for prx in m.get_restraints_manager().geometry.pair_proxies(
      ).bond_proxies.simple:
    i, j = prx.i_seqs
    result.setdefault(i, []).append(j)
    result.setdefault(j, []).append(i)
  return result

def is_hd(a):
  return a.element.strip().upper() in ('H', 'D')

def first_rg(h):
  return h.models()[0].chains()[0].residue_groups()[0]

def h_on_n(h, altloc=''):
  """names of the H/D within 1.3 A of the N of residue 1, conformer altloc
  (and the blank group)"""
  rg = first_rg(h)
  ns = [a for a in rg.atoms() if a.name.strip() == 'N' and
        a.parent().altloc in ('', altloc)]
  return sorted(a.name.strip() for ag in rg.atom_groups()
    if ag.altloc == altloc for a in ag.atoms()
    if is_hd(a) and [n for n in ns if a.distance(n) < 1.3])

def worst_angle_at_n(m):
  """largest |model - ideal| of the angles at N of residue 1 with an H"""
  atoms = m.get_hierarchy().atoms()
  worst = 0
  for prx in m.get_restraints_manager().geometry.angle_proxies:
    x, y, z = [atoms[i] for i in prx.i_seqs]
    if y.name.strip() != 'N' or y.parent().parent().resseq.strip() != '1':
      continue
    if not (is_hd(x) or is_hd(z)): continue
    d = math.degrees((matrix.col(x.xyz) - matrix.col(y.xyz)).angle(
      matrix.col(z.xyz) - matrix.col(y.xyz)))
    worst = max(worst, abs(d - prx.angle_ideal))
  return worst

def unrestrained_h(m):
  bonded = partners(m)
  return [a.id_str() for a in first_rg(m.get_hierarchy()).atoms()
    if is_hd(a) and a.i_seq not in bonded]

def test_converted_termini():
  # names on N, the carbon H on its carbon, every H restrained in the
  # output read again, with NH3-<comp>; N-CA as mod_NH3 for an sp3 CA
  for comp, pdb_str, names, carbon_h, n_ca in CONVERTED:
    h, _, out = run(comp, pdb_str)
    assert h_on_n(h) == names, (comp, h_on_n(h))
    if carbon_h is not None:
      ag = first_rg(h).atom_groups()[0]
      heavy = [a for a in ag.atoms() if not is_hd(a)]
      x = ag.get_atom(carbon_h[0])
      near = sorted(heavy, key=lambda a: a.distance(x))[0]
      assert near.name.strip() == carbon_h[1] and near.distance(x) < 1.2, \
        (comp, near.name, near.distance(x))
    m, residue_names = reread(out)
    assert '%s%%NH3-%s' % (comp, comp) in residue_names, (comp,
      residue_names)
    assert unrestrained_h(m) == [], (comp, unrestrained_h(m))
    atoms = m.get_hierarchy().atoms()
    ideal = [prx.distance_ideal for prx in m.get_restraints_manager(
      ).geometry.pair_proxies().bond_proxies.simple
      if sorted(atoms[i].name.strip() for i in prx.i_seqs) == ['CA', 'N']
      and atoms[prx.i_seqs[0]].parent().resname.strip() == comp]
    assert len(ideal) == 1 and abs(ideal[0] - n_ca) < 1.e-3, (comp, ideal)
    assert worst_angle_at_n(m) < 1.0, (comp, worst_angle_at_n(m))

def test_standard_terminus():
  # ALA: H1-H3 and mod_NH3, as before
  h, _, out = run('ALA', tcs.pdb_ala_gly)
  assert h_on_n(h) == ['H1', 'H2', 'H3'], h_on_n(h)
  m, residue_names = reread(out)
  assert 'ALA%NH3' in residue_names, residue_names
  assert unrestrained_h(m) == []

def turned_b(pdb_str):
  """blank N, every other atom of residue 1 in A and B (B turned 35 deg
  about N-C)"""
  lines = pdb_str.splitlines()
  xyz = lambda l: matrix.col((float(l[30:38]), float(l[38:46]),
    float(l[46:54])))
  first = [l for l in lines if l[:6] in ('ATOM  ', 'HETATM') and
           l[22:26].strip() == '1']
  n = xyz([l for l in first if l[12:16] == ' N  '][0])
  c = xyz([l for l in first if l[12:16] == ' C  '][0])
  rot = (c - n).axis_and_angle_as_r3_rotation_matrix(35, deg=True)
  result = []
  for l in lines:
    if l not in first or l[12:16] == ' N  ':
      result.append(l)
      continue
    result.append(l[:16] + 'A' + l[17:54] + '  0.50' + l[60:])
    result.append(l[:16] + 'B' + l[17:30] + '%8.3f%8.3f%8.3f' %
      tuple(rot * (xyz(l) - n) + n) + '  0.50' + l[60:])
  return '\n'.join(result) + '\n'

def with_ala(pdb_str, own):
  """residue 1 as conformer own, ALA (N CA C O CB) as the other"""
  other = 'B' if own == 'A' else 'A'
  result = []
  for l in pdb_str.splitlines():
    if not (l[:6] in ('ATOM  ', 'HETATM') and l[22:26].strip() == '1'):
      result.append(l)
      continue
    mine = [l[:16] + own + l[17:54] + '  0.50' + l[60:]]
    ala = []
    if l[12:16].strip() in ('N', 'CA', 'C', 'O', 'CB'):
      ala = [l[:16] + other + 'ALA' + l[20:54] + '  0.50' + l[60:]]
    result += mine + ala if own == 'A' else ala + mine
  return '\n'.join(result) + '\n'

def test_split_ca():
  # a blank N with CA in A and B: all three N-H in each conformer, none in
  # the blank group, angles at N at their ideal values
  for comp, pdb_str, names, _, _ in CONVERTED[:1] + CONVERTED[3:4]:
    h, _, out = run(comp + '_split', turned_b(pdb_str))
    assert h_on_n(h) == [], (comp, h_on_n(h))
    for alt in 'AB':
      assert h_on_n(h, alt) == names, (comp, alt, h_on_n(h, alt))
    m, _ = reread(out)
    assert unrestrained_h(m) == [], (comp, unrestrained_h(m))
    assert worst_angle_at_n(m) < 1.0, (comp, worst_angle_at_n(m))
    # the three N-H supplied in the blank group (keep_existing_H): they
    # follow the split CA too, computed names included
    h = h.deep_copy()
    rg = first_rg(h)
    blank = [ag for ag in rg.atom_groups() if ag.altloc == ''][0]
    for ag in rg.atom_groups():
      for name in names:
        a = ag.get_atom(name)
        if a is None: continue
        if ag.altloc == 'A':
          x = a.detached_copy()
          x.set_occ(1.0)
          blank.append_atom(x)
        ag.remove_atom(a)
    h, _, out = run(comp + '_split_blank_h', h.as_pdb_string(
      crystal_symmetry=iotbx.pdb.input(source_info=None,
      lines=pdb_str).crystal_symmetry()), 'keep_existing_H=True')
    assert h_on_n(h) == [], (comp, h_on_n(h))
    for alt in 'AB':
      assert h_on_n(h, alt) == names, (comp, alt, h_on_n(h, alt))
    m, _ = reread(out)
    assert unrestrained_h(m) == [], (comp, unrestrained_h(m))
    assert worst_angle_at_n(m) < 1.0, (comp, worst_angle_at_n(m))

def test_mixed_alternates():
  # 0A0 or 2KZ with ALA as the other conformer, in both orders: each by its
  # own dictionary (ALA H1-H3, as for a standard residue)
  for comp, pdb_str, names, _, _ in CONVERTED[:1] + CONVERTED[3:4]:
    for own in 'AB':
      h, _, out = run('%s_ala_%s' % (comp, own), with_ala(pdb_str, own))
      other = 'B' if own == 'A' else 'A'
      assert h_on_n(h, own) == names, (comp, own, h_on_n(h, own))
      assert h_on_n(h, other) == ['H1', 'H2', 'H3'], (comp, own,
        h_on_n(h, other))
      m, _ = reread(out)
      assert unrestrained_h(m) == [], (comp, own, unrestrained_h(m))
      assert worst_angle_at_n(m) < 1.0, (comp, own, worst_angle_at_n(m))

def supplied(comp):
  """the converted output with the N-H named H1, H2, H3 (as master wrote
  2KZ: H2 on N, none on CB)"""
  h, _, _ = run(comp, [c[1] for c in CONVERTED if c[0] == comp][0])
  h = h.deep_copy()
  ag = first_rg(h).atom_groups()[0]
  names = h_on_n(h)
  for a in list(ag.atoms()):
    if is_hd(a) and a.name.strip() not in names:
      if a.name.strip() == 'H2': ag.remove_atom(a)
  for old, new in zip(names, ('H1', 'H2', 'H3')):
    ag.get_atom(old).set_name((' ' + new).ljust(4))
  return h.as_pdb_string(crystal_symmetry=iotbx.pdb.input(source_info=None,
    lines=[c[1] for c in CONVERTED if c[0] == comp][0]).crystal_symmetry())

def test_supplied_h1_h3():
  # keep_existing_H: supplied H1-H3 the dictionary lacks are never parented
  # by name: listed as unknown, not renamed, with the message; the terminus
  # is built with the dictionary's names. 2KZ's supplied H2 (CB's name) is
  # the name-collision line; H2 goes back to CB
  for comp, unknown, names in (('0A0', 'H1, H2, H3', 'HN1, HN2, HN3'),
                               ('2KZ', 'H1, H3', 'H, H4, H5')):
    h, description, _ = run(comp + '_supplied', supplied(comp),
      'keep_existing_H=True')
    line = ('%s %s A 1: on N, not in the %s dictionary (its terminal N-H: '
      '%s); listed as unknown atoms, not renamed; reduce2 without '
      'keep_existing_H rebuilds them' % (unknown, comp, comp, names))
    assert line in [l.strip() for l in description.splitlines()], \
      description
    assert h_on_n(h) == names.split(', '), (comp, h_on_n(h))
  assert ('H2 2KZ A 1: within bond length of N, but the 2KZ dictionary '
    'bonds H2 to CB; not renamed; name collision unresolved' in
    [l.strip() for l in description.splitlines()]), description
  # a free 0A0 (no next residue: not converted) with a supplied H1: no
  # promise of a terminus it does not get
  h, description, _ = run('0A0_free_supplied', tnp.pdb_0a0_supplied_h1(),
    'keep_existing_H=True')
  assert 'on N, not in the 0A0 dictionary' not in description, description
  assert 'terminal conversion not applied' in description, description

def test_master_written_file():
  # 2KZ with H1-H3 on N read as is: mod_NH3 is not applied (it would take
  # CB's H2), H1 and H3 are unknown atoms; never the apply_mod() error
  mon_lib_srv = mmtbx.monomer_library.server.server()
  nh3 = mon_lib_srv.mod_mod_id_dict['NH3']
  assert terminal_nh3.mod_takes_name(
    mon_lib_srv.get_comp_comp_id_direct('2KZ'), nh3)
  assert not terminal_nh3.mod_takes_name(
    mon_lib_srv.get_comp_comp_id_direct('SCY'), nh3)
  m = mmtbx.model.manager(model_input=iotbx.pdb.input(source_info=None,
    lines=supplied('2KZ')), log=null_out())
  try:
    m.process(make_restraints=True)
  except Sorry as e:
    assert 'apply_mod' not in str(e), str(e)
    # H1 and H3
    assert 'unknown nonbonded energy type symbols: 2' in str(e), str(e)
  else:
    raise AssertionError('no stop on unknown H1/H3')

NEUTRON = ('all_hd', 'hd_and_h', 'hd_and_d', 'all_d')

def kept_pair(od1=None):
  """2KZ with a blank N, CA etc. in A and B (0.5/0.5) and a supplied H/D
  pair on N (0.3/0.7, B 21/31); od1: OD1's occupancy in both conformers"""
  lines = turned_b(tnp.with_next_gly(pdb_2kz)).splitlines()
  n = [l for l in lines if l[12:16] == ' N  ' and '2KZ' in l][0]
  h, _, _ = run('2KZ_split', turned_b(tnp.with_next_gly(pdb_2kz)))
  # the site: the conformer A H of the run without it
  x = [a for a in first_rg(h).atoms() if a.name.strip() == 'H' and
       a.parent().altloc == 'A'][0].xyz
  if od1 is not None:
    lines = [l[:54] + '%6.2f' % od1 + l[60:] if l[12:16] == ' OD1' else l
      for l in lines]
  i = lines.index(n) + 1
  pair = [
    'HETATM   16  H  C2KZ A   1    %8.3f%8.3f%8.3f  0.30 21.00           H'
      % x,
    'HETATM   17  D  D2KZ A   1    %8.3f%8.3f%8.3f  0.70 31.00           D'
      % x]
  return '\n'.join(lines[:i] + pair + lines[i:]) + '\n'

def test_neutron():
  # 0A0: each option read again with every H/D bonded, the three N-H sites
  # under the dictionary's names. A supplied H/D pair on a blank N under a
  # split CA: one isotope per conformer at its occupancy as placed (0.5 +
  # 0.5 = the input 1.0; with OD1 at 0.4 the conformer has no single value:
  # 0.49, as without the pair), the supplied B, every record used
  for option in NEUTRON:
    h, description, out = run('0A0_' + option,
      tnp.with_next_gly(tnp.pdb_0a0), 'neutron_option=' + option)
    rg = first_rg(h)
    n = [a for a in rg.atoms() if a.name.strip() == 'N'][0]
    on_n = [a.name.strip() for a in rg.atoms() if is_hd(a) and
      a.distance(n) < 1.3]
    assert set('H' + x[1:] for x in on_n) == set(['HN1', 'HN2', 'HN3']), \
      (option, on_n)
    m, _ = reread(out)
    assert unrestrained_h(m) == [], (option, unrestrained_h(m))
    h, description, out = run('pair_' + option, kept_pair(),
      'keep_existing_H=True', 'neutron_option=' + option)
    site = [a for a in first_rg(h).atoms() if a.name.strip() in ('H', 'D')]
    assert sorted((a.parent().altloc, a.occ) for a in site) == [
      ('A', 0.5), ('B', 0.5)], (option, [(a.id_str(), a.occ) for a in site])
    assert [a.b for a in site if a.name.strip() == 'D'] in ([], [31.0]*2), \
      (option, [(a.id_str(), a.b) for a in site])
    assert 'Input values not used' not in description, description
    m, _ = reread(out)
    assert unrestrained_h(m) == [], (option, unrestrained_h(m))
    h, description, out = run('pair_od1_' + option, kept_pair(od1=0.4),
      'keep_existing_H=True', 'neutron_option=' + option)
    site = [a for a in first_rg(h).atoms() if a.name.strip() in ('H', 'D')]
    assert sorted((a.parent().altloc, a.occ) for a in site) == [
      ('A', 0.49), ('B', 0.49)], (option, [(a.id_str(), a.occ) for a in site])
    assert 'Input values not used' not in description, description

def test_metal_bound_start():
  # 0A0 at a chain start with its N 2.1 A from a Zn (default
  # metal_bound_amines_as_nh2): NH2, listed as for a standard residue
  lines = tnp.with_next_gly(tnp.pdb_0a0).splitlines()
  xyz = lambda name: matrix.col([float(x) for x in [l for l in lines
    if l[12:16] == name and '0A0' in l][0][30:54].split()])
  n, ca, c = xyz(' N  '), xyz(' CA '), xyz(' C  ')
  a = (ca - n).normalize()
  q = (c - ca)
  q = (q - a * q.dot(a)).normalize() * -1
  t = math.radians(109.5)
  zn = n + (a * math.cos(t) + q * math.sin(t)) * 2.1
  i = max(k for k, l in enumerate(lines) if l[:6] in ('ATOM  ', 'HETATM'))
  lines.insert(i + 1, 'HETATM   99 ZN    ZN B   1    %8.3f%8.3f%8.3f  1.00 '
    '20.00          ZN' % tuple(zn))
  h, description, _ = run('0A0_zn', '\n'.join(lines) + '\n')
  assert h_on_n(h) == ['HN2', 'HN3'], h_on_n(h)
  assert 'N 0A0 A 1: metal-bound chain start, placed as NH2' in \
    [l.strip() for l in description.splitlines()], description

if __name__ == '__main__':
  test_converted_termini()
  test_standard_terminus()
  test_split_ca()
  test_mixed_alternates()
  test_supplied_h1_h3()
  test_master_written_file()
  test_neutron()
  test_metal_bound_start()
  print('OK')
