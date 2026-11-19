"""validate_ligands with reduce2's flips, symmetry copies and metal-bound
amines as NH2 (add_flip_movers, include_symmetric_copies,
metal_bound_amines_as_nh2; all True by default), and their report in the log.

- 1avd fragment (Asn A 12 and the residues within 5 A of its side chain, with
  BTN A 400): the Asn A 12 amide flips next to the biotin with the default
  and not with add_flip_movers=False; listed and marked in the log.
- Ethanol EOH A 1 and Asp B 1, P 1 2 1: the EOH O is 2.65 A from a
  carboxylate O of the Asp mate -x+1,y,-z+1. With the symmetry copies its HO
  H-bonds to that O; without (include_symmetric_copies=False), or with map
  input, it does not.
- 7T5E His A 1 on Cu A 301 with its glycan: the N-terminal amine is placed as
  NH2 with the default and as NH3+ with metal_bound_amines_as_nh2=False.
- Flips are read from the final model: a blank His flipped only in the first
  conformer pass is not listed; split conformers once each, labelled; with
  add_flip_movers=False a His whose Zn binds only the flipped ring is not
  flipped.
- Map input: no metal contact across the map box (a Cu one box edge away).
- 1avd with its data and the defaults: the log and a few values pinned.
"""
from __future__ import absolute_import, division, print_function
import os
from six.moves import cStringIO as StringIO
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out
from iotbx.cli_parser import run_program
from mmtbx.programs import validate_ligands as val_lig
from mmtbx.regression.tst_reduce2_symmetry_copies import nearest_mate
from mmtbx.regression.reduce2_metal_amine_utils import pdb_his_cu_glycan
from mmtbx.regression import tst_reduce_switch as trs

from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL) # Only show critical errors

# 1avd: Asn A 12 and the residues within 5 A of its side chain, BTN A 400
pdb_1avd_asn12 = """\
CRYST1   80.150   80.150   85.280  90.00  90.00  90.00 P 42 21 2
ATOM     64  N   THR A  11      36.411  22.423  25.750  1.00 21.35           N
ATOM     65  CA  THR A  11      37.347  21.971  24.719  1.00 32.99           C
ATOM     66  C   THR A  11      38.326  20.880  25.194  1.00 42.31           C
ATOM     67  O   THR A  11      38.787  20.896  26.345  1.00 55.37           O
ATOM     68  CB  THR A  11      38.141  23.157  24.137  1.00 42.06           C
ATOM     69  OG1 THR A  11      37.327  23.723  23.143  1.00 63.77           O
ATOM     70  CG2 THR A  11      39.463  22.723  23.500  1.00 48.19           C
ATOM     71  N   ASN A  12      38.706  19.997  24.247  1.00 29.82           N
ATOM     72  CA  ASN A  12      39.648  18.873  24.392  1.00 20.35           C
ATOM     73  C   ASN A  12      41.065  19.115  23.959  1.00 23.66           C
ATOM     74  O   ASN A  12      41.363  19.999  23.132  1.00 28.59           O
ATOM     75  CB  ASN A  12      39.286  17.800  23.347  1.00 24.45           C
ATOM     76  CG  ASN A  12      38.324  16.943  24.008  1.00 37.16           C
ATOM     77  OD1 ASN A  12      37.894  17.367  25.091  1.00 35.61           O
ATOM     78  ND2 ASN A  12      38.107  15.746  23.474  1.00 26.09           N
ATOM     79  N   ASP A  13      41.907  18.165  24.408  1.00 10.60           N
ATOM     80  CA  ASP A  13      43.276  18.088  23.985  1.00 12.67           C
ATOM     81  C   ASP A  13      43.351  17.565  22.499  1.00 29.38           C
ATOM     82  O   ASP A  13      44.421  17.380  21.876  1.00 34.89           O
ATOM     83  CB  ASP A  13      44.005  17.136  24.934  1.00 15.79           C
ATOM     84  CG  ASP A  13      43.311  15.817  25.034  1.00  6.48           C
ATOM     85  OD1 ASP A  13      42.120  15.743  25.231  1.00 16.24           O
ATOM     86  OD2 ASP A  13      44.121  14.793  24.959  1.00 21.34           O
ATOM     87  N   LEU A  14      42.195  17.247  21.936  1.00 31.59           N
ATOM     88  CA  LEU A  14      42.146  16.835  20.551  1.00 42.19           C
ATOM     89  C   LEU A  14      41.603  18.010  19.720  1.00 36.47           C
ATOM     90  O   LEU A  14      41.780  18.066  18.513  1.00 37.31           O
ATOM     91  CB  LEU A  14      41.301  15.576  20.253  1.00 36.87           C
ATOM     92  CG  LEU A  14      41.754  14.252  20.833  1.00 13.98           C
ATOM     93  CD1 LEU A  14      40.755  13.241  20.332  1.00 15.46           C
ATOM     94  CD2 LEU A  14      43.144  13.764  20.386  1.00 10.28           C
ATOM     95  N   GLY A  15      40.869  18.907  20.381  1.00 24.84           N
ATOM     96  CA  GLY A  15      40.308  20.083  19.738  1.00  5.23           C
ATOM     97  C   GLY A  15      38.812  20.070  19.666  1.00 15.74           C
ATOM     98  O   GLY A  15      38.244  20.991  19.132  1.00 33.27           O
ATOM     99  N   SER A  16      38.183  19.040  20.220  1.00 21.52           N
ATOM    100  CA  SER A  16      36.722  18.900  20.236  1.00 22.88           C
ATOM    101  C   SER A  16      36.000  19.960  21.048  1.00 22.99           C
ATOM    102  O   SER A  16      36.568  20.578  21.932  1.00 43.38           O
ATOM    103  CB  SER A  16      36.273  17.493  20.604  1.00 42.89           C
ATOM    104  OG  SER A  16      36.674  16.542  19.609  1.00 23.04           O
ATOM    219  N   TYR A  33      28.188  17.767  20.747  1.00 23.41           N
ATOM    220  CA  TYR A  33      29.593  17.508  20.582  1.00 33.40           C
ATOM    221  C   TYR A  33      30.017  17.825  19.206  1.00 27.52           C
ATOM    222  O   TYR A  33      29.236  17.568  18.306  1.00 26.44           O
ATOM    223  CB  TYR A  33      29.937  16.028  20.866  1.00 40.37           C
ATOM    224  CG  TYR A  33      31.355  15.789  21.344  1.00 20.90           C
ATOM    225  CD1 TYR A  33      31.642  15.815  22.704  1.00 12.90           C
ATOM    226  CD2 TYR A  33      32.401  15.536  20.451  1.00 18.25           C
ATOM    227  CE1 TYR A  33      32.919  15.514  23.173  1.00  7.61           C
ATOM    228  CE2 TYR A  33      33.700  15.295  20.899  1.00  6.81           C
ATOM    229  CZ  TYR A  33      33.956  15.329  22.271  1.00 12.27           C
ATOM    230  OH  TYR A  33      35.210  15.145  22.777  1.00 27.87           O
ATOM    894  N   ASN A 118      40.641  12.072  27.338  1.00  4.64           N
ATOM    895  CA  ASN A 118      40.350  13.427  26.970  1.00 16.67           C
ATOM    896  C   ASN A 118      40.309  14.392  28.142  1.00 28.56           C
ATOM    897  O   ASN A 118      39.731  14.170  29.182  1.00 31.95           O
ATOM    898  CB  ASN A 118      39.048  13.575  26.200  1.00 22.48           C
ATOM    899  CG  ASN A 118      39.271  13.495  24.720  1.00 35.04           C
ATOM    900  OD1 ASN A 118      38.553  12.753  24.079  1.00 20.67           O
ATOM    901  ND2 ASN A 118      40.351  14.097  24.213  1.00  4.00           N
ATOM    902  N   ILE A 119      40.847  15.543  27.924  1.00 36.37           N
ATOM    903  CA  ILE A 119      40.812  16.561  28.930  1.00 29.02           C
ATOM    904  C   ILE A 119      39.941  17.724  28.495  1.00 22.86           C
ATOM    905  O   ILE A 119      40.052  18.303  27.416  1.00 28.10           O
ATOM    906  CB  ILE A 119      42.247  16.961  29.262  1.00 22.81           C
ATOM    907  CG1 ILE A 119      42.986  15.656  29.439  1.00 23.78           C
ATOM    908  CG2 ILE A 119      42.361  17.727  30.578  1.00 16.63           C
ATOM    909  CD1 ILE A 119      44.383  15.827  30.019  1.00 13.11           C
ATOM    910  N   PHE A 120      39.027  18.087  29.314  1.00 21.80           N
ATOM    911  CA  PHE A 120      38.293  19.241  28.887  1.00 28.97           C
ATOM    912  C   PHE A 120      38.579  20.407  29.777  1.00 39.83           C
ATOM    913  O   PHE A 120      38.847  20.298  30.971  1.00 34.58           O
ATOM    914  CB  PHE A 120      36.785  19.083  28.693  1.00 24.13           C
ATOM    915  CG  PHE A 120      36.444  17.877  27.889  1.00 29.81           C
ATOM    916  CD1 PHE A 120      36.887  16.601  28.255  1.00 24.01           C
ATOM    917  CD2 PHE A 120      35.720  18.050  26.705  1.00 31.22           C
ATOM    918  CE1 PHE A 120      36.560  15.501  27.456  1.00 20.73           C
ATOM    919  CE2 PHE A 120      35.364  16.955  25.915  1.00 22.78           C
ATOM    920  CZ  PHE A 120      35.819  15.681  26.281  1.00 12.91           C
TER
HETATM 1954  C10 BTN A 400      33.064   8.812  15.357  1.00  5.67           C
HETATM 1955  C11 BTN A 400      33.347   7.551  14.580  1.00 90.00           C
HETATM 1956  C2  BTN A 400      36.012  10.285  18.933  1.00 12.79           C
HETATM 1957  C3  BTN A 400      37.005  12.897  20.510  1.00 60.52           C
HETATM 1958  C4  BTN A 400      37.338  10.937  19.271  1.00 45.96           C
HETATM 1959  C5  BTN A 400      37.595  10.645  20.766  1.00 31.24           C
HETATM 1960  C6  BTN A 400      36.511   9.733  21.328  1.00  4.00           C
HETATM 1961  C7  BTN A 400      35.408  11.018  17.703  1.00 40.70           C
HETATM 1962  C8  BTN A 400      34.100  10.397  17.155  1.00 13.99           C
HETATM 1963  C9  BTN A 400      34.337   9.270  16.140  1.00 24.45           C
HETATM 1964  N1  BTN A 400      37.432  11.964  21.340  1.00 10.17           N
HETATM 1965  N2  BTN A 400      36.917  12.344  19.280  1.00 35.54           N
HETATM 1966  O11 BTN A 400      32.463   6.780  14.309  1.00 26.97           O
HETATM 1967  O12 BTN A 400      34.631   7.358  14.257  1.00 22.93           O
HETATM 1968  O3  BTN A 400      36.748  14.066  20.803  1.00 27.90           O
HETATM 1969  S1  BTN A 400      35.049  10.359  20.469  1.00 46.22           S
END
"""

# ethanol in place of Ser A 1 (C1 at CB, O at OG) of the Ser + Asp pair of
# tst_reduce2_symmetry_copies; C2 is Ser CA turned 120 deg about C1-O, so
# that HO anti to C2, its orientation without an acceptor, points away from
# the Asp mate
pdb_eoh_asp = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1 2 1
HETATM    1  C2  EOH A   1      23.236  27.462  15.276  1.00 20.00           C
HETATM    2  C1  EOH A   1      22.432  26.665  14.248  1.00 20.00           C
HETATM    3  O   EOH A   1      22.850  25.299  14.272  1.00 20.00           O
ATOM      4  N   ASP B   1       5.600  20.940  11.664  1.00 20.00           N
ATOM      5  CA  ASP B   1       7.082  21.049  11.550  1.00 20.00           C
ATOM      6  C   ASP B   1       7.764  20.210  12.632  1.00 20.00           C
ATOM      7  O   ASP B   1       8.982  20.391  12.840  1.00 20.00           O
ATOM      8  CB  ASP B   1       7.499  22.520  11.683  1.00 20.00           C
ATOM      9  CG  ASP B   1       7.251  23.076  13.088  1.00 20.00           C
ATOM     10  OD1 ASP B   1       7.554  24.272  13.319  1.00 20.00           O
ATOM     11  OD2 ASP B   1       6.766  22.312  13.954  1.00 20.00           O
ATOM     12  OXT ASP B   1       7.061  19.382  13.265  1.00 20.00           O
END
"""

OFF = ['add_flip_movers=False', 'include_symmetric_copies=False',
  'metal_bound_amines_as_nh2=False']

# ------------------------------------------------------------------------------

def run_vl(prefix, pdb_str, args=(), map_from=None):
  """validate_ligands on pdb_str; returns (result, log text, working model
  hierarchy). map_from: also an Fcalc map of that model text."""
  fn = prefix + '.pdb'
  with open(fn, 'w') as fh:
    fh.write(pdb_str)
  files = [fn]
  if map_from is not None:
    from mmtbx.regression.tst_validate_ligands_4 import _simulated_map_manager
    inp = iotbx.pdb.input(lines=map_from.split('\n'), source_info=None)
    model = mmtbx.model.manager(model_input=inp, log=null_out())
    model.process(make_restraints=True)
    model.setup_scattering_dictionaries(scattering_table='electron')
    map_fn = prefix + '.map'
    _simulated_map_manager(model, d_min=2.0).write_map(map_fn)
    files.append(map_fn)
  log = StringIO()
  result = run_program(program_class=val_lig.Program,
    args=files + list(args), logger=log)
  return result, log.getvalue(), result.ligand_manager.model.get_hierarchy()

def section(log):
  """The lines of the log section 'Changes by reduce2 besides adding H'."""
  lines = log.splitlines()
  i = [k for k, l in enumerate(lines) if 'Changes by reduce2 besides' in l]
  assert len(i) == 1, log
  out = []
  for l in lines[i[0] + 1:]:
    if '----------' in l: break
    if l.strip(): out.append(l.rstrip())
  return out

def atom(h, chain, resseq, name):
  r = [a for a in h.atoms() if a.parent().parent().parent().id == chain and
    a.parent().parent().resseq.strip() == resseq and a.name.strip() == name]
  assert len(r) == 1, (chain, resseq, name, len(r))
  return r[0]

# ------------------------------------------------------------------------------

def test_asn_flip():
  """The Asn A 12 amide flips next to the biotin with the default, listed and
  marked; not with add_flip_movers=False."""
  inp = iotbx.pdb.input(lines=pdb_1avd_asn12.split('\n'), source_info=None)
  od1 = atom(inp.construct_hierarchy(), 'A', '12', 'OD1').xyz
  _, log, h = run_vl('tst_validate_ligands_5_flip', pdb_1avd_asn12)
  lines = section(log)
  assert 'Flipped Asn/Gln/His: 1' in lines, lines
  assert '  ASN A 12 (amide)  * within 5 A of ligand: BTN A 400 3.5 A' in \
    lines, lines
  # the validated model carries the flipped heavy atoms: ND2 where OD1 was
  # (within 0.1 A here: the flip rotates about CB-CG)
  assert approx_xyz(atom(h, 'A', '12', 'ND2').xyz, od1), od1
  _, log, h = run_vl('tst_validate_ligands_5_noflip', pdb_1avd_asn12,
    ['add_flip_movers=False'])
  lines = section(log)
  assert 'Flips: off (add_flip_movers=False), none made' in lines, lines
  assert not [l for l in lines if 'ASN A 12' in l], lines
  assert approx_xyz(atom(h, 'A', '12', 'OD1').xyz, od1), od1
  print('test_asn_flip OK')

def approx_xyz(a, b, eps=0.3):
  return sum((x - y) ** 2 for x, y in zip(a, b)) ** 0.5 < eps

def ho_to_asp(h):
  """EOH HO to the nearest mate of an Asp carboxylate O."""
  cs = iotbx.pdb.input(lines=pdb_eoh_asp.split('\n'),
    source_info=None).crystal_symmetry()
  ho = atom(h, 'A', '1', 'HO').xyz
  return min(nearest_mate(cs, ho, atom(h, 'B', '1', n).xyz)
    for n in ('OD1', 'OD2'))

def ho_torsion(h):
  from cctbx import geometry_restraints
  sites = [atom(h, 'A', '1', n).xyz for n in ('C2', 'C1', 'O', 'HO')]
  return geometry_restraints.dihedral(sites=sites, angle_ideal=0,
    weight=1).angle_model

def sym_hbonds(result):
  """The H-bond entries of the EOH profile across symmetry."""
  lig = result.ligand_manager.interactions_as_dict()['ligands'][0]
  return [(e['labels'][1:], e['symop']) for e in lig['entries']
    if e['type'] == 'hbond' and e['symop']]

def test_symmetry_contact():
  """The EOH HO H-bonds to the Asp mate with the copies (default), also in
  the interaction profile; not with include_symmetric_copies=False, and not
  with map input, where no copies are used: there the HO keeps the
  orientation it has without copies."""
  r, log, h = run_vl('tst_validate_ligands_5_copies', pdb_eoh_asp,
    ['interactions=True'])
  lines = section(log)
  assert 'Symmetry copies (crystal contacts): used' in lines, lines
  assert [l for l in lines if l.startswith(
    '  Symmetry copies (fixed environment): model 1,')], lines
  assert ho_to_asp(h) < 2.0, ho_to_asp(h)
  assert sym_hbonds(r) == [(['A EOH 1 HO', 'B ASP 1 OD2 (-x+1,y,-z+1)'],
    '-x+1,y,-z+1')], sym_hbonds(r)
  r, log, h_off = run_vl('tst_validate_ligands_5_nocopies', pdb_eoh_asp,
    ['include_symmetric_copies=False', 'interactions=True'])
  lines = section(log)
  assert ('Symmetry copies (crystal contacts): off '
    '(include_symmetric_copies=False)') in lines, lines
  assert ho_to_asp(h_off) > 2.5, ho_to_asp(h_off)
  assert sym_hbonds(r) == [], sym_hbonds(r)
  _, log, h_map = run_vl('tst_validate_ligands_5_map', pdb_eoh_asp,
    map_from=pdb_eoh_asp)
  lines = section(log)
  assert ('Symmetry copies (crystal contacts): none: map input (the map box '
    'is not a lattice)') in lines, lines
  assert abs(ho_torsion(h_map) - ho_torsion(h_off)) < 1, (
    ho_torsion(h_map), ho_torsion(h_off))
  assert abs(ho_torsion(h) - ho_torsion(h_off)) > 10, (
    ho_torsion(h), ho_torsion(h_off))
  print('test_symmetry_contact OK')

def test_metal_amine():
  """His A 1 on Cu A 301: NH2 with the default, listed; NH3+ with
  metal_bound_amines_as_nh2=False."""
  def n_h(h):
    return sorted(a.name.strip() for a in h.atoms()
      if a.parent().resname == 'HIS' and a.name.strip() in ('H', 'H1', 'H2',
      'H3'))
  _, log, h = run_vl('tst_validate_ligands_5_nh2', pdb_his_cu_glycan)
  lines = section(log)
  assert 'Metal-bound amines placed as NH2: 1' in lines, lines
  assert [l for l in lines if l.startswith('  N   HIS A   1: ')], lines
  assert len(n_h(h)) == 2, n_h(h)
  _, log, h = run_vl('tst_validate_ligands_5_nh3', pdb_his_cu_glycan,
    ['metal_bound_amines_as_nh2=False'])
  lines = section(log)
  assert ('Metal-bound amines as NH2: off (metal_bound_amines_as_nh2=False)'
    in lines), lines
  assert n_h(h) == ['H1', 'H2', 'H3'], n_h(h)
  print('test_metal_amine OK')

def test_all_off():
  """All three False: no flip, no copies, no NH2, each said so."""
  _, log, _ = run_vl('tst_validate_ligands_5_off', pdb_1avd_asn12, OFF)
  assert section(log) == [
    'Symmetry copies (crystal contacts): off (include_symmetric_copies=False)',
    'Flips: off (add_flip_movers=False), none made',
    'Metal-bound amines as NH2: off (metal_bound_amines_as_nh2=False)'], \
    section(log)
  print('test_all_off OK')

def test_flip_label():
  """The label of a flip key: insertion code, conformer, MODEL id."""
  label = val_lig.Program.flip_label
  assert label(('', 'A', '  12', ' ', 'ASN', '', 'amide')) == \
    'ASN A 12 (amide)'
  assert label(('', 'B', '  14', 'A', 'GLN', 'A', 'amide')) == \
    'GLN B 14A conformer A (amide)'
  assert label(('   2', 'C', '   5', ' ', 'HIS', 'B', 'His')) == \
    'model 2 HIS C 5 conformer B (His)'
  print('test_flip_label OK')

def test_flips_final_model():
  """His D next to two SO4 (a ring flip gives both N-H an acceptor). All
  blank: flipped and listed, marked from the ring atoms. SO4 split, B at the
  site and A 30 A away: the last pass (A) leaves the blank His unflipped, so
  it is not listed. His split: each conformer once."""
  mark = '  * within 5 A of ligand: SO4 D 101 3.0 A; SO4 D 102 2.9 A'
  _, log, _ = run_vl('tst_validate_ligands_5_his', trs.his_d_str('blank'))
  lines = section(log)
  assert lines[1:3] == ['Flipped Asn/Gln/His: 1', '  HIS D 1 (His)' + mark], \
    lines
  _, log, _ = run_vl('tst_validate_ligands_5_his_so4ab',
    trs.his_d_str('so4_split'))
  lines = section(log)
  assert 'Flipped Asn/Gln/His: 0' in lines, lines
  assert not [l for l in lines if 'HIS' in l], lines
  _, log, _ = run_vl('tst_validate_ligands_5_hisab', trs.his_d_str('his_split'))
  lines = section(log)
  assert lines[1:4] == ['Flipped Asn/Gln/His: 2',
    '  HIS D 1 conformer A (His)' + mark,
    '  HIS D 1 conformer B (His)' + mark], lines
  print('test_flips_final_model OK')

def test_ion_lock_flip():
  """A His whose flipped ring binds a Zn: flipped by the ion lock-down with
  the default, not with add_flip_movers=False (no heavy atom moves)."""
  _, log, _ = run_vl('tst_validate_ligands_5_zn_flip', trs.his_zn_str())
  lines = section(log)
  assert lines[1:3] == ['Flipped Asn/Gln/His: 1',
    '  HIS A 1 (His)  * within 5 A of ligand: EOH A 201 3.3 A'], lines
  _, log, _ = run_vl('tst_validate_ligands_5_zn', trs.his_zn_str(),
    ['add_flip_movers=False'])
  lines = section(log)
  assert lines[1:] == [
    'Flips: off (add_flip_movers=False), none made',
    'Metal-bound amines placed as NH2: 0'], lines
  print('test_ion_lock_flip OK')

def test_map_box_contacts():
  """His A 1 and a Cu one edge away in a P 1 box: across the model's own
  symmetry the amine binds the Cu (NH2); with a map of that box, whose box
  is not a lattice, it does not."""
  pdb_str = trs.his_cu_box_str()
  _, log, _ = run_vl('tst_validate_ligands_5_box', pdb_str)
  lines = section(log)
  assert 'Metal-bound amines placed as NH2: 1' in lines, lines
  assert [l for l in lines if l.startswith('  N   HIS A   1: ') and
    ' sym x-1,y,z ' in l], lines
  _, log, _ = run_vl('tst_validate_ligands_5_box_map', pdb_str,
    map_from=pdb_str)
  lines = section(log)
  assert ('Symmetry copies (crystal contacts): none: map input (the map box '
    'is not a lattice)') in lines, lines
  assert 'Metal-bound amines placed as NH2: 0' in lines, lines
  print('test_map_box_contacts OK')

def test_1avd_defaults():
  """1avd with its data and the defaults: the log section and values of
  the ligands near the flips."""
  import libtbx.load_env
  from libtbx.test_utils import approx_equal
  from mmtbx.regression.tst_validate_ligands import find_lr
  pdb = libtbx.env.find_in_repositories(
    relative_path='phenix_regression/pdb/pdb1avd.ent.gz', test=os.path.isfile)
  mtz = libtbx.env.find_in_repositories(
    relative_path='phenix_regression/reflection_files/1avd.mtz',
    test=os.path.isfile)
  if pdb is None or mtz is None:
    print('test_1avd_defaults skipped: no phenix_regression')
    return
  log = StringIO()
  result = run_program(program_class=val_lig.Program, args=[pdb, mtz],
    logger=log)
  lines = section(log.getvalue())
  assert lines[0] == 'Symmetry copies (crystal contacts): used', lines
  assert lines[1].startswith('  Symmetry copies (fixed environment): model 1,'
    ' 2731 copy atoms'), lines
  assert lines[2:] == ['Flipped Asn/Gln/His: 4',
    '  ASN A 12 (amide)  * within 5 A of ligand: BTN A 400 3.5 A',
    '  ASN A 24 (amide)', '  ASN A 69 (amide)', '  ASN B 54 (amide)',
    'Metal-bound amines placed as NH2: 0'], lines
  btn = find_lr(result.ligand_manager, 'chain A and resseq 400')
  # Asn A 12 flipped: other atoms within the radius than in the input (23.49)
  assert approx_equal(btn.get_adps().b_mean_within, 23.23, eps=0.01)
  assert approx_equal(btn.get_ccs().rscc, 0.969, eps=0.002)
  print('test_1avd_defaults OK')

# ------------------------------------------------------------------------------

if __name__ == '__main__':
  test_asn_flip()
  test_symmetry_contact()
  test_metal_amine()
  test_all_off()
  test_flip_label()
  test_flips_final_model()
  test_ion_lock_flip()
  test_map_box_contacts()
  test_1avd_defaults()
  print('OK')
