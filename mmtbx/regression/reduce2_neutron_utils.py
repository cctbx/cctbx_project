"""Inputs, runners, checks and shared runs for the tst_reduce2_neutron_*.py
tests (not a test file).

reduce2 with neutron distances: every step uses the same setting.
use_neutron_distances=True (and qmready=True, which switches it on) must give
H at the neutron restraint targets, water O-H from the HOH restraints, and an
Optimizer run with neutron radii (main run and flipkins). Bond lengths are
checked against each restraint's own target; bonds whose X-ray and neutron
targets differ by more than 0.04 A must also be far from the wrong target.
"""
from __future__ import absolute_import, division, print_function
import os
from six.moves import cStringIO as StringIO
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out, Sorry
from libtbx.test_utils import approx_equal
from scitbx import matrix
from cctbx import geometry_restraints
from cctbx.geometry_restraints.auto_linking_types import iterate_covalent
from cctbx.geometry_restraints.linking_class import linking_class
from iotbx.cli_parser import run_program
from mmtbx.programs import reduce2
from mmtbx.hydrogens import reduce_hydrogen
from mmtbx.reduce import Optimizers

# File-name prefix of what the helpers write (cached_run(),
# protein_ligand_src(), dipeptide_water_out()): run_tests runs the test files
# in parallel in one directory, so each file sets its own.
PREFIX = "tst_reduce2_neutron"

def set_prefix(prefix):
  global PREFIX
  PREFIX = prefix

# 1yjp: Asn/Gln for flipkins, Tyr OH, waters, P21
pdb_1yjp = """\
CRYST1   21.937    4.866   23.477  90.00 107.08  90.00 P 1 21 1
ATOM      1  N   GLY A   1      -9.009   4.612   6.102  1.00 16.77           N
ATOM      2  CA  GLY A   1      -9.052   4.207   4.651  1.00 16.57           C
ATOM      3  C   GLY A   1      -8.015   3.140   4.419  1.00 16.16           C
ATOM      4  O   GLY A   1      -7.523   2.521   5.381  1.00 16.78           O
ATOM      5  N   ASN A   2      -7.656   2.923   3.155  1.00 15.02           N
ATOM      6  CA  ASN A   2      -6.522   2.038   2.831  1.00 14.10           C
ATOM      7  C   ASN A   2      -5.241   2.537   3.427  1.00 13.13           C
ATOM      8  O   ASN A   2      -4.978   3.742   3.426  1.00 11.91           O
ATOM      9  CB  ASN A   2      -6.346   1.881   1.341  1.00 15.38           C
ATOM     10  CG  ASN A   2      -7.584   1.342   0.692  1.00 14.08           C
ATOM     11  OD1 ASN A   2      -8.025   0.227   1.016  1.00 17.46           O
ATOM     12  ND2 ASN A   2      -8.204   2.155  -0.169  1.00 11.72           N
ATOM     13  N   ASN A   3      -4.438   1.590   3.905  1.00 12.26           N
ATOM     14  CA  ASN A   3      -3.193   1.904   4.589  1.00 11.74           C
ATOM     15  C   ASN A   3      -1.955   1.332   3.895  1.00 11.10           C
ATOM     16  O   ASN A   3      -1.872   0.119   3.648  1.00 10.42           O
ATOM     17  CB  ASN A   3      -3.259   1.378   6.042  1.00 12.15           C
ATOM     18  CG  ASN A   3      -2.006   1.739   6.861  1.00 12.82           C
ATOM     19  OD1 ASN A   3      -1.702   2.925   7.072  1.00 15.05           O
ATOM     20  ND2 ASN A   3      -1.271   0.715   7.306  1.00 13.48           N
ATOM     21  N   GLN A   4      -1.005   2.228   3.598  1.00 10.29           N
ATOM     22  CA  GLN A   4       0.384   1.888   3.199  1.00 10.53           C
ATOM     23  C   GLN A   4       1.435   2.606   4.088  1.00 10.24           C
ATOM     24  O   GLN A   4       1.547   3.843   4.115  1.00  8.86           O
ATOM     25  CB  GLN A   4       0.656   2.148   1.711  1.00  9.80           C
ATOM     26  CG  GLN A   4       1.944   1.458   1.213  1.00 10.25           C
ATOM     27  CD  GLN A   4       2.504   2.044  -0.089  1.00 12.43           C
ATOM     28  OE1 GLN A   4       2.744   3.268  -0.190  1.00 14.62           O
ATOM     29  NE2 GLN A   4       2.750   1.161  -1.091  1.00  9.05           N
ATOM     30  N   GLN A   5       2.154   1.821   4.871  1.00 10.38           N
ATOM     31  CA  GLN A   5       3.270   2.361   5.640  1.00 11.39           C
ATOM     32  C   GLN A   5       4.594   1.768   5.172  1.00 11.52           C
ATOM     33  O   GLN A   5       4.768   0.546   5.054  1.00 12.05           O
ATOM     34  CB  GLN A   5       3.056   2.183   7.147  1.00 11.96           C
ATOM     35  CG  GLN A   5       1.829   2.950   7.647  1.00 10.81           C
ATOM     36  CD  GLN A   5       1.344   2.414   8.954  1.00 13.10           C
ATOM     37  OE1 GLN A   5       0.774   1.325   9.002  1.00 10.65           O
ATOM     38  NE2 GLN A   5       1.549   3.187  10.039  1.00 12.30           N
ATOM     39  N   ASN A   6       5.514   2.664   4.856  1.00 11.99           N
ATOM     40  CA  ASN A   6       6.831   2.310   4.318  1.00 12.30           C
ATOM     41  C   ASN A   6       7.854   2.761   5.324  1.00 13.40           C
ATOM     42  O   ASN A   6       8.219   3.943   5.374  1.00 13.92           O
ATOM     43  CB  ASN A   6       7.065   3.016   2.993  1.00 12.13           C
ATOM     44  CG  ASN A   6       5.961   2.735   2.003  1.00 12.77           C
ATOM     45  OD1 ASN A   6       5.798   1.604   1.551  1.00 14.27           O
ATOM     46  ND2 ASN A   6       5.195   3.747   1.679  1.00 10.07           N
ATOM     47  N   TYR A   7       8.292   1.817   6.147  1.00 14.70           N
ATOM     48  CA  TYR A   7       9.159   2.144   7.299  1.00 15.18           C
ATOM     49  C   TYR A   7      10.603   2.331   6.885  1.00 15.91           C
ATOM     50  O   TYR A   7      11.041   1.811   5.855  1.00 15.76           O
ATOM     51  CB  TYR A   7       9.061   1.065   8.369  1.00 15.35           C
ATOM     52  CG  TYR A   7       7.665   0.929   8.902  1.00 14.45           C
ATOM     53  CD1 TYR A   7       6.771   0.021   8.327  1.00 15.68           C
ATOM     54  CD2 TYR A   7       7.210   1.756   9.920  1.00 14.80           C
ATOM     55  CE1 TYR A   7       5.480  -0.094   8.796  1.00 13.46           C
ATOM     56  CE2 TYR A   7       5.904   1.649  10.416  1.00 14.33           C
ATOM     57  CZ  TYR A   7       5.047   0.729   9.831  1.00 15.09           C
ATOM     58  OH  TYR A   7       3.766   0.589  10.291  1.00 14.39           O
ATOM     59  OXT TYR A   7      11.358   2.999   7.612  1.00 17.49           O
HETATM   61  O   HOH A   8      -6.471   5.227   7.124  1.00 22.62           O
HETATM   62  O   HOH A   9      10.431   1.858   3.216  1.00 19.71           O
HETATM   63  O   HOH A  10     -11.286   1.756  -1.468  1.00 17.08           O
HETATM   64  O   HOH A  11      11.808   4.179   9.970  1.00 23.99           O
HETATM   65  O   HOH A  12      13.605   1.327   9.198  1.00 26.17           O
HETATM   66  O   HOH A  13      -2.749   3.429  10.024  1.00 39.15           O
HETATM   67  O   HOH A  14      -1.500   0.682  10.967  1.00 43.49           O
END
"""

# Lys 3 stops at CD: its HD2/HD3 cannot be parameterized, so placement drops
# them with model.select()
pdb_truncated_lys = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
ATOM      1  N   ALA A   1      11.204   2.315  10.312  1.00 20.00           N
ATOM      2  CA  ALA A   1      12.555   2.835  10.528  1.00 20.00           C
ATOM      3  C   ALA A   1      13.588   1.726  10.381  1.00 20.00           C
ATOM      4  O   ALA A   1      13.325   0.567  10.721  1.00 20.00           O
ATOM      5  CB  ALA A   1      12.866   3.965   9.556  1.00 20.00           C
ATOM      6  N   SER A   2      14.767   2.090   9.884  1.00 20.00           N
ATOM      7  CA  SER A   2      15.853   1.130   9.701  1.00 20.00           C
ATOM      8  C   SER A   2      16.968   1.370  10.713  1.00 20.00           C
ATOM      9  O   SER A   2      17.281   2.516  11.041  1.00 20.00           O
ATOM     10  CB  SER A   2      16.411   1.224   8.278  1.00 20.00           C
ATOM     11  OG  SER A   2      15.456   0.826   7.311  1.00 20.00           O
ATOM     12  N   LYS A   3      17.560   0.282  11.199  1.00 20.00           N
ATOM     13  CA  LYS A   3      18.643   0.370  12.172  1.00 20.00           C
ATOM     14  C   LYS A   3      19.945  -0.148  11.570  1.00 20.00           C
ATOM     15  O   LYS A   3      19.964  -1.176  10.892  1.00 20.00           O
ATOM     16  CB  LYS A   3      18.290  -0.424  13.432  1.00 20.00           C
ATOM     17  CG  LYS A   3      19.308  -0.301  14.556  1.00 20.00           C
ATOM     18  CD  LYS A   3      18.913  -1.137  15.762  1.00 20.00           C
ATOM     19  OXT LYS A   3      20.995   0.451  11.763  1.00 20.00           O
HETATM   20  O   HOH A   4      24.000   6.000  20.000  1.00 20.00           O
END
"""

INFO_TRUE = "useNeutronDistances = True"
INFO_FALSE = "useNeutronDistances = False"

def run_reduce2(prefix, pdb_str, args=()):
  fn_in, fn_out = prefix + ".pdb", prefix + "_out.pdb"
  fn_txt = prefix + ".txt"
  with open(fn_in, "w") as f:
    f.write(pdb_str)
  for fn in (fn_out, fn_txt):
    if os.path.exists(fn): os.remove(fn)
  log = StringIO()
  run_program(program_class=reduce2.Program, logger=log,
    args=[fn_in, "output.filename=" + fn_out, "overwrite=True",
          "output.description_file_name=" + fn_txt] + list(args))
  with open(fn_txt) as f:
    description = f.read()
  return fn_out, log.getvalue(), description

def interpret(file_name, use_neutron_distances):
  inp = iotbx.pdb.input(file_name)
  model = mmtbx.model.manager(model_input=inp, log=null_out())
  p = reduce_hydrogen.get_reduce_pdb_interpretation_params(
    use_neutron_distances)
  model.process(make_restraints=True, pdb_interpretation_params=p)
  return model

def h_bond_lengths(model):
  """(i_seq, j_seq, j_sym) -> (ideal, actual) for covalent bonds with an H/D
  end; asu bonds with their symmetry operator"""
  metal = linking_class().get_origin_id('metal coordination')
  covalent = set(iterate_covalent()) - set([metal])
  sites_cart = model.get_sites_cart()
  hd = model.get_hd_selection()
  bond_proxies = model.get_restraints_manager().geometry.pair_proxies(
    sites_cart=sites_cart).bond_proxies
  result = {}
  for p in bond_proxies.simple:
    i, j = p.i_seqs
    if p.origin_id not in covalent or not (hd[i] or hd[j]): continue
    b = geometry_restraints.bond(sites_cart=sites_cart, proxy=p)
    result[(i, j, 0)] = (p.distance_ideal, b.distance_model)
  asu_mappings = bond_proxies.asu_mappings()
  for p in bond_proxies.asu:
    if p.origin_id not in covalent or not (hd[p.i_seq] or hd[p.j_seq]):
      continue
    b = geometry_restraints.bond(sites_cart=sites_cart,
      asu_mappings=asu_mappings, proxy=p)
    result[(p.i_seq, p.j_seq, p.j_sym)] = (p.distance_ideal, b.distance_model)
  return result

def check_bond_lengths(file_name, use_neutron_distances):
  right = h_bond_lengths(interpret(file_name, use_neutron_distances))
  wrong = h_bond_lengths(interpret(file_name, not use_neutron_distances))
  assert len(right) > 0
  n_discriminating = 0
  for key, (ideal, actual) in right.items():
    assert abs(actual - ideal) < 0.02, (file_name, key, ideal, actual)
    if key not in wrong: continue
    ideal_wrong = wrong[key][0]
    if abs(ideal - ideal_wrong) > 0.04:
      assert abs(actual - ideal_wrong) > 0.02, (file_name, key, ideal_wrong)
      n_discriminating += 1
  assert n_discriminating > 0, file_name
  return n_discriminating

def water_geometry(file_name):
  """per water O and conformer: ([O-H/D lengths], H-O-H angle or None);
  conformer residues combine blank-altloc O with altloc H/D"""
  h = iotbx.pdb.input(file_name).construct_hierarchy()
  seen = set()
  result = []
  for chain in h.only_model().chains():
    for conformer in chain.conformers():
      for residue in conformer.residues():
        if residue.resname not in ("HOH", "DOD"): continue
        ra = residue.atoms()
        o = [a for a in ra if a.element.strip() == "O"]
        hs = [a for a in ra if a.element.strip() in ("H", "D")]
        if len(o) != 1 or not hs: continue
        key = tuple(sorted(a.i_seq for a in [o[0]] + hs))
        if key in seen: continue
        seen.add(key)
        vs = [matrix.col(a.xyz) - matrix.col(o[0].xyz) for a in hs]
        angle = None
        if len(vs) == 2: angle = vs[0].angle(vs[1], deg=True)
        result.append(([v.length() for v in vs], angle))
  return result

def check_water(file_name, d_ideal, n_expected):
  waters = water_geometry(file_name)
  assert len(waters) == n_expected, (file_name, len(waters))
  for lengths, angle in waters:
    for d in lengths:
      assert approx_equal(d, d_ideal, eps=0.005, out=None), (file_name, d)
    if angle is not None:
      assert approx_equal(angle, 103.91, eps=0.1, out=None), (file_name, angle)

def check_no_d(fn):
  h = iotbx.pdb.input(fn).construct_hierarchy()
  assert "D" not in [a.element.strip() for a in h.atoms()]
  assert "DOD" not in [ag.resname for ag in h.atom_groups()]

FLOATING = "H/D without a covalent partner, removed before optimization:"

def h_input_without_tyr_oh(**params):
  """1yjp with H from reduce2, then Tyr OH deleted: its HH has no partner"""
  lines = [l for l in cached_run(pdb_1yjp, **params).text.splitlines()
           if " OH  TYR" not in l]
  return "\n".join(lines) + "\n"

def atom_labels(file_name):
  h = iotbx.pdb.input(file_name).construct_hierarchy()
  return sorted((a.parent().parent().resseq, a.parent().resname, a.name,
                 a.parent().altloc) for a in h.atoms())

# ------------------------------------------------------------------------------
# step 3.1: phil, validate(), file names, _RebuildModel()

OPTIONS = ("none", "all_h", "all_d", "hd_and_h", "hd_and_d", "all_hd")

def make_program(prefix, pdb_str, cif_str=None, **params):
  """reduce2.Program as run_program builds it, validate() not yet called;
  params: name=value, with __ for the scope separator"""
  from iotbx.cli_parser import CCTBXParser
  with open(prefix + ".pdb", "w") as f:
    f.write(pdb_str)
  args = [prefix + ".pdb"]
  cif_strs = [] if cif_str is None else (
    [cif_str] if isinstance(cif_str, str) else list(cif_str))
  for i, text in enumerate(cif_strs):
    fn = prefix + ("_%d" % i if i else "") + ".cif"
    with open(fn, "w") as f:
      f.write(text)
    args.append(fn)
  args += ["%s=%s" % (name.replace("__", "."), value)
           for name, value in params.items()]
  parser = CCTBXParser(program_class=reduce2.Program, logger=null_out())
  parser.parse_args(args)
  log = StringIO()
  prog = reduce2.Program(parser.data_manager, parser.working_phil.extract(),
    logger=log)
  return prog, log

def check_sorry(prefix, pdb_str, text, **params):
  prog, log = make_program(prefix, pdb_str, **params)
  try:
    prog.validate()
  except Sorry as e:
    assert text in str(e), str(e)
  else:
    raise AssertionError("no Sorry for %s" % params)

pdb_ser_dg = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
ATOM      1  N   SER A   2      14.767   2.090   9.884  1.00 20.00           N
ATOM      2  CA  SER A   2      15.853   1.130   9.701  1.00 20.00           C
ATOM      3  C   SER A   2      16.968   1.370  10.713  1.00 20.00           C
ATOM      4  O   SER A   2      17.281   2.516  11.041  1.00 20.00           O
ATOM      5  CB  SER A   2      16.411   1.224   8.278  1.00 20.00           C
ATOM      6  OG  SER A   2      15.456   0.826   7.311  1.00 20.00           O
ATOM      7  DG  SER A   2      15.856   0.874   6.441  1.00 20.00           D
END
"""

# ligand unknown to the CCD and the monomer library; C1-O1 1.4321 A
pdb_lig = """\
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1
HETATM    1  C1   01 A   1      10.000  10.000  10.000  1.00 20.00           C
HETATM    2  O1   01 A   1      11.432  10.000  10.000  1.00 20.00           O
HETATM    3  H11  01 A   1       9.640  11.020  10.000  1.00 20.00           H
HETATM    4  H12  01 A   1       9.640   9.490  10.880  1.00 20.00           H
HETATM    5  H13  01 A   1       9.640   9.490   9.120  1.00 20.00           H
HETATM    6  HO1  01 A   1      11.750  10.900  10.000  1.00 20.00           H
END
"""
cif_lig = """\
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
 01  01  'test methanol'  ligand  6  2  .

data_comp_01
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
 01  C1   C  CH3  0   0.000   0.000   0.000
 01  O1   O  OH1  0   1.432   0.000   0.000
 01  H11  H  HCH3 0  -0.360   1.020   0.000
 01  H12  H  HCH3 0  -0.360  -0.510   0.880
 01  H13  H  HCH3 0  -0.360  -0.510  -0.880
 01  HO1  H  HOH1 0   1.750   0.900   0.000

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
_chem_comp_bond.value_dist_neutron
 01  C1  O1   single  1.4321  0.020  1.4321
 01  C1  H11  single  0.970   0.020  1.090
 01  C1  H12  single  0.970   0.020  1.090
 01  C1  H13  single  0.970   0.020  1.090
 01  O1  HO1  single  0.850   0.020  0.980

loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
 01  H11  C1  H12  109.5  3.0
 01  H11  C1  H13  109.5  3.0
 01  H12  C1  H13  109.5  3.0
 01  O1   C1  H11  109.5  3.0
 01  O1   C1  H12  109.5  3.0
 01  O1   C1  H13  109.5  3.0
 01  C1   O1  HO1  109.5  3.0
"""

# ------------------------------------------------------------------------------
# step 3.2: output step for add and place

NEUTRON = ("all_h", "all_d", "hd_and_h", "hd_and_d", "all_hd")

def run_prog(prefix, pdb_str, cif_str=None, **params):
  """validate() and run() a reduce2.Program; returns the program"""
  params.setdefault("output__filename", prefix + "_out.pdb")
  params.setdefault("output__description_file_name", prefix + ".txt")
  params.setdefault("output__overwrite", True)
  prog, log = make_program(prefix, pdb_str, cif_str, **params)
  prog.validate()
  prog.run()
  prog.log_text = log.getvalue()
  return prog

_cached_runs = {}
def cached_run(pdb_str, cif_str=None, **params):
  """run_prog() once per input and parameters, shared by tests that only
  read it: the program (with .text, .description) and its files
  <PREFIX>_cache_<n>*"""
  key = (pdb_str, cif_str if cif_str is None or isinstance(cif_str, str)
    else tuple(cif_str), tuple(sorted(params.items())))
  if key not in _cached_runs:
    prog = run_prog(PREFIX + "_cache_%d" % len(_cached_runs),
      pdb_str, cif_str, **params)
    with open(prog.params.output.filename) as f: prog.text = f.read()
    with open(prog.params.output.description_file_name) as f:
      prog.description = f.read()
    _cached_runs[key] = prog
  return _cached_runs[key]

def hd_pairs(hierarchy):
  """[(H, D)] from exchangeable_hd_selections()"""
  atoms = hierarchy.atoms()
  result = []
  for x in hierarchy.exchangeable_hd_selections():
    a, b = atoms[x[0][0]], atoms[x[1][0]]
    if a.element.strip() == "D": a, b = b, a
    result.append((a, b))
  return result

def heavy_atoms(hierarchy):
  from mmtbx.hydrogens import hd_exchange
  return dict((hd_exchange.atom_key(a), (a.name, a.parent().altloc, a.xyz,
    a.occ, a.b, tuple(a.uij) if a.uij_is_defined() else None))
    for a in hierarchy.atoms() if a.element.strip() not in ("H", "D"))

def pdb_1yjp_h():
  return cached_run(pdb_1yjp).text

def make_model_from_str(pdb_str):
  model = mmtbx.model.manager(model_input=iotbx.pdb.input(
    lines=pdb_str.split("\n"), source_info=None), log=null_out())
  model.process(make_restraints=True)
  return model

pdb_1yjp_dod = pdb_1yjp.replace("HOH A", "DOD A")

RENAMED = "water(s) renamed DOD -> HOH (their D were replaced by H)"

def water_groups(fn):
  h = iotbx.pdb.input(fn).construct_hierarchy()
  return [(ag.resname, sorted(set(a.element.strip() for a in ag.atoms())))
          for ag in h.atom_groups() if ag.resname in ("HOH", "DOD")]

pdb_1yjp_dod_d = pdb_1yjp_dod.replace(
  "HETATM   61  O   DOD A   8      -6.471   5.227   7.124  1.00 22.62"
  "           O\n",
  "HETATM   61  O   DOD A   8      -6.471   5.227   7.124  1.00 22.62"
  "           O\n"
  "HETATM   61  D1  DOD A   8      -5.621   5.227   7.124  1.00 22.62"
  "           D\n"
  "HETATM   61  D2  DOD A   8      -6.675   6.052   7.124  1.00 22.62"
  "           D\n")

# ------------------------------------------------------------------------------
# step 3.3: per-mode labels, His, ligands, glycan, mmCIF, description

def atom_label_set(src, is_file=True):
  """{(resseq, name, altloc, element, occupancy)}"""
  inp = iotbx.pdb.input(file_name=src) if is_file else iotbx.pdb.input(
    source_info=None, lines=src)
  h = inp.construct_hierarchy()
  return set((a.parent().parent().resseq.strip(), a.name.strip(),
              a.parent().altloc, a.element.strip(), round(a.occ, 2))
             for a in h.atoms())

def hd_table(text):
  """{(resseq, name, altloc, element, occupancy)} from lines
  'resseq name altloc element occupancy', altloc '-' for blank"""
  result = set()
  for line in text.strip().splitlines():
    resseq, name, altloc, element, occ = line.split()
    result.add((resseq, name, "" if altloc == "-" else altloc, element,
                float(occ)))
  return result

def hd_labels(fn):
  return set(x for x in atom_label_set(fn) if x[3] in ("H", "D"))

# Ser-Leu dipeptide; Leu stops after CA
dipeptide = """\
CRYST1   14.228   15.246   13.577  90.00  90.00  90.00 P 1
ATOM    577  N   SER A  37       8.872   5.000   5.995  1.00 25.33           N
ATOM    578  CA  SER A  37       7.827   5.915   6.440  1.00 31.94           C
ATOM    579  CB  SER A  37       7.709   5.892   7.965  1.00 42.16           C
ATOM    580  OG  SER A  37       8.926   6.283   8.577  1.00 54.74           O
ATOM    581  C   SER A  37       8.095   7.336   5.956  1.00 32.18           C
ATOM    582  O   SER A  37       9.228   7.814   6.002  1.00 31.99           O
ATOM    588  N   LEU A  38       7.045   8.005   5.492  1.00 34.04           N
ATOM    589  CA  LEU A  38       7.165   9.372   5.000  1.00 36.58           C
ATOM    594  C   LEU A  38       6.026  10.246   5.515  1.00 37.03           C
ATOM    595  O   LEU A  38       5.000   9.741   5.969  1.00 35.00           O
END
"""
# expected H/D per mode: the chain-start N-H (Ser 37 H, not resseq 1, so no
# NH3) and HA on the truncated Leu are placed like any other H
dipeptide_hd = {
"all_h": """
37 H - H 1.0
37 HA - H 1.0
37 HB2 - H 1.0
37 HB3 - H 1.0
37 HG - H 1.0
38 H - H 1.0
38 HA - H 1.0""",
"all_d": """
37 D - D 1.0
37 DA - D 1.0
37 DB2 - D 1.0
37 DB3 - D 1.0
37 DG - D 1.0
38 D - D 1.0
38 DA - D 1.0""",
"hd_and_h": """
37 H A H 0.5
37 D B D 0.5
37 HG A H 0.5
37 DG B D 0.5
37 HA - H 1.0
37 HB2 - H 1.0
37 HB3 - H 1.0
38 H A H 0.5
38 D B D 0.5
38 HA - H 1.0""",
"hd_and_d": """
37 H A H 0.5
37 D B D 0.5
37 HG A H 0.5
37 DG B D 0.5
37 DA - D 1.0
37 DB2 - D 1.0
37 DB3 - D 1.0
38 H A H 0.5
38 D B D 0.5
38 DA - D 1.0""",
"all_hd": """
37 H A H 0.5
37 D B D 0.5
37 HA A H 0.5
37 DA B D 0.5
37 HB2 A H 0.5
37 DB2 B D 0.5
37 HB3 A H 0.5
37 DB3 B D 0.5
37 HG A H 0.5
37 DG B D 0.5
38 H A H 0.5
38 D B D 0.5
38 HA A H 0.5
38 DA B D 0.5""",
}

his_no_h = """\
CRYST1   12.670   14.862   15.186  90.00  90.00  90.00 P 1
ATOM    780  N   HIS A  45       5.209   7.219   5.800  1.00  9.38           N
ATOM    781  CA  HIS A  45       6.450   6.586   6.229  1.00 10.30           C
ATOM    782  C   HIS A  45       6.629   5.221   5.572  1.00 11.06           C
ATOM    783  O   HIS A  45       7.669   4.940   4.976  1.00 12.76           O
ATOM    784  CB  HIS A  45       6.482   6.444   7.752  1.00 11.71           C
ATOM    785  CG  HIS A  45       6.466   7.752   8.482  1.00 17.20           C
ATOM    786  ND1 HIS A  45       6.374   7.837   9.854  1.00 28.53           N
ATOM    787  CD2 HIS A  45       6.532   9.027   8.030  1.00 22.78           C
ATOM    788  CE1 HIS A  45       6.383   9.107  10.217  1.00 26.38           C
ATOM    789  NE2 HIS A  45       6.478   9.850   9.129  1.00 20.23           N
END
"""
his_hd = """
45 H A H 0.5
45 D B D 0.5
45 HD1 A H 0.5
45 DD1 B D 0.5
45 HE2 A H 0.5
45 DE2 B D 0.5
45 HA - H 1.0
45 HB2 - H 1.0
45 HB3 - H 1.0
45 HD2 - H 1.0
45 HE1 - H 1.0"""

hc4_no_h = """\
CRYST1   67.145   67.145   41.119  90.00  90.00 120.00 P 63          6
HETATM 1987  C1  HC4 A 169      12.270  -2.851   7.690  1.00  2.51           C
HETATM 1988  O1  HC4 A 169      12.412  -3.051   8.915  1.00  5.58           O
HETATM 1989  C2  HC4 A 169      13.073  -1.969   6.902  1.00  3.66           C
HETATM 1990  C3  HC4 A 169      13.948  -1.115   7.480  1.00  4.10           C
HETATM 1991  C1' HC4 A 169      14.960  -0.357   6.862  1.00  3.50           C
HETATM 1992  C2' HC4 A 169      15.734   0.515   7.642  1.00  7.24           C
HETATM 1993  C3' HC4 A 169      16.796   1.251   7.105  1.00  7.06           C
HETATM 1994  C4' HC4 A 169      17.093   1.090   5.743  1.00  6.80           C
HETATM 1995  C5' HC4 A 169      16.317   0.232   4.947  1.00  7.54           C
HETATM 1996  C6' HC4 A 169      15.255  -0.474   5.508  1.00  4.60           C
HETATM 1997  O4' HC4 A 169      18.097   1.748   5.166  1.00  7.92           O
END
"""
hc4_all_d = """
169 D2 - D 1.0
169 D3 - D 1.0
169 D2' - D 1.0
169 D3' - D 1.0
169 D5' - D 1.0
169 D6' - D 1.0
169 DO4' - D 1.0"""

def make_cif(code, atoms, bonds, angles):
  """monomer CIF: atoms (name, element, energy type, x, y, z), bonds
  (a, b, x-ray, neutron), angles (a, b, c, value)"""
  lines = ["data_comp_list", "loop_", "_chem_comp.id",
    "_chem_comp.three_letter_code", "_chem_comp.name", "_chem_comp.group",
    "_chem_comp.number_atoms_all", "_chem_comp.number_atoms_nh",
    "_chem_comp.desc_level",
    " %s  %s  'test ligand'  ligand  %d  %d  ." % (code, code, len(atoms),
      len([x for x in atoms if x[1] != "H"])), "", "data_comp_" + code,
    "loop_", "_chem_comp_atom.comp_id", "_chem_comp_atom.atom_id",
    "_chem_comp_atom.type_symbol", "_chem_comp_atom.type_energy",
    "_chem_comp_atom.charge", "_chem_comp_atom.x", "_chem_comp_atom.y",
    "_chem_comp_atom.z"]
  lines += [" %s  %-4s  %s  %-4s  0  %.3f  %.3f  %.3f" % ((code,) + x)
            for x in atoms]
  lines += ["", "loop_", "_chem_comp_bond.comp_id",
    "_chem_comp_bond.atom_id_1", "_chem_comp_bond.atom_id_2",
    "_chem_comp_bond.type", "_chem_comp_bond.value_dist",
    "_chem_comp_bond.value_dist_esd", "_chem_comp_bond.value_dist_neutron"]
  lines += [" %s  %-4s  %-4s  single  %.3f  0.020  %.3f" % ((code,) + x)
            for x in bonds]
  lines += ["", "loop_", "_chem_comp_angle.comp_id",
    "_chem_comp_angle.atom_id_1", "_chem_comp_angle.atom_id_2",
    "_chem_comp_angle.atom_id_3", "_chem_comp_angle.value_angle",
    "_chem_comp_angle.value_angle_esd"]
  lines += [" %s  %-4s  %-4s  %-4s  %.1f  3.0" % ((code,) + x)
            for x in angles]
  return "\n".join(lines) + "\n"

def heavy_pdb(code, atoms, offset=(10.0, 10.0, 10.0), cryst1=True):
  lines = []
  if cryst1:
    lines.append("CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1")
  for i, (name, element, energy, x, y, z) in enumerate(atoms):
    if element == "H": continue
    lines.append("HETATM%5d %-4s %3s A   1    %8.3f%8.3f%8.3f  1.00 20.00"
      "          %2s" % (i + 1, name if len(name) == 4 else " " + name, code,
      x + offset[0], y + offset[1], z + offset[2], element))
  return "\n".join(lines) + "\n"

def _tetra(center, a, b):
  """angles for the H on a CH2 between a and b"""
  return [(a, center, "H%s1" % center[-1]), (a, center, "H%s2" % center[-1]),
          (b, center, "H%s1" % center[-1]), (b, center, "H%s2" % center[-1]),
          ("H%s1" % center[-1], center, "H%s2" % center[-1])]

# 03: ethylene glycol; one hydroxyl H named 1HO2 (digit first)
zq3_atoms = [
  ("O2", "O", "OH1", -0.700, 1.210, 0.0), ("C1", "C", "CH2", 0.0, 0.0, 0.0),
  ("C2", "C", "CH2", 1.520, 0.0, 0.0), ("O3", "O", "OH1", 2.220, 1.210, 0.0),
  ("1HO2", "H", "HOH1", -1.600, 1.000, 0.0),
  ("H11", "H", "HCH2", -0.330, -0.550, 0.890),
  ("H12", "H", "HCH2", -0.330, -0.550, -0.890),
  ("H21", "H", "HCH2", 1.850, -0.550, 0.890),
  ("H22", "H", "HCH2", 1.850, -0.550, -0.890),
  ("HO3", "H", "HOH1", 3.120, 1.000, 0.0)]
zq3_cif = make_cif("03", zq3_atoms,
  [("C1", "O2", 1.43, 1.43), ("C1", "C2", 1.52, 1.52),
   ("C2", "O3", 1.43, 1.43), ("O2", "1HO2", 0.85, 0.98),
   ("O3", "HO3", 0.85, 0.98), ("C1", "H11", 0.97, 1.09),
   ("C1", "H12", 0.97, 1.09), ("C2", "H21", 0.97, 1.09),
   ("C2", "H22", 0.97, 1.09)],
  [(x[0], x[1], x[2], 109.5) for x in
    [("O2", "C1", "C2"), ("C1", "C2", "O3"), ("C1", "O2", "1HO2"),
     ("C2", "O3", "HO3")] + _tetra("C1", "O2", "C2") +
    _tetra("C2", "C1", "O3")])

# 04: N-methyl ethanolamine, NH and OH (a terminal NH2 has no dihedral
# anchor for its riding H and is not placed: one heavy neighbour)
zq4_atoms = [
  ("C3", "C", "CH3", -1.430, 1.800, 0.0), ("N1", "N", "NH1", -0.700, 1.210, 0.0),
  ("C1", "C", "CH2", 0.0, 0.0, 0.0), ("C2", "C", "CH2", 1.520, 0.0, 0.0),
  ("O1", "O", "OH1", 2.220, 1.210, 0.0),
  ("HN1", "H", "HNH1", -0.200, 1.900, 0.500),
  ("H31", "H", "HCH3", -2.400, 1.400, 0.0),
  ("H32", "H", "HCH3", -1.300, 2.400, 0.880),
  ("H33", "H", "HCH3", -1.300, 2.400, -0.880),
  ("H11", "H", "HCH2", -0.330, -0.550, 0.890),
  ("H12", "H", "HCH2", -0.330, -0.550, -0.890),
  ("H21", "H", "HCH2", 1.850, -0.550, 0.890),
  ("H22", "H", "HCH2", 1.850, -0.550, -0.890),
  ("HO1", "H", "HOH1", 3.120, 1.000, 0.0)]
zq4_cif = make_cif("04", zq4_atoms,
  [("C3", "N1", 1.47, 1.47), ("C1", "N1", 1.47, 1.47),
   ("C1", "C2", 1.52, 1.52), ("C2", "O1", 1.43, 1.43),
   ("N1", "HN1", 0.89, 1.01), ("O1", "HO1", 0.85, 0.98),
   ("C3", "H31", 0.96, 1.09), ("C3", "H32", 0.96, 1.09),
   ("C3", "H33", 0.96, 1.09), ("C1", "H11", 0.97, 1.09),
   ("C1", "H12", 0.97, 1.09), ("C2", "H21", 0.97, 1.09),
   ("C2", "H22", 0.97, 1.09)],
  [(x[0], x[1], x[2], 109.5) for x in
    [("C3", "N1", "C1"), ("N1", "C1", "C2"), ("C1", "C2", "O1"),
     ("C3", "N1", "HN1"), ("C1", "N1", "HN1"), ("C2", "O1", "HO1"),
     ("N1", "C3", "H31"), ("N1", "C3", "H32"), ("N1", "C3", "H33"),
     ("H31", "C3", "H32"), ("H31", "C3", "H33"), ("H32", "C3", "H33")] +
    _tetra("C1", "N1", "C2") + _tetra("C2", "C1", "O1")])

MISSING_RESTRAINTS = "Restraints were not found for the following residues:"
UNCHANGED = "Residues without restraints, left unchanged:"

# ------------------------------------------------------------------------------
# step 4.1: input check and normalization

def hd_values(fn):
  """{(resseq, name, altloc, element): (occupancy, B)} of the H/D"""
  h = iotbx.pdb.input(fn).construct_hierarchy()
  return dict(((a.parent().parent().resseq.strip(), a.name.strip(),
                a.parent().altloc, a.element.strip()),
               (round(a.occ, 2), round(a.b, 2)))
              for a in h.atoms() if a.element.strip() in ("H", "D"))

def edit_atoms(pdb_str, edit):
  """pdb_str with edit(atom) applied to every atom (CRYST1 kept)"""
  inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  h = inp.construct_hierarchy()
  for a in h.atoms():
    edit(a)
  return h.as_pdb_string(crystal_symmetry=inp.crystal_symmetry())

def dipeptide_out(option):
  """written output of the Ser-Leu dipeptide in mode option"""
  return cached_run(dipeptide, neutron_option=option).text

# Ser-Leu dipeptide plus a water, for the reruns and the cross-mode runs
dipeptide_water = dipeptide.replace("END\n",
  "HETATM  600  O   HOH A 101      11.500  11.000   9.500  1.00 30.00"
  "           O\nEND\n")
# exchangeable sites (H name) of the dipeptide; the water is separate
dipeptide_exchangeable = set([("37", "H"), ("37", "HG"), ("38", "H")])
# the Optimizer rotates these (Ser 37 N has one heavy neighbour)
dipeptide_rotatable = set([("37", "H"), ("37", "HG")])

_dipeptide_water_out = {}
def dipeptide_water_out(option):
  """written output of dipeptide_water in mode option, exclude_water=False"""
  if option not in _dipeptide_water_out:
    fn, _, _ = run_reduce2(PREFIX + "_dwo_" + option,
      dipeptide_water, ["neutron_option=" + option, "exclude_water=False"])
    with open(fn) as f: _dipeptide_water_out[option] = f.read()
  return _dipeptide_water_out[option]

def h_name(name, element):
  """the H name of an H or D (first D -> H for D)"""
  return name.replace("D", "H", 1) if element == "D" else name

def all_atoms(src, is_file=True):
  """{(resseq, name, altloc, element): (occupancy, B, xyz)}"""
  inp = iotbx.pdb.input(file_name=src) if is_file else iotbx.pdb.input(
    source_info=None, lines=src)
  h = inp.construct_hierarchy()
  return dict(((a.parent().parent().resseq.strip(), a.name.strip(),
                a.parent().altloc, a.element.strip()),
               (round(a.occ, 2), round(a.b, 2), matrix.col(a.xyz)))
              for a in h.atoms())

def pairs_03_07():
  """hd_and_h dipeptide, each pair H 0.3 B 15 / D 0.7 B 25 (parents 1.0)"""
  def edit(a):
    if a.parent().altloc == "A":
      a.set_occ(0.3); a.set_b(15.0)
    elif a.parent().altloc == "B":
      a.set_occ(0.7); a.set_b(25.0)
  return edit_atoms(dipeptide_out("hd_and_h"), edit)

def dg_moved():
  """hd_and_h dipeptide with Ser DG 0.5 A away from HG"""
  def edit(a):
    if a.name.strip() == "DG":
      a.set_xyz((a.xyz[0] + 0.5, a.xyz[1], a.xyz[2]))
  src = edit_atoms(dipeptide_out("hd_and_h"), edit)
  assert src != dipeptide_out("hd_and_h")
  return src

def ser_cif_hg_cb():
  """geostd SER with an extra CB-HG bond (below the covalent limit)"""
  import libtbx.load_env
  fn = libtbx.env.find_in_repositories(
    relative_path="chem_data/geostd/s/data_SER.cif", test=os.path.isfile)
  with open(fn) as f: cif = f.read()
  old = " SER  OG  HG   coval  0.840  0.020  0.980\n"
  assert cif.count(old) == 1
  return cif.replace(old, old + " SER  CB  HG   coval  1.100  0.020  1.100\n")

def atom_records(fn):
  """[(resseq, resname, name, altloc, element, xyz, occ, b)] in file order"""
  h = iotbx.pdb.input(fn).construct_hierarchy()
  return [(a.parent().parent().resseq.strip(), a.parent().resname,
           a.name.strip(), a.parent().altloc, a.element.strip(),
           tuple(round(x, 3) for x in a.xyz), round(a.occ, 2), round(a.b, 2))
          for a in h.atoms()]

def check_removed(fn_in, fn_out):
  """no H/D, no empty or H/D-only atom_groups, no DOD, heavy atoms as in
  the input (DOD read as HOH)"""
  h = iotbx.pdb.input(fn_out).construct_hierarchy()
  assert not [a for a in h.atoms() if a.element.strip() in ("H", "D")]
  for ag in h.atom_groups():
    assert ag.atoms_size() > 0 and ag.resname != "DOD", ag.id_str()
  heavy = lambda records: sorted(
    (r[0], "HOH" if r[1] == "DOD" else r[1]) + r[2:]
    for r in records if r[4] not in ("H", "D"))
  assert heavy(atom_records(fn_out)) == heavy(atom_records(fn_in))

# ------------------------------------------------------------------------------
# step 4.2: the broad tests

def written_kind(option, resseq, name):
  """what option writes at a dipeptide_water site ("What each mode writes"):
  "split", "H" or "D"; all its parents are in the blank group"""
  if resseq == "101":
    return "H" if option == "all_h" else "D"
  if option in ("all_h", "all_d"):
    return option[-1].upper()
  if option == "all_hd" or (resseq, name) in dipeptide_exchangeable:
    return "split"
  return "H" if option == "hd_and_h" else "D"

def expected_hd(src, option):
  """{(resseq, name, altloc, element): (occupancy, B)} that option writes
  for the H/D of src, per "Restoration policy" (all parents at 1.0):
  recorded pair as pair -> both restored; recorded pair as one isotope ->
  parent occupancy, that isotope's B; new pair -> 0.5 each, B of the input
  atom; single without record -> unchanged"""
  sites = {}
  for (resseq, name, altloc, element), (occ, b, xyz) in \
      all_atoms(src, is_file=False).items():
    if element not in ("H", "D"): continue
    sites.setdefault((resseq, h_name(name, element)), {})[element] = (occ, b)
  result = {}
  for (resseq, name), values in sites.items():
    kind = written_kind(option, resseq, name)
    d_name = name.replace("H", "D", 1)
    recorded = len(values) == 2
    if kind == "split":
      if recorded:
        h, d = values["H"], values["D"]
      else:
        (occ, b), = values.values()
        h = d = (0.5, b)
      result[(resseq, name, "A", "H")] = h
      result[(resseq, d_name, "B", "D")] = d
    else:
      if recorded:
        v = (1.0, values[kind][1])
      else:
        (v,) = values.values()
      result[(resseq, name if kind == "H" else d_name, "", kind)] = v
  return result

def hg_pair_03_07(src):
  """src with the Ser HG pair at H 0.3 B 15 / D 0.7 B 25"""
  def edit(a):
    if a.parent().parent().resseq.strip() != "37": return
    if a.name.strip() == "HG" and a.parent().altloc == "A":
      a.set_occ(0.3); a.set_b(15.0)
    if a.name.strip() == "DG" and a.parent().altloc == "B":
      a.set_occ(0.7); a.set_b(25.0)
  return edit_atoms(src, edit)

# ------------------------------------------------------------------------------
# optimize and the input check: residues without restraints, CCD-only residues

_ligand_src = {}
def protein_ligand_src(option):
  """1yjp plus ligand 04 (chain B) from add with the CIF, for option"""
  if option not in _ligand_src:
    lig = heavy_pdb("04", zq4_atoms, offset=(15.0, 10.0, 15.0),
      cryst1=False).replace(" A   1 ", " B   1 ")
    prog = run_prog(PREFIX + "_lig_src_" + option,
      pdb_1yjp.replace("END\n", "") + lig + "END\n", zq4_cif,
      neutron_option=option)
    with open(prog.params.output.filename) as f:
      _ligand_src[option] = f.read()
  return _ligand_src[option]

def ligand_split_records(pdb_str, ligand):
  """(label, element, xyz, occ, b) of the ligand's atoms or the others"""
  h = iotbx.pdb.input(source_info=None, lines=pdb_str).construct_hierarchy()
  return sorted((a.id_str(), a.element.strip(), tuple(round(x, 3)
    for x in a.xyz), round(a.occ, 2), round(a.b, 2)) for a in h.atoms()
    if (a.parent().resname.strip() == "04") == ligand)

def bgm_atoms():
  """BGM (no monomer library entry, CCD only): heavy atoms at the CCD ideal
  coordinates, as heavy_pdb() atoms"""
  from mmtbx.chemical_components import get_cif_dictionary
  cc = get_cif_dictionary("BGM")
  return [(a.atom_id, a.type_symbol, "", float(a.pdbx_model_Cartn_x_ideal),
           float(a.pdbx_model_Cartn_y_ideal), float(a.pdbx_model_Cartn_z_ideal))
          for a in cc["_chem_comp_atom"] if a.type_symbol != "H"]

UNKNOWN_ATOMS = ("Atom names not in the restraint dictionary of their "
  "residue (dictionary older\nthan the CCD?) - ")
UNKNOWN_ATOMS_DO = ("Rename these atoms to match the dictionary, or supply a "
  "restraint CIF for ")

# ------------------------------------------------------------------------------
# flipkins: one input with every feature the flipkin passes must handle

def combined_flipkin_model():
  """1yjp (Asn/Gln amides, waters), His 45 (chain C), BGM (CCD only, chain
  B), 04 (chain D; its CIF is not given to the flipkin runs) and 03 (chain
  E, user CIF): heavy atoms in a P1 box"""
  def atom_lines(text):
    return "".join(l + "\n" for l in text.splitlines()
                   if l.startswith(("ATOM", "HETATM")))
  his = iotbx.pdb.input(source_info=None, lines=his_no_h).construct_hierarchy()
  his.atoms().set_xyz(his.atoms().extract_xyz() + (25.0, 0.0, 0.0))
  his.only_chain().id = "C"
  def ligand(code, atoms, offset, chain):
    return heavy_pdb(code, atoms, offset=offset, cryst1=False).replace(
      " A   1 ", " %s   1 " % chain)
  return ("CRYST1   80.000   80.000   80.000  90.00  90.00  90.00 P 1\n" +
    atom_lines(pdb_1yjp) + atom_lines(his.as_pdb_string()) +
    ligand("BGM", bgm_atoms(), (15.0, 10.0, 15.0), "B") +
    ligand("04", zq4_atoms, (-10.0, 10.0, -8.0), "D") +
    ligand("03", zq3_atoms, (-10.0, -10.0, 20.0), "E") + "END\n")

def run_prog_on_model(prefix, model, **params):
  """validate() and run() a reduce2.Program on a model in memory (no file)"""
  from iotbx.cli_parser import CCTBXParser
  params.setdefault("output__filename", prefix + "_out.pdb")
  params.setdefault("output__description_file_name", prefix + ".txt")
  params.setdefault("output__overwrite", True)
  parser = CCTBXParser(program_class=reduce2.Program, logger=null_out())
  parser.parse_args(["%s=%s" % (name.replace("__", "."), value)
                     for name, value in params.items()])
  parser.data_manager.add_model(prefix, model)
  log = StringIO()
  prog = reduce2.Program(parser.data_manager, parser.working_phil.extract(),
    logger=log)
  prog.validate()
  prog.run()
  prog.log_text = log.getvalue()
  return prog

def bgm_xh(file_name):
  """[(parent element, H or D, X-H length, CCD ideal X-H length)] of BGM"""
  from mmtbx.chemical_components import get_cif_dictionary
  cc = get_cif_dictionary("BGM")
  element = dict((a.atom_id, a.type_symbol) for a in cc["_chem_comp_atom"])
  ideal_xyz = dict((a.atom_id, matrix.col((float(a.pdbx_model_Cartn_x_ideal),
    float(a.pdbx_model_Cartn_y_ideal), float(a.pdbx_model_Cartn_z_ideal))))
    for a in cc["_chem_comp_atom"])
  parent = {}
  for b in cc["_chem_comp_bond"]:
    for x, y in ((b.atom_id_1, b.atom_id_2), (b.atom_id_2, b.atom_id_1)):
      if element[x] == "H": parent[x] = y
  h = iotbx.pdb.input(file_name).construct_hierarchy()
  atoms = [a for a in h.atoms() if a.parent().resname == "BGM"]
  heavy = dict((a.name.strip(), a) for a in atoms
               if a.element.strip() not in ("H", "D"))
  result = []
  for a in atoms:
    e = a.element.strip()
    if e not in ("H", "D"): continue
    name = a.name.strip()
    if e == "D": name = "H" + name[1:]
    x = parent[name]
    result.append((element[x], e, (matrix.col(a.xyz) -
      matrix.col(heavy[x].xyz)).length(), (ideal_xyz[name] -
      ideal_xyz[x]).length()))
  return result

def placed_text(pdb_str):
  """pdb_str with H placed by the library (no reduce2 normalization)"""
  model = mmtbx.model.manager(model_input=iotbx.pdb.input(source_info=None,
    lines=pdb_str), log=null_out())
  ph = reduce_hydrogen.place_hydrogens(model=model)
  ph.run()
  return ph.get_model().model_as_pdb()

def h_to_d(text, names=None):
  """H -> D (name and element), for all H or the given names"""
  result = []
  for l in text.splitlines(True):
    if l.startswith(("ATOM", "HETATM")) and l[76:78] == " H" and (
        names is None or l[12:16].strip() in names):
      i = l[12:16].index("H")
      l = l[:12+i] + "D" + l[13+i:76] + " D" + l[78:]
    result.append(l)
  return "".join(result)

def processed_model(text):
  model = mmtbx.model.manager(model_input=iotbx.pdb.input(source_info=None,
    lines=text), log=null_out())
  model.process(make_restraints=True, pdb_interpretation_params=
    reduce_hydrogen.get_reduce_pdb_interpretation_params(False))
  return model

def probe_setup(model):
  """atoms, bonded neighbor lists and extra atom info, as the Optimizer"""
  from mmtbx.probe import Helpers
  atoms = model.get_atoms()
  bond_proxies, asu_proxies = model.get_restraints_manager().geometry.\
    get_all_bond_proxies(sites_cart=model.get_sites_cart())
  bonded = Helpers.getBondedNeighborLists(atoms, bond_proxies, asu_proxies)
  extra = Helpers.getExtraAtomInfo(model=model, bondedNeighborLists=bonded,
    useNeutronDistances=False, probePhil=Optimizers._philLike()).extraAtomInfo
  return atoms, bonded, extra
