"""reduce2 reports the IUPAC renaming of symmetric side-chain atoms
(flip_symmetric_amino_acids in pdb_interpretation): the description lists every
renamed residue and atom pair, the log a capped list. Placement renames (add,
place, the flipkin passes); optimize and remove keep the input names, also
with a restraint file: reduce2's processing of the input for a restraint file
runs without the flip (B23), so a restraint file never changes the names.
"""
from __future__ import absolute_import, division, print_function
import os
from libtbx.utils import null_out
from six.moves import cStringIO as StringIO
from iotbx.cli_parser import CCTBXParser
from iotbx.pdb import hierarchy
from mmtbx.monomer_library import idealized_aa
from mmtbx.programs import reduce2
from mmtbx.regression.reduce2_neutron_utils import zq4_cif
from iotbx.pdb.tst_hierarchy_flip_symmetric import pdb_1us0_phe_276
import iotbx.pdb
import mmtbx.model
from mmtbx.hydrogens import reduce_hydrogen
from scitbx import matrix

HEADER = "IUPAC renaming of symmetric side-chain atoms"

# atoms the 180 deg turn about (first, second) moves
TURN = {
  "ASP": ("CB", "CG", ("OD1", "OD2", "HD2")),
  "GLU": ("CG", "CD", ("OE1", "OE2", "HE2")),
  "ARG": ("NE", "CZ", ("NH1", "NH2", "HH11", "HH12", "HH21", "HH22")),
  "PHE": ("CB", "CG", ("CD1", "CD2", "CE1", "CE2", "CZ", "HD1", "HD2", "HE1",
                       "HE2", "HZ")),
  "TYR": ("CB", "CG", ("CD1", "CD2", "CE1", "CE2", "CZ", "OH", "HD1", "HD2",
                       "HE1", "HE2", "HH")),
}

def build(spec, altlocs=False):
  """spec: (idealized_aa key, resseq, turn the side-chain end 180 deg), 25 A
  apart. altlocs: each residue as conformers A and B, the second 0.1 A away."""
  residues = idealized_aa.residue_dict()
  root = hierarchy.root()
  model = hierarchy.model()
  root.append_model(model)
  chain = hierarchy.chain(id="A")
  model.append_chain(chain)
  for i, (key, resseq, turned) in enumerate(spec):
    rg = residues[key].only_model().only_chain().only_residue_group(
      ).detached_copy()
    rg.resseq = "%4d" % resseq
    ag = rg.only_atom_group()
    if turned:
      first, second, names = TURN[ag.resname]
      origin = matrix.col(ag.get_atom(first).xyz)
      axis = (matrix.col(ag.get_atom(second).xyz) - origin).normalize()
      r = axis.axis_and_angle_as_r3_rotation_matrix(180, deg=True)
      for n in names:
        a = ag.get_atom(n)
        if a is not None: a.xyz = (r * (matrix.col(a.xyz) - origin) + origin).elems
    for a in ag.atoms():
      a.xyz = (a.xyz[0] + 10 + 25 * i, a.xyz[1] + 10, a.xyz[2] + 10)
      a.b = 20
      a.occ = 1
    if altlocs:
      other = ag.detached_copy()
      ag.altloc = "A"
      other.altloc = "B"
      for a in other.atoms():
        a.xyz = (a.xyz[0] + 0.1, a.xyz[1], a.xyz[2])
      for a in list(ag.atoms()) + list(other.atoms()):
        a.occ = 0.5
      rg.append_atom_group(other)
    chain.append_residue_group(rg)
  root.atoms().reset_serial()
  return root

def expected_renamed(root):
  """resseq of the residues flip_symmetric_amino_acids() itself renames"""
  flipped = root.deep_copy()
  flipped.flip_symmetric_amino_acids()
  result = set()
  for rg, other in zip(root.residue_groups(), flipped.residue_groups()):
    if list(rg.atoms().extract_xyz()) != list(other.atoms().extract_xyz()):
      result.add(rg.resseq.strip())
  return result

# heavy atoms; Glu 11 and Phe 51 get renamed, Glu 10 and Phe 50 do not
HEAVY = [("glu", 10, True), ("glu", 11, False), ("phe", 50, False),
         ("phe", 51, True)]

def write(prefix, root):
  fn = prefix + "_in.pdb"
  with open(fn, "w") as f: f.write(root.as_pdb_string())
  return fn

def run(prefix, fn, args, extra=()):
  log = StringIO()
  parser = CCTBXParser(program_class=reduce2.Program, logger=null_out())
  parser.parse_args([fn] + list(extra) + ["output.filename=%s.pdb" % prefix,
    "output.description_file_name=%s.txt" % prefix, "output.overwrite=True"] +
    args)
  task = reduce2.Program(parser.data_manager, parser.working_phil.extract(),
    master_phil=parser.master_phil, logger=log)
  task.validate()
  task.run()
  with open(prefix + ".txt") as f: description = f.read()
  return task, description, log.getvalue()

def listed(text):
  """The lines of the first renaming section: {'GLU A 11': 'OE1/OE2'}; None if
  the text has no such section."""
  lines = text.splitlines()
  start = [i for i, l in enumerate(lines) if l.startswith(HEADER)]
  if not start: return None
  result = {}
  for l in lines[start[0] + 1:]:
    if not l.startswith("  ") or ":" not in l: break
    label, pairs = l.strip().split(": ")
    result[label] = pairs
  return result

def sites(fn):
  h = iotbx.pdb.input(fn).construct_hierarchy()
  return dict(((rg.resseq.strip(), ag.altloc, a.name.strip()), a.xyz)
    for rg in h.residue_groups() for ag in rg.atom_groups() for a in ag.atoms())

def distance(a, b):
  return (matrix.col(a) - matrix.col(b)).length()

GLU_PHE = {"GLU A 11": "OE1/OE2", "PHE A 51": "CD1/CD2, CE1/CE2"}

def without_acid_h(root):
  for ag in root.atom_groups():
    for name in ("HD2", "HE2"):
      atom = ag.get_atom(name)
      if atom is not None and ag.resname in ("ASP", "GLU"):
        ag.remove_atom(atom)
  root.atoms().reset_serial()
  return root

def swapped_resseqs(renames):
  """resseq of each OD1/OD2 or OE1/OE2 swap in a list of rename records"""
  result = []
  for r in renames:
    if ("OD1", "OD2") in r.pairs or ("OD2", "OD1") in r.pairs:
      result.append(r.label.split()[-1])
    if ("OE1", "OE2") in r.pairs or ("OE2", "OE1") in r.pairs:
      result.append(r.label.split()[-1])
  return result

def test_add_and_place():
  # exactly the renamed ones, with the atom pairs; the atoms stay where they
  # were, under the partner's name
  prefix = "tst_reduce2_report_renaming_add"
  root = build(HEAVY)
  assert expected_renamed(root) == set(["11", "51"])
  fn = write(prefix, root)
  for approach in ("add", "place"):
    out = prefix + "_" + approach
    task, description, log = run(out, fn, ["approach=" + approach])
    assert listed(description) == GLU_PHE, description
    assert "%s: 2 residue(s)" % HEADER in description
    assert listed(log) == GLU_PHE, log
    before, after = sites(fn), sites(out + ".pdb")
    assert before["11", "", "OE1"] == after["11", "", "OE2"]
    assert before["11", "", "OE2"] == after["11", "", "OE1"]
    assert before["51", "", "CE1"] == after["51", "", "CE2"]
    assert before["10", "", "OE1"] == after["10", "", "OE1"]
    assert before["50", "", "CE1"] == after["50", "", "CE1"]

def test_add_nothing_renamed():
  prefix = "tst_reduce2_report_renaming_none"
  fn = write(prefix, build([("glu", 10, True), ("phe", 50, False)]))
  task, description, log = run(prefix, fn, ["approach=add"])
  assert listed(description) == {}, description
  assert "%s: none" % HEADER in description
  assert "%s: none" % HEADER in log

# all H supplied. Renamed: Tyr 31, Arg 40, Phe 50. Glu 11 would be, but its HE2
# has no partner: an acid with one H stays
SUPPLIED = [("glu_h", 10, True), ("glu_h", 11, False), ("tyr_h", 30, True),
            ("tyr_h", 31, False), ("arg_h", 40, True), ("arg_h", 41, False),
            ("phe_h", 50, True), ("phe_h", 51, False)]
SUPPLIED_LIST = {"TYR A 31": "CD1/CD2, CE1/CE2, HD1/HD2, HE1/HE2",
  "ARG A 40": "NH1/NH2, HH11/HH21, HH12/HH22",
  "PHE A 50": "CD1/CD2, CE1/CE2, HD1/HD2, HE1/HE2"}

def test_optimize():
  # the input stays as it is: no renaming, and the description says so
  prefix = "tst_reduce2_report_renaming_opt"
  root = build(SUPPLIED)
  assert expected_renamed(root) == set(["31", "40", "50"])
  fn = write(prefix, root)
  task, description, log = run(prefix, fn, ["approach=optimize"])
  assert listed(description) == {}, description
  assert "%s: none" % HEADER in description
  before, after = sites(fn), sites(prefix + ".pdb")
  assert before["40", "", "NH1"] == after["40", "", "NH1"]

def turned(spec):
  """build() with the side-chain end of a residue turned by any angle:
  spec (idealized_aa key, resseq, degrees about the TURN axis or 0)."""
  root = build([(key, resseq, False) for key, resseq, degrees in spec])
  for rg, (key, resseq, degrees) in zip(root.residue_groups(), spec):
    if not degrees: continue
    ag = rg.only_atom_group()
    first, second, names = TURN[ag.resname]
    origin = matrix.col(ag.get_atom(first).xyz)
    axis = (matrix.col(ag.get_atom(second).xyz) - origin).normalize()
    r = axis.axis_and_angle_as_r3_rotation_matrix(degrees, deg=True)
    for n in names:
      a = ag.get_atom(n)
      if a is not None: a.xyz = (r * (matrix.col(a.xyz) - origin) + origin).elems
  return root

# Phe 60 turned 220 deg (chi2 133), Asp 61 by 60 deg (chi2 150), both named
# against IUPAC; Phe 62 as IUPAC wants. All H supplied but the acid H.
NON_IUPAC = [("phe_h", 60, 220), ("asp_h", 61, 60), ("phe_h", 62, 0)]
NON_IUPAC_LIST = {"PHE A 60": "CD1/CD2, CE1/CE2", "ASP A 61": "OD1/OD2"}

def test_restraint_file():
  # B23: an unrelated restraint file gives the output of the run without it,
  # for every approach. add and place rename once (placement), optimize and
  # remove keep the input names.
  prefix = "tst_reduce2_report_renaming_cif"
  root = without_acid_h(turned(NON_IUPAC))
  assert expected_renamed(root) == set(["60", "61"])
  fn = write(prefix, root)
  cif = prefix + ".cif"
  with open(cif, "w") as f: f.write(zq4_cif)
  for approach, renamed in (("add", True), ("place", True),
                            ("optimize", False), ("remove", False)):
    out = prefix + "_" + approach
    task, with_cif, log = run(out + "_cif", fn, ["approach=" + approach],
      extra=[cif])
    task, without, log = run(out, fn, ["approach=" + approach])
    with open(out + "_cif.pdb") as f: a = f.read()
    with open(out + ".pdb") as f: b = f.read()
    assert a == b, approach
    assert listed(with_cif) == listed(without), (approach, with_cif)
    assert task._input_renames == []
    before, after = sites(fn), sites(out + "_cif.pdb")
    if renamed:
      assert listed(with_cif) == NON_IUPAC_LIST, (approach, with_cif)
      assert before["61", "", "OD1"] == after["61", "", "OD2"]
      assert before["60", "", "CE1"] == after["60", "", "CE2"]
    else:
      assert HEADER not in with_cif or listed(with_cif) == {}, with_cif
      assert before["61", "", "OD1"] == after["61", "", "OD1"]
      assert before["60", "", "CE1"] == after["60", "", "CE1"]
    assert before["62", "", "CE1"] == after["62", "", "CE1"]

# 8afw Phe C 298 (chi2 -90.8 deg) and every residue within 4 A of it: B20's
# double flip with a restraint file (renamed in the input's processing, back
# in placement's)
pdb_8afw_phe_298 = """\
CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1
ATOM   6321  N   ALA C 252     326.798 341.856 233.595  1.00 76.44           N
ATOM   6322  CA  ALA C 252     325.720 342.833 233.539  1.00 76.44           C
ATOM   6323  C   ALA C 252     325.984 343.836 232.425  1.00 76.44           C
ATOM   6324  O   ALA C 252     327.129 344.052 232.023  1.00 76.44           O
ATOM   6325  CB  ALA C 252     325.563 343.572 234.869  1.00 76.44           C
ATOM   6326  N   ILE C 253     324.909 344.456 231.937  1.00 76.78           N
ATOM   6327  CA  ILE C 253     324.962 345.454 230.876  1.00 76.78           C
ATOM   6328  C   ILE C 253     324.284 346.726 231.366  1.00 76.78           C
ATOM   6329  O   ILE C 253     323.256 346.662 232.048  1.00 76.78           O
ATOM   6330  CB  ILE C 253     324.288 344.960 229.581  1.00 76.78           C
ATOM   6331  CG1 ILE C 253     324.937 343.674 229.087  1.00 76.78           C
ATOM   6332  CG2 ILE C 253     324.368 346.010 228.495  1.00 76.78           C
ATOM   6333  CD1 ILE C 253     324.188 343.040 227.958  1.00 76.78           C
ATOM   6334  N   SER C 254     324.855 347.874 231.014  1.00 80.14           N
ATOM   6335  CA  SER C 254     324.256 349.147 231.381  1.00 80.14           C
ATOM   6336  C   SER C 254     322.879 349.296 230.745  1.00 80.14           C
ATOM   6337  O   SER C 254     322.483 348.529 229.864  1.00 80.14           O
ATOM   6338  CB  SER C 254     325.151 350.310 230.957  1.00 80.14           C
ATOM   6339  OG  SER C 254     325.626 350.125 229.639  1.00 80.14           O
ATOM   6386  N   ALA C 261     330.157 344.442 232.133  1.00 76.55           N
ATOM   6387  CA  ALA C 261     331.072 344.140 233.224  1.00 76.55           C
ATOM   6388  C   ALA C 261     330.959 342.670 233.599  1.00 76.55           C
ATOM   6389  O   ALA C 261     329.855 342.153 233.791  1.00 76.55           O
ATOM   6390  CB  ALA C 261     330.788 345.015 234.443  1.00 76.55           C
ATOM   6677  N   SER C 297     325.387 344.941 238.941  1.00 74.44           N
ATOM   6678  CA  SER C 297     324.351 345.804 238.386  1.00 74.44           C
ATOM   6679  C   SER C 297     324.886 347.217 238.222  1.00 74.44           C
ATOM   6680  O   SER C 297     325.585 347.725 239.103  1.00 74.44           O
ATOM   6681  CB  SER C 297     323.108 345.809 239.278  1.00 74.44           C
ATOM   6682  OG  SER C 297     323.384 346.423 240.524  1.00 74.44           O
ATOM   6683  N   PHE C 298     324.549 347.840 237.097  1.00 75.41           N
ATOM   6684  CA  PHE C 298     324.922 349.224 236.842  1.00 75.41           C
ATOM   6685  C   PHE C 298     323.887 350.168 237.432  1.00 75.41           C
ATOM   6686  O   PHE C 298     322.783 349.755 237.796  1.00 75.41           O
ATOM   6687  CB  PHE C 298     325.038 349.482 235.340  1.00 75.41           C
ATOM   6688  CG  PHE C 298     326.280 348.925 234.715  1.00 75.41           C
ATOM   6689  CD1 PHE C 298     327.427 349.692 234.626  1.00 75.41           C
ATOM   6690  CD2 PHE C 298     326.294 347.644 234.196  1.00 75.41           C
ATOM   6691  CE1 PHE C 298     328.565 349.187 234.045  1.00 75.41           C
ATOM   6692  CE2 PHE C 298     327.430 347.136 233.617  1.00 75.41           C
ATOM   6693  CZ  PHE C 298     328.567 347.908 233.541  1.00 75.41           C
ATOM   6694  N   SER C 299     324.249 351.443 237.523  1.00 80.83           N
ATOM   6695  CA  SER C 299     323.255 352.462 237.811  1.00 80.83           C
ATOM   6696  C   SER C 299     322.485 352.808 236.540  1.00 80.83           C
ATOM   6697  O   SER C 299     322.878 352.448 235.428  1.00 80.83           O
ATOM   6698  CB  SER C 299     323.915 353.710 238.388  1.00 80.83           C
ATOM   6699  OG  SER C 299     324.670 354.387 237.400  1.00 80.83           O
ATOM   6732  N   TYR C 304     327.810 354.379 238.427  1.00 73.40           N
ATOM   6733  CA  TYR C 304     328.401 353.480 239.403  1.00 73.40           C
ATOM   6734  C   TYR C 304     328.043 352.038 239.079  1.00 73.40           C
ATOM   6735  O   TYR C 304     327.061 351.760 238.387  1.00 73.40           O
ATOM   6736  CB  TYR C 304     327.945 353.814 240.831  1.00 73.40           C
ATOM   6737  CG  TYR C 304     327.938 355.289 241.147  1.00 73.40           C
ATOM   6738  CD1 TYR C 304     329.113 355.959 241.450  1.00 73.40           C
ATOM   6739  CD2 TYR C 304     326.756 356.009 241.158  1.00 73.40           C
ATOM   6740  CE1 TYR C 304     329.109 357.305 241.741  1.00 73.40           C
ATOM   6741  CE2 TYR C 304     326.742 357.356 241.451  1.00 73.40           C
ATOM   6742  CZ  TYR C 304     327.922 357.998 241.741  1.00 73.40           C
ATOM   6743  OH  TYR C 304     327.916 359.342 242.033  1.00 73.40           O
ATOM   6744  N   LEU C 305     328.858 351.122 239.599  1.00 73.40           N
ATOM   6745  CA  LEU C 305     328.657 349.690 239.435  1.00 73.40           C
ATOM   6746  C   LEU C 305     328.623 349.035 240.808  1.00 73.40           C
ATOM   6747  O   LEU C 305     329.450 349.352 241.667  1.00 73.40           O
ATOM   6748  CB  LEU C 305     329.767 349.072 238.583  1.00 73.40           C
ATOM   6749  CG  LEU C 305     329.619 347.584 238.286  1.00 73.40           C
ATOM   6750  CD1 LEU C 305     328.369 347.354 237.473  1.00 73.40           C
ATOM   6751  CD2 LEU C 305     330.840 347.062 237.556  1.00 73.40           C
ATOM   6848  N   LEU C 319     332.273 354.582 237.959  1.00 76.75           N
ATOM   6849  CA  LEU C 319     332.630 354.513 236.547  1.00 76.75           C
ATOM   6850  C   LEU C 319     332.719 355.894 235.908  1.00 76.75           C
ATOM   6851  O   LEU C 319     332.726 356.005 234.680  1.00 76.75           O
ATOM   6852  CB  LEU C 319     331.634 353.640 235.788  1.00 76.75           C
ATOM   6853  CG  LEU C 319     331.571 352.173 236.208  1.00 76.75           C
ATOM   6854  CD1 LEU C 319     330.172 351.627 236.018  1.00 76.75           C
ATOM   6855  CD2 LEU C 319     332.564 351.362 235.408  1.00 76.75           C
TER
"""

def test_8afw_double_flip():
  # B20: with a restraint file the flip ran twice and near chi2 +-90 the second
  # call swapped back. Now (B23, and an idempotent flip since #1286) one
  # renaming, in placement, the same with and without the restraint file, and
  # a further flip renames nothing.
  prefix = "tst_reduce2_report_renaming_8afw"
  fn = prefix + "_in.pdb"
  with open(fn, "w") as f: f.write(pdb_8afw_phe_298)
  cif = prefix + ".cif"
  with open(cif, "w") as f: f.write(zq4_cif)
  task, with_cif, log = run(prefix + "_cif", fn, ["approach=place"], extra=[cif])
  task, without, log = run(prefix, fn, ["approach=place"])
  assert listed(with_cif) == listed(without), with_cif
  assert listed(with_cif).get("PHE C 298") == "CD1/CD2, CE1/CE2", with_cif
  for out in (prefix + "_cif.pdb", prefix + ".pdb"):
    h = iotbx.pdb.input(out).construct_hierarchy()
    before = list(h.atoms().extract_xyz())
    h.flip_symmetric_amino_acids()
    assert list(h.atoms().extract_xyz()) == before, out
  with open(prefix + "_cif.pdb") as f: a = f.read()
  with open(prefix + ".pdb") as f: b = f.read()
  assert a == b

def test_second_processing_renames_nothing():
  # Asp from idealized_aa has chi2 near 90 through both O, without acid H: the
  # flip used to swap it in each processing; now in neither. The output has the
  # input names.
  prefix = "tst_reduce2_report_renaming_cancel"
  root = without_acid_h(build([("asp_h", 21, False), ("glu_h", 11, False)]))
  fn = write(prefix, root)
  assert expected_renamed(root) == set(["11"])
  cif = prefix + ".cif"
  with open(cif, "w") as f: f.write(zq4_cif)
  task, description, log = run(prefix, fn,
    ["approach=place", "keep_existing_H=True"], extra=[cif])
  assert swapped_resseqs(task._input_renames).count("21") == 0
  assert swapped_resseqs(task._placement_renames).count("21") == 0
  assert listed(description) == {"GLU A 11": "OE1/OE2"}, description
  before, after = sites(fn), sites(prefix + ".pdb")
  assert before["21", "", "OD1"] == after["21", "", "OD1"]

def phe_276_matched(fn):
  # conformers A and B renamed together: each CD1, CE1 nearest its namesake
  h = sites(fn)
  for name, partner in (("CD1", "CD2"), ("CE1", "CE2")):
    same = distance(h["276", "A", name], h["276", "B", name])
    crossed = distance(h["276", "A", name], h["276", "B", partner])
    assert same < 0.5 and crossed > 2, (fn, name, same, crossed)

def test_corresponding_conformers_in_placement():
  # H sit on placeholder sites during placement; correspondence must not see
  # them. As the direct call: A and B flip together.
  prefix = "tst_reduce2_report_renaming_corresponding"
  fn = prefix + "_in.pdb"
  with open(fn, "w") as f: f.write(pdb_1us0_phe_276)
  model = mmtbx.model.manager(
    model_input = iotbx.pdb.input(file_name=fn), log = null_out())
  placed = reduce_hydrogen.place_hydrogens(model = model)
  placed.run()
  out = prefix + "_placed.pdb"
  with open(out, "w") as f: f.write(placed.get_model().model_as_pdb())
  phe_276_matched(out)
  task, description, log = run(prefix, fn, ["approach=add"])
  assert listed(description) == {
    "PHE A 276 altloc A": "CD1/CD2, CE1/CE2",
    "PHE A 276 altloc B": "CD1/CD2, CE1/CE2"}, description
  phe_276_matched(prefix + ".pdb")

def test_two_renamings_cancel():
  # B has one H. The input's processing for a restraint file used to rename A
  # alone (B, incomplete, out of the set) and placement took A back. Without
  # the flip there (B23) neither processing renames: the input names, nothing
  # listed, as without the restraint file.
  prefix = "tst_reduce2_report_renaming_two_cancel"
  occ = {"A": "0.38", "B": "0.62"}
  lines = []
  for l in pdb_1us0_phe_276.splitlines():
    lines.append(l[:54] + "  " + occ[l[16]] + l[60:])
    if l[12:16] == " CD1" and l[16] == "B":
      lines.append("ATOM     99  HD1BPHE A 276    %8.3f" % (float(l[30:38]) + 1) +
        l[38:54] + "  0.62  8.00           H")
  fn = prefix + "_in.pdb"
  with open(fn, "w") as f: f.write("\n".join(lines) + "\n")
  cif = prefix + ".cif"
  with open(cif, "w") as f: f.write(zq4_cif)
  task, description, log = run(prefix, fn, ["approach=place"], extra=[cif])
  assert task._input_renames == [] and task._placement_renames == []
  assert listed(description) == {}, description
  before, after = sites(fn), sites(prefix + ".pdb")
  assert before["276", "A", "CD1"] == after["276", "A", "CD1"]
  task, description, log = run(prefix + "_no_cif", fn, ["approach=place"])
  with open(prefix + ".pdb") as f: a = f.read()
  with open(prefix + "_no_cif.pdb") as f: b = f.read()
  assert a == b

def test_neutron():
  # the same list as the run without a neutron_option
  prefix = "tst_reduce2_report_renaming_neutron"
  fn = write(prefix, build(HEAVY))
  for option in ("all_d", "hd_and_h"):
    out = prefix + "_" + option
    task, description, log = run(out, fn,
      ["approach=add", "neutron_option=" + option])
    assert listed(description) == GLU_PHE, (option, description)
  # optimize keeps the input names
  fn2 = write(prefix + "_supplied", build(SUPPLIED))
  out = prefix + "_optimize"
  task, description, log = run(out, fn2,
    ["approach=optimize", "neutron_option=all_d"])
  assert listed(description) == {}, description

def test_flipkin_passes():
  # each flipkin pass places H again and says so, in the log; optimize itself
  # renames nothing but its flipkin passes do
  prefix = "tst_reduce2_report_renaming_flipkin"
  spec = HEAVY + [("asn", 70, False)]
  fn = write(prefix, build(spec))
  folder = prefix + "_kin"
  if not os.path.isdir(folder): os.mkdir(folder)
  flip = ["add_flip_movers=True", "output.flipkin_directory=" + folder]
  task, description, log = run(prefix, fn, ["approach=add"] + flip)
  assert listed(description) == GLU_PHE, description
  assert log.count("%s: 2 residue(s)" % HEADER) == 3, log   # main + 2 passes
  out = prefix + "_opt"
  fn2 = write(out, build([("glu_h", 11, False), ("asn_h", 70, False)]))
  task, description, log = run(out, fn2, ["approach=optimize"] + flip)
  assert listed(description) == {}, description
  assert log.count("%s: none" % HEADER) == 1, log
  assert log.count("%s: 1 residue(s)" % HEADER) == 2, log

def test_altlocs_and_models():
  # one line per conformer, one residue in the count
  prefix = "tst_reduce2_report_renaming_alt"
  fn = write(prefix, build([("glu", 11, False), ("glu", 12, True)], altlocs=True))
  task, description, log = run(prefix, fn, ["approach=add"])
  assert listed(description) == {"GLU A 11 altloc A": "OE1/OE2",
    "GLU A 11 altloc B": "OE1/OE2"}, description
  assert "%s: 1 residue(s)" % HEADER in description
  # two models; with model_id only that one
  prefix = "tst_reduce2_report_renaming_models"
  root = build([("glu", 11, False)])
  second = root.only_model().detached_copy()
  root.only_model().id = "1"
  second.id = "2"
  rg2 = second.only_chain().only_residue_group()
  rg2.resseq = "  12"
  root.append_model(second)
  fn = write(prefix, root)
  task, description, log = run(prefix, fn, ["approach=place"])
  assert listed(description) == {"GLU A 11 (model 1)": "OE1/OE2",
    "GLU A 12 (model 2)": "OE1/OE2"}, description
  # with model_id only the selected model, also with a restraint file
  cif = prefix + ".cif"
  with open(cif, "w") as f: f.write(zq4_cif)
  out = prefix + "_id2"
  task, description, log = run(out, fn, ["approach=place", "model_id=2"],
    extra=[cif])
  assert listed(description) == {"GLU A 12": "OE1/OE2"}, description
  # remove keeps the names of every model (B23)
  out = prefix + "_remove"
  task, description, log = run(out, fn, ["approach=remove", "model_id=2"],
    extra=[cif])
  assert HEADER not in description, description

def test_insertion_codes():
  # same resseq, same pair: two residues, neither cancels the other
  prefix = "tst_reduce2_report_renaming_icode"
  root = build([("glu", 11, False), ("glu", 12, False)])
  rg = root.only_model().only_chain().residue_groups()[1]
  rg.resseq = "  11"
  rg.icode = "A"
  fn = write(prefix, root)
  task, description, log = run(prefix, fn, ["approach=place"])
  assert listed(description) == {"GLU A 11": "OE1/OE2",
    "GLU A 11A": "OE1/OE2"}, description
  assert "%s: 2 residue(s)" % HEADER in description
  assert listed(log) == listed(description), log

def test_model_number():
  # model_id is a number: MODEL 0002 is model_id=2
  prefix = "tst_reduce2_report_renaming_model_number"
  root = build([("glu", 11, False)])
  second = root.only_model().detached_copy()
  root.only_model().id = "1"
  second.id = "0002"
  second.only_chain().only_residue_group().resseq = "  12"
  root.append_model(second)
  fn = write(prefix, root)
  with open(fn) as fh:
    assert "MODEL     0002" in fh.read()
  cif = prefix + ".cif"
  with open(cif, "w") as f: f.write(zq4_cif)
  task, description, log = run(prefix, fn, ["approach=place", "model_id=2"],
    extra=[cif])
  assert listed(description) == {"GLU A 12": "OE1/OE2"}, description
  # Asp 22 is not renamed, in the model selected from the second place
  out = prefix + "_cancel"
  root = without_acid_h(build([("asp_h", 21, False), ("glu_h", 11, False)]))
  second = root.only_model().detached_copy()
  root.only_model().id = "1"
  second.id = "0002"
  second.only_chain().residue_groups()[0].resseq = "  22"
  second.only_chain().residue_groups()[1].resseq = "  12"
  root.append_model(second)
  fn = write(out, root)
  task, description, log = run(out, fn,
    ["approach=place", "keep_existing_H=True", "model_id=2"], extra=[cif])
  assert swapped_resseqs(task._input_renames).count("22") == 0
  assert swapped_resseqs(task._placement_renames).count("22") == 0
  assert listed(description) == {"GLU A 12": "OE1/OE2"}, description

def test_log_cap():
  # the log caps the list at 20 unless verbosity > 2; the description is full
  prefix = "tst_reduce2_report_renaming_cap"
  spec = [("glu", 100 + i, False) for i in range(22)]
  fn = write(prefix, build(spec))
  task, description, log = run(prefix, fn, ["approach=place"])
  assert len(listed(description)) == 22, description
  assert len(listed(log)) == 20, log
  assert "(+2 more)" in log
  assert "(+2 more)" not in description
  out = prefix + "_v"
  task, description, log = run(out, fn, ["approach=place", "verbosity=3"])
  assert len(listed(log)) == 22, log
  assert "(+2 more)" not in log

def test_supplied_h():
  # keep_existing_H: H named in a pair go with their heavy atom: same site,
  # partner's name, bonded to the atom of that name. Glu 11, whose test asks for
  # the renaming, has one acid H: not renamed, the H stays on its O.
  prefix = "tst_reduce2_report_renaming_supplied"
  fn = write(prefix, build(SUPPLIED))
  task, description, log = run(prefix, fn,
    ["approach=place", "keep_existing_H=True"])
  assert listed(description) == SUPPLIED_LIST, description
  before, after = sites(fn), sites(prefix + ".pdb")
  for resseq, h1, h2 in (("40", "HH11", "HH21"), ("40", "HH12", "HH22"),
                         ("50", "HD1", "HD2"), ("50", "HE1", "HE2"),
                         ("31", "HD1", "HD2"), ("31", "HE1", "HE2")):
    assert distance(before[resseq, "", h1], after[resseq, "", h2]) < 0.1
    assert distance(before[resseq, "", h2], after[resseq, "", h1]) < 0.1
  assert before["11", "", "OE1"] == after["11", "", "OE1"]
  assert before["11", "", "OE2"] == after["11", "", "OE2"]
  assert distance(before["11", "", "HE2"], after["11", "", "HE2"]) < 0.1
  for resseq, h, parent, ideal in (("40", "HH11", "NH1", 0.86),
      ("40", "HH22", "NH2", 0.86), ("50", "HD1", "CD1", 0.93),
      ("31", "HE2", "CE2", 0.93), ("31", "HH", "OH", 0.84),
      ("11", "HE2", "OE2", 0.95), ("10", "HE2", "OE2", 0.95)):
    d = distance(after[resseq, "", h], after[resseq, "", parent])
    assert abs(d - ideal) < 0.05, (resseq, h, parent, d)
  assert distance(before["31", "", "HH"], after["31", "", "HH"]) < 0.1

if __name__ == "__main__":
  test_add_and_place()
  test_add_nothing_renamed()
  test_optimize()
  test_restraint_file()
  test_8afw_double_flip()
  test_second_processing_renames_nothing()
  test_corresponding_conformers_in_placement()
  test_two_renamings_cancel()
  test_neutron()
  test_flipkin_passes()
  test_altlocs_and_models()
  test_insertion_codes()
  test_model_number()
  test_log_cap()
  test_supplied_h()
  print("OK")
