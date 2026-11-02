"""reduce2 reports the IUPAC renaming of symmetric side-chain atoms
(flip_symmetric_amino_acids in pdb_interpretation): the description lists every
renamed residue and atom pair, the log a capped list. Placement renames (add,
place, the flipkin passes), optimize does not unless a restraint file makes
reduce2 process the input; two renamings of a residue can cancel.
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
from scitbx import matrix

HEADER = "IUPAC renaming of symmetric side-chain atoms"

# atoms the 180 deg turn about (first, second) moves
TURN = {
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
  import iotbx.pdb
  h = iotbx.pdb.input(fn).construct_hierarchy()
  return dict(((rg.resseq.strip(), ag.altloc, a.name.strip()), a.xyz)
    for rg in h.residue_groups() for ag in rg.atom_groups() for a in ag.atoms())

def distance(a, b):
  return (matrix.col(a) - matrix.col(b)).length()

GLU_PHE = {"GLU A 11": "OE1/OE2", "PHE A 51": "CD1/CD2, CE1/CE2"}

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

# all H supplied. Renamed: Glu 11, Tyr 31, Arg 40, Phe 50
SUPPLIED = [("glu_h", 10, True), ("glu_h", 11, False), ("tyr_h", 30, True),
            ("tyr_h", 31, False), ("arg_h", 40, True), ("arg_h", 41, False),
            ("phe_h", 50, True), ("phe_h", 51, False)]
SUPPLIED_LIST = {"GLU A 11": "OE1/OE2",
  "TYR A 31": "CD1/CD2, CE1/CE2, HD1/HD2, HE1/HE2",
  "ARG A 40": "NH1/NH2, HH11/HH21, HH12/HH22",
  "PHE A 50": "CD1/CD2, CE1/CE2, HD1/HD2, HE1/HE2"}

def test_optimize():
  # the input stays as it is: no renaming, and the description says so
  prefix = "tst_reduce2_report_renaming_opt"
  root = build(SUPPLIED)
  assert expected_renamed(root) == set(["11", "31", "40", "50"])
  fn = write(prefix, root)
  task, description, log = run(prefix, fn, ["approach=optimize"])
  assert listed(description) == {}, description
  assert "%s: none" % HEADER in description
  before, after = sites(fn), sites(prefix + ".pdb")
  assert before["40", "", "NH1"] == after["40", "", "NH1"]

def test_restraint_file():
  # a restraint file: reduce2 processes the input first, which renames, also
  # for optimize and remove
  prefix = "tst_reduce2_report_renaming_cif"
  fn = write(prefix, build(SUPPLIED))
  cif = prefix + ".cif"
  with open(cif, "w") as f: f.write(zq4_cif)
  task, description, log = run(prefix, fn, ["approach=optimize"], extra=[cif])
  assert listed(description) == SUPPLIED_LIST, description
  before, after = sites(fn), sites(prefix + ".pdb")
  assert before["40", "", "NH1"] == after["40", "", "NH2"]
  out = prefix + "_remove"
  task, description, log = run(out, fn, ["approach=remove"], extra=[cif])
  assert listed(description) == {"GLU A 11": "OE1/OE2",
    "TYR A 31": "CD1/CD2, CE1/CE2, HD1/HD2, HE1/HE2",
    "ARG A 40": "NH1/NH2, HH11/HH21, HH12/HH22",
    "PHE A 50": "CD1/CD2, CE1/CE2, HD1/HD2, HE1/HE2"}, description
  # without a restraint file remove does not mention it
  out = prefix + "_remove_no_cif"
  task, description, log = run(out, fn, ["approach=remove"])
  assert HEADER not in description, description

def test_two_renamings_cancel():
  # Asp from idealized_aa has chi2 near 90: the second call of
  # flip_symmetric_amino_acids() swaps back; the output has the input names
  prefix = "tst_reduce2_report_renaming_cancel"
  fn = write(prefix, build([("asp_h", 21, False), ("glu_h", 11, False)]))
  assert expected_renamed(build([("asp_h", 21, False)])) == set(["21"])
  cif = prefix + ".cif"
  with open(cif, "w") as f: f.write(zq4_cif)
  task, description, log = run(prefix, fn,
    ["approach=place", "keep_existing_H=True"], extra=[cif])
  assert listed(description) == {"GLU A 11": "OE1/OE2"}, description
  before, after = sites(fn), sites(prefix + ".pdb")
  assert before["21", "", "OD1"] == after["21", "", "OD1"]

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
  # a restraint file renames all models before model_id selects one
  cif = prefix + ".cif"
  with open(cif, "w") as f: f.write(zq4_cif)
  out = prefix + "_id2"
  task, description, log = run(out, fn, ["approach=remove", "model_id=2"],
    extra=[cif])
  assert listed(description) == {"GLU A 12": "OE1/OE2"}, description

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
  assert "MODEL     0002" in open(fn).read()
  cif = prefix + ".cif"
  with open(cif, "w") as f: f.write(zq4_cif)
  task, description, log = run(prefix, fn, ["approach=remove", "model_id=2"],
    extra=[cif])
  assert listed(description) == {"GLU A 12": "OE1/OE2"}, description
  # the input's rename and placement's second swap of Asp 22 net to nothing,
  # though selection moved the model from the second place to the first
  out = prefix + "_cancel"
  root = build([("asp_h", 21, False), ("glu_h", 11, False)])
  second = root.only_model().detached_copy()
  root.only_model().id = "1"
  second.id = "0002"
  second.only_chain().residue_groups()[0].resseq = "  22"
  second.only_chain().residue_groups()[1].resseq = "  12"
  root.append_model(second)
  fn = write(out, root)
  task, description, log = run(out, fn,
    ["approach=place", "keep_existing_H=True", "model_id=2"], extra=[cif])
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
  # partner's name, bonded to the atom of that name. The renaming has no pair
  # for the acid H: the H keeps its name and ends on the O of that name.
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
  test_two_renamings_cancel()
  test_neutron()
  test_flipkin_passes()
  test_altlocs_and_models()
  test_insertion_codes()
  test_model_number()
  test_log_cap()
  test_supplied_h()
  print("OK")
