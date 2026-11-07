"""reduce2 add_h_to_water=True: the water H are placed by
mmtbx.hydrogens.water_protonation (virvel) once the other H are final, each
water oriented toward its own H-bond acceptors.
"""
import os
import iotbx.pdb
from libtbx.utils import null_out, Sorry
from libtbx.test_utils import approx_equal, Exception_expected
from iotbx.cli_parser import run_program
from mmtbx.programs import reduce2
from mmtbx.hydrogens import reduce_hydrogen
from mmtbx.hydrogens.water_protonation import place_water_hydrogens
from scitbx import matrix

# Each water is 2.8 A from one carbonyl O, along its C=O axis, and has no
# other acceptor within 3.5 A.
pdb_str = """\
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
ATOM     19  CE  LYS A   3      19.919  -1.002  16.894  1.00 20.00           C
ATOM     20  NZ  LYS A   3      19.531  -1.812  18.059  1.00 20.00           N
ATOM     21  OXT LYS A   3      20.995   0.451  11.763  1.00 20.00           O
HETATM   22  O   HOH A   4      12.729  -2.058  11.491  1.00 20.00           O
HETATM   23  O   HOH A   5      20.007  -3.513   9.351  1.00 20.00           O
END
"""

# Only waters and a Zn: placement adds no H, so the Optimizer has nothing to do.
pdb_waters_only = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
HETATM    1 ZN    ZN A   1      10.000  10.000  10.000  1.00 20.00          ZN
HETATM    2  O   HOH A   2      12.100  10.000  10.000  1.00 20.00           O
HETATM    3  O   HOH A   3      10.000  12.100  10.000  1.00 20.00           O
HETATM    4  O   HOH A   4      14.900  10.000  10.000  1.00 20.00           O
END
"""

def run_reduce2(prefix, extra_args=(), model_str=pdb_str):
  """Run reduce2 on model_str; returns the output hierarchy."""
  fn_in, fn_out = prefix + ".pdb", prefix + "_out.pdb"
  with open(fn_in, "w") as f:
    f.write(model_str)
  if os.path.exists(fn_out): os.remove(fn_out)
  run_program(program_class=reduce2.Program, logger=null_out(),
    args=[fn_in, "output.filename=" + fn_out, "overwrite=True",
          "output.description_file=" + prefix + ".txt"] + list(extra_args))
  return iotbx.pdb.input(fn_out).construct_hierarchy()

def xyz_by_id(hierarchy, water):
  """Coordinates by (resseq, atom name), of the water or of the other atoms."""
  sel = hierarchy.atom_selection_cache().selection("water")
  if not water: sel = ~sel
  return {(a.parent().parent().resseq.strip(), a.name.strip()): a.xyz
          for a in hierarchy.select(sel).atoms()}

def exercise_parity():
  """
  The water H of add_h_to_water=True are those that place_water_hydrogens
  gives on reduce2's output without them, and nothing else changes: the
  Optimizer sees the same bare waters either way.
  """
  with_h = run_reduce2("tst_reduce2_water_h_parity_on",
    ["add_h_to_water=True"])
  without = run_reduce2("tst_reduce2_water_h_parity_off")
  other_on, other_off = xyz_by_id(with_h, False), xyz_by_id(without, False)
  assert sorted(other_on) == sorted(other_off), set(other_on) ^ set(other_off)
  for key in other_off:
    assert approx_equal(other_on[key], other_off[key], eps=1e-3), key
  cs = iotbx.pdb.input(source_info=None,
    lines=pdb_str.split("\n")).crystal_symmetry()
  place_water_hydrogens(without, geometry="xray", element="H",
    existing_h="reorient",
    crystal_symmetry=reduce_hydrogen.contact_symmetry(cs))
  got, expected = xyz_by_id(with_h, True), xyz_by_id(without, True)
  assert sorted(got) == sorted(expected), (sorted(got), sorted(expected))
  for key in expected:
    assert approx_equal(got[key], expected[key], eps=2e-3), key

def check_toward_acceptors(h):
  """
  Each water has an O-H within 25 deg of its acceptor; returns the H1
  directions.
  """
  acceptor = {"4": ("1", "O"), "5": ("3", "O")}
  xyz = xyz_by_id(h, False)
  xyz.update(xyz_by_id(h, True))
  directions = []
  for resseq, acc in sorted(acceptor.items()):
    o = matrix.col(xyz[(resseq, "O")])
    to_acc = (matrix.col(xyz[acc]) - o).normalize()
    angles = []
    for name in ("H1", "H2"):
      oh = (matrix.col(xyz[(resseq, name)]) - o).normalize()
      angles.append(oh.angle(to_acc, deg=True))
    assert min(angles) < 25, (resseq, angles)
    directions.append(
      (matrix.col(xyz[(resseq, "H1")]) - o).normalize())
  return directions

def exercise_oriented_to_acceptors():
  """
  #1150: two waters with their acceptors in different directions get
  different orientations, each with an O-H within 25 deg of its acceptor.
  """
  h = run_reduce2("tst_reduce2_water_h_oriented", ["add_h_to_water=True"])
  directions = check_toward_acceptors(h)
  assert directions[0].angle(directions[1], deg=True) > 25, directions

def exercise_waters_only():
  """
  A model of waters and a Zn gets its water H with approach=add (no H for
  the Optimizer: its section is skipped), also with the atom dump and with a
  neutron_option, and with approach=place.
  """
  for i, args in enumerate((["approach=add"],
      ["approach=add", "output.print_atom_info=True"],
      ["approach=add", "neutron_option=all_d"],
      ["approach=place"])):
    prefix = f"tst_reduce2_water_h_only_{i}"
    h = run_reduce2(prefix, args + ["add_h_to_water=True"], pdb_waters_only)
    n = h.atom_selection_cache().selection(
      "water and (element H or element D)").count(True)
    assert n == 6, (args, n)
    if "approach=add" in args:
      with open(prefix + ".txt") as f:
        text = f.read()
      assert "No H to optimize" in text, (args, text)

def water_h_of(hierarchy):
  """Water H/D coordinates by (model index, resseq, name)."""
  sel = hierarchy.atom_selection_cache().selection(
    "water and (element H or element D)")
  result = {}
  for i, m in enumerate(hierarchy.select(sel).models()):
    for a in m.atoms():
      result[(i, a.parent().parent().resseq.strip(), a.name.strip())] = a.xyz
  return result

def exercise_models():
  """
  Two identical MODELs: the placer takes one model at a time, and both get
  the same water H.
  """
  inp = iotbx.pdb.input(source_info=None, lines=pdb_str.split("\n"))
  h = inp.construct_hierarchy()
  copy = h.only_model().detached_copy()
  h.only_model().id = "1"
  copy.id = "2"
  h.append_model(copy)
  two = h.as_pdb_string(crystal_symmetry=inp.crystal_symmetry())
  h = run_reduce2("tst_reduce2_water_h_models", ["add_h_to_water=True"], two)
  assert h.models_size() == 2, h.models_size()
  xyz = water_h_of(h)
  assert len(xyz) == 8, sorted(xyz)
  for (i, resseq, name), site in xyz.items():
    if i == 1:
      assert approx_equal(site, xyz[(0, resseq, name)], eps=1e-3), (resseq,
        name)

def exercise_kept_water_h():
  """
  keep_existing_H=True with water H in the input, both turned away from the
  acceptor: they are placed again (existing_h="reorient"), toward it, with
  their names.
  """
  with_h = pdb_str.replace("END\n",
    "HETATM   24  H1  HOH A   4      12.729  -2.058  12.341  1.00 20.00"
    "           H\n"
    "HETATM   25  H2  HOH A   4      13.553  -2.058  11.286  1.00 20.00"
    "           H\nEND\n")
  h = run_reduce2("tst_reduce2_water_h_kept",
    ["add_h_to_water=True", "keep_existing_H=True"], with_h)
  sel = h.atom_selection_cache().selection("resseq 4")
  names = sorted(a.name.strip() for a in h.select(sel).atoms())
  assert names == ["H1", "H2", "O"], names
  xyz = xyz_by_id(h, True)
  xyz.update(xyz_by_id(h, False))
  o = matrix.col(xyz[("4", "O")])
  to_acc = (matrix.col(xyz[("1", "O")]) - o).normalize()
  angles = [(matrix.col(xyz[("4", n)]) - o).normalize().angle(to_acc,
    deg=True) for n in ("H1", "H2")]
  assert min(angles) < 25, angles

def exercise_no_crystal_symmetry():
  """Without CRYST1 the water H are placed without symmetry, and say so."""
  no_cryst1 = "\n".join(l for l in pdb_str.splitlines()
    if not l.startswith("CRYST1")) + "\n"
  h = run_reduce2("tst_reduce2_water_h_nocs", ["add_h_to_water=True"],
    no_cryst1)
  n = len(water_h_of(h))
  assert n == 4, n
  with open("tst_reduce2_water_h_nocs.txt") as f:
    text = f.read()
  assert "without crystal symmetry" in text, text

def exercise_qmready():
  """qmready=True places water H too, each water toward its acceptor."""
  h = run_reduce2("tst_reduce2_water_h_qm", ["qmready=True"])
  n = len(water_h_of(h))
  assert n == 4, n
  check_toward_acceptors(h)

def exercise_optimize():
  """
  approach=optimize on a model with H everywhere but on the waters:
  add_h_to_water=True places the water H after the optimization, toward
  their acceptors; without it, and with qmready=True alone, the waters stay
  without H.
  """
  run_reduce2("tst_reduce2_water_h_opt_in")
  with open("tst_reduce2_water_h_opt_in_out.pdb") as f:
    with_h = f.read()
  h = run_reduce2("tst_reduce2_water_h_opt_on",
    ["approach=optimize", "add_h_to_water=True"], with_h)
  n = len(water_h_of(h))
  assert n == 4, n
  check_toward_acceptors(h)
  for args in (["approach=optimize"], ["approach=optimize", "qmready=True"]):
    h = run_reduce2("tst_reduce2_water_h_opt_off", args, with_h)
    assert water_h_of(h) == {}, args

def exercise_remove_sorry():
  """add_h_to_water=True with approach=remove is a Sorry."""
  try:
    run_reduce2("tst_reduce2_water_h_remove",
      ["approach=remove", "add_h_to_water=True"])
  except Sorry as e:
    assert "add_h_to_water=True cannot be used with approach=remove" in \
      str(e), str(e)
  else:
    raise Exception_expected

if __name__ == "__main__":
  exercise_parity()
  exercise_oriented_to_acceptors()
  exercise_waters_only()
  exercise_models()
  exercise_kept_water_h()
  exercise_no_crystal_symmetry()
  exercise_qmready()
  exercise_optimize()
  exercise_remove_sorry()
  print("OK")
