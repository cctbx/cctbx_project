"""reduce2 approach=place adds H atoms without optimizing them: the output is
that of reduce_hydrogen.place_hydrogens. Options that need the optimizer stop
with a Sorry.
"""
from __future__ import absolute_import, division, print_function
import os
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out, Sorry
from libtbx.test_utils import approx_equal
from iotbx.cli_parser import run_program
from mmtbx.programs import reduce2
from mmtbx.hydrogens import reduce_hydrogen

# Ser OG-H and Lys NZ-H3 are rotatable: approach=add moves them.
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
END
"""

def xyz_by_name(hierarchy):
  return {(a.parent().parent().resseq, a.name): a.xyz
          for a in hierarchy.atoms()}

def run_reduce2(prefix, extra_args=()):
  fn_in, fn_out = prefix + ".pdb", prefix + "_out.pdb"
  with open(fn_in, "w") as f:
    f.write(pdb_str)
  if os.path.exists(fn_out): os.remove(fn_out)
  run_program(program_class=reduce2.Program, logger=null_out(),
    args=[fn_in, "output.filename=" + fn_out, "overwrite=True",
          "output.description_file=" + prefix + ".txt"] + list(extra_args))
  return xyz_by_name(iotbx.pdb.input(fn_out).construct_hierarchy())

def placed_by_library():
  inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  model = mmtbx.model.manager(model_input=inp, log=null_out())
  o = reduce_hydrogen.place_hydrogens(model=model)
  o.run()
  return xyz_by_name(o.get_model().get_hierarchy())

_placed = []
def placed_by_reduce2():
  """the approach=place run, shared by the first two tests (read only)"""
  if not _placed:
    _placed.append(run_reduce2("tst_reduce2_place_1", ["approach=place"]))
  return _placed[0]

def test_place_is_place_hydrogens():
  placed = placed_by_reduce2()
  ref = placed_by_library()
  assert sorted(placed.keys()) == sorted(ref.keys())
  for key in ref:
    assert approx_equal(placed[key], ref[key], eps=2e-3), key

def test_add_optimizes():
  placed = placed_by_reduce2()
  added = run_reduce2("tst_reduce2_place_3", ["approach=add"])
  moved = [key for key in placed if key in added and
    max(abs(a-b) for a, b in zip(placed[key], added[key])) > 0.01]
  assert moved, "approach=add moved no H"

def test_place_rejects_flips():
  try:
    run_reduce2("tst_reduce2_place_4", ["approach=place", "add_flip_movers=True"])
  except Sorry as e:
    assert "add_flip_movers requires optimization" in str(e), str(e)
  else:
    raise AssertionError("approach=place add_flip_movers=True did not stop")

if __name__ == "__main__":
  test_place_is_place_hydrogens()
  test_add_optimizes()
  test_place_rejects_flips()
  print("OK")
