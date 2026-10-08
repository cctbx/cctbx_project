"""reduce2 and hydrogenate on a model without crystal symmetry.

A P1 box is made for the calculation only. It must not be written out, and
the library path (place_hydrogens) must not move the coordinates. A cryo-EM
placeholder cell (CRYST1 1 1 1 P 1) in the input is written back.
"""
from __future__ import absolute_import, division, print_function
import os
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out
from libtbx.test_utils import approx_equal
from iotbx.cli_parser import run_program
from mmtbx.programs import reduce2, hydrogenate
from mmtbx.hydrogens import reduce_hydrogen

# Far from the origin, as cryo-EM models usually are.
pdb_str = """\
ATOM      1  N   ALA A   1     111.204 102.315 150.312  1.00 20.00           N
ATOM      2  CA  ALA A   1     112.555 102.835 150.528  1.00 20.00           C
ATOM      3  C   ALA A   1     113.588 101.726 150.381  1.00 20.00           C
ATOM      4  O   ALA A   1     113.325 100.567 150.721  1.00 20.00           O
ATOM      5  CB  ALA A   1     112.866 103.965 149.556  1.00 20.00           C
ATOM      6  N   SER A   2     114.767 102.090 149.884  1.00 20.00           N
ATOM      7  CA  SER A   2     115.853 101.130 149.701  1.00 20.00           C
ATOM      8  C   SER A   2     116.968 101.370 150.713  1.00 20.00           C
ATOM      9  O   SER A   2     117.281 102.516 151.041  1.00 20.00           O
ATOM     10  CB  SER A   2     116.411 101.224 148.278  1.00 20.00           C
ATOM     11  OG  SER A   2     115.456 100.826 147.311  1.00 20.00           O
ATOM     12  N   LYS A   3     117.560 100.282 151.199  1.00 20.00           N
ATOM     13  CA  LYS A   3     118.643 100.370 152.172  1.00 20.00           C
ATOM     14  C   LYS A   3     119.945  99.852 151.570  1.00 20.00           C
ATOM     15  O   LYS A   3     119.964  98.824 150.892  1.00 20.00           O
ATOM     16  CB  LYS A   3     118.290  99.576 153.432  1.00 20.00           C
ATOM     17  CG  LYS A   3     119.308  99.699 154.556  1.00 20.00           C
ATOM     18  CD  LYS A   3     118.913  98.863 155.762  1.00 20.00           C
ATOM     19  CE  LYS A   3     119.919  98.998 156.894  1.00 20.00           C
ATOM     20  NZ  LYS A   3     119.531  98.188 158.059  1.00 20.00           N
ATOM     21  OXT LYS A   3     120.995 100.451 151.763  1.00 20.00           O
END
"""

cryst1 = "CRYST1   40.000   50.000   60.000  90.00  90.00  90.00 P 21 21 21\n"

# cryo-EM placeholder; cctbx reads it as no symmetry.
cryst1_placeholder = \
  "CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1\n"

# Acetate ~5 A from the peptide; with a restraint file the programs call
# process() before H placement, and process() boxes a model without symmetry.
act_str = """\
HETATM   22  C   ACT B   1     125.000 104.000 152.000  1.00 20.00           C
HETATM   23  O   ACT B   1     125.266 103.768 153.198  1.00 20.00           O
HETATM   24  OXT ACT B   1     125.812 104.059 151.055  1.00 20.00           O
HETATM   25  CH3 ACT B   1     123.517 104.263 151.678  1.00 20.00           C
"""

def act_cif():
  '''Path to the geostd ACT restraints, or None.'''
  import libtbx.load_env
  return libtbx.env.find_in_repositories(
    relative_path='chem_data/geostd/a/data_ACT.cif', test=os.path.isfile)

def pdb_with_act():
  return pdb_str.replace("END\n", "") + act_str + "END\n"

def placeholder_cif_str():
  from cctbx import crystal
  return iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None
    ).construct_hierarchy().as_mmcif_string(
      crystal_symmetry=crystal.symmetry((1, 1, 1, 90, 90, 90), "P 1"))

def assert_placeholder_pdb(txt):
  lines = [l for l in txt.split("\n") if l.startswith("CRYST1")]
  assert len(lines) == 1, lines
  assert lines[0].split()[1:8] == \
    ["1.000", "1.000", "1.000", "90.00", "90.00", "90.00", "P"], lines[0]

def assert_placeholder_cif(txt):
  cell = [l.split()[1] for l in txt.split("\n")
          if l.startswith("_cell.length_")]
  assert cell == ["1.000", "1.000", "1.000"], cell

def heavy_xyz(hierarchy):
  return {(a.parent().parent().resseq, a.name): a.xyz
          for a in hierarchy.atoms() if a.element.strip() not in ("H", "D")}

def all_xyz(hierarchy):
  return {(a.parent().parent().resseq, a.name): a.xyz
          for a in hierarchy.atoms()}

def n_h(hierarchy):
  return sum(1 for a in hierarchy.atoms() if a.element.strip() in ("H", "D"))

def run_cli(pdb_text, prefix, suffix, in_suffix=".pdb", extra_args=()):
  fn_in  = prefix + in_suffix
  fn_out = prefix + "_out" + suffix
  with open(fn_in, "w") as f:
    f.write(pdb_text)
  if os.path.exists(fn_out): os.remove(fn_out)
  args = [fn_in, "output.filename=" + fn_out, "overwrite=True",
          "output.description_file=" + prefix + ".txt"] + list(extra_args)
  results = run_program(program_class=reduce2.Program, logger=null_out(),
    args=args)
  with open(fn_out) as f:
    txt = f.read()
  return txt, results.model

def run_hydrogenate(model_text, prefix, suffix, extra_args=()):
  fn_in  = prefix + suffix
  fn_out = prefix + "_hydrogenate" + suffix
  with open(fn_in, "w") as f:
    f.write(model_text)
  if os.path.exists(fn_out): os.remove(fn_out)
  run_program(program_class=hydrogenate.Program, logger=null_out(),
    args=[fn_in, "output.prefix=" + prefix] + list(extra_args))
  with open(fn_out) as f:
    return f.read()

def place_lib(pdb_text):
  inp = iotbx.pdb.input(lines=pdb_text.split("\n"), source_info=None)
  model = mmtbx.model.manager(model_input=inp, log=null_out())
  assert model.crystal_symmetry() is None
  o = reduce_hydrogen.place_hydrogens(model=model)
  o.run()
  return o.get_model()

def test_cli_pdb_no_cryst1():
  txt, _ = run_cli(pdb_str, "tst_reduce2_no_symmetry_1", ".pdb")
  assert "CRYST1" not in txt, txt[:200]
  h = iotbx.pdb.input(lines=txt.split("\n"), source_info=None
    ).construct_hierarchy()
  assert n_h(h) > 0
  ref = heavy_xyz(iotbx.pdb.input(lines=pdb_str.split("\n"),
    source_info=None).construct_hierarchy())
  out = heavy_xyz(h)
  assert set(ref) == set(out)
  for k in ref:
    assert approx_equal(ref[k], out[k], eps=1.e-3), k

def test_cli_cif_no_cell():
  txt, _ = run_cli(pdb_str, "tst_reduce2_no_symmetry_2", ".cif")
  assert "_cell.length_a" not in txt
  assert "_symmetry.space_group_name_H-M" not in txt
  assert "_space_group.name_H-M_alt" not in txt

def test_cli_keeps_real_cryst1():
  txt, _ = run_cli(cryst1 + pdb_str, "tst_reduce2_no_symmetry_3", ".pdb")
  cs = iotbx.pdb.input(lines=txt.split("\n"), source_info=None
    ).crystal_symmetry()
  assert cs is not None
  assert approx_equal(cs.unit_cell().parameters(), (40, 50, 60, 90, 90, 90))
  assert str(cs.space_group_info()) == "P 21 21 21"

def test_cli_box_cushion():
  # Working box = heavy-atom extent + 2*5 A.
  _, model = run_cli(pdb_str, "tst_reduce2_no_symmetry_4", ".pdb")
  xyz = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None
    ).atoms().extract_xyz()
  extent = [mx - mn for mx, mn in zip(xyz.max(), xyz.min())]
  abc = model.crystal_symmetry().unit_cell().parameters()[:3]
  assert approx_equal(abc, [e + 10 for e in extent], eps=1.e-3), abc

def test_hydrogenate_pdb_no_cryst1():
  txt = run_hydrogenate(pdb_str, "tst_reduce2_no_symmetry_5", ".pdb")
  assert "CRYST1" not in txt, txt[:200]
  h = iotbx.pdb.input(lines=txt.split("\n"), source_info=None
    ).construct_hierarchy()
  assert n_h(h) > 0
  ref = heavy_xyz(iotbx.pdb.input(lines=pdb_str.split("\n"),
    source_info=None).construct_hierarchy())
  out = heavy_xyz(h)
  assert set(ref) == set(out)
  for k in ref:
    assert approx_equal(ref[k], out[k], eps=1.e-3), k

def test_hydrogenate_cif_no_cell():
  cif_str = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None
    ).construct_hierarchy().as_mmcif_string()
  txt = run_hydrogenate(cif_str, "tst_reduce2_no_symmetry_6", ".cif")
  assert "_atom_site" in txt
  assert "_cell.length_a" not in txt
  assert "_symmetry.space_group_name_H-M" not in txt
  assert "_space_group.name_H-M_alt" not in txt

def test_hydrogenate_keeps_real_cryst1():
  txt = run_hydrogenate(cryst1 + pdb_str, "tst_reduce2_no_symmetry_7", ".pdb")
  cs = iotbx.pdb.input(lines=txt.split("\n"), source_info=None
    ).crystal_symmetry()
  assert cs is not None
  assert approx_equal(cs.unit_cell().parameters(), (40, 50, 60, 90, 90, 90))
  assert str(cs.space_group_info()) == "P 21 21 21"

def test_cli_restraints_no_cryst1():
  fn_cif = act_cif()
  if fn_cif is None:
    print("  skipping test_cli_restraints_no_cryst1: geostd not available")
    return
  txt, _ = run_cli(pdb_with_act(), "tst_reduce2_no_symmetry_8", ".pdb",
    extra_args=[fn_cif])
  assert "CRYST1" not in txt, txt[:200]
  h = iotbx.pdb.input(lines=txt.split("\n"), source_info=None
    ).construct_hierarchy()
  assert "ACT" in h.overall_counts().resnames

def test_hydrogenate_restraints_no_cryst1():
  fn_cif = act_cif()
  if fn_cif is None:
    print("  skipping test_hydrogenate_restraints_no_cryst1: "
          "geostd not available")
    return
  txt = run_hydrogenate(pdb_with_act(), "tst_reduce2_no_symmetry_9", ".pdb",
    extra_args=[fn_cif])
  assert "CRYST1" not in txt, txt[:200]

def test_cli_keeps_placeholder_pdb():
  txt, _ = run_cli(cryst1_placeholder + pdb_str,
    "tst_reduce2_no_symmetry_10", ".pdb")
  assert_placeholder_pdb(txt)
  h = iotbx.pdb.input(lines=txt.split("\n"), source_info=None
    ).construct_hierarchy()
  ref = heavy_xyz(iotbx.pdb.input(lines=pdb_str.split("\n"),
    source_info=None).construct_hierarchy())
  out = heavy_xyz(h)
  for k in ref:
    assert approx_equal(ref[k], out[k], eps=1.e-3), k

def test_cli_keeps_placeholder_cif():
  txt, _ = run_cli(placeholder_cif_str(), "tst_reduce2_no_symmetry_11",
    ".cif", in_suffix=".cif")
  assert_placeholder_cif(txt)

def test_hydrogenate_keeps_placeholder():
  txt = run_hydrogenate(cryst1_placeholder + pdb_str,
    "tst_reduce2_no_symmetry_12", ".pdb")
  assert_placeholder_pdb(txt)
  txt = run_hydrogenate(placeholder_cif_str(),
    "tst_reduce2_no_symmetry_13", ".cif")
  assert_placeholder_cif(txt)

def test_placeholder_with_restraints():
  # Both the placeholder and the restraint-file path at once.
  fn_cif = act_cif()
  if fn_cif is None:
    print("  skipping test_placeholder_with_restraints: geostd not available")
    return
  txt, _ = run_cli(cryst1_placeholder + pdb_with_act(),
    "tst_reduce2_no_symmetry_14", ".pdb", extra_args=[fn_cif])
  assert_placeholder_pdb(txt)
  txt = run_hydrogenate(cryst1_placeholder + pdb_with_act(),
    "tst_reduce2_no_symmetry_15", ".pdb", extra_args=[fn_cif])
  assert_placeholder_pdb(txt)

def test_place_hydrogens_does_not_move_model():
  model = place_lib(pdb_str)
  assert n_h(model.get_hierarchy()) > 0
  ref = heavy_xyz(iotbx.pdb.input(lines=pdb_str.split("\n"),
    source_info=None).construct_hierarchy())
  out = heavy_xyz(model.get_hierarchy())
  assert set(ref) == set(out)
  for k in ref:
    assert approx_equal(ref[k], out[k], eps=1.e-6), k

def test_place_hydrogens_translation_invariant():
  # The unshifted P1 box has lattice planes through the model; placement must
  # not depend on where they fall.
  h = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None
    ).construct_hierarchy()
  xyz = h.atoms().extract_xyz()
  extent = [mx - mn for mx, mn in zip(xyz.max(), xyz.min())]
  mid = [(mx + mn) / 2 for mx, mn in zip(xyz.max(), xyz.min())]
  t = [3 * (e + 10) - m for e, m in zip(extent, mid)]
  h.atoms().set_xyz(xyz + tuple(t))
  a = all_xyz(place_lib(pdb_str).get_hierarchy())
  b = all_xyz(place_lib(h.as_pdb_string()).get_hierarchy())
  assert set(a) == set(b)
  assert len(a) > 21
  for k in a:
    back = [x - s for x, s in zip(b[k], t)]
    assert approx_equal(a[k], back, eps=1.e-3), k

if __name__ == "__main__":
  test_cli_pdb_no_cryst1()
  test_cli_cif_no_cell()
  test_cli_keeps_real_cryst1()
  test_cli_box_cushion()
  test_hydrogenate_pdb_no_cryst1()
  test_hydrogenate_cif_no_cell()
  test_hydrogenate_keeps_real_cryst1()
  test_cli_restraints_no_cryst1()
  test_hydrogenate_restraints_no_cryst1()
  test_cli_keeps_placeholder_pdb()
  test_cli_keeps_placeholder_cif()
  test_hydrogenate_keeps_placeholder()
  test_placeholder_with_restraints()
  test_place_hydrogens_does_not_move_model()
  test_place_hydrogens_translation_invariant()
  print("OK")
