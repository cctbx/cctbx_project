"""reduce2 water_hydrogens: the settings of the water H placement
(mmtbx.hydrogens.water_protonation) and when they can take effect.
"""
import os
import iotbx.pdb
import iotbx.phil
from libtbx.utils import null_out, Sorry
from libtbx.test_utils import approx_equal, Exception_expected
from iotbx.cli_parser import CCTBXParser, run_program
from iotbx.data_manager import DataManager
from mmtbx.programs import reduce2

# Water 4 has two H turned away from its acceptor (Ala 1 O), water 5 one H.
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
HETATM   23  H1  HOH A   4      12.729  -2.058  12.341  1.00 20.00           H
HETATM   24  H2  HOH A   4      13.553  -2.058  11.286  1.00 20.00           H
HETATM   25  O   HOH A   5      20.007  -3.513   9.351  1.00 20.00           O
HETATM   26  H1  HOH A   5      20.007  -3.513  10.201  1.00 20.00           H
END
"""

# A water on the two-fold of P 1 2 1 (along b through x = z = 0), bare and
# with one H at 1.0, which with its equivalent makes the water's two.
two_fold_str = """\
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1 2 1
HETATM    1  O   HOH W   1       0.000   5.000   0.000  1.00 10.00           O
END
"""
two_fold_h_str = two_fold_str.replace("END", (
  "HETATM    2  H1  HOH W   1       0.757   5.586   0.000  1.00 10.00"
  "           H\nEND"))

def run_reduce2(prefix, extra_args=(), model_str=pdb_str):
  """Run reduce2 on model_str; returns the output hierarchy and description."""
  fn_in, fn_out = prefix + ".pdb", prefix + "_out.pdb"
  with open(fn_in, "w") as f:
    f.write(model_str)
  if os.path.exists(fn_out): os.remove(fn_out)
  run_program(program_class=reduce2.Program, logger=null_out(),
    args=[fn_in, "output.filename=" + fn_out, "overwrite=True",
          "output.description_file=" + prefix + ".txt"] + list(extra_args))
  with open(prefix + ".txt") as f:
    description = f.read()
  return iotbx.pdb.input(fn_out).construct_hierarchy(), description

def validate(args):
  """Program.validate() of reduce2 on pdb_str with args, without a run."""
  dm = DataManager(reduce2.Program.datatypes)
  dm.process_model_str("tst_reduce2_water_h_settings.pdb", pdb_str)
  master = CCTBXParser(program_class=reduce2.Program,
    logger=null_out()).master_phil
  working, unused = master.fetch(source=iotbx.phil.parse("\n".join(args)),
    track_unused_definitions=True)
  assert not unused, [u.path for u in unused]
  reduce2.Program(dm, working.extract(), logger=null_out()).validate()

def water_h(hierarchy):
  """Water H coordinates by (resseq, atom name)."""
  sel = hierarchy.atom_selection_cache().selection("water and element H")
  return {(a.parent().parent().resseq.strip(), a.name.strip()): a.xyz
          for a in hierarchy.select(sel).atoms()}

def placer_arguments(prefix, args):
  """Run reduce2; returns the keyword arguments of each placer call."""
  from mmtbx.hydrogens import water_protonation
  calls = []
  original = water_protonation.place_water_hydrogens
  def recording(hier, **kwargs):
    calls.append(kwargs)
    return original(hier, **kwargs)
  water_protonation.place_water_hydrogens = recording
  try:
    run_reduce2(prefix, args)
  finally:
    water_protonation.place_water_hydrogens = original
  return calls

def exercise_defaults():
  """
  reduce2's water_hydrogens defaults, from its master PHIL: the search
  settings of mmtbx.development.virvel's program, existing_h=reorient and
  symmetry=True.
  """
  from mmtbx.programs import water_protonation as virvel
  r = iotbx.phil.parse(reduce2.Program.master_phil_str,
    process_includes=True).extract().water_hydrogens
  v = iotbx.phil.parse(virvel.Program.master_phil_str).extract()
  for name, got, want in (
      ("refine.max_sweeps", r.refine.max_sweeps, v.refine.max_sweeps),
      ("refine.tolerance", r.refine.tolerance, v.refine.tolerance),
      ("basin.rounds", r.basin.rounds, v.basin.rounds),
      ("cluster_search", r.cluster_search, v.cluster_search)):
    assert got == want, (name, got, want)
  assert r.existing_h == "reorient", r.existing_h
  assert r.symmetry is True, r.symmetry

def exercise_settings_reach_placer():
  """
  existing_h, refine.max_sweeps, refine.tolerance, basin.rounds and
  cluster_search reach place_water_hydrogens as given, and symmetry=False
  as symmetry_mates=False with the crystal symmetry kept.
  """
  calls = placer_arguments("tst_reduce2_water_h_settings_placer", [
    "add_h_to_water=True", "keep_existing_H=True",
    "water_hydrogens.existing_h=keep", "water_hydrogens.symmetry=False",
    "water_hydrogens.refine.max_sweeps=0",
    "water_hydrogens.refine.tolerance=0", "water_hydrogens.basin.rounds=2",
    "water_hydrogens.cluster_search=True"])
  assert len(calls) == 1, len(calls)
  expected = dict(existing_h="keep", symmetry_mates=False, n_refine=0,
    refine_tol=0, n_basin=2, cluster_search=True)
  for key, value in expected.items():
    assert calls[0][key] == value, (key, calls[0][key], value)
  assert calls[0]["crystal_symmetry"] is not None, "crystal symmetry dropped"

def exercise_existing_h():
  """
  Water H in the input, kept by keep_existing_H=True: existing_h=keep
  leaves both waters as they are, places no H and reports water 5 as a
  single-H water left as is; existing_h=complete adds water 5's second H,
  leaves the input H where they are, places one H and reports water 5 as
  completed. The default, existing_h=reorient, places all four water H again
  and reports water 5 as placed again.
  """
  inp = water_h(iotbx.pdb.input(source_info=None,
    lines=pdb_str.split("\n")).construct_hierarchy())
  for mode, n_placed, report in (("keep", 0, "left as is"),
                                 ("complete", 1, "completed")):
    h, description = run_reduce2("tst_reduce2_water_h_settings_" + mode,
      ["add_h_to_water=True", "keep_existing_H=True",
       "water_hydrogens.existing_h=" + mode])
    out = water_h(h)
    assert len(out) == len(inp) + n_placed, (mode, sorted(out))
    for key, site in inp.items():
      assert approx_equal(out[key], site, eps=0.01), (mode, key)
    for line in (f"Water H placed: {n_placed} ",
                 f"1 single-H water {report}:"):
      assert line in description, (mode, line, description)
  _, description = run_reduce2("tst_reduce2_water_h_settings_reorient",
    ["add_h_to_water=True", "keep_existing_H=True"])
  for line in ("Water H placed: 4 ",
               "Water H in the input: 3 (existing_h=reorient)",
               "1 single-H water placed again:"):
    assert line in description, ("reorient", line, description)

def exercise_report_names():
  """
  The report names the waters and their H as the written model has them.
  Two models with neutron_option=all_d and existing_h=keep: the placed
  water H are D, the input ones H/D, and water 5 is listed once per model,
  as DOD with "(model N)" after it, as in reduce2's other reports.
  approach=optimize with the waters named DOD: water 5 carries H, so it is
  listed as HOH, its name in the written model. Waters 5 and 5A, each with
  one H, are listed apart.
  """
  inp = iotbx.pdb.input(source_info=None, lines=pdb_str.split("\n"))
  h = inp.construct_hierarchy()
  copy = h.only_model().detached_copy()
  h.only_model().id = "1"
  copy.id = "2"
  h.append_model(copy)
  out, description = run_reduce2("tst_reduce2_water_h_settings_models",
    ["add_h_to_water=True", "keep_existing_H=True", "neutron_option=all_d",
     "water_hydrogens.existing_h=keep"],
    model_str=h.as_pdb_string(crystal_symmetry=inp.crystal_symmetry()))
  assert out.models_size() == 2, out.models_size()
  names = set(ag.resname for ag in out.atom_groups()
              if ag.parent().resseq.strip() == "5")
  assert names == {"DOD"}, names
  for line in ("Water D placed: 0 ", "Water H/D in the input: 6 ",
               "2 single-D waters left as is:", "  DOD A 5 (model 1)\n",
               "  DOD A 5 (model 2)\n"):
    assert line in description, (line, description)
  out, description = run_reduce2("tst_reduce2_water_h_settings_dod",
    ["approach=optimize", "add_h_to_water=True",
     "water_hydrogens.existing_h=keep"],
    model_str=pdb_str.replace("HOH", "DOD"))
  names = set(ag.resname for ag in out.atom_groups()
              if ag.parent().resseq.strip() == "5")
  assert names == {"HOH"}, names
  assert "  HOH A 5\n" in description, description
  # waters 5 and 5A, each with one H, are listed apart
  water_5a = (
    "HETATM   27  O   HOH A   5A     24.007  -3.513   9.351  1.00 20.00"
    "           O\n"
    "HETATM   28  H1  HOH A   5A     24.007  -3.513  10.201  1.00 20.00"
    "           H\nEND")
  _, description = run_reduce2("tst_reduce2_water_h_settings_icode",
    ["add_h_to_water=True", "keep_existing_H=True",
     "water_hydrogens.existing_h=keep"],
    model_str=pdb_str.replace("END", water_5a))
  for line in ("2 single-H waters left as is:", "  HOH A 5\n",
               "  HOH A 5A\n"):
    assert line in description, (line, description)

def exercise_symmetry_off_site_symmetry():
  """
  water_hydrogens.symmetry=False leaves out the symmetry mates, not the
  site symmetry: a bare water on a two-fold, placed again or completed,
  gets two H at half occupancy, and a water whose H at 1.0 already makes
  its two gets nothing added when completed.
  """
  for name, model_str, existing_h, expected in (
      ("bare", two_fold_str, "reorient", [0.5, 0.5]),
      ("bare", two_fold_str, "complete", [0.5, 0.5]),
      ("h", two_fold_h_str, "reorient", [0.5, 0.5]),
      ("h", two_fold_h_str, "complete", [1.0])):
    out, _ = run_reduce2(
      f"tst_reduce2_water_h_settings_two_fold_{name}_{existing_h}",
      ["add_h_to_water=True", "keep_existing_H=True",
       "water_hydrogens.symmetry=False",
       f"water_hydrogens.existing_h={existing_h}"], model_str=model_str)
    sel = out.atom_selection_cache().selection("water and element H")
    occ = sorted(a.occ for a in out.select(sel).atoms())
    assert approx_equal(occ, expected), (name, existing_h, occ)

def exercise_settings_accepted():
  """
  validate() accepts the settings where they can take effect: qmready=True
  (approach=add) with a water_hydrogens setting and with existing_h=keep,
  and approach=optimize with add_h_to_water=True and existing_h=complete.
  """
  for args in (["qmready=True", "water_hydrogens.basin.rounds=2"],
               ["qmready=True", "water_hydrogens.existing_h=keep"],
               ["approach=optimize", "add_h_to_water=True",
                "water_hydrogens.existing_h=complete"]):
    validate(args)

def exercise_settings_without_water_h():
  """
  A water_hydrogens setting when reduce2 places no water H is a Sorry: with
  approach=add and approach=optimize without add_h_to_water=True, and with
  approach=remove. So is existing_h=keep or complete when approach=add or
  place drops the input H (keep_existing_H=False); the message names the
  settings that keep them for that approach.
  """
  for args, text in (
      (["water_hydrogens.basin.rounds=2"],
       "water_hydrogens requires add_h_to_water=True"),
      (["approach=optimize", "water_hydrogens.existing_h=keep"],
       "water_hydrogens requires add_h_to_water=True"),
      (["approach=remove", "water_hydrogens.basin.rounds=2"],
       "water_hydrogens cannot be used with approach=remove"),
      (["add_h_to_water=True", "water_hydrogens.existing_h=complete"],
       "existing_h=complete needs the input H kept: keep_existing_H=True or "
       "qmready=True."),
      (["approach=place", "add_h_to_water=True",
        "water_hydrogens.existing_h=keep"],
       "existing_h=keep needs the input H kept: keep_existing_H=True.")):
    try:
      run_reduce2("tst_reduce2_water_h_settings_without", args)
    except Sorry as e:
      assert text in str(e), (args, str(e))
    else:
      raise Exception_expected

def exercise_settings_need_values():
  """A water_hydrogens setting given as None is a Sorry naming it."""
  for path in ("existing_h", "symmetry", "refine.max_sweeps",
               "refine.tolerance", "basin.rounds", "cluster_search"):
    try:
      run_reduce2("tst_reduce2_water_h_settings_none",
        ["add_h_to_water=True", f"water_hydrogens.{path}=None"])
    except Sorry as e:
      assert f"water_hydrogens.{path} needs a value" in str(e), (path,
        str(e))
    else:
      raise Exception_expected

if __name__ == "__main__":
  exercise_defaults()
  exercise_settings_reach_placer()
  exercise_existing_h()
  exercise_report_names()
  exercise_symmetry_off_site_symmetry()
  exercise_settings_accepted()
  exercise_settings_without_water_h()
  exercise_settings_need_values()
  print("OK")
