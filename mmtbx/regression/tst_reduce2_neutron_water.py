"""reduce2 neutron tests: waters: layouts, waters without O or with altloc H,
DOD input"""
from __future__ import absolute_import, division, print_function
import iotbx.pdb
from libtbx.utils import Sorry
from libtbx.test_utils import approx_equal
from scitbx import matrix
from mmtbx.regression.model.tst_model_idealize_h_riding_water import layouts
from mmtbx.regression import reduce2_neutron_utils
from mmtbx.regression.reduce2_neutron_utils import (
  RENAMED, cached_run, check_water, pdb_1yjp, pdb_1yjp_dod, pdb_1yjp_dod_d,
  run_reduce2, water_geometry, water_groups)

reduce2_neutron_utils.set_prefix("tst_r2n_water")

def test_water_layouts():
  # 9l: waters the atom_group rule skips, from 0.90 A inputs; per O and
  # conformer: (number of H/D, expected number of groups)
  cryst1 = "CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1\n"
  for key, n_groups in (("a", 1), ("b", 1), ("c", 2), ("d", 2), ("e", 1)):
    pdb_str = cryst1 + layouts[key] + "END\n"
    prefix = "tst_r2n_water_layout_" + key
    with open(prefix + "_in.pdb", "w") as f:
      f.write(pdb_str)
    for lengths, angle in water_geometry(prefix + "_in.pdb"):
      for d in lengths: assert approx_equal(d, 0.90, eps=0.001, out=None)
    if key == "e":
      # DOD with D, kept by optimize: needs a neutron_option since step 4
      # (the hd_and_h run of this layout is test 9l(e))
      try:
        run_reduce2(prefix, pdb_str,
          ["approach=optimize", "use_neutron_distances=True"])
      except Sorry as e:
        assert "input contains deuterium" in str(e), str(e)
      else:
        raise AssertionError("D input with none ran")
      continue
    fn, log, description = run_reduce2(prefix, pdb_str,
      ["approach=optimize", "use_neutron_distances=True"])
    check_water(fn, 0.98, n_groups)
    if key == "c":
      # A and B keep distinct H positions
      h = iotbx.pdb.input(fn).construct_hierarchy()
      xyz = dict(((a.name.strip(), a.parent().altloc), matrix.col(a.xyz))
                 for a in h.atoms())
      for name in ("H1", "H2"):
        assert (xyz[(name, "A")] - xyz[(name, "B")]).length() > 0.5, name

def test_water_without_o_add():
  # approach=add keep_existing_H=True on 1yjp with water H and water 8's O
  # deleted: no crash in workaround_002(); the two H are removed and reported
  # as missing their heavy atom, not as not parameterized
  pdb_str = "\n".join(l for l in cached_run(pdb_1yjp,
    exclude_water=False).text.splitlines()
    if " O   HOH A   8" not in l) + "\n"
  fn, log, description = run_reduce2("tst_r2n_water_no_o", pdb_str,
    ["approach=add", "keep_existing_H=True", "exclude_water=False"])
  h = iotbx.pdb.input(fn).construct_hierarchy()
  assert not [a for a in h.atoms() if a.parent().resname == "HOH" and
              a.parent().parent().resseq.strip() == "8"]
  report = log.split("were not placed because a neighbouring heavy atom is")
  assert len(report) == 2 and "HOH A 8" in report[1].split("\n\n")[0], log
  for block in log.split("could not be parameterized")[1:]:
    assert "HOH A   8" not in block.split("\n\n")[0]

def test_water_single_altloc_h():
  # water 8 with a blank O and H1/H2 in altloc A, approach=add
  # keep_existing_H=True: one pair, still in A, at the neutron distance; no
  # blank pair added next to it, nothing listed as not parameterized
  o_line = ("HETATM   61  O   HOH A   8      -6.471   5.227   7.124  1.00 22.62"
            "           O\n")
  assert pdb_1yjp.count(o_line) == 1
  pdb_str = pdb_1yjp.replace(o_line, o_line +
    "HETATM   61  H1 AHOH A   8      -5.571   5.327   7.124  1.00 22.62"
    "           H\n"
    "HETATM   61  H2 AHOH A   8      -6.771   6.127   7.124  1.00 22.62"
    "           H\n")
  fn, log, description = run_reduce2("tst_r2n_water_w8a", pdb_str,
    ["approach=add", "keep_existing_H=True", "exclude_water=False",
     "use_neutron_distances=True"])
  h = iotbx.pdb.input(fn).construct_hierarchy()
  w8 = sorted((a.name.strip(), a.parent().altloc) for a in h.atoms()
              if a.parent().parent().resseq.strip() == "8")
  assert w8 == [("H1", "A"), ("H2", "A"), ("O", "")], w8
  check_water(fn, 0.98, 7)
  for block in log.split("could not be parameterized")[1:]:
    assert "HOH A   8" not in block.split("\n\n")[0], log

def test_dod_input_all_h():
  # input waters named DOD (O only), all_h: placement adds H, the written
  # waters are HOH with H and there is no DOD anywhere
  assert pdb_1yjp_dod.count("DOD A") == 7
  for approach in ("add", "place"):
    fn, log, description = run_reduce2("tst_r2n_water_dod_" + approach,
      pdb_1yjp_dod, ["neutron_option=all_h", "approach=" + approach,
      "exclude_water=False", "keep_existing_H=False"])
    h = iotbx.pdb.input(fn).construct_hierarchy()
    waters = [ag for ag in h.atom_groups() if ag.resname in ("HOH", "DOD")]
    assert len(waters) == 7
    assert set(ag.resname for ag in waters) == set(["HOH"])
    assert set(a.element.strip() for ag in waters for a in ag.atoms()) == \
      set(["O", "H"])
    assert "DOD" not in open(fn).read()
    check_water(fn, 0.98, 7)

def test_output_check_water_names():
  # the water-name check: renaming left out, a water with H named DOD
  from mmtbx.hydrogens import hd_exchange
  original = hd_exchange.apply_neutron_option
  def no_rename(hierarchy, *args, **kwargs):
    r = original(hierarchy, *args, **kwargs)
    for ag in hierarchy.atom_groups():
      if ag.resname == "HOH": ag.resname = "DOD"
    return r
  hd_exchange.apply_neutron_option = no_rename
  try:
    run_reduce2("tst_r2n_water_wname", pdb_1yjp_dod,
      ["neutron_option=all_h", "exclude_water=False"])
  except Sorry as e:
    assert "Internal error: H/D output check failed" in str(e), str(e)
    assert "water name DOD does not match its H/D" in str(e), str(e)
  else:
    raise AssertionError("water-name check did not fire")
  finally:
    hd_exchange.apply_neutron_option = original

def test_dod_input_none():
  # none, add and place: placement replaced the D, so the waters are HOH
  for approach in ("add", "place"):
    fn, log, description = run_reduce2("tst_r2n_water_dodn_" + approach,
      pdb_1yjp_dod, ["approach=" + approach, "exclude_water=False"])
    assert water_groups(fn) == [("HOH", ["H", "O"])] * 7, water_groups(fn)
    with open(fn) as f:
      assert "DOD" not in f.read()
    assert "7 " + RENAMED in log and "7 " + RENAMED in description

def test_dod_input_none_exclude_water():
  # no H added to water: the waters stay DOD with O only, nothing logged
  fn, log, description = run_reduce2("tst_r2n_water_dodx", pdb_1yjp_dod,
    [])
  assert water_groups(fn) == [("DOD", ["O"])] * 7, water_groups(fn)
  assert RENAMED not in log and RENAMED not in description

def test_dod_input_with_d_kept():
  # keep_existing_H=True and a DOD water with D1/D2, neutron_option=none:
  # the D would be kept, so the run stops (before step 4 placement kept
  # D1/D2 and added H1/H2)
  assert pdb_1yjp_dod_d.count("D1  DOD") == 1
  try:
    run_reduce2("tst_r2n_water_dodd", pdb_1yjp_dod_d,
      ["keep_existing_H=True", "exclude_water=False"])
  except Sorry as e:
    assert "input contains deuterium; set neutron_option" in str(e), str(e)
    assert str(e).endswith("rebuild H: keep_existing_H=False"), str(e)
  else:
    raise AssertionError("D input with none ran")

def test_dod_input_with_d_normalized():
  # the same water with hd_and_h: D1/D2 normalized to H1/H2 before
  # placement, no extra H, written as DOD with D1/D2 at 0.98 A
  fn, log, description = run_reduce2("tst_r2n_water_doddn",
    pdb_1yjp_dod_d, ["keep_existing_H=True", "exclude_water=False",
                     "neutron_option=hd_and_h"])
  assert water_groups(fn) == [("DOD", ["D", "O"])] * 7, water_groups(fn)
  h = iotbx.pdb.input(fn).construct_hierarchy()
  w8 = sorted(a.name.strip() for a in h.atoms()
              if a.parent().parent().resseq.strip() == "8")
  assert w8 == ["D1", "D2", "O"], w8
  check_water(fn, 0.98, 7)

def test_dod_input_hd_and_h():
  # the earlier rename to HOH does not break the neutron path or the
  # resname invariant: written as DOD with D, output check passes
  fn, log, description = run_reduce2("tst_r2n_water_dodhd",
    pdb_1yjp_dod, ["neutron_option=hd_and_h", "exclude_water=False"])
  assert water_groups(fn) == [("DOD", ["D", "O"])] * 7, water_groups(fn)
  check_water(fn, 0.98, 7)

def test_water_layout_e_neutron():
  # 9l(e): a DOD water with D at 0.90 A, optimize hd_and_h: normalized and
  # written as DOD with D at 0.98 A / 103.91 deg
  pdb_str = ("CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1\n" +
    layouts["e"] + "END\n")
  with open("tst_r2n_water_l9e_in.pdb", "w") as f:
    f.write(pdb_str)
  for lengths, angle in water_geometry("tst_r2n_water_l9e_in.pdb"):
    for d in lengths: assert approx_equal(d, 0.90, eps=0.001, out=None)
  fn, log, description = run_reduce2("tst_r2n_water_l9e", pdb_str,
    ["approach=optimize", "neutron_option=hd_and_h"])
  assert water_groups(fn) == [("DOD", ["D", "O"])], water_groups(fn)
  check_water(fn, 0.98, 1)

if __name__ == "__main__":
  test_water_layouts()
  test_water_without_o_add()
  test_water_single_altloc_h()
  test_dod_input_all_h()
  test_output_check_water_names()
  test_dod_input_none()
  test_dod_input_none_exclude_water()
  test_dod_input_with_d_kept()
  test_dod_input_with_d_normalized()
  test_dod_input_hd_and_h()
  test_water_layout_e_neutron()
  print("OK")
