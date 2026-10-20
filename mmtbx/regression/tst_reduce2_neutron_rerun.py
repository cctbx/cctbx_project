"""reduce2 neutron tests: each mode's output rerun in the same and in every
other mode"""
from __future__ import absolute_import, division, print_function
from mmtbx.regression import reduce2_neutron_utils
from mmtbx.regression.reduce2_neutron_utils import (
  NEUTRON, all_atoms, check_water, dipeptide_rotatable, dipeptide_water_out,
  expected_hd, h_name, hg_pair_03_07, run_reduce2, water_groups)

reduce2_neutron_utils.set_prefix("tst_r2n_rerun")

def test_rerun_same_mode():
  # 2: each mode's written output (dipeptide + water), run again in the same
  # mode with optimize and with add keep_existing_H=True -> the same labels,
  # occupancy and B per atom, and the same positions except for the H the
  # Optimizer rotates. all_d goes through the rename path without pairs
  restored = {"all_h": 0, "all_d": 0, "hd_and_h": 3, "hd_and_d": 3,
              "all_hd": 7}
  for option in NEUTRON:
    src = dipeptide_water_out(option)
    first = all_atoms(src, is_file=False)
    for args in (["approach=optimize"],
                 ["approach=add", "keep_existing_H=True"]):
      fn, log, description = run_reduce2("tst_r2n_rerun_rr_%s_%s" % (
        option, args[0][9:]), src,
        args + ["neutron_option=" + option, "exclude_water=False"])
      second = all_atoms(fn)
      assert sorted(second) == sorted(first), (option, args,
        sorted(set(second) ^ set(first)))
      for key, (occ, b, xyz) in first.items():
        assert second[key][:2] == (occ, b), (option, args, key)
        if (key[0], h_name(key[1], key[3])) in dipeptide_rotatable: continue
        assert abs(second[key][2] - xyz) < 0.01, (option, args, key)
      assert "Values restored from the input: %d" % restored[option] in \
        description, (option, args)

def test_cross_mode():
  # 2b: each mode's output (dipeptide + water) as input for every other mode,
  # through optimize; in the paired inputs the Ser HG pair is at 0.3 / 0.7
  # with different B. Labels per "What each mode writes", values per
  # "Restoration policy", water O-D/H at 0.98 A
  inputs = {}
  for option in NEUTRON:
    src = dipeptide_water_out(option)
    if option in ("hd_and_h", "hd_and_d", "all_hd"):
      src = hg_pair_03_07(src)
      assert src.count("  0.30 15.00") == 1 and src.count("  0.70 25.00") == 1
    inputs[option] = src
  failed = []
  for source in NEUTRON:
    for option in NEUTRON:
      if option == source: continue
      fn, log, description = run_reduce2("tst_r2n_rerun_x_%s_%s" % (
        source, option), inputs[source], ["approach=optimize",
        "neutron_option=" + option, "exclude_water=False"])
      got = dict((k, v[:2]) for k, v in all_atoms(fn).items()
                 if k[3] in ("H", "D"))
      expected = expected_hd(inputs[source], option)
      if got != expected:
        failed.append("%s -> %s: %s" % (source, option, sorted(
          (k, expected.get(k), got.get(k)) for k in set(got) | set(expected)
          if got.get(k) != expected.get(k))))
      check_water(fn, 0.98, 1)
      assert water_groups(fn) == [("HOH" if option == "all_h" else "DOD",
        ["H" if option == "all_h" else "D", "O"])], (source, option)
  assert not failed, "\n".join(failed)

if __name__ == "__main__":
  test_rerun_same_mode()
  test_cross_mode()
  print("OK")
