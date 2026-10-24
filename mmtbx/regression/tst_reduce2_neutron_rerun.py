"""reduce2 neutron tests: each mode's output rerun in the same and in every
other mode"""
from __future__ import absolute_import, division, print_function
from mmtbx.regression import reduce2_neutron_utils
from mmtbx.regression.reduce2_neutron_utils import (
  NEUTRON, all_atoms, check_water, dipeptide, dipeptide_out,
  dipeptide_rotatable, dipeptide_water_out, edit_atoms, expected_hd, h_name,
  hg_pair_03_07, run_reduce2, water_groups)

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

def hd_by_site(src, is_file=True, names=("HB2", "HB3")):
  """[(xyz, {element: (occupancy, B)})] of the Ser 37 H/D names, by site"""
  sites = []
  for key, (occ, b, xyz) in all_atoms(src, is_file).items():
    if key[0] != "37" or h_name(key[1], key[3]) not in names:
      continue
    for site_xyz, values in sites:
      if abs(site_xyz - xyz) < 0.2:
        values[key[3]] = (occ, b)
        break
    else:
      sites.append((xyz, {key[3]: (occ, b)}))
  return sites

def test_reversed_ch2_names():
  # 3: Ser 37 CB pairs with HB2/DB2 and HB3/DB3 names reversed (placement
  # swaps them back) and unequal values per site, through place
  # keep_existing_H=True: each physical site keeps its input values. all_hd
  # restores both, hd_and_h / hd_and_d write the parent occupancy with that
  # isotope's B (CB H are not exchangeable)
  src = dipeptide_out("all_hd")
  swap = {"HB2": "HB3", "HB3": "HB2", "DB2": "DB3", "DB3": "DB2"}
  for first, second in (((0.2, 15.), (0.7, 18.)), ((0.7, 18.), (0.2, 15.))):
    # (H occupancy, H B) by input name; D: 1 - occupancy, B + 10
    values = {"HB2": first, "HB3": second}
    def edit(a):
      if a.parent().parent().resseq.strip() != "37": return
      name = a.name.strip()
      if name not in swap: return
      occ, b = values[h_name(name, a.element.strip())]
      if a.element.strip() == "D": occ, b = round(1 - occ, 2), b + 10
      a.set_occ(occ); a.set_b(b)
      a.set_name(a.name.replace(name, swap[name]))
    src_reversed = edit_atoms(src, edit)
    sites_in = hd_by_site(src_reversed, is_file=False)
    assert len(sites_in) == 2 and all(len(v) == 2 for _, v in sites_in)
    for option in ("all_hd", "hd_and_h", "hd_and_d"):
      fn, log, description = run_reduce2("tst_r2n_rerun_ch2_%s_%s" % (
        option, first[0]), src_reversed, ["approach=place",
        "keep_existing_H=True", "neutron_option=" + option])
      sites_out = hd_by_site(fn)
      assert len(sites_out) == 2, (option, sites_out)
      for xyz, v_in in sites_in:
        v_out = [v for x, v in sites_out if abs(x - xyz) < 0.2]
        assert len(v_out) == 1, (option, xyz)
        if option == "all_hd":
          expected = v_in
        else:
          element = "H" if option == "hd_and_h" else "D"
          expected = {element: (1.0, v_in[element][1])}
        assert v_out[0] == expected, (option, first, v_in, v_out[0])

_nh3_out = []
def dipeptide_nh3_out():
  """all_hd output of the dipeptide with an N-terminal NH3 on Ser 37"""
  if not _nh3_out:
    fn, _, _ = run_reduce2("tst_r2n_rerun_nh3_src", dipeptide,
      ["neutron_option=all_hd", "n_terminal_charge=first_in_chain"])
    with open(fn) as f: _nh3_out.append(f.read())
  return _nh3_out[0]

def test_reversed_propeller_names():
  # 3: Ser 37 NH3 pairs with H2/D2 and H3/D3 names swapped and unequal
  # values per site, through place keep_existing_H=True: each physical site
  # keeps its input values (riding idealization used to turn the kept
  # propeller by one site). NH3 H are exchangeable: pairs in all three modes
  names = ("H1", "H2", "H3")
  swap = {"H2": "H3", "H3": "H2", "D2": "D3", "D3": "D2"}
  values = {"H1": (0.5, 20.), "H2": (0.2, 15.), "H3": (0.7, 18.)}
  def edit(a):
    name = a.name.strip()
    if a.parent().parent().resseq.strip() != "37" or \
       h_name(name, a.element.strip()) not in names: return
    occ, b = values[h_name(name, a.element.strip())]
    if a.element.strip() == "D": occ, b = round(1 - occ, 2), b + 10
    a.set_occ(occ); a.set_b(b)
    if name in swap: a.set_name(a.name.replace(name, swap[name]))
  src = edit_atoms(dipeptide_nh3_out(), edit)
  sites_in = hd_by_site(src, is_file=False, names=names)
  assert len(sites_in) == 3 and all(len(v) == 2 for _, v in sites_in)
  for option in ("all_hd", "hd_and_h", "hd_and_d"):
    fn, log, description = run_reduce2("tst_r2n_rerun_nh3_%s" % option, src,
      ["approach=place", "keep_existing_H=True", "neutron_option=" + option,
       "n_terminal_charge=first_in_chain"])
    sites_out = hd_by_site(fn, names=names)
    assert len(sites_out) == 3, (option, sites_out)
    for xyz, v_in in sites_in:
      v_out = [v for x, v in sites_out if abs(x - xyz) < 0.2]
      assert v_out == [v_in], (option, v_in, v_out)

def test_own_output_fed_back():
  # 3: reduce2's own all_hd output (NH3 on Ser 37) through place
  # keep_existing_H=True: every H/D keeps its label, values and site; the
  # rotatable NH3 and HG return to the ideal dihedral (the Optimizer turned
  # them), each H on its own site
  src = dipeptide_nh3_out()
  fn, log, description = run_reduce2("tst_r2n_rerun_own_nh3", src,
    ["approach=place", "keep_existing_H=True", "neutron_option=all_hd",
     "n_terminal_charge=first_in_chain"])
  first, second = all_atoms(src, is_file=False), all_atoms(fn)
  assert sorted(first) == sorted(second), sorted(set(first) ^ set(second))
  for key, (occ, b, xyz) in second.items():
    # H and D of a pair share the site: nearest of the same element
    nearest = min([k for k in first if k[3] == key[3]],
      key=lambda k: abs(first[k][2] - xyz))
    assert nearest == key, (key, nearest)
    assert first[key][:2] == (occ, b), key
    limit = 0.35 if h_name(key[1], key[3]) in ("H1", "H2", "H3", "HG") \
      else 0.01
    assert abs(first[key][2] - xyz) < limit, (key, abs(first[key][2] - xyz))

if __name__ == "__main__":
  test_rerun_same_mode()
  test_cross_mode()
  test_reversed_ch2_names()
  test_reversed_propeller_names()
  test_own_output_fed_back()
  print("OK")
