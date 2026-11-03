"""reduce2 neutron tests: neutron and X-ray X-H and water distances, per
approach and mode"""
from __future__ import absolute_import, division, print_function
import iotbx.pdb
from libtbx.test_utils import approx_equal
from mmtbx.reduce import Optimizers
from mmtbx.regression.model.tst_model_idealize_h_riding_water import \
  pdb_water, hoh_cif_096, make_model, check_oh
from mmtbx.regression import reduce2_neutron_utils
from mmtbx.regression.reduce2_neutron_utils import (
  INFO_FALSE, INFO_TRUE, NEUTRON, cached_run, check_bond_lengths, check_no_d,
  check_water, hd_pairs, pdb_1yjp, pdb_truncated_lys, run_reduce2)

reduce2_neutron_utils.set_prefix("tst_r2n_distances")

def test_xray_default():
  # control: the default run stays at X-ray distances
  prog = cached_run(pdb_1yjp, exclude_water=False)
  fn, description = prog.params.output.filename, prog.description
  check_bond_lengths(fn, use_neutron_distances=False)
  check_water(fn, 0.85, 7)
  assert INFO_FALSE in description and INFO_TRUE not in description

def test_neutron_add():
  # 9b: also explicit neutron_option=none (none keeps today's distances)
  fn, log, description = run_reduce2("tst_r2n_distances_add", pdb_1yjp,
    ["use_neutron_distances=True", "exclude_water=False",
     "neutron_option=none"])
  check_bond_lengths(fn, use_neutron_distances=True)
  check_water(fn, 0.98, 7)
  assert INFO_TRUE in description and INFO_FALSE not in description
  check_no_d(fn)

def test_neutron_place():
  fn, log, description = run_reduce2("tst_r2n_distances_place", pdb_1yjp,
    ["use_neutron_distances=True", "exclude_water=False", "approach=place"])
  check_bond_lengths(fn, use_neutron_distances=True)
  check_water(fn, 0.98, 7)

def test_neutron_optimize():
  # input water H at 0.85 A (X-ray run), optimized with neutron distances
  prog = cached_run(pdb_1yjp, exclude_water=False)
  check_water(prog.params.output.filename, 0.85, 7)
  pdb_str = prog.text
  fn, log, description = run_reduce2("tst_r2n_distances_opt", pdb_str,
    ["use_neutron_distances=True", "exclude_water=False",
     "approach=optimize"])
  check_water(fn, 0.98, 7)
  assert INFO_TRUE in description and INFO_FALSE not in description

def test_qmready():
  # 9b: also explicit neutron_option=none
  fn, log, description = run_reduce2("tst_r2n_distances_qmready", pdb_1yjp,
    ["qmready=True", "neutron_option=none"])
  assert "qmready=True: use_neutron_distances set to True" in log
  check_bond_lengths(fn, use_neutron_distances=True)
  check_water(fn, 0.98, 7)
  assert INFO_TRUE in description and INFO_FALSE not in description
  check_no_d(fn)

def test_neutron_after_select():
  # placement drops the unparameterizable Lys HD2/HD3 with model.select();
  # the water must still get the neutron target afterwards
  fn, log, description = run_reduce2("tst_r2n_distances_select",
    pdb_truncated_lys, ["use_neutron_distances=True", "exclude_water=False"])
  h = iotbx.pdb.input(fn).construct_hierarchy()
  names = [a.name.strip() for a in h.atoms() if a.parent().resname == "LYS"]
  assert "HD2" not in names and "HD3" not in names, names
  check_water(fn, 0.98, 1)

def test_optimizer_passes_flag():
  # the Optimizer's own setup_riding_h_manager() call idealizes the water
  for flag, d in ((True, 0.96), (False, 0.85)):
    model = make_model(pdb_water, True, "n_gaussian", cif=hoh_cif_096)
    Optimizers.Optimizer(Optimizers._philLike(), False, model,
      useNeutronDistances=flag, verbosity=0)
    check_oh(model, d)
  # default useNeutronDistances: scattering-table rule as before
  model = make_model(pdb_water, False, "electron")
  Optimizers.Optimizer(Optimizers._philLike(), False, model, verbosity=0)
  check_oh(model, 0.98)

def test_distances_all_modes():
  # neutron targets for every mode (add), and with place for a split mode
  # and a one-isotope mode (the output step is the same for add and place);
  # water 0.98 A; H and D of a pair at the same position. all_h switches
  # neutron distances on and writes H only
  for option in NEUTRON:
    for approach in (("add", "place") if option in ("hd_and_h", "all_d")
                     else ("add",)):
      prefix = "tst_r2n_distances_32_%s_%s" % (option, approach)
      fn, log, description = run_reduce2(prefix, pdb_1yjp,
        ["neutron_option=" + option, "approach=" + approach,
         "exclude_water=False"])
      check_bond_lengths(fn, use_neutron_distances=True)
      check_water(fn, 0.98, 7)
      h = iotbx.pdb.input(fn).construct_hierarchy()
      waters = set(ag.resname for ag in h.atom_groups()
                   if ag.resname in ("HOH", "DOD"))
      assert waters == set(["HOH" if option == "all_h" else "DOD"]), waters
      pairs = hd_pairs(h)
      assert (len(pairs) > 0) == (option in ("hd_and_h", "hd_and_d",
                                             "all_hd")), (option, len(pairs))
      for hh, dd in pairs:
        assert approx_equal(hh.xyz, dd.xyz, eps=1e-3), (hh.id_str())
      if option == "all_h":
        assert "neutron_option=all_h: use_neutron_distances set to True" in log
        assert "D" not in [a.element.strip() for a in h.atoms()]

def test_water_after_select_hd():
  # 9k(a): placement drops the unparameterizable Lys HD2/HD3 (select());
  # the water O-D stays 0.98 A
  fn, log, description = run_reduce2("tst_r2n_distances_9ka",
    pdb_truncated_lys, ["neutron_option=hd_and_h", "exclude_water=False"])
  h = iotbx.pdb.input(fn).construct_hierarchy()
  names = [a.name.strip() for a in h.atoms() if a.parent().resname == "LYS"]
  assert "HD2" not in names and "HD3" not in names, names
  assert [ag.resname for ag in h.atom_groups()
          if ag.resname in ("HOH", "DOD")] == ["DOD"]
  check_water(fn, 0.98, 1)

if __name__ == "__main__":
  test_xray_default()
  test_neutron_add()
  test_neutron_place()
  test_neutron_optimize()
  test_qmready()
  test_neutron_after_select()
  test_optimizer_passes_flag()
  test_distances_all_modes()
  test_water_after_select_hd()
  print("OK")
