"""reduce2 neutron tests: H/D without a covalent partner (floating) and the
Optimizer and probe with D"""
from __future__ import absolute_import, division, print_function
import os
import iotbx.pdb
from libtbx.utils import null_out, Sorry
from scitbx import matrix
from mmtbx.programs import reduce2
from mmtbx.reduce import Optimizers
from mmtbx.regression import reduce2_neutron_utils
from mmtbx.regression.reduce2_neutron_utils import (
  FLOATING, atom_labels, cached_run, check_water, h_input_without_tyr_oh,
  h_to_d, make_program, pdb_1yjp, pdb_1yjp_h, pdb_ser_dg, placed_text,
  probe_setup, processed_model, run_reduce2, water_groups)

reduce2_neutron_utils.set_prefix("tst_r2n_optimizer")

def test_floating_h_removed():
  # 5a (none): approach=optimize on input with an H whose parent heavy atom
  # is missing: the run completes, the H is gone and listed, all other atoms
  # unchanged in count and names
  pdb_str = h_input_without_tyr_oh()
  with open("tst_r2n_optimizer_float_in.pdb", "w") as f:
    f.write(pdb_str)
  fn, log, description = run_reduce2("tst_r2n_optimizer_float", pdb_str,
    ["approach=optimize"])
  assert FLOATING in description and FLOATING in log
  listed = description.split(FLOATING)[1].strip().splitlines()[0]
  assert listed.split() == ["HH", "TYR", "A", "7"], listed
  before = [l for l in atom_labels("tst_r2n_optimizer_float_in.pdb")
            if l[1:3] != ("TYR", " HH ")]
  assert atom_labels(fn) == before

def test_optimize_without_h():
  # approach=optimize adds no H: input without H, or whose only H is floating
  # (removed before optimization), stops with a Sorry. Before: AttributeError
  # in the Optimizer (no riding H manager)
  text = "The model has no H or D atoms to optimize"
  only_hh = [l for l in h_input_without_tyr_oh().splitlines()
             if l[76:78].strip() != "H" or l[12:16] == " HH "]
  assert sum(1 for l in only_hh if l[76:78].strip() == "H") == 1
  for prefix, pdb_str in (("tst_r2n_optimizer_noh", pdb_1yjp),
      ("tst_r2n_optimizer_noh_float", "\n".join(only_hh) + "\n")):
    try:
      run_reduce2(prefix, pdb_str, ["approach=optimize"])
    except Sorry as e:
      assert text in str(e), str(e)
    else:
      raise AssertionError("no Sorry for %s" % prefix)

def test_floating_h_water_098():
  # 9k(b): water at 0.85 A in the input, a floating H removed with select();
  # with neutron distances the written water O-H is 0.98 A
  pdb_str = h_input_without_tyr_oh(exclude_water=False)
  fn, log, description = run_reduce2("tst_r2n_optimizer_float_w", pdb_str,
    ["approach=optimize", "use_neutron_distances=True",
     "exclude_water=False"])
  assert FLOATING in description
  assert not [l for l in atom_labels(fn) if l[1:3] == ("TYR", " HH ")]
  check_water(fn, 0.98, 7)

def test_no_floating_h_unchanged():
  # input without floating H/D: no removal, output as without the step
  pdb_str = pdb_1yjp_h()
  fn, log, description = run_reduce2("tst_r2n_optimizer_nofloat", pdb_str,
    ["approach=optimize"])
  assert FLOATING not in description
  original = reduce2.Program._RemoveFloatingHD
  reduce2.Program._RemoveFloatingHD = lambda self: ("", [])
  try:
    fn_ref, _, _ = run_reduce2("tst_r2n_optimizer_nofloat_ref", pdb_str,
      ["approach=optimize"])
  finally:
    reduce2.Program._RemoveFloatingHD = original
  def atom_lines(file_name):
    with open(file_name) as f:
      return [l for l in f if l.startswith(("ATOM", "HETATM", "ANISOU"))]
  assert atom_lines(fn) == atom_lines(fn_ref)

def test_floating_h_flipkins():
  # orphan water H (water 8 without its O), keep_existing_H=False: removed
  # in the main run and not brought back by the flipkin passes (the
  # keep_existing_H=True case is in test_flipkins_combined)
  pdb_str = "\n".join(l for l in cached_run(pdb_1yjp,
    exclude_water=False).text.splitlines()
    if " O   HOH A   8" not in l) + "\n"
  assert pdb_str.count("HOH A   8") == 2       # H1, H2
  def orphans(hierarchy):
    return [a for a in hierarchy.atoms() if a.parent().resname == "HOH" and
            a.parent().parent().resseq.strip() == "8"]
  seen = []
  original = Optimizers.Optimizer.__init__
  def init(self, probePhil, addFlipMovers, model, *args, **kwargs):
    seen.append(len(orphans(model.get_hierarchy())))
    return original(self, probePhil, addFlipMovers, model, *args, **kwargs)
  if not os.path.isdir("tst_r2n_optimizer_orphan_fk"):
    os.mkdir("tst_r2n_optimizer_orphan_fk")
  for keep in ("False",):
    del seen[:]
    Optimizers.Optimizer.__init__ = init
    try:
      fn, log, description = run_reduce2(
        "tst_r2n_optimizer_orphan_" + keep, pdb_str,
        ["approach=optimize", "keep_existing_H=" + keep, "exclude_water=False",
         "add_flip_movers=True",
         "output.flipkin_directory=tst_r2n_optimizer_orphan_fk"])
    finally:
      Optimizers.Optimizer.__init__ = original
    # main run and both flipkin configurations, none with the orphan H
    assert seen == [0, 0, 0], (keep, seen)
    assert description.count(FLOATING) == 1
    listed = description.split(FLOATING)[1].strip().splitlines()[:2]
    assert sorted(l.split()[0] for l in listed) == ["H1", "H2"], listed
    assert all(l.split()[1:] == ["HOH", "A", "8"] for l in listed), listed
    assert orphans(iotbx.pdb.input(fn).construct_hierarchy()) == []

def test_floating_h_water_hd():
  # 9k(b), hd_and_h: Tyr HH without its OH is dropped by the removal step;
  # the written water O-D stays at 0.98 A
  pdb_str = h_input_without_tyr_oh(exclude_water=False)
  fn, log, description = run_reduce2("tst_r2n_optimizer_float_wd", pdb_str,
    ["approach=optimize", "neutron_option=hd_and_h", "exclude_water=False"])
  assert FLOATING in description
  assert not [l for l in atom_labels(fn) if l[1] == "TYR" and
              l[2].strip() in ("HH", "DH")]
  assert water_groups(fn) == [("DOD", ["D", "O"])] * 7, water_groups(fn)
  check_water(fn, 0.98, 7)

def test_atom_characteristics_unbonded_h():
  # _GetAtomCharacteristics(): an H without bonded neighbour raises, unless
  # its residue is left unchanged (no restraints)
  from mmtbx.probe import Helpers
  from mmtbx.hydrogens import hd_exchange
  from scitbx.array_family import flex
  from libtbx import group_args
  pdb_str = h_input_without_tyr_oh()
  prog, log = make_program("tst_r2n_optimizer_gac", pdb_str)
  prog.validate()
  prog.model = prog.data_manager.get_model()
  prog._ReinterpretModel()
  carts = flex.vec3_double([a.xyz for a in prog.model.get_atoms()])
  bonds, asu = prog.model.get_restraints_manager().geometry.\
    get_all_bond_proxies(sites_cart=carts)
  neighbors = Helpers.getBondedNeighborLists(prog.model.get_atoms(), bonds,
    asu)
  try:
    prog._GetAtomCharacteristics(neighbors)
  except Sorry as e:
    assert str(e) == "Found Hydrogen with no neigbors.", str(e)
  else:
    raise AssertionError("unbonded HH accepted")
  keys = hd_exchange.atom_keys(prog.model.get_hierarchy())
  hh = [k for k, a in zip(keys, prog.model.get_atoms())
        if a.name.strip() == "HH"]
  assert len(hh) == 1
  prog._unrestrained = group_args(keys=set(hh), labels=["TYR A 7"])
  in_water, in_het, in_main, in_side = prog._GetAtomCharacteristics(neighbors)
  a = [a for a in prog.model.get_atoms() if a.name.strip() == "HH"][0]
  assert in_main[a] is False and not in_water[a]

def test_optimizer_single_rotator_d():
  # Ser with DG (no reduce2 normalization) gets a MoverSingleHydrogenRotator,
  # as with HG
  from mmtbx.reduce import Movers
  text = placed_text("".join(l for l in pdb_ser_dg.splitlines(True)
    if l[12:16].strip() != "DG"))
  for name, model_text in (("HG", text), ("DG", h_to_d(text, ["HG"]))):
    model = processed_model(model_text)
    model.setup_riding_h_manager()
    atoms, bonded, extra = probe_setup(model)
    atom = [a for a in atoms if a.name.strip() == name][0]
    mover = Movers.MoverSingleHydrogenRotator(atom, bonded, extra,
      model.get_riding_h_manager().h_parameterization)
    assert len(mover.CoarsePositions().positions) == 7, name

def test_phantom_targets_skip_d():
  # getPhantomHydrogensFor(): a D next to the water O is skipped as an
  # orientation target exactly as an H is (Optimizer call, acceptorOnly=False)
  from mmtbx.probe import Helpers
  import mmtbx_probe_ext as probe_ext
  text = placed_text("".join(l for l in pdb_ser_dg.splitlines(True)
    if l[12:16].strip() != "DG"))
  model = processed_model(text)
  og, hg = [[a for a in model.get_atoms() if a.name.strip() == n][0]
            for n in ("OG", "HG")]
  # water O 2.8 A from OG, on the OG-HG line
  o = matrix.col(og.xyz) + 2.8 * (matrix.col(hg.xyz) -
    matrix.col(og.xyz)).normalize()
  text = text.replace("END", "HETATM  999  O   HOH B   1    %8.3f%8.3f%8.3f"
    "  1.00 20.00           O\nEND" % tuple(o))
  phantoms = {}
  for name, model_text in (("HG", text), ("DG", h_to_d(text, ["HG"]))):
    model = processed_model(model_text)
    atoms, bonded, extra = probe_setup(model)
    water = [a for a in atoms if a.parent().resname == "HOH"][0]
    result = Helpers.getPhantomHydrogensFor(Helpers.getMaxISeq(model), water,
      probe_ext.SpatialQuery(atoms), extra, 0.0, False)
    phantoms[name] = sorted(tuple(round(x, 3) for x in p.xyz) for p in result)
  assert len(phantoms["HG"]) > 0
  assert phantoms["DG"] == phantoms["HG"], phantoms

def test_probe2_deuterated():
  # probe2 on the fully deuterated model gives the contacts of the H model
  import iotbx.cli_parser
  from iotbx.data_manager import DataManager
  from mmtbx.programs import probe2
  def run_probe2(model, output_format):
    parser = iotbx.cli_parser.CCTBXParser(program_class=probe2.Program,
      logger=null_out())
    parser.parse_args(["output.write_files=False",
      "output.format=" + output_format])
    dm = DataManager()
    dm.add_model("model", model)
    p2 = probe2.Program(dm, parser.working_phil.extract(),
      master_phil=parser.master_phil, logger=null_out())
    p2.overrideModel(model)
    return p2.run()[1]
  def as_h(line):
    # raw line with D atom names written as H (the element columns give the
    # parent's element for H and D)
    f = line.split(":")
    for i in (3, 4):
      if f[i][-5:].lstrip().startswith("D"):
        f[i] = f[i][:-5] + f[i][-5:].replace("D", "H", 1)
    return ":".join(f)
  text = placed_text(pdb_1yjp)
  h_model, d_model = processed_model(text), processed_model(h_to_d(text))
  assert set(a.element for a in d_model.get_atoms()
             if a.element_is_hydrogen()) == set(["D"])
  assert run_probe2(h_model, "oneline") == run_probe2(d_model, "oneline")
  h_raw = run_probe2(h_model, "raw").splitlines()
  d_raw = run_probe2(d_model, "raw").splitlines()
  assert len(h_raw) > 100, len(h_raw)
  assert d_raw != h_raw  # D names in the D contacts
  assert [as_h(l) for l in d_raw] == h_raw

if __name__ == "__main__":
  test_floating_h_removed()
  test_optimize_without_h()
  test_floating_h_water_098()
  test_no_floating_h_unchanged()
  test_floating_h_flipkins()
  test_floating_h_water_hd()
  test_atom_characteristics_unbonded_h()
  test_optimizer_single_rotator_d()
  test_phantom_targets_skip_d()
  test_probe2_deuterated()
  print("OK")
