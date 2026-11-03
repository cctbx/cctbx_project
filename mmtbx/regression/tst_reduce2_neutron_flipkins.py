"""reduce2 neutron tests: flipkins (amide and His) on one input with every
feature the flipkin passes must handle"""
from __future__ import absolute_import, division, print_function
import os
import iotbx.pdb
from mmtbx.reduce import Optimizers
from mmtbx.regression import reduce2_neutron_utils
from mmtbx.regression.reduce2_neutron_utils import (
  FLOATING, INFO_FALSE, INFO_TRUE, UNCHANGED, cached_run,
  combined_flipkin_model, ligand_split_records, make_program, water_groups,
  zq3_cif, zq4_cif)

reduce2_neutron_utils.set_prefix("tst_r2n_flipkins")

def test_flipkins_combined():
  # Flipkins (amide and His) on one input, add and optimize, none and
  # hd_and_h, keep_existing_H=True, ignore_missing_restraints=True:
  # - the flipkin Optimizer runs use the neutron flag (none + neutron
  #   distances: every configuration logs it)
  # - flipkins are made from the H-only self.model: hd_and_h and none with
  #   neutron distances give the same flipkins
  # - an O-only water named DOD gets H in every pass (renamed HOH)
  # - an orphan water H removed as floating stays out of every flipkin pass;
  #   optimize: a water H the dictionary does not know (H3), which placement
  #   in the passes keeps, too: removed from the snapshot
  #   (_RemoveFromSnapshot())
  # - 04 left unchanged: its unbonded H/D are classified by their residue
  #   (_GetAtomCharacteristics()) and kept out of probe2 (_ProbeModel());
  #   written exactly as in the input, listed
  # - BGM (CCD only) is bonded in the passes (restraints attached once);
  #   amide and His flipkins written; BGM H/D in the output
  # - no auto restraint objects in the written files, the results model and
  #   self.model; the data manager's input model keeps exactly its own
  # - D input (optimize hd_and_h) and add hd_and_h: the flipkin passes do not
  #   change the written model or the H/D section (run without flipkins)
  model = combined_flipkin_model()
  src = {}
  for option in ("none", "hd_and_h"):
    src[option] = cached_run(model, (zq3_cif, zq4_cif), approach="place",
      neutron_option=option, exclude_water=False).text
  # H input: water 8 without its O (orphan H1/H2), water 9 O only, as DOD
  lines = []
  for l in src["none"].splitlines():
    if l[17:26] == "HOH A   8" and l[12:16] == " O  ": continue
    if l[17:26] == "HOH A   9":
      if l[76:78].strip() == "H": continue
      l = l[:17] + "DOD" + l[20:]
    lines.append(l)
  h_src = "\n".join(lines) + "\n"
  assert h_src.count("HOH A   8") == 2 and h_src.count("DOD A   9") == 1
  # both His ring N-H (a His with one gets no flip mover on optimize)
  assert " HD1 HIS" in h_src and " HE2 HIS" in h_src
  # optimize: also an H3 on water 10 (no restraint, floating)
  o = [l for l in lines if l[17:26] == "HOH A  10" and l[12:16] == " O  "][0]
  h3 = o[:12] + " H3 " + o[16:30] + "%8.3f" % (float(o[30:38]) + 0.96) + \
    o[38:76] + " H"
  opt_src = h_src.replace(o + "\n", o + "\n" + h3 + "\n")
  assert opt_src.count(" H3  HOH A  10") == 1
  d_src = src["hd_and_h"]
  assert "DD21BASN" in d_src
  def orphans(hierarchy):
    return [a for a in hierarchy.atoms() if a.parent().resname == "HOH" and
            a.parent().parent().resseq.strip() == "8"]
  seen = []
  original_init = Optimizers.Optimizer.__init__
  def init(self, probePhil, addFlipMovers, model, *args, **kwargs):
    seen.append(len(orphans(model.get_hierarchy())))
    return original_init(self, probePhil, addFlipMovers, model, *args,
      **kwargs)
  names = lambda m: [n for n, o in (m.get_restraint_objects() or [])]
  def run(prefix, source, ext, directory, **params):
    params.update(keep_existing_H=True, ignore_missing_restraints=True,
      exclude_water=False, add_flip_movers=True,
      output__filename=prefix + "_out" + ext,
      output__description_file_name=prefix + ".txt", output__overwrite=True)
    if directory is not None:
      if not os.path.isdir(directory): os.mkdir(directory)
      for f in os.listdir(directory): os.remove(os.path.join(directory, f))
      params["output__flipkin_directory"] = directory
    prog, log = make_program(prefix, source, zq3_cif, **params)
    input_model = prog.data_manager.get_model()
    before = names(input_model)
    assert len(before) == 1, before
    prog.validate()
    prog.run()
    prog.log_text = log.getvalue()
    assert names(prog.data_manager.get_model()) == before
    assert prog.data_manager.get_model() is input_model
    for m in (prog.get_results().model, prog.model):
      assert not [n for n in names(m) if n.startswith("auto_")], names(m)
    with open(prog.params.output.filename) as f: prog.text = f.read()
    with open(prog.params.output.description_file_name) as f:
      prog.description = f.read()
    assert "auto_" not in prog.text and "data_comp_BGM" not in prog.text
    assert "BGM" in prog.text
    return prog
  def hd_section(description):
    return description[description.index("H/D model (neutron_option="):]
  kins = {}
  for key, approach, option, source, ext, extra in (
      ("add_none", "add", "none", h_src, ".pdb",
       dict(use_neutron_distances=True)),
      ("add_hd", "add", "hd_and_h", h_src, ".cif", {}),
      ("opt_none", "optimize", "none", opt_src, ".cif", {}),
      ("opt_hd", "optimize", "hd_and_h", d_src, ".pdb", {})):
    prefix = "tst_r2n_flipkins_fkc_" + key
    directory = prefix + "_fk"
    del seen[:]
    Optimizers.Optimizer.__init__ = init
    try:
      prog = run(prefix, source, ext, directory, approach=approach,
        neutron_option=option, **extra)
    finally:
      Optimizers.Optimizer.__init__ = original_init
    label = (approach, option)
    h = iotbx.pdb.input(source_info=None, lines=prog.text).construct_hierarchy()
    assert sorted(os.listdir(directory)) == [prefix + "FH-fliphis.kin",
      prefix + "FH-flipnq.kin"], (label, os.listdir(directory))
    kins[key] = {}
    for f in sorted(os.listdir(directory)):
      with open(os.path.join(directory, f)) as fh:
        text = fh.read().replace(prefix + "_out" + ext, "BASE_out").replace(
          prefix, "BASE")
      kins[key][f.replace(prefix, "BASE")] = [l for l in text.splitlines()
        if not l.startswith("@caption")]
    # BGM bonded in the passes, its H/D written
    assert len([a for a in h.atoms() if a.parent().resname == "BGM" and
      a.element.strip() in ("H", "D")]) == (13 if option == "none" else 19)
    # 04 left unchanged and listed
    assert ligand_split_records(prog.text, True) == \
      ligand_split_records(source, True), label
    for text in (prog.description, prog.log_text):
      assert UNCHANGED + "\n04 D 1\n" in text, label
    has_d = "D" in [a.element.strip() for a in h.atoms()]
    assert has_d == (option != "none"), label
    # no orphan water H in any Optimizer run, nor written
    assert orphans(h) == [], label
    if key == "add_none":
      assert INFO_TRUE in prog.description
      assert prog.log_text.count(INFO_TRUE) >= 2, prog.log_text.count(INFO_TRUE)
      assert INFO_FALSE not in prog.log_text
      assert water_groups(prog.params.output.filename) == \
        [("HOH", ["H", "O"])] * 6, water_groups(prog.params.output.filename)
    if key == "add_hd":
      assert kins["add_hd"] == kins["add_none"]
    if key == "opt_none":
      # main run and the four flipkin configurations (amide, His)
      assert seen == [0] * 5, seen
      assert prog.description.count(FLOATING) == 1
      listed = prog.description.split(FLOATING)[1].strip().splitlines()[:2]
      assert sorted(l.split()[0] for l in listed) == ["H1", "H2"], listed
      assert all(l.split()[1:] == ["HOH", "A", "8"] for l in listed), listed
      listed = prog.description.split(FLOATING)[1].strip().splitlines()[2]
      assert listed.split() == ["H3", "HOH", "A", "10"], listed
    if option == "hd_and_h":
      # the same run without flipkins writes the same model and H/D section
      off = run(prefix + "_off", source, ext, None, approach=approach,
        neutron_option=option, **extra)
      assert off.text.replace(prefix + "_off", prefix) == prog.text, label
      assert hd_section(off.description) == hd_section(prog.description)
    if key == "opt_hd":
      assert "Values restored from the input: 29" in prog.description, \
        hd_section(prog.description)

if __name__ == "__main__":
  test_flipkins_combined()
  print("OK")
