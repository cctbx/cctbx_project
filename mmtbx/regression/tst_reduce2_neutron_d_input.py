"""reduce2 neutron tests: input H/D: normalization, recorded values, Sorrys,
approach=remove"""
from __future__ import absolute_import, division, print_function
import iotbx.pdb
from libtbx.utils import Sorry
from mmtbx.regression import reduce2_neutron_utils
from mmtbx.regression.reduce2_neutron_utils import (
  atom_records, cached_run, check_removed, cif_lig, dg_moved, dipeptide,
  dipeptide_hd, dipeptide_out, edit_atoms, hd_labels, hd_table, hd_values,
  his_no_h, make_model_from_str, make_program, pairs_03_07, pdb_1yjp,
  pdb_1yjp_dod_d, pdb_1yjp_h, pdb_lig, pdb_ser_dg, run_prog, run_reduce2,
  ser_cif_hg_cb)

reduce2_neutron_utils.set_prefix("tst_r2n_d_input")

def test_ambiguous_parent():
  # 5b: an H/D with two covalent partners stops the run before the Optimizer.
  # reduce2 takes no geometry_restraints.edits, and placement drops an H
  # that a CIF bonds to two heavy atoms, so the check is called directly on
  # a model whose edits make Ser HG ambiguous
  from mmtbx.hydrogens import tst_hd_exchange as t
  model = t.make_model(t.pdb_ser_hd.replace(
    "ATOM     11  DG BSER A   2      15.856   0.874   6.441  0.50 20.00"
    "           D\n", ""), [t.edits_bond("name HG", "name CB", 1.1)])
  prog, log = make_program("tst_r2n_d_input_amb", pdb_1yjp,
    neutron_option="hd_and_h")
  prog.validate()
  prog.model = model
  try:
    prog._CheckAmbiguousParents()
  except Sorry as e:
    text = str(e)
    assert "more than one covalent bond partner" in text
    assert 'HG ASER' in text and 'OG  SER' in text and 'CB  SER' in text
    assert "covalent geometry" in text and "edits" in text
  else:
    raise AssertionError("ambiguous HG not reported")
  # a normal model passes the check
  prog.model = make_model_from_str(pdb_1yjp_h())
  prog._CheckAmbiguousParents()
  # end to end, a CIF bond from HO1 to a second heavy atom: placement drops
  # HO1 (not parameterizable), so the run completes with none and hd_and_h
  cif = cif_lig.replace(" 01  O1  HO1  single  0.850   0.020  0.980\n",
    " 01  O1  HO1  single  0.850   0.020  0.980\n"
    " 01  C1  HO1  single  1.050   0.020  1.050\n")
  assert cif != cif_lig
  for option in ("none", "hd_and_h"):
    prog = run_prog("tst_r2n_d_input_amb_" + option, pdb_lig, cif,
      neutron_option=option, keep_existing_H=True)
    h = iotbx.pdb.input(prog.params.output.filename).construct_hierarchy()
    assert "HO1" not in [a.name.strip() for a in h.atoms()]

def test_record_restore():
  # 3: refined values of input H/D pairs are kept, not reset to 0.5/0.5;
  # written as one isotope: parent occupancy, that isotope's B
  src = pairs_03_07()
  exchangeable = [("37", "H"), ("37", "HG"), ("38", "H")]
  fn, log, description = run_reduce2("tst_r2n_d_input_rec_hd", src,
    ["approach=optimize", "neutron_option=hd_and_h"])
  v = hd_values(fn)
  for resseq, name in exchangeable:
    d_name = name.replace("H", "D", 1)
    assert v[(resseq, name, "A", "H")] == (0.3, 15.0), v
    assert v[(resseq, d_name, "B", "D")] == (0.7, 25.0), v
  assert "Values restored from the input: 3" in description
  for option, element, b in (("all_d", "D", 25.0), ("all_h", "H", 15.0)):
    fn, log, description = run_reduce2("tst_r2n_d_input_rec_" + option,
      src, ["approach=optimize", "neutron_option=" + option])
    v = hd_values(fn)
    for resseq, name in exchangeable:
      if element == "D": name = name.replace("H", "D", 1)
      assert v[(resseq, name, "", element)] == (1.0, b), (option, v)
    assert not [k for k in v if k[2] != ""], v
    assert "Values restored from the input: 3" in description

def test_pair_sum_note():
  # 4: report, don't repair. Ser at 0.6 whose pairs have H + D = 1 (0.5 +
  # 0.5): the run completes, the values are unchanged and the sum is noted.
  # This tests that reduce2 passes such input through with a note; it does
  # not say that H + D = 1 on a parent at 0.6 is right
  def edit(a):
    if a.parent().parent().resseq.strip() == "37" and \
       a.parent().altloc == "":
      a.set_occ(0.6)
  src = edit_atoms(dipeptide_out("hd_and_h"), edit)
  fn, log, description = run_reduce2("tst_r2n_d_input_sum", src,
    ["approach=optimize", "neutron_option=hd_and_h"])
  v = hd_values(fn)
  for name in ("H", "HG"):
    assert v[("37", name, "A", "H")][0] == 0.5, v
    assert v[("37", name.replace("H", "D", 1), "B", "D")][0] == 0.5, v
  notes = description.split("Input H/D notes (values not changed)")[1]
  assert notes.count("H + D occupancy 1.00, parent 0.60 (sum 1)") == 2, notes
  # listed as split sites on a partial parent, with the remark that restored
  # pairs keep their input values (0.5/0.5 above, not 0.3/0.3)
  partial = description.split(
    "Split sites on a partially occupied parent (2):")[1].splitlines()
  assert "restored pairs keep their input values" in partial[1], partial
  for line, name in zip(partial[2:4], ("H  ASER A  37", "HG ASER A  37")):
    assert line.strip().startswith(name), partial
    assert line.endswith("parent occupancy 0.60"), partial
  assert "Values restored from the input: 3" in description

def test_invalid_input():
  # 5: a matcher problem stops the run when the input H/D are kept, with
  # the parameters that rebuild them for this run
  src = dg_moved()
  for args, hint in (
      (["approach=optimize"], "approach=add keep_existing_H=False"),
      (["approach=optimize", "keep_existing_H=False"],
       "approach=add keep_existing_H=False"),
      (["approach=add", "keep_existing_H=True"], "keep_existing_H=False")):
    try:
      run_reduce2("tst_r2n_d_input_inv", src,
        args + ["neutron_option=hd_and_h"])
    except Sorry as e:
      text = str(e)
      assert "H/D pairs that cannot be collapsed (1)" in text, text
      assert "HG ASER A  37" in text and "DG BSER A  37" in text, text
      assert "H and D 0.50 A apart" in text, text
      assert text.endswith("Rebuild H: " + hint), text
    else:
      raise AssertionError("ran: %s" % args)
  # discarding the input H/D: runs, Ser HG/DG split as new
  fn, log, description = run_reduce2("tst_r2n_d_input_inv_ok", src,
    ["approach=add", "keep_existing_H=False", "neutron_option=hd_and_h"])
  assert hd_labels(fn) == hd_table(dipeptide_hd["hd_and_h"])

def test_floating_h_hd_and_h():
  # 5a, hd_and_h with D input: Ser OG missing, so HG and DG have no covalent
  # partner; they are dropped and listed, the other pairs are kept as H/D
  src = "\n".join(l for l in dipeptide_out("hd_and_h").splitlines()
                  if l[12:16] != " OG ") + "\n"
  fn, log, description = run_reduce2("tst_r2n_d_input_flhd", src,
    ["approach=optimize", "neutron_option=hd_and_h"])
  dropped = description.split(
    "H/D without a covalent partner, dropped at the input check (2):")[1]
  assert "HG ASER A  37" in dropped and "DG BSER A  37" in dropped, dropped
  expected = set(x for x in hd_table(dipeptide_hd["hd_and_h"])
                 if x[1] not in ("HG", "DG"))
  assert hd_labels(fn) == expected, sorted(hd_labels(fn) ^ expected)

def test_ambiguous_parent_optimize():
  # 5b, optimize end to end: H-only input (no D, so no input check), a user
  # SER CIF bonds HG to OG and CB; hd_and_h stops before the Optimizer,
  # naming HG and both restraints; none runs
  fn, _, _ = run_reduce2("tst_r2n_d_input_ambo_in", dipeptide, [])
  with open(fn) as f: src = f.read()
  cif = ser_cif_hg_cb()
  try:
    run_prog("tst_r2n_d_input_ambo", src, cif, approach="optimize",
      neutron_option="hd_and_h")
  except Sorry as e:
    text = str(e)
    assert "more than one covalent bond partner" in text, text
    assert "HG  SER A  37" in text, text
    assert "OG  SER A  37" in text and "CB  SER A  37" in text, text
    assert text.count("covalent geometry") == 2, text
  else:
    raise AssertionError("ambiguous HG not reported")
  prog = run_prog("tst_r2n_d_input_ambo_none", src, cif,
    approach="optimize", neutron_option="none")
  assert "Optimizing" in prog.log_text

def test_d_input_none():
  # D input with none: optimize and add keep_existing_H=True would keep the
  # D -> Sorry with the hint for the run; add keep_existing_H=False discards
  # them and runs
  for args, hint in (
      (["approach=optimize"], "approach=add keep_existing_H=False"),
      (["approach=add", "keep_existing_H=True"], "keep_existing_H=False")):
    try:
      run_reduce2("tst_r2n_d_input_dnone", pdb_ser_dg, args)
    except Sorry as e:
      assert str(e) == ("input contains deuterium; set neutron_option, or "
        "rebuild H: " + hint), str(e)
    else:
      raise AssertionError("ran: %s" % args)
  fn, log, description = run_reduce2("tst_r2n_d_input_dnone_ok",
    pdb_ser_dg, ["approach=add", "keep_existing_H=False"])
  assert not [k for k in hd_values(fn) if k[3] == "D"]

def test_d_input_qmready():
  # qmready forces keep_existing_H=True and allows only none: D input stops,
  # with the hint to remove the H/D first
  try:
    run_reduce2("tst_r2n_d_input_dqm", pdb_ser_dg,
      ["approach=add", "qmready=True"])
  except Sorry as e:
    assert str(e) == ("input contains deuterium; set neutron_option, or "
      "rebuild H: remove the H/D first (approach=remove), then run "
      "qmready=True"), str(e)
  else:
    raise AssertionError("D input with qmready ran")

def test_unused_records():
  # an input pair on His HE2: the Optimizer flips the His and deletes HE2,
  # so its recorded values are not used and listed by site
  src = cached_run(his_no_h, neutron_option="hd_and_h").text
  fn, log, description = run_reduce2("tst_r2n_d_input_unused", src,
    ["approach=optimize", "neutron_option=hd_and_h", "add_flip_movers=True"])
  assert "Deleting chain A HIS 45 HE2" in description, description
  unused = description.split(
    "Input values not used (site no longer in the model) (1):")[1]
  assert unused.splitlines()[1].strip() == "HE2 HIS A 45", unused
  assert not [k for k in hd_values(fn) if k[1] in ("HE2", "DE2")]
  assert "Values restored from the input: 2" in description

def test_remove_d_input():
  # 9a: approach=remove is its own path: D input, an H/D layout the matcher
  # rejects (Leu HA and DA in one atom_group) and a DOD water with D1/D2,
  # also with keep_existing_H=True -> no Sorry, every H/D removed
  da = [l for l in dipeptide_out("hd_and_h").splitlines()
        if l[12:16] == " HA " and l[17:26] == "LEU A  38"]
  assert len(da) == 1
  src = dipeptide_out("hd_and_h").replace(da[0] + "\n",
    da[0] + "\n" + da[0][:12] + " DA " + da[0][16:76] + " D\n")
  src = src.replace("END\n", "") + "\n".join(l for l in
    pdb_1yjp_dod_d.splitlines() if " DOD A   8" in l) + "\nEND\n"
  assert src.count("DOD A   8") == 3 and src.count(" DA  LEU") == 1
  with open("tst_r2n_d_input_rmd_in.pdb", "w") as f:
    f.write(src)
  # the layout is rejected when the H/D are kept
  try:
    run_reduce2("tst_r2n_d_input_rmd_kept", src,
      ["approach=add", "keep_existing_H=True", "neutron_option=hd_and_h"])
  except Sorry as e:
    assert "H and D of one site in one conformer" in str(e), str(e)
  else:
    raise AssertionError("HA/DA in one atom_group accepted")
  for args in ([], ["keep_existing_H=True"]):
    fn, log, description = run_reduce2("tst_r2n_d_input_rmd", src,
      ["approach=remove"] + args)
    check_removed("tst_r2n_d_input_rmd_in.pdb", fn)
    h = iotbx.pdb.input(fn).construct_hierarchy()
    assert [ag.resname for ag in h.atom_groups()
            if ag.parent().resseq.strip() == "8"] == ["HOH"]
    # Ser and Leu: only the blank group is left
    assert sorted(set(ag.altloc for ag in h.atom_groups())) == [""]

def test_remove_then_qmready():
  # the qmready hint: hd_and_h output -> approach=remove -> qmready=True runs
  # free His: atom-complete, as qmready requires (the dipeptide is not)
  src = cached_run(his_no_h, neutron_option="hd_and_h").text
  assert " DD1BHIS" in src
  fn, log, description = run_reduce2("tst_r2n_d_input_rmq", src,
    ["approach=remove"])
  with open(fn) as f: src = f.read()
  fn, log, description = run_reduce2("tst_r2n_d_input_rmq2", src,
    ["approach=add", "qmready=True"])
  assert not [k for k in hd_values(fn) if k[3] == "D"]
  assert hd_values(fn)

def test_remove_h_input():
  # H-only input: the output is the input without its H, atom for atom
  from mmtbx.hydrogens import tst_hd_exchange as th
  fn_conf, _, _ = run_reduce2("tst_r2n_d_input_rmh_conf",
    th.pdb_pep_og_ab, [])
  with open(fn_conf) as f: conf = f.read()
  for key, src in (("1yjp", pdb_1yjp_h()), ("conf", conf)):
    fn_in = "tst_r2n_d_input_rmh_%s_in.pdb" % key
    with open(fn_in, "w") as f: f.write(src)
    fn, log, description = run_reduce2("tst_r2n_d_input_rmh_" + key, src,
      ["approach=remove"])
    assert atom_records(fn) == [r for r in atom_records(fn_in)
                                if r[4] != "H"], key

if __name__ == "__main__":
  test_ambiguous_parent()
  test_record_restore()
  test_pair_sum_note()
  test_invalid_input()
  test_floating_h_hd_and_h()
  test_ambiguous_parent_optimize()
  test_d_input_none()
  test_d_input_qmready()
  test_unused_records()
  test_remove_d_input()
  test_remove_then_qmready()
  test_remove_h_input()
  print("OK")
