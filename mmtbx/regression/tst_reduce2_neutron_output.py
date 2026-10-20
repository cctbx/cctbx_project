"""reduce2 neutron tests: neutron_option output: names, Phil, per-site rules,
labels per mode, output checks, mmCIF, description file"""
from __future__ import absolute_import, division, print_function
import os
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out, Sorry
from libtbx.test_utils import approx_equal
from mmtbx.regression import reduce2_neutron_utils
from mmtbx.regression.reduce2_neutron_utils import (
  NEUTRON, OPTIONS, atom_label_set, cached_run, check_bond_lengths, check_no_d,
  check_sorry, cif_lig, dipeptide, dipeptide_hd, hc4_all_d, hc4_no_h,
  hd_labels, hd_table, heavy_atoms, heavy_pdb, his_hd, his_no_h, make_program,
  pdb_1yjp, pdb_lig, pdb_ser_dg, run_prog, run_reduce2, zq3_atoms, zq3_cif)

reduce2_neutron_utils.set_prefix("tst_r2n_output")

def test_output_names():
  # 13: default output and description file names
  for option in OPTIONS:
    for flips, pad in ((False, "H"), (True, "FH")):
      prog, log = make_program("tst_r2n_output_names", pdb_1yjp,
        neutron_option=option, add_flip_movers=flips)
      prog.validate()
      suffix = "" if option == "none" else "_" + option
      expected = "tst_r2n_output_names" + pad + suffix
      assert prog.params.output.filename == expected + ".pdb", (
        option, flips, prog.params.output.filename)
      assert prog.params.output.description_file_name == expected + ".txt"

def test_validate_sorrys():
  prefix = "tst_r2n_output_validate"
  for option in OPTIONS[1:]:
    check_sorry(prefix, pdb_1yjp, "cannot be used with approach=remove",
      neutron_option=option, approach="remove")
    check_sorry(prefix, pdb_1yjp, "cannot be used with qmready=True",
      neutron_option=option, qmready=True)
    # optimize, and D input with keep_existing_H=True, pass validate(); the
    # input check in run() decides (step 4)
    for pdb_str, params in ((pdb_1yjp, dict(approach="optimize")),
                            (pdb_1yjp, dict(keep_existing_H=True)),
                            (pdb_ser_dg, dict(keep_existing_H=True))):
      prog, log = make_program(prefix, pdb_str, neutron_option=option,
        **params)
      prog.validate()
  # none: the same runs pass validate()
  for params in (dict(approach="remove"), dict(qmready=True),
                 dict(approach="optimize"), dict(keep_existing_H=True)):
    prog, log = make_program(prefix, pdb_ser_dg, neutron_option="none",
      **params)
    prog.validate()

def test_rebuild_model():
  # _RebuildModel(): user CIF restraints, crystal symmetry, atom count
  prog, log = make_program("tst_r2n_output_rebuild", pdb_lig, cif_lig,
    stop_on_any_missing_hydrogen=True)
  prog.validate()
  prog.model = prog.data_manager.get_model()
  hierarchy = prog.model.get_hierarchy().deep_copy()
  for make_restraints in (True, False):
    rebuilt = prog._RebuildModel(hierarchy.deep_copy(), make_restraints)
    assert rebuilt.get_number_of_atoms() == prog.model.get_number_of_atoms()
    assert rebuilt.crystal_symmetry().is_similar_symmetry(
      prog.model.crystal_symmetry())
    assert (rebuilt.get_restraints_manager() is not None) == make_restraints
  # the CIF's C1-O1 target is used: the restraints come from the user CIF
  simple, asu = prog._RebuildModel(hierarchy.deep_copy(),
    True).get_restraints_manager().geometry.get_all_bond_proxies()
  ideal = sorted(round(p.distance_ideal, 4) for p in simple)
  assert 1.4321 in ideal, ideal
  # without the CIF the ligand cannot be interpreted
  prog, log = make_program("tst_r2n_output_rebuild_nocif", pdb_lig,
    stop_on_any_missing_hydrogen=True)
  prog.validate()
  prog.model = prog.data_manager.get_model()
  try:
    prog._RebuildModel(prog.model.get_hierarchy().deep_copy(), True)
  except Sorry:
    pass
  else:
    raise AssertionError("01 interpreted without its CIF")

def test_n_terminal_nh3():
  # 1: Gly 1 H1, H2, H3 all split in hd_and_h
  fn = cached_run(pdb_1yjp, neutron_option="hd_and_h").params.output.filename
  h = iotbx.pdb.input(fn).construct_hierarchy()
  gly = [rg for rg in h.residue_groups() if rg.resseq.strip() == "1"][0]
  by_alt = dict((ag.altloc, sorted(a.name.strip() for a in ag.atoms()))
                for ag in gly.atom_groups())
  assert by_alt["A"] == ["H1", "H2", "H3"], by_alt
  assert by_alt["B"] == ["D1", "D2", "D3"], by_alt

def test_conformers():
  # 8: Ser with OG in A and B: backbone H on the blank N split into C/D,
  # HG in A and B written as D and listed as conformer sites
  from mmtbx.hydrogens.tst_hd_exchange import pdb_pep_og_ab
  for approach in ("add", "place"):
    prog = run_prog("tst_r2n_output_conf_" + approach, pdb_pep_og_ab,
      neutron_option="hd_and_h", approach=approach)
    h = iotbx.pdb.input(prog.params.output.filename).construct_hierarchy()
    ser = [rg for rg in h.residue_groups() if rg.resseq.strip() == "2"][0]
    hd = dict((ag.altloc, sorted(a.name.strip() for a in ag.atoms()
               if a.element.strip() in ("H", "D"))) for ag in ser.atom_groups())
    assert hd["C"] == ["H"] and hd["D"] == ["D"], hd
    assert "DG" in hd["A"] and "DG" in hd["B"], hd
    sites = prog._neutron_result.conformer_sites
    assert sorted(x.site[5:9].strip() for x in sites
                  if x.reason == "conformer") == ["HG", "HG"], sites
    # 8 split sites and the two conformer D (water excluded)
    assert prog._neutron_result.n_exchangeable == 10
    with open(prog.params.output.description_file_name) as f:
      assert "Exchangeable H sites (not water): 10" in f.read()

def test_output_symmetry():
  # 9d: none in -> none out; cryo-EM placeholder -> 1 1 1 P 1; cell -> cell
  no_cs = "\n".join(l for l in pdb_1yjp.splitlines()
                    if not l.startswith("CRYST1"))
  placeholder = pdb_1yjp.replace(
    "CRYST1   21.937    4.866   23.477  90.00 107.08  90.00 P 1 21 1",
    "CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1")
  for label, pdb_str, cryst1 in (("nocs", no_cs, None),
      ("placeholder", placeholder, "CRYST1    1.000    1.000    1.000"),
      ("cell", pdb_1yjp, "CRYST1   21.937    4.866   23.477")):
    if label == "cell":
      fn = cached_run(pdb_str,
        neutron_option="hd_and_h").params.output.filename
    else:
      fn, log, description = run_reduce2("tst_r2n_output_cs_" + label,
        pdb_str, ["neutron_option=hd_and_h"])
    with open(fn) as f:
      lines = [l for l in f.read().splitlines() if l.startswith("CRYST1")]
    if cryst1 is None:
      assert lines == [], lines
    else:
      assert len(lines) == 1 and lines[0].startswith(cryst1), lines

def test_output_model():
  # 9e, 9f, 9h: results and heavy atoms; 9g: no files, and the results model
  # is the H/D model
  prefix = "tst_r2n_output_out"
  for ext in ("_out.pdb", ".txt"):
    if os.path.exists(prefix + ext): os.remove(prefix + ext)
  prog = run_prog(prefix, pdb_1yjp, neutron_option="hd_and_h",
    exclude_water=False, output__write_files=False)
  assert not os.path.exists(prefix + "_out.pdb")
  assert not os.path.exists(prefix + ".txt")
  results = prog.get_results()
  assert results.model is prog._output_model
  assert hasattr(results, "total_charge")
  out_h = results.model.get_hierarchy()
  assert "D" in [a.element.strip() for a in out_h.atoms()]
  # self.model stays H-only
  assert "D" not in [a.element.strip()
                     for a in prog.model.get_hierarchy().atoms()]
  # 9h: heavy atoms identical (xyz to round-off from rebuilding the model)
  out_heavy = heavy_atoms(out_h)
  in_heavy = heavy_atoms(prog.model.get_hierarchy())
  assert sorted(out_heavy) == sorted(in_heavy)
  for key, (name, altloc, xyz, occ, b, uij) in in_heavy.items():
    o = out_heavy[key]
    assert (o[0], o[1], o[3], o[4], o[5]) == (name, altloc, occ, b, uij), key
    assert approx_equal(o[2], xyz, eps=1e-9, out=None), key
  # 9e: TYR ring names as in self.model
  def ring(h):
    return dict((a.name.strip(), a.xyz) for a in h.atoms()
      if a.parent().resname == "TYR" and a.name.strip() in
      ("CD1", "CD2", "CE1", "CE2"))
  assert ring(out_h) == ring(prog.model.get_hierarchy())

def test_occupancy_groups():
  # 10: one exchangeable pair per split site; phenix.refine's
  # occupancy_selections() puts each split site in one two-atom group
  from mmtbx.refinement.occupancies import occupancy_selections
  prog = cached_run(pdb_1yjp, neutron_option="hd_and_h")
  fn = prog.params.output.filename
  model = mmtbx.model.manager(model_input=iotbx.pdb.input(fn), log=null_out())
  model.process(make_restraints=True)
  h = model.get_hierarchy()
  pairs = h.exchangeable_hd_selections()
  assert len(pairs) == prog._neutron_result.n_split
  groups = occupancy_selections(model=model, as_flex_arrays=False)
  for x in pairs:
    i_h, i_d = x[0][0], x[1][0]
    found = [g for g in groups if any(i_h in s for s in g)]
    assert len(found) == 1, (h.atoms()[i_h].id_str(), found)
    members = sorted(i for s in found[0] for i in s)
    assert members == sorted([i_h, i_d]), (h.atoms()[i_h].id_str(), members)

def test_partial_parent():
  # 12: Tyr OH at 0.8 -> HH 0.4 / DH 0.4, listed
  pdb_str = pdb_1yjp.replace(
    "OH  TYR A   7       3.766   0.589  10.291  1.00",
    "OH  TYR A   7       3.766   0.589  10.291  0.80")
  prog = run_prog("tst_r2n_output_partial", pdb_str,
    neutron_option="hd_and_h")
  h = iotbx.pdb.input(prog.params.output.filename).construct_hierarchy()
  occ = dict((a.name.strip(), a.occ) for a in h.atoms()
             if a.parent().resname == "TYR" and a.name.strip() in ("HH", "DH"))
  assert approx_equal(occ["HH"], 0.4) and approx_equal(occ["DH"], 0.4), occ
  sites = prog._neutron_result.partial_parent_sites
  assert [(x.site[5:9].strip(), x.parent_occ) for x in sites] == [
    ("HH", 0.8)], sites

def test_output_check_fires():
  # a provenance map missing an entry is an internal error
  from mmtbx.hydrogens import hd_exchange
  original = hd_exchange.apply_neutron_option
  def broken(*args, **kwargs):
    r = original(*args, **kwargs)
    r.provenance.pop(sorted(r.provenance)[0])
    return r
  hd_exchange.apply_neutron_option = broken
  try:
    run_reduce2("tst_r2n_output_broken", pdb_1yjp,
      ["neutron_option=hd_and_h"])
  except Sorry as e:
    assert "Internal error: H/D output check failed" in str(e), str(e)
    assert "no provenance" in str(e)
  else:
    raise AssertionError("output check did not fire")
  finally:
    hd_exchange.apply_neutron_option = original

def test_dipeptide_modes():
  # checked against phenix.ready_set on 2026-10-17; differences explained in
  # the neutron_option plan
  for option in NEUTRON:
    fn = cached_run(dipeptide, neutron_option=option).params.output.filename
    ours, expected = hd_labels(fn), hd_table(dipeptide_hd[option])
    assert ours == expected, (option, sorted(ours ^ expected))

def test_his_split():
  # His in hd_and_h: HD1 and HE2 split into A/B with DD1/DE2. approach=place
  # keeps both ring N-H (no Optimizer, so no His flip or protonation choice).
  # checked against phenix.ready_set on 2026-10-17; differences explained in
  # the neutron_option plan
  fn, log, description = run_reduce2("tst_r2n_output_his", his_no_h,
    ["neutron_option=hd_and_h", "approach=place"])
  ours, expected = hd_labels(fn), hd_table(his_hd)
  assert ours == expected, sorted(ours ^ expected)

def test_hc4_all_d():
  # HC4 ligand, all_d: 7 D and no H
  # checked against phenix.ready_set on 2026-10-17; differences explained in
  # the neutron_option plan
  fn, log, description = run_reduce2("tst_r2n_output_hc4_d", hc4_no_h,
    ["neutron_option=all_d"])
  ours, expected = hd_labels(fn), hd_table(hc4_all_d)
  assert ours == expected, sorted(ours ^ expected)

def test_hc4_hd_and_h():
  # 6: ligand OH: HO4' split, the six C-H stay H
  fn, log, description = run_reduce2("tst_r2n_output_hc4_hd",
    hc4_no_h, ["neutron_option=hd_and_h"])
  labels = atom_label_set(fn)
  hd = sorted((x[1], x[2]) for x in labels if x[3] in ("H", "D"))
  assert ("HO4'", "A") in hd and ("DO4'", "B") in hd, hd
  assert len([x for x in hd if x[1] == ""]) == 6, hd

def test_glycan_link():
  # 7: NAG-ASN (fixture of tst_add_hydrogen_16, 1a7s ASN 114 - NAG 401):
  # ASN ND2 H split, ring OH split, no H on the O1 the link replaces
  from mmtbx.hydrogens.tst_add_hydrogen_16 import pdb_str_019_control
  fn, log, description = run_reduce2("tst_r2n_output_nag",
    pdb_str_019_control, ["neutron_option=hd_and_h"])
  h = iotbx.pdb.input(fn).construct_hierarchy()
  def hd(resname):
    return sorted((a.name.strip(), a.parent().altloc) for a in h.atoms()
      if a.parent().resname == resname and a.element.strip() in ("H", "D"))
  asn, nag = hd("ASN"), hd("NAG")
  assert ("HD21", "A") in asn and ("DD21", "B") in asn, asn
  assert not [x for x in asn if x[0][1:] == "D22"], asn
  for name in ("O3", "O4", "O6"):
    assert ("H" + name, "A") in nag and ("D" + name, "B") in nag, nag
  assert not [x for x in nag if x[0][1:] == "O1"], nag

def test_name_not_pairable_end_to_end():
  # 9i: 1HO2 stays H (never D) and is listed; the other OH is split
  # (hd_and_h) or D (all_d); the written model reinterprets with restraints
  for option, ho3 in (("hd_and_h", "split"), ("all_d", "D")):
    prog = run_prog("tst_r2n_output_zq3_" + option,
      heavy_pdb("03", zq3_atoms), zq3_cif, neutron_option=option)
    fn = prog.params.output.filename
    h = iotbx.pdb.input(fn).construct_hierarchy()
    hd = sorted((a.name.strip(), a.parent().altloc, a.element.strip())
                for a in h.atoms() if a.element.strip() in ("H", "D"))
    assert ("1HO2", "", "H") in hd, (option, hd)
    assert not [x for x in hd if x[0] == "1DO2"], hd
    if ho3 == "split":
      assert ("HO3", "A", "H") in hd and ("DO3", "B", "D") in hd, hd
    else:
      assert ("DO3", "", "D") in hd, hd
    with open(prog.params.output.description_file_name) as f:
      description = f.read()
    line = [l for l in description.splitlines() if "1HO2  03" in l]
    assert len(line) == 1 and "name not pairable by phenix.refine" in line[0]
    split = sorted(tuple(sorted(h.atoms()[i[0]].name.strip() for i in x))
                   for x in h.exchangeable_hd_selections())
    assert split == ([("DO3", "HO3")] if ho3 == "split" else []), split
    # reads back and reinterprets with restraints (and the CIF)
    from iotbx.cif import reader as cif_reader
    model = mmtbx.model.manager(model_input=iotbx.pdb.input(fn),
      restraint_objects=[("zq3.cif", cif_reader(
        input_string=zq3_cif).model())], log=null_out())
    model.process(make_restraints=True)
    from mmtbx.hydrogens import hd_exchange
    table = hd_exchange.h_parent_table(model)
    assert set(e.status for e in table.values()) == set(["ok"])

def test_section_in_log():
  # the H/D section also goes to the log, so it is there without files
  prefix = "tst_r2n_output_zq3_log"
  for ext in ("_out.pdb", ".txt"):
    if os.path.exists(prefix + ext): os.remove(prefix + ext)
  prog = run_prog(prefix, heavy_pdb("03", zq3_atoms), zq3_cif,
    neutron_option="all_d", output__write_files=False)
  assert not os.path.exists(prefix + "_out.pdb")
  log = prog.log_text
  assert log.count("H/D model (neutron_option=all_d):") == 1
  line = [l for l in log.splitlines() if "1HO2  03" in l]
  assert len(line) == 1 and "kept as H" in line[0], line
  # the counts are printed once
  assert log.count("Split into H/D pairs") == 1

def test_mmcif_output():
  # 11: .cif output reads back with the PDB output's H/D labels and
  # occupancies and the same exchangeable pairs
  results = {}
  for ext in ("pdb", "cif"):
    if ext == "pdb":
      fn = cached_run(pdb_1yjp,
        neutron_option="hd_and_h").params.output.filename
    else:
      fn = "tst_r2n_output_mmcif_out." + ext
      run_reduce2("tst_r2n_output_mmcif_" + ext, pdb_1yjp,
        ["neutron_option=hd_and_h", "output.filename=" + fn])
    h = iotbx.pdb.input(fn).construct_hierarchy()
    labels = set((a.parent().parent().resseq.strip(), a.parent().resname,
      a.name.strip(), a.parent().altloc, a.element.strip(), round(a.occ, 2))
      for a in h.atoms())
    atoms = h.atoms()
    pairs = set(frozenset(atoms[i[0]].parent().parent().resseq.strip() + ":" +
      atoms[i[0]].name.strip() for i in x)
      for x in h.exchangeable_hd_selections())
    results[ext] = (labels, pairs)
  assert results["cif"][0] == results["pdb"][0], sorted(
    results["cif"][0] ^ results["pdb"][0])
  assert results["cif"][1] == results["pdb"][1]
  assert len(results["pdb"][1]) == 20

def test_description_file():
  # every category, by written name and with its reason: conformer site,
  # merged copies, partial parent, name not pairable, water renamed
  from mmtbx.hydrogens.tst_hd_exchange import pdb_pep_og_ab
  pdb_str = pdb_pep_og_ab.replace(
    "NZ  LYS A   3      19.531  -1.812  18.059  1.00",
    "NZ  LYS A   3      19.531  -1.812  18.059  0.80").replace("END\n", "")
  assert "0.80" in pdb_str
  pdb_str += heavy_pdb("03", zq3_atoms, offset=(5.0, 15.0, 5.0),
    cryst1=False) + "END\n"
  prog = run_prog("tst_r2n_output_desc", pdb_str, zq3_cif,
    neutron_option="hd_and_h", exclude_water=False)
  with open(prog.params.output.description_file_name) as f:
    description = f.read()
  section = description[description.index("H/D model (neutron_option="
    "hd_and_h)"):]
  def entry(label):
    found = [l for l in section.splitlines() if l.strip().startswith(label)]
    assert len(found) == 1, (label, found)
    return found[0]
  assert "Neutron distances used (neutron_option=hd_and_h)" in section
  assert "Waters renamed DOD: 1" in section
  assert entry("DG ASER A   2").endswith("conformer")
  assert entry("DG BSER A   2").endswith("conformer")
  assert "merged conformer copies" in entry("HA  SER A   2")
  for name in ("HZ1", "HZ2", "HZ3"):
    assert "parent occupancy 0.80" in entry("%sALYS A   3" % name)
  assert ("Split sites on a partially occupied parent (3):\n"
    "    New pairs get half the parent occupancy each, restored pairs keep "
    "their input values; phenix.refine currently refines H + D = 1.") in \
    section, section
  assert "name not pairable by phenix.refine" in entry("1HO2  03 A   1")
  # none: no H/D section; 9b: explicit none keeps today's distances
  fn, log, description = run_reduce2("tst_r2n_output_desc_none",
    pdb_1yjp, ["neutron_option=none"])
  assert "H/D model" not in description
  check_bond_lengths(fn, use_neutron_distances=False)
  check_no_d(fn)

if __name__ == "__main__":
  test_output_names()
  test_validate_sorrys()
  test_rebuild_model()
  test_n_terminal_nh3()
  test_conformers()
  test_output_symmetry()
  test_output_model()
  test_occupancy_groups()
  test_partial_parent()
  test_output_check_fires()
  test_dipeptide_modes()
  test_his_split()
  test_hc4_all_d()
  test_hc4_hd_and_h()
  test_glycan_link()
  test_name_not_pairable_end_to_end()
  test_section_in_log()
  test_mmcif_output()
  test_description_file()
  print("OK")
