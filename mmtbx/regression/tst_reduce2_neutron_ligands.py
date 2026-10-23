"""reduce2 neutron tests: ligands: user CIF, missing restraints, name mismatch,
CCD-only"""
from __future__ import absolute_import, division, print_function
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out, Sorry
from mmtbx.hydrogens import reduce_hydrogen
from mmtbx.regression import reduce2_neutron_utils
from mmtbx.regression.reduce2_neutron_utils import (
  FLOATING, MISSING_RESTRAINTS, UNCHANGED, UNKNOWN_ATOMS, UNKNOWN_ATOMS_DO,
  bgm_atoms, bgm_xh, cached_run, check_bond_lengths, heavy_pdb,
  ligand_split_records, pdb_1yjp, protein_ligand_src, run_prog, zq4_atoms,
  zq4_cif)

reduce2_neutron_utils.set_prefix("tst_r2n_ligands")

def test_user_cif_ligand():
  # 9c (add): invented ligand with OH and NH2 from a CIF written here: the
  # polar H are split, no missing restraints; without the CIF the run fails
  prog = cached_run(heavy_pdb("04", zq4_atoms), zq4_cif,
    neutron_option="hd_and_h")
  h = iotbx.pdb.input(prog.params.output.filename).construct_hierarchy()
  hd = sorted((a.name.strip(), a.parent().altloc) for a in h.atoms()
              if a.element.strip() in ("H", "D"))
  for name in ("HN1", "HO1"):
    assert (name, "A") in hd and ("D" + name[1:], "B") in hd, hd
  assert len([x for x in hd if x[1] == ""]) == 7, hd
  with open(prog.params.output.description_file_name) as f:
    description = f.read()
  for text in ("Restraints were not found", "Insufficient restraints",
               "not placed", "Could not find atom info"):
    assert text not in description and text not in prog.log_text, text
  try:
    run_prog("tst_r2n_ligands_zq4_nocif", heavy_pdb("04", zq4_atoms),
      neutron_option="hd_and_h")
  except Sorry as e:
    assert str(e) == MISSING_RESTRAINTS + " 04", str(e)
  else:
    raise AssertionError("04 ran without its CIF")

def test_user_cif_ligand_optimize():
  # 9c (optimize): the 04 hd_and_h output (D input) with its CIF, optimize
  # hd_and_h: runs, the polar H are split, no missing restraints
  src = cached_run(heavy_pdb("04", zq4_atoms), zq4_cif,
    neutron_option="hd_and_h").text
  assert " DN1B 04" in src and " DO1B 04" in src
  prog = run_prog("tst_r2n_ligands_zq4o", src, zq4_cif,
    approach="optimize", neutron_option="hd_and_h")
  h = iotbx.pdb.input(prog.params.output.filename).construct_hierarchy()
  hd = sorted((a.name.strip(), a.parent().altloc) for a in h.atoms()
              if a.element.strip() in ("H", "D"))
  for name in ("HN1", "HO1"):
    assert (name, "A") in hd and ("D" + name[1:], "B") in hd, hd
  assert len([x for x in hd if x[1] == ""]) == 7, hd
  with open(prog.params.output.description_file_name) as f:
    description = f.read()
  for text in ("Restraints were not found", "Insufficient restraints",
               "without a covalent partner", "Could not find atom info"):
    assert text not in description and text not in prog.log_text, text
  assert "Values restored from the input: 2" in description
  # without the CIF: the missing-restraints Sorry
  try:
    run_prog("tst_r2n_ligands_zq4o_nocif", src, approach="optimize",
      neutron_option="hd_and_h")
  except Sorry as e:
    assert str(e) == MISSING_RESTRAINTS + " 04", str(e)
  else:
    raise AssertionError("04 ran without its CIF")

def test_missing_restraints_optimize():
  # ligand 04 without its CIF stops optimize with the missing-restraints
  # Sorry, H input and D input: checked once at the start of the run,
  # before the input check
  import traceback
  for option in ("none", "hd_and_h"):
    try:
      run_prog("tst_r2n_ligands_mr_" + option,
        protein_ligand_src(option), approach="optimize",
        neutron_option=option)
    except Sorry as e:
      assert str(e) == MISSING_RESTRAINTS + " 04", str(e)
      functions = [f[2] for f in traceback.extract_tb(e.__traceback__)]
      assert "_PrepareRestraints" in functions and \
        "_NormalizeInput" not in functions, functions
    else:
      raise AssertionError("04 ran without its CIF (%s)" % option)

def test_missing_restraints_ignored():
  # ignore_missing_restraints=True: ligand 04 (no CIF) is written exactly as
  # in the input (D stay D in all_h) and listed; the protein is converted
  # as usual; nothing of the ligand is floating
  for option, source in (("none", "none"), ("hd_and_h", "hd_and_h"),
                         ("all_h", "hd_and_h")):
    src = protein_ligand_src(source)
    prog = run_prog("tst_r2n_ligands_mri_" + option, src,
      approach="optimize", neutron_option=option,
      ignore_missing_restraints=True)
    with open(prog.params.output.filename) as f: out = f.read()
    lig_in = ligand_split_records(src, True)
    assert len(lig_in) == (16 if source == "hd_and_h" else 14), len(lig_in)
    assert ligand_split_records(out, True) == lig_in, option
    with open(prog.params.output.description_file_name) as f:
      description = f.read()
    for text in (description, prog.log_text):
      assert UNCHANGED + "\n04 B 1\n" in text, option
    assert FLOATING not in description
    # the protein: H/D labels as from the same mode without the ligand
    def protein(pdb_str):
      return sorted((r[0], r[1])
                    for r in ligand_split_records(pdb_str, False))
    if option == "all_h":
      assert "D" not in [r[1] for r in ligand_split_records(out, False)]
      assert sorted(r[1] for r in lig_in).count("D") == 2
    else:
      assert protein(out) == protein(src), option

def test_ccd_only_ligand_optimize():
  # BGM has no monomer library entry: placement builds its restraints from
  # the CCD; optimize and the input check get the same objects, so its H
  # are kept (before: removed as floating, then the Optimizer crashed when
  # BGM was the only residue)
  prog = run_prog("tst_r2n_ligands_bgm", heavy_pdb("BGM", bgm_atoms()),
    neutron_option="none")
  with open(prog.params.output.filename) as f: src = f.read()
  def hd(pdb_str):
    h = iotbx.pdb.input(source_info=None, lines=pdb_str).construct_hierarchy()
    return sorted((a.name.strip(), a.parent().altloc, a.element.strip())
                  for a in h.atoms() if a.element.strip() in ("H", "D"))
  placed = hd(src)
  assert len(placed) == 13, placed
  polar = ["H1", "H21", "H22", "HO3'", "HOP2", "HOP3"]
  for option in ("none", "hd_and_h"):
    prog = run_prog("tst_r2n_ligands_bgm_" + option, src,
      approach="optimize", neutron_option=option)
    with open(prog.params.output.filename) as f: out = f.read()
    with open(prog.params.output.description_file_name) as f:
      description = f.read()
    assert FLOATING not in description and UNCHANGED not in description
    if option == "none":
      assert hd(out) == placed
    else:
      expected = sorted([x for x in placed if x[0] not in polar] +
        [(n, "A", "H") for n in polar] +
        [("D" + n[1:], "B", "D") for n in polar])
      assert hd(out) == expected, hd(out)
    # the OH (O-P-OH) H are now Optimizer movers
    movers = [l for l in description.splitlines()
              if l.strip().startswith("Added MoverSingleHydrogenRotator")
              and "BGM" in l]
    assert len(movers) == 3, movers
  # input check (D input) gets the CCD restraints too
  with open(prog.params.output.filename) as f: d_src = f.read()
  prog = run_prog("tst_r2n_ligands_bgm_d", d_src, approach="optimize",
    neutron_option="hd_and_h")
  with open(prog.params.output.filename) as f: out = f.read()
  assert hd(out) == hd(d_src)
  # the CCD objects are dropped again
  assert [n for n, o in (prog.model.get_restraint_objects() or [])
          if n.startswith("auto_")] == []

def test_missing_restraints_keep_existing_h():
  # add/place/optimize keep_existing_H=True ignore_missing_restraints=True:
  # ligand 04 (no CIF) is written exactly as in the input, D and their 0.5
  # occupancy and B too, and listed (with flipkins: test_flipkins_combined).
  # add keep_existing_H=False: no H
  for approach in ("add", "place", "optimize"):
    for option in ("none", "hd_and_h"):
      prefix = "tst_r2n_ligands_mrk_%s_%s" % (approach, option)
      src = protein_ligand_src(option)
      prog = run_prog(prefix, src, approach=approach, neutron_option=option,
        keep_existing_H=True, ignore_missing_restraints=True)
      with open(prog.params.output.filename) as f: out = f.read()
      lig_in = ligand_split_records(src, True)
      assert len(lig_in) == (16 if option == "hd_and_h" else 14), len(lig_in)
      assert ligand_split_records(out, True) == lig_in, (approach, option)
      if option == "hd_and_h":
        d_in = [r for r in lig_in if r[1] == "D"]
        assert len(d_in) == 2 and [r[3] for r in d_in] == [0.5, 0.5], d_in
        assert [r for r in ligand_split_records(out, True)
                if r[1] == "D"] == d_in
      with open(prog.params.output.description_file_name) as f:
        description = f.read()
      for text in (description, prog.log_text):
        assert UNCHANGED + "\n04 B 1\n" in text, (approach, option)
  prog = run_prog("tst_r2n_ligands_mrk_nokeep", protein_ligand_src("none"),
    approach="add", ignore_missing_restraints=True)
  with open(prog.params.output.filename) as f: out = f.read()
  assert [r[1] for r in ligand_split_records(out, True)] == \
    ["C", "C", "C", "N", "O"]
  assert UNCHANGED not in prog.log_text

def test_name_mismatch_reported():
  # ligand 04's atoms named ZQ4 (a GeoStd entry without N1), no CIF: optimize
  # lists the residue, its atom and the H/D removed there under that cause
  # (not as floating); add has the report in the description file too, and
  # optimize on its output (no H/D affected) without the H/D clause; a run
  # without such a residue has no such text
  for option in ("none", "hd_and_h"):
    src = protein_ligand_src(option).replace(" 04 B", "ZQ4 B")
    prog = run_prog("tst_r2n_ligands_nm_opt_" + option, src,
      approach="optimize", neutron_option=option)
    with open(prog.params.output.description_file_name) as f:
      description = f.read()
    removed = ["H11 ZQ4 B   1", "H12 ZQ4 B   1", "H21 ZQ4 B   1",
      "H31 ZQ4 B   1", "H32 ZQ4 B   1", "H33 ZQ4 B   1"]
    if option == "none":
      removed += ["HN1 ZQ4 B   1", "HO1 ZQ4 B   1"]
    else:
      removed += ["HN1AZQ4 B   1", "HO1AZQ4 B   1", "DN1BZQ4 B   1",
                  "DO1BZQ4 B   1"]
    block = (UNKNOWN_ATOMS + "1 residue; %d H/D on them were removed:\n"
      " ZQ4 B 1: N1\n" % len(removed) +
      "".join("   %s\n" % l for l in removed) + UNKNOWN_ATOMS_DO + "ZQ4.\n")
    for text in (description, prog.log_text):
      assert block in text, (option, text)
      # not also as floating (optimize, input check)
      assert "without a covalent partner" not in text, option
  prog = run_prog("tst_r2n_ligands_nm_add",
    protein_ligand_src("none").replace(" 04 B", "ZQ4 B"), approach="add")
  with open(prog.params.output.description_file_name) as f:
    description = f.read()
  block = (UNKNOWN_ATOMS + "1 residue; 4 H atoms on them were not placed:\n"
    " ZQ4 B 1: N1\n" + UNKNOWN_ATOMS_DO + "ZQ4.\n")
  for text in (description, prog.log_text):
    assert block in text, text
  # optimize on that output: no H/D affected, no H/D clause
  with open(prog.params.output.filename) as f: added = f.read()
  prog = run_prog("tst_r2n_ligands_nm_opt_added", added,
    approach="optimize")
  with open(prog.params.output.description_file_name) as f:
    description = f.read()
  block = (UNKNOWN_ATOMS + "1 residue:\n ZQ4 B 1: N1\n" + UNKNOWN_ATOMS_DO +
    "ZQ4.\n")
  for text in (description, prog.log_text):
    assert block in text, text
  for approach in ("add", "optimize"):
    if approach == "add":
      prog = cached_run(pdb_1yjp)
    else:
      prog = run_prog("tst_r2n_ligands_nm_no_" + approach,
        protein_ligand_src("none"), zq4_cif, approach=approach)
    with open(prog.params.output.description_file_name) as f:
      description = f.read()
    for text in (description, prog.log_text):
      assert "Atom names not in" not in text and "Rename these" not in text

def test_ccd_only_neutron_targets():
  # CCD-only ligand (BGM): with neutron distances its X-H targets are the
  # library's (as the protein's in the same model: C-H, N-H, O-H; aromatic
  # C-H for TPT); X-ray targets stay 0.9 x the CCD length
  bgm = heavy_pdb("BGM", bgm_atoms(), offset=(15.0, 10.0, 15.0),
    cryst1=False).replace(" A   1 ", " B   1 ")
  h_model = cached_run(pdb_1yjp.replace("END\n", "") + bgm + "END\n").text
  for neutron in (False, True):
    model = mmtbx.model.manager(model_input=iotbx.pdb.input(
      source_info=None, lines=h_model), log=null_out())
    reduce_hydrogen.add_ccd_restraint_objects(model)
    model.process(make_restraints=True, plain_pairs_radius=None,
      pdb_interpretation_params=
        reduce_hydrogen.get_reduce_pdb_interpretation_params(neutron))
    atoms = model.get_hierarchy().atoms()
    simple, asu = model.get_restraints_manager().geometry.\
      get_all_bond_proxies(sites_cart=model.get_sites_cart())
    targets = {}
    for proxy in simple:
      a, b = [atoms[i] for i in proxy.i_seqs]
      if (a.element.strip() == "H") == (b.element.strip() == "H"): continue
      hh, x = (a, b) if a.element.strip() == "H" else (b, a)
      kind = "BGM" if hh.parent().resname == "BGM" else "protein"
      if kind == "protein" and hh.name.strip() not in (
          "HA", "H", "HH", "HD1"): continue
      key = (kind, x.element.strip(), x.name.strip() == "CD1")
      targets.setdefault(key, set()).add(round(proxy.distance_ideal, 3))
    bgm_targets = dict((k[1], v) for k, v in targets.items()
                       if k[0] == "BGM")
    assert sorted(bgm_targets) == ["C", "N", "O"], bgm_targets
    if neutron:
      for element in ("C", "N", "O"):
        assert bgm_targets[element] == \
          targets[("protein", element, False)], (element, targets)
      aromatic = targets[("protein", "C", True)]
    else:
      # 0.9 x the CCD ideal lengths (C-H 1.089-1.090, N-H/O-H ~0.97)
      assert bgm_targets["C"] == set([0.98, 0.981]), bgm_targets
      assert bgm_targets["N"] == set([0.873]), bgm_targets
      assert bgm_targets["O"] == set([0.87, 0.871]), bgm_targets
  # aromatic C-H of a CCD-only ligand (TPT): the library's aromatic value
  md = reduce_hydrogen.get_h_restraints("TPT", strict=False)
  values = set(round(float(b.value_dist_neutron), 3) for b in md.bond_list
               if b.value_dist_neutron is not None and "C" in (
               b.atom_id_1[:1], b.atom_id_2[:1]))
  assert aromatic <= values, (aromatic, values)

def test_ccd_only_neutron_lengths():
  # 1yjp plus BGM: hd_and_h (add, and optimize on the add output) writes
  # BGM's X-H/X-D within 0.02 A of the library's neutron targets, as the
  # protein's; add none keeps 0.9 x the CCD length (X-ray unchanged)
  library = reduce_hydrogen.get_library_neutron_xh()
  bgm = heavy_pdb("BGM", bgm_atoms(), offset=(15.0, 10.0, 15.0),
    cryst1=False).replace(" A   1 ", " B   1 ")
  model = pdb_1yjp.replace("END\n", "") + bgm + "END\n"
  prog = cached_run(model)
  for x, e, d, d_ideal in bgm_xh(prog.params.output.filename):
    assert abs(d - 0.9 * d_ideal) < 0.005, (x, e, d, d_ideal)
  prog = run_prog("tst_r2n_ligands_ccdl_add", model,
    neutron_option="hd_and_h")
  fn = prog.params.output.filename
  with open(fn) as f: src = f.read()
  prog = run_prog("tst_r2n_ligands_ccdl_opt", src, approach="optimize",
    neutron_option="hd_and_h")
  for fn in (fn, prog.params.output.filename):
    found = bgm_xh(fn)
    assert len(found) == 19, len(found)
    for x, e, d, d_ideal in found:
      assert abs(d - library[x]) < 0.02, (fn, x, e, d)
    # the protein, by its own restraints (interpreted without BGM, which
    # has no restraints outside reduce2)
    with open(fn) as f:
      protein = [l for l in f if l[17:20] != "BGM"]
    with open(fn[:-4] + "_protein.pdb", "w") as f:
      f.write("".join(protein))
    check_bond_lengths(fn[:-4] + "_protein.pdb", use_neutron_distances=True)

def test_ccd_only_heavy_atom_bonds():
  # CCD placeholder (BGM): heavy-atom bonds at the CCD ideal length, X-H at
  # 0.9 x; the grm that place_hydrogens leaves on the model has the same
  from mmtbx.chemical_components import get_cif_dictionary
  from scitbx import matrix
  cc = get_cif_dictionary("BGM")
  xyz = dict((a.atom_id, matrix.col((float(a.pdbx_model_Cartn_x_ideal),
    float(a.pdbx_model_Cartn_y_ideal), float(a.pdbx_model_Cartn_z_ideal))))
    for a in cc["_chem_comp_atom"])
  is_h = dict((a.atom_id, a.type_symbol in ("H", "D"))
              for a in cc["_chem_comp_atom"])
  md = reduce_hydrogen.get_h_restraints("BGM", strict=False)
  n_heavy = n_xh = 0
  for b in md.bond_list:
    d = (xyz[b.atom_id_1] - xyz[b.atom_id_2]).length()
    if is_h[b.atom_id_1] or is_h[b.atom_id_2]:
      assert abs(float(b.value_dist) - 0.9 * d) < 0.001, (b.atom_id_1,
        b.atom_id_2, b.value_dist, d)
      n_xh += 1
    else:
      assert abs(float(b.value_dist) - d) < 0.001, (b.atom_id_1,
        b.atom_id_2, b.value_dist, d)
      n_heavy += 1
  assert (n_heavy, n_xh) == (26, 13), (n_heavy, n_xh)
  bgm = heavy_pdb("BGM", bgm_atoms(), offset=(15.0, 10.0, 15.0),
    cryst1=False).replace(" A   1 ", " B   1 ")
  model = mmtbx.model.manager(model_input=iotbx.pdb.input(source_info=None,
    lines=pdb_1yjp.replace("END\n", "") + bgm + "END\n"), log=null_out())
  ph = reduce_hydrogen.place_hydrogens(model=model)
  ph.run()
  model = ph.get_model()
  atoms = model.get_hierarchy().atoms()
  simple, asu = model.get_restraints_manager().geometry.\
    get_all_bond_proxies(sites_cart=model.get_sites_cart())
  n = 0
  for proxy in simple:
    a, b = [atoms[i] for i in proxy.i_seqs]
    if a.parent().resname != "BGM" or b.parent().resname != "BGM": continue
    if "H" in (a.element.strip(), b.element.strip()): continue
    d = (xyz[a.name.strip()] - xyz[b.name.strip()]).length()
    assert abs(proxy.distance_ideal - d) < 0.001, (a.name, b.name,
      proxy.distance_ideal, d)
    n += 1
  assert n == n_heavy, (n, n_heavy)

if __name__ == "__main__":
  test_user_cif_ligand()
  test_user_cif_ligand_optimize()
  test_missing_restraints_optimize()
  test_missing_restraints_ignored()
  test_ccd_only_ligand_optimize()
  test_missing_restraints_keep_existing_h()
  test_name_mismatch_reported()
  test_ccd_only_neutron_targets()
  test_ccd_only_neutron_lengths()
  test_ccd_only_heavy_atom_bonds()
  print("OK")
