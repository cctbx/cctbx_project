"""reduce2 neutron tests: deposited H/D layouts repaired at the input check
(B16), approach=optimize and approach=add keep_existing_H=True:
(a) a blank H/D on a parent split into conformer copies, (b) names permuted
between the H set and the D set, (c) one isotope per conformer on the altloc
copies of one parent, (d) an H/D pair in heavy-atom conformer groups on a
blank parent, (e) H and D apart on a one-H parent. Each repair case passes
with its values kept at their sites (by isotope and conformer), and a rerun
on the output is stable; each control stops with the reason it was not
repaired. The fixtures are cuts of 7T5E, 6TAE, 7KCU and 6KK8
(tst_hd_exchange)."""
from __future__ import absolute_import, division, print_function
import os
import iotbx.pdb
from libtbx import group_args
from libtbx.utils import Sorry
from mmtbx.hydrogens import hd_exchange
from mmtbx.hydrogens import tst_hd_exchange as t
from mmtbx.regression import reduce2_neutron_utils
from mmtbx.regression.reduce2_neutron_utils import run_prog

reduce2_neutron_utils.set_prefix("tst_r2n_repair")

# tier E of the deposited-model audit
SETTINGS = dict(neutron_option="hd_and_h", use_neutron_distances=True,
                add_flip_movers=True)
APPROACHES = (("opt", dict(approach="optimize")),
              ("keep", dict(approach="add", keep_existing_H=True)))

def run(prefix, pdb_str, cif_str=None, **params):
  p = dict(SETTINGS)
  p.update(params)
  return run_prog(prefix, pdb_str, cif_str, **p)

def hierarchy(text):
  return iotbx.pdb.input(source_info=None, lines=text).construct_hierarchy()

def is_hd(a):
  return a.element.strip() in ("H", "D")

def site_of(a):
  """(chain, resseq, resname, name as H, conformer); the conformer of a blank
  or H/D-only (layout) atom_group is '' (a pair site holds both isotopes)"""
  ag = a.parent()
  rg = ag.parent()
  conformer = "" if hd_exchange._hd_only(ag) else ag.altloc
  return (rg.parent().id, rg.resseq.strip(), ag.resname,
          hd_exchange.normalized_name(a).strip(), conformer)

def value(a):
  return (round(a.occ, 2), round(a.b, 2))

def written(h):
  """{site_of(): sorted [(isotope, occ, b)]} of every H/D"""
  result = {}
  for a in h.atoms():
    if is_hd(a):
      result.setdefault(site_of(a), []).append((a.element.strip(),) + value(a))
  return dict((k, sorted(v)) for k, v in result.items())

def expected_input(pdb_str, cif_str=None):
  """The input with the planned repairs applied, matched as reduce2 matches
  it (accepted pairs included): the pairs and singles on the parents of the
  repaired atoms (the sites checked; other sites follow reduce2's usual
  rules, e.g. merged conformer copies). pairs: {site: {isotope: (occ, b)}};
  singles: {site: (occ, b)}; atoms: those H/D, for the check by position."""
  def make(text):
    return t.make_model(text) if cif_str is None else \
      t.make_model_cif(text, cif_str)
  model = make(pdb_str)
  plan, c, m = t.plan_of(model)
  h = model.get_hierarchy().deep_copy()
  origin = hd_exchange.apply_repairs(h, plan)
  m2 = make(h.as_pdb_string(crystal_symmetry=model.crystal_symmetry()))
  accepted = set((x.h_key, x.d_key) for x in plan.accept + plan.conformer)
  table = hd_exchange.h_parent_table(m2)
  c2 = hd_exchange.find_hd_candidates(m2.get_hierarchy(), table, model=m2)
  match = hd_exchange.collapse_check(c2, accepted=accepted or None)
  t.check_clean(match)
  atoms = m2.get_hierarchy().atoms()
  keys = hd_exchange.atom_keys(m2.get_hierarchy())
  repaired = set(origin) | accepted_keys(accepted)
  parents = set(table[k].parent_id for k in keys if k in repaired)
  assert parents
  on_parent = lambda i: table[keys[i]].status == "ok" and \
    table[keys[i]].parent_id in parents
  result = group_args(pairs={}, singles={}, atoms=[])
  for p in match.pairs:
    if not on_parent(p.i_h): continue
    h_atom, d_atom = atoms[p.i_h], atoms[p.i_d]
    site = site_of(h_atom)[:4] + ("",)
    assert site not in result.pairs, site
    result.pairs[site] = {"H": value(h_atom), "D": value(d_atom)}
    result.atoms.extend([h_atom, d_atom])
  for s in c2.singles:
    if not on_parent(s.i_seq): continue
    a = atoms[s.i_seq]
    assert site_of(a) not in result.singles, site_of(a)
    result.singles[site_of(a)] = value(a)
    result.atoms.append(a)
  result.hierarchy = m2.get_hierarchy()   # the atoms keep their hierarchy
  return result

def accepted_keys(accepted):
  return set(k for pair in accepted for k in pair)

def check_kept(expected, out_text, label, keep):
  """Each pair site holds exactly its H and its D with their own values, in
  both approaches; each single is in its conformer, with its values unless
  keep (a supplied single gets placement's occupancy and B there, as any
  supplied H). No H/D at an expected position under another name."""
  out = hierarchy(out_text)
  have = written(out)
  for site, v in sorted(expected.pairs.items()):
    want = sorted((iso,) + v[iso] for iso in ("H", "D"))
    assert have.get(site) == want, (label, site, want, have.get(site))
  for site, v in sorted(expected.singles.items()):
    got = have.get(site, [])
    assert len(got) == 1, (label, site, got)
    if not keep:
      assert got[0][1:] == v, (label, site, v, got)
  out_atoms = [a for a in out.atoms() if is_hd(a)]
  for a in expected.atoms:
    for b in out_atoms:
      if b.element.strip() != a.element.strip(): continue
      if b.distance(a) > 0.05: continue
      assert site_of(b)[:4] == site_of(a)[:4], (label, a.id_str(), b.id_str())
  return out

def check_repair(name, pdb_str, n_repairs, title, cif_str=None):
  expected = expected_input(pdb_str, cif_str)
  outputs = {}
  for tag, params in APPROACHES:
    prefix = "tst_r2n_repair_%s_%s" % (name, tag)
    prog = run(prefix, pdb_str, cif_str, **params)
    with open(prog.params.output.description_file_name) as f:
      description = f.read()
    assert "  Input H/D repaired (%d):\n" % n_repairs in description, (
      name, tag, description)
    assert "    %s (%d):\n" % (title, n_repairs) in description, (name, tag)
    with open(prog.params.output.filename) as f:
      out_text = f.read()
    check_kept(expected, out_text, (name, tag), tag == "keep")
    # rerun on the output: nothing to repair, every H/D and value kept
    again = run(prefix + "_rerun", out_text, cif_str, **params)
    with open(again.params.output.description_file_name) as f:
      assert "Input H/D repaired" not in f.read(), (name, tag)
    with open(again.params.output.filename) as f:
      rerun_text = f.read()
    before = written(hierarchy(out_text))
    after = written(hierarchy(rerun_text))
    assert before == after, (name, tag, [(k, before.get(k), after.get(k))
      for k in set(before) | set(after) if before.get(k) != after.get(k)])
    outputs[tag] = (out_text, description)
  return outputs

def check_refused(name, pdb_str, reason):
  for tag, params in APPROACHES:
    try:
      run("tst_r2n_repair_%s_%s" % (name, tag), pdb_str, **params)
    except Sorry as e:
      assert "(sites not repaired)" in str(e), str(e)
      assert "(not repaired: %s" % reason in str(e), (name, tag, str(e))
    else:
      raise AssertionError("no Sorry: %s %s" % (name, tag))

SPLIT = "Blank H/D on a split parent, split per conformer"
RENAME = "H/D names permuted between H and D, renamed by position"
RELABEL = "H/D pair in conformer groups, moved to H/D-only groups"
ACCEPT = "H/D pair on a one-H parent beyond 0.20 A, paired"
CONFORMER = "H and D on altloc copies of one parent, one isotope per conformer"

def test_split():
  # (a): CA copies 0.006 A (7T5E Glu A 140), CA A/C with an H/D-only B
  # (Val B 210), the blank D on N the only blank atom (6TAE Gln A 13),
  # constructed copies 0.08 A apart
  for name, pdb_str in (("a140", t.pdb_a_7t5e_glu_a140),
                        ("b210", t.pdb_a_7t5e_val_b210),
                        ("a13", t.pdb_a_6tae_gln_a13),
                        ("ala008", t.ala_fixture(shift=0.08))):
    check_repair(name, pdb_str, 1, SPLIT)

def test_split_user_cif():
  # review fix 1: with a CIF of CA-HA 1.20 A, a DA at 1.19 A is accepted
  import libtbx.load_env
  fn = libtbx.env.find_in_repositories(
    relative_path=os.path.join("chem_data", "geostd", "a", "data_ALA.cif"),
    test=os.path.isfile)
  if fn is None:
    print("Skipping test_split_user_cif: no chem_data/geostd")
    return
  with open(fn) as f:
    cif = f.read().replace(" ALA  CA  HA   single  0.970  0.020  1.090",
                           " ALA  CA  HA   single  1.200  0.020  1.200")
  pdb_str = t.ala_fixture(d_length=1.19)
  check_refused("ala119", pdb_str, "bond 1.190 A to CA AALA A   1 "
    "(dictionary 0.97/1.09 A)")
  check_repair("ala119_cif", pdb_str, 1, SPLIT, cif_str=cif)

def test_split_refused():
  for name, pdb_str, reason in (
      ("rot14", t.ala_fixture(rot=14),
       "ideal positions on the copies 0.527 A apart (limit 0.10)"),
      ("shift015", t.ala_fixture(shift=0.15),
       "ideal positions on the copies 0.150 A apart (limit 0.10)"),
      ("occ0", t.ala_fixture(occ_a=0, occ_b=0),
       "parent copies with total occupancy 0 (no basis for sharing it)"),
      # step B review fix 3: on the parent's position, no TypeError
      ("d_on_ca", t.ala_d_on_parent(),
       "bond 0.000 A to CA AALA A   1 (dictionary 0.97/1.09 A)")):
    check_refused(name, pdb_str, reason)

def test_rename():
  # (b): amide (7T5E Gln A 14), Lys NH3 shifted (7KCU Lys A 3), with one H
  # completed to the name left free (Lys A 41)
  for name, pdb_str, n in (("a14", t.pdb_b_7t5e_gln_a14, 2),
                           ("a3", t.pdb_b_7kcu_lys_a3, 5),
                           ("a41", t.pdb_b_7kcu_lys_a41, 3)):
    check_repair(name, pdb_str, n, RENAME)

def test_rename_refused():
  check_refused("b125", t.pdb_ctl_7kcu_lys_b125, "1 leftover H, 0 free names")
  check_refused("b4", t.pdb_ctl_7kcu_gln_b4, "no partner within 0.20 A")

def test_relabel():
  # (d): 6KK8 Met A 41, H in conformer A, D in B on the blank N
  check_repair("met41", t.pdb_d_6kk8_met_a41, 1, RELABEL)

def test_accept():
  # (e): 7T5E Arg A 21 HE/DE 0.23 A apart on NE
  check_repair("arg21", t.pdb_e_7t5e_arg_a21, 1, ACCEPT)

def test_accept_floating_chain():
  # step B review fix 1: a chain A holding only a floating H before the
  # fixture; removing it shifts the Arg's chain occurrence, and the accepted
  # pair must follow
  check_repair("arg21_floating", t.floating_chain(t.pdb_e_7t5e_arg_a21), 1,
    ACCEPT)

def test_unrestrained_floating_chain():
  # step B review fix 1, the keys of residues without restraints: ligand 04
  # (no CIF, ignore_missing_restraints=True) in a chain B after a chain B
  # that holds only a floating H. Before, the ligand was converted
  # (optimize) or stopped the run (add: "has parent status floating")
  src = reduce2_neutron_utils.protein_ligand_src("hd_and_h")
  pdb_str = t.floating_chain(src, chain_id="B")
  for tag, params in APPROACHES:
    prog = run("tst_r2n_repair_unrestrained_%s" % tag, pdb_str,
      ignore_missing_restraints=True, **params)
    with open(prog.params.output.filename) as f:
      out = f.read()
    assert reduce2_neutron_utils.ligand_split_records(out, True) == \
      reduce2_neutron_utils.ligand_split_records(src, True), tag
    with open(prog.params.output.description_file_name) as f:
      description = f.read()
    assert reduce2_neutron_utils.UNCHANGED + "\n04 B 1\n" in description
    assert "H   ALA B 900" in description.split(
      "dropped at the input check")[1], tag

def test_accept_refused():
  # (e) on a two-H parent (7KCU Gln A 18 NE2, D out of the amide plane)
  check_refused("gln18", t.pdb_ctl_7kcu_gln_a18,
    "NE2 carries 2 H in the dictionary (one H site needed); D ")

def test_conformer():
  # (c): 7T5E Tyr B 25, HH on OH A at 0.65, DH on OH B at 0.35: two singles,
  # each in its conformer; written as D in both, listed with the conversion
  # and among the conformer sites
  outputs = check_repair("tyr25", t.pdb_c_7t5e_tyr_b25, 1, CONFORMER)
  for tag, (out_text, description) in outputs.items():
    assert "written as D in both conformers (the H converted to D)" in \
      description, (tag, description)
    sites = description.split("Conformer sites written as one isotope")[1]
    for label in ("DH ATYR B  25", "DH BTYR B  25"):
      assert label in sites, (tag, label)
    have = written(hierarchy(out_text))
    for alt, occ in (("A", 0.65), ("B", 0.35)):
      got = have[("B", "25", "TYR", "HH", alt)]
      assert [g[:2] for g in got] == [("D", occ)], (tag, alt, got)

def test_conformer_carbon():
  # (c) on a carbon-bound site (step C review P3): Ala HA in conformer A, DA
  # in B, each at its conformer's occupancy. hd_and_h and all_hd write H in
  # both conformers, and the line says so
  pdb_str = t.ala_conformer_ha()
  for option in ("hd_and_h", "all_hd"):
    for tag, params in APPROACHES:
      prog = run("tst_r2n_repair_ala_c_%s_%s" % (option, tag), pdb_str,
        neutron_option=option, **params)
      with open(prog.params.output.description_file_name) as f:
        description = f.read()
      assert "    %s (1):\n" % CONFORMER in description, (option, tag)
      assert "written as H in both conformers (the D converted to H)" in \
        description, (option, tag, description)
      assert "written as D in both" not in description, (option, tag)
      with open(prog.params.output.filename) as f:
        have = written(hierarchy(f.read()))
      for alt, occ in (("A", 0.6), ("B", 0.4)):
        got = have[("A", "1", "ALA", "HA", alt)]
        assert [g[:2] for g in got] == [("H", occ)], (option, tag, alt, got)

def test_conformer_refused():
  # Ser B 82: H 0.15 on N A, D 0.85 on N B, H + D = 1
  check_refused("ser82", t.pdb_ctl_7t5e_ser_b82, t.SER_B82_REASON)

def test_different_parents_refused():
  # (c) only for each H/D at its conformer's occupancy: HH at 0.50 in
  # conformer A (0.65) is neither that nor an exchange pair
  pdb_str = t.pdb_c_7t5e_tyr_b25.replace("  0.65 24.15           H",
                                         "  0.50 24.15           H")
  check_refused("tyr25_occ", pdb_str, "HH 0.50 on OH ATYR B  25 (conformer "
    "0.65), DH 0.35 on OH BTYR B  25 (conformer 0.35): neither one isotope "
    "per conformer nor an exchange pair")

if __name__ == "__main__":
  test_split()
  test_split_user_cif()
  test_split_refused()
  test_rename()
  test_rename_refused()
  test_relabel()
  test_accept()
  test_accept_floating_chain()
  test_unrestrained_floating_chain()
  test_accept_refused()
  test_conformer()
  test_conformer_carbon()
  test_conformer_refused()
  test_different_parents_refused()
  print("OK")
