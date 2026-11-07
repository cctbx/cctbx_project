"""reduce2 neutron tests: deposited H/D layouts repaired at the input check
(B16), approach=optimize and approach=add keep_existing_H=True:
(a) a blank H/D on a parent split into conformer copies, (b) names permuted
between the H set and the D set, (d) an H/D pair in heavy-atom conformer
groups on a blank parent, (e) H and D apart on a one-H parent. Each repair
case passes with its values kept at their sites, and a rerun on the output
is stable; each control stops with the reason it was not repaired. The
fixtures are cuts of 7T5E, 6TAE, 7KCU and 6KK8 (tst_hd_exchange)."""
from __future__ import absolute_import, division, print_function
import os
import iotbx.pdb
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

def hd_atoms(h):
  return [a for a in h.atoms() if a.element.strip() in ("H", "D")]

def residue_name(a):
  """(chain, resseq, resname, name as H): the written isotope follows the
  output rules (hd_and_h: D on carbon as H, an exchangeable single in an
  altloc group as D), the values must not"""
  rg = a.parent().parent()
  return (rg.parent().id, rg.resseq.strip(), a.parent().resname,
          hd_exchange.normalized_name(a).strip())

def values(atoms):
  """{residue_name(): sorted [(occ, b)]}"""
  result = {}
  for a in atoms:
    result.setdefault(residue_name(a), []).append((round(a.occ, 2),
      round(a.b, 1)))
  return dict((k, sorted(v)) for k, v in result.items())

def expected_input(pdb_str, cif_str=None):
  """the input with the planned repairs applied (before normalization) and
  the H/D on the parents of the repaired atoms (the sites checked: other
  sites follow reduce2's usual rules, e.g. merged conformer copies)"""
  model = t.make_model(pdb_str) if cif_str is None else \
    t.make_model_cif(pdb_str, cif_str)
  plan, c, m = t.plan_of(model)
  h = model.get_hierarchy().deep_copy()
  origin = hd_exchange.apply_repairs(h, plan)
  repaired = set(origin) | set(k for x in plan.accept for k in (x.h_key,
    x.d_key))
  atoms = list(h.atoms())
  heavy = [a for a in atoms if a.element.strip() not in ("H", "D")]
  def parent(a):
    return min(heavy, key=lambda b: b.distance(a))
  parents = set(parent(a).memory_id() for a in atoms
                if hd_exchange.atom_key(a) in repaired)
  assert parents
  expected = [a for a in hd_atoms(h) if parent(a).memory_id() in parents]
  # pair sites: an H and a D within 0.20 A (their values are recorded)
  pair = set()
  for a in expected:
    for b in expected:
      if a.element.strip() != b.element.strip() and b.distance(a) <= 0.2:
        pair.add(a.memory_id())
  return expected, pair, h

def check_kept(expected, pair, out_text, label, keep):
  """every expected H/D value in the output, and no H/D at an expected
  site under another name (values swapped between sites). keep: values
  only at pair sites (recorded and restored); a supplied single gets
  placement's values there, as any supplied H (occupancy and B reset)"""
  out = hierarchy(out_text)
  if keep: expected = [a for a in expected if a.memory_id() in pair]
  want, have = values(expected), values(hd_atoms(out))
  for k, v in want.items():
    got = list(have.get(k, []))
    for x in v:
      assert x in got, (label, k, v, have.get(k))
      got.remove(x)
  out_atoms = hd_atoms(out)
  for a in expected:
    for b in out_atoms:
      if b.element.strip() != a.element.strip(): continue
      if b.distance(a) > 0.05: continue
      assert residue_name(b) == residue_name(a), (label, a.id_str(),
        b.id_str())

def check_repair(name, pdb_str, n_repairs, title, cif_str=None):
  expected, pair, h = expected_input(pdb_str, cif_str)
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
    check_kept(expected, pair, out_text, (name, tag), tag == "keep")
    # rerun on the output: nothing to repair, the H/D and values kept
    again = run(prefix + "_rerun", out_text, cif_str, **params)
    with open(again.params.output.description_file_name) as f:
      assert "Input H/D repaired" not in f.read(), (name, tag)
    with open(again.params.output.filename) as f:
      rerun_text = f.read()
    # the atoms keep their hierarchies while in use
    h_out, h_rerun = hierarchy(out_text), hierarchy(rerun_text)
    before, after = values(hd_atoms(h_out)), values(hd_atoms(h_rerun))
    assert before == after, (name, tag, [(k, before.get(k), after.get(k))
      for k in set(before) | set(after) if before.get(k) != after.get(k)])

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
       "parent copies with total occupancy 0 (no basis for sharing it)")):
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

def test_accept_refused():
  # (e) on a two-H parent (7KCU Gln A 18 NE2, D out of the amide plane)
  check_refused("gln18", t.pdb_ctl_7kcu_gln_a18,
    "NE2 carries 2 H in the dictionary (one H site needed); D ")

def test_different_parents_refused():
  # (c): 7T5E Tyr B 25, HH on OH A, DH on OH B
  check_refused("tyr25", t.pdb_ctl_7t5e_tyr_b25,
    "no repair for H and D on different parents")

if __name__ == "__main__":
  test_split()
  test_split_user_cif()
  test_split_refused()
  test_rename()
  test_rename_refused()
  test_relabel()
  test_accept()
  test_accept_refused()
  test_different_parents_refused()
  print("OK")
