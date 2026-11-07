"""hd_exchange: atom/site identity, H parent table, repairs of deposited H/D
layouts (no reduce2)."""
from __future__ import absolute_import, division, print_function
import math
import iotbx.pdb
import iotbx.phil
import mmtbx.model
from libtbx import group_args
from libtbx.utils import null_out, Sorry
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex
from scitbx import matrix
from cctbx import geometry_restraints, sgtbx
from cctbx.geometry_restraints.linking_class import linking_class
from mmtbx.monomer_library.pdb_interpretation import grand_master_phil_str
from mmtbx.hydrogens import reduce_hydrogen
from mmtbx.hydrogens import hd_exchange

# Ser with H HG in A and D DG in B
pdb_ser_hd = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
ATOM      1  N   SER A   2      14.767   2.090   9.884  1.00 20.00           N
ATOM      2  CA  SER A   2      15.853   1.130   9.701  1.00 20.00           C
ATOM      3  C   SER A   2      16.968   1.370  10.713  1.00 20.00           C
ATOM      4  O   SER A   2      17.281   2.516  11.041  1.00 20.00           O
ATOM      5  CB  SER A   2      16.411   1.224   8.278  1.00 20.00           C
ATOM      6  OG  SER A   2      15.456   0.826   7.311  1.00 20.00           O
ATOM      7  HA  SER A   2      15.444   0.123   9.808  1.00 20.00           H
ATOM      8  HB2 SER A   2      17.304   0.601   8.200  1.00 20.00           H
ATOM      9  HB3 SER A   2      16.665   2.256   8.039  1.00 20.00           H
ATOM     10  HG ASER A   2      15.856   0.874   6.441  0.50 20.00           H
ATOM     11  DG BSER A   2      15.856   0.874   6.441  0.50 20.00           D
END
"""

# 30F: Se-H from the library (1.444 A)
pdb_30f = """\
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1
HETATM    1  C10 30F A   1      12.392   9.810   9.963  1.00 20.00           C
HETATM    2  C08 30F A   1      10.909  10.206  10.137  1.00 20.00           C
HETATM    3  C09 30F A   1       9.982   9.241  10.307  1.00 20.00           C
HETATM    4  O2  30F A   1      12.678   8.610  10.088  1.00 20.00           O
HETATM    5 SE2  30F A   1       8.140   9.609  10.598  1.00 20.00          SE
HETATM    6  N04 30F A   1      10.682  11.544  10.125  1.00 20.00           N
HETATM    7  H13 30F A   1      10.297   8.208  10.345  1.00 20.00           H
HETATM    8  H2  30F A   1       7.642   9.275   9.267  1.00 20.00           H
HETATM    9  H1  30F A   1      11.429  12.105   9.752  1.00 20.00           H
HETATM   10  OXT 30F A   1      13.164  10.751   9.700  1.00 20.00           O
HETATM   11  H3  30F A   1       9.748  11.882   9.956  1.00 20.00           H
END
"""
# the library does not know HX, so its only partner comes from an edits bond
pdb_30f_hx = pdb_30f.replace(" H2  30F", " HX  30F")
# HX one cell away along x: bonded to Se through x-1,y,z
pdb_30f_hx_sym = pdb_30f_hx.replace("HX  30F A   1       7.642",
                                    "HX  30F A   1      27.642")

# the same with HX in its own residue before 30F, so that the H is the i
# endpoint of the asu proxy (above it is j)
pdb_30f_hx_sym_h_first = "\n".join(
  [pdb_30f_hx_sym.split("\n")[0]] +
  [l.replace("30F A   1", "HYD A   0") for l in pdb_30f_hx_sym.split("\n")
   if " HX " in l] +
  [l for l in pdb_30f_hx_sym.split("\n")[1:] if " HX " not in l])

# H1 on the 2-fold of P2: 0.85 A from O and from O's symmetry copy
pdb_water_special = """\
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1 2 1
HETATM    1  O   HOH A   1       0.527   5.000   0.668  1.00 20.00           O
HETATM    2  H1  HOH A   1       0.000   5.000   0.000  1.00 20.00           H
HETATM    3  H2  HOH A   1       1.360   5.000   0.850  1.00 20.00           H
END
"""

def make_model(pdb_str, edits=None, process=True):
  inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  model = mmtbx.model.manager(model_input=inp, log=null_out())
  if not process: return model
  model.set_stop_for_unknowns(False)
  p = reduce_hydrogen.get_reduce_pdb_interpretation_params(False)
  if edits:
    master = iotbx.phil.parse(grand_master_phil_str, process_includes=True)
    e = master.fetch(iotbx.phil.parse("\n".join(edits))).extract()
    p.geometry_restraints.edits.bond = e.geometry_restraints.edits.bond
  model.process(make_restraints=True, pdb_interpretation_params=p)
  return model

def edits_bond(sel1, sel2, distance_ideal, symmetry_operation=None):
  sym = ""
  if symmetry_operation:
    sym = "symmetry_operation = %s" % symmetry_operation
  return """geometry_restraints.edits.bond {
    atom_selection_1 = %s
    atom_selection_2 = %s
    distance_ideal = %s
    sigma = 0.02
    %s
  }""" % (sel1, sel2, distance_ideal, sym)

def entry(table, name):
  found = [v for k, v in table.items() if k[-2].strip() == name]
  assert len(found) == 1, (name, len(found))
  return found[0]

def test_parent_table_30f():
  table = hd_exchange.h_parent_table(make_model(pdb_30f))
  assert len(table) == 4
  for name, parent, element, exchangeable in (
      ("H13", "C09", "C", False), ("H2", "SE2", "SE", True),
      ("H1", "N04", "N", True), ("H3", "N04", "N", True)):
    e = entry(table, name)
    assert e.status == "ok", (name, e.status)
    assert parent in e.parent_id, (name, e.parent_id)
    assert e.parent_element == element
    assert e.exchangeable is exchangeable
    assert e.parent_occ == 1.0 and e.parent_altloc == ""
  # Se-H at 1.444 A from the library
  assert entry(table, "H2").source[0][1:3] == ("covalent geometry", 1.444)

def test_edits_covalent_limit():
  # 9: without restraints an H is floating; an edits Se-H at 1.444 A is
  # covalent (limit 1.20 + 0.31 + 0.4 = 1.91 A)
  e = entry(hd_exchange.h_parent_table(make_model(pdb_30f_hx)), "HX")
  assert e.status == "floating" and e.parent_id is None and e.source == []
  se_h = edits_bond("name SE2", "name HX", 1.444)
  e = entry(hd_exchange.h_parent_table(make_model(pdb_30f_hx, [se_h])), "HX")
  assert e.status == "ok" and "SE2" in e.parent_id and e.exchangeable
  assert e.source[0][1] == "edits"
  # an edits H...O at 1.9 A is above O's limit (0.66 + 0.31 + 0.4 = 1.37 A)
  h_o = edits_bond("name O2", "name HX", 1.9)
  e = entry(hd_exchange.h_parent_table(
    make_model(pdb_30f_hx, [se_h, h_o])), "HX")
  assert e.status == "ok" and "SE2" in e.parent_id
  assert len(e.source) == 1
  # below the limit it is a second partner
  h_o = edits_bond("name O2", "name HX", 1.3)
  e = entry(hd_exchange.h_parent_table(
    make_model(pdb_30f_hx, [se_h, h_o])), "HX")
  assert e.status == "ambiguous" and e.parent_id is None
  assert sorted(s[2] for s in e.source) == [1.3, 1.444], e.source
  assert e.exchangeable is None

def test_hydrogen_bond_origin_is_not_a_parent():
  # 9: a hydrogen-bond restraint (secondary-structure type) on an H
  model = make_model(pdb_ser_hd)
  atoms = model.get_hierarchy().atoms()
  i_hg = [a.i_seq for a in atoms if a.name.strip() == "HG"][0]
  i_o = [a.i_seq for a in atoms if a.name.strip() == "O"][0]
  proxy = geometry_restraints.bond_simple_proxy(i_seqs=(i_hg, i_o),
    distance_ideal=1.95, weight=1,
    origin_id=linking_class().get_origin_id("hydrogen bonds"))
  model.get_restraints_manager().geometry.add_new_bond_restraints_in_place(
    proxies=[proxy], sites_cart=model.get_sites_cart())
  e = entry(hd_exchange.h_parent_table(model), "HG")
  assert e.status == "ok" and "OG" in e.parent_id and len(e.source) == 1

def check_sym_op(model, h_name, parent_name):
  """each stored sym_op, applied to the parent, gives the image bonded to the
  H: distance within 0.05 A of the restraint target"""
  atoms = model.get_hierarchy().atoms()
  uc = model.crystal_symmetry().unit_cell()
  h = [a for a in atoms if a.name.strip() == h_name][0]
  parent = [a for a in atoms if a.name.strip() == parent_name][0]
  e = entry(hd_exchange.h_parent_table(model), h_name)
  for partner_id, origin, distance_ideal, sym_op in e.source:
    site_frac = sgtbx.rt_mx(sym_op) * uc.fractionalize(parent.xyz)
    d = (matrix.col(uc.orthogonalize(site_frac)) - matrix.col(h.xyz)).length()
    assert abs(d - distance_ideal) < 0.05, (sym_op, d, distance_ideal)
  return e

def test_symmetry():
  # 9: the only bond is to a symmetry mate (asu proxy), with the H as the i
  # and as the j endpoint; the wrong operator puts Se ~39.5 A away
  se_h = edits_bond("name SE2", "name HX", 1.444, "x-1,y,z")
  for pdb_str in (pdb_30f_hx_sym, pdb_30f_hx_sym_h_first):
    model = make_model(pdb_str, [se_h])
    e = check_sym_op(model, "HX", "SE2")
    assert e.status == "ok" and "SE2" in e.parent_id, e.status
    assert len(e.source) == 1 and e.source[0][1] == "edits"
    assert e.source[0][3] != "x,y,z"
  # an H bonded to its parent and to the parent's symmetry copy: one parent
  model = make_model(pdb_water_special)
  e = check_sym_op(model, "H1", "O")
  assert e.status == "ok" and " O  " in e.parent_id
  assert len(e.source) == 2
  assert sorted(s[3] == "x,y,z" for s in e.source) == [False, True], e.source

def test_hd_pair_keys():
  # 18a: H HG in A, D DG in B: no duplicate key, different atom keys, one
  # site key
  model = make_model(pdb_ser_hd)
  table = hd_exchange.h_parent_table(model)
  hg, dg = entry(table, "HG"), entry(table, "DG")
  assert hg.status == dg.status == "ok"
  assert hg.parent_id == dg.parent_id and hg.exchangeable and dg.exchangeable
  atoms = dict((a.name.strip(), a) for a in model.get_hierarchy().atoms())
  assert hd_exchange.atom_key(atoms["HG"]) != hd_exchange.atom_key(atoms["DG"])
  assert hd_exchange.site_key(atoms["HG"], "") == \
         hd_exchange.site_key(atoms["DG"], "")

def test_rekey_records():
  # records follow placement's renames of input H, in order: a swap of two
  # recorded names swaps the records, one recorded name moves to the other
  part = (0, "A", 0, "  37", " ", 0, "SER")
  hb2, hb3, h1, h = [part + ("", n) for n in (" HB2", " HB3", " H1 ", " H  ")]
  records = {hb2: "x", hb3: "y", h1: "z", part + ("A", " HB2"): "w"}
  renames = [(part, "", " HB2", " HB3"), (part, "", " H1 ", " H  ")]
  result = hd_exchange.rekey_records(records, renames)
  assert result == {hb2: "y", hb3: "x", h: "z", part + ("A", " HB2"): "w"}
  assert records[hb2] == "x"   # input unchanged
  # renamed twice: back where it was
  assert hd_exchange.rekey_records(records, renames[:1] * 2) == records

def test_name_normalization():
  # 18c: only D atoms are renamed, first D only
  h = iotbx.pdb.input(source_info=None, lines="""\
ATOM      1 HD21 ASN A   1       1.000   1.000   1.000  0.50 20.00           H
ATOM      2 DD21 ASN A   1       1.000   1.000   1.000  0.50 20.00           D
ATOM      3 1HB  ASN A   1       2.000   1.000   1.000  0.50 20.00           H
ATOM      4 1DB  ASN A   1       2.000   1.000   1.000  0.50 20.00           D
ATOM      5 HD22 ASN A   1       3.000   1.000   1.000  1.00 20.00           H
""".split("\n")).construct_hierarchy()
  atoms = dict((a.name, a) for a in h.atoms())
  names = dict((n, hd_exchange.normalized_name(a)) for n, a in atoms.items())
  assert names == {"HD21": "HD21", "DD21": "HD21", "1HB ": "1HB ",
    "1DB ": "1HB ", "HD22": "HD22"}, names
  sk = dict((n, hd_exchange.site_key(a, "")) for n, a in atoms.items())
  assert sk["HD21"] == sk["DD21"] and sk["1HB "] == sk["1DB "]
  assert len(set(sk.values())) == 3

def test_repeated_labels():
  # 16 (identity): two models with identical labels
  h = iotbx.pdb.input(source_info=None, lines="""\
MODEL        1
ATOM      1  OG  SER A   1       1.000   1.000   1.000  1.00 20.00           O
ATOM      2  HG  SER A   1       1.900   1.000   1.000  1.00 20.00           H
ENDMDL
MODEL        2
ATOM      1  OG  SER A   1       1.000   1.000   1.000  1.00 20.00           O
ATOM      2  HG  SER A   1       1.900   1.000   1.000  1.00 20.00           H
ENDMDL
""".split("\n")).construct_hierarchy()
  keys = hd_exchange.atom_keys(h)
  assert len(set(keys)) == 4 and keys[1][0] == 0 and keys[3][0] == 1
  # 17 (identity): repeated chain id (two separate chain A blocks with the
  # same resseq) and a duplicated residue number within one chain
  h = iotbx.pdb.input(source_info=None, lines="""\
ATOM      1  OG  SER A   1       1.000   1.000   1.000  1.00 20.00           O
ATOM      2  HG  SER A   1       1.900   1.000   1.000  1.00 20.00           H
ATOM      3  OG  SER B   1       5.000   1.000   1.000  1.00 20.00           O
ATOM      4  HG  SER B   1       5.900   1.000   1.000  1.00 20.00           H
ATOM      5  OG  SER A   1       9.000   1.000   1.000  1.00 20.00           O
ATOM      6  HG  SER A   1       9.900   1.000   1.000  1.00 20.00           H
ATOM      7  OG  SER B   1      12.000   1.000   1.000  1.00 20.00           O
ATOM      8  HG  SER B   1      12.900   1.000   1.000  1.00 20.00           H
ATOM      9  OG  SER B   2      15.000   1.000   1.000  1.00 20.00           O
ATOM     10  OG  SER B   1      18.000   1.000   1.000  1.00 20.00           O
ATOM     11  HG  SER B   1      18.900   1.000   1.000  1.00 20.00           H
""".split("\n")).construct_hierarchy()
  keys = hd_exchange.atom_keys(h)
  assert len(set(keys)) == 11
  # atom_keys() and atom_key() agree
  assert keys == [hd_exchange.atom_key(a) for a in h.atoms()]
  hg = [k for k in keys if k[-2] == " HG "]
  assert [(k[1], k[2], k[5]) for k in hg] == [
    ("A", 0, 0), ("B", 0, 0), ("A", 1, 0), ("B", 1, 0), ("B", 1, 1)], hg

def test_duplicate_is_an_error():
  h = iotbx.pdb.input(source_info=None, lines="""\
ATOM      1  OG  SER A   1       1.000   1.000   1.000  1.00 20.00           O
ATOM      2  HG  SER A   1       1.900   1.000   1.000  1.00 20.00           H
ATOM      3  HG  SER A   1       2.900   1.000   1.000  1.00 20.00           H
""".split("\n")).construct_hierarchy()
  try:
    hd_exchange.atom_keys(h)
  except Sorry as e:
    assert "duplicate atom identity" in str(e)
    assert str(e).count('HG  SER A   1') == 2, str(e)
  else:
    raise AssertionError("duplicate atom not reported")

def test_water_canonical_resname():
  for resname in ("HOH", "DOD", "WAT"):
    assert hd_exchange.canonical_resname(resname) == "HOH"
  assert hd_exchange.canonical_resname("SER") == "SER"

def test_keys_survive_copies():
  # 18: atom_key() is unchanged by deep_copy(), sort_atoms_in_place() and
  # select() of unrelated atoms; site_key() stays when H moves to blank
  model = make_model(pdb_ser_hd)
  keys = dict((a.name, hd_exchange.atom_key(a))
              for a in model.get_hierarchy().atoms())
  copy = model.deep_copy()
  assert dict((a.name, hd_exchange.atom_key(a))
              for a in copy.get_hierarchy().atoms()) == keys
  h = model.get_hierarchy().deep_copy()
  h.sort_atoms_in_place()
  assert dict((a.name, hd_exchange.atom_key(a)) for a in h.atoms()) == keys
  sel = flex.bool([a.name.strip() != "HA" for a in
                   model.get_hierarchy().atoms()])
  sub = model.select(sel)
  sub_keys = dict((a.name, hd_exchange.atom_key(a))
                  for a in sub.get_hierarchy().atoms())
  assert sub_keys == dict((k, v) for k, v in keys.items() if k != " HA ")
  # move HG from altloc A into the blank atom_group
  h = model.get_hierarchy().deep_copy()
  rg = h.only_residue_group()
  hg = [a for a in h.atoms() if a.name == " HG "][0]
  before = hd_exchange.site_key(hg, "")
  ag_a = hg.parent()
  blank = [ag for ag in rg.atom_groups() if ag.altloc == ""][0]
  moved = hg.detached_copy()
  ag_a.remove_atom(hg)
  blank.append_atom(moved)
  assert hd_exchange.site_key(moved, "") == before
  assert hd_exchange.atom_key(moved) != keys[" HG "]

# ------------------------------------------------------------------------------
# 2.2: candidates and collapse check

# Asn with amide H in A and D in B; positions P1 (21) and P2 (22)
pdb_asn_hd = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
ATOM      1  N   ASN A   2      12.344  12.923  13.155  1.00 15.02           N
ATOM      2  CA  ASN A   2      13.478  12.038  12.831  1.00 14.10           C
ATOM      3  C   ASN A   2      14.759  12.537  13.427  1.00 13.13           C
ATOM      4  O   ASN A   2      15.022  13.742  13.426  1.00 11.91           O
ATOM      5  CB  ASN A   2      13.654  11.881  11.341  1.00 15.38           C
ATOM      6  CG  ASN A   2      12.416  11.342  10.692  1.00 14.08           C
ATOM      7  OD1 ASN A   2      11.975  10.227  11.016  1.00 17.46           O
ATOM      8  ND2 ASN A   2      11.796  12.155   9.831  1.00 11.72           N
ATOM      9 HD21AASN A   2      10.985  11.856   9.379  0.50 15.00           H
ATOM     10 HD22AASN A   2      12.141  13.050   9.648  0.50 15.00           H
ATOM     11 DD21BASN A   2      10.985  11.856   9.379  0.50 15.00           D
ATOM     12 DD22BASN A   2      12.141  13.050   9.648  0.50 15.00           D
END
"""

def check_indices(hierarchy, c):
  """each index field names the atom of its label"""
  atoms = hierarchy.atoms()
  for p in c.pairs:
    assert atoms[p.i_h].id_str() == p.h_id and atoms[p.i_d].id_str() == p.d_id
    if p.geometric_partner is None:
      assert p.i_geometric_partner is None
    else:
      assert atoms[p.i_geometric_partner].id_str() == p.geometric_partner
  for x in c.singles + c.unresolved + c.same_conformer:
    assert atoms[x.i_seq].id_str() == x.atom_id

def candidates_of(model, xyz_tolerance=0.2):
  table = hd_exchange.h_parent_table(model)
  c = hd_exchange.find_hd_candidates(model.get_hierarchy(), table, model=model)
  check_indices(model.get_hierarchy(), c)
  return c, hd_exchange.collapse_check(c, xyz_tolerance=xyz_tolerance)

def names(atom_id):
  return atom_id[5:9].strip()

def pair_names(pairs):
  return sorted((names(p.h_id), names(p.d_id)) for p in pairs)

def test_basic_pair():
  c, r = candidates_of(make_model(pdb_ser_hd))
  assert pair_names(r.pairs) == [("HG", "DG")]
  p = r.pairs[0]
  assert p.distance < 1e-3 and not p.name_mismatch and p.n_competing == 0
  assert p.layout_ok and p.parent_altloc == "" and p.rotatable
  assert r.problems == [] and r.unresolved == [] and r.notes == []
  # the C-H are singles
  assert sorted(names(x.atom_id) for x in c.singles) == ["HA", "HB2", "HB3"]

def test_pair_0p5_apart():
  # 1a: kept as a candidate with its distance; collapse_check: problem
  pdb = pdb_ser_hd.replace(
    "DG BSER A   2      15.856   0.874   6.441",
    "DG BSER A   2      16.356   0.874   6.441")
  c, r = candidates_of(make_model(pdb))
  assert pair_names(c.pairs) == [("HG", "DG")]
  assert abs(c.pairs[0].distance - 0.5) < 1e-3
  assert r.pairs == [] and len(r.problems) == 1
  assert "0.50 A apart" in r.problems[0].reasons[0], r.problems[0].reasons

def test_swapped_names():
  # 1b: control, then H HD21 on D DD22's position and HD22 on DD21's
  c, r = candidates_of(make_model(pdb_asn_hd))
  assert pair_names(r.pairs) == [("HD21", "DD21"), ("HD22", "DD22")]
  assert not any(p.name_mismatch for p in c.pairs)
  pdb = pdb_asn_hd.replace("DD21BASN", "DDXXBASN").replace(
    "DD22BASN", "DD21BASN").replace("DDXXBASN", "DD22BASN")
  c, r = candidates_of(make_model(pdb))
  assert pair_names(c.pairs) == [("HD21", "DD21"), ("HD22", "DD22")]
  for p in c.pairs:
    assert p.name_mismatch and p.geometric_partner is not None
    other = {"HD21": "DD22", "HD22": "DD21"}[names(p.h_id)]
    assert names(p.geometric_partner) == other, p.geometric_partner
    assert p.i_geometric_partner not in (None, p.i_h, p.i_d)
  assert r.pairs == [] and len(r.problems) == 2
  for prob in r.problems:
    assert any(x.startswith("names swapped") for x in prob.reasons)
  # control: no name mismatch, no geometric partner
  c, r = candidates_of(make_model(pdb_asn_hd))
  assert [p.i_geometric_partner for p in c.pairs] == [None, None]

def test_one_h_two_d():
  # 2: one H with two D candidates
  pdb = pdb_ser_hd.replace("END",
    "ATOM     12  DG CSER A   2      15.856   0.874   6.441  0.50 20.00"
    "           D\nEND")
  c, r = candidates_of(make_model(pdb))
  assert len(c.pairs) == 1 and c.pairs[0].n_competing == 1
  assert r.pairs == []
  assert any("competing" in x for x in r.problems[0].reasons)

def test_v2_names():
  # 3: 1HB / 1DB
  pdb = pdb_ser_hd.replace(
    " HB2 SER A   2      17.304   0.601   8.200  1.00",
    "1HB ASER A   2      17.304   0.601   8.200  0.50").replace("END",
    "ATOM     12 1DB BSER A   2      17.304   0.601   8.200  0.50 20.00"
    "           D\nEND")
  c, r = candidates_of(make_model(pdb))
  assert pair_names(r.pairs) == [("1HB", "1DB"), ("HG", "DG")]

def test_pair_in_one_conformer():
  # 4: an H and a D of one site in compatible atom_groups (same altloc, or
  # one blank): pdb_interpretation restrains only one of them; the other is
  # floating. That is a stop, not a removal. First: OG in A and B, HG and DG
  # both in A.
  pdb = pdb_ser_hd.replace(
    "ATOM      6  OG  SER A   2      15.456   0.826   7.311  1.00",
    "ATOM      6  OG ASER A   2      15.456   0.826   7.311  0.50").replace(
    "DG BSER", "DG ASER").replace("END",
    "ATOM     12  OG BSER A   2      16.000   0.300   7.100  0.50 20.00"
    "           O\nEND")
  c, r = candidates_of(make_model(pdb))
  assert c.pairs == [] and len(c.same_conformer) == 1
  assert r.pairs == [] and r.unresolved == []
  assert [p.kind for p in r.problems] == ["same conformer"]
  assert "H/D pair in one conformer not supported" in \
    r.problems[0].reasons[0]
  removed = [names(u.atom_id) for u in r.unresolved if u.action == "remove"]
  assert "HG" not in removed and "DG" not in removed
  # compatible atom_groups: one isotope in the blank group, the other in B,
  # both orders; the one in B is floating and must not be removed
  hg_blank = pdb_ser_hd.replace(
    "HG ASER A   2      15.856   0.874   6.441  0.50",
    "HG  SER A   2      15.856   0.874   6.441  1.00")
  dg_blank = pdb_ser_hd.replace(
    "HG ASER A   2      15.856   0.874   6.441  0.50",
    "HG BSER A   2      15.856   0.874   6.441  0.50").replace(
    "DG BSER A   2      15.856   0.874   6.441  0.50",
    "DG  SER A   2      15.856   0.874   6.441  1.00")
  for pdb, floating in ((hg_blank, "DG"), (dg_blank, "HG")):
    c, r = candidates_of(make_model(pdb))
    assert c.pairs == [] and r.pairs == [] and r.unresolved == []
    assert [names(x.atom_id) for x in c.same_conformer] == [floating]
    assert c.same_conformer[0].i_seq is not None
    assert [p.kind for p in r.problems] == ["same conformer"]
    assert "H/D pair in one conformer not supported" in \
      r.problems[0].reasons[0]
  # a floating H/D without a same-name partner is still removed
  pdb = pdb_30f_hx.replace("END",
    "HETATM   12  HY  30F A   1       8.000   9.000   9.000  1.00 20.00"
    "           H\nEND")
  c, r = candidates_of(make_model(pdb))
  assert c.same_conformer == [] and r.problems == []
  assert sorted((names(u.atom_id), u.action) for u in r.unresolved) == [
    ("HX", "remove"), ("HY", "remove")]
  assert sorted(names(u.atom_id) for u in c.unresolved
                if u.i_seq is not None) == ["HX", "HY"]

def test_layout_needs_two_atom_groups():
  # layout_ok: H and D in different atom_groups. Real restraints cannot give
  # both a parent inside one atom_group, so the parent table is built by hand.
  h = iotbx.pdb.input(source_info=None, lines="""\
ATOM      1  OG  SER A   2      15.456   0.826   7.311  1.00 20.00           O
ATOM      2  HG ASER A   2      15.856   0.874   6.441  0.50 20.00           H
ATOM      3  DG ASER A   2      15.856   0.874   6.441  0.50 20.00           D
ATOM      4  HG BSER A   2      15.856   0.874   6.441  0.50 20.00           H
ATOM      5  DG CSER A   2      15.856   0.874   6.441  0.50 20.00           D
""".split("\n")).construct_hierarchy()
  keys = hd_exchange.atom_keys(h)
  og = keys[0]
  table = {}
  for key in keys[1:]:
    table[key] = group_args(status="ok", parent_key=og,
      parent_sym_ops=("x,y,z",), parent_id='pdb=" OG  SER A   2 "',
      parent_occ=1.0, parent_altloc="", source=[])
  c = hd_exchange.find_hd_candidates(h, table)
  layout = dict(((p.h_id[9], p.d_id[9]), p.layout_ok) for p in c.pairs)
  assert layout == {("A", "A"): False, ("B", "C"): True}, layout

def test_renamed_d_collides():
  # 5 with real restraints: HG and DG in A (DG 1 A away), another DG in B on
  # HG. HG is floating next to the bonded DG in its atom_group, so this is
  # the same-conformer problem; the A DG is not reported again as a
  # collision.
  pdb = pdb_ser_hd.replace(
    "DG BSER A   2      15.856   0.874   6.441",
    "DG ASER A   2      15.856   1.874   6.441").replace("END",
    "ATOM     12  DG BSER A   2      15.856   0.874   6.441  0.50 20.00"
    "           D\nEND")
  c, r = candidates_of(make_model(pdb))
  assert [x.kind for x in r.problems] == ["same conformer"], \
    [x.kind for x in r.problems]
  assert r.unresolved == []
  # 5 as such: a single D whose H name already exists in its atom_group.
  # Real restraints leave one of the two floating (above), so the parent
  # table is built by hand: the H is paired with a D in B, the A D is single.
  h = iotbx.pdb.input(source_info=None, lines="""\
ATOM      1  OG  SER A   2      15.456   0.826   7.311  1.00 20.00           O
ATOM      2  HG ASER A   2      15.856   0.874   6.441  0.50 20.00           H
ATOM      3  DG ASER A   2      15.856   1.874   6.441  0.50 20.00           D
ATOM      4  DG BSER A   2      15.856   0.874   6.441  0.50 20.00           D
""".split("\n")).construct_hierarchy()
  keys = hd_exchange.atom_keys(h)
  table = dict((key, group_args(status="ok", parent_key=keys[0],
    parent_sym_ops=("x,y,z",), parent_id='pdb=" OG  SER A   2 "',
    parent_occ=1.0, parent_altloc="", source=[])) for key in keys[1:])
  c = hd_exchange.find_hd_candidates(h, table)
  r = hd_exchange.collapse_check(c)
  dg = sorted((x.atom_id[9], x.name_collision) for x in c.singles)
  assert dg == [("A", True)], dg
  assert "name collision" in [x.kind for x in r.problems]

def test_notes():
  # 6: occupancy sum and B differences are notes, never problems
  pdb = pdb_ser_hd.replace(
    "HG ASER A   2      15.856   0.874   6.441  0.50 20.00",
    "HG ASER A   2      15.856   0.874   6.441  0.30 20.00").replace(
    "DG BSER A   2      15.856   0.874   6.441  0.50 20.00",
    "DG BSER A   2      15.856   0.874   6.441  0.60 30.00")
  c, r = candidates_of(make_model(pdb))
  assert pair_names(r.pairs) == [("HG", "DG")] and r.problems == []
  kinds = dict((n.kind, n) for n in r.notes)
  assert set(kinds) == set(["occupancy", "b"])
  assert not kinds["occupancy"].sum_is_one
  # parent at 0.6, H + D = 1 (phenix.refine convention)
  pdb = pdb_ser_hd.replace(
    "OG  SER A   2      15.456   0.826   7.311  1.00",
    "OG  SER A   2      15.456   0.826   7.311  0.60")
  c, r = candidates_of(make_model(pdb))
  assert pair_names(r.pairs) == [("HG", "DG")] and r.problems == []
  assert [n.kind for n in r.notes] == ["occupancy"]
  assert r.notes[0].sum_is_one

def test_microheterogeneity():
  # 7: SER in A with HG, CYS in B with DG: no cross-resname pair
  pdb = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
ATOM      1  N  ASER A   2      14.767   2.090   9.884  0.50 20.00           N
ATOM      2  CA ASER A   2      15.853   1.130   9.701  0.50 20.00           C
ATOM      3  C  ASER A   2      16.968   1.370  10.713  0.50 20.00           C
ATOM      4  O  ASER A   2      17.281   2.516  11.041  0.50 20.00           O
ATOM      5  CB ASER A   2      16.411   1.224   8.278  0.50 20.00           C
ATOM      6  OG ASER A   2      15.456   0.826   7.311  0.50 20.00           O
ATOM      7  HG ASER A   2      15.856   0.874   6.441  0.50 20.00           H
ATOM      8  N  BCYS A   2      14.767   2.090   9.884  0.50 20.00           N
ATOM      9  CA BCYS A   2      15.853   1.130   9.701  0.50 20.00           C
ATOM     10  C  BCYS A   2      16.968   1.370  10.713  0.50 20.00           C
ATOM     11  O  BCYS A   2      17.281   2.516  11.041  0.50 20.00           O
ATOM     12  CB BCYS A   2      16.411   1.224   8.278  0.50 20.00           C
ATOM     13  SG BCYS A   2      15.210   0.700   6.900  0.50 20.00           S
ATOM     14  DG BCYS A   2      15.856   0.874   6.441  0.50 20.00           D
END
"""
  c, r = candidates_of(make_model(pdb))
  assert c.pairs == [] and c.different_parents == []
  assert sorted(names(x.atom_id) for x in c.singles) == ["DG", "HG"]

def test_one_isotope_per_site():
  # 7a: sharing a parent does not make a pair. Water O/H1/D2 and Asn HD21
  # (H) + DD22 (D) at their own positions are singles, no problem.
  water = """\
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1
HETATM    1  O   HOH A   1      10.000  10.000  10.000  1.00 20.00           O
HETATM    2  H1  HOH A   1      10.850  10.000  10.000  1.00 20.00           H
HETATM    3  D2  HOH A   1       9.795  10.825  10.000  1.00 20.00           D
END
"""
  asn = pdb_asn_hd.replace(
    "ATOM     10 HD22AASN A   2      12.141  13.050   9.648  0.50 15.00"
    "           H\n", "").replace(
    "ATOM     11 DD21BASN A   2      10.985  11.856   9.379  0.50 15.00"
    "           D\n", "")
  for pdb, expected in ((water, ["D2", "H1"]), (asn, ["DD22", "HD21"])):
    c, r = candidates_of(make_model(pdb))
    assert c.pairs == [] and r.problems == [] and r.unresolved == []
    assert sorted(names(x.atom_id) for x in c.singles
                  if names(x.atom_id) in expected) == expected
  # control: HD21 (H) and DD22 (D) within 0.5 A of each other
  asn_close = asn.replace(
    "DD22BASN A   2      12.141  13.050   9.648",
    "DD22BASN A   2      11.185  11.856   9.379")
  c, r = candidates_of(make_model(asn_close))
  assert pair_names(c.pairs) == [("HD21", "DD22")]
  assert c.pairs[0].name_mismatch and abs(c.pairs[0].distance - 0.2) < 1e-3
  # matched by position: the geometric partner is the pair's own D
  assert c.pairs[0].i_geometric_partner == c.pairs[0].i_d
  assert r.pairs == [] and len(r.problems) == 1
  assert any("names" in x for x in r.problems[0].reasons)

def test_ambiguous_is_unresolved():
  # 8b: a pair whose H has an ambiguous parent is unresolved, with the
  # restraints named, not an invalid pair
  pdb = pdb_30f_hx.replace(
    "HX  30F A   1       7.642   9.275   9.267  1.00",
    "HX A30F A   1       7.642   9.275   9.267  0.50").replace("END",
    "HETATM   12  DX B30F A   1       7.642   9.275   9.267  0.50 20.00"
    "           D\nEND")
  edits = [edits_bond("name SE2", "name HX", 1.444),
           edits_bond("name SE2", "name DX", 1.444),
           edits_bond("name O2", "name HX", 1.3)]
  c, r = candidates_of(make_model(pdb, edits))
  assert c.pairs == [] and r.problems == []
  assert [(names(u.atom_id), u.status, u.action) for u in r.unresolved] == [
    ("HX", "ambiguous", "stop")]
  assert sorted(x[1] for x in r.unresolved[0].source) == ["edits", "edits"]
  assert sorted(x[2] for x in r.unresolved[0].source) == [1.3, 1.444]
  assert [names(u.atom_id) for u in c.unresolved] == ["HX"]
  assert c.unresolved[0].i_seq is not None

def test_symmetry_parent():
  # 8a (symmetry part): H and D both bonded to Se across the boundary
  pdb = pdb_30f_hx.replace(
    "HX  30F A   1       7.642   9.275   9.267  1.00",
    "HX A30F A   1      27.642   9.275   9.267  0.50").replace("END",
    "HETATM   12  DX B30F A   1      27.642   9.275   9.267  0.50 20.00"
    "           D\nEND")
  edits = [edits_bond("name SE2", "name HX", 1.444, "x-1,y,z"),
           edits_bond("name SE2", "name DX", 1.444, "x-1,y,z")]
  c, r = candidates_of(make_model(pdb, edits))
  assert pair_names(r.pairs) == [("HX", "DX")] and r.problems == []
  # H bonded to Se, D bonded to Se's symmetry copy: different parents
  pdb = pdb.replace("HX A30F A   1      27.642", "HX A30F A   1       7.642")
  edits = [edits_bond("name SE2", "name HX", 1.444),
           edits_bond("name SE2", "name DX", 1.444, "x-1,y,z")]
  c, r = candidates_of(make_model(pdb, edits))
  assert c.pairs == [] and len(c.different_parents) == 1
  assert [p.kind for p in r.problems] == ["different parents"]

def test_different_parents_5A():
  # 1: same-name H/D on different parents, 5 A apart: listed, not a pair
  pdb = pdb_30f_hx.replace(
    "HX  30F A   1       7.642   9.275   9.267  1.00",
    "HX A30F A   1       7.642   9.275   9.267  0.50").replace("END",
    "HETATM   12  DX B30F A   1      11.200  12.400  10.500  0.50 20.00"
    "           D\nEND")
  edits = [edits_bond("name SE2", "name HX", 1.444),
           edits_bond("name N04", "name DX", 1.0)]
  c, r = candidates_of(make_model(pdb, edits))
  assert len(c.different_parents) == 1
  dp = c.different_parents[0]
  assert (names(dp.h_id), names(dp.d_id)) == ("HX", "DX")
  assert abs(dp.distance - 5.0) < 0.3
  # HX and DX are no pair; DX is left over on N04 0.83 A from H1 there:
  # different sites, so both are singles
  assert c.pairs == [] and r.pairs == []
  assert [p.kind for p in r.problems] == ["different parents"]
  assert sorted(names(x.atom_id) for x in c.singles
                if names(x.atom_id) in ("DX", "H1")) == ["DX", "H1"]

def test_reversed_order():
  # 11: D in A (its atom_group first), H in B
  pdb = pdb_ser_hd.replace(
    "ATOM     10  HG ASER A   2      15.856   0.874   6.441  0.50 20.00"
    "           H\n"
    "ATOM     11  DG BSER A   2      15.856   0.874   6.441  0.50 20.00"
    "           D\n",
    "ATOM     10  DG ASER A   2      15.856   0.874   6.441  0.50 20.00"
    "           D\n"
    "ATOM     11  HG BSER A   2      15.856   0.874   6.441  0.50 20.00"
    "           H\n")
  model = make_model(pdb)
  ags = [ag.altloc for ag in model.get_hierarchy().atom_groups()]
  assert ags.index("A") < ags.index("B")
  c, r = candidates_of(model)
  assert pair_names(r.pairs) == [("HG", "DG")] and r.problems == []

def test_alternative_h_orientations():
  # 18b (candidate part): HG in A and another HG in B, both H
  pdb = pdb_ser_hd.replace(
    "DG BSER A   2      15.856   0.874   6.441  0.50 20.00           D",
    "HG BSER A   2      14.956   0.300   6.900  0.50 20.00           H")
  c, r = candidates_of(make_model(pdb))
  assert c.pairs == [] and r.problems == []
  hg = sorted((x.atom_id[9], names(x.atom_id)) for x in c.singles
              if names(x.atom_id) == "HG")
  assert hg == [("A", "HG"), ("B", "HG")], hg

# ------------------------------------------------------------------------------
# 2.3: normalization, output, restoration

# Ala-Ser-Lys and a water; H are placed by reduce_hydrogen
pdb_pep = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
ATOM      1  N   ALA A   1      11.204   2.315  10.312  1.00 20.00           N
ATOM      2  CA  ALA A   1      12.555   2.835  10.528  1.00 20.00           C
ATOM      3  C   ALA A   1      13.588   1.726  10.381  1.00 20.00           C
ATOM      4  O   ALA A   1      13.325   0.567  10.721  1.00 20.00           O
ATOM      5  CB  ALA A   1      12.866   3.965   9.556  1.00 20.00           C
ATOM      6  N   SER A   2      14.767   2.090   9.884  1.00 20.00           N
ATOM      7  CA  SER A   2      15.853   1.130   9.701  1.00 20.00           C
ATOM      8  C   SER A   2      16.968   1.370  10.713  1.00 20.00           C
ATOM      9  O   SER A   2      17.281   2.516  11.041  1.00 20.00           O
ATOM     10  CB  SER A   2      16.411   1.224   8.278  1.00 20.00           C
ATOM     11  OG  SER A   2      15.456   0.826   7.311  1.00 20.00           O
ATOM     12  N   LYS A   3      17.560   0.282  11.199  1.00 20.00           N
ATOM     13  CA  LYS A   3      18.643   0.370  12.172  1.00 20.00           C
ATOM     14  C   LYS A   3      19.945  -0.148  11.570  1.00 20.00           C
ATOM     15  O   LYS A   3      19.964  -1.176  10.892  1.00 20.00           O
ATOM     16  CB  LYS A   3      18.290  -0.424  13.432  1.00 20.00           C
ATOM     17  CG  LYS A   3      19.308  -0.301  14.556  1.00 20.00           C
ATOM     18  CD  LYS A   3      18.913  -1.137  15.762  1.00 20.00           C
ATOM     19  CE  LYS A   3      19.919  -1.002  16.894  1.00 20.00           C
ATOM     20  NZ  LYS A   3      19.531  -1.812  18.059  1.00 20.00           N
ATOM     21  OXT LYS A   3      20.995   0.451  11.763  1.00 20.00           O
HETATM   22  O   HOH A   4      24.000   6.000  20.000  1.00 20.00           O
END
"""
EXCHANGEABLE_PEP = 9   # Ala H1 H2 H3, Ser H HG, Lys H HZ1 HZ2 HZ3
ALL_H_PEP = 25         # without the water

def model_from_str(pdb_str):
  return make_model(pdb_str)

def placed(pdb_str):
  model = make_model(pdb_str, process=False)
  o = reduce_hydrogen.place_hydrogens(model=model, exclude_water=False)
  o.run()
  return o.get_model()

def model_from_hierarchy(hierarchy, crystal_symmetry):
  return make_model(hierarchy.as_pdb_string(crystal_symmetry=crystal_symmetry))

def atom_values(hierarchy):
  return dict((hd_exchange.atom_key(a), (a.name, a.xyz, a.occ, a.b,
    tuple(a.uij) if a.uij_is_defined() else None))
    for a in hierarchy.atoms())

def heavy(hierarchy):
  return dict((k, v) for k, v in atom_values(hierarchy).items()
              if k[-1] not in ("H", "D"))

def check_output(h_out, table, r, h_in):
  """output check, exchangeable_hd_selections(), heavy atoms, numbering"""
  c = hd_exchange.find_hd_candidates(h_out, table, provenance=r.provenance)
  k = hd_exchange.collapse_check(c)
  assert k.problems == [], [(p.kind, p.reasons) for p in k.problems]
  assert k.unresolved == []
  assert len(k.pairs) == r.n_split, (len(k.pairs), r.n_split)
  pairs = set(frozenset((p.i_h, p.i_d)) for p in k.pairs)
  exchangeable = set(frozenset((x[0][0], x[1][0]))
                     for x in h_out.exchangeable_hd_selections())
  assert exchangeable == pairs
  # provenance: a bijection onto the written H/D
  hd = [hd_exchange.atom_key(a) for a in h_out.atoms()
        if a.element.strip().upper() in ("H", "D")]
  assert sorted(hd) == sorted(r.provenance)
  # heavy atoms untouched (canonical water name in the key)
  assert heavy(h_out) == heavy(h_in)
  atoms = h_out.atoms()
  assert list(atoms.extract_i_seq()) == list(range(atoms.size()))
  assert [int(a.serial) for a in atoms] == list(range(1, atoms.size() + 1))
  return k

def run_option(model, option, records=None):
  h_in = model.get_hierarchy()
  h = h_in.deep_copy()
  table = hd_exchange.h_parent_table(model)
  r = hd_exchange.apply_neutron_option(h, option, table, records=records)
  check_output(h, table, r, h_in)
  return h, r

def normalize(pdb_str):
  """H/D input -> (H-only model, records)"""
  model = model_from_str(pdb_str)
  h = model.get_hierarchy()
  table = hd_exchange.h_parent_table(model)
  match = hd_exchange.collapse_check(hd_exchange.find_hd_candidates(h, table))
  assert match.problems == [] and match.unresolved == []
  heavy_in = heavy(h)
  result = hd_exchange.normalize_to_h(h, match)
  assert heavy(h) == heavy_in
  assert not any(a.element.strip().upper() == "D" for a in h.atoms())
  atoms = h.atoms()
  assert list(atoms.extract_i_seq()) == list(range(atoms.size()))
  assert [int(a.serial) for a in atoms] == list(range(1, atoms.size() + 1))
  return model_from_hierarchy(h, model.crystal_symmetry()), result.records

def find(hierarchy, resname, name, altloc=None):
  found = [a for a in hierarchy.atoms() if a.parent().resname == resname
           and a.name.strip() == name and
           (altloc is None or a.parent().altloc == altloc)]
  assert len(found) == 1, (resname, name, altloc, len(found))
  return found[0]

def values(atom):
  return (round(atom.occ, 3), round(atom.b, 3),
          tuple(atom.uij) if atom.uij_is_defined() else None)

_pep_h = []
def pep_h():
  if not _pep_h: _pep_h.append(placed(pdb_pep))
  return _pep_h[0]

def hd_input(option, edit=None):
  """pdb string of the tripeptide written in option, with edit(hierarchy)"""
  h, r = run_option(pep_h(), option)
  if edit is not None: edit(h)
  return h.as_pdb_string(crystal_symmetry=pep_h().crystal_symmetry())

def test_modes():
  # "What each mode writes" on H-only input (with the water)
  model = pep_h()
  # none: hierarchy unchanged, empty results
  h = model.get_hierarchy().deep_copy()
  before = h.as_pdb_string()
  r = hd_exchange.apply_neutron_option(h, "none",
    hd_exchange.h_parent_table(model))
  assert h.as_pdb_string() == before
  assert (r.n_split, r.n_to_d, r.n_water, r.n_restored) == (0, 0, 0, 0)
  assert r.provenance == {} and r.merged_sites == [] and r.notes == []
  for option, n_split, n_d, water in (
      ("all_h", 0, 0, "HOH"),
      ("all_d", 0, ALL_H_PEP + 2, "DOD"),
      ("hd_and_h", EXCHANGEABLE_PEP, 2, "DOD"),
      ("hd_and_d", EXCHANGEABLE_PEP, ALL_H_PEP - EXCHANGEABLE_PEP + 2, "DOD"),
      ("all_hd", ALL_H_PEP, 2, "DOD")):
    h, r = run_option(model, option)
    assert (r.n_split, r.n_to_d) == (n_split, n_d), (option, r.n_split,
      r.n_to_d)
    w = [ag for ag in h.atom_groups()
         if iotbx.pdb.common_residue_names_get_class(ag.resname) ==
            "common_water"]
    assert [ag.resname for ag in w] == [water]
    e = set(a.element.strip() for a in w[0].atoms()) - set(["O"])
    assert e == set(["H" if water == "HOH" else "D"])
    # new pairs: 0.5 x parent occupancy, H's B
    for a in h.atoms():
      if a.parent().altloc != "":
        assert abs(a.occ - 0.5) < 1e-6 and abs(a.b - 20.0) < 1e-6
    assert r.conformer_sites == [] and r.unclassified_sites == []

def test_altloc_groups():
  # 12: Ser H and HG in A, D and DG in B; exactly two new groups per residue
  h, r = run_option(pep_h(), "hd_and_h")
  for resseq, h_names in (("   1", ["H1", "H2", "H3"]), ("   2", ["H", "HG"]),
                          ("   3", ["H", "HZ1", "HZ2", "HZ3"])):
    rg = [x for x in h.only_chain().residue_groups()
          if x.resseq == resseq][0]
    by_alt = dict((ag.altloc, sorted(a.name.strip() for a in ag.atoms()))
                  for ag in rg.atom_groups())
    assert sorted(by_alt) == ["", "A", "B"], by_alt
    assert by_alt["A"] == sorted(h_names)
    assert by_alt["B"] == sorted("D" + n[1:] for n in h_names)

def test_altloc_groups_conformers():
  # 12: OG/HG in A and B, backbone H blank: HG is a conformer site (D in A
  # and in B), the backbone H of that residue is split into C/D.
  # (place_hydrogens would put every H of the residue into A and B.)
  h = pep_h().get_hierarchy().deep_copy()
  og, hg = find(h, "SER", "OG"), find(h, "SER", "HG")
  rg = og.parent().parent()
  for altloc, shift in (("A", 0.0), ("B", 0.6)):
    ag = iotbx.pdb.hierarchy.atom_group(altloc=altloc, resname="SER")
    for atom in (og, hg):
      a = atom.detached_copy()
      a.set_xyz((atom.xyz[0] + shift, atom.xyz[1], atom.xyz[2]))
      a.set_occ(0.5)
      ag.append_atom(a)
    rg.append_atom_group(ag)
  og.parent().remove_atom(og)
  hg.parent().remove_atom(hg)
  model = model_from_hierarchy(h, pep_h().crystal_symmetry())
  h, r = run_option(model, "hd_and_h")
  rg = [x for x in h.only_chain().residue_groups() if x.resseq == "   2"][0]
  hd = dict((ag.altloc, sorted(a.name.strip() for a in ag.atoms()
             if a.element.strip() in ("H", "D"))) for ag in rg.atom_groups())
  assert hd["A"] == ["DG"] and hd["B"] == ["DG"], hd
  assert hd["C"] == ["H"] and hd["D"] == ["D"], hd
  assert sorted(x.reason for x in r.conformer_sites) == ["conformer"] * 2

pdb_pep_og_ab = pdb_pep.replace(
  "ATOM     11  OG  SER A   2      15.456   0.826   7.311  1.00",
  "ATOM     11  OG ASER A   2      15.456   0.826   7.311  0.50").replace(
  "ATOM     12  N   LYS",
  "ATOM     11  OG BSER A   2      16.300   0.000   7.300  0.50 20.00"
  "           O\nATOM     12  N   LYS")

def ser_h_d(h):
  rg = [x for x in h.only_chain().residue_groups() if x.resseq == "   2"][0]
  return dict((ag.altloc, sorted(a.name.strip() for a in ag.atoms()
    if a.element.strip() in ("H", "D"))) for ag in rg.atom_groups())

def test_merge_conformer_copies():
  # placement puts the amide H once into the blank group (blank N, CA and
  # Ala 1 C) and the other Ser H into A and B; HA (blank parent, coincident
  # copies) is merged into the blank group, HB2/HB3/HG differ
  model = placed(pdb_pep_og_ab)
  hd = ser_h_d(model.get_hierarchy())
  assert hd[""] == ["H"] and hd["A"] == hd["B"] == ["HA", "HB2", "HB3", "HG"]
  h, r = run_option(model, "hd_and_h")
  hd = ser_h_d(h)
  assert hd[""] == ["HA"], hd                      # merged, not exchangeable
  assert hd["C"] == ["H"] and hd["D"] == ["D"], hd # merged and split
  assert hd["A"] == hd["B"] == ["DG", "HB2", "HB3"], hd
  assert sorted(names(x.site) for x in r.merged_sites) == ["HA"]
  # exchangeable output sites (not water): Ala H1 H2 H3, Ser H once, HG in A
  # and in B, Lys H HZ1 HZ2 HZ3
  assert r.n_exchangeable == 10, r.n_exchangeable
  assert all(x.reason == "merged conformer copies" for x in r.merged_sites)
  assert all(len(x.copies) == 2 for x in r.merged_sites)
  assert sorted((names(x.site), x.reason) for x in r.conformer_sites) == [
    ("HG", "conformer"), ("HG", "conformer")]
  # merged H at the parent occupancy; the split pair 0.5 each
  assert find(h, "SER", "HA").occ == 1.0
  assert find(h, "SER", "H").occ == 0.5 and find(h, "SER", "D").occ == 0.5
  # the split pair maps to the blank H
  src = [r.provenance[hd_exchange.atom_key(find(h, "SER", n))]
         for n in ("H", "D")]
  assert src[0] == src[1] and src[0][-3:] == ("", " H  ", "H"), src

def test_conformer_copies_differ():
  # a backbone H whose copies are 0.1 A apart is not merged: conformer site
  # placement gives one blank amide H here: copy it into A and B
  lines = placed(pdb_pep_og_ab).model_as_pdb().splitlines()
  h = [l for l in lines if l[12:20] == " H   SER"]
  assert len(h) == 1, h
  lines = sum([[l[:16] + alt + l[17:54] + "  0.50" + l[60:]
                for alt in "AB"] if l == h[0] else [l] for l in lines], [])
  model = make_model("\n".join(lines))
  copies = [a for a in model.get_hierarchy().atoms()
            if a.parent().resname == "SER" and a.name == " H  "]
  assert sorted(a.parent().altloc for a in copies) == ["A", "B"]
  b = [a for a in copies if a.parent().altloc == "B"][0]
  b.set_xyz((b.xyz[0] + 0.1, b.xyz[1], b.xyz[2]))
  h, r = run_option(model, "hd_and_h")
  hd = ser_h_d(h)
  assert "C" not in hd and hd["A"] == hd["B"] == ["D", "DG", "HB2", "HB3"], hd
  assert [names(x.site) for x in r.merged_sites] == ["HA"]
  assert sorted(names(x.site) for x in r.conformer_sites) == [
    "H", "H", "HG", "HG"]

def test_no_free_altloc_letters():
  # 12: all letters but one in use: no site of that residue is split
  lines = ["ATOM      1  N   SER A   2      14.767   2.090   9.884  1.00 20.00"
           "           N",
           "ATOM      2  OG  SER A   2      15.456   0.826   7.311  1.00 20.00"
           "           O",
           "ATOM      3  H   SER A   2      14.965   2.892   9.644  1.00 20.00"
           "           H",
           "ATOM      4  HG  SER A   2      14.778   1.319   7.357  1.00 20.00"
           "           H",
           "ATOM      5  CA  SER A   2      15.853   1.130   9.701  1.00 20.00"
           "           C",
           "ATOM      6  HA  SER A   2      15.510   0.233   9.835  1.00 20.00"
           "           H"]
  for i, alt in enumerate("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxy"):
    lines.append("ATOM  %5d  CB %sSER A   2      16.411   1.224   8.278  0.02"
                 " 20.00           C" % (i + 7, alt))
  h = iotbx.pdb.input(source_info=None, lines=lines).construct_hierarchy()
  keys = hd_exchange.atom_keys(h)
  table = {}
  for key in keys:
    if key[-1] != "H": continue
    parent = {" H  ": "N", " HG ": "OG", " HA ": "CA"}[key[-2]]
    table[key] = group_args(status="ok", exchangeable=parent != "CA",
      parent_occ=1.0, parent_altloc="",
      parent_id='pdb=" %-3s SER A   2 "' % parent)
  hc = h.deep_copy()
  r = hd_exchange.apply_neutron_option(hc, "hd_and_h", table)
  assert r.n_split == 0 and r.n_to_d == 2
  assert sorted(x.reason for x in r.conformer_sites) == \
    ["no free altloc letters"] * 2
  assert sorted(a.name.strip() for a in hc.atoms()
                if a.element.strip() == "D") == ["D", "DG"]
  # all_hd: written as conformer sites are; exchangeable -> D, the
  # non-exchangeable C-H -> a single H in its own (blank) group, no D
  hc = h.deep_copy()
  r = hd_exchange.apply_neutron_option(hc, "all_hd", table)
  assert r.n_split == 0
  ha = [a for a in hc.atoms() if a.name.strip() in ("HA", "DA")]
  assert [(a.name.strip(), a.element.strip(), a.parent().altloc)
          for a in ha] == [("HA", "H", "")], ha
  # sites are listed by their input (H) name
  assert sorted((names(x.site), x.reason) for x in r.conformer_sites) == [
    ("H", "no free altloc letters"), ("HA", "no free altloc letters"),
    ("HG", "no free altloc letters")]
  assert sorted(a.name.strip() for a in hc.atoms()
                if a.element.strip() == "D") == ["D", "DG"]

def test_name_not_pairable():
  # 1HO2 (exchangeable) and 1HC1 (not): never split or D; HO3 follows the
  # mode
  h = iotbx.pdb.input(source_info=None, lines="""\
HETATM    1  O2  LIG A   1       1.000   1.000   1.000  1.00 20.00           O
HETATM    2 1HO2 LIG A   1       1.900   1.000   1.000  1.00 20.00           H
HETATM    3  O3  LIG A   1       4.000   1.000   1.000  1.00 20.00           O
HETATM    4  HO3 LIG A   1       4.900   1.000   1.000  1.00 20.00           H
HETATM    5  C1  LIG A   1       7.000   1.000   1.000  1.00 20.00           C
HETATM    6 1HC1 LIG A   1       8.000   1.000   1.000  1.00 20.00           H
""".split("\n")).construct_hierarchy()
  def table_for(h):
    return dict((key, group_args(status="ok", exchangeable=key[-2] != "1HC1",
      parent_occ=1.0, parent_altloc="", parent_id="O")) for key in
      hd_exchange.atom_keys(h) if key[-1] == "H")
  # an H whose name does not start with H is never split and never D, in
  # any mode; listed wherever the mode would split it or write D
  for option, listed, ho3 in (
      ("hd_and_h", ["1HO2"], "split"), ("hd_and_d", ["1HC1", "1HO2"], "split"),
      ("all_hd", ["1HC1", "1HO2"], "split"), ("all_d", ["1HC1", "1HO2"], "D"),
      ("all_h", [], "H")):
    hc = h.deep_copy()
    r = hd_exchange.apply_neutron_option(hc, option, table_for(h))
    digit = sorted((a.name, a.element.strip(), a.parent().altloc)
                   for a in hc.atoms() if a.name[0] == "1")
    assert digit == [("1HC1", "H", ""), ("1HO2", "H", "")], (option, digit)
    assert sorted(x.site[5:9] for x in r.unclassified_sites) == listed, option
    for x in r.unclassified_sites:
      assert x.reason == hd_exchange.NOT_PAIRABLE
      assert "kept as H" in x.reason
      # written as it was: output label equals the input label
      assert x.output == x.site.replace("pdb=", "").replace('"', "")
    o3 = sorted((a.name.strip(), a.parent().altloc) for a in hc.atoms()
                if a.name.strip() in ("HO3", "DO3"))
    pairs = hc.exchangeable_hd_selections()
    if ho3 == "split":
      assert o3 == [("DO3", "B"), ("HO3", "A")], (option, o3)
      assert r.n_split == 1 and len(pairs) == 1
    else:
      assert o3 == [("%sO3" % ho3, "")], (option, o3)
      assert r.n_split == 0 and pairs == []

def test_parent_occupancy_from_table():
  # 10: H written at 1.0 on a parent at 0.6 -> 0.3 / 0.3, partial parent
  pdb = pdb_pep.replace(
    "OG  SER A   2      15.456   0.826   7.311  1.00",
    "OG  SER A   2      15.456   0.826   7.311  0.60")
  model = placed(pdb)
  hg = find(model.get_hierarchy(), "SER", "HG")
  hg.set_occ(1.0)
  h, r = run_option(model, "hd_and_h")
  assert values(find(h, "SER", "HG"))[0] == 0.3
  assert values(find(h, "SER", "DG"))[0] == 0.3
  assert [x.parent_occ for x in r.partial_parent_sites] == [0.6]

def test_anisotropic():
  # 13: split keeps uij on H and copies it to D
  uij = (0.20, 0.21, 0.22, 0.01, 0.02, 0.03)
  model = pep_h()
  hg = find(model.get_hierarchy(), "SER", "HG")
  hg.set_uij(uij)
  try:
    h, r = run_option(model, "hd_and_h")
  finally:
    hg.uij_erase()
  for name in ("HG", "DG"):
    assert approx_equal(find(h, "SER", name).uij, uij)
  # record/restore keeps uij of both exactly
  uij_d = (0.30, 0.31, 0.32, 0.04, 0.05, 0.06)
  def edit(h):
    find(h, "SER", "HG").set_uij(uij)
    find(h, "SER", "DG").set_uij(uij_d)
  model, records = normalize(hd_input("hd_and_h", edit))
  h, r = run_option(model, "hd_and_h", records)
  assert approx_equal(find(h, "SER", "HG").uij, uij, eps=1e-4)
  assert approx_equal(find(h, "SER", "DG").uij, uij_d, eps=1e-4)

def set_pair(h, resname, h_name, d_name, h_vals, d_vals):
  for name, (occ, b) in ((h_name, h_vals), (d_name, d_vals)):
    a = find(h, resname, name)
    a.set_occ(occ)
    a.set_b(b)

def test_round_trip():
  # 8: normalize_to_h() + apply_neutron_option() restore occupancy and B
  def edit(h):
    set_pair(h, "SER", "HG", "DG", (0.3, 15.0), (0.7, 25.0))
    set_pair(h, "LYS", "HZ1", "DZ1", (0.45, 18.0), (0.55, 22.0))
  model, records = normalize(hd_input("hd_and_h", edit))
  assert len(records) == EXCHANGEABLE_PEP
  # normalized: H at the parent occupancy, in the blank group
  hg = find(model.get_hierarchy(), "SER", "HG")
  assert hg.parent().altloc == "" and hg.occ == 1.0
  for option in ("hd_and_h", "hd_and_d", "all_hd"):
    h, r = run_option(model, option, records)
    assert values(find(h, "SER", "HG"))[:2] == (0.3, 15.0)
    assert values(find(h, "SER", "DG"))[:2] == (0.7, 25.0)
    assert values(find(h, "LYS", "HZ1"))[:2] == (0.45, 18.0)
    assert values(find(h, "LYS", "DZ1"))[:2] == (0.55, 22.0)
    assert r.n_restored >= EXCHANGEABLE_PEP

def test_30f_split():
  # 8a: Se-H (1.444 A) is split and the output check accepts it
  model = make_model(pdb_30f)
  h, r = run_option(model, "hd_and_h")
  assert r.n_split == 3   # H2 on Se, H1 and H3 on N
  assert find(h, "30F", "H2").parent().altloc == "A"
  assert find(h, "30F", "D2").parent().altloc == "B"

def test_alternative_orientations_output():
  # 18b: Ser HG in A and another HG in B (both H): single isotope with their
  # altlocs, listed; the backbone H of the residue is split into C/D
  h = pep_h().get_hierarchy().deep_copy()
  hg = find(h, "SER", "HG")
  rg = hg.parent().parent()
  for altloc, xyz in (("A", hg.xyz), ("B", (15.300, 0.200, 6.600))):
    ag = iotbx.pdb.hierarchy.atom_group(altloc=altloc, resname="SER")
    a = hg.detached_copy()
    a.set_xyz(xyz)
    a.set_occ(0.5)
    ag.append_atom(a)
    rg.append_atom_group(ag)
  hg.parent().remove_atom(hg)
  model = model_from_hierarchy(h, pep_h().crystal_symmetry())
  h, r = run_option(model, "hd_and_h")
  dg = sorted(a.parent().altloc for a in h.atoms() if a.name.strip() == "DG")
  assert dg == ["A", "B"], dg
  assert [x.reason for x in r.conformer_sites] == ["conformer"] * 2
  rg = [x for x in h.only_chain().residue_groups() if x.resseq == "   2"][0]
  hd = dict((ag.altloc, sorted(a.name.strip() for a in ag.atoms()
             if a.element.strip() in ("H", "D"))) for ag in rg.atom_groups())
  assert hd["C"] == ["H"] and hd["D"] == ["D"], hd

def test_cross_mode():
  # 19, fully occupied parent, input pair H 0.3 / D 0.7, different B and uij
  uij_h = (0.20, 0.21, 0.22, 0.01, 0.02, 0.03)
  uij_d = (0.30, 0.31, 0.32, 0.04, 0.05, 0.06)
  def edit(h):
    set_pair(h, "SER", "HG", "DG", (0.3, 15.0), (0.7, 25.0))
    find(h, "SER", "HG").set_uij(uij_h)
    find(h, "SER", "DG").set_uij(uij_d)
  model, records = normalize(hd_input("hd_and_h", edit))
  # with ANISOU, B is the isotropic equivalent of uij: compare with records
  rec = records[hd_exchange.site_key(find(model.get_hierarchy(), "SER", "HG"),
                                     "")]
  assert rec["H"][0] == 0.3 and rec["D"][0] == 0.7
  assert abs(rec["H"][1] - rec["D"][1]) > 1
  assert approx_equal(rec["H"][2], uij_h) and approx_equal(rec["D"][2], uij_d)
  def same(atom, occ, record):
    assert approx_equal(atom.occ, occ) and approx_equal(atom.b, record[1])
    assert approx_equal(atom.uij, record[2])
  for option in ("hd_and_h", "hd_and_d", "all_hd"):
    h, r = run_option(model, option, records)
    same(find(h, "SER", "HG"), 0.3, rec["H"])
    same(find(h, "SER", "DG"), 0.7, rec["D"])
  h, r = run_option(model, "all_d", records)
  dg = find(h, "SER", "DG")
  same(dg, 1.0, rec["D"])
  assert dg.parent().altloc == ""
  h, r = run_option(model, "all_h", records)
  hg = find(h, "SER", "HG")
  same(hg, 1.0, rec["H"])
  assert hg.parent().altloc == ""
  # parent at 0.6, input sum 1: a single isotope gets the parent's 0.6
  pdb = hd_input("hd_and_h", edit).replace(
    "OG  SER A   2      15.456   0.826   7.311  1.00",
    "OG  SER A   2      15.456   0.826   7.311  0.60")
  model6, records6 = normalize(pdb)
  for option, name in (("all_d", "DG"), ("all_h", "HG")):
    h, r = run_option(model6, option, records6)
    assert values(find(h, "SER", name))[0] == 0.6
  # input all_d (single D, no record) -> hd_and_h: new pair from the parent,
  # both with the former D's B and uij
  def edit_d(h):
    dg = find(h, "SER", "DG")
    dg.set_b(33.0)
    dg.set_uij(uij_d)
  model_d, records_d = normalize(hd_input("all_d", edit_d))
  assert records_d == {}
  former_d = find(model_d.get_hierarchy(), "SER", "HG")   # renamed to H
  b_d = former_d.b
  h, r = run_option(model_d, "hd_and_h", records_d)
  for name in ("HG", "DG"):
    a = find(h, "SER", name)
    assert approx_equal(a.occ, 0.5) and approx_equal(a.b, b_d)
    assert approx_equal(a.uij, uij_d, eps=1e-4)
  # single D at 0.7 on a parent at 1.0, no record: kept and noted as one
  # isotope; split from the parent
  def edit_d7(h):
    find(h, "SER", "DG").set_occ(0.7)
  model7, records7 = normalize(hd_input("all_d", edit_d7))
  for option, name in (("all_d", "DG"), ("all_h", "HG")):
    h, r = run_option(model7, option, records7)
    assert values(find(h, "SER", name))[0] == 0.7
    assert [x.site[5:9].strip() for x in r.notes] == [name]
  h, r = run_option(model7, "hd_and_h", records7)
  assert values(find(h, "SER", "HG"))[0] == 0.5
  assert values(find(h, "SER", "DG"))[0] == 0.5

pdb_water_hd = """\
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1
HETATM    1  O   %(res)s A   1      10.000  10.000  10.000  0.80 20.00           O
HETATM    2  H1 A%(res)s A   1      10.850  10.000  10.000  0.30 14.00           H
HETATM    3  H2 A%(res)s A   1       9.795  10.825  10.000  0.30 14.00           H
HETATM    4  D1 B%(res)s A   1      10.850  10.000  10.000  0.70 26.00           D
HETATM    5  D2 B%(res)s A   1       9.795  10.825  10.000  0.70 26.00           D
END
"""

def test_water_cross_mode():
  # 19: an input water H/D pair becomes one isotope at the O's occupancy
  # with that isotope's B: DOD with D in every neutron mode except all_h,
  # HOH with H in all_h; never a pair. Input named HOH, DOD and WAT.
  for res in ("HOH", "DOD", "WAT"):
    model, records = normalize(pdb_water_hd % {"res": res})
    assert len(records) == 2
    w = [ag for ag in model.get_hierarchy().atom_groups()]
    assert len(w) == 1 and w[0].altloc == ""
    assert w[0].resname == ("HOH" if res == "DOD" else res)
    for option in ("all_d", "hd_and_h", "hd_and_d", "all_hd", "all_h"):
      h, r = run_option(model, option, records)
      assert r.n_split == 0 and r.n_restored == 2
      ags = list(h.atom_groups())
      assert len(ags) == 1 and ags[0].altloc == ""
      element = "H" if option == "all_h" else "D"
      expected = "DOD" if element == "D" else (
        "HOH" if res == "DOD" else res)
      assert ags[0].resname == expected, (res, option, ags[0].resname)
      for a in ags[0].atoms():
        if a.name.strip() == "O": continue
        assert a.element.strip() == element
        assert values(a)[:2] == (0.8, 14.0 if element == "H" else 26.0)
  # input HOH (H only) -> DOD
  h, r = run_option(pep_h(), "hd_and_h")
  assert "DOD" in [ag.resname for ag in h.atom_groups()]

pdb_water_alt = """\
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1
HETATM    1  O   %(res)s A   1      10.000  10.000  10.000  1.00 20.00           O
HETATM    2 %(x)s1 A%(res)s A   1      10.850  10.000  10.000  0.50 20.00          %(x)s
HETATM    3 %(x)s2 A%(res)s A   1       9.795  10.825  10.000  0.50 20.00          %(x)s
HETATM    4 %(x)s1 B%(res)s A   1       9.150  10.000  10.000  0.50 20.00          %(x)s
HETATM    5 %(x)s2 B%(res)s A   1      10.205   9.175  10.000  0.50 20.00          %(x)s
END
"""

def check_water_reinterpreted(h, crystal_symmetry, element):
  """fresh restraints: every H/D bonded to the O, nothing floating"""
  assert set(ag.resname for ag in h.atom_groups()) == set(
    ["DOD" if element == "D" else "HOH"])
  model = model_from_hierarchy(h, crystal_symmetry)
  table = hd_exchange.h_parent_table(model)
  assert len(table) == 4
  for e in table.values():
    assert e.status == "ok" and " O  " in e.parent_id, (e.h_id, e.status)
    assert e.h_id[5:7] == " " + element

def test_water_alternative_orientations():
  # 19: blank O, alternative H1/H2 in A and B: in hd_and_h every atom_group
  # of the residue_group is DOD, and reinterpretation bonds every D to the O
  model = model_from_str(pdb_water_alt % {"res": "HOH", "x": " H"})
  h, r = run_option(model, "hd_and_h")
  assert r.n_water == 1 and r.n_split == 0
  check_water_reinterpreted(h, model.crystal_symmetry(), "D")
  # the normalize_to_h() direction: DOD with alternative D orientations
  model, records = normalize(pdb_water_alt % {"res": "DOD", "x": " D"})
  assert records == {}
  check_water_reinterpreted(model.get_hierarchy(), model.crystal_symmetry(),
                            "H")

# an O in each altloc, each with its own H
pdb_water_alt_o = """\
HETATM    6  O  AHOH A   2       5.000   5.000   5.000  0.60 20.00           O
HETATM    7  H1 AHOH A   2       5.850   5.000   5.000  0.60 20.00           H
HETATM    8  H2 AHOH A   2       4.795   5.825   5.000  0.60 20.00           H
HETATM    9  O  BHOH A   2       5.300   5.000   5.000  0.40 20.00           O
HETATM   10  H1 BHOH A   2       6.150   5.000   5.000  0.40 20.00           H
HETATM   11  H2 BHOH A   2       5.095   5.825   5.000  0.40 20.00           H
END
"""

def test_add_water_h_parents():
  """
  add_water_h_parents() gives water H the entries h_parent_table() gives
  them from restraints, for a blank O with H in altlocs A and B and for an O
  in each altloc; source names the O. Entries already in the table stay.
  """
  pdb = (pdb_water_alt % {"res": "HOH", "x": " H"}).replace(
    "END\n", pdb_water_alt_o)
  model = model_from_str(pdb)
  expected = hd_exchange.h_parent_table(model)
  assert len(expected) == 8, len(expected)
  table = {}
  n = hd_exchange.add_water_h_parents(table, model.get_hierarchy(), 0.85)
  assert n == 8, n
  assert sorted(table) == sorted(expected), set(table) ^ set(expected)
  for key, e in expected.items():
    for field in ("h_id", "parent_id", "parent_element", "parent_occ",
                  "parent_altloc", "exchangeable", "status", "parent_key",
                  "parent_sym_ops"):
      got, want = getattr(table[key], field), getattr(e, field)
      assert got == want, (key, field, got, want)
    source = table[key].source
    assert len(source) == 1 and len(source[0]) == 4, source
    assert source[0][0] == e.parent_id and source[0][3] == "x,y,z", source
  key = sorted(expected)[0]
  kept = group_args(status="kept")
  table = {key: kept}
  n = hd_exchange.add_water_h_parents(table, model.get_hierarchy(), 0.85)
  assert n == 7 and table[key] is kept, n

def test_lone_alternate_h():
  # blank OG, HG only in A at 0.30: never merged or split; one isotope in A
  pdb = pdb_ser_hd.replace(
    "ATOM     10  HG ASER A   2      15.856   0.874   6.441  0.50",
    "ATOM     10  HG ASER A   2      15.856   0.874   6.441  0.30").replace(
    "ATOM     11  DG BSER A   2      15.856   0.874   6.441  0.50 20.00"
    "           D\n", "")
  model = model_from_str(pdb)
  h, r = run_option(model, "all_h")
  hg = find(h, "SER", "HG")
  assert hg.parent().altloc == "A" and values(hg)[0] == 0.3
  assert r.merged_sites == []
  assert [(names(x.site), x.reason) for x in r.notes] == [
    ("HG", "occupancy differs from parent")]
  h, r = run_option(model, "hd_and_h")
  assert r.n_split == 0 and r.merged_sites == []
  dg = find(h, "SER", "DG")
  assert dg.parent().altloc == "A" and values(dg)[0] == 0.3
  assert [(names(x.site), x.reason) for x in r.conformer_sites] == [
    ("HG", "conformer")]
  assert not [a for a in h.atoms() if a.name.strip() == "HG"]

pdb_two_models = """\
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1
MODEL        1
%s
ENDMDL
MODEL        2
%s
ENDMDL
END
"""

def water_block(chain, resseq, h_vals, d_vals, x=10.0):
  lines = ["HETATM    1  O   HOH %s%4d    %8.3f  10.000  10.000  1.00 20.00"
           "           O" % (chain, resseq, x)]
  for name, alt, element, (occ, b), dx, dy in (
      ("H1", "A", "H", h_vals, 0.850, 0.0), ("H2", "A", "H", h_vals, -0.205,
       0.825), ("D1", "B", "D", d_vals, 0.850, 0.0), ("D2", "B", "D", d_vals,
       -0.205, 0.825)):
    lines.append("HETATM    1  %-3s%sHOH %s%4d    %8.3f%8.3f  10.000%6.2f%6.2f"
      "           %s" % (name, alt, chain, resseq, x + dx, 10.0 + dy, occ, b,
      element))
  return "\n".join(lines)

def test_two_models():
  # 16: identical labels in two models, different values: each keeps its own
  pdb = pdb_two_models % (
    water_block("A", 1, (0.3, 11.0), (0.7, 21.0)),
    water_block("A", 1, (0.8, 12.0), (0.2, 22.0)))
  model, records = normalize(pdb)
  assert len(records) == 4
  for option, b in (("hd_and_h", (21.0, 22.0)), ("all_h", (11.0, 12.0))):
    h, r = run_option(model, option, records)
    for m, b_expected in zip(h.models(), b):
      for a in m.atoms():
        if a.name.strip() != "O":
          assert values(a)[:2] == (1.0, b_expected), (option, values(a))

def hand_table(hierarchy):
  """parent table for waters by hand: every H/D on the O of its residue"""
  table = {}
  for key, atom in zip(hd_exchange.atom_keys(hierarchy), hierarchy.atoms()):
    if key[-1] not in ("H", "D"): continue
    rg = atom.parent().parent()
    o = [a for a in rg.atoms() if a.name.strip() == "O"][0]
    table[key] = group_args(status="ok", exchangeable=True, parent_occ=o.occ,
      parent_altloc="", parent_key=hd_exchange.atom_key(o),
      parent_sym_ops=("x,y,z",), parent_id=o.id_str(), source=[])
  return table

def test_repeated_labels_restore():
  # 17: repeated chain id (two chain A blocks) and a duplicated residue
  # number within one chain: each site keeps its own values.
  # pdb_interpretation rejects duplicate labels, so the tables are by hand.
  pdb = "\n".join([
    water_block("A", 1, (0.3, 11.0), (0.7, 21.0), x=5.0),
    water_block("A", 2, (0.3, 12.0), (0.7, 22.0), x=10.0),
    water_block("A", 1, (0.3, 13.0), (0.7, 23.0), x=15.0),
    water_block("B", 1, (0.3, 14.0), (0.7, 24.0), x=20.0),
    water_block("A", 1, (0.3, 15.0), (0.7, 25.0), x=25.0)])
  h = iotbx.pdb.input(source_info=None, lines=pdb.split("\n")
    ).construct_hierarchy()
  assert [c.id for c in h.only_model().chains()] == ["A", "B", "A"]
  match = hd_exchange.collapse_check(
    hd_exchange.find_hd_candidates(h, hand_table(h)))
  assert match.problems == [] and len(match.pairs) == 10
  records = hd_exchange.normalize_to_h(h, match).records
  assert len(records) == 10
  table = hand_table(h)
  r = hd_exchange.apply_neutron_option(h, "hd_and_h", table, records=records)
  assert r.n_restored == 10
  b = [round(a.b, 1) for a in h.atoms() if a.name.strip() == "D1"]
  assert b == [21.0, 22.0, 23.0, 24.0, 25.0], b
  c = hd_exchange.find_hd_candidates(h, table, provenance=r.provenance)
  k = hd_exchange.collapse_check(c)
  assert k.problems == [] and k.unresolved == []

def test_water_names_follow_content():
  # DOD with H -> HOH in all_h (every atom_group), DOD with D in hd_and_h;
  # DOD with only O and WAT with H keep their names
  def water(resname, h_altloc=" ", h=True):
    lines = ["HETATM    1  O   %s A   1      10.000  10.000  10.000  1.00 20.00"
             "           O" % resname]
    if h:
      for name, x, y in (("H1", 10.850, 10.0), ("H2", 9.795, 10.825)):
        lines.append("HETATM    1  %-3s%s%s A   1    %8.3f%8.3f  10.000  1.00"
          " 20.00           H" % (name, h_altloc, resname, x, y))
    return iotbx.pdb.input(source_info=None, lines=lines).construct_hierarchy()
  for label, h, option, names, element in (
      ("DOD+H", water("DOD"), "all_h", ["HOH"], "H"),
      ("DOD+H", water("DOD"), "hd_and_h", ["DOD"], "D"),
      ("DOD+H in A", water("DOD", "A"), "all_h", ["HOH", "HOH"], "H"),
      ("DOD+H in A", water("DOD", "A"), "hd_and_h", ["DOD", "DOD"], "D"),
      ("DOD O only", water("DOD", h=False), "all_h", ["DOD"], None),
      ("DOD O only", water("DOD", h=False), "hd_and_h", ["DOD"], None),
      ("WAT+H", water("WAT"), "all_h", ["WAT"], "H")):
    r = hd_exchange.apply_neutron_option(h, option, hand_table(h))
    assert [ag.resname for ag in h.atom_groups()] == names, (label, option,
      [ag.resname for ag in h.atom_groups()])
    elements = set(a.element.strip() for a in h.atoms()) - set(["O"])
    assert elements == (set([element]) if element else set()), (label, option)
    assert r.n_water == (1 if element == "D" else 0)
    assert r.n_water_to_hoh == (1 if label.startswith("DOD+H") and
                                option == "all_h" else 0), (label, option)

def test_conformer_copy_occupancy_note():
  # H copied into each conformer (blank parent CB) at the conformer's
  # occupancy: no note; a copy that differs from its conformer: one note
  model = placed(pdb_pep_og_ab)
  h, r = run_option(model, "hd_and_h")
  hb = [a for a in h.atoms() if a.parent().resname == "SER" and
        a.name.strip() in ("HB2", "HB3")]
  assert sorted((a.name.strip(), a.parent().altloc, a.occ) for a in hb) == [
    ("HB2", "A", 0.5), ("HB2", "B", 0.5), ("HB3", "A", 0.5),
    ("HB3", "B", 0.5)]
  assert r.notes == [], [(x.output, x.text) for x in r.notes]
  hb2 = [a for a in model.get_hierarchy().atoms() if a.parent().resname ==
         "SER" and a.name.strip() == "HB2" and a.parent().altloc == "A"][0]
  hb2.set_occ(0.3)
  try:
    h, r = run_option(model, "hd_and_h")
  finally:
    hb2.set_occ(0.5)
  assert [(names(x.site), x.reason, x.text) for x in r.notes] == [
    ("HB2", "occupancy differs from conformer",
     "occupancy 0.30, conformer 0.50")]

# B16 repair fixtures, cut from deposited models (atoms with their H/D as
# deposited; refused by reduce2 before the repairs)

# (a) 7T5E A 139-141: blank HA on CA A/B (copies 0.006 A)
pdb_a_7t5e_glu_a140 = """\
CRYST1   68.296   42.269   70.413  90.00  98.47  90.00 P 1 21 1
ATOM      1  N   PRO A 139     -10.887 -13.427 -40.812  1.00 20.62           N
ATOM      2  CA  PRO A 139      -9.986 -13.938 -41.852  1.00 18.23           C
ATOM      3  C   PRO A 139     -10.737 -14.683 -42.956  1.00 16.15           C
ATOM      4  O   PRO A 139     -11.730 -15.370 -42.711  1.00 15.86           O
ATOM      5  CB  PRO A 139      -9.066 -14.897 -41.080  1.00 23.92           C
ATOM      6  CG  PRO A 139      -9.195 -14.559 -39.732  1.00 19.78           C
ATOM      7  CD  PRO A 139     -10.552 -14.028 -39.509  1.00 19.48           C
ATOM      8  HA  PRO A 139      -9.382 -13.127 -42.271  1.00 23.23           H
ATOM      9  HB2 PRO A 139      -9.378 -15.930 -41.229  1.00 30.07           H
ATOM     10  HB3 PRO A 139      -8.035 -14.724 -41.403  1.00 24.85           H
ATOM     11  HG2 PRO A 139      -9.025 -15.456 -39.153  1.00 19.95           H
ATOM     12  HG3 PRO A 139      -8.441 -13.818 -39.485  1.00 19.40           H
ATOM     13  HD2 PRO A 139     -11.231 -14.832 -39.265  1.00 17.27           H
ATOM     14  HD3 PRO A 139     -10.547 -13.289 -38.717  1.00 17.21           H
ATOM     15  N   GLU A 140     -10.221 -14.536 -44.178  1.00 18.62           N
ATOM     16  C   GLU A 140     -10.444 -16.633 -45.486  1.00 24.06           C
ATOM     17  O   GLU A 140     -11.187 -17.379 -46.148  1.00 26.07           O
ATOM     18  HA  GLU A 140     -11.880 -15.061 -45.321  1.00 25.82           H
ATOM     19  CA AGLU A 140     -10.796 -15.153 -45.362  0.46 23.38           C
ATOM     20  CB AGLU A 140     -10.316 -14.402 -46.611  0.46 30.70           C
ATOM     21  CG AGLU A 140      -8.891 -14.723 -47.041  0.46 31.17           C
ATOM     22  CD AGLU A 140      -8.054 -13.488 -47.367  0.46 48.28           C
ATOM     23  OE1AGLU A 140      -6.841 -13.493 -47.043  0.46 31.57           O
ATOM     24  OE2AGLU A 140      -8.600 -12.511 -47.923  0.46 47.29           O
ATOM     25  D  AGLU A 140      -9.388 -13.992 -44.381  0.46 36.01           D
ATOM     26  CA BGLU A 140     -10.792 -15.154 -45.366  0.54 23.34           C
ATOM     27  CB BGLU A 140     -10.295 -14.415 -46.618  0.54 30.68           C
ATOM     28  CG BGLU A 140     -10.873 -14.921 -47.948  0.54 40.69           C
ATOM     29  CD BGLU A 140     -10.689 -13.924 -49.098  0.54 50.00           C
ATOM     30  OE1BGLU A 140      -9.674 -13.187 -49.112  0.54 46.66           O
ATOM     31  OE2BGLU A 140     -11.573 -13.866 -49.982  0.54 56.50           O
ATOM     32  D  BGLU A 140      -9.389 -13.994 -44.381  0.54 36.11           D
ATOM     33  N   CYS A 141      -9.352 -17.071 -44.862  1.00 18.75           N
ATOM     34  CA  CYS A 141      -8.802 -18.402 -45.065  1.00 23.36           C
ATOM     35  C   CYS A 141      -9.215 -19.418 -44.021  1.00 24.23           C
ATOM     36  O   CYS A 141      -8.716 -20.547 -44.055  1.00 21.50           O
ATOM     37  CB  CYS A 141      -7.272 -18.332 -45.115  1.00 26.50           C
ATOM     38  SG  CYS A 141      -6.473 -17.757 -43.539  1.00 21.70           S
ATOM     39  HA  CYS A 141      -9.162 -18.738 -46.030  1.00 21.61           H
ATOM     40  HB2 CYS A 141      -6.875 -19.325 -45.345  1.00 31.54           H
ATOM     41  HB3 CYS A 141      -7.017 -17.619 -45.902  1.00 21.32           H
ATOM     42  D  BCYS A 141      -8.828 -16.516 -44.198  1.00 25.27           D
END
"""

# (a) 7T5E B 209-211: blank HA on CA A/C, an H/D-only group B
pdb_a_7t5e_val_b210 = """\
CRYST1   68.296   42.269   70.413  90.00  98.47  90.00 P 1 21 1
ATOM      1  N   PRO B 209       7.559  -1.211   1.234  1.00 25.58           N
ATOM      2  CA  PRO B 209       6.584  -1.728   2.200  1.00 31.86           C
ATOM      3  C   PRO B 209       6.183  -0.661   3.204  1.00 31.10           C
ATOM      4  O   PRO B 209       6.997   0.127   3.670  1.00 28.33           O
ATOM      5  CB  PRO B 209       7.317  -2.900   2.872  1.00 27.20           C
ATOM      6  CG  PRO B 209       8.510  -3.177   2.072  1.00 32.84           C
ATOM      7  CD  PRO B 209       8.847  -1.917   1.325  1.00 43.87           C
ATOM      8  HA  PRO B 209       5.710  -2.105   1.666  1.00 29.73           H
ATOM      9  HB2 PRO B 209       7.589  -2.631   3.896  1.00 35.05           H
ATOM     10  HB3 PRO B 209       6.652  -3.766   2.849  1.00 30.17           H
ATOM     11  HG2 PRO B 209       9.337  -3.434   2.740  1.00 32.59           H
ATOM     12  HG3 PRO B 209       8.284  -3.988   1.374  1.00 32.77           H
ATOM     13  HD2 PRO B 209       9.559  -1.307   1.878  1.00 35.57           H
ATOM     14  HD3 PRO B 209       9.230  -2.171   0.340  1.00 37.38           H
ATOM     15  N   VAL B 210       4.905  -0.625   3.496  1.00 24.63           N
ATOM     16  C   VAL B 210       4.441   0.020   5.792  1.00 24.93           C
ATOM     17  O   VAL B 210       4.215  -1.128   6.157  1.00 25.97           O
ATOM     18  HA  VAL B 210       4.849   1.350   4.166  1.00 31.78           H
ATOM     19  CA AVAL B 210       4.309   0.419   4.321  0.79 25.69           C
ATOM     20  CB AVAL B 210       2.845   0.659   3.874  0.79 25.16           C
ATOM     21  CG1AVAL B 210       1.911   1.043   4.959  0.79 23.77           C
ATOM     22  CG2AVAL B 210       2.858   1.728   2.812  0.79 33.60           C
ATOM     23  HB AVAL B 210       2.429  -0.245   3.425  0.79 23.75           H
ATOM     24 HG11AVAL B 210       2.434   1.646   5.704  0.79 25.73           H
ATOM     25 HG12AVAL B 210       1.093   1.609   4.513  0.79 24.77           H
ATOM     26 HG13AVAL B 210       1.527   0.128   5.410  0.79 25.88           H
ATOM     27 HG21AVAL B 210       1.838   1.877   2.470  0.79 21.17           H
ATOM     28 HG22AVAL B 210       3.277   2.650   3.223  0.79 24.37           H
ATOM     29 HG23AVAL B 210       3.489   1.378   1.992  0.79 25.53           H
ATOM     30  CA CVAL B 210       4.307   0.418   4.326  0.21 25.86           C
ATOM     31  CB CVAL B 210       2.831   0.655   3.933  0.21 25.25           C
ATOM     32  CG1CVAL B 210       2.115  -0.645   3.570  0.21 24.85           C
ATOM     33  CG2CVAL B 210       2.066   1.372   5.028  0.21 25.00           C
ATOM     34  HB CVAL B 210       2.804   1.297   3.053  0.21 27.44           H
ATOM     35 HG11CVAL B 210       2.778  -1.501   3.677  0.21 27.98           H
ATOM     36 HG12CVAL B 210       1.251  -0.776   4.219  0.21 25.33           H
ATOM     37 HG13CVAL B 210       1.776  -0.564   2.541  0.21 23.38           H
ATOM     38 HG21CVAL B 210       1.571   0.634   5.659  0.21 25.49           H
ATOM     39 HG22CVAL B 210       2.754   1.968   5.627  0.21 25.83           H
ATOM     40 HG23CVAL B 210       1.326   2.019   4.559  0.21 24.45           H
ATOM     41  D  BVAL B 210       4.205  -0.264   2.864  1.00 29.44           D
ATOM     42  N   THR B 211       4.893   0.956   6.626  1.00 25.79           N
ATOM     43  CA  THR B 211       4.963   0.728   8.073  1.00 21.85           C
ATOM     44  C   THR B 211       3.884   1.478   8.849  1.00 24.60           C
ATOM     45  O   THR B 211       3.597   1.117  10.001  1.00 28.61           O
ATOM     46  CB  THR B 211       6.351   1.144   8.611  1.00 28.38           C
ATOM     47  OG1 THR B 211       6.582   2.537   8.332  1.00 36.29           O
ATOM     48  CG2 THR B 211       7.444   0.302   7.956  1.00 38.77           C
ATOM     49  HA  THR B 211       4.835  -0.337   8.275  1.00 34.01           H
ATOM     50  HB  THR B 211       6.392   0.970   9.687  1.00 32.03           H
ATOM     51 HG21 THR B 211       8.371   0.872   7.874  1.00 48.05           H
ATOM     52 HG22 THR B 211       7.636  -0.585   8.567  1.00 28.85           H
ATOM     53 HG23 THR B 211       7.115  -0.020   6.968  1.00 27.33           H
ATOM     54  H  ATHR B 211       5.211   1.881   6.341  0.26 26.62           H
ATOM     55  HG1ATHR B 211       6.119   3.100   8.990  0.47 28.61           H
ATOM     56  D  BTHR B 211       5.252   1.853   6.297  0.74 26.61           D
ATOM     57  DG1BTHR B 211       6.134   3.085   9.006  0.53 28.62           D
END
"""

# (a) 6TAE A 12-14: blank D on N A/B, the only atom of its blank group
pdb_a_6tae_gln_a13 = """\
CRYST1   82.904   82.904   75.868  90.00  90.00  90.00 P 42 21 2
ATOM      1  N   TYR A  12     100.801  -4.172  27.264  1.00 25.27           N
ATOM      2  CA  TYR A  12     100.581  -2.844  26.673  1.00 22.65           C
ATOM      3  C   TYR A  12      99.761  -1.952  27.636  1.00 24.90           C
ATOM      4  O   TYR A  12     100.095  -0.788  27.859  1.00 23.01           O
ATOM      5  CB  TYR A  12      99.878  -2.985  25.320  1.00 21.33           C
ATOM      6  CG  TYR A  12     100.012  -1.781  24.449  1.00 20.56           C
ATOM      7  CD1 TYR A  12      99.309  -0.617  24.686  1.00 21.69           C
ATOM      8  CD2 TYR A  12     100.862  -1.799  23.385  1.00 20.24           C
ATOM      9  CE1 TYR A  12      99.494   0.523  23.872  1.00 21.51           C
ATOM     10  CE2 TYR A  12     101.052  -0.693  22.573  1.00 24.49           C
ATOM     11  CZ  TYR A  12     100.360   0.461  22.802  1.00 22.85           C
ATOM     12  OH  TYR A  12     100.533   1.565  21.976  1.00 24.00           O
ATOM     13  D   TYR A  12     100.486  -4.807  26.773  1.00 25.81           D
ATOM     14  DA  TYR A  12     101.450  -2.441  26.480  1.00 19.28           D
ATOM     15  DB2 TYR A  12     100.335  -3.696  24.849  1.00 16.09           D
ATOM     16  DB3 TYR A  12      98.957  -3.228  25.421  1.00 23.28           D
ATOM     17  DD1 TYR A  12      98.725  -0.584  25.407  1.00 29.19           D
ATOM     18  DD2 TYR A  12     101.348  -2.575  23.234  1.00 26.19           D
ATOM     19  DE1 TYR A  12      99.004   1.299  24.032  1.00 31.31           D
ATOM     20  DE2 TYR A  12     101.654  -0.754  21.875  1.00 27.01           D
ATOM     21  DH  TYR A  12     100.994   2.120  22.372  1.00 33.18           D
ATOM     22  D   GLN A  13      98.409  -3.302  28.054  1.00 28.51           D
ATOM     23  N  AGLN A  13      98.683  -2.508  28.189  0.58 24.72           N
ATOM     24  CA AGLN A  13      97.818  -1.802  29.118  0.58 24.21           C
ATOM     25  C  AGLN A  13      98.576  -1.343  30.370  0.58 22.95           C
ATOM     26  O  AGLN A  13      98.369  -0.234  30.835  0.58 26.94           O
ATOM     27  CB AGLN A  13      96.597  -2.686  29.435  0.58 26.05           C
ATOM     28  CG AGLN A  13      95.662  -2.841  28.227  0.58 28.35           C
ATOM     29  CD AGLN A  13      94.469  -3.786  28.438  0.58 40.85           C
ATOM     30  OE1AGLN A  13      94.372  -4.872  27.830  0.58 29.84           O
ATOM     31  NE2AGLN A  13      93.534  -3.345  29.251  0.58 36.97           N
ATOM     32  DA AGLN A  13      97.482  -0.992  28.687  0.58 27.36           D
ATOM     33  DB2AGLN A  13      96.914  -3.574  29.661  0.58 27.18           D
ATOM     34  DB3AGLN A  13      96.100  -2.302  30.169  0.58 29.66           D
ATOM     35  DG2AGLN A  13      95.298  -1.971  28.002  0.58 29.84           D
ATOM     36  DG3AGLN A  13      96.179  -3.194  27.487  0.58 26.84           D
ATOM     37 DE21AGLN A  13      93.665  -2.625  29.702  0.58 32.90           D
ATOM     38 DE22AGLN A  13      92.788  -3.764  29.316  0.58 37.09           D
ATOM     39  N  BGLN A  13      98.681  -2.490  28.191  0.42 24.76           N
ATOM     40  CA BGLN A  13      97.866  -1.677  29.093  0.42 24.29           C
ATOM     41  C  BGLN A  13      98.617  -1.297  30.372  0.42 22.99           C
ATOM     42  O  BGLN A  13      98.433  -0.198  30.878  0.42 27.01           O
ATOM     43  CB BGLN A  13      96.525  -2.340  29.406  0.42 26.14           C
ATOM     44  CG BGLN A  13      95.652  -1.530  30.381  0.42 33.27           C
ATOM     45  CD BGLN A  13      95.330  -0.100  29.916  0.42 36.68           C
ATOM     46  OE1BGLN A  13      94.969   0.134  28.757  0.42 31.67           O
ATOM     47  NE2BGLN A  13      95.448   0.857  30.839  0.42 25.08           N
ATOM     48  DA BGLN A  13      97.662  -0.832  28.645  0.42 27.20           D
ATOM     49  DB2BGLN A  13      96.034  -2.447  28.577  0.42 28.48           D
ATOM     50  DB3BGLN A  13      96.695  -3.213  29.790  0.42 27.73           D
ATOM     51  DG2BGLN A  13      94.806  -1.995  30.490  0.42 33.81           D
ATOM     52  DG3BGLN A  13      96.108  -1.471  31.234  0.42 31.50           D
ATOM     53 DE21BGLN A  13      95.522   0.644  31.668  0.42 28.67           D
ATOM     54 DE22BGLN A  13      95.442   1.683  30.604  0.42 27.64           D
ATOM     55  N   LYS A  14      99.497  -2.166  30.861  1.00 21.77           N
ATOM     56  CA  LYS A  14     100.220  -1.884  32.084  1.00 25.70           C
ATOM     57  C   LYS A  14     101.170  -0.780  31.775  1.00 25.88           C
ATOM     58  O   LYS A  14     101.459   0.050  32.629  1.00 22.65           O
ATOM     59  CB  LYS A  14     101.084  -3.059  32.530  1.00 26.98           C
ATOM     60  CG  LYS A  14     100.374  -4.009  33.493  1.00 52.63           C
ATOM     61  CD  LYS A  14     101.404  -4.810  34.354  1.00 73.04           C
ATOM     62  CE  LYS A  14     100.727  -5.980  35.154  1.00 81.79           C
ATOM     63  NZ  LYS A  14     100.050  -7.029  34.279  1.00 71.80           N
ATOM     64  D   LYS A  14      99.738  -2.891  30.481  1.00 30.84           D
ATOM     65  DA  LYS A  14      99.615  -1.612  32.805  1.00 24.13           D
ATOM     66  DB2 LYS A  14     101.379  -3.549  31.745  1.00 28.60           D
ATOM     67  DB3 LYS A  14     101.864  -2.717  32.995  1.00 30.46           D
ATOM     68  DD2 LYS A  14     102.079  -5.186  33.767  1.00 58.25           D
ATOM     69  DD3 LYS A  14     101.843  -4.203  34.982  1.00 46.82           D
ATOM     70  DE2 LYS A  14     101.405  -6.427  35.687  1.00 64.32           D
ATOM     71  DE3 LYS A  14     100.039  -5.616  35.740  1.00 60.96           D
ATOM     72  DG2 LYS A  14      99.825  -3.470  34.085  1.00 40.59           D
ATOM     73  DG3 LYS A  14      99.810  -4.632  33.005  1.00 37.06           D
ATOM     74  DZ1 LYS A  14      99.670  -7.651  34.789  1.00 62.87           D
ATOM     75  DZ2 LYS A  14      99.416  -6.638  33.789  1.00 55.35           D
ATOM     76  DZ3 LYS A  14     100.639  -7.426  33.738  1.00 57.61           D
END
"""

# (b) 7T5E A 13-15: amide NH2 names permuted between H and D
pdb_b_7t5e_gln_a14 = """\
CRYST1   68.296   42.269   70.413  90.00  98.47  90.00 P 1 21 1
ATOM      1  N   ASN A  13     -26.283  -2.044 -37.427  1.00 23.34           N
ATOM      2  CA  ASN A  13     -26.280  -2.568 -36.058  1.00 16.40           C
ATOM      3  C   ASN A  13     -27.315  -3.664 -35.928  1.00 23.60           C
ATOM      4  O   ASN A  13     -28.476  -3.496 -36.313  1.00 21.20           O
ATOM      5  CB  ASN A  13     -26.555  -1.479 -35.029  1.00 19.12           C
ATOM      6  CG  ASN A  13     -26.211  -1.941 -33.623  1.00 21.24           C
ATOM      7  OD1 ASN A  13     -26.932  -2.725 -33.009  1.00 21.09           O
ATOM      8  ND2 ASN A  13     -25.100  -1.456 -33.119  1.00 22.95           N
ATOM      9  HA  ASN A  13     -25.318  -3.021 -35.830  1.00 21.26           H
ATOM     10  HB2 ASN A  13     -25.948  -0.598 -35.237  1.00 20.31           H
ATOM     11  HB3 ASN A  13     -27.619  -1.240 -35.067  1.00 24.70           H
ATOM     12  H  AASN A  13     -26.642  -1.090 -37.506  0.25 25.94           H
ATOM     13 HD21AASN A  13     -24.808  -1.725 -32.180  0.12 23.76           H
ATOM     14 HD22AASN A  13     -24.531  -0.808 -33.665  0.07 19.85           H
ATOM     15  D  BASN A  13     -26.649  -1.104 -37.524  0.75 25.85           D
ATOM     16 DD21BASN A  13     -24.813  -1.727 -32.173  0.88 24.20           D
ATOM     17 DD22BASN A  13     -24.536  -0.805 -33.666  0.93 19.39           D
ATOM     18  N   GLN A  14     -26.871  -4.802 -35.424  1.00 17.10           N
ATOM     19  CA  GLN A  14     -27.707  -5.973 -35.275  1.00 16.62           C
ATOM     20  C   GLN A  14     -28.612  -5.945 -34.057  1.00 19.80           C
ATOM     21  O   GLN A  14     -29.414  -6.869 -33.897  1.00 27.40           O
ATOM     22  CB  GLN A  14     -26.809  -7.198 -35.224  1.00 17.87           C
ATOM     23  CG  GLN A  14     -26.012  -7.405 -36.488  1.00 16.93           C
ATOM     24  CD  GLN A  14     -26.915  -7.582 -37.760  1.00 28.34           C
ATOM     25  OE1 GLN A  14     -27.607  -8.581 -37.893  1.00 34.38           O
ATOM     26  NE2 GLN A  14     -26.907  -6.593 -38.658  1.00 26.16           N
ATOM     27  HA  GLN A  14     -28.354  -6.057 -36.151  1.00 26.25           H
ATOM     28  HB2 GLN A  14     -26.101  -7.058 -34.408  1.00 20.19           H
ATOM     29  HB3 GLN A  14     -27.422  -8.084 -35.050  1.00 18.17           H
ATOM     30  HG2 GLN A  14     -25.339  -6.555 -36.628  1.00 23.06           H
ATOM     31  HG3 GLN A  14     -25.389  -8.297 -36.371  1.00 19.91           H
ATOM     32  H  AGLN A  14     -25.914  -4.940 -35.100  0.29 18.71           H
ATOM     33 HE21AGLN A  14     -26.319  -5.776 -38.504  0.51 23.46           H
ATOM     34 HE22AGLN A  14     -27.482  -6.657 -39.501  0.23 26.23           H
ATOM     35  D  BGLN A  14     -25.914  -4.970 -35.134  0.71 18.20           D
ATOM     36 DE21BGLN A  14     -27.488  -6.682 -39.485  0.49 26.50           D
ATOM     37 DE22BGLN A  14     -26.330  -5.756 -38.532  0.77 23.75           D
ATOM     38  N   GLY A  15     -28.482  -4.964 -33.179  1.00 19.89           N
ATOM     39  CA  GLY A  15     -29.389  -4.856 -32.060  1.00 18.11           C
ATOM     40  C   GLY A  15     -28.701  -5.079 -30.725  1.00 21.75           C
ATOM     41  O   GLY A  15     -27.773  -5.889 -30.588  1.00 13.78           O
ATOM     42  HA2 GLY A  15     -29.835  -3.859 -32.055  1.00 19.75           H
ATOM     43  HA3 GLY A  15     -30.198  -5.576 -32.179  1.00 32.09           H
ATOM     44  H  AGLY A  15     -27.773  -4.233 -33.223  0.13 18.82           H
ATOM     45  D  BGLY A  15     -27.773  -4.229 -33.189  0.87 18.50           D
END
"""

# (b) 7KCU A 2-4: Lys NH3 names shifted cyclically
pdb_b_7kcu_lys_a3 = """\
CRYST1   60.830   67.100   69.040  90.00  90.00  90.00 P 21 21 21
ATOM      1  N   PHE A   2     -14.930  19.755 -12.199  1.00 23.11           N
ATOM      2  CA  PHE A   2     -14.616  18.377 -11.841  1.00 23.36           C
ATOM      3  C   PHE A   2     -15.861  17.627 -11.398  1.00 23.01           C
ATOM      4  O   PHE A   2     -16.097  16.525 -11.861  1.00 25.31           O
ATOM      5  CB  PHE A   2     -13.559  18.374 -10.753  1.00 27.51           C
ATOM      6  CG  PHE A   2     -13.002  17.034 -10.383  1.00 21.34           C
ATOM      7  CD1 PHE A   2     -12.045  16.371 -11.180  1.00 23.50           C
ATOM      8  CD2 PHE A   2     -13.285  16.545  -9.144  1.00 22.53           C
ATOM      9  CE1 PHE A   2     -11.465  15.090 -10.754  1.00 22.58           C
ATOM     10  CE2 PHE A   2     -12.723  15.379  -8.697  1.00 20.00           C
ATOM     11  CZ  PHE A   2     -11.791  14.621  -9.513  1.00 24.97           C
ATOM     12  HA  PHE A   2     -14.230  17.930 -12.613  1.00 29.42           H
ATOM     13  HB2 PHE A   2     -12.835  18.936 -11.047  1.00 32.00           H
ATOM     14  HB3 PHE A   2     -13.958  18.752  -9.952  1.00 17.13           H
ATOM     15  HD1 PHE A   2     -11.824  16.711 -12.017  1.00 30.27           H
ATOM     16  HD2 PHE A   2     -13.886  16.991  -8.594  1.00 38.10           H
ATOM     17  HE1 PHE A   2     -10.847  14.634 -11.296  1.00 30.19           H
ATOM     18  HE2 PHE A   2     -12.997  15.063  -7.865  1.00 28.10           H
ATOM     19  HZ  PHE A   2     -11.456  13.819  -9.198  1.00 30.72           H
ATOM     20  H  APHE A   2     -14.546  20.344 -11.705  0.01 28.28           H
ATOM     21  D  BPHE A   2     -14.612  20.362 -11.697  0.99 28.25           D
ATOM     22  N   LYS A   3     -16.717  18.226 -10.577  1.00 27.11           N
ATOM     23  CA  LYS A   3     -17.921  17.508 -10.173  1.00 24.87           C
ATOM     24  C   LYS A   3     -18.879  17.239 -11.319  1.00 25.40           C
ATOM     25  O   LYS A   3     -19.470  16.162 -11.367  1.00 22.87           O
ATOM     26  CB  LYS A   3     -18.681  18.218  -9.071  1.00 36.77           C
ATOM     27  CG  LYS A   3     -17.995  18.036  -7.734  1.00 40.44           C
ATOM     28  CD  LYS A   3     -18.331  16.755  -7.078  1.00 39.56           C
ATOM     29  CE  LYS A   3     -17.496  16.724  -5.798  1.00 44.43           C
ATOM     30  NZ  LYS A   3     -17.643  15.457  -5.098  1.00 42.99           N
ATOM     31  HA  LYS A   3     -17.617  16.646  -9.848  1.00 35.12           H
ATOM     32  HB2 LYS A   3     -18.655  19.150  -9.313  1.00 30.45           H
ATOM     33  HB3 LYS A   3     -19.601  17.879  -8.962  1.00 26.26           H
ATOM     34  HG2 LYS A   3     -17.038  18.065  -7.878  1.00 33.99           H
ATOM     35  HG3 LYS A   3     -18.267  18.746  -7.119  1.00 31.75           H
ATOM     36  HD2 LYS A   3     -19.264  16.740  -6.799  1.00 32.01           H
ATOM     37  HD3 LYS A   3     -18.138  15.989  -7.641  1.00 34.70           H
ATOM     38  HE2 LYS A   3     -16.555  16.867  -5.994  1.00 38.59           H
ATOM     39  HE3 LYS A   3     -17.819  17.429  -5.213  1.00 43.84           H
ATOM     40  H  ALYS A   3     -16.641  19.020 -10.249  0.12 29.23           H
ATOM     41  HZ1ALYS A   3     -18.224  15.529  -4.426  0.47 37.79           H
ATOM     42  HZ2ALYS A   3     -17.929  14.838  -5.671  0.85 36.85           H
ATOM     43  HZ3ALYS A   3     -16.849  15.203  -4.782  0.33 35.53           H
ATOM     44  D  BLYS A   3     -16.589  18.998 -10.224  0.88 29.16           D
ATOM     45  DZ1BLYS A   3     -16.831  15.198  -4.846  0.53 35.40           D
ATOM     46  DZ2BLYS A   3     -18.151  15.559  -4.372  0.15 37.78           D
ATOM     47  DZ3BLYS A   3     -18.002  14.830  -5.624  0.67 36.59           D
ATOM     48  N   GLN A   4     -19.028  18.187 -12.251  1.00 26.44           N
ATOM     49  CA  GLN A   4     -19.849  17.928 -13.416  1.00 20.67           C
ATOM     50  C   GLN A   4     -19.244  16.847 -14.307  1.00 18.56           C
ATOM     51  O   GLN A   4     -19.976  16.059 -14.895  1.00 21.84           O
ATOM     52  CB  GLN A   4     -20.074  19.191 -14.236  1.00 22.80           C
ATOM     53  CG  GLN A   4     -20.864  18.991 -15.490  1.00 24.52           C
ATOM     54  CD  GLN A   4     -22.227  18.348 -15.256  1.00 36.10           C
ATOM     55  OE1 GLN A   4     -22.909  18.580 -14.252  1.00 24.22           O
ATOM     56  NE2 GLN A   4     -22.619  17.493 -16.217  1.00 35.95           N
ATOM     57  HA  GLN A   4     -20.709  17.630 -13.099  1.00 26.66           H
ATOM     58  HB2 GLN A   4     -20.592  19.778 -13.662  1.00 29.26           H
ATOM     59  HB3 GLN A   4     -19.238  19.619 -14.491  1.00 28.51           H
ATOM     60  HG2 GLN A   4     -21.031  19.853 -15.878  1.00 31.71           H
ATOM     61  HG3 GLN A   4     -20.373  18.447 -16.128  1.00 18.72           H
ATOM     62  D  BGLN A   4     -18.568  18.912 -12.278  1.00 21.94           D
ATOM     63 DE21BGLN A   4     -23.404  17.149 -16.180  0.46 29.99           D
ATOM     64 DE22BGLN A   4     -22.075  17.287 -16.855  0.86 31.94           D
ATOM     65 HE21AGLN A   4     -22.122  17.328 -16.908  0.54 31.41           H
ATOM     66 HE22AGLN A   4     -23.378  17.089 -16.139  0.14 30.04           H
END
"""

# (b) 7KCU A 40-42: Lys NH3, one H completed to the name left free
pdb_b_7kcu_lys_a41 = """\
CRYST1   60.830   67.100   69.040  90.00  90.00  90.00 P 21 21 21
ATOM      1  N   GLY A  40      -3.498  -9.297 -28.698  1.00 34.09           N
ATOM      2  CA  GLY A  40      -3.772  -9.907 -29.972  1.00 28.36           C
ATOM      3  C   GLY A  40      -5.112  -9.569 -30.541  1.00 36.23           C
ATOM      4  O   GLY A  40      -5.539 -10.227 -31.481  1.00 31.67           O
ATOM      5  HA2 GLY A  40      -3.095  -9.660 -30.613  1.00 30.39           H
ATOM      6  HA3 GLY A  40      -3.751 -10.860 -29.846  1.00 33.22           H
ATOM      7  D  BGLY A  40      -2.982  -8.604 -28.726  1.00 36.40           D
ATOM      8  N   LYS A  41      -5.722  -8.482 -30.080  1.00 37.86           N
ATOM      9  CA  LYS A  41      -7.087  -8.159 -30.414  1.00 23.78           C
ATOM     10  C   LYS A  41      -7.212  -6.675 -30.783  1.00 23.88           C
ATOM     11  O   LYS A  41      -6.529  -5.841 -30.211  1.00 28.30           O
ATOM     12  CB  LYS A  41      -7.968  -8.575 -29.246  1.00 35.48           C
ATOM     13  CG  LYS A  41      -8.356 -10.033 -29.363  1.00 46.19           C
ATOM     14  CD  LYS A  41      -9.114 -10.503 -28.119  1.00 49.69           C
ATOM     15  CE  LYS A  41      -8.687 -11.967 -27.658  1.00 58.61           C
ATOM     16  NZ  LYS A  41      -8.904 -12.942 -28.813  1.00 66.90           N
ATOM     17  HA  LYS A  41      -7.338  -8.669 -31.175  1.00 28.81           H
ATOM     18  HB2 LYS A  41      -7.451  -8.487 -28.434  1.00 38.40           H
ATOM     19  HB3 LYS A  41      -8.754  -8.011 -29.183  1.00 39.06           H
ATOM     20  HG2 LYS A  41      -8.907 -10.190 -30.157  1.00 49.15           H
ATOM     21  HG3 LYS A  41      -7.539 -10.556 -29.421  1.00 39.07           H
ATOM     22  HD2 LYS A  41      -8.939  -9.868 -27.408  1.00 42.80           H
ATOM     23  HD3 LYS A  41     -10.064 -10.522 -28.309  1.00 43.46           H
ATOM     24  HE2 LYS A  41      -7.745 -11.973 -27.411  1.00 50.91           H
ATOM     25  HE3 LYS A  41      -9.231 -12.261 -26.917  1.00 50.52           H
ATOM     26  H  ALYS A  41      -5.348  -7.908 -29.560  0.02 35.86           H
ATOM     27  HZ1ALYS A  41      -9.472 -12.612 -29.414  0.23 50.70           H
ATOM     28  HZ2ALYS A  41      -8.128 -13.112 -29.215  0.37 52.04           H
ATOM     29  HZ3ALYS A  41      -9.224 -13.713 -28.493  1.00 59.31           H
ATOM     30  D  BLYS A  41      -5.331  -7.874 -29.604  0.98 36.01           D
ATOM     31  DZ1BLYS A  41      -9.229 -13.711 -28.497  0.77 59.20           D
ATOM     32  DZ2BLYS A  41      -9.495 -12.605 -29.395  0.63 50.36           D
ATOM     33  N   SER A  42      -8.017  -6.351 -31.788  1.00 30.08           N
ATOM     34  CA  SER A  42      -8.217  -4.966 -32.188  1.00 22.47           C
ATOM     35  C   SER A  42      -9.212  -4.304 -31.234  1.00 22.75           C
ATOM     36  O   SER A  42      -9.871  -4.971 -30.432  1.00 26.30           O
ATOM     37  CB  SER A  42      -8.757  -4.892 -33.615  1.00 23.14           C
ATOM     38  OG  SER A  42      -9.964  -5.647 -33.762  1.00 27.58           O
ATOM     39  HA  SER A  42      -7.371  -4.487 -32.176  1.00 32.60           H
ATOM     40  HB2 SER A  42      -8.898  -3.967 -33.888  1.00 25.63           H
ATOM     41  HB3 SER A  42      -8.053  -5.239 -34.177  1.00 22.37           H
ATOM     42  H  ASER A  42      -8.485  -6.913 -32.238  0.18 32.70           H
ATOM     43  HG ASER A  42     -10.607  -5.242 -33.405  0.26 29.40           H
ATOM     44  D  BSER A  42      -8.427  -6.927 -32.279  0.82 33.02           D
ATOM     45  DG BSER A  42     -10.605  -5.243 -33.415  0.74 29.27           D
END
"""

# (d) 6KK8 A 40-42: Met A 41 backbone H in conformer A, D in B on a blank N
pdb_d_6kk8_met_a41 = """\
CRYST1  133.396  133.396  133.396  90.00  90.00  90.00 P 21 3
ATOM      1  N   LEU A  40      -2.758  24.217  21.522  1.00  4.58           N
ATOM      2  CA  LEU A  40      -3.461  24.736  20.364  1.00  4.17           C
ATOM      3  C   LEU A  40      -3.829  23.634  19.388  1.00  5.01           C
ATOM      4  O   LEU A  40      -4.997  23.486  19.001  1.00  5.39           O
ATOM      5  CB  LEU A  40      -2.599  25.783  19.661  1.00  4.81           C
ATOM      6  CG  LEU A  40      -3.074  26.193  18.260  1.00  5.85           C
ATOM      7  CD1 LEU A  40      -4.458  26.825  18.293  1.00  5.87           C
ATOM      8  CD2 LEU A  40      -2.047  27.099  17.601  1.00  6.93           C
ATOM      9  HA  LEU A  40      -4.397  25.194  20.704  1.00  4.72           H
ATOM     10  HB2 LEU A  40      -2.541  26.671  20.295  1.00  5.57           H
ATOM     11  HB3 LEU A  40      -1.615  25.356  19.496  1.00  4.41           H
ATOM     12  HG  LEU A  40      -3.131  25.281  17.670  1.00  5.34           H
ATOM     13 HD11 LEU A  40      -4.699  27.226  17.311  1.00 10.39           H
ATOM     14 HD12 LEU A  40      -5.189  26.069  18.591  1.00  7.07           H
ATOM     15 HD13 LEU A  40      -4.438  27.664  19.000  1.00  5.73           H
ATOM     16 HD21 LEU A  40      -2.455  27.530  16.682  1.00  8.05           H
ATOM     17 HD22 LEU A  40      -1.762  27.911  18.275  1.00  4.92           H
ATOM     18 HD23 LEU A  40      -1.167  26.513  17.347  1.00  8.05           H
ATOM     19  H  ALEU A  40      -1.859  24.668  21.704  0.41  6.20           H
ATOM     20  D  BLEU A  40      -1.874  24.677  21.712  0.59  6.31           D
ATOM     21  N   MET A  41      -2.837  22.885  18.907  1.00  4.34           N
ATOM     22  C   MET A  41      -3.903  20.770  18.226  1.00  5.42           C
ATOM     23  O   MET A  41      -4.681  20.231  17.425  1.00  5.22           O
ATOM     24  CA AMET A  41      -3.089  21.977  17.798  0.35  4.37           C
ATOM     25  CB AMET A  41      -1.771  21.552  17.160  0.35  5.38           C
ATOM     26  CG AMET A  41      -1.047  22.686  16.464  0.35  6.33           C
ATOM     27  SD AMET A  41      -1.608  22.916  14.759  0.35  8.97           S
ATOM     28  CE AMET A  41      -3.188  23.730  14.941  0.35 11.42           C
ATOM     29  H  AMET A  41      -1.868  22.896  19.230  0.90  6.45           H
ATOM     30  HA AMET A  41      -3.680  22.500  17.050  0.35  6.91           H
ATOM     31  HB2AMET A  41      -1.119  21.165  17.944  0.35  7.32           H
ATOM     32  HB3AMET A  41      -1.977  20.774  16.422  0.35  6.78           H
ATOM     33  HG2AMET A  41      -1.237  23.611  17.007  0.35  7.26           H
ATOM     34  HG3AMET A  41       0.024  22.493  16.450  0.35  8.45           H
ATOM     35  HE1AMET A  41      -3.953  22.961  15.028  0.35  9.56           H
ATOM     36  HE2AMET A  41      -3.180  24.342  15.839  0.35  7.86           H
ATOM     37  HE3AMET A  41      -3.388  24.354  14.070  0.35  9.87           H
ATOM     38  CA BMET A  41      -3.108  21.993  17.793  0.65  4.22           C
ATOM     39  CB BMET A  41      -1.814  21.632  17.061  0.65  5.22           C
ATOM     40  CG BMET A  41      -1.231  22.799  16.271  0.65  6.83           C
ATOM     41  SD BMET A  41      -2.361  23.279  14.917  0.65 11.52           S
ATOM     42  CE BMET A  41      -1.741  24.907  14.505  0.65 16.10           C
ATOM     43  D  BMET A  41      -1.879  22.877  19.256  0.10  6.33           D
ATOM     44  HA BMET A  41      -3.750  22.518  17.083  0.65  6.95           H
ATOM     45  HB2BMET A  41      -1.079  21.322  17.805  0.65  7.44           H
ATOM     46  HB3BMET A  41      -2.030  20.818  16.364  0.65  6.77           H
ATOM     47  HG2BMET A  41      -1.089  23.650  16.934  0.65  7.26           H
ATOM     48  HG3BMET A  41      -0.280  22.530  15.812  0.65  8.90           H
ATOM     49  HE1BMET A  41      -1.765  25.560  15.382  0.65 13.24           H
ATOM     50  HE2BMET A  41      -0.722  24.815  14.129  0.65 10.36           H
ATOM     51  HE3BMET A  41      -2.381  25.333  13.734  0.65 10.95           H
ATOM     52  N   ASN A  42      -3.729  20.314  19.468  1.00  4.19           N
ATOM     53  CA  ASN A  42      -4.597  19.267  19.981  1.00  4.30           C
ATOM     54  C   ASN A  42      -6.049  19.696  19.916  1.00  5.90           C
ATOM     55  O   ASN A  42      -6.893  18.977  19.361  1.00  6.20           O
ATOM     56  CB  ASN A  42      -4.203  18.876  21.406  1.00  4.79           C
ATOM     57  CG  ASN A  42      -3.176  17.788  21.419  1.00  6.64           C
ATOM     58  OD1 ASN A  42      -3.514  16.638  21.273  1.00  8.22           O
ATOM     59  ND2 ASN A  42      -1.927  18.144  21.572  1.00  6.74           N
ATOM     60  HA  ASN A  42      -4.519  18.390  19.336  1.00  7.77           H
ATOM     61  HB2 ASN A  42      -3.798  19.734  21.935  1.00  5.69           H
ATOM     62  HB3 ASN A  42      -5.071  18.478  21.926  1.00  4.49           H
ATOM     63  H  AASN A  42      -3.055  20.684  20.139  0.92  6.71           H
ATOM     64 HD21AASN A  42      -1.689  19.129  21.688  0.02  8.69           H
ATOM     65 HD22AASN A  42      -1.191  17.436  21.586  0.20  9.81           H
ATOM     66  D  BASN A  42      -3.013  20.637  20.119  0.08  6.54           D
ATOM     67 DD21BASN A  42      -1.672  19.119  21.659  0.98  8.89           D
ATOM     68 DD22BASN A  42      -1.214  17.423  21.560  0.80 10.00           D
END
"""

# (e) 7T5E A 20-22: Arg A 21 HE/DE 0.23 A apart on NE (one H)
pdb_e_7t5e_arg_a21 = """\
CRYST1   68.296   42.269   70.413  90.00  98.47  90.00 P 1 21 1
ATOM      1  N   VAL A  20     -24.884 -11.104 -29.077  1.00 15.04           N
ATOM      2  CA  VAL A  20     -23.817 -10.784 -28.126  1.00 12.33           C
ATOM      3  C   VAL A  20     -24.336 -11.136 -26.748  1.00 11.49           C
ATOM      4  O   VAL A  20     -25.421 -10.683 -26.376  1.00 13.44           O
ATOM      5  CB  VAL A  20     -23.418  -9.301 -28.187  1.00 13.55           C
ATOM      6  CG1 VAL A  20     -22.219  -9.069 -27.245  1.00 16.66           C
ATOM      7  CG2 VAL A  20     -23.141  -8.837 -29.624  1.00 13.09           C
ATOM      8  HA  VAL A  20     -22.940 -11.386 -28.345  1.00 10.37           H
ATOM      9  HB  VAL A  20     -24.247  -8.714 -27.786  1.00 10.91           H
ATOM     10 HG11 VAL A  20     -21.540  -8.370 -27.719  1.00  9.92           H
ATOM     11 HG12 VAL A  20     -22.596  -8.633 -26.328  1.00 13.84           H
ATOM     12 HG13 VAL A  20     -21.754 -10.017 -27.009  1.00 16.77           H
ATOM     13 HG21 VAL A  20     -22.473  -7.965 -29.616  1.00 18.08           H
ATOM     14 HG22 VAL A  20     -22.682  -9.677 -30.157  1.00 22.69           H
ATOM     15 HG23 VAL A  20     -24.087  -8.566 -30.110  1.00 17.37           H
ATOM     16  H  AVAL A  20     -25.740 -10.592 -28.901  1.00 11.54           H
ATOM     17  N   ARG A  21     -23.602 -11.977 -26.022  1.00 11.30           N
ATOM     18  CA  ARG A  21     -24.003 -12.423 -24.678  1.00 13.46           C
ATOM     19  C   ARG A  21     -23.605 -11.325 -23.690  1.00 15.41           C
ATOM     20  O   ARG A  21     -22.439 -11.167 -23.379  1.00 14.25           O
ATOM     21  CB  ARG A  21     -23.328 -13.737 -24.318  1.00 15.63           C
ATOM     22  CG  ARG A  21     -23.587 -14.915 -25.236  1.00 17.35           C
ATOM     23  CD  ARG A  21     -22.718 -16.037 -24.869  1.00 17.39           C
ATOM     24  NE  ARG A  21     -22.872 -17.143 -25.819  1.00 12.95           N
ATOM     25  CZ  ARG A  21     -22.267 -18.323 -25.714  1.00 22.59           C
ATOM     26  NH1 ARG A  21     -21.455 -18.584 -24.694  1.00 21.39           N
ATOM     27  NH2 ARG A  21     -22.486 -19.251 -26.628  1.00 16.88           N
ATOM     28  HA  ARG A  21     -25.090 -12.566 -24.649  1.00 20.38           H
ATOM     29  HB2 ARG A  21     -22.252 -13.553 -24.294  1.00  6.75           H
ATOM     30  HB3 ARG A  21     -23.698 -14.039 -23.335  1.00 10.60           H
ATOM     31  HG2 ARG A  21     -24.606 -15.268 -25.083  1.00 25.89           H
ATOM     32  HG3 ARG A  21     -23.454 -14.651 -26.295  1.00 15.71           H
ATOM     33  HD2 ARG A  21     -21.665 -15.739 -24.838  1.00 11.46           H
ATOM     34  HD3 ARG A  21     -23.006 -16.394 -23.891  1.00 16.80           H
ATOM     35  H  AARG A  21     -22.683 -12.312 -26.321  1.00 10.59           H
ATOM     36  HE AARG A  21     -23.508 -17.007 -26.606  0.01 17.49           H
ATOM     37 HH11AARG A  21     -20.999 -19.495 -24.631  0.01 18.85           H
ATOM     38 HH21AARG A  21     -22.023 -20.156 -26.553  0.13 17.05           H
ATOM     39  DE BARG A  21     -23.658 -17.080 -26.450  0.99 17.39           D
ATOM     40 DH11BARG A  21     -20.991 -19.486 -24.606  1.00 18.52           D
ATOM     41 DH12BARG A  21     -21.302 -17.861 -23.994  0.99 20.50           D
ATOM     42 DH21BARG A  21     -22.040 -20.175 -26.575  1.00 16.60           D
ATOM     43 DH22BARG A  21     -23.120 -19.018 -27.396  0.87 14.35           D
ATOM     44  N   VAL A  22     -24.552 -10.524 -23.232  1.00 12.34           N
ATOM     45  CA  VAL A  22     -24.183  -9.319 -22.491  1.00 16.31           C
ATOM     46  C   VAL A  22     -24.424  -9.498 -20.998  1.00 16.43           C
ATOM     47  O   VAL A  22     -25.463 -10.061 -20.602  1.00 15.20           O
ATOM     48  CB  VAL A  22     -24.939  -8.071 -22.978  1.00 19.16           C
ATOM     49  CG1 VAL A  22     -24.441  -7.617 -24.338  1.00 16.61           C
ATOM     50  CG2 VAL A  22     -26.458  -8.290 -23.012  1.00 20.27           C
ATOM     51  HA  VAL A  22     -23.117  -9.126 -22.640  1.00 15.33           H
ATOM     52  HB  VAL A  22     -24.742  -7.298 -22.239  1.00 14.80           H
ATOM     53 HG11 VAL A  22     -23.361  -7.459 -24.322  1.00 10.25           H
ATOM     54 HG12 VAL A  22     -24.716  -8.394 -25.055  1.00 14.60           H
ATOM     55 HG13 VAL A  22     -24.932  -6.677 -24.600  1.00 11.07           H
ATOM     56 HG21 VAL A  22     -26.724  -9.181 -23.574  1.00 21.36           H
ATOM     57 HG22 VAL A  22     -26.805  -8.336 -21.977  1.00 16.17           H
ATOM     58 HG23 VAL A  22     -26.894  -7.440 -23.533  1.00 19.63           H
ATOM     59  H  AVAL A  22     -25.550 -10.698 -23.348  0.95 16.45           H
ATOM     60  D  BVAL A  22     -25.554 -10.675 -23.346  0.05 16.54           D
END
"""

# control (b) 7KCU B 124-126: the free name is carried by a D
pdb_ctl_7kcu_lys_b125 = """\
CRYST1   60.830   67.100   69.040  90.00  90.00  90.00 P 21 21 21
ATOM      1  N   GLY B 124     -39.286   6.594  -4.689  1.00 28.46           N
ATOM      2  CA  GLY B 124     -40.158   6.473  -5.857  1.00 23.44           C
ATOM      3  C   GLY B 124     -41.388   7.340  -5.762  1.00 19.38           C
ATOM      4  O   GLY B 124     -41.780   8.001  -6.732  1.00 21.66           O
ATOM      5  HA2 GLY B 124     -39.704   6.753  -6.668  1.00 24.81           H
ATOM      6  HA3 GLY B 124     -40.443   5.546  -5.941  1.00 23.28           H
ATOM      7  H  AGLY B 124     -39.190   5.835  -4.275  1.00 21.94           H
ATOM      8  N   LYS B 125     -41.962   7.423  -4.597  1.00 22.91           N
ATOM      9  CA  LYS B 125     -43.083   8.288  -4.442  1.00 26.59           C
ATOM     10  C   LYS B 125     -42.699   9.758  -4.535  1.00 32.68           C
ATOM     11  O   LYS B 125     -43.379  10.502  -5.152  1.00 23.60           O
ATOM     12  CB  LYS B 125     -43.778   8.028  -3.148  1.00 24.88           C
ATOM     13  CG  LYS B 125     -44.571   6.771  -3.178  1.00 52.11           C
ATOM     14  CD  LYS B 125     -44.627   6.175  -1.797  1.00 54.41           C
ATOM     15  CE  LYS B 125     -45.926   5.473  -1.541  1.00 58.21           C
ATOM     16  NZ  LYS B 125     -45.909   4.096  -2.056  1.00 59.08           N
ATOM     17  HA  LYS B 125     -43.686   8.094  -5.176  1.00 20.96           H
ATOM     18  HB2 LYS B 125     -43.122   7.978  -2.435  1.00 34.32           H
ATOM     19  HB3 LYS B 125     -44.364   8.773  -2.942  1.00 36.18           H
ATOM     20  HG2 LYS B 125     -45.468   6.952  -3.499  1.00 38.11           H
ATOM     21  HG3 LYS B 125     -44.171   6.140  -3.796  1.00 39.98           H
ATOM     22  HD2 LYS B 125     -43.895   5.549  -1.686  1.00 48.10           H
ATOM     23  HD3 LYS B 125     -44.504   6.876  -1.138  1.00 43.89           H
ATOM     24  HE2 LYS B 125     -46.105   5.460  -0.588  1.00 51.25           H
ATOM     25  HE3 LYS B 125     -46.649   5.966  -1.959  1.00 38.39           H
ATOM     26  H  ALYS B 125     -41.722   6.993  -3.892  0.02 25.69           H
ATOM     27  HZ1ALYS B 125     -45.101   3.910  -2.378  0.46 44.34           H
ATOM     28  HZ2ALYS B 125     -46.103   3.530  -1.397  1.00 51.46           H
ATOM     29  D  BLYS B 125     -41.722   6.993  -3.892  0.98 22.88           D
ATOM     30  DZ1BLYS B 125     -46.088   3.529  -1.394  1.00 52.68           D
ATOM     31  DZ2BLYS B 125     -46.523   4.010  -2.694  0.54 48.35           D
ATOM     32  N   ASN B 126     -41.581  10.134  -3.951  1.00 23.24           N
ATOM     33  CA  ASN B 126     -41.156  11.506  -3.990  1.00 17.83           C
ATOM     34  C   ASN B 126     -40.778  11.905  -5.375  1.00 22.58           C
ATOM     35  O   ASN B 126     -40.993  13.056  -5.755  1.00 22.28           O
ATOM     36  CB  ASN B 126     -40.051  11.711  -2.995  1.00 19.48           C
ATOM     37  CG  ASN B 126     -40.567  11.525  -1.615  1.00 26.39           C
ATOM     38  OD1 ASN B 126     -41.782  11.658  -1.400  1.00 24.70           O
ATOM     39  ND2 ASN B 126     -39.684  11.230  -0.671  1.00 25.69           N
ATOM     40  HA  ASN B 126     -41.868  12.093  -3.739  1.00 14.42           H
ATOM     41  HB2 ASN B 126     -39.330  11.077  -3.141  1.00 16.40           H
ATOM     42  HB3 ASN B 126     -39.693  12.606  -3.105  1.00 30.38           H
ATOM     43  H  AASN B 126     -41.060   9.583  -3.533  0.88 23.96           H
ATOM     44 HD21AASN B 126     -39.946  11.114   0.139  0.01 24.88           H
ATOM     45 HD22AASN B 126     -38.852  11.150  -0.871  0.09 24.35           H
ATOM     46  D  BASN B 126     -41.065   9.630  -3.485  0.12 24.01           D
ATOM     47 DD21BASN B 126     -38.831  11.217  -0.835  0.99 24.17           D
ATOM     48 DD22BASN B 126     -39.981  11.007   0.110  0.91 24.37           D
END
"""

# control (b) 7KCU B 3-5: no D within 0.20 A of the H
pdb_ctl_7kcu_gln_b4 = """\
CRYST1   60.830   67.100   69.040  90.00  90.00  90.00 P 21 21 21
ATOM      1  N   LYS B   3     -27.953   5.326   0.770  1.00 23.50           N
ATOM      2  CA  LYS B   3     -26.735   5.753   0.129  1.00 19.96           C
ATOM      3  C   LYS B   3     -25.810   4.594  -0.197  1.00 22.78           C
ATOM      4  O   LYS B   3     -25.226   4.543  -1.293  1.00 21.89           O
ATOM      5  CB  LYS B   3     -26.079   6.800   1.016  1.00 26.82           C
ATOM      6  CG  LYS B   3     -25.266   7.779   0.286  1.00 42.46           C
ATOM      7  CD  LYS B   3     -26.188   8.793  -0.478  1.00 51.10           C
ATOM      8  CE  LYS B   3     -25.428  10.057  -0.558  1.00 53.71           C
ATOM      9  NZ  LYS B   3     -25.803  10.968  -1.681  1.00 41.57           N
ATOM     10  HA  LYS B   3     -26.989   6.149  -0.744  1.00 21.62           H
ATOM     11  HB2 LYS B   3     -26.737   7.361   1.444  1.00 26.19           H
ATOM     12  HB3 LYS B   3     -25.493   6.371   1.657  1.00 23.56           H
ATOM     13  HG2 LYS B   3     -24.784   8.268   0.963  1.00 31.29           H
ATOM     14  HG3 LYS B   3     -24.647   7.324  -0.309  1.00 28.97           H
ATOM     15  HD2 LYS B   3     -26.356   8.507  -1.395  1.00 34.75           H
ATOM     16  HD3 LYS B   3     -27.006   8.961   0.018  1.00 38.68           H
ATOM     17  HE2 LYS B   3     -25.546  10.523   0.288  1.00 34.63           H
ATOM     18  HE3 LYS B   3     -24.498   9.812  -0.698  1.00 39.99           H
ATOM     19  H  ALYS B   3     -28.051   5.530   1.599  0.61 25.45           H
ATOM     20  HZ1ALYS B   3     -25.625  10.594  -2.469  0.52 34.08           H
ATOM     21  HZ2ALYS B   3     -26.671  11.157  -1.633  0.29 33.77           H
ATOM     22  HZ3ALYS B   3     -25.344  11.728  -1.617  0.24 36.56           H
ATOM     23  D  BLYS B   3     -28.062   5.577   1.579  0.39 25.49           D
ATOM     24  DZ1BLYS B   3     -25.632  10.553  -2.453  0.48 33.73           D
ATOM     25  DZ2BLYS B   3     -26.678  11.136  -1.623  0.71 33.90           D
ATOM     26  DZ3BLYS B   3     -25.346  11.736  -1.634  0.76 36.50           D
ATOM     27  N   GLN B   4     -25.781   3.581   0.624  1.00 27.57           N
ATOM     28  CA  GLN B   4     -24.946   2.438   0.313  1.00 20.91           C
ATOM     29  C   GLN B   4     -25.573   1.559  -0.801  1.00 19.19           C
ATOM     30  O   GLN B   4     -24.870   0.960  -1.628  1.00 19.88           O
ATOM     31  CB  GLN B   4     -24.714   1.612   1.593  1.00 30.13           C
ATOM     32  CG  GLN B   4     -24.021   0.259   1.416  1.00 26.38           C
ATOM     33  CD  GLN B   4     -22.680   0.430   0.756  1.00 31.34           C
ATOM     34  OE1 GLN B   4     -21.943   1.407   0.999  1.00 28.63           O
ATOM     35  NE2 GLN B   4     -22.368  -0.500  -0.110  1.00 34.77           N
ATOM     36  HA  GLN B   4     -24.077   2.757   0.005  1.00 28.57           H
ATOM     37  HB2 GLN B   4     -24.099   2.098   2.159  1.00 17.16           H
ATOM     38  HB3 GLN B   4     -25.551   1.475   2.069  1.00 27.79           H
ATOM     39  HG2 GLN B   4     -23.857  -0.107   2.304  1.00 31.87           H
ATOM     40  HG3 GLN B   4     -24.560  -0.344   0.870  1.00 34.98           H
ATOM     41  H  AGLN B   4     -26.217   3.525   1.361  0.29 26.82           H
ATOM     42 HE21AGLN B   4     -22.922  -1.132  -0.249  0.44 30.50           H
ATOM     43 HE22AGLN B   4     -21.631  -0.456  -0.533  0.19 29.32           H
ATOM     44  D  BGLN B   4     -26.275   3.486   1.315  0.71 26.91           D
ATOM     45 DE21BGLN B   4     -21.803  -0.336  -0.745  0.56 24.59           D
ATOM     46 DE22BGLN B   4     -22.719  -1.291  -0.048  0.81 33.60           D
ATOM     47  N   ASP B   5     -26.881   1.502  -0.840  1.00 20.17           N
ATOM     48  CA  ASP B   5     -27.557   0.770  -1.900  1.00 20.07           C
ATOM     49  C   ASP B   5     -27.159   1.303  -3.281  1.00 29.00           C
ATOM     50  O   ASP B   5     -27.020   0.543  -4.228  1.00 21.26           O
ATOM     51  CB  ASP B   5     -29.079   0.885  -1.729  1.00 23.95           C
ATOM     52  CG  ASP B   5     -29.631   0.049  -0.602  1.00 24.38           C
ATOM     53  OD1 ASP B   5     -30.859   0.274  -0.345  1.00 24.53           O
ATOM     54  OD2 ASP B   5     -28.864  -0.720   0.033  1.00 26.81           O
ATOM     55  HA  ASP B   5     -27.328  -0.161  -1.822  1.00 17.96           H
ATOM     56  HB2 ASP B   5     -29.317   1.813  -1.515  1.00 22.22           H
ATOM     57  HB3 ASP B   5     -29.535   0.575  -2.525  1.00 29.38           H
ATOM     58  H  AASP B   5     -27.408   1.893  -0.280  0.46 23.67           H
ATOM     59  D  BASP B   5     -27.391   1.881  -0.264  0.54 23.73           D
END
"""

# control (e) 7KCU A 17-19: Gln A 18 NE2 (two H), D out of the amide plane
pdb_ctl_7kcu_gln_a18 = """\
CRYST1   60.830   67.100   69.040  90.00  90.00  90.00 P 21 21 21
ATOM      1  N   ALA A  17     -15.752   4.347 -20.308  1.00 21.20           N
ATOM      2  CA  ALA A  17     -14.751   3.486 -19.682  1.00 20.61           C
ATOM      3  C   ALA A  17     -13.986   2.686 -20.748  1.00 28.92           C
ATOM      4  O   ALA A  17     -12.757   2.649 -20.735  1.00 22.33           O
ATOM      5  CB  ALA A  17     -15.330   2.512 -18.686  1.00 24.95           C
ATOM      6  HA  ALA A  17     -14.166   4.086 -19.186  1.00 26.50           H
ATOM      7  HB1 ALA A  17     -14.638   1.856 -18.526  1.00 22.23           H
ATOM      8  HB2 ALA A  17     -15.608   2.917 -17.851  1.00 26.46           H
ATOM      9  HB3 ALA A  17     -16.112   2.103 -19.085  1.00 26.07           H
ATOM     10  H  AALA A  17     -16.572   4.209 -20.097  0.80 24.92           H
ATOM     11  D  BALA A  17     -16.562   4.203 -20.062  0.20 25.12           D
ATOM     12  N   GLN A  18     -14.676   2.090 -21.696  1.00 28.68           N
ATOM     13  CA  GLN A  18     -13.962   1.377 -22.720  1.00 28.68           C
ATOM     14  C   GLN A  18     -13.172   2.256 -23.722  1.00 26.53           C
ATOM     15  O   GLN A  18     -12.144   1.871 -24.111  1.00 22.06           O
ATOM     16  CB  GLN A  18     -14.697   0.144 -23.290  1.00 30.79           C
ATOM     17  CG  GLN A  18     -16.086   0.325 -23.807  1.00 46.60           C
ATOM     18  CD  GLN A  18     -16.845  -0.964 -24.110  1.00 38.21           C
ATOM     19  OE1 GLN A  18     -16.300  -2.003 -24.298  1.00 44.16           O
ATOM     20  NE2 GLN A  18     -18.108  -0.854 -24.181  1.00 46.89           N
ATOM     21  HA  GLN A  18     -13.230   0.958 -22.240  1.00 28.85           H
ATOM     22  HB2 GLN A  18     -14.157  -0.215 -24.011  1.00 41.35           H
ATOM     23  HB3 GLN A  18     -14.726  -0.531 -22.594  1.00 27.96           H
ATOM     24  HG2 GLN A  18     -16.593   0.836 -23.156  1.00 40.60           H
ATOM     25  HG3 GLN A  18     -16.047   0.857 -24.617  1.00 37.09           H
ATOM     26  H  AGLN A  18     -15.533   2.086 -21.763  0.16 28.26           H
ATOM     27 HE21AGLN A  18     -18.482  -0.092 -24.041  0.24 36.06           H
ATOM     28 HE22AGLN A  18     -18.590  -1.541 -24.369  0.55 40.89           H
ATOM     29  D  BGLN A  18     -15.533   2.086 -21.763  0.84 28.45           D
ATOM     30 DE21BGLN A  18     -18.472  -0.521 -24.886  0.76 36.12           D
ATOM     31 DE22BGLN A  18     -18.600  -1.112 -23.525  0.45 39.76           D
ATOM     32  N   ASP A  19     -13.645   3.443 -24.055  1.00 20.56           N
ATOM     33  CA  ASP A  19     -12.924   4.338 -24.922  1.00 18.96           C
ATOM     34  C   ASP A  19     -11.663   4.812 -24.279  1.00 27.17           C
ATOM     35  O   ASP A  19     -10.642   4.944 -24.946  1.00 23.56           O
ATOM     36  CB  ASP A  19     -13.785   5.549 -25.228  1.00 31.19           C
ATOM     37  CG  ASP A  19     -14.662   5.350 -26.485  1.00 35.98           C
ATOM     38  OD1 ASP A  19     -14.684   4.220 -27.003  1.00 48.28           O
ATOM     39  OD2 ASP A  19     -15.378   6.300 -26.859  1.00 54.99           O
ATOM     40  HA  ASP A  19     -12.692   3.920 -25.762  1.00 25.69           H
ATOM     41  HB2 ASP A  19     -14.410   5.633 -24.497  1.00 26.93           H
ATOM     42  HB3 ASP A  19     -13.235   6.342 -25.367  1.00 23.51           H
ATOM     43  D  BASP A  19     -14.364   3.791 -23.721  1.00 29.74           D
END
"""

# (c) 7T5E B 24-26: Tyr B 25 HH on OH A (0.65, conformer A 0.65), DH on
# OH B (0.35, conformer B 0.35): one isotope per conformer
pdb_c_7t5e_tyr_b25 = """\
CRYST1   68.296   42.269   70.413  90.00  98.47  90.00 P 1 21 1
ATOM      1  N   THR B  24     -26.961  -5.468  -6.207  1.00 20.38           N
ATOM      2  CA  THR B  24     -27.028  -6.786  -6.835  1.00 20.32           C
ATOM      3  C   THR B  24     -25.625  -7.338  -7.070  1.00 14.52           C
ATOM      4  O   THR B  24     -25.370  -8.534  -6.878  1.00 22.10           O
ATOM      5  CB  THR B  24     -27.780  -6.670  -8.168  1.00 24.71           C
ATOM      6  OG1 THR B  24     -28.974  -5.858  -8.029  1.00 27.68           O
ATOM      7  CG2 THR B  24     -28.110  -8.058  -8.706  1.00 31.92           C
ATOM      8  HA  THR B  24     -27.555  -7.493  -6.195  1.00 21.82           H
ATOM      9  HB  THR B  24     -27.102  -6.145  -8.830  1.00 22.88           H
ATOM     10 HG21 THR B  24     -28.214  -8.745  -7.863  1.00 29.11           H
ATOM     11 HG22 THR B  24     -29.051  -8.025  -9.251  1.00 31.80           H
ATOM     12 HG23 THR B  24     -27.328  -8.425  -9.373  1.00 29.32           H
ATOM     13  D  BTHR B  24     -27.120  -4.681  -6.844  1.00 14.69           D
ATOM     14  DG1BTHR B  24     -29.364  -5.905  -7.128  0.78 23.47           D
ATOM     15  HG1ATHR B  24     -29.324  -5.928  -7.117  0.22 24.01           H
ATOM     16  N   TYR B  25     -24.693  -6.450  -7.428  1.00 11.52           N
ATOM     17  C   TYR B  25     -22.324  -6.412  -6.782  1.00 12.31           C
ATOM     18  O   TYR B  25     -22.506  -5.436  -6.072  1.00 12.70           O
ATOM     19  D   TYR B  25     -24.862  -5.441  -7.456  1.00 17.07           D
ATOM     20  CA ATYR B  25     -23.346  -6.866  -7.806  0.65 10.53           C
ATOM     21  CB ATYR B  25     -22.968  -6.351  -9.196  0.65 13.23           C
ATOM     22  CG ATYR B  25     -23.969  -6.828 -10.212  0.65 12.62           C
ATOM     23  CD1ATYR B  25     -23.985  -8.150 -10.618  0.65 22.24           C
ATOM     24  CD2ATYR B  25     -24.948  -5.981 -10.702  0.65 22.27           C
ATOM     25  CE1ATYR B  25     -24.906  -8.602 -11.516  0.65 24.58           C
ATOM     26  CE2ATYR B  25     -25.878  -6.426 -11.611  0.65 22.27           C
ATOM     27  CZ ATYR B  25     -25.855  -7.750 -11.997  0.65 22.01           C
ATOM     28  OH ATYR B  25     -26.761  -8.203 -12.922  0.65 28.19           O
ATOM     29  HA ATYR B  25     -23.308  -7.949  -7.827  0.65 12.38           H
ATOM     30  HB2ATYR B  25     -22.992  -5.260  -9.177  0.65 15.05           H
ATOM     31  HB3ATYR B  25     -21.972  -6.712  -9.484  0.65 16.28           H
ATOM     32  HD1ATYR B  25     -23.249  -8.836 -10.238  0.65 21.90           H
ATOM     33  HD2ATYR B  25     -24.962  -4.952 -10.392  0.65 19.16           H
ATOM     34  HE1ATYR B  25     -24.881  -9.635 -11.826  0.65 19.65           H
ATOM     35  HE2ATYR B  25     -26.628  -5.761 -12.011  0.65 27.24           H
ATOM     36  HH ATYR B  25     -26.258  -8.472 -13.729  0.65 24.15           H
ATOM     37  CA BTYR B  25     -23.340  -6.785  -7.856  0.35 10.74           C
ATOM     38  CB BTYR B  25     -23.058  -6.025  -9.146  0.35 14.88           C
ATOM     39  CG BTYR B  25     -21.756  -6.257  -9.832  0.35 14.10           C
ATOM     40  CD1BTYR B  25     -21.452  -7.490 -10.390  0.35 16.43           C
ATOM     41  CD2BTYR B  25     -20.868  -5.209 -10.034  0.35 15.32           C
ATOM     42  CE1BTYR B  25     -20.279  -7.685 -11.079  0.35 15.77           C
ATOM     43  CE2BTYR B  25     -19.697  -5.400 -10.723  0.35 14.70           C
ATOM     44  CZ BTYR B  25     -19.415  -6.628 -11.253  0.35 13.72           C
ATOM     45  OH BTYR B  25     -18.220  -6.810 -11.895  0.35 10.58           O
ATOM     46  DH BTYR B  25     -17.694  -7.460 -11.353  0.35 13.48           D
ATOM     47  HA BTYR B  25     -23.263  -7.855  -8.035  0.35 12.90           H
ATOM     48  HB2BTYR B  25     -23.850  -6.262  -9.860  0.35 16.04           H
ATOM     49  HB3BTYR B  25     -23.098  -4.964  -8.903  0.35 15.10           H
ATOM     50  HD1BTYR B  25     -22.139  -8.308 -10.268  0.35 17.64           H
ATOM     51  HD2BTYR B  25     -21.093  -4.234  -9.634  0.35 16.40           H
ATOM     52  HE1BTYR B  25     -20.049  -8.660 -11.489  0.35 17.67           H
ATOM     53  HE2BTYR B  25     -19.003  -4.585 -10.855  0.35 15.66           H
ATOM     54  N   ASN B  26     -21.248  -7.180  -6.695  1.00  9.86           N
ATOM     55  CA  ASN B  26     -20.267  -7.014  -5.626  1.00 14.46           C
ATOM     56  C   ASN B  26     -19.040  -6.173  -6.017  1.00 12.39           C
ATOM     57  O   ASN B  26     -18.050  -6.147  -5.269  1.00 15.56           O
ATOM     58  CB  ASN B  26     -19.844  -8.387  -5.105  1.00 11.05           C
ATOM     59  CG  ASN B  26     -19.123  -9.230  -6.137  1.00  9.06           C
ATOM     60  OD1 ASN B  26     -18.894  -8.838  -7.284  1.00 13.07           O
ATOM     61  ND2 ASN B  26     -18.771 -10.399  -5.725  1.00 14.50           N
ATOM     62  HA  ASN B  26     -20.734  -6.513  -4.778  1.00 11.13           H
ATOM     63  HB2 ASN B  26     -19.151  -8.252  -4.275  1.00 16.54           H
ATOM     64  HB3 ASN B  26     -20.705  -8.955  -4.758  1.00 20.91           H
ATOM     65  D  BASN B  26     -21.097  -7.972  -7.315  1.00 10.28           D
ATOM     66 DD21BASN B  26     -18.998 -10.676  -4.771  1.00 15.49           D
ATOM     67 DD22BASN B  26     -18.297 -11.030  -6.347  0.95 20.65           D
ATOM     68 HD21AASN B  26     -18.980 -10.687  -4.769  0.05 15.72           H
END
"""

# control (c) 7T5E B 81-83: Ser B 82 H on N A (0.15, conformer 0.56), D on
# N B (0.85, conformer 0.44): H + D = 1, consistent with an exchange pair
# over the two copies, which the occupancies do not distribute
pdb_ctl_7t5e_ser_b82 = """\
CRYST1   68.296   42.269   70.413  90.00  98.47  90.00 P 1 21 1
ATOM      1  N   ASP B  81     -11.576 -11.702  -1.277  1.00 18.52           N
ATOM      2  CA  ASP B  81     -10.253 -12.182  -1.673  1.00 22.89           C
ATOM      3  C   ASP B  81      -9.153 -11.635  -0.763  1.00 18.26           C
ATOM      4  O   ASP B  81      -9.161 -10.465  -0.377  1.00 18.25           O
ATOM      5  CB  ASP B  81      -9.946 -11.796  -3.119  1.00 22.45           C
ATOM      6  CG  ASP B  81      -8.762 -12.585  -3.690  1.00 30.22           C
ATOM      7  OD1 ASP B  81      -9.009 -13.673  -4.226  1.00 35.50           O
ATOM      8  OD2 ASP B  81      -7.594 -12.159  -3.580  1.00 24.83           O
ATOM      9  HA  ASP B  81     -10.216 -13.265  -1.618  1.00 19.86           H
ATOM     10  HB2 ASP B  81     -10.798 -12.031  -3.759  1.00 30.20           H
ATOM     11  HB3 ASP B  81      -9.718 -10.737  -3.187  1.00 17.62           H
ATOM     12  H  AASP B  81     -12.081 -11.198  -2.005  0.31 17.56           H
ATOM     13  D  BASP B  81     -12.058 -11.208  -2.025  0.69 17.73           D
ATOM     14  N  ASER B  82      -8.190 -12.497  -0.404  0.56 19.69           N
ATOM     15  CA ASER B  82      -7.139 -12.096   0.535  0.56 20.47           C
ATOM     16  C  ASER B  82      -6.167 -11.056  -0.035  0.56 17.24           C
ATOM     17  O  ASER B  82      -5.431 -10.452   0.743  0.56 18.12           O
ATOM     18  CB ASER B  82      -6.354 -13.336   1.009  0.56 26.37           C
ATOM     19  OG ASER B  82      -7.094 -14.060   1.996  0.56 28.69           O
ATOM     20  DG ASER B  82      -7.951 -14.366   1.633  0.56 24.69           D
ATOM     21  H  ASER B  82      -8.111 -13.457  -0.738  0.15 20.86           H
ATOM     22  HA ASER B  82      -7.612 -11.642   1.406  0.56 23.40           H
ATOM     23  HB2ASER B  82      -6.163 -13.985   0.160  0.56 24.54           H
ATOM     24  HB3ASER B  82      -5.408 -12.998   1.432  0.56 22.13           H
ATOM     25  N  BSER B  82      -8.182 -12.519  -0.468  0.44 19.76           N
ATOM     26  CA BSER B  82      -7.051 -12.227   0.424  0.44 21.62           C
ATOM     27  C  BSER B  82      -6.233 -11.030  -0.041  0.44 17.42           C
ATOM     28  O  BSER B  82      -5.646 -10.319   0.780  0.44 18.26           O
ATOM     29  CB BSER B  82      -6.119 -13.452   0.521  0.44 19.15           C
ATOM     30  OG BSER B  82      -5.800 -13.975  -0.771  0.44 24.96           O
ATOM     31  D  BSER B  82      -8.184 -13.456  -0.864  0.85 19.97           D
ATOM     32  DG BSER B  82      -5.982 -13.331  -1.480  0.44 22.97           D
ATOM     33  HA BSER B  82      -7.431 -12.005   1.420  0.44 23.61           H
ATOM     34  HB2BSER B  82      -5.197 -13.136   1.005  0.44 20.95           H
ATOM     35  HB3BSER B  82      -6.607 -14.225   1.116  0.44 23.29           H
ATOM     36  N   SER B  83      -6.173 -10.806  -1.349  1.00 11.89           N
ATOM     37  CA  SER B  83      -5.411  -9.695  -1.911  1.00 11.37           C
ATOM     38  C   SER B  83      -5.903  -8.334  -1.430  1.00 17.87           C
ATOM     39  O   SER B  83      -5.145  -7.360  -1.428  1.00 13.17           O
ATOM     40  CB  SER B  83      -5.528  -9.723  -3.434  1.00 18.80           C
ATOM     41  OG  SER B  83      -6.876  -9.521  -3.867  1.00 17.48           O
ATOM     42  HA  SER B  83      -4.365  -9.833  -1.634  1.00 19.36           H
ATOM     43  HB2 SER B  83      -4.912  -8.910  -3.820  1.00 23.88           H
ATOM     44  HB3 SER B  83      -5.158 -10.680  -3.813  1.00 15.26           H
ATOM     45  D  BSER B  83      -6.624 -11.353  -2.062  1.00 20.57           D
ATOM     46  DG DSER B  83      -7.377 -10.365  -3.873  1.00 22.34           D
END
"""

# A47 (open, not used by a test yet), 7T5E A 130-132: Lys A 131 HZ1-3 in A,
# DZ1/DZ2 in B; after the input normalization HZ1/HZ2 are blank, HZ3 is
# alone in A, the amide H in B. add keep_existing_H=True stops at the output
# H/D check ("HZ3C/DZ3D: 1 competing partner")
pdb_a47_7t5e_lys_a131 = """\
CRYST1   68.296   42.269   70.413  90.00  98.47  90.00 P 1 21 1
ATOM      1  N   GLY A 130     -13.378   7.004 -30.766  1.00 26.79           N
ATOM      2  CA  GLY A 130     -13.919   6.573 -32.040  1.00 22.20           C
ATOM      3  C   GLY A 130     -15.335   6.071 -31.961  1.00 28.48           C
ATOM      4  O   GLY A 130     -15.721   5.148 -32.691  1.00 28.25           O
ATOM      5  HA2 GLY A 130     -13.915   7.381 -32.767  1.00 23.99           H
ATOM      6  HA3 GLY A 130     -13.297   5.766 -32.435  1.00 30.22           H
ATOM      7  H  AGLY A 130     -13.185   6.225 -30.138  0.03 27.78           H
ATOM      8  D  BGLY A 130     -13.190   6.219 -30.151  0.97 28.12           D
ATOM      9  N   LYS A 131     -16.141   6.717 -31.124  1.00 23.22           N
ATOM     10  CA  LYS A 131     -17.548   6.357 -30.962  1.00 17.59           C
ATOM     11  C   LYS A 131     -17.756   4.897 -30.524  1.00 16.38           C
ATOM     12  O   LYS A 131     -18.812   4.306 -30.776  1.00 19.42           O
ATOM     13  CB  LYS A 131     -18.326   6.623 -32.256  1.00 35.34           C
ATOM     14  CG  LYS A 131     -18.337   8.088 -32.643  1.00 43.86           C
ATOM     15  CD  LYS A 131     -19.046   8.339 -33.976  1.00 48.35           C
ATOM     16  CE  LYS A 131     -19.374   9.821 -34.151  1.00 47.46           C
ATOM     17  NZ  LYS A 131     -20.670  10.246 -33.479  1.00 70.24           N
ATOM     18  HA  LYS A 131     -17.999   6.995 -30.208  1.00 25.58           H
ATOM     19  HB2 LYS A 131     -17.874   6.065 -33.075  1.00 24.06           H
ATOM     20  HB3 LYS A 131     -19.361   6.309 -32.090  1.00 30.72           H
ATOM     21  HG2 LYS A 131     -18.885   8.648 -31.882  1.00 38.19           H
ATOM     22  HG3 LYS A 131     -17.313   8.449 -32.732  1.00 41.48           H
ATOM     23  HD2 LYS A 131     -18.368   8.038 -34.773  1.00 36.31           H
ATOM     24  HD3 LYS A 131     -19.973   7.770 -34.017  1.00 37.34           H
ATOM     25  HE2 LYS A 131     -18.559  10.385 -33.695  1.00 42.17           H
ATOM     26  HE3 LYS A 131     -19.443  10.059 -35.213  1.00 54.16           H
ATOM     27  D  BLYS A 131     -15.836   7.525 -30.576  1.00 24.23           D
ATOM     28  DZ1BLYS A 131     -20.704   9.907 -32.496  0.14 45.83           D
ATOM     29  DZ2BLYS A 131     -20.735  11.285 -33.473  0.53 52.21           D
ATOM     30  HZ1ALYS A 131     -20.696   9.876 -32.507  1.00 45.46           H
ATOM     31  HZ2ALYS A 131     -20.749  11.283 -33.454  0.86 52.34           H
ATOM     32  HZ3ALYS A 131     -21.477   9.854 -34.001  0.47 51.02           H
ATOM     33  N   GLY A 132     -16.783   4.304 -29.835  1.00 19.85           N
ATOM     34  CA  GLY A 132     -16.974   2.967 -29.318  1.00 18.69           C
ATOM     35  C   GLY A 132     -16.764   1.878 -30.346  1.00 24.96           C
ATOM     36  O   GLY A 132     -17.070   0.711 -30.048  1.00 17.33           O
ATOM     37  HA2 GLY A 132     -16.284   2.772 -28.498  1.00 23.60           H
ATOM     38  HA3 GLY A 132     -17.993   2.873 -28.943  1.00 23.32           H
ATOM     39  H  AGLY A 132     -15.878   4.727 -29.642  0.68 20.26           H
ATOM     40  D  BGLY A 132     -15.859   4.696 -29.637  0.32 20.43           D
END
"""

# ------------------------------------------------------------------------------
# B16: repair_plan() and apply_repairs()

# an Ala with H, from ideal geometry; ala_fixture() builds the constructed
# cases from it
_ala = """\
ATOM      1  N   ALA A   1      -0.966   0.493   1.500  1.00 10.00           N
ATOM      2  CA  ALA A   1       0.257   0.418   0.692  1.00 10.00           C
ATOM      3  C   ALA A   1      -0.094   0.017  -0.716  1.00 10.00           C
ATOM      4  O   ALA A   1      -1.056  -0.682  -0.923  1.00 10.00           O
ATOM      5  CB  ALA A   1       1.204  -0.620   1.296  1.00 10.00           C
ATOM      6  OXT ALA A   1       0.661   0.439  -1.742  1.00 10.00           O
ATOM      7  H   ALA A   1      -1.383  -0.425   1.482  1.00 10.00           H
ATOM      8  H2  ALA A   1      -0.676   0.661   2.452  1.00 10.00           H
ATOM      9  HA  ALA A   1       0.746   1.392   0.682  1.00 10.00           H
ATOM     10  HB1 ALA A   1       1.459  -0.330   2.316  1.00 10.00           H
ATOM     11  HB2 ALA A   1       0.715  -1.594   1.307  1.00 10.00           H
ATOM     12  HB3 ALA A   1       2.111  -0.676   0.693  1.00 10.00           H
ATOM     13  HXT ALA A   1       0.435   0.182  -2.647  1.00 10.00           H
"""

def ala_fixture(rot=0.0, shift=0.01, occ_a=0.5, occ_b=0.5, d_length=1.09):
  """One Ala in conformers A/B, a blank DA d_length from the CA copies'
  mean along the ideal direction; the copies' other atoms rotated by
  -rot/+rot about an axis through CA perpendicular to that direction, or
  shifted by -shift/2 and +shift/2 along it. The copies' ideal D positions
  are then 2 * 1.09 * sin(rot) or shift apart."""
  h = iotbx.pdb.input(source_info=None, lines=_ala).construct_hierarchy()
  xyz = dict((a.name.strip(), matrix.col(a.xyz)) for a in h.atoms())
  ca = xyz["CA"]
  s = matrix.col((0, 0, 0))
  for n in ("N", "C", "CB"): s += (xyz[n] - ca).normalize()
  u = -s.normalize()
  axis = u.cross(xyz["N"] - ca).normalize()
  lines = ["CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1"]
  def atom(name, alt, p, occ, element):
    p = p + matrix.col((15, 15, 15))
    lines.append("ATOM  %5d %-4s%1sALA A   1    %8.3f%8.3f%8.3f%6.2f%6.2f"
      "          %2s" % (len(lines), (" " + name) if len(name) < 4 else name,
      alt, p[0], p[1], p[2], occ, 10.0, element))
  for alt, sign, occ in (("A", -1, occ_a), ("B", 1, occ_b)):
    r = axis.axis_and_angle_as_r3_rotation_matrix(sign * rot, deg=True)
    for a in h.atoms():
      n = a.name.strip()
      if n == "HA": continue
      atom(n, alt, ca + r * (xyz[n] - ca) + axis * (sign * shift / 2.0),
           occ, a.element.strip())
  atom("DA", " ", ca + u * d_length, 1.0, "D")
  return "\n".join(lines) + "\nEND\n"

def ala_d_on_parent():
  """ala_fixture() with the DA on CA A (zero bond length)"""
  lines = ala_fixture().splitlines()
  ca = [l for l in lines if l[12:17] == " CA A"][0]
  return "\n".join(l[:30] + ca[30:54] + l[54:] if l[12:16] == " DA " else l
                   for l in lines) + "\n"

def ala_conformer_ha(occ_a=0.6, occ_b=0.4, shift=0.3):
  """ala_fixture() in conformers A/B shift apart, HA on CA A at occ_a and
  DA on CA B at occ_b (each at its conformer's occupancy): a carbon-bound
  (c) site, no blank atom"""
  lines = ala_fixture(shift=shift, occ_a=occ_a, occ_b=occ_b).splitlines()
  xyz = {}
  for l in lines:
    if l.startswith("ATOM"):
      xyz[l[12:17]] = matrix.col((float(l[30:38]), float(l[38:46]),
                                  float(l[46:54])))
  mean = (xyz[" CA A"] + xyz[" CA B"]) / 2.0
  result = []
  for l in lines:
    if l[12:17] != " DA  ":
      result.append(l)
      continue
    for name, alt, occ, element in (("HA", "A", occ_a, "H"),
                                    ("DA", "B", occ_b, "D")):
      p = xyz[" DA  "] + xyz[" CA " + alt] - mean
      result.append("ATOM  %5d  %-3s%1sALA A   1    %8.3f%8.3f%8.3f%6.2f"
        "%6.2f          %2s" % (len(result), name, alt, p[0], p[1], p[2],
        occ, 10.0, element))
  return "\n".join(result) + "\n"

def floating_chain(pdb_str, chain_id="A"):
  """pdb_str after a chain with only a floating H (Ala 900), so that the
  chains chain_id of pdb_str are their second occurrence until it is
  removed"""
  lines = pdb_str.splitlines()
  i = [k for k, l in enumerate(lines) if l.startswith(("ATOM", "HETATM"))][0]
  return "\n".join(lines[:i] + [
    "ATOM      1  H   ALA %s 900      40.000  40.000  40.000  1.00 20.00"
    "           H" % chain_id, "TER"] + lines[i:]) + "\n"

def make_model_cif(pdb_str, cif_str):
  """make_model() with a restraint CIF"""
  import iotbx.cif
  inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  model = mmtbx.model.manager(model_input=inp, log=null_out(),
    restraint_objects=[("user.cif",
      iotbx.cif.reader(input_string=cif_str).model())])
  model.set_stop_for_unknowns(False)
  model.process(make_restraints=True,
    pdb_interpretation_params=reduce_hydrogen.
      get_reduce_pdb_interpretation_params(False))
  return model

def plan_of(model):
  c, m = candidates_of(model)
  return hd_exchange.repair_plan(model, c, m), c, m

def refused_reasons(plan):
  return sorted(set(plan.refused.values()))

def repaired_match(model, plan):
  """the check on the repaired hierarchy: (match, repaired hierarchy)"""
  h = model.get_hierarchy().deep_copy()
  origin = hd_exchange.apply_repairs(h, plan)
  m2 = model_from_hierarchy(h, model.crystal_symmetry())
  c = hd_exchange.find_hd_candidates(m2.get_hierarchy(),
    hd_exchange.h_parent_table(m2), model=m2)
  accepted = set((x.h_key, x.d_key) for x in plan.accept + plan.conformer) \
    or None
  return hd_exchange.collapse_check(c, accepted=accepted), h, origin

def check_clean(match):
  assert not match.problems, [p.reasons for p in match.problems]
  assert not [u for u in match.unresolved if u.action == "stop"]

def test_repair_split():
  # (a) accepted: one copy per conformer, occupancy shared, values copied
  for pdb_str in (ala_fixture(), ala_fixture(shift=0.08),
                  pdb_a_7t5e_glu_a140, pdb_a_7t5e_val_b210,
                  pdb_a_6tae_gln_a13):
    model = make_model(pdb_str)
    plan, c, m = plan_of(model)
    assert [u for u in m.unresolved if u.action == "stop"]
    assert len(plan.split) == 1 and not plan.refused, (plan.split,
      plan.refused)
    assert "ideal positions" in plan.split[0].reason
    before = [a for a in model.get_hierarchy().atoms()
              if hd_exchange.atom_key(a) == plan.split[0].key][0]
    match, h, origin = repaired_match(model, plan)
    check_clean(match)
    copies = [a for a in h.atoms() if a.name == before.name and
              a.parent().altloc != "" and
              a.element.strip() == before.element.strip()]
    assert len(copies) == len(plan.split[0].parent_keys), len(copies)
    assert approx_equal(sum(a.occ for a in copies), before.occ, eps=1.e-3)
    for a in copies:
      assert approx_equal(a.xyz, before.xyz, eps=1.e-3)
      assert approx_equal(a.b, before.b)
      assert origin[hd_exchange.atom_key(a)] == plan.split[0].key
  # the D was the only atom of its blank group: the group is gone
  model = make_model(pdb_a_6tae_gln_a13)
  plan, c, m = plan_of(model)
  match, h, origin = repaired_match(model, plan)
  gln = [rg for rg in h.residue_groups() if rg.resseq_as_int() == 13][0]
  assert "" not in [ag.altloc for ag in gln.atom_groups()]

def test_repair_split_refused():
  # (a) refused, each with its reason; the site stays a stop
  for pdb_str, reason in (
      (ala_fixture(rot=14), "ideal positions on the copies 0.527 A apart "
                            "(limit 0.10)"),
      (ala_fixture(shift=0.15), "ideal positions on the copies 0.150 A apart "
                                "(limit 0.10)"),
      (ala_fixture(occ_a=0, occ_b=0), "parent copies with total occupancy 0 "
                                      "(no basis for sharing it)"),
      (ala_fixture(rot=16, shift=0), "16.0 deg from the ideal direction")):
    model = make_model(pdb_str)
    plan, c, m = plan_of(model)
    assert not plan.split
    assert len(plan.refused) == 1, plan.refused
    assert list(plan.refused.values())[0].startswith(reason), plan.refused

def test_repair_split_user_cif():
  # review fix 1: the lengths of the model's own dictionary. GeoStd CA-HA
  # 0.97/1.09 A refuses a DA at 1.19 A; a CIF with 1.20 A accepts it
  import libtbx.load_env, os
  fn = libtbx.env.find_in_repositories(
    relative_path=os.path.join("chem_data", "geostd", "a", "data_ALA.cif"),
    test=os.path.isfile)
  if fn is None:
    print("Skipping test_repair_split_user_cif: no chem_data/geostd")
    return
  with open(fn) as f:
    cif = f.read()
  old = " ALA  CA  HA   single  0.970  0.020  1.090"
  assert cif.count(old) == 1
  cif = cif.replace(old, " ALA  CA  HA   single  1.200  0.020  1.200")
  pdb_str = ala_fixture(d_length=1.19)
  plan, c, m = plan_of(make_model(pdb_str))
  assert not plan.split
  assert refused_reasons(plan)[0].startswith("bond 1.190 A to CA AALA A   1"
    " (dictionary 0.97/1.09 A)"), refused_reasons(plan)
  model = make_model_cif(pdb_str, cif)
  plan, c, m = plan_of(model)
  assert len(plan.split) == 1 and not plan.refused, plan.refused
  check_clean(repaired_match(model, plan)[0])

def test_repair_rename():
  # (b) each H takes the name of the D at its site; a leftover H the one
  # name left free
  for pdb_str, n_renames, completed in ((pdb_b_7t5e_gln_a14, 2, 0),
      (pdb_b_7kcu_lys_a3, 5, 0), (pdb_b_7kcu_lys_a41, 3, 1)):
    model = make_model(pdb_str)
    plan, c, m = plan_of(model)
    assert len(plan.rename) == n_renames and not plan.refused, (
      [x.reason for x in plan.rename], plan.refused)
    reasons = [x.reason for x in plan.rename]
    assert len([r for r in reasons if "the one name left free" in r]) == \
      completed, reasons
    assert len([r for r in reasons if "the name of the D at its site" in r]) \
      == n_renames - completed, reasons
    match, h, origin = repaired_match(model, plan)
    check_clean(match)
    # names, not sites: every H/D at its input position and values
    before = sorted((a.xyz, a.occ, a.b, a.element) for a in
      model.get_hierarchy().atoms() if a.element.strip() in ("H", "D"))
    after = sorted((a.xyz, a.occ, a.b, a.element) for a in h.atoms()
      if a.element.strip() in ("H", "D"))
    assert before == after

def test_repair_rename_refused():
  # (the Lys B 125 cut also holds Asn B 126, whose amide is renamed)
  for pdb_str, resseq, reason in (
      (pdb_ctl_7kcu_lys_b125, 125, "1 leftover H, 0 free names"),
      (pdb_ctl_7kcu_gln_b4, 4, "no partner within 0.20 A")):
    plan, c, m = plan_of(make_model(pdb_str))
    assert not [x for x in plan.rename if x.key[3].strip() == str(resseq)]
    assert refused_reasons(plan) == [reason], refused_reasons(plan)

def test_repair_relabel():
  # (d) H and D moved into two free H/D-only letters, values unchanged
  model = make_model(pdb_d_6kk8_met_a41)
  plan, c, m = plan_of(model)
  assert len(plan.relabel) == 1 and not plan.refused, plan.refused
  assert "read as an exchange pair" in plan.relabel[0].reason
  assert "an inference" in plan.relabel[0].reason
  match, h, origin = repaired_match(model, plan)
  check_clean(match)
  rg = [rg for rg in h.residue_groups() if rg.resseq_as_int() == 41][0]
  hd = dict((a.element.strip(), a) for a in rg.atoms() if a.name.strip() in
            ("H", "D"))
  old = dict((a.element.strip(), a) for a in model.get_hierarchy().atoms()
    if a.parent().parent().resseq_as_int() == 41 and a.name.strip() in
    ("H", "D"))
  conformers = set(ag.altloc for ag in rg.atom_groups()
                   if not hd_exchange._hd_only(ag))
  for e in ("H", "D"):
    a = hd[e]
    assert a.parent().altloc not in conformers
    assert hd_exchange._hd_only(a.parent())
    assert (a.xyz, a.occ, a.b) == (old[e].xyz, old[e].occ, old[e].b)
  assert hd["H"].parent().altloc != hd["D"].parent().altloc

def test_repair_relabel_refused():
  # (d) guard: H and D each at its conformer's occupancy read as one isotope
  # per conformer as well
  h = iotbx.pdb.input(source_info=None,
    lines=pdb_d_6kk8_met_a41).construct_hierarchy()
  rg = [rg for rg in h.residue_groups() if rg.resseq_as_int() == 41][0]
  occ = dict((ag.altloc, ag.atoms()[0].occ) for ag in rg.atom_groups()
             if ag.altloc and not hd_exchange._hd_only(ag))
  for a in rg.atoms():
    if a.name.strip() in ("H", "D"): a.set_occ(occ[a.parent().altloc])
  pdb_str = h.as_pdb_string(crystal_symmetry=iotbx.pdb.input(
    source_info=None, lines=pdb_d_6kk8_met_a41).crystal_symmetry())
  plan, c, m = plan_of(make_model(pdb_str))
  assert not plan.relabel
  assert refused_reasons(plan)[0].startswith(
    "H and D occupancies equal their conformers'"), refused_reasons(plan)

def test_repair_accept():
  # (e) one-H parent: the pair collapses with accepted, the H kept
  model = make_model(pdb_e_7t5e_arg_a21)
  plan, c, m = plan_of(model)
  assert len(plan.accept) == 1 and not plan.refused, plan.refused
  assert "an empirical repair limit" in plan.accept[0].reason
  assert not hd_exchange.repairs_made(plan)
  accepted = set((x.h_key, x.d_key) for x in plan.accept)
  check_clean(hd_exchange.collapse_check(c, accepted=accepted))
  # beyond the limit: refused
  h = model.get_hierarchy().deep_copy()
  d = [a for a in h.atoms() if a.name.strip() == "DE" and
       a.parent().parent().resseq_as_int() == 21][0]
  ne = [a for a in h.atoms() if a.name.strip() == "NE" and
        a.parent().parent().resseq_as_int() == 21][0]
  hh = [a for a in h.atoms() if a.name.strip() == "HE" and
        a.parent().parent().resseq_as_int() == 21][0]
  # D rotated about NE away from the H until 0.7 A apart, bond kept
  axis = (matrix.col(hh.xyz) - matrix.col(ne.xyz)).cross(
    matrix.col(d.xyz) - matrix.col(ne.xyz)).normalize()
  v = matrix.col(hh.xyz) - matrix.col(ne.xyz)
  L = (matrix.col(d.xyz) - matrix.col(ne.xyz)).length()
  angle = 2 * math.degrees(math.asin(0.35 / L))
  d.xyz = (matrix.col(ne.xyz) + axis.axis_and_angle_as_r3_rotation_matrix(
    angle, deg=True) * v.normalize() * L).elems
  plan, c, m = plan_of(model_from_hierarchy(h, model.crystal_symmetry()))
  assert not plan.accept
  assert refused_reasons(plan) == [
    "H and D 0.70 A apart (limit 0.60 A, empirical)"], refused_reasons(plan)

def test_repair_accept_refused_multi_h():
  # (e) on a parent with two H sites: refused, the D's twist named
  plan, c, m = plan_of(make_model(pdb_ctl_7kcu_gln_a18))
  assert not plan.accept
  reasons = refused_reasons(plan)
  assert reasons, plan.refused
  for r in reasons:
    assert r.startswith("NE2 carries 2 H in the dictionary (one H site "
      "needed); D ") and r.endswith(" deg out of the CG-CD-NE2 plane"), r

def test_repair_conformer():
  # (c) one isotope per conformer (7T5E Tyr B 25): kept as two singles, no
  # atom changed; the different-parents problem goes with accepted
  model = make_model(pdb_c_7t5e_tyr_b25)
  plan, c, m = plan_of(model)
  assert [x.kind for x in m.problems] == ["different parents"]
  assert len(plan.conformer) == 1 and not plan.refused, plan.refused
  assert not (plan.split or plan.rename or plan.relabel or plan.accept)
  assert not hd_exchange.repairs_made(plan)
  x = plan.conformer[0]
  assert x.atom == "HH ATYR B  25 / DH BTYR B  25", x.atom
  assert x.reason == ("kept as two singles: HH 0.65 on OH ATYR B  25 "
    "(conformer 0.65), DH 0.35 on OH BTYR B  25 (conformer 0.35); read as "
    "one isotope per conformer (each at its conformer's occupancy), an "
    "inference"), x.reason
  match = hd_exchange.collapse_check(c, accepted=set([(x.h_key, x.d_key)]))
  check_clean(match)
  singles = set(s.atom_key for s in c.singles)
  assert x.h_key in singles and x.d_key in singles
  # the isotope the output writes in both conformers (step C review P3): an
  # exchangeable site D (H for all_h); a carbon-bound one (Ala HA in A, DA
  # in B) per neutron option
  for option, kind in (("hd_and_h", "D"), ("all_hd", "D"), ("all_h", "H")):
    assert hd_exchange.conformer_site_isotope(option, x) == kind, option
  plan, c, m = plan_of(make_model(ala_conformer_ha()))
  assert len(plan.conformer) == 1 and not plan.refused, plan.refused
  x = plan.conformer[0]
  assert x.atom == "HA AALA A   1 / DA BALA A   1", x.atom
  assert not x.exchangeable and not x.water
  for option, kind in (("hd_and_h", "H"), ("all_hd", "H"), ("hd_and_d", "D"),
                       ("all_d", "D"), ("all_h", "H")):
    assert hd_exchange.conformer_site_isotope(option, x) == kind, option
  # an H/D at another occupancy than its conformer's: refused
  for old, new, reason in (
      ("  0.65 24.15           H", "  0.50 24.15           H",
       "HH 0.50 on OH ATYR B  25 (conformer 0.65), DH 0.35 on OH BTYR B  25 "
       "(conformer 0.35): neither one isotope per conformer nor an exchange "
       "pair"),
      ("  0.65 24.15           H", "  0.15 24.15           H",
       "HH 0.15 on OH ATYR B  25 (conformer 0.65), DH 0.35 on OH BTYR B  25 "
       "(conformer 0.35): neither one isotope per conformer nor an exchange "
       "pair")):
    assert pdb_c_7t5e_tyr_b25.count(old) == 1
    plan, c, m = plan_of(make_model(pdb_c_7t5e_tyr_b25.replace(old, new)))
    assert not plan.conformer
    assert refused_reasons(plan) == [reason], refused_reasons(plan)

SER_B82_REASON = ("H 0.15 on N  ASER B  82 (conformer 0.56), D 0.85 on N  "
  "BSER B  82 (conformer 0.44), H + D = 1.00: consistent with an exchange "
  "pair over the two copies, but the occupancies do not determine how H and "
  "D are distributed over the conformers (0.56/0.44), so the site is not "
  "repaired; remove the H and D at this site and rerun with approach=add "
  "keep_existing_H=True, which places them anew (approach=optimize adds no "
  "H)")

def test_repair_conformer_refused():
  # Ser B 82: H + D = 1 over the two copies of the split N: consistent with
  # an exchange pair, but how H and D sit in the conformers is not known
  plan, c, m = plan_of(make_model(pdb_ctl_7t5e_ser_b82))
  assert not plan.conformer
  assert refused_reasons(plan) == [SER_B82_REASON], refused_reasons(plan)

def test_repair_split_on_parent():
  # step B review fix 3: a DA on CA A has no direction; refused with the
  # bond length, no TypeError
  plan, c, m = plan_of(make_model(ala_d_on_parent()))
  assert not plan.split
  assert refused_reasons(plan) == ["bond 0.000 A to CA AALA A   1 "
    "(dictionary 0.97/1.09 A)"], refused_reasons(plan)

def test_collapse_check_accepted_default():
  # accepted=None (validate_H, the output check) is collapse_check() as
  # before: same pairs and problems
  model = make_model(pdb_e_7t5e_arg_a21)
  c, m = candidates_of(model)
  m0 = hd_exchange.collapse_check(c)
  m1 = hd_exchange.collapse_check(c, accepted=None)
  m2 = hd_exchange.collapse_check(c, accepted=set())
  for x in (m1, m2):
    assert [p.h_id for p in x.pairs] == [p.h_id for p in m0.pairs]
    assert [p.reasons for p in x.problems] == [p.reasons for p in m0.problems]
  assert [p.reasons for p in m0.problems] == [["H and D 0.23 A apart"]]

if __name__ == "__main__":
  test_parent_table_30f()
  test_edits_covalent_limit()
  test_hydrogen_bond_origin_is_not_a_parent()
  test_symmetry()
  test_hd_pair_keys()
  test_name_normalization()
  test_repeated_labels()
  test_duplicate_is_an_error()
  test_water_canonical_resname()
  test_keys_survive_copies()
  test_basic_pair()
  test_pair_0p5_apart()
  test_swapped_names()
  test_one_h_two_d()
  test_v2_names()
  test_pair_in_one_conformer()
  test_layout_needs_two_atom_groups()
  test_renamed_d_collides()
  test_notes()
  test_microheterogeneity()
  test_one_isotope_per_site()
  test_ambiguous_is_unresolved()
  test_symmetry_parent()
  test_different_parents_5A()
  test_reversed_order()
  test_alternative_h_orientations()
  test_modes()
  test_altloc_groups()
  test_altloc_groups_conformers()
  test_merge_conformer_copies()
  test_conformer_copies_differ()
  test_no_free_altloc_letters()
  test_name_not_pairable()
  test_parent_occupancy_from_table()
  test_anisotropic()
  test_round_trip()
  test_add_water_h_parents()
  test_30f_split()
  test_alternative_orientations_output()
  test_cross_mode()
  test_water_cross_mode()
  test_water_alternative_orientations()
  test_lone_alternate_h()
  test_two_models()
  test_repeated_labels_restore()
  test_water_names_follow_content()
  test_conformer_copy_occupancy_note()
  test_rekey_records()
  test_repair_split()
  test_repair_split_refused()
  test_repair_split_user_cif()
  test_repair_rename()
  test_repair_rename_refused()
  test_repair_relabel()
  test_repair_relabel_refused()
  test_repair_accept()
  test_repair_accept_refused_multi_h()
  test_repair_conformer()
  test_repair_conformer_refused()
  test_repair_split_on_parent()
  test_collapse_check_accepted_default()
  print("OK")
