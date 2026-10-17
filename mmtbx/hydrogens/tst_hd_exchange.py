"""hd_exchange: atom/site identity and H parent table (no reduce2)."""
from __future__ import absolute_import, division, print_function
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

def candidates_of(model, xyz_tolerance=0.2):
  table = hd_exchange.h_parent_table(model)
  c = hd_exchange.find_hd_candidates(model.get_hierarchy(), table, model=model)
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
  assert r.pairs == [] and len(r.problems) == 2
  for prob in r.problems:
    assert any(x.startswith("names swapped") for x in prob.reasons)

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
  # placement copies every H of Ser into A and B; H and HA (blank parents,
  # coincident copies) are merged into the blank group, HB2/HB3/HG differ
  model = placed(pdb_pep_og_ab)
  hd = ser_h_d(model.get_hierarchy())
  assert hd[""] == [] and hd["A"] == hd["B"] == ["H", "HA", "HB2", "HB3", "HG"]
  h, r = run_option(model, "hd_and_h")
  hd = ser_h_d(h)
  assert hd[""] == ["HA"], hd                      # merged, not exchangeable
  assert hd["C"] == ["H"] and hd["D"] == ["D"], hd # merged and split
  assert hd["A"] == hd["B"] == ["DG", "HB2", "HB3"], hd
  assert sorted(names(x.site) for x in r.merged_sites) == ["H", "HA"]
  # exchangeable output sites (not water): Ala H1 H2 H3, Ser H once (merged
  # from two copies), HG in A and in B, Lys H HZ1 HZ2 HZ3
  assert r.n_exchangeable == 10, r.n_exchangeable
  assert all(x.reason == "merged conformer copies" for x in r.merged_sites)
  assert all(len(x.copies) == 2 for x in r.merged_sites)
  assert sorted((names(x.site), x.reason) for x in r.conformer_sites) == [
    ("HG", "conformer"), ("HG", "conformer")]
  # merged H at the parent occupancy; the split pair 0.5 each
  assert find(h, "SER", "HA").occ == 1.0
  assert find(h, "SER", "H").occ == 0.5 and find(h, "SER", "D").occ == 0.5
  # the split pair maps to the first (A) copy
  src = [r.provenance[hd_exchange.atom_key(find(h, "SER", n))]
         for n in ("H", "D")]
  assert src[0] == src[1] and src[0][-3:] == ("A", " H  ", "H"), src

def test_conformer_copies_differ():
  # a backbone H whose copies are 0.1 A apart is not merged: conformer site
  model = placed(pdb_pep_og_ab)
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
  print("OK")
