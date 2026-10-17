"""hd_exchange: atom/site identity and H parent table (no reduce2)."""
from __future__ import absolute_import, division, print_function
import iotbx.pdb
import iotbx.phil
import mmtbx.model
from libtbx import group_args
from libtbx.utils import null_out, Sorry
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
  print("OK")
