from __future__ import absolute_import, division, print_function
import iotbx.pdb
import mmtbx.model
from libtbx.test_utils import approx_equal
from mmtbx.monomer_library import server, pdb_interpretation

# Ideal poly-Ala helix with riding hydrogens; H(i+4)...O(i) = 2.1 A.
pdb_str = """\
CRYST1   18.325   16.888   20.596  90.00  90.00  90.00 P 1
ATOM      1  N   ALA A   1      -5.606  -2.251 -12.878  1.00  0.00           N
ATOM      2  CA  ALA A   1      -5.850  -1.194 -13.852  1.00  0.00           C
ATOM      3  C   ALA A   1      -5.186  -1.524 -15.184  1.00  0.00           C
ATOM      4  O   ALA A   1      -5.744  -1.260 -16.249  1.00  0.00           O
ATOM      5  CB  ALA A   1      -5.339   0.137 -13.324  1.00  0.00           C
ATOM      6  H   ALA A   1      -6.119  -2.935 -12.969  1.00  0.00           H
ATOM      7  HA  ALA A   1      -6.805  -1.113 -14.000  1.00  0.00           H
ATOM      8  HB1 ALA A   1      -5.803   0.348 -12.499  1.00  0.00           H
ATOM      9  HB3 ALA A   1      -5.510   0.826 -13.985  1.00  0.00           H
ATOM     10  HB2 ALA A   1      -4.386   0.066 -13.159  1.00  0.00           H
ATOM     11  N   ALA A   2      -3.992  -2.102 -15.115  1.00  0.00           N
ATOM     12  CA  ALA A   2      -3.261  -2.499 -16.313  1.00  0.00           C
ATOM     13  C   ALA A   2      -3.961  -3.660 -17.011  1.00  0.00           C
ATOM     14  O   ALA A   2      -4.016  -3.716 -18.240  1.00  0.00           O
ATOM     15  CB  ALA A   2      -1.833  -2.881 -15.957  1.00  0.00           C
ATOM     16  H   ALA A   2      -3.580  -2.277 -14.381  1.00  0.00           H
ATOM     17  HA  ALA A   2      -3.228  -1.750 -16.928  1.00  0.00           H
ATOM     18  HB1 ALA A   2      -1.365  -3.141 -16.766  1.00  0.00           H
ATOM     19  HB3 ALA A   2      -1.851  -3.623 -15.332  1.00  0.00           H
ATOM     20  HB2 ALA A   2      -1.393  -2.117 -15.552  1.00  0.00           H
ATOM     21  N   ALA A   3      -4.492  -4.585 -16.219  1.00  0.00           N
ATOM     22  CA  ALA A   3      -5.216  -5.731 -16.755  1.00  0.00           C
ATOM     23  C   ALA A   3      -6.531  -5.289 -17.389  1.00  0.00           C
ATOM     24  O   ALA A   3      -6.939  -5.814 -18.425  1.00  0.00           O
ATOM     25  CB  ALA A   3      -5.476  -6.751 -15.658  1.00  0.00           C
ATOM     26  H   ALA A   3      -4.446  -4.571 -15.360  1.00  0.00           H
ATOM     27  HA  ALA A   3      -4.676  -6.156 -17.439  1.00  0.00           H
ATOM     28  HB1 ALA A   3      -4.627  -7.049 -15.296  1.00  0.00           H
ATOM     29  HB3 ALA A   3      -5.958  -7.504 -16.035  1.00  0.00           H
ATOM     30  HB2 ALA A   3      -6.006  -6.336 -14.960  1.00  0.00           H
ATOM     31  N   ALA A   4      -7.189  -4.323 -16.758  1.00  0.00           N
ATOM     32  CA  ALA A   4      -8.442  -3.785 -17.273  1.00  0.00           C
ATOM     33  C   ALA A   4      -8.205  -3.003 -18.561  1.00  0.00           C
ATOM     34  O   ALA A   4      -9.007  -3.065 -19.492  1.00  0.00           O
ATOM     35  CB  ALA A   4      -9.102  -2.897 -16.231  1.00  0.00           C
ATOM     36  H   ALA A   4      -6.928  -3.959 -16.024  1.00  0.00           H
ATOM     37  HA  ALA A   4      -9.045  -4.519 -17.470  1.00  0.00           H
ATOM     38  HB1 ALA A   4      -9.282  -3.423 -15.436  1.00  0.00           H
ATOM     39  HB3 ALA A   4      -9.932  -2.550 -16.592  1.00  0.00           H
ATOM     40  HB2 ALA A   4      -8.503  -2.165 -16.015  1.00  0.00           H
ATOM     41  N   ALA A   5      -7.099  -2.269 -18.604  1.00  0.00           N
ATOM     42  CA  ALA A   5      -6.735  -1.498 -19.787  1.00  0.00           C
ATOM     43  C   ALA A   5      -6.358  -2.423 -20.939  1.00  0.00           C
ATOM     44  O   ALA A   5      -6.687  -2.157 -22.094  1.00  0.00           O
ATOM     45  CB  ALA A   5      -5.585  -0.556 -19.469  1.00  0.00           C
ATOM     46  H   ALA A   5      -6.538  -2.200 -17.956  1.00  0.00           H
ATOM     47  HA  ALA A   5      -7.495  -0.963 -20.063  1.00  0.00           H
ATOM     48  HB1 ALA A   5      -5.859   0.050 -18.763  1.00  0.00           H
ATOM     49  HB3 ALA A   5      -5.359  -0.053 -20.267  1.00  0.00           H
ATOM     50  HB2 ALA A   5      -4.821  -1.078 -19.179  1.00  0.00           H
ATOM     51  N   ALA A   6      -5.665  -3.509 -20.614  1.00  0.00           N
ATOM     52  CA  ALA A   6      -5.268  -4.493 -21.614  1.00  0.00           C
ATOM     53  C   ALA A   6      -6.485  -5.236 -22.153  1.00  0.00           C
ATOM     54  O   ALA A   6      -6.565  -5.533 -23.345  1.00  0.00           O
ATOM     55  CB  ALA A   6      -4.272  -5.476 -21.019  1.00  0.00           C
ATOM     56  H   ALA A   6      -5.411  -3.700 -19.815  1.00  0.00           H
ATOM     57  HA  ALA A   6      -4.837  -4.037 -22.354  1.00  0.00           H
ATOM     58  HB1 ALA A   6      -4.021  -6.121 -21.698  1.00  0.00           H
ATOM     59  HB3 ALA A   6      -4.686  -5.930 -20.269  1.00  0.00           H
ATOM     60  HB2 ALA A   6      -3.488  -4.989 -20.719  1.00  0.00           H
ATOM     61  N   ALA A   7      -7.430  -5.532 -21.267  1.00  0.00           N
ATOM     62  CA  ALA A   7      -8.660  -6.212 -21.655  1.00  0.00           C
ATOM     63  C   ALA A   7      -9.529  -5.303 -22.518  1.00  0.00           C
ATOM     64  O   ALA A   7     -10.158  -5.756 -23.474  1.00  0.00           O
ATOM     65  CB  ALA A   7      -9.427  -6.658 -20.421  1.00  0.00           C
ATOM     66  H   ALA A   7      -7.382  -5.348 -20.428  1.00  0.00           H
ATOM     67  HA  ALA A   7      -8.437  -7.001 -22.173  1.00  0.00           H
ATOM     68  HB1 ALA A   7      -8.873  -7.266 -19.907  1.00  0.00           H
ATOM     69  HB3 ALA A   7     -10.240  -7.107 -20.700  1.00  0.00           H
ATOM     70  HB2 ALA A   7      -9.646  -5.879 -19.887  1.00  0.00           H
TER
END
"""

def make_params(shrink_h, shrink_da=0, h_missing=False):
  return pdb_interpretation.ener_lib_as_nonbonded_params(
    ener_lib                     = server.ener_lib(),
    assume_hydrogens_all_missing = h_missing,
    factor_1_4_interactions      = 0.8,
    default_distance             = 1,
    minimum_distance             = 1,
    const_shrink_donor_acceptor  = shrink_da,
    const_shrink_hydrogen_acceptor = shrink_h)

def table_lookup(params, t1, t2):
  for a, b in [(t1, t2), (t2, t1)]:
    if a in params.distance_table and b in params.distance_table[a]:
      return params.distance_table[a][b]
  return None

def exercise_params_table():
  p = make_params(shrink_h=0.6)
  # amide H ... carbonyl O: 1.05 + 1.40 - 0.6
  assert approx_equal(table_lookup(p, "HNH1", "O"), 1.85)
  # water H ... water O (hb_type B)
  assert approx_equal(table_lookup(p, "HOH2", "OH2"), 1.85)
  # hydroxyl H ... ring N acceptor: 1.05 + 1.55 - 0.6
  assert approx_equal(table_lookup(p, "HOH1", "NR56"), 2.00)
  # heavy donor ... acceptor untouched (still from radius sum, no table entry)
  assert table_lookup(p, "NH1", "O") is None
  assert table_lookup(p, "OH1", "O") is None
  # nonpolar H untouched
  assert table_lookup(p, "HCH1", "O") is None
  # polar H against non-acceptor untouched
  assert table_lookup(p, "HNH1", "CH1") is None
  assert table_lookup(p, "HNH1", "NH1") is None
  assert table_lookup(p, "HNH1", "HNH1") is None
  # shrink 0 restores the old behaviour: no table entries at all
  p0 = make_params(shrink_h=0)
  assert len(list(p0.distance_table.keys())) == 0
  # no-H model (vdwh radii): entries are harmless but use the vdwh radii
  ph = make_params(shrink_h=0.6, h_missing=True)
  assert approx_equal(table_lookup(ph, "HNH1", "O"), 1.05 + 1.52 - 0.6)

def get_vdw(model, name_i, resseq_i, name_j, resseq_j):
  grm = model.get_restraints_manager().geometry
  sites_cart = model.get_sites_cart()
  pp = grm.pair_proxies(sites_cart=sites_cart)
  atoms = model.get_hierarchy().atoms()
  want = set()
  for a in atoms:
    if a.name.strip() == name_i and a.parent().parent().resseq_as_int() == resseq_i:
      want.add(("i", a.i_seq))
    if a.name.strip() == name_j and a.parent().parent().resseq_as_int() == resseq_j:
      want.add(("j", a.i_seq))
  i_seq = [s for k, s in want if k == "i"][0]
  j_seq = [s for k, s in want if k == "j"][0]
  result = []
  for proxy in pp.nonbonded_proxies.simple:
    if set(proxy.i_seqs) == set([i_seq, j_seq]):
      result.append(proxy.vdw_distance)
  for proxy in pp.nonbonded_proxies.asu:
    if set([proxy.i_seq, proxy.j_seq]) == set([i_seq, j_seq]):
      result.append(proxy.vdw_distance)
  assert len(result) == 1, (name_i, resseq_i, name_j, resseq_j, result)
  return result[0]

def get_model(pdb_str, shrink_h=None, strip_h=False):
  params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  if shrink_h is not None:
    params.pdb_interpretation.const_shrink_hydrogen_acceptor = shrink_h
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  model = mmtbx.model.manager(model_input=pdb_inp, log=None)
  if strip_h:
    model = model.select(~model.get_hd_selection())
  model.process(make_restraints=True, pdb_interpretation_params=params)
  return model

def exercise_model_default():
  # default: polar H ... acceptor pairs are allowed at H-bond distance
  model = get_model(pdb_str)
  assert approx_equal(get_vdw(model, "H", 5, "O", 1), 1.85)
  # heavy atoms of the same H-bond unchanged (1.55 + 1.40)
  assert approx_equal(get_vdw(model, "N", 5, "O", 1), 2.95)
  # nonpolar H against O unchanged (1.22 + 1.40)
  assert approx_equal(get_vdw(model, "HB1", 5, "O", 1), 2.62)
  # no nonbonded violation in the ideal helix H-bonds
  grm = model.get_restraints_manager().geometry
  sites_cart = model.get_sites_cart()
  pp = grm.pair_proxies(sites_cart=sites_cart)
  atoms = model.get_hierarchy().atoms()
  for proxy in pp.nonbonded_proxies.simple:
    i, j = proxy.i_seqs
    names = set([atoms[i].name.strip(), atoms[j].name.strip()])
    if names == set(["H", "O"]):
      d = atoms[i].distance(atoms[j])
      assert d > proxy.vdw_distance, (atoms[i].id_str(), atoms[j].id_str(), d, proxy.vdw_distance)

def exercise_model_shrink_zero():
  # old behaviour on request: H-bond H...O scored against the full radius sum
  model = get_model(pdb_str, shrink_h=0)
  assert approx_equal(get_vdw(model, "H", 5, "O", 1), 2.45)
  assert approx_equal(get_vdw(model, "N", 5, "O", 1), 2.95)

def exercise_model_no_h():
  # a model without H is not affected (united-atom radii 1.60 + 1.52)
  for shrink_h in [None, 0]:
    model = get_model(pdb_str, shrink_h=shrink_h, strip_h=True)
    assert approx_equal(get_vdw(model, "N", 5, "O", 1), 3.12)

def run():
  exercise_params_table()
  exercise_model_default()
  exercise_model_shrink_zero()
  exercise_model_no_h()
  print("OK")

if (__name__ == "__main__"):
  run()
