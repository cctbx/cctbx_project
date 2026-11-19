import iotbx.pdb
import mmtbx.model
import mmtbx.f_model
from mmtbx.validation import geo_vs_map
from scitbx.array_family import flex
from scitbx.matrix import col
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out
from libtbx import easy_pickle

# Water goes first and there is a H atom so that atom indices of the atoms
# used differ from indices in the model.
pdb_str = """
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1
HETATM    1  O   HOH S   1       4.000   4.000   4.000  1.00 20.00           O
ATOM      2  N   GLY A   1       6.011   8.357   7.697  1.00 10.00           N
ATOM      3  CA  GLY A   1       7.303   7.697   7.846  1.00 10.00           C
ATOM      4  C   GLY A   1       8.456   8.672   7.697  1.00 10.00           C
ATOM      5  O   GLY A   1       8.272   9.884   7.580  1.00 10.00           O
ATOM      6  N   ALA A   2       9.660   8.117   7.703  1.00 10.00           N
ATOM      7  H   ALA A   2       9.737   7.123   7.800  1.00 10.00           H
ATOM      8  CA  ALA A   2      10.874   8.925   7.571  1.00 15.00           C
ATOM      9  C   ALA A   2      12.103   8.031   7.694  1.00 15.00           C
ATOM     10  O   ALA A   2      12.000   6.811   7.818  1.00 15.00           O
ATOM     11  CB  ALA A   2      10.897   9.670   6.244  1.00 25.00           C
ATOM     12  N   GLY A   3      13.262   8.672   7.656  1.00 20.00           N
ATOM     13  CA  GLY A   3      14.531   7.950   7.759  1.00 20.00           C
ATOM     14  C   GLY A   3      15.697   8.920   7.697  1.00 20.00           C
ATOM     15  O   GLY A   3      15.500  10.130   7.600  1.00 30.00           O
ATOM     16  OXT GLY A   3      16.850   8.450   7.750  1.00 30.00           O
TER
END
"""

def get_model_and_fmodel(lines = pdb_str):
  pdb_inp = iotbx.pdb.input(source_info=None, lines=lines)
  model = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  model.process(make_restraints=True)
  xrs = model.get_xray_structure()
  f_obs = abs(xrs.structure_factors(d_min=1.0).f_calc())
  fmodel = mmtbx.f_model.manager(
    f_obs          = f_obs,
    r_free_flags   = f_obs.generate_r_free_flags(),
    xray_structure = xrs)
  return model, fmodel

def exercise_rank():
  assert geo_vs_map.rank_0_1(flex.double()).size() == 0
  assert approx_equal(geo_vs_map.rank_0_1(flex.double([5])), [1])
  assert approx_equal(
    geo_vs_map.rank_0_1(flex.double([3, -1, 10, 7])), [1/3, 0, 1, 2/3])
  # equal values: average rank
  assert approx_equal(
    geo_vs_map.rank_0_1(flex.double([2, 1, 2, 3])), [0.5, 0, 0.5, 1])

def exercise_map_values():
  model, fmodel = get_model_and_fmodel()
  sites_frac = model.get_sites_frac()
  v1 = geo_vs_map.map_values(fmodel = fmodel, sites_frac = sites_frac)
  v2 = geo_vs_map.map_values(fmodel = fmodel, sites_frac = sites_frac,
    sigma_scale = True)
  assert v1.size() == v2.size() == model.size()
  # atoms are in density
  not_h = ~model.get_hd_selection()
  assert v2.select(not_h).all_gt(3)
  # sigma scaling does not change the order
  assert approx_equal(geo_vs_map.rank_0_1(v1), geo_vs_map.rank_0_1(v2))
  # model is the source of data: difference map is flat compared to the atoms
  v3 = geo_vs_map.map_values(fmodel = fmodel, sites_frac = sites_frac,
    map_type = "mFo-DFmodel")
  assert flex.max(flex.abs(v3)) < flex.min(v1.select(not_h)) / 100

def check_restraints(restraints, n_atoms, model, atoms, value):
  sites_cart = model.get_sites_cart()
  used = list(atoms.i_seqs)
  n = len(restraints.i_seqs)
  assert n > 0
  for key in ["ideal", "model", "delta", "weight", "origin_id",
              "map_2mfo_dfc", "map_mfo_dfc", "rank"]:
    assert len(getattr(restraints, key)) == n, key
  assert approx_equal(restraints.delta, restraints.model - restraints.ideal)
  assert restraints.weight.all_gt(0)
  for i in range(n):
    i_seqs = restraints.i_seqs[i]
    assert len(i_seqs) == n_atoms
    # indices are those of the model given, and of atoms used only
    for i_seq in i_seqs: assert i_seq in used
    sites = [col(sites_cart[i_seq]) for i_seq in i_seqs]
    assert approx_equal(value(sites), restraints.model[i], 1.e-4)
    # map values are those of the atoms of the restraint
    for key in ["map_2mfo_dfc", "map_mfo_dfc", "rank"]:
      expected = [getattr(atoms, key)[used.index(i_seq)] for i_seq in i_seqs]
      assert approx_equal(getattr(restraints, key)[i], expected)

def check_per_atom(result, atoms):
  used = list(atoms.i_seqs)
  # mean absolute deviation per atom, atoms in order of first appearance
  deltas = {}
  for i_seqs, delta in zip(result.restraints.i_seqs, result.restraints.delta):
    for i_seq in i_seqs: deltas.setdefault(i_seq, []).append(abs(delta))
  assert result.deltas.size() == result.maps.size() == len(deltas)
  for i, i_seq in enumerate(deltas.keys()):
    assert approx_equal(result.deltas[i], flex.mean(flex.double(deltas[i_seq])))
    assert approx_equal(result.maps[i], atoms.rank[used.index(i_seq)])
  # bins
  assert len(result.bins) == 5
  assert approx_equal(result.bins[0], [0, 0.2])
  assert approx_equal(result.bins[-1], [0.8, 1])
  assert result.x.size() == result.y.size() > 0
  for x, y in zip(result.x, result.y):
    sel = (result.maps >= x - 0.2) & (result.maps <= x)
    assert approx_equal(y, flex.mean(result.deltas.select(sel)))
  assert result.d_min == "1.0"

def exercise_run():
  model, fmodel = get_model_and_fmodel()
  sites_cart = model.get_sites_cart().deep_copy()
  result = geo_vs_map.run(model = model, fmodel = fmodel)
  # model is not changed
  assert model.size() == 16
  assert approx_equal(model.get_sites_cart(), sites_cart)
  # atoms used: no water, no H
  atoms = result.atoms
  # (H is placed after CB of its residue when the model is made)
  assert list(atoms.i_seqs) == [1,2,3,4,5,6,7,8,9,11,12,13,14,15]
  model_atoms = model.get_hierarchy().atoms()
  for i_seq, id_str in zip(atoms.i_seqs, atoms.id_strs):
    assert model_atoms[i_seq].id_str() == id_str
  for key in ["map_2mfo_dfc", "map_mfo_dfc", "rank"]:
    assert getattr(atoms, key).size() == 14
  assert approx_equal(atoms.rank, geo_vs_map.rank_0_1(atoms.map_2mfo_dfc))
  assert approx_equal(flex.min(atoms.rank), 0)
  assert approx_equal(flex.max(atoms.rank), 1)
  #
  check_restraints(result.bonds.restraints, 2, model, atoms,
    value = lambda s: abs(s[0]-s[1]))
  check_restraints(result.angles.restraints, 3, model, atoms,
    value = lambda s: (s[0]-s[1]).angle(s[2]-s[1], deg=True))
  check_per_atom(result.bonds, atoms)
  check_per_atom(result.angles, atoms)
  # no bonds to H, all 13 bonds of 14 atoms
  assert len(result.bonds.restraints.i_seqs) == 13
  #
  assert approx_equal(result.r_work, fmodel.r_work())
  assert approx_equal(result.r_free, fmodel.r_free())
  assert result.g is not None
  # can be stored
  easy_pickle.dump("tst_geo_vs_map.pkl", result)
  loaded = easy_pickle.load("tst_geo_vs_map.pkl")
  assert approx_equal(loaded.bonds.deltas, result.bonds.deltas)
  assert loaded.bonds.restraints.i_seqs == result.bonds.restraints.i_seqs

def exercise_no_macromolecule():
  lines = "\n".join([l for l in pdb_str.splitlines()
    if l.startswith("CRYST1") or l.startswith("HETATM")])
  model, fmodel = get_model_and_fmodel(lines = lines)
  assert geo_vs_map.run(model = model, fmodel = fmodel) is None

if (__name__ == "__main__"):
  exercise_rank()
  exercise_map_values()
  exercise_run()
  exercise_no_macromolecule()
  print("OK")
