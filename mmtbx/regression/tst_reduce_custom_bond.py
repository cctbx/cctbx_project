"""Regression tests for issue cctbx/cctbx_project#1199.

reduce2's Optimizer decides Asn/Gln/His flips (and the ion "lock-down" that
suppresses them) from local geometry and element type only -- a spatial
neighbour search plus element_is_positive_ion().  It did not consult the model's
bond restraints, so a user-defined (custom) bond on a flipping atom was invisible
to the flip logic and could be silently broken.

The Histidine ring is the exposed case.  MoverHisFlip validates the ring
nitrogens ND1/NE2 by counting hydrogens (== 1) and carbons (== 2), not their
total number of bonded neighbours, so an extra bond -- to a metal, or here to a
Cys SG -- slips through, the Mover is built and the ring can be flipped into an
orientation that makes the bond geometrically impossible.  (The ring carbons use
strict neighbour-count checks and are protected; only the two ring nitrogens
slip through.  The Asn/Gln amide is protected because MoverAmideFlip requires an
exact neighbour count and raises otherwise.)

Tests:
- test_his_is_flippable_without_custom_bond: control; the His is flippable.
- test_custom_bond_on_his_ring_n_suppresses_flip: a user-defined bond on ND1
  must suppress the flip (the bug fixed in _PlaceMovers by checking the bonded
  neighbours before adding a flip Mover).
- test_metal_coordinated_his_stays_locked: guards the fix -- a His whose ND1 is
  genuinely bonded to a nearby Zn must still be locked in place by the ion
  lock-down, not skipped.  This fails if the bonded-neighbour check is placed
  before the lock-down instead of after it.
- test_deleted_h_in_model_order: getHydrogensToDelete() returns the H in model
  order. It returned a set's order, which follows memory addresses, so
  reduce2's "Deleting ..." lines changed order from run to run.
"""
from __future__ import absolute_import, division, print_function
from iotbx.data_manager import DataManager
from mmtbx.hydrogens import reduce_hydrogen
from mmtbx.reduce import Optimizers
from mmtbx.probe import Helpers
from cctbx import geometry_restraints
from libtbx.utils import null_out
import mmtbx.model

# His A 41 from 1xso (a metal-free, flippable histidine) plus a cysteine whose
# side chain has been translated so that its SG sits ~2.3 A from His 41 ND1 -- a
# plausible distance for a user-defined bond.  Only CRYST1 is needed; the
# DataManager fills in the rest.
his_cys_pdb = """\
CRYST1   73.450   68.940   58.760  90.00  90.00  90.00 P 21 21 21    8
ATOM      1  N   HIS A  41      15.717  35.160  -4.931  1.00 10.07           N
ATOM      2  CA  HIS A  41      15.914  34.644  -3.612  1.00 10.05           C
ATOM      3  C   HIS A  41      17.386  34.282  -3.367  1.00  8.36           C
ATOM      4  O   HIS A  41      18.015  33.572  -4.158  1.00 10.78           O
ATOM      5  CB  HIS A  41      15.089  33.373  -3.411  1.00 11.44           C
ATOM      6  CG  HIS A  41      13.587  33.619  -3.505  1.00 11.15           C
ATOM      7  ND1 HIS A  41      12.813  33.864  -2.376  1.00 13.52           N
ATOM      8  CD2 HIS A  41      12.741  33.657  -4.566  1.00 13.29           C
ATOM      9  CE1 HIS A  41      11.563  34.039  -2.768  1.00 16.05           C
ATOM     10  NE2 HIS A  41      11.504  33.918  -4.071  1.00 14.57           N
ATOM     11  N   CYS B   1      13.821  33.913   2.939  1.00 12.92           N
ATOM     12  CA  CYS B   1      14.913  34.541   2.199  1.00 11.46           C
ATOM     13  C   CYS B   1      15.676  35.564   3.004  1.00 13.18           C
ATOM     14  O   CYS B   1      16.758  35.999   2.616  1.00 13.08           O
ATOM     15  CB  CYS B   1      14.388  35.150   0.916  1.00 12.38           C
ATOM     16  SG  CYS B   1      13.558  33.954  -0.202  1.00 13.46           S
TER
END
"""

# His A 69 and its Zn from 1xso: ND1 coordinates the Zn at ~2.07 A, so the ion
# lock-down should fire and leave the ring in place (state 0).
his_zn_pdb = """\
CRYST1   73.450   68.940   58.760  90.00  90.00  90.00 P 21 21 21    8
ATOM      1  N   HIS A  69      30.438  39.998   3.098  1.00  9.85           N
ATOM      2  CA  HIS A  69      29.461  39.439   2.195  1.00  8.06           C
ATOM      3  C   HIS A  69      29.530  40.119   0.822  1.00 10.28           C
ATOM      4  O   HIS A  69      29.676  41.356   0.727  1.00 11.99           O
ATOM      5  CB  HIS A  69      28.065  39.593   2.817  1.00 11.09           C
ATOM      6  CG  HIS A  69      26.961  39.000   1.988  1.00  7.48           C
ATOM      7  ND1 HIS A  69      26.646  37.666   2.008  1.00  8.20           N
ATOM      8  CD2 HIS A  69      26.183  39.603   1.064  1.00  9.47           C
ATOM      9  CE1 HIS A  69      25.670  37.467   1.090  1.00  9.09           C
ATOM     10  NE2 HIS A  69      25.385  38.635   0.523  1.00 10.94           N
TER
HETATM   11 ZN    ZN A 152      27.539  36.010   2.881  1.00  9.34          ZN
END
"""

def _i_seq(model, resname, name):
  for a in model.get_atoms():
    if a.parent().resname.strip() == resname and a.name.strip() == name:
      return a.i_seq
  raise RuntimeError("atom %s/%s not found" % (resname, name))

def _nd1_neighbor_names(model):
  """Names of the atoms bonded to His ND1 in the restraints, per get_all_bond_proxies."""
  atoms = model.get_atoms()
  bp, ap = model.get_restraints_manager().geometry.get_all_bond_proxies(
    sites_cart=model.get_sites_cart())
  bnl = Helpers.getBondedNeighborLists(atoms, bp, ap)
  for a in atoms:
    if a.parent().resname.strip() == "HIS" and a.name.strip() == "ND1":
      return set(n.name.strip() for n in bnl[a])
  return set()

def _optimizer_for(pdb_str, custom_bond=None):
  """Place H, make restraints, optionally add a custom bond, build the Optimizer.

  :param custom_bond: optional ((resname1, name1), (resname2, name2)) to bond.
  :return: (optimizer, model).  addFlipMovers=True so flip Movers are considered.
  """
  dm = DataManager(['model'])
  dm.process_model_str("tst_reduce_custom_bond.pdb", pdb_str)
  model = dm.get_model()
  model.set_log(null_out())
  add_h = reduce_hydrogen.place_hydrogens(model=model)
  add_h.run()
  model = add_h.get_model()
  model.set_log(null_out())
  p = mmtbx.model.manager.get_default_pdb_interpretation_params()
  p.pdb_interpretation.allow_polymer_cross_special_position = True
  p.pdb_interpretation.clash_guard.nonbonded_distance_threshold = None
  p.pdb_interpretation.proceed_with_excessive_length_bonds = True
  p.pdb_interpretation.disable_uc_volume_vs_n_atoms_check = True
  model.process(make_restraints=True, pdb_interpretation_params=p)
  if custom_bond is not None:
    (r1, n1), (r2, n2) = custom_bond
    proxy = geometry_restraints.bond_simple_proxy(
      i_seqs=(_i_seq(model, r1, n1), _i_seq(model, r2, n2)),
      distance_ideal=2.3, weight=1/0.02**2)
    model.get_restraints_manager().geometry.add_new_bond_restraints_in_place(
      proxies=[proxy], sites_cart=model.get_sites_cart())
  opt = Optimizers.Optimizer(Optimizers._philLike(), True, model,
                             bondedNeighborDepth=4)
  return opt, model

def test_his_is_flippable_without_custom_bond():
  """Control: the histidine gets a flip Mover when nothing unusual is bonded.

  Guards the fixture -- if this His ever stopped being flippable the bug test
  below would pass for the wrong reason.
  """
  opt, model = _optimizer_for(his_cys_pdb)
  assert _nd1_neighbor_names(model) == set(["CG", "CE1", "HD1"]), \
    _nd1_neighbor_names(model)
  assert "Added MoverHisFlip" in opt.getInfo(), opt.getInfo()
  print("test_his_is_flippable_without_custom_bond OK")

def test_custom_bond_on_his_ring_n_suppresses_flip():
  """Issue #1199: a user-defined bond on His ND1 must suppress the flip.

  MoverHisFlip only counts H/C neighbours, so without the _PlaceMovers check the
  SG bond is ignored and a free-to-flip His flip Mover is added, which can break
  the bond.
  """
  opt, model = _optimizer_for(his_cys_pdb, custom_bond=(("HIS", "ND1"), ("CYS", "SG")))
  # Sanity: the user-defined ND1-SG bond really is in the restraints the
  # optimizer sees (else the assertion below would pass for the wrong reason).
  assert "SG" in _nd1_neighbor_names(model), _nd1_neighbor_names(model)
  info = opt.getInfo()
  assert "Added MoverHisFlip" not in info, (
    "issue #1199: a His flip Mover was added for a residue whose ND1 carries a "
    "user-defined bond to Cys SG; flipping the ring would break that bond.\n"
    + info)
  # reduce2 protonates both ring nitrogens, so the Optimizer must delete the
  # hydrogen that conflicts with the bond (HD1 on the bonded ND1) rather than
  # leaving it dangling. This holds regardless of when/how the H was placed.
  deleted = set((h.parent().resname.strip(), h.name.strip())
                for h in opt.getHydrogensToDelete())
  assert ("HIS", "HD1") in deleted, deleted
  print("test_custom_bond_on_his_ring_n_suppresses_flip OK")

def test_metal_coordinated_his_stays_locked():
  """A His whose ND1 is genuinely bonded to a nearby Zn must stay locked.

  The ion lock-down handles this correctly (state 0, hydrogen removed).  The
  #1199 bonded-neighbour check must run AFTER the lock-down, so it must not turn
  this into a skipped residue -- otherwise the metal-coordinating hydrogen is
  left in place.
  """
  opt, model = _optimizer_for(his_zn_pdb, custom_bond=(("HIS", "ND1"), ("ZN", "ZN")))
  assert "ZN" in _nd1_neighbor_names(model), _nd1_neighbor_names(model)
  info = opt.getInfo()
  assert "Set MoverHisFlip" in info, (
    "the ion lock-down no longer handles a genuinely metal-coordinated His; the "
    "bonded-neighbour check must run AFTER the lock-down, not before it.\n" + info)
  assert "Added MoverHisFlip" not in info, info
  print("test_metal_coordinated_his_stays_locked OK")

# 4KW9 His 139, 147, 169, 202, 204 (Thr 203) and Zn 301-304: the ion
# lock-down removes one ring H from each His (5 H in 5 residues, no scoring).
his_zn4_pdb = """\
CRYST1   72.798   72.798  114.942  90.00  90.00  90.00 P 41 21 2     8
ATOM   2047  N   HIS A 139      12.238 -10.703  -4.125  1.00 45.45           N
ATOM   2048  CA  HIS A 139      13.628 -11.030  -4.461  1.00 43.59           C
ATOM   2049  C   HIS A 139      13.827 -11.172  -5.971  1.00 49.65           C
ATOM   2050  O   HIS A 139      13.993 -12.279  -6.486  1.00 53.96           O
ATOM   2051  CB  HIS A 139      14.042 -12.345  -3.792  1.00 46.80           C
ATOM   2052  CG  HIS A 139      14.837 -12.177  -2.536  1.00 46.84           C
ATOM   2053  ND1 HIS A 139      16.212 -12.046  -2.537  1.00 45.73           N
ATOM   2054  CD2 HIS A 139      14.455 -12.164  -1.236  1.00 45.15           C
ATOM   2055  CE1 HIS A 139      16.639 -11.937  -1.291  1.00 49.67           C
ATOM   2056  NE2 HIS A 139      15.593 -12.004  -0.483  1.00 47.69           N
ATOM   2176  N   HIS A 147      26.829   2.625 -10.023  1.00 15.30           N
ATOM   2177  CA  HIS A 147      27.480   3.904  -9.776  1.00 14.14           C
ATOM   2178  C   HIS A 147      26.421   4.969  -9.582  1.00 13.47           C
ATOM   2179  O   HIS A 147      25.523   5.118 -10.417  1.00 14.09           O
ATOM   2180  CB  HIS A 147      28.419   4.273 -10.938  1.00 14.62           C
ATOM   2181  CG  HIS A 147      27.774   4.175 -12.266  1.00 15.69           C
ATOM   2182  ND1 HIS A 147      27.251   5.269 -12.918  1.00 17.34           N
ATOM   2183  CD2 HIS A 147      27.531   3.106 -13.060  1.00 17.15           C
ATOM   2184  CE1 HIS A 147      26.715   4.876 -14.056  1.00 16.84           C
ATOM   2185  NE2 HIS A 147      26.881   3.570 -14.165  1.00 17.85           N
ATOM   2508  N   HIS A 169      18.101  -0.488  -9.449  1.00 35.45           N
ATOM   2509  CA  HIS A 169      18.572  -1.866  -9.572  1.00 35.70           C
ATOM   2510  C   HIS A 169      17.760  -2.552 -10.637  1.00 35.59           C
ATOM   2511  O   HIS A 169      16.546  -2.375 -10.711  1.00 35.81           O
ATOM   2512  CB  HIS A 169      18.477  -2.626  -8.236  1.00 33.64           C
ATOM   2513  CG  HIS A 169      19.486  -2.178  -7.226  1.00 33.72           C
ATOM   2514  ND1 HIS A 169      19.147  -1.456  -6.100  1.00 36.24           N
ATOM   2515  CD2 HIS A 169      20.833  -2.326  -7.187  1.00 36.16           C
ATOM   2516  CE1 HIS A 169      20.240  -1.185  -5.407  1.00 36.38           C
ATOM   2517  NE2 HIS A 169      21.278  -1.702  -6.042  1.00 34.35           N
ATOM   2951  N   HIS A 202      31.151  10.385  -5.909  1.00 12.17           N
ATOM   2952  CA  HIS A 202      31.040   9.189  -6.766  1.00 11.65           C
ATOM   2953  C   HIS A 202      31.046   7.970  -5.873  1.00 11.68           C
ATOM   2954  O   HIS A 202      31.973   7.779  -5.079  1.00 12.62           O
ATOM   2955  CB  HIS A 202      32.243   9.035  -7.699  1.00 13.00           C
ATOM   2956  CG  HIS A 202      32.300  10.055  -8.784  1.00 12.37           C
ATOM   2957  ND1 HIS A 202      31.710   9.866 -10.022  1.00 12.23           N
ATOM   2958  CD2 HIS A 202      32.900  11.270  -8.824  1.00 14.43           C
ATOM   2959  CE1 HIS A 202      31.910  10.953 -10.756  1.00 12.94           C
ATOM   2960  NE2 HIS A 202      32.638  11.808 -10.057  1.00 14.43           N
ATOM   2968  N   THR A 203      30.024   7.141  -5.990  1.00 11.93           N
ATOM   2969  CA  THR A 203      30.007   5.897  -5.233  1.00 11.18           C
ATOM   2970  C   THR A 203      29.849   4.756  -6.237  1.00 12.68           C
ATOM   2971  O   THR A 203      28.786   4.594  -6.845  1.00 12.31           O
ATOM   2972  CB  THR A 203      28.873   5.891  -4.205  1.00 17.21           C
ATOM   2973  OG1 THR A 203      28.972   7.055  -3.357  1.00 19.29           O
ATOM   2974  CG2 THR A 203      28.967   4.647  -3.328  1.00 16.95           C
ATOM   2982  N   HIS A 204      30.921   3.994  -6.412  1.00 12.55           N
ATOM   2983  CA  HIS A 204      30.984   2.935  -7.400  1.00 13.68           C
ATOM   2984  C   HIS A 204      31.145   1.619  -6.668  1.00 18.40           C
ATOM   2985  O   HIS A 204      32.121   1.428  -5.921  1.00 18.45           O
ATOM   2986  CB  HIS A 204      32.183   3.128  -8.320  1.00 14.32           C
ATOM   2987  CG  HIS A 204      32.186   4.437  -9.051  1.00 13.89           C
ATOM   2988  ND1 HIS A 204      33.288   4.902  -9.734  1.00 14.75           N
ATOM   2989  CD2 HIS A 204      31.227   5.380  -9.205  1.00 13.80           C
ATOM   2990  CE1 HIS A 204      33.009   6.071 -10.287  1.00 14.01           C
ATOM   2991  NE2 HIS A 204      31.759   6.377  -9.994  1.00 13.29           N
HETATM 3384 ZN    ZN A 301      30.927   8.097 -10.572  0.91 15.37          ZN
HETATM 3385 ZN    ZN A 302      25.968   2.617 -15.780  0.70 40.90          ZN
HETATM 3386 ZN    ZN A 303      17.515  -0.626  -5.210  0.64 32.42          ZN
HETATM 3387 ZN    ZN A 304      17.894 -11.737  -3.667  0.58 46.77          ZN
END
"""

def test_deleted_h_in_model_order():
  opt, model = _optimizer_for(his_zn4_pdb)
  deleted = opt.getHydrogensToDelete()
  names = [(h.parent().parent().resseq.strip(), h.name.strip()) for h in deleted]
  assert names == [("139", "HD1"), ("147", "HE2"), ("169", "HD1"), ("202", "HD1"),
    ("204", "HE2")], names
  i_seqs = [h.i_seq for h in deleted]
  assert i_seqs == sorted(i_seqs), i_seqs
  print("test_deleted_h_in_model_order OK")

if __name__ == "__main__":
  test_his_is_flippable_without_custom_bond()
  test_custom_bond_on_his_ring_n_suppresses_flip()
  test_metal_coordinated_his_stays_locked()
  test_deleted_h_in_model_order()
  print("OK")
