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

if __name__ == "__main__":
  test_his_is_flippable_without_custom_bond()
  test_custom_bond_on_his_ring_n_suppresses_flip()
  test_metal_coordinated_his_stays_locked()
  print("OK")
