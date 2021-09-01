from __future__ import absolute_import, division, print_function
import mmtbx.model
import libtbx.load_env
from libtbx.utils import format_cpu_times
from libtbx.test_utils import approx_equal
import iotbx.pdb
from mmtbx.regression.model import tst_model_cart_ref_restraints


def exercise_adopting_ref_tors_restraints_h():
  params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  params.pdb_interpretation.flip_symmetric_amino_acids=False
  inp_1 = iotbx.pdb.input(lines=tst_model_cart_ref_restraints.pdb_str_h, source_info=None)
  h_model = mmtbx.model.manager(model_input = inp_1)
  h_model.process(pdb_interpretation_params=params, make_restraints=True)

  inp_2 = iotbx.pdb.input(lines=tst_model_cart_ref_restraints.pdb_str, source_info=None)
  model = mmtbx.model.manager(model_input = inp_2)
  model.process(pdb_interpretation_params=params, make_restraints=True)
  params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  params.pdb_interpretation.flip_symmetric_amino_acids=False
  params.reference_model.enabled=True
  params.reference_model.sigma = 2

  # case 1: big model adopting small model.
  h_model.set_reference_torsion_restraints(
      ref_model = model,
      params=params)
  assert h_model.get_restraints_manager().geometry.get_n_reference_dihedral_proxies() == 34
  reference_dihedral_proxies = h_model.get_restraints_manager().geometry.reference_dihedral_manager.reference_dihedral_proxies
  answer = [
[(1,   2,  9, 10), 171.515643141],
[(9,  10, 13, 14), -57.7746022999],
[(10, 13, 14, 15), -60.8369797925],
[(10, 11, 23, 24), -176.236177141],
[(23, 24, 27, 28), -176.840947675],
[(24, 27, 28, 29), 61.8839679236],
[(24, 25, 37, 38), 166.212120415],
[(37, 38, 41, 42), -164.6572256],
[(38, 41, 42, 43), -160.727030714],
[(41, 42, 43, 44), 54.0842898531],
[(38, 39, 54, 55), 173.187552548],
[(54, 55, 58, 59), -62.8213683054],
[(55, 58, 59, 60), 159.926433285],
[(58, 59, 60, 61), -72.3924878566],
[(55, 56, 71, 72), -176.726828257],
[(71, 72, 75, 76), -54.9041387187],
[(72, 75, 76, 77), -66.5778205786],
[(72, 73, 85, 86), -171.777645364],
[(85, 86, 89, 90), -60.9936732039],
[(86, 89, 90, 91), 91.6226963647],
[(91, 93, 95, 96), -179.330180577],
[(85, 86, 87, 97), 157.562627012],
[(0,   1,  2,  9), -162.221584042],
[(2,   9, 10, 11), -60.5820017501],
[(9,  10, 11, 23), 141.187189652],
[(11, 23, 24, 25), -119.250599331],
[(23, 24, 25, 37), 125.146913149],
[(25, 37, 38, 39), -126.159284697],
[(37, 38, 39, 54), 112.807868032],
[(39, 54, 55, 56), -114.975951407],
[(54, 55, 56, 71), 126.756425626],
[(56, 71, 72, 73), -116.419938986],
[(71, 72, 73, 85), 97.6949671859],
[(73, 85, 86, 87), -80.0636781782],
  ]

  for i, dp in enumerate(reference_dihedral_proxies):
    # print(dir(dp))
    # print(dp.i_seqs, dp.angle_ideal, answer[i][1])
    assert dp.i_seqs == answer[i][0]
    assert approx_equal(dp.angle_ideal, answer[i][1])
    assert approx_equal(dp.weight, 0.25)

  answer = [
[(1,   2,  4,  5), 171.515643141],
[(4,   5,  8,  9), -57.7746022999],
[(5,   8,  9, 10), -60.8369797925],
[(5,   6, 12, 13), -176.236177141],
[(12, 13, 16, 17), -176.840947675],
[(13, 16, 17, 18), 61.8839679236],
[(13, 14, 20, 21), 166.212120415],
[(20, 21, 24, 25), -164.6572256],
[(21, 24, 25, 26), -160.727030714],
[(24, 25, 26, 27), 54.0842898531],
[(21, 22, 29, 30), 173.187552548],
[(29, 30, 33, 34), -62.8213683054],
[(30, 33, 34, 35), 159.926433285],
[(33, 34, 35, 36), -72.3924878566],
[(30, 31, 38, 39), -176.726828257],
[(38, 39, 42, 43), -54.9041387187],
[(39, 42, 43, 44), -66.5778205786],
[(39, 40, 46, 47), -171.777645364],
[(46, 47, 50, 51), -60.9936732039],
[(47, 50, 51, 52), 91.6226963647],
[(52, 54, 56, 57), -179.330180577],
[(46, 47, 48, 58), 157.562627012],
[(0,   1,  2,  4), -162.221584042],
[(2,   4,  5,  6), -60.5820017501],
[(4,   5,  6, 12), 141.187189652],
[(6,  12, 13, 14), -119.250599331],
[(12, 13, 14, 20), 125.146913149],
[(14, 20, 21, 22), -126.159284697],
[(20, 21, 22, 29), 112.807868032],
[(22, 29, 30, 31), -114.975951407],
[(29, 30, 31, 38), 126.756425626],
[(31, 38, 39, 40), -116.419938986],
[(38, 39, 40, 46), 97.6949671859],
[(40, 46, 47, 48), -80.0636781782],
  ]
  params.reference_model.sigma = 4
  model.set_reference_torsion_restraints(
      ref_model = h_model,
      params=params)
  assert model.get_restraints_manager().geometry.get_n_reference_dihedral_proxies() == 34
  reference_dihedral_proxies = model.get_restraints_manager().geometry.reference_dihedral_manager.reference_dihedral_proxies
  for i, dp in enumerate(reference_dihedral_proxies):
    # print dp.i_seqs, dp.angle_ideal
    assert dp.i_seqs == answer[i][0]
    assert approx_equal(dp.angle_ideal, answer[i][1])
    assert approx_equal(dp.weight, 0.0625)

  # changing coords for every other atoms to modify torsion angles
  f = True
  for a in h_model.get_hierarchy().atoms():
    if f:
      a.xyz = (a.xyz[0]+0.5, a.xyz[1], a.xyz[2])
      f = False
    else:
      f = True
  new_targets = [
173.233628557,
-65.911305977,
-14.0162554531,
162.629939942,
-177.848862887,
84.054568441,
150.088700474,
-169.113940647,
-162.410889324,
98.509473131,
176.247405316,
-77.6303938113,
161.640293376,
-109.88499213,
-179.297617725,
-73.9436357254,
-5.54761865294,
154.179231454,
-58.3916258056,
89.6379144637,
-167.538405775,
155.271991752,
-168.005387354,
-38.998355844,
140.733342841,
-116.817012423,
131.119009065,
-122.898932369,
133.945203506,
-136.386148677,
113.946903271,
-102.742029233,
104.376893428,
-65.9529100865,
  ]
  h_model.set_sites_cart_from_hierarchy()
  model.set_reference_torsion_restraints(
      ref_model = h_model,
      params=params)
  assert model.get_restraints_manager().geometry.get_n_reference_dihedral_proxies() == 34
  reference_dihedral_proxies = model.get_restraints_manager().geometry.reference_dihedral_manager.reference_dihedral_proxies
  for i, dp in enumerate(reference_dihedral_proxies):
    # print dp.i_seqs, dp.angle_ideal
    assert dp.i_seqs == answer[i][0]
    assert approx_equal(dp.angle_ideal, new_targets[i])
    assert approx_equal(dp.weight, 0.0625)



def run():
  if (not libtbx.env.has_module("reduce")):
    print("Reduce not installed")
    return
  exercise_adopting_ref_tors_restraints_h()

  print(format_cpu_times())

if (__name__ == "__main__"):
  run()
