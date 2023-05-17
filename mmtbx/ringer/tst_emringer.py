
from __future__ import absolute_import, division, print_function
from mmtbx.ringer import em_scoring as score
from mmtbx.programs import emringer
import iotbx.pdb
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal, Exception_expected
from libtbx.utils import null_out, Sorry
import libtbx.load_env
import warnings
import os.path
from iotbx.cli_parser import run_program
from six.moves import range

def exercise_emringer_residue_scan():
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/mmtbx/em_ringer/tst_emringer_model.pdb",
    test=os.path.isfile)
  map_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/mmtbx/em_ringer/tst_emringer_map.ccp4",
    test=os.path.isfile)
  assert (not None in [pdb_file, map_file])
  emringer_results = run_program(program_class=emringer.Program, args=[pdb_file, map_file, 'quiet=True'])
  results = emringer_results.ringer_result
  scoring = emringer_results.scoring_result
  rolling = emringer_results.rolling_result
  #results, scoring, rolling = emringer.run([pdb_file, map_file], out=null_out())
  # Make sure the right number of residues (22 out of 28) get scanned
  assert len(results)==22
  modelled_list = [290.742121792,192.844056257,45.4781110306,294.247825632,303.618891108,58.7694040824,331.70068496,46.7136045049,290.167261226,304.261231829,282.651244586,268.729721112,195.972333785,305.321933311,314.81066224,286.028424514,311.180807466,313.004918133,296.67781565,296.949191638,169.644245088,192.496265164]
  peak_list = [270,180,260,75,305,30,310,90,265,270,270,240,280,260,310,285,295,100,260,165,155,200]
  peak_rhos = [0.175600306502,0.351591946536,0.206238983746,0.3269057296,0.68375562882,0.251143527693,0.29106077218,0.199922124642,0.298461589197,0.563313760047,0.412696803251,0.511080434089,0.310001828446,0.228239176285,0.563148497472,0.490755919184,0.200978032127,0.274929619102,0.299229846335,0.179215798655,0.150783734124,0.210869945593]
  for i in range(22):
    # Make sure the modelled angle is correctly read
    assert approx_equal(results[i]._angles[1].angle_current, modelled_list[i])
    # Make sure the peak is chosen correctly
    assert approx_equal(results[i]._angles[1].peak_chi, peak_list[i])
    # Make sure the peak rhos are correct
    assert approx_equal(results[i]._angles[1].peak_rho, peak_rhos[i])

  emringer_results2 = run_program(program_class=emringer.Program, args=[pdb_file, map_file, "rolling_window_threshold=0.5", 'quiet=True'])
  #results, scoring2, rolling2 = emringer.run([pdb_file, map_file, "rolling_window_threshold=0.5"], out=null_out())
  results2 = emringer_results2.ringer_result
  scoring2 = emringer_results2.scoring_result
  rolling2 = emringer_results2.rolling_result
  assert rolling.threshold == 0
  assert rolling2.threshold == 0.5
  #print rolling.results_a[0]
  #print rolling2.results_a[0]
  # just making sure this doesn't break!
  #results, scoring2, rolling = emringer.run([pdb_file, map_file, "sampling_angle=2"], out=null_out())


def exercise_emringer_insertion_codes():
  """
  Checks that emringer doesn't crash when there are insertion codes.
  The correctness of output is not checked.
  """
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/mmtbx/em_ringer/tst_emringer_insertion_codes_model.pdb",
    test=os.path.isfile)
  map_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/mmtbx/em_ringer/tst_emringer_map.ccp4",
    test=os.path.isfile)
  assert (not None in [pdb_file, map_file])
  emringer_results = run_program(program_class=emringer.Program, args=[pdb_file, map_file, 'quiet=True'])

#  results, scoring, rolling = emringer.run([pdb_file, map_file], out=null_out())

# FIXME this will fail right now, which is deliberate
def exercise_emringer_out_of_bounds():
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/mmtbx/em_ringer/tst_emringer_model.pdb",
    test=os.path.isfile)
  map_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/mmtbx/em_ringer/tst_emringer_map.ccp4",
    test=os.path.isfile)
  assert (not None in [pdb_file, map_file])
  pdb_in = iotbx.pdb.input(pdb_file)
  hierarchy = pdb_in.construct_hierarchy()
  xyz = hierarchy.atoms().extract_xyz()
  xyz += flex.vec3_double(xyz.size(), (200., 0.0, 0.0))
  hierarchy.atoms().set_xyz(xyz)
  with open("tst_emringer_shifted.pdb", "w") as f:
    f.write(hierarchy.as_pdb_string(
      crystal_symmetry=pdb_in.crystal_symmetry()))
  args = ["tst_emringer_shifted.pdb", map_file]
  try:
    results, s, r = emringer.run(args, out=null_out())
  except Sorry as e:
    pass
  else:
    raise Exception_expected

def exercise_emringer_pickle_loading():
  pkl_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/mmtbx/em_ringer/tst_emringer_pickle.pkl",
    test=os.path.isfile)
  waves,thresholds = score.parse_pickle(pkl_file)
  assert approx_equal(thresholds,[0.15145128496323776, 0.17806650215606262, 0.20468171934888746, 0.23129693654171235, 0.25791215373453719, 0.28452737092736202, 0.31114258812018691, 0.33775780531301181, 0.36437302250583659, 0.39098823969866148, 0.41760345689148637, 0.44421867408431126, 0.47083389127713604, 0.49744910846996093, 0.52406432566278582, 0.5506795428556106, 0.5772947600484355, 0.60390997724126039, 0.63052519443408517, 0.65714041162691006])
  return waves, thresholds

def exercise_emringer_peakfinding(waves):
  list = score.Peaklist()
  for i in waves:
    list.append_lists(score.calculate_peaks(i,0.4))
  assert len(list) == 6
  print(list)
  assert [i.chi_value*5 for i in list.get_peaks()] == [305, 270, 270, 240, 310, 285]


def exercise_emringer_statistics():
  # Test statistic calculation
  peak_list = [5]+[0]*5+[50]+[0]*12+[25]+[0]*10+[10]+[0]*12+[5]+[0]*10+[10]+[0]*17
  binned_peaks = score.calculate_binned_counts(peak_list)
  assert binned_peaks==[50,25,10,5,10,5]

  zscore, rotamer_ratio = score.statistic(binned_peaks)
  assert approx_equal(rotamer_ratio, 2.0/3)
  assert approx_equal(zscore, 2.570679)

  new_rotamer_ratio, new_zscore  = score.calc_ratio(peak_list)
  print(zscore)
  print(new_zscore)
  assert approx_equal(new_rotamer_ratio, rotamer_ratio)
  assert approx_equal(new_zscore,zscore)

if __name__=='__main__':
  keep_going=True
  try:
    import wx # special import
  except ImportError:
    print("Required cctbx irrelevant dependencies are missing, skipping test.")
    keep_going=False
  tstdir = libtbx.env.find_in_repositories("phenix_regression/mmtbx/em_ringer")
  if (tstdir is None):
    warnings.warn("phenix_regression not available, skipping test")
  else :
    if(keep_going):
      exercise_emringer_residue_scan()
      exercise_emringer_insertion_codes()
      # FIXME
      #exercise_emringer_out_of_bounds()
      #w, t = exercise_emringer_pickle_loading()
      #exercise_emringer_peakfinding(w)
      print("OK")
