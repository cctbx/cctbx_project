from __future__ import absolute_import, division, print_function
import os
import copy
import shutil
from pathlib import Path

from iotbx.data_manager import DataManager
import libtbx
from libtbx import phil
from cctbx.maptbx.box import shift_and_box_model
from iotbx.map_model_manager import map_model_manager
from libtbx.utils import null_out
from cctbx.programs.qscore import Program as QscoreProgram

import numpy as np



from cctbx.maptbx.qscore import (
  generate_probes_np,
  shell_probes_precalculate,
  calc_qscore,

)

def isclose_or_nan(a, b, atol=1e-3):
  # perform isclose comparison, treating nans as equal
  return np.isclose(a, b, atol=atol) | (np.isnan(a) & np.isnan(b))


################################################################################
#### Test probe generation
################################################################################

def test_probe_generation():
  # test the primary points generation function against expected data
  sites_cart = np.array([[ 5.276, 12.488, 16.069],
                    [ 5.649, 13.947, 16.076]])

  probes_expected = np.array([[
          [ 5.276 , 12.488 , 15.969 ],
          [ 5.2588, 12.5558, 15.9976],
          [ 5.186 , 12.4803, 16.0261],
          [ 5.2564, 12.391 , 16.0547],
          [ 5.3636, 12.442 , 16.0833],
          [ 5.3304, 12.5601, 16.1119],
          [ 5.2115, 12.5151, 16.1404],
          [ 5.276 , 12.488 , 16.169 ]],

        [[ 5.649 , 13.947 , 15.976 ],
          [ 5.6318, 14.0148, 16.0046],
          [ 5.559 , 13.9393, 16.0331],
          [ 5.6294, 13.85  , 16.0617],
          [ 5.7366, 13.901 , 16.0903],
          [ 5.7034, 14.0191, 16.1189],
          [ 5.5845, 13.9741, 16.1474],
          [ 5.649 , 13.947 , 16.176 ]]])
  probes_xyz = generate_probes_np(sites_cart,0.1,8)
  assert np.all(np.isclose(probes_xyz,probes_expected,atol=1e-3))





def test_shell_probes():
  sites_cart = np.array([[ 5.276, 12.488, 16.069],
                [ 5.649, 13.947, 16.076]])


  # test precalculate (numpy)

  expected_probes = np.array(

                      [[[5.2760, 12.4880, 14.5690],
                        [5.0515, 13.4037, 14.9023],
                        [4.0297, 12.4397, 15.2357],
                        [4.8250, 11.1476, 15.5690],
                        [6.3669, 11.4721, 15.9023],
                        [6.6515, 13.0626, 16.2357],
                        [5.2302, 13.9015, 16.5690],
                        [4.0466, 12.6981, 16.9023],
                        [5.3430, 11.5476, 17.2357],
                        [5.2760, 12.4880, 17.5690]],

                        [[5.6490, 13.9470, 14.5760],
                        [5.4245, 14.8627, 14.9093],
                        [4.4027, 13.8987, 15.2427],
                        [5.1980, 12.6066, 15.5760],
                        [6.7399, 12.9311, 15.9093],
                        [7.0245, 14.5216, 16.2427],
                        [5.6032, 15.3605, 16.5760],
                        [4.4196, 14.1571, 16.9093],
                        [5.7160, 13.0066, 17.2427],
                        [5.6490, 13.9470, 17.5760]]])

  shell_func = shell_probes_precalculate
  probe_xyz, probe_mask = shell_func(
                      sites_cart=sites_cart,
                      atoms_tree = None,
                      selection_bool=None,
                      n_probes=10,
                      RAD=1.5,
                      rtol=0.9,
                      log = null_out())

  assert np.all(isclose_or_nan(probe_xyz,expected_probes,atol=1e-3))


  # test at single shell
  probe_xyz = np.array(probe_xyz)
  probe_xyz = probe_xyz.reshape(expected_probes.shape)
  assert np.all(isclose_or_nan(probe_xyz,expected_probes,atol=1e-3))



################################################################################
#### Test templates for real data
################################################################################

def convert_dict_to_group_args(d):
  # Recursive conversion of nested dictionary to nested group_args
  def convert_func(d):
    if isinstance(d, dict):
      converted_items = {k: convert_func(v) for k, v in d.items()}
      return libtbx.group_args(**converted_items)
    else:
      return d

  return convert_func(d)

def convert_group_args_to_dict(g):
  # Recursive conversion of nested group args to nested dictionary
  def convert_func(g):
    if isinstance(g, libtbx.group_args):
      converted_items = {k: convert_func(v) for k, v in g.__dict__.items() if not k.startswith('_')}
      return converted_items
    else:
      return g

  return convert_func(g)


# a template to store configuration for a single test
test_template ={

         "data":{
          "name":None,
          "model_str":None,
          "model_file":None,
          "map_file":None,
          "test_dir":None,
          },
        "results":{
           "expected":{},
           "calc":{},
        },
         "params":{
            "selection":None, # Just calculate q score for a sub-selection
            "shells": [0.0,0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.,
                       1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2. ],
            "n_probes":8,
            "nproc":4,
            "debug":True,
            "rtol":0.9,
         }
}
def run_test(test):

  if isinstance(test,dict):
    test = convert_dict_to_group_args(test)

  print()
  print("#"*79)
  print("Running test "+test.data.name)
  print("\tModel File:",test.data.model_file)
  print("\tMap:",test.data.map_file)
  if test.data.model_str is not None:
    print("\tModel Str:\n",test.data.map_file)

  dm = DataManager()
  assert [test.data.model_file,test.data.model_str].count(None)==1
  if test.data.model_file is not None:
    dm.process_model_file(test.data.model_file)
  elif test.data.model_str is not None:
    dm.process_model_str(test.data.name,test.data.model_str)

  # get model
  model = dm.get_model()

  # rebox as p1
  #model = shift_and_box_model(model)


  # get data from file or calculated
  if test.data.map_file is not None:
    dm.process_real_map_file(test.data.map_file)
    mm = dm.get_real_map()
    mmm = map_model_manager(model=model,map_manager=mm)
  else:
    mmm = map_model_manager(model=model)
    mmm.generate_map(d_min=2)

  params = test.params
  q_func = calc_qscore

  params = convert_group_args_to_dict(params)
  result = q_func(mmm,
                  **params,
                  log=null_out())
  test.results.calc = convert_dict_to_group_args(result)
  for key,value_expected in test.results.expected.__dict__.items():
    value_calc= np.array(test.results.calc.__dict__[key])
    assert np.all(isclose_or_nan(value_expected,value_calc)), (
f"Failed test: {test.data.name}\nExpected:\n{value_expected}\nCalc:\n{value_calc}")
  print("Done.\n")
  return convert_group_args_to_dict(test)

################################################################################
#### Define and run tests
################################################################################

def build_tests(test_dir="qscore_tst_dir"):

  # make test dir
  test_dir = Path("qscore_tst_dir")
  if test_dir.exists():
    shutil.rmtree(test_dir)
  test_dir.mkdir()


  tests = {}


  #4. map/model from regression dir (tst2) numpy/precalculate
  test = copy.deepcopy(test_template)
  tst2_model_file = libtbx.env.find_in_repositories(
  relative_path=\
  f"cctbx_project/iotbx/regression/data/non_zero_origin_model_split.pdb",
  test=os.path.isfile)
  tst2_map_file = libtbx.env.find_in_repositories(
  relative_path=\
  f"cctbx_project/iotbx/regression/data/non_zero_origin_map.ccp4",
  test=os.path.isfile)
  test["data"]["model_file"] = tst2_model_file
  test["data"]["map_file"] = tst2_map_file
  test["data"]["name"] = "tst2_precalc_numpy"
  test["data"]["test_dir"] = test_dir
  test["params"]["n_probes"] = 32
  expected_qscore_per_atom= np.array([
    0.79852,  0.80609,  0.80133,  0.72103,  0.75883,  0.81456,  0.82049,
    0.77932,  0.77675,  0.78246,  0.84899,  0.71687,  0.77178,  0.82013,
    0.82884,  0.82516,  0.88237,  0.79658,  0.84557,  0.80930,  0.78813,
    0.80125,  0.66020,  0.77447,  0.82493,  0.80411,  0.81775,  0.75656,
    0.79463,  0.84372,  0.66545,  0.78583,  0.80262,  0.88354,  0.80967,
    0.83346,  0.70284,  0.83037,  0.79684,  0.69537,  0.74061,  0.79931,
    0.88578,  0.86837,  0.79809,  0.70589,  0.76989,  0.82633,  0.60241,
    0.77312,  0.67250,  0.79516,  0.78177,  0.72057,  0.77841,  0.84436,
    0.81322,  0.77128,  0.79111,  0.83423,  0.82586,  0.81274,  0.70658,
    0.76125,  0.83175,  0.76124,  0.78605,  0.83369,  0.76933,  0.81073,
    0.67650,  0.73330,  0.83762,  0.76627,  0.77703,  0.77652,  0.78147,
    0.85242,  0.70684,  0.81237,  0.80886,  0.87575,  0.74811,  0.74911,
    0.87672,  0.76362,
  ])
  test["results"]["expected"]["qscore_per_atom"] = expected_qscore_per_atom
  tests[test["data"]["name"]] = test



  #5. 1yjp with simulated density
  # precalculated with numpy
  test = copy.deepcopy(test_template)
  yjp_model_file = libtbx.env.find_in_repositories(
    relative_path=\
    f"cctbx_project/iotbx/regression/data/1yjp.pdb",
    test=os.path.isfile)
  test["data"]["model_file"] = yjp_model_file
  test["data"]["name"] = "1yjp_precalc_numpy"
  test["data"]["test_dir"] = test_dir
  test["data"]["map_file"] = None
  test["params"]["n_probes"] = 32
  expected_qscore_per_atom= np.array([
    0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,
    0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,
    -0.10505,  0.03062,  0.00000,  0.67700,  0.95088,  0.94763,  0.75051,
    0.89511,  0.87489,  0.97221,  0.91645,  0.91834,  0.13628,  0.22108,
    0.09771,  0.94044,  0.87528,  0.92334,  0.95104,  0.92881,  0.91000,
    0.88136,  0.97055,  0.95707,  0.93397,  0.86285,  0.87010,  0.97105,
    0.90625,  0.88578,  0.96754,  0.96147,  0.93298,  0.86254,  0.89024,
    0.95584,  0.90034,  0.91163,  0.62957,  0.91963,  0.23041,  0.91662,
    0.89745,  0.95870,  0.94299,  0.00000,  0.97742,  0.00000,  0.95322,
    0.96996,  0.89275,  0.91926,
  ])
  test["results"]["expected"]["qscore_per_atom"] = expected_qscore_per_atom
  tests[test["data"]["name"]] = test


  #6. 1yjp with simulated density and shift and box
  dm = DataManager()
  dm.process_model_file(yjp_model_file)
  model = dm.get_model()
  model = shift_and_box_model(model)
  boxed_1yjp_path = Path(test_dir,"1yjp_boxed.pdb")
  dm.write_model_file(model.model_as_pdb(),filename=str(boxed_1yjp_path),overwrite=True)

  # precalculated with numpy
  test = copy.deepcopy(test_template)
  test["data"]["name"] = "1yjp_boxed_precalc_numpy"
  test["data"]["test_dir"] = test_dir
  #test["data"]["model_file"] = str(boxed_1yjp_path)
  test["data"]["model_str"] = model.model_as_pdb()
  test["data"]["map_file"] = None
  test["params"]["n_probes"] = 32
  expected_qscore_per_atom = np.array([
    0.96035,  0.91590,  0.92283,  0.96654,  0.94164,  0.89421,  0.91464,
    0.97026,  0.91649,  0.91408,  0.97073,  0.95491,  0.94386,  0.88791,
    0.90950,  0.97556,  0.92062,  0.92067,  0.96957,  0.95537,  0.94550,
    0.87780,  0.89461,  0.97148,  0.93788,  0.92123,  0.90062,  0.96383,
    0.96746,  0.94327,  0.88461,  0.90352,  0.97037,  0.93105,  0.92621,
    0.90266,  0.96954,  0.96153,  0.94760,  0.86363,  0.88383,  0.96882,
    0.92100,  0.91835,  0.96474,  0.95877,  0.92953,  0.88268,  0.89251,
    0.95668,  0.91766,  0.91907,  0.92488,  0.92805,  0.93958,  0.93576,
    0.90368,  0.96675,  0.95459,  0.97762,  0.97334,  0.98168,  0.97580,
    0.97306,  0.96717,  0.96175,
  ])
  test["results"]["expected"]["qscore_per_atom"] = expected_qscore_per_atom
  tests[test["data"]["name"]] = test


  #7 tst2 with shift_and_box (numpy,precalculate)

  tst2_model_file = libtbx.env.find_in_repositories(
  relative_path=\
  f"cctbx_project/iotbx/regression/data/non_zero_origin_model_split.pdb",
  test=os.path.isfile)
  tst2_map_file = libtbx.env.find_in_repositories(
  relative_path=\
  f"cctbx_project/iotbx/regression/data/non_zero_origin_map.ccp4",
  test=os.path.isfile)

  dm = DataManager()
  dm.process_model_file(tst2_model_file)
  model = dm.get_model()
  #model = shift_and_box_model(model)
  boxed_tst2_path = Path(test_dir,"tst2_boxed.pdb")
  dm.write_model_file(model.model_as_pdb(),filename=str(boxed_tst2_path),overwrite=True)

  test = copy.deepcopy(test_template)
  test["data"]["model_file"] = str(boxed_tst2_path)
  test["data"]["map_file"] = tst2_map_file
  test["data"]["name"] = "tst2_boxed_precalc_numpy"
  test["data"]["test_dir"] = test_dir
  test["params"]["n_probes"] = 32
  expected_qscore_per_atom = np.array([
  0.79852,  0.80609,  0.80133,  0.72103,  0.75883,  0.81456,  0.82049,
  0.77932,  0.77675,  0.78246,  0.84899,  0.71687,  0.77178,  0.82013,
  0.82884,  0.82516,  0.88237,  0.79658,  0.84557,  0.80930,  0.78813,
  0.80125,  0.66020,  0.77447,  0.82493,  0.80411,  0.81775,  0.75656,
  0.79463,  0.84372,  0.66545,  0.78583,  0.80262,  0.88354,  0.80967,
  0.83346,  0.70284,  0.83037,  0.79684,  0.69537,  0.74061,  0.79931,
  0.88578,  0.86837,  0.79809,  0.70589,  0.76989,  0.82633,  0.60241,
  0.77312,  0.67250,  0.79516,  0.78177,  0.72057,  0.77841,  0.84436,
  0.81322,  0.77128,  0.79111,  0.83423,  0.82586,  0.81274,  0.70658,
  0.76125,  0.83175,  0.76124,  0.78605,  0.83369,  0.76933,  0.81073,
  0.67650,  0.73330,  0.83762,  0.76627,  0.77703,  0.77652,  0.78147,
  0.85242,  0.70684,  0.81237,  0.80886,  0.87575,  0.74811,  0.74911,
  0.87672,  0.76362,
  ])
  test["results"]["expected"]["qscore_per_atom"] = expected_qscore_per_atom
  tests[test["data"]["name"]] = test



  return tests


def test_program_template(test):
  print("Running prgram template test: ",test["data"]["name"])
  model_file = test["data"]["model_file"]
  map_file = test["data"]["map_file"]
  dm = DataManager()
  dm.process_model_file(str(model_file))
  dm.process_real_map_file(str(map_file))
  params = phil.parse(QscoreProgram.master_phil_str,process_includes=True).extract()
  task = QscoreProgram(dm,params)
  task.run()
  results = task.get_results()


if (__name__ == "__main__"):


  # test probe generation
  test_probe_generation()

  # test single shell probe generation
  test_shell_probes()



  # Test on some real models
  tests = build_tests()
  for test_name,test in tests.items():
    test = run_test(test)


  # Test if program template runs, with default phil
  test_program_template(list(tests.values())[-1])
