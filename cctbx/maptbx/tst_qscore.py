from __future__ import absolute_import, division, print_function
import os
<<<<<<< Updated upstream
from cctbx.array_family import flex
from iotbx.cli_parser import run_program
import libtbx
from libtbx.utils import null_out
from libtbx.test_utils import approx_equal
from cctbx.programs import qscore


expected_results = {
    17: [
        0.35535516324938227,
        0.24550687752861167,
        0.42717508108780494,
        0.2630431196992697,
        0.2512972739404068,
        0.0748338849699907,
        0.3951188616344501,
        0.017261865581353543,
        0.47942130214316503,
        0.2927000401754054,
        0.08603369430277283,
        0.30219565421918,
        0.39417849168365854,
        0.6189297826154754,
        0.36598113428758416,
        0.07821030825399476,
        0.5275589560551526,
        0.19580632885270882,
        0.27077417168391377,
        0.31471829887209857,
        0.3158786931123718,
        0.5096481626657788,
        0.4609767956764144,
        0.28235509392736674,
        0.4709590709747482,
        0.2939374853294539,
        0.10241906482206845,
        0.4015066315039059,
        0.35209047260985793,
        0.35059709887741874,
        0.3265106583648,
        0.2806311093313803,
        0.37458818716858017,
        0.385126780724512,
        0.06154713873358951,
        0.10991564773911004,
        0.49364138455921214,
        0.2544807372533646,
        0.5054205573973468,
        0.2637558217581824,
        0.3140624661999957,
        0.2365647784766615,
        0.28262674351648415,
        0.1969617894777505,
        0.3197038813888359,
        0.004434119579554559,
        0.49852526974799255,
        0.4284389897113832,
        0.32266761533489585,
        0.5064439221443267,
        0.4922792436039985,
        0.6303751667061026,
        0.1542182251535256,
        0.17137072347023904,
        0.5974219741230471,
        0.3533762354980602,
        0.27000081206540955,
        0.2007048933666327,
        0.30627931992625224,
        0.4333678839790156,
        0.5452464031846214,
        0.39074700178563093,
        0.44160485208870487,
        0.3012758760420372,
        0.3250739032545345,
        0.43362607699544053,
        0.3625137062692818,
        0.3737674696637951,
        0.42857066342176636,
        0.35808570562399467,
        0.22410204647321708,
        0.3767186160303378,
        0.28651652046017084,
        -0.049724839252815733,
        0.19356789252833098,
        0.09094302118136462,
        0.42016933552912894,
        0.13060556101358295,
        0.22990514778100643,
        0.24778269855730797,
        0.3321627503379433,
        0.2793544145496639,
        0.36436645199153783,
        0.1990021248328906,
        0.5972764932446221,
        0.38459252021474777,
        0.15376046116368622,
        0.4350468866909971,
        0.34446835591022873,
        0.36878927260006605,
        0.05794444507950947,
        0.1368028652112549,
        0.47365762715495435,
        0.40950604509893035,
        0.21697233895086965,
        0.24302703044362328,
        0.3419794880192022,
        0.5585194498374837,
        0.459308252423094,
        0.3230117486740559,
        0.3148210060317072,
        0.13707459657834334,
        0.35686464668359086,
        0.3070311263154991,
        0.3055250982319954,
        0.5000373760197506,
        0.2828409443923543,
        0.22088517997784513,
        0.34473673110006714,
        0.39401167930098646,
        0.1409604644258536,
        0.17904487505641237,
        0.41512532606270575,
        0.051477839079248605,
        0.430233223946139,
        0.10701871986502859,
        0.3554414516983775,
        0.15227178169975683,
        0.20902369264131399,
        0.27523031292228356,
        0.5014761466259059,
        0.07273646238705293,
        0.5228074093643269,
        0.28213307002470817,
        0.32834396845308833,
        0.4546359850826332,
        0.2640898256601553,
        0.4042889388009938,
        -0.058819305685534916,
        0.3320355211739059,
        0.5368734032575433,
        0.38723919567921916,
        0.27920683488494874,
        0.21293610438730504,
        0.28101278857543577,
        0.5005434542705948,
        0.3471817144906104,
        0.21127407819640626,
        0.32067511222833683,
        0.11079335636572182,
        0.17036780464871323,
        0.41955675987187646,
        0.3905657146499568,
        0.3401149986039323,
        0.36035223718270976,
        0.23917228790998196,
        0.36381317879482006,
        0.39077694083716824,
        0.03118681161607613,
        0.04641661292929102,
        0.3724988531432701,
        0.2643036127826901,
        0.30537188074803256,
        0.16202886247768936,
        0.43568155630265554,
        0.3084634599307411,
        0.387489124753754,
        0.2666390040428182,
        0.46279277413454706,
        0.02002279771517451,
        0.4827489588800846,
        0.445680062618974,
        0.14784923408302839,
        0.4184692111286882,
        0.3933408944964084,
        0.4550840595198136,
        -0.06497059743774909,
        0.20326928453476295,
        0.5410683745517895,
        0.3934791271999539,
        0.24834218260971685,
        0.2284335972515956,
        0.2749252359472737,
        0.5811000263463834,
        0.5886279113579868,
        0.10014876331174488,
        0.3558056816456358,
        0.23931433448674996,
        0.28193700267688226,
        0.30923904707665817,
        0.3937787470567507,
        0.16583160065281066,
        0.41937141588607457,
        0.11413599799171258,
        0.3671212694643467,
        0.31069640919716385,
        0.05525220154724908,
        0.247288738983681,
        0.39689271509298957,
        0.20134818886327846,
        0.5060624692935539,
        0.23340144114502143,
        0.1047045872234278,
        0.26762911598483546,
        0.39263746442670566,
        0.12692533917487203,
        0.3980087754375635,
        0.03756336993253084,
        0.46101152409796214,
        0.42446162573971874,
        0.3643343019533042,
        0.5130778316649313,
        0.3063184055501381,
        0.47635554117162227,
        0.2022659217974044,
        0.2554080063069391,
        0.5501695893864158,
        0.2846628653680698,
        0.3242176784889862,
        0.24811871246412262,
        0.3255252546682871,
        0.5453376573614725,
        0.41301426891418463,
        0.07307555071072433,
        0.4536301022017574,
        0.24638397492302436,
        0.19814077792624293,
        0.3602238859424958,
        0.3170086525644145,
        0.2972960450133945,
        0.4145260546137287,
        0.2661909058919775,
        0.49307051785060885,
        0.42049092767916485,
        0.010562206055517411,
        0.05167120683170355,
        0.48650651743605644,
        0.07522293396734439,
        0.4990564081521944,
        0.3334757193364281,
        0.14185320419595457,
        0.12828315708419813,
        0.11089188858896354,
        0.14452747870340907,
        0.3699940929152232,
        -0.041235846823679784,
        0.513812878496174,
        0.3694188603001491,
        0.28901185662505796,
        0.3803523709807769,
        0.3765605158390536,
        0.4527272022777316,
        -0.08048732044100083,
        0.18539904739806815,
        0.5232189411542995,
        0.3022522764188308,
        0.35572393347262615,
        0.27116261764046895,
        0.3626082715676051,
        0.5619703843998924,
        0.4944482273684381,
        0.16248076555523522,
        0.3206340194265775,
        0.14751443872046044,
        0.09595366681825422,
        0.40639326593721653,
        0.31519352252160526,
        0.2403981456993558,
        0.36692833148492576,
        0.15957380970815432,
        0.49190636939034976,
        0.3189327704420888,
        0.10792266773848919,
        0.04212798503994418,
        0.39647128156913286,
        0.19546149360621035
    ],
    42: [

    ]
=======
import copy
import shutil
from pathlib import Path

from iotbx.data_manager import DataManager
from cctbx.array_family import flex
import libtbx
from cctbx.maptbx.box import shift_and_box_model
from iotbx.map_model_manager import map_model_manager
from libtbx.utils import null_out
import numpy as np
from scipy.spatial import KDTree



from cctbx.maptbx.qscore import (
  generate_probes_np,
  get_probe_mask,
  shell_probes_progressive,
  shell_probes_precalculate,
  shell_probes_precalculate_flex,
  calc_qscore,
  calc_qscore_flex

)

def isclose_or_nan(a, b, atol=1e-3):
  # perform isclose comparison, treating nans as equal
  return np.isclose(a, b, atol=atol) | (np.isnan(a) & np.isnan(b))


################################################################################
#### Test probe generation
################################################################################

def test_probe_generation():
  # test the primary points generation function against expected data
  atoms_xyz = np.array([[ 5.276, 12.488, 16.069],
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

  probes_xyz = generate_probes_np(atoms_xyz,0.1,8)
  assert np.all(np.isclose(probes_xyz,probes_expected,atol=1e-3))

def test_probe_masking():
  # test the progressive probe masking function against test data
  atoms_xyz = np.array([
      [0,0,-1],
      [0,0,1],
  ])

  # probes_xyz shape (2,4,3), (n_atoms,n_probes,3)
  probes_xyz = np.array([
    [[0,0,-2],
    [0,0,-0.5],
    [0,0,0],
    [0,0,0.5]],

    [[0,0,-2],
    [0,0,-0.5],
    [0,0,0],
    [0,0,0.5]]])

  atom_tree = KDTree(atoms_xyz)

  calculated_result = get_probe_mask(atom_tree,probes_xyz,r=1.4)
  manual_result = np.array([[ True,  True, False, False],
                            [False, False, False,  True]])


  assert np.all(calculated_result==manual_result)



def test_shell_probes():
  # Test full progressive probe generation for a single shell
  atoms_xyz = np.array([[ 5.276, 12.488, 16.069],
                [ 5.649, 13.947, 16.076]])

  # Test progressive

  expected_probes = np.array([[[ 5.276 , 12.488 , 14.569 ],
                          [ 5.0515, 13.4037, 14.9023],
                          [ 4.0297, 12.4397, 15.2357],
                          [ 4.825 , 11.1476, 15.569 ],
                          [ 6.3669, 11.4721, 15.9023],
                          [ 4.0466, 12.6981, 16.9023],
                          [ 5.343 , 11.5476, 17.2357],
                          [ 5.276 , 12.488 , 17.569 ],
                          [ np.nan,  np.nan,  np.nan],
                          [ np.nan,  np.nan,  np.nan]],

                          [[ 5.649 , 13.947 , 14.576 ],
                          [ 5.4245, 14.8627, 14.9093],
                          [ 4.4027, 13.8987, 15.2427],
                          [ 6.7399, 12.9311, 15.9093],
                          [ 7.0245, 14.5216, 16.2427],
                          [ 5.6032, 15.3605, 16.576 ],
                          [ 4.4196, 14.1571, 16.9093],
                          [ 5.716 , 13.0066, 17.2427],
                          [ 5.649 , 13.947 , 17.576 ],
                          [ np.nan,  np.nan,  np.nan]]])
  shell_func = shell_probes_progressive
  probe_xyz, probe_mask = shell_func(
                      atoms_xyz=atoms_xyz,
                      atoms_tree = None,
                      selection_bool=None,
                      n_probes_target=8,
                      n_probes_max=10,
                      RAD=1.5,
                      rtol=0.9,
                      log = null_out())
  assert np.all(isclose_or_nan(probe_xyz,expected_probes,atol=1e-3))

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
                      atoms_xyz=atoms_xyz,
                      atoms_tree = None,
                      selection_bool=None,
                      n_probes_target=8,
                      n_probes_max=10,
                      RAD=1.5,
                      rtol=0.9,
                      log = null_out())

  assert np.all(isclose_or_nan(probe_xyz,expected_probes,atol=1e-3))

  # test precalculate (flex)
  shell_func = shell_probes_precalculate_flex
  probe_xyz,probe_mask  = shell_func(
                      atoms_xyz=flex.vec3_double(atoms_xyz),
                      atoms_tree = None,
                      selection_bool=None,
                      n_probes_target=8,
                      n_probes_max=10,
                      RAD=1.5,
                      rtol=0.9,
                      log = null_out())

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
            "selection_str":None, # Just calculate q score for a sub-selection
            "iselection":None,
            "shells": [0.0,0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.,
                       1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2. ],
            "n_probes_target":8,
            "n_probes_max":16,
            "n_probes_min":4,
            "nproc":4,
            "probe_allocation_method":None,
            "backend":None,
            "debug":True,
            "rtol":0.9,
         }
>>>>>>> Stashed changes
}
# make flex arrays
expected_results = {key: flex.double(val)
                    for key, val in list(expected_results.items())}


def exercise(test_name):
    pdb_file = libtbx.env.find_in_repositories(
        relative_path=f"phenix_regression/real_space_refine/data/tst_{test_name}.pdb",
        test=os.path.isfile)

    map_file = libtbx.env.find_in_repositories(
        relative_path=f"phenix_regression/real_space_refine/data/tst_{test_name}.ccp4",
        test=os.path.isfile)

    result = run_program(
        program_class=qscore.Program,
        args=[pdb_file, map_file,
              "probe_allocation_method=progressive shell_radius_num=5 nproc=1 n_probes=8"],
        logger=null_out(),
    )

    try:
        expected_result = expected_results[test_name]
        assert approx_equal(expected_result, result.qscore, eps=1.e-2)
    except Exception:
        for val in result.qscore:
            print(str(val)+",")
        raise


<<<<<<< Updated upstream
if (__name__ == "__main__"):
    """
    Test random files to verify basic functionality remains unchanged
    Data from phenix_regression/real_space_refine/data
    """
    for test_name in [17]:  # [17,42,48]:
        exercise(test_name)
    print("OK")
=======
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
  if params.backend == "numpy":
     q_func = calc_qscore
  elif params.backend == "flex":
     q_func = calc_qscore_flex

  result = q_func(mmm,params,log=null_out(),debug=True)
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
  test["params"]["n_probes_max"] = 32
  test["params"]["probe_allocation_method"] = "precalculate"
  test["params"]["backend"] = "numpy"
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

  #tst2 flex precalculate
  test = copy.deepcopy(test_template)
  test["data"]["model_file"] = tst2_model_file
  test["data"]["map_file"] = tst2_map_file
  test["data"]["name"] = "tst2_precalc_flex"
  test["data"]["test_dir"] = test_dir
  test["params"]["n_probes_max"] = 32
  test["params"]["probe_allocation_method"] = "precalculate"
  test["params"]["backend"] = "flex"
  test["results"]["expected"]["qscore_per_atom"] = expected_qscore_per_atom
  tests[test["data"]["name"]] = test

  # tst2 progressive (numpy)
  test = copy.deepcopy(test_template)
  test["data"]["model_file"] = tst2_model_file
  test["data"]["map_file"] = tst2_map_file
  test["data"]["name"] = "tst2_progressive_numpy"
  test["data"]["test_dir"] = test_dir
  test["params"]["n_probes_max"] = 16
  test["params"]["n_probes_target"] = 8
  test["params"]["probe_allocation_method"] = "progressive"
  test["params"]["backend"] = "numpy"
  expected_qscore_per_atom= np.array([
    0.81621,  0.79426,  0.83739,  0.67616,  0.75113,  0.81278,  0.75789,
    0.75623,  0.77865,  0.80018,  0.83847,  0.67525,  0.78909,  0.81843,
    0.81285,  0.80816,  0.89982,  0.73878,  0.81402,  0.79254,  0.81353,
    0.81543,  0.64347,  0.75470,  0.84479,  0.80404,  0.77552,  0.75578,
    0.80234,  0.84508,  0.66298,  0.76894,  0.76924,  0.90342,  0.78303,
    0.79723,  0.73253,  0.83709,  0.84759,  0.60567,  0.75056,  0.77734,
    0.89625,  0.85381,  0.74985,  0.69442,  0.80130,  0.82290,  0.57231,
    0.70518,  0.72775,  0.83440,  0.82770,  0.73152,  0.76289,  0.84503,
    0.79984,  0.75651,  0.79504,  0.82964,  0.83653,  0.83177,  0.68769,
    0.79369,  0.83867,  0.67165,  0.78174,  0.85420,  0.73200,  0.82028,
    0.73856,  0.79636,  0.83043,  0.69672,  0.79881,  0.75317,  0.78247,
    0.83621,  0.74694,  0.81975,  0.79633,  0.87402,  0.74882,  0.72080,
    0.87380,  0.74778,
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
  test["params"]["n_probes_max"] = 32
  test["params"]["probe_allocation_method"] = "precalculate"
  test["params"]["backend"] = "numpy"
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


  # precalcualted with flex
  test = copy.deepcopy(test_template)
  test["data"]["model_file"] = yjp_model_file
  test["data"]["name"] = "1yjp_precalc_flex"
  test["data"]["test_dir"] = test_dir
  test["data"]["map_file"] = None
  test["params"]["n_probes_max"] = 32
  test["params"]["probe_allocation_method"] = "precalculate"
  test["params"]["backend"] = "flex"
  test["results"]["expected"]["qscore_per_atom"] = expected_qscore_per_atom
  tests[test["data"]["name"]] = test

  # progressive (numpy)
  test = copy.deepcopy(test_template)
  test["data"]["model_file"] = yjp_model_file
  test["data"]["name"] = "1yjp_progressive_numpy"
  test["data"]["test_dir"] = test_dir
  test["data"]["map_file"] = None
  test["params"]["n_probes_max"] = 16
  test["params"]["n_probes_target"] = 8
  test["params"]["probe_allocation_method"] = "progressive"
  test["params"]["backend"] = "numpy"
  expected_qscore_per_atom= np.array([
    0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,
    0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.07689,
    0.19082,  0.02124,  0.00000,  0.74434,  0.93599,  0.93209,  0.77202,
    0.91962,  0.91928,  0.96171,  0.91251,  0.91642,  0.36156,  0.31734,
    0.00000,  0.94018,  0.90381,  0.93852,  0.94966,  0.93188,  0.91575,
    0.91731,  0.97200,  0.95272,  0.94102,  0.89297,  0.91385,  0.96626,
    0.89274,  0.91538,  0.96871,  0.96205,  0.93226,  0.89114,  0.92745,
    0.94155,  0.90361,  0.92891,  0.74236,  0.92431,  0.28466,  0.91468,
    0.93081,  0.96238,  0.93960,  0.00000,  0.97779,  0.00000,  0.94679,
    0.97174,  0.89276,  0.90919,
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
  test["params"]["n_probes_max"] = 32
  test["params"]["probe_allocation_method"] = "precalculate"
  test["params"]["backend"] = "numpy"
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

  # precalculated with flex
  test = copy.deepcopy(test_template)
  test["data"]["model_file"] = str(boxed_1yjp_path)
  test["data"]["name"] = "1yjp_boxed_precalc_flex"
  test["data"]["test_dir"] = test_dir
  test["data"]["map_file"] = None
  test["params"]["n_probes_max"] = 32
  test["params"]["probe_allocation_method"] = "precalculate"
  test["params"]["backend"] = "flex"
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
  test["params"]["n_probes_max"] = 32
  test["params"]["probe_allocation_method"] = "precalculate"
  test["params"]["backend"] = "numpy"
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

  # tst2, shift_and_box, precalculate, flex
  test = copy.deepcopy(test_template)
  test["data"]["model_file"] = str(boxed_tst2_path)
  test["data"]["map_file"] = tst2_map_file
  test["data"]["name"] = "tst2_boxed_precalc_flex"
  test["data"]["test_dir"] = test_dir
  test["params"]["n_probes_max"] = 32
  test["params"]["probe_allocation_method"] = "precalculate"
  test["params"]["backend"] = "flex"
  test["results"]["expected"]["qscore_per_atom"] = expected_qscore_per_atom
  tests[test["data"]["name"]] = test

  # tst2, shift_and_box, progressive, numpy
  test = copy.deepcopy(test_template)
  test["data"]["model_file"] = str(boxed_tst2_path)
  test["data"]["map_file"] = tst2_map_file
  test["data"]["name"] = "tst2_boxed_progressive_numpy"
  test["data"]["test_dir"] = test_dir
  test["params"]["n_probes_max"] = 16
  test["params"]["n_probes_target"] = 8
  test["params"]["probe_allocation_method"] = "progressive"
  test["params"]["backend"] = "numpy"
  expected_qscore_per_atom = np.array([
    0.81621,  0.79426,  0.83739,  0.67616,  0.75113,  0.81278,  0.75789,
    0.75623,  0.77865,  0.80018,  0.83847,  0.67525,  0.78909,  0.81843,
    0.81285,  0.80816,  0.89982,  0.73878,  0.81402,  0.79254,  0.81353,
    0.81543,  0.64347,  0.75470,  0.84479,  0.80404,  0.77552,  0.75578,
    0.80234,  0.84508,  0.66298,  0.76894,  0.76924,  0.90342,  0.78303,
    0.79723,  0.73253,  0.83709,  0.84759,  0.60567,  0.75056,  0.77734,
    0.89625,  0.85381,  0.74985,  0.69442,  0.80130,  0.82290,  0.57231,
    0.70518,  0.72775,  0.83440,  0.82770,  0.73152,  0.76289,  0.84503,
    0.79984,  0.75651,  0.79504,  0.82964,  0.83653,  0.83177,  0.68769,
    0.79369,  0.83867,  0.67165,  0.78174,  0.85420,  0.73200,  0.82028,
    0.73856,  0.79636,  0.83043,  0.69672,  0.79881,  0.75317,  0.78247,
    0.83621,  0.74694,  0.81975,  0.79633,  0.87402,  0.74882,  0.72080,
    0.87380,  0.74778,

  ])
  test["results"]["expected"]["qscore_per_atom"] = expected_qscore_per_atom
  tests[test["data"]["name"]] = test


  return tests




if (__name__ == "__main__"):


  #1. test probe generation
  test_probe_generation()

  #2. test probe masking (for progressive)
  test_probe_masking()

  #3. test single shell probe generation
  test_shell_probes()



  # Test on some real models
  tests = build_tests()
  for test_name,test in tests.items():
     test = run_test(test)
>>>>>>> Stashed changes
