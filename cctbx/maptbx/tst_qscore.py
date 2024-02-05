from __future__ import absolute_import, division, print_function
import os
import copy
import subprocess
import shutil
from pathlib import Path
import json
import argparse

from iotbx.data_manager import DataManager
from cctbx.array_family import flex
import libtbx
from libtbx import group_args
from libtbx.utils import null_out
from iotbx.cli_parser import run_program
import numpy as np
from scipy.spatial import KDTree



from cctbx.maptbx.qscore import (
  generate_probes_np,
  SpherePtsVectorized,
  get_probe_mask,
  get_probes,
  shell_probes_progressive,
  shell_probes_precalculate,
  _shell_probes_progressive_wrapper,
  _shell_probes_precalculate_wrapper
)
from cctbx.programs.qscore import Program as QscoreProgram




def isclose_or_nan(a, b, atol=1e-3):
  # perform isclose comparison, treating nans as equal
  return np.isclose(a, b, atol=atol) | (np.isnan(a) & np.isnan(b))

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


  # test that our points generator functions don't diverge on
  # a large amount of points
  points = np.random.random((1000,3))
  rads = np.linspace(0,2.0,20)
  Ns = [2]
  # test point by point
  for point in points:
    for rad in rads:
      for N in Ns:
        mapq_values = SpherePtsVectorized(point[None,:],rad,N)
        cctbx_values = generate_probes_np(point[None,:],rad,N)
        assert np.all(np.isclose(mapq_values,cctbx_values,atol=1e-3))

  # test vectorized over points
  for rad in rads:
    for N in Ns:
      mapq_values = SpherePtsVectorized(points,rad,N)
      cctbx_values = generate_probes_np(points,rad,N)
      assert np.all(np.isclose(mapq_values,cctbx_values,atol=1e-3))


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

def test_shell_probes(probe_allocation_method="progressive"):
  # Test full progressive probe generation for a single shell
  atoms_xyz = np.array([[ 5.276, 12.488, 16.069],
                [ 5.649, 13.947, 16.076]])

  if probe_allocation_method=="progressive":
    shell_func = shell_probes_progressive
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

  elif probe_allocation_method == "precalculate":
    shell_func = shell_probes_precalculate
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

  probes_xyz, _ = shell_func(atoms_xyz=atoms_xyz,
                      atoms_tree = None,
                      selection=None,
                      n_probes_target=8,
                      n_probes_max=10,
                      RAD=1.5,
                      rtol=0.9,
                      log = null_out())

  if probe_allocation_method == "precalculate":
    print(probes_xyz)
  assert np.all(isclose_or_nan(probes_xyz,expected_probes,atol=1e-3))

def test_get_probes(probe_allocation_method="progressive"):
  # Test the full progressive probe generation for multiple shells

  atoms_xyz = np.array([[ 5.276, 12.488, 16.069],
                    [ 5.649, 13.947, 16.076]])

  params = group_args(
    selection=None,
    shells = np.array([0.0,0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3,
          1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2. ]),
    n_probes_target=8,
    n_probes_max=16,
    n_probes_min=4,
    rtol=0.9,
    nproc=1,
    probe_allocation_method=probe_allocation_method,
    log = null_out()
    )

  if probe_allocation_method == "progressive":
    worker_func = _shell_probes_progressive_wrapper
    sum_expected = 12343.8878
  elif probe_allocation_method == "precalculate":
    worker_func = _shell_probes_precalculate_wrapper
    sum_expected = 23337.6710

  probes_xyz,probe_mask = get_probes(atoms_xyz=atoms_xyz,
                      params=params,
                      worker_func=worker_func,
                      log = null_out())


  sum_calc = probes_xyz[~np.isnan(np.around(probes_xyz,3))].sum()


  print(probe_allocation_method,sum_calc)
  assert np.all(np.isclose(sum_calc,sum_expected)), (
    "The sum of generated probes do not match previous values obtained from mapq. Debug probe generation."
  )


################################
# Tests with real data
################################

def prepare_test_data(templates):
  """
  Prepare folders with data files for each test
  """

  for i,template in enumerate(templates):
    test_dir_name = template["data"]["name"]

    is_fragment = False
    if template["data"]["fragment_iselection"] not in [None,[]]:
      i_sel = template["data"]["fragment_iselection"]
      is_fragment = True
    test_dir = Path(test_dir_name)

    # make test directory and copy data
    test_dir.mkdir(exist_ok=True)
    templates[i]["data"]["test_dir"] = str(test_dir.absolute())
    model_path_base = Path(template["data"]["model_file_base"])
    map_path_base = Path(template["data"]["map_file_base"])
    model_path = test_dir / Path(f"model.pdb")
    map_path = test_dir / Path(f"map.ccp4")
    shutil.copyfile(model_path_base,model_path)
    shutil.copyfile(map_path_base,map_path)

    # possibly  make fragment files
    if is_fragment:
      dm = DataManager()
      dm.process_model_file(str(model_path))
      model = dm.get_model()
      sel = np.full(model.get_number_of_atoms(),False)
      sel[i_sel] = True
      model_sel = model.select(flex.bool(sel))
      dm.write_model_file(model_sel.model_as_pdb(),str(model_path),overwrite=True)

    # record in template the data in test directory
    templates[i]["data"]["model_file"] = str(model_path)
    templates[i]["data"]["map_file"] = str(map_path)

  return templates

def run_test_template_mapq(template,
                          mapq_location=None,
                          mapq_debug_data_filename="debug_data_mapq.json"):
  """
  Run mapq via Chimera on the command line. Load debug results
  NOTE: The debug results rely on a modified version of mapq. The releaseed version does not write
  all the intermediate data.
  """
  assert mapq_location is not None
  mapq_location = Path(mapq_location)
  model_path = Path(template["data"]["model_file"])
  map_path = Path(template["data"]["map_file"])
  test_dir = Path(template["data"]["test_dir"])

  # run program
  mapq_executable = mapq_location / Path("mapq_cmd.py")
  chimera_path = mapq_location / Path("../../../../../Chimera.app")
  mapq_command = f"python {mapq_executable.absolute()} {chimera_path.absolute()}  map={map_path.absolute()} pdb={model_path.absolute()} bfactor=1.0"
  # q is stored in bfactor column of pdb, with:  bfactor = f * (1.0-Qscore). The scale factor f is provided with the 'bfactor' arg
  # Bfactor 'f' should not affect q result
  print(f"Running mapq  with command:")
  print(mapq_command)
  print("\n\n")
  subprocess.run(mapq_command.split())
  debug_data = load_mapq_debug_data(test_dir)
  return debug_data



def load_mapq_debug_data(test_dir):
    """
    Load the mapq intermediate results. Looks for a file 'debug_data_mapq.json' in
    the test directory.
    """

    # anticipate output data file path
    data_file = Path(test_dir,Path("debug_data_mapq.json")).absolute()
    # load debug data
    with open(data_file,"r") as fh:
      debug_data = json.load(fh)
      def _is_ragged(a):
        # don't force arrays for ragged data
        if isinstance(a, list):
          # Check if all elements are lists and have the same length
          if all(isinstance(i, list) for i in a):
            length = len(a[0])
            return any(len(i) != length for i in a)
          else:
            # It's a list, but not a list of lists
            return False
        else:
          # Not a list, so it's not a ragged array in the typical sense
          return False

      debug_data = {key:np.array(value) if not _is_ragged(value) else value for key,value in debug_data.items()}
    return debug_data


def run_test_template_cctbx(template):
  """
  Run a test using the progressive method added to cctbx
  """


  params = group_args(
    selection=None,
    shells = template["params"]["shells"],
    n_probes_target=template["params"]["n_probes_target"],
    n_probes_max=template["params"]["n_probes_max"],
    n_probes_min=template["params"]["n_probes_min"],
    rtol=template["params"]["rtol"],
    nproc=template["params"]["nproc"],
    probe_allocation_method = template["params"]["probe_allocation_method"],
    log = null_out(),
  )
  model_filename = template["data"]["model_file"]
  map_filename = template["data"]["map_file"]

  param_args = [f"{key}={getattr(params,key)}" for key in params.keys() if key not in  ["shells","log"]]
  for shell in params.shells:
    param_args.append(f"shells={shell}")
  args = [model_filename,map_filename, "debug=True"] + param_args
  print(args)
  result = run_program(program_class=QscoreProgram,args=args)
  result = {key:getattr(result,key) for key in result.keys()} # group_args to dict
  return result

def run_template(template,mapq_location=None):
  print("Template")
  print(json.dumps(template,indent=2))

  # get data
  debug_data = run_test_template_cctbx(template)
  probe_xyz = debug_data["probe_xyz"]
  probe_mask = debug_data["probe_mask"]
  d_vals = debug_data["d_vals"]
  g_vals = debug_data["g_vals"]
  qscore_per_atom = debug_data["qscore_per_atom"]

  # Record sums
  template["results"]["probe_sum"] = probe_xyz[~np.isnan(probe_xyz)].sum()
  template["results"]["q_sum"] = qscore_per_atom.sum()

  if template["results"]["probe_sum_expected"] is not None:
    assert np.isclose(template["results"]["probe_sum"],template["results"]["probe_sum_expected"],atol=1e-2)
  if template["results"]["q_sum_expected"]:
    assert np.isclose(template["results"]["q_sum"],template["results"]["q_sum_expected"],atol=1e-2)

  # Record all data
  template["results"]["cctbx"] = debug_data

  # run mapq, check results with progressive
  if mapq_location is not None:
    if template["params"]["probe_allocation_method"] == "progressive":
      debug_data_mapq = run_test_template_mapq(template,mapq_location = mapq_location)

      probe_xyz_mapq = debug_data_mapq["probe_xyz"]
      probe_mask_mapq = debug_data_mapq["probe_mask"]
      d_vals_mapq = debug_data_mapq["d_vals"]
      g_vals_mapq = debug_data_mapq["g_vals"]
      qscore_per_atom_mapq = debug_data_mapq["qscore_per_atom"]

      # Check probes
      assert np.all(isclose_or_nan(probe_xyz,probe_xyz_mapq))
      assert np.all(isclose_or_nan(probe_mask,probe_mask_mapq))

      # Check d and g
      assert np.all(isclose_or_nan(d_vals,d_vals_mapq))
      assert np.all(isclose_or_nan(g_vals,g_vals_mapq))

      # Check q
      assert np.all(isclose_or_nan(qscore_per_atom,qscore_per_atom_mapq))

      # Check q from actual pdb output file
      test_dir = Path(template["data"]["test_dir"])
      model_path = Path(template["data"]["model_file"])
      map_path = Path(template["data"]["map_file"])
      dm = DataManager()
      mapq_output_file = Path(test_dir, f"{model_path.stem}.pdb__Q__{map_path.stem}.ccp4.pdb")
      _ = dm.process_model_file(str(mapq_output_file))
      model = dm.get_model()
      pseudo_b = model.get_b_iso().as_numpy_array()
      q_test = pseudo_b
      # Round heavily to match bfactor
      q_calc = np.around(qscore_per_atom,decimals=2)
      assert set(q_test)==set(q_calc)

      # record sums
      template["results"]["probe_sum"] = probe_xyz_mapq[~np.isnan(probe_xyz)].sum()
      template["results"]["q_sum"] = qscore_per_atom_mapq.sum()
      template["results"]["mapq"] = debug_data_mapq
  return template


# a template to store configuration for a single test
test_template ={

         "data":{
          "name":None,
          "model_file":None,
          "model_file_base":None,
          "map_file":None,
          "map_file_base":None,
          "test_dir":None,
          "fragment_iselection":None, # actually reduce the file
         },
         "results":{
          "probe_sum":None,
          "q_sum":None,
          "probe_sum_expected":None,
          "q_sum_expected":None,
         },
         "params":{
            "selection":None, # Just calculate q score for a sub-selection
            "iselection":None,
            "shells": [0.0,0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2. ],
            "n_probes_target":8,
            "n_probes_max":16,
            "n_probes_min":4,
            "nproc":4,
            "probe_allocation_method":"progressive",
            "rtol":0.9,
         }
}

if (__name__ == "__main__"):
  parser = argparse.ArgumentParser(description="Run qscore tests")

  # Figure out if using mapq
  parser.add_argument('--mapq_location',
                      type=str,
                      help='Compare to mapq results. Example: /Users/user/Desktop/Chimera.app/Contents/Resources/share/mapq')
  args = parser.parse_args()
  mapq_location = args.mapq_location
  if mapq_location is None:
    mapq = False
  else:
    mapq = True


  # Start tests...

  #1. test probe generation
  test_probe_generation()

  #2. test probe masking (for progressive)
  test_probe_masking()

  #3. test single shell probe generation
  test_shell_probes(probe_allocation_method="progressive")
  test_shell_probes(probe_allocation_method="precalculate")

  #4. test multi-shell probe generation
  test_get_probes(probe_allocation_method="progressive")
  test_get_probes(probe_allocation_method="precalculate")


  #5. test some real files
  templates = []


  def get_base_map_model(test_name):
    # model
    pdb_file_repo = libtbx.env.find_in_repositories(
      relative_path=f"phenix_regression/real_space_refine/data/tst_{test_name}.pdb",
      test=os.path.isfile)

    # map
    map_file_repo = libtbx.env.find_in_repositories(
      relative_path=f"phenix_regression/real_space_refine/data/tst_{test_name}.ccp4",
      test=os.path.isfile)
    return pdb_file_repo,map_file_repo

  # Do some small fragments (progressive)
  base_model_name, base_map_name = get_base_map_model(17)
  test_name = Path(base_model_name).stem

  # 1 atom fragment
  test = copy.deepcopy(test_template)
  test["data"]["name"] = f"{test_name}_0"
  test["data"]["model_file_base"] = base_model_name
  test["data"]["map_file_base"] = base_map_name
  test["data"]["fragment_iselection"] = [0]
  templates.append(test)

  # 2 atom fragment
  test = copy.deepcopy(test_template)
  test["data"]["name"] = f"{test_name}_0-1"
  test["data"]["model_file_base"] = base_model_name
  test["data"]["map_file_base"] = base_map_name
  test["data"]["fragment_iselection"] = [0,1]
  templates.append(test)

  # 3 atom fragment
  test = copy.deepcopy(test_template)
  test["data"]["name"] = f"{test_name}_0-2"
  test["data"]["model_file_base"] = base_model_name
  test["data"]["map_file_base"] = base_map_name
  test["data"]["fragment_iselection"] = [0,1,2]
  templates.append(test)

  # 4 atom fragment
  test = copy.deepcopy(test_template)
  test["data"]["name"] = f"{test_name}_0-3"
  test["data"]["model_file_base"] = base_model_name
  test["data"]["map_file_base"] = base_map_name
  test["data"]["fragment_iselection"] = [0,1,2,3]
  templates.append(test)

  # full molecules

  # 17
  # progressive
  base_model_name, base_map_name = get_base_map_model(17)
  test_name = Path(base_model_name).stem
  test = copy.deepcopy(test_template)
  test["data"]["name"] = test_name+"_progressive"
  test["data"]["model_file_base"] = base_model_name
  test["data"]["map_file_base"] = base_map_name
  test["params"]["probe_allocation_method"] = "progressive"
  test["results"]["q_sum_expected"] = 85.74585520218474
  test["results"]["probe_sum_expected"] = 7033351.9964000005
  templates.append(test)

  # precalculate
  test = copy.deepcopy(test_template)
  test["data"]["name"] = test_name+"_precalc"
  test["data"]["model_file_base"] = base_model_name
  test["data"]["map_file_base"] = base_map_name
  test["params"]["n_probes_max"] = 128
  test["params"]["probe_allocation_method"] = "precalculate"
  test["results"]["q_sum_expected"] = 87.87610388144196
  test["results"]["probe_sum_expected"] = 108903540.12935309
  templates.append(test)

  # 42
  # progressive
  base_model_name, base_map_name = get_base_map_model(42)
  test_name = Path(base_model_name).stem
  test = copy.deepcopy(test_template)
  test["data"]["name"] = test_name+"_progressive"
  test["data"]["model_file_base"] = base_model_name
  test["data"]["map_file_base"] = base_map_name
  test["results"]["q_sum_expected"] = 65.18062030153557
  test["results"]["probe_sum_expected"] = 479151.29789999995
  templates.append(test)

  # precalculate
  test = copy.deepcopy(test_template)
  test["data"]["name"] = test_name+"_precalc"
  test["data"]["model_file_base"] = base_model_name
  test["data"]["map_file_base"] = base_map_name
  test["params"]["n_probes_max"] = 128
  test["params"]["probe_allocation_method"] = "precalculate"
  test["results"]["q_sum_expected"] = 65.98930853051141
  test["results"]["probe_sum_expected"] = 7393520.8060661685
  templates.append(test)

  # 48
  # progressive
  base_model_name, base_map_name = get_base_map_model(48)
  test_name = Path(base_model_name).stem
  test = copy.deepcopy(test_template)
  test["data"]["name"] = test_name+"_progressive"
  test["data"]["model_file_base"] = base_model_name
  test["data"]["map_file_base"] = base_map_name
  test["results"]["q_sum_expected"] = 74.57025665508365
  test["results"]["probe_sum_expected"] = 904759.264
  templates.append(test)

  # precalculate
  test = copy.deepcopy(test_template)
  test["data"]["name"] = test_name+"_precalc"
  test["data"]["model_file_base"] = base_model_name
  test["data"]["map_file_base"] = base_map_name
  test["params"]["n_probes_max"] = 128
  test["params"]["probe_allocation_method"] = "precalculate"
  test["results"]["q_sum_expected"] = 75.99978032145603
  test["results"]["probe_sum_expected"] = 13979341.58057615
  templates.append(test)

  # prepare files
  templates = prepare_test_data(templates)

  # run templates
  templates = [run_template(template,mapq_location=mapq_location) for template in templates]


  # Compare progressive and precalculate
  for template in templates:
    name = template["data"]["name"]
    for other_template in templates:
      other_name = other_template["data"]["name"]

      if "progressive" in name:
        if name.strip("progressive") == other_name.strip("precalc"):
          debug_data_cctbx = template["results"]["cctbx"]
          debug_data_cctbx_mp = other_template["results"]["cctbx"]
          probe_xyz = debug_data_cctbx["probe_xyz"]
          probe_xyz_mp = debug_data_cctbx_mp["probe_xyz"]
          qscore_per_atom = debug_data_cctbx["qscore_per_atom"]
          qscore_per_atom_mp = debug_data_cctbx_mp["qscore_per_atom"]

          rmsd = np.sqrt(np.mean((qscore_per_atom-qscore_per_atom_mp) ** 2))
          cc = np.corrcoef(qscore_per_atom, qscore_per_atom_mp)[0][1]
          print(template["data"]["name"],"rmsd",rmsd)
          print(template["data"]["name"],"CC",cc)
          assert rmsd<0.1, rmsd
          assert cc>0.9, cc

  print("OK")
