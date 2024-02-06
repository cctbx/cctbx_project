from __future__ import absolute_import, division, print_function
import os
<<<<<<< Updated upstream
from cctbx.array_family import flex
import libtbx
from libtbx import group_args
from libtbx.utils import null_out
from iotbx.cli_parser import run_program
import numpy as np
from scipy.spatial import KDTree



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
