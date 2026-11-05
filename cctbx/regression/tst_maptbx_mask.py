from __future__ import absolute_import, division, print_function
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex
from iotbx.map_model_manager import map_model_manager
from libtbx import group_args
import libtbx.load_env

def get_random_structure_and_map(
   use_static_structure = False,
   random_seed = 171413,
  ):

  if use_static_structure:
    mmm = map_model_manager()
    mmm.generate_map()
    return group_args(model = mmm.model(), mm = mmm.map_manager())
  import random
  random.seed(random_seed)
  i = random.randint(1, 714717)
  flex.set_random_seed(i)

  xrs = random_structure.xray_structure(
    space_group_info = space_group_info(19),
    volume_per_atom  = 25.,
    elements         = ('C', 'N', 'O', 'H')*10,
    min_distance     = 1.5)
  fc = xrs.structure_factors(d_min = 2).f_calc()
  fft_map = fc.fft_map(resolution_factor = 0.25)
  fft_map.apply_volume_scaling()
  ph = iotbx.pdb.input(
    source_info = None, lines = xrs.as_pdb_file()).construct_hierarchy()
  ph.atoms().set_xyz(xrs.sites_cart())
  map_data = fft_map.real_map_unpadded()
  mm = map_manager(
    unit_cell_grid             = map_data.accessor().all(),
    unit_cell_crystal_symmetry = fc.crystal_symmetry(),
    origin_shift_grid_units    = (0, 0, 0),
    map_data                   = map_data)
  model = mmtbx.model.manager(
    model_input = None, pdb_hierarchy = ph, crystal_symmetry = fc.crystal_symmetry())
  return group_args(model = model, mm = mm)

def check_masked_map(new_map, orig_map, mask, outside_value):
  """New map must blend orig_map and outside_value by the mask value"""
  expected = mask * orig_map + (1 - mask) * outside_value
  assert approx_equal(flex.max(flex.abs((new_map - expected).as_1d())), 0), \
    "masked map is not mask*map + (1-mask)*outside"

def exercise_mask():
  mam = get_random_structure_and_map(use_static_structure = True)
  mm_orig = mam.mm.deep_copy()
  map_data_orig   = mam.mm.map_data().deep_copy()
  assert map_data_orig.origin() == (0, 0, 0)

  # Create mask with gridding supplied
  from cctbx.maptbx.mask import create_mask_around_atoms
  cm = create_mask_around_atoms(model = mam.model,
     mask_atoms_atom_radius = 3,
     n_real = mam.mm.map_data().all(),
     wrapping = mam.mm.wrapping())

  new_mask = cm.mask()
  assert (new_mask.count(0), new_mask.count(1), new_mask.size()) == (31852, 6548, 38400)

  # Soft mask: smooth, positive everywhere, rescaled to a maximum of 1
  cm.soft_mask(soft_mask_radius = 5)
  new_mask = cm.mask()
  assert flex.min(new_mask.as_1d()) > 0
  assert approx_equal(flex.max(new_mask.as_1d()), 1), \
    "soft mask is not rescaled to a maximum of 1"

  # Create mask starting with a map_manager object
  cm = create_mask_around_atoms(model = mam.model,
     mask_atoms_atom_radius = 3,
     map_manager = mam.mm)
  assert mm_orig.map_data().origin() == (0, 0, 0)
  mask = cm.mask()
  assert (mask.count(0), mask.count(1)) == (31852, 6548)
  inside = (mask > 0.5).as_1d()
  mean_inside = flex.mean(map_data_orig.as_1d().select(inside))

  new_mm = cm.apply_mask_to_other_map_manager(mm_orig)
  check_masked_map(new_mm.map_data(), map_data_orig, mask, 0)

  new_mm = cm.apply_mask_to_other_map_manager(mm_orig,
     set_outside_to_mean_inside = True)
  check_masked_map(new_mm.map_data(), map_data_orig, mask, mean_inside)

  # Create mask using n_real and check applying mask
  cm = create_mask_around_atoms(model = mam.model,
     mask_atoms_atom_radius = 3,
     n_real = mam.mm.map_data().all(),
     wrapping = mam.mm.wrapping())
  assert mm_orig.map_data().origin() == (0, 0, 0)

  new_mm = cm.apply_mask_to_other_map_manager(mm_orig,
     set_outside_to_mean_inside = True)
  check_masked_map(new_mm.map_data(), map_data_orig, cm.mask(), mean_inside)

  # Create mask using map_manager object
  mm = mm_orig.deep_copy()
  mm.create_mask_around_atoms(
    model = mam.model,
    mask_atoms_atom_radius = 3)
  hard_mask = mm.get_mask_as_map_manager().map_data()
  assert (hard_mask.count(0), hard_mask.count(1)) == (31852, 6548)
  mm.soft_mask(soft_mask_radius = 5)
  soft_mask = mm.get_mask_as_map_manager().map_data().deep_copy()
  mm.apply_mask()
  check_masked_map(mm.map_data(), map_data_orig, soft_mask, 0)

  # create mask around edges
  mm = mm_orig.deep_copy()
  mm.create_mask_around_edges(boundary_radius = 5)
  edge_mask = mm.get_mask_as_map_manager().map_data()
  assert (edge_mask.count(0), edge_mask.count(1)) == (17184, 21216)
  mm.soft_mask(soft_mask_radius = 5)
  soft_mask = mm.get_mask_as_map_manager().map_data().deep_copy()
  mm.apply_mask()
  check_masked_map(mm.map_data(), map_data_orig, soft_mask, 0)

  # The rest uses create_mask_around_density with a sequence, which needs
  # phenix.loop_lib to estimate the solvent content
  if not libtbx.env.has_module("phenix"):
    print("phenix is not available, skipping create_mask_around_density")
    return

  # auto-generate mask
  mm = mm_orig.deep_copy()
  mm.create_mask_around_density(
    resolution = 3,
    molecular_mass = 2000,
    sequence = 'AVAGS',
    solvent_content = None)
  mask_mm = mm.get_mask_as_map_manager()
  mask_data = mask_mm.map_data().deep_copy()
  print ("Mask zero/one: ",
     mask_data.count(0), mask_data.count(1), mask_data.size())
  assert (mask_data.count(0), mask_data.count(1), mask_data.size()) == \
    (36705, 1695, 38400)

  mm.apply_mask()
  check_masked_map(mm.map_data(), map_data_orig, mask_data, 0)


  # create zero_boundary_mask

  map_data=mm.map_data()
  mean=map_data.as_1d().min_max_mean().mean
  sd=map_data.as_1d().standard_deviation_of_the_sample()

  from cctbx.maptbx import binary_filter
  bf=binary_filter(map_data,mean+sd).result()
  direct_mask_mm=mm.customized_copy(map_data=bf)

  map_data=mm.map_data()
  mean=map_data.as_1d().min_max_mean().mean
  sd=map_data.as_1d().standard_deviation_of_the_sample()

  manager_mask_mm=mm.deep_copy()
  manager_mask_mm.binary_filter(mean+sd)  # filters in place
  assert approx_equal(direct_mask_mm.map_data(), manager_mask_mm.map_data())
  mam=map_model_manager(map_manager_1=direct_mask_mm,map_manager_2=manager_mask_mm)
  assert approx_equal(mam.map_map_cc(map_id='map_manager_1',other_map_id='map_manager_2'),1)


if (__name__  ==  "__main__"):
  exercise_mask()
