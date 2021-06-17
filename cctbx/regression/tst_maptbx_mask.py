from __future__ import absolute_import, division, print_function
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex
from iotbx.map_model_manager import map_model_manager
from libtbx import group_args

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

def exercise_mask():
  mam = get_random_structure_and_map(use_static_structure = True)
  mm_orig = mam.mm.deep_copy()
  map_data_orig   = mam.mm.map_data().deep_copy()
  assert map_data_orig.origin() == (0, 0, 0)
  sites_frac_orig = mam.model.get_sites_frac().deep_copy()
  sites_cart_orig = mam.model.get_sites_cart().deep_copy()
  cs_orig         = mam.model.crystal_symmetry()

  # Create mask with gridding supplied
  from cctbx.maptbx.mask import create_mask_around_atoms
  cm = create_mask_around_atoms(model = mam.model,
     mask_atoms_atom_radius = 3,
     n_real = mam.mm.map_data().all(),
     wrapping = mam.mm.wrapping())

  new_mask = cm.mask()
  assert (new_mask.count(0), new_mask.count(1), new_mask.size()) == (31852, 6548, 38400)

  cm.soft_mask(soft_mask_radius = 5)
  new_mask = cm.mask()
  assert (new_mask.count(0), new_mask.count(1), new_mask.size()) == (0, 0, 38400)
  assert approx_equal( (new_mask[0], new_mask[100]), (0.0232240949563, 0.0216530808725))

  # Create mask starting with a map_manager object
  cm = create_mask_around_atoms(model = mam.model,
     mask_atoms_atom_radius = 3,
     map_manager = mam.mm)
  assert mm_orig.map_data().origin() == (0, 0, 0)

  new_mm = cm.apply_mask_to_other_map_manager(mm_orig)
  new_map_data = new_mm.map_data()
  assert approx_equal( (new_map_data[5215], new_map_data[8432]), (0.0, -0.0514926330954))

  new_mm = cm.apply_mask_to_other_map_manager(mm_orig,
     set_outside_to_mean_inside = True)
  new_map_data = new_mm.map_data()
  assert approx_equal( (new_map_data[5215], new_map_data[8432]), (0.0646441178021, -0.0514926330954))

  # Create mask using n_real and check applying mask
  cm = create_mask_around_atoms(model = mam.model,
     mask_atoms_atom_radius = 3,
     n_real = mam.mm.map_data().all(),
     wrapping = mam.mm.wrapping())
  assert mm_orig.map_data().origin() == (0, 0, 0)

  new_mm = cm.apply_mask_to_other_map_manager(mm_orig,
     set_outside_to_mean_inside = True)
  new_map_data = new_mm.map_data()
  assert approx_equal( (new_map_data[5215], new_map_data[8432]), (0.0646441178021, -0.0514926330954))

  # Create mask using map_manager object
  mm = mm_orig.deep_copy()
  orig_map = mm.map_data().deep_copy()
  mm.create_mask_around_atoms(
    model = mam.model,
    mask_atoms_atom_radius = 3)
  mm.soft_mask(soft_mask_radius = 5)
  mm.apply_mask()
  new_map = mm.map_data()
  assert approx_equal(
    ( orig_map[5215], orig_map[8432], new_map[5215], orig_map[8432]),
    (-0.012385224860014447, -0.05149263309541403, -0.0004143453062615644,
        -0.05149263309541403)
    )

  # create mask around edges
  mm = mm_orig.deep_copy()
  orig_map = mm.map_data().deep_copy()
  mm.create_mask_around_edges(
    soft_mask_radius = 5)
  mm.soft_mask(soft_mask_radius = 5)
  mm.apply_mask()
  new_map = mm.map_data()
  assert approx_equal( (orig_map[0], new_map[0]),
     (-0.0071372973593,-0.00204624804668))
  assert approx_equal(
   ( orig_map[5215], orig_map[8432], new_map[5215], orig_map[8432]),
  (-0.01238522486, -0.0514926330954, -0.00401527420775, -0.0514926330954))


  # auto-generate mask
  mm = mm_orig.deep_copy()
  orig_map = mm.map_data().deep_copy()
  mm.create_mask_around_density(
    resolution = 3,
    molecular_mass = 2000,
    sequence = 'AVAGS',
    solvent_content = None)
  mask_mm = mm.get_mask_as_map_manager()
  mask_data = mask_mm.map_data()
  print ("Mask zero/one: ",
     mask_data.count(0), mask_data.count(1), mask_data.size())
  assert approx_equal ((mask_data.count(0), mask_data.count(1), mask_data.size()),
    (36705, 1695, 38400))

  mm.apply_mask()
  new_map = mm.map_data()
  print( "Map values before/after mask:",
       orig_map[4322], orig_map[9680], new_map[4322], new_map[9680]),
  assert approx_equal(
     (orig_map[4322], orig_map[9680], new_map[4322], new_map[9680]),
    (-0.0214128152846, -0.0249896972752, 0.0, -0.0249896972752))


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

  manager_mask_mm=mm.binary_filter(mean+sd)
  mam=map_model_manager(map_manager_1=direct_mask_mm,map_manager_2=manager_mask_mm)
  assert approx_equal(mam.map_map_cc(map_id='map_manager_1',other_map_id='map_manager_2'),1)


if (__name__  ==  "__main__"):
  exercise_mask()
