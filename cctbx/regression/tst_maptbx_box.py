from __future__ import absolute_import, division, print_function
import os
import iotbx.pdb
from libtbx.test_utils import approx_equal
from cctbx.sgtbx import space_group_info
from cctbx.development import random_structure
import cctbx.maptbx.box
from libtbx import group_args
import iotbx.pdb
from iotbx.map_manager import map_manager
from iotbx.map_model_manager import map_model_manager
from iotbx.data_manager import DataManager
import mmtbx.model
from scitbx.array_family import flex
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

def exercise_around_model():

  from cctbx.maptbx.box import make_list_symmetric
  a=[3,4,5,3,9,1,6,3,2,5,6,6]
  new_a=make_list_symmetric(a)
  from scitbx.array_family import flex
  aa=flex.double(a)
  new_aa=flex.double(new_a)
  assert (aa.size(),new_aa.size())== (12, 12)
  assert aa.min_max_mean().mean == new_aa.min_max_mean().mean
  print (a,new_a)

  a=[3,4,5,3,8,1,6,7,3,2,5,6,6]
  new_a=make_list_symmetric(a)
  from scitbx.array_family import flex
  aa=flex.double(a)
  new_aa=flex.double(new_a)
  print (a,new_a)
  assert (aa.size(),new_aa.size())== (13, 13)
  assert aa.min_max_mean().mean == new_aa.min_max_mean().mean

  mam = get_random_structure_and_map(use_static_structure = True)

  map_data_orig   = mam.mm.map_data().deep_copy()
  sites_frac_orig = mam.model.get_sites_frac().deep_copy()
  sites_cart_orig = mam.model.get_sites_cart().deep_copy()
  cs_orig         = mam.model.crystal_symmetry()

  box = cctbx.maptbx.box.around_model(
    map_manager = mam.mm,
    model       = mam.model.deep_copy(),
    box_cushion     = 10,
    wrapping    = True)
  new_mm1 = box.map_manager()
  new_mm2 = box.apply_to_map(map_manager = mam.mm.deep_copy())
  assert approx_equal(new_mm1.map_data(), new_mm2.map_data())

  new_model1 = box.model()
  new_model2 = box.apply_to_model(model = mam.model.deep_copy())
  assert new_model1.crystal_symmetry().is_similar_symmetry(
         new_model2.crystal_symmetry())
  assert new_model1.crystal_symmetry().is_similar_symmetry(
         box.crystal_symmetry)

  assert approx_equal(new_model1.get_sites_cart()[0], (19.705233333333336, 15.631525, 13.5040625))
  # make sure things did change
  assert new_mm2.map_data().size() !=  map_data_orig.size()

  # make sure things are changed in-place and are therefore different from start
  assert box.map_manager().map_data().size() !=  map_data_orig.size()
  assert box.model().get_sites_frac() !=  sites_frac_orig
  assert box.model().get_sites_cart() !=   sites_cart_orig
  assert (not cs_orig.is_similar_symmetry(box.model().crystal_symmetry()))

  # make sure box, model and map_manager remember original crystal symmetry
  assert cs_orig.is_similar_symmetry(box.map_manager().unit_cell_crystal_symmetry())
  assert cs_orig.is_similar_symmetry(
    box.map_manager().unit_cell_crystal_symmetry())

  assert approx_equal (box.model().shift_cart(),
     [5.229233333333334, 5.061524999999999, 5.162062499999999])

  assert box.model().unit_cell_crystal_symmetry().is_similar_symmetry(cs_orig)
  assert (not box.model().crystal_symmetry().is_similar_symmetry(cs_orig))

  assert approx_equal(
     box.model()._figure_out_hierarchy_to_output(do_not_shift_back = False
       ).atoms().extract_xyz()[0],
        (14.476, 10.57, 8.342))

  # make sure we can stack shifts
  sel = box.model().selection("resseq 219:219")
  m_small = box.model().select(selection = sel)

  assert approx_equal(box.model().shift_cart(),
     m_small.shift_cart())

  # Now box again:
  small_box = cctbx.maptbx.box.around_model(
    map_manager = mam.mm,
    model       = mam.model.deep_copy(),
    box_cushion     = 5,
    wrapping    = True)

  # Make sure nothing was zeroed out in this map (wrapping = True)
  assert new_mm1.map_data().as_1d().count(0) == 0

  # Now without wrapping...
  box = cctbx.maptbx.box.around_model(
    map_manager = mam.mm,
    model       = mam.model.deep_copy(),
    box_cushion     = 10,
    wrapping    = False)

  # make sure things are changed in-place and are therefore different from start
  assert box.map_manager().map_data().size() !=  map_data_orig.size()
  assert box.model().get_sites_frac() !=  sites_frac_orig
  assert box.model().get_sites_cart() !=   sites_cart_orig
  assert (not cs_orig.is_similar_symmetry(box.model().crystal_symmetry()))

  # make sure box, model and map_manager remember original crystal symmetry
  assert cs_orig.is_similar_symmetry(box.model().unit_cell_crystal_symmetry())
  assert cs_orig.is_similar_symmetry(
    box.map_manager().unit_cell_crystal_symmetry())

  assert box.map_manager().map_data().as_1d().count(0) == 81264

  # Now specify bounds directly
  new_box = cctbx.maptbx.box.with_bounds(
    map_manager = mam.mm.deep_copy(),
    lower_bounds =  (-7, -7, -7),
    upper_bounds =  (37, 47, 39),
    wrapping    = False)

  new_model = new_box.apply_to_model(mam.model.deep_copy())
  # make sure things are changed in-place and are therefore different from start
  assert new_box.map_manager().map_data().size() !=  map_data_orig.size()
  assert new_model.get_sites_frac() !=  sites_frac_orig
  assert new_model.get_sites_cart() !=   sites_cart_orig
  assert (not cs_orig.is_similar_symmetry(new_model.crystal_symmetry()))

  # make sure box, model and map_manager remember original crystal symmetry
  assert cs_orig.is_similar_symmetry(box.model().unit_cell_crystal_symmetry())
  assert cs_orig.is_similar_symmetry(
    box.map_manager().unit_cell_crystal_symmetry())

  assert box.map_manager().map_data().as_1d().count(0) == 81264

  # Now specify bounds directly and init with model
  box = cctbx.maptbx.box.with_bounds(
    map_manager = mam.mm.deep_copy(),
    lower_bounds =  (-7, -7, -7),
    upper_bounds =  (37, 47, 39),
    wrapping    = False,
    model = mam.model.deep_copy())

  new_model = box.model()
  # make sure things are changed in-place and are therefore different from start
  assert box.map_manager().map_data().size() !=  map_data_orig.size()
  assert new_model.get_sites_frac() !=  sites_frac_orig
  assert new_model.get_sites_cart() !=   sites_cart_orig
  assert (not cs_orig.is_similar_symmetry(new_model.crystal_symmetry()))

  # make sure box, model and map_manager remember original crystal symmetry
  assert cs_orig.is_similar_symmetry(box.model().unit_cell_crystal_symmetry())
  assert cs_orig.is_similar_symmetry(
    box.map_manager().unit_cell_crystal_symmetry())

  assert box.map_manager().map_data().as_1d().count(0) == 81264

  # Extract using around_unique

  data_dir = os.path.dirname(os.path.abspath(__file__))
  data_ccp4 = os.path.join(data_dir, 'data', 'D7.ccp4')
  data_ncs =  os.path.join(data_dir, 'data', 'D7.ncs_spec')
  data_seq =  os.path.join(data_dir, 'data', 'D7.seq')

  dm = DataManager(['real_map', 'phil', 'ncs_spec', 'sequence'])
  dm.process_real_map_file(data_ccp4)
  mm = dm.get_real_map(data_ccp4)

  dm.process_ncs_spec_file(data_ncs)
  ncs_obj = dm.get_ncs_spec(data_ncs)

  dm.process_sequence_file(data_seq)
  sequence = dm.get_sequence(data_seq)
  sequence_as_text = sequence[0].sequence

  map_model_mgr=map_model_manager(map_manager=mm,ncs_object=ncs_obj)
  mm=map_model_mgr.map_manager()
  mm.show_summary()

  box = cctbx.maptbx.box.around_unique(
    map_manager = mm.deep_copy(),
    resolution = 3,
    box_cushion = 1,
    sequence = sequence_as_text,
    soft_mask = True,
    wrapping    = False,
   )

  box.map_manager().write_map('new_box.ccp4')

  # run again from map_manager

  map_model_mgr.box_all_maps_around_unique_and_shift_origin(
    resolution = 3,
    box_cushion= 1,
    sequence = sequence_as_text,
    soft_mask = True,
   )

  # Get bounds around density
  box = cctbx.maptbx.box.around_density(
    map_manager = mam.mm.deep_copy(),
    wrapping = False)

  # Create a mask

  mm = mam.mm.deep_copy()

  mm.create_mask_around_density(
        resolution = 3,
        molecular_mass = 2100,
        sequence = "GAVAGA",
        solvent_content = 0.5,
        )
  mask_mm = mm.get_mask_as_map_manager()
  assert approx_equal(
     (mask_mm.map_data().count(0), mask_mm.map_data().count(1),
      mask_mm.map_data().size()),
     (19184, 19216, 38400))

  # Box around the mask
  box = cctbx.maptbx.box.around_mask(
    map_manager = mam.mm.deep_copy(),
    mask_as_map_manager = mask_mm,
    wrapping = False,
    )

  assert (box.gridding_first, box.gridding_last) == ([0, 0, 0] , [29, 39, 31])

  # Box around the mask with cubic box
  box = cctbx.maptbx.box.around_mask(
    map_manager = mam.mm.deep_copy(),
    mask_as_map_manager = mask_mm,
    use_cubic_boxing = True,
    wrapping = False,
    )

  assert (box.gridding_first, box.gridding_last) == ([1, 6, 2], [30, 35, 31])

  #
  # IF you are about to change this - THINK TWICE!
  #
  from libtbx.introspection import getfullargspec
  r = getfullargspec(cctbx.maptbx.box.around_model.__init__)
  assert sorted(r.args)  == sorted(
        ['self', 'map_manager', 'model', 'box_cushion',
        'wrapping', 'model_can_be_outside_bounds', 'stay_inside_current_map',
        'use_cubic_boxing', 'require_match_unit_cell_crystal_symmetry',
        'log']), r.args
  r = getfullargspec(cctbx.maptbx.box.with_bounds.__init__)
  assert sorted(r.args)  ==  sorted(
       ['self', 'map_manager', 'lower_bounds', 'upper_bounds',
        'model', 'wrapping', 'model_can_be_outside_bounds',
        'stay_inside_current_map', 'use_cubic_boxing',
        'require_match_unit_cell_crystal_symmetry', 'log']), r.args

  print ("OK")

if (__name__  ==  "__main__"):
  if libtbx.env.has_module('phenix'):
    exercise_around_model()
  else:
    print('phenix is not available, skipping test')
