from __future__ import absolute_import, division, print_function
import os
import libtbx.load_env
from iotbx.data_manager import DataManager
from iotbx.map_model_manager import _match_map_model_ncs, map_model_manager
from iotbx.phil import parse
from libtbx.program_template import ProgramTemplate
from libtbx.test_utils import approx_equal
from libtbx.utils import Sorry
from libtbx.test_utils import show_diff

def test_01():

  # Source data

  data_dir = os.path.dirname(os.path.abspath(__file__))
  data_ccp4 = os.path.join(data_dir, 'data',
                          'non_zero_origin_map.ccp4')
  data_pdb = os.path.join(data_dir, 'data',
                          'non_zero_origin_model.cif')
  data_ncs_spec = os.path.join(data_dir, 'data',
                          'non_zero_origin_ncs_spec.ncs_spec')

  # DataManager

  dm = DataManager(['ncs_spec','model', 'real_map', 'phil'])
  dm.set_overwrite(True)

  # Read in map and model and ncs

  map_file=data_ccp4
  dm.process_real_map_file(map_file)
  mm = dm.get_real_map(map_file)

  model_file=data_pdb
  dm.process_model_file(model_file)
  model = dm.get_model(model_file)

  ncs_file=data_ncs_spec
  dm.process_ncs_spec_file(ncs_file)
  ncs = dm.get_ncs_spec(ncs_file)

  ncs_dc = ncs.deep_copy()

  mmmn = _match_map_model_ncs()
  mmmn.add_map_manager(mm)
  mmmn.add_model(model)
  mmmn.add_ncs_object(ncs)

  # Save it
  mmmn_dc=mmmn.deep_copy()

  # Test creating mmm from model:
  mmm_from_model = model.as_map_model_manager(create_model_map = False)
  mmm_from_model = model.as_map_model_manager(create_model_map = True,
    resolution = 5)
  assert mmm_from_model.map_manager() is not None

  # Make sure we can add an ncs object that is either shifted or not
  mmmn_dcdc=mmmn.deep_copy()
  new_mmmn = _match_map_model_ncs()
  new_mmmn.add_map_manager(mmmn_dcdc.map_manager())
  new_mmmn.add_model(mmmn_dcdc.model())
  new_mmmn.add_ncs_object(mmmn_dcdc.ncs_object())
  assert new_mmmn.ncs_object().shift_cart() == new_mmmn.map_manager().shift_cart()

  mmmn_dcdc=mmmn.deep_copy()
  new_mmmn = _match_map_model_ncs()
  new_mmmn.add_map_manager(mmmn_dcdc.map_manager())
  new_mmmn.add_model(mmmn_dcdc.model())
  new_mmmn.add_ncs_object(ncs_dc)
  assert new_mmmn.ncs_object().shift_cart() == new_mmmn.map_manager().shift_cart()


  original_ncs=mmmn.ncs_object()
  assert approx_equal((24.0528, 11.5833, 20.0004),
     tuple(original_ncs.ncs_groups()[0].translations_orth()[-1]),
     eps=0.1)

  assert tuple(mmmn._map_manager.origin_shift_grid_units) == (0,0,0)

  # Shift origin to (0,0,0)
  mmmn=mmmn_dc.deep_copy()  # fresh version of _match_map_model_ncs
  mmmn.shift_origin()
  new_ncs=mmmn.ncs_object()
  assert tuple(mmmn._map_manager.origin_shift_grid_units) == (100,100,100)

  mmmn.write_model('s.pdb')
  mmmn.write_map('s.mrc')

  shifted_ncs=mmmn.ncs_object()
  assert approx_equal((-153.758, -74.044, -127.487),
      tuple(shifted_ncs.ncs_groups()[0].translations_orth()[-1]),eps=0.1)


  # Shift a model and shift it back

  mmmn=mmmn_dc.deep_copy()  # fresh version of _match_map_model_ncs
  model=mmmn.model()
  shifted_model=mmmn.shift_model_to_match_working_map(model=model)
  model_in_original_position=mmmn.shift_model_to_match_original_map(
      model=shifted_model)
  assert (approx_equal(model.get_sites_cart(), # not a copy
                      shifted_model.get_sites_cart()))
  assert approx_equal(model.get_sites_cart(),
                      model_in_original_position.get_sites_cart())

  # test data_manager map_model_manager
  generated_mmm = dm.get_map_model_manager()
  print(generated_mmm)
  assert(isinstance(generated_mmm,map_model_manager))


  # Generate a map and model

  import sys
  mmm=map_model_manager(log=sys.stdout)
  mmm.generate_map()
  model=mmm.model()
  mm=mmm.map_manager()
  assert approx_equal(
     model.get_sites_cart()[0], (14.476, 10.57, 8.34) ,eps=0.01)
  assert approx_equal(mm.map_data()[10,10,10],-0.0506,eps=0.001)
  # Save it
  mmm_dc=mmm.deep_copy()

  # Create model from sites
  mmm_sites = mmm_dc.deep_copy()
  from scitbx.array_family import flex
  sites_cart = flex.vec3_double()
  sites_cart.append((3,4,5))
  mmm_sites.model_from_sites_cart(sites_cart = sites_cart,
    model_id = 'new_model')
  assert mmm_sites.get_model_by_id('new_model').get_sites_cart()[0] == (3,4,5)
  ph_sites = mmm_sites.get_model_by_id('new_model').get_hierarchy()
  text_sites = mmm_sites.get_model_by_id('new_model').model_as_pdb()

  # Create model from hierarchy
  mmm_sites = mmm_dc.deep_copy()
  mmm_sites.model_from_hierarchy(hierarchy = ph_sites,
    model_id = 'new_model')
  assert mmm_sites.get_model_by_id('new_model').get_sites_cart()[0] == (3,4,5)

  # Create model from text
  mmm_sites = mmm_dc.deep_copy()
  mmm_sites.model_from_text(text = text_sites,
    model_id = 'new_model')
  assert mmm_sites.get_model_by_id('new_model').get_sites_cart()[0] == (3,4,5)



  # Set crystal_symmetry and unit_cell_crystal_symmetry and shift_cart
  # Box and shift the map_model_manager so we have new coordinate system
  mmm_sites.box_all_maps_around_model_and_shift_origin(box_cushion=4.5)
  new_model = mmm_sites.get_model_by_id('new_model')
  assert approx_equal((3., 4., 5.0),
     mmm_sites.get_model_by_id('new_model').get_sites_cart()[0])

  # arbitrarily set unit_cell crystal symmetry of model to
  #  match crystal_symmetry. First have to set shift_cart to None
  new_model.set_shift_cart(shift_cart = None)
  new_model.set_unit_cell_crystal_symmetry_and_shift_cart()
  assert new_model.crystal_symmetry() != mmm_sites.crystal_symmetry()

  # now set crystal symmetries and shift cart of model to match the manager
  mmm_sites.set_model_symmetries_and_shift_cart_to_match_map(new_model)
  assert new_model.crystal_symmetry().is_similar_symmetry(
    mmm_sites.crystal_symmetry())
  assert new_model.unit_cell_crystal_symmetry().is_similar_symmetry(
    mmm_sites.unit_cell_crystal_symmetry())
  assert new_model.shift_cart() == mmm_sites.shift_cart()

  # Import hierarchy into a model and set symmetries and shift to match
  mmm_sites.model_from_hierarchy(
   hierarchy = mmm_sites.model().get_hierarchy(),
   model_id='model_from_hierarchy')
  assert not show_diff(mmm_sites.get_model_by_id('model_from_hierarchy').model_as_pdb(),
     mmm_sites.get_model_by_id('model').model_as_pdb())


  # Check on wrapping
  assert not mm.wrapping()  # this one should not wrap because it is zero at edges

  # Make a new one with no buffer so it is not zero at edges
  mmm=map_model_manager()
  mmm.generate_map(box_cushion=0)
  mm=mmm.map_manager()
  # check its compatibility with wrapping
  assert mm.is_consistent_with_wrapping()
  mmm.show_summary()

  # now box it
  sel=mmm.model().selection("resseq 221:221")
  new_model=mmm.model().deep_copy().select(sel)
  new_mmm=map_model_manager(model=new_model,map_manager=mm.deep_copy())
  new_mmm.box_all_maps_around_model_and_shift_origin()
  new_mm=new_mmm.map_manager()

  assert not new_mm.wrapping()
  assert not new_mm.is_consistent_with_wrapping()

  # now box it with selection
  new_mmm_1=map_model_manager(
      model=mmm.model().deep_copy(),map_manager=mm.deep_copy())
  new_mmm_1.box_all_maps_around_model_and_shift_origin(
      selection_string="resseq 221:221")
  new_mm_1=new_mmm_1.map_manager()

  assert not new_mm_1.wrapping()
  assert not new_mm_1.is_consistent_with_wrapping()
  assert new_mm_1.map_data().all()== new_mm.map_data().all()

  # create map_model_manager with just half-maps
  mm1=mm.deep_copy()
  mm2=mm.deep_copy()
  map_data=mm2.map_data()
  map_data+=1.
  new_mmm=map_model_manager(model=mmm.model().deep_copy(),
    map_manager_1=mm1,
    map_manager_2=mm2)
  assert new_mmm._map_dict.get('map_manager') is None # should not be any yet
  assert approx_equal(new_mmm.map_manager().map_data()[232],
     mm.deep_copy().map_data()[232]+0.5)
  assert new_mmm._map_dict.get('map_manager') is not None # now should be there

  # generate map data from a model
  mm1=mm.deep_copy()
  mm2=mm.deep_copy()
  new_mmm=map_model_manager(model=mmm.model().deep_copy(), map_manager=mm1)
  mmm.generate_map(model=mmm.model())
  mm=mmm.map_manager()
  mmm.show_summary()

  # check get_map_model_manager function
  dm = DataManager(['model'])
  assert not hasattr(dm, 'get_map_model_manager')
  dm = DataManager(['real_map'])
  assert hasattr(dm, 'get_map_model_manager')
  dm = DataManager(['sequence'])
  assert not hasattr(dm, 'get_map_model_manager')
  dm = DataManager(['model', 'real_map'])
  assert hasattr(dm, 'get_map_model_manager')

  # usage
  dm.get_map_model_manager(model_file=data_pdb, map_files=data_ccp4)
  dm.get_map_model_manager(model_file=data_pdb, map_files=[data_ccp4])
  dm.get_map_model_manager(model_file=data_pdb, map_files=[data_ccp4, data_ccp4, data_ccp4])
  dm.get_map_model_manager(model_file=data_pdb, map_files=data_ccp4, ignore_symmetry_conflicts=True)

  # errors
  try:
    dm.get_map_model_manager(model_file=data_pdb, map_files=data_ccp4, from_phil=True)
  except Sorry as e:
    assert 'from_phil is set to True' in str(e)
  try:
    dm.get_map_model_manager(model_file=data_pdb, map_files=data_ccp4, abc=123)
  except TypeError as e:
    assert 'unexpected keyword argument' in str(e)
  try:
    dm.get_map_model_manager(model_file=data_pdb, map_files=[data_ccp4, data_ccp4])
  except Sorry as e:
    assert '1 full map and 2 half maps' in str(e)

  # PHIL
  class test_program(ProgramTemplate):
    master_phil_str = '''
include scope iotbx.map_model_manager.map_model_phil_str
'''
  working_phil_str = '''
  map_model {
    full_map = %s
    half_map = %s
    half_map = s.mrc
    model = %s
  }
''' % (data_ccp4, data_ccp4, data_pdb)

  master_phil = parse(test_program.master_phil_str, process_includes=True)
  working_phil = master_phil.fetch(parse(working_phil_str))
  tp = test_program(dm, working_phil.extract())

  try:
    dm.get_map_model_manager(from_phil=True)
  except Exception as e:
    assert 'ignore_symmetry_conflicts' in str(e)
  try:
    dm.get_map_model_manager(from_phil=True, ignore_symmetry_conflicts=True)
  except AssertionError:
    pass

# ----------------------------------------------------------------------------

if (__name__ == '__main__'):
  if libtbx.env.find_in_repositories(relative_path='chem_data') is not None:
    test_01()
  else:
    print('Skip test_01, chem_data not available')
  print ("OK")
