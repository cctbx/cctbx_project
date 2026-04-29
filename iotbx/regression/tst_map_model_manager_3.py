from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
import os, sys
from libtbx.utils import Sorry
from libtbx.test_utils import approx_equal
from mmtbx.model import manager as model_manager
from iotbx.data_manager import DataManager

def exercise(file_name, out = sys.stdout):

  # Set up source data
  if not os.path.isfile(file_name):
    raise Sorry("Missing the file: %s" %(file_name)+"\n")

  print ("Reading from %s" %(file_name))
  from iotbx.map_manager import map_manager

  m = map_manager(file_name)

  # make a little model
  sites_cart = flex.vec3_double( ((8, 10, 12), (14, 15, 16)))
  model = model_manager.from_sites_cart(
         atom_name = ' CA ',
         resname = 'ALA',
         chain_id = 'A',
         b_iso = 30.,
         occ = 1.,
         scatterer = 'C',
         sites_cart = sites_cart,
         crystal_symmetry = m.crystal_symmetry())

  # make a map_model_manager with lots of maps and model and ncs
  from iotbx.map_model_manager import map_model_manager

  from mmtbx.ncs.ncs import ncs
  ncs_object=ncs()
  ncs_object.set_unit_ncs()
  mask_mm=m.deep_copy()
  mask_mm.set_is_mask(True)
  mam = map_model_manager(
          map_manager =  m,
          ncs_object =  ncs_object,
          map_manager_1 =  m.deep_copy(),
          map_manager_2 =  m.deep_copy(),
          extra_map_manager_list =  [m.deep_copy(),m.deep_copy(),m.deep_copy()],
          extra_map_manager_id_list = ["extra_1","extra_2","map_manager_mask"],
          model     = model.deep_copy(),)
  print (mam.map_manager())
  print (mam.model())
  print (mam.map_manager_1())
  print (mam.map_manager_2())
  print (mam.map_manager_mask())
  print (mam.map_manager().ncs_object())
  all_map_names=mam.map_id_list()
  for id in all_map_names:
    print("Map_manager %s: %s " %(id,mam.get_map_manager_by_id(id)))

  # Check the method all_objects_have_same_symmetry_and_shift_cart()
  # It should apply to changes in crystal_symmetry, unit_cell_crystal_symmetry,
  #   and shift_cart in model, map_model_manager, map_manager, and ncs_object

  print("\nChecking crystal_symmetry, unit_cell_crystal_symmetry, shift_cart"+
     " consistency")

  from cctbx import uctbx
  from cctbx import crystal
  from iotbx.map_model_manager import all_objects_have_same_symmetry_and_shift_cart

  # Change gridding on a map
  new_mam = mam.deep_copy()
  assert new_mam.check_consistency(stop_on_errors=False,print_errors=True)
  new_mm = new_mam.map_manager().resample_on_different_grid(n_real=(15,16,17))
  new_mam._map_dict['test'] = new_mm
  assert not new_mam.check_consistency(stop_on_errors=False,print_errors=True)

  # Check consistencies of cs, unit_cell cs and shift_cart

  for delta in [0.00001,0.15]:
    # Check unit_cell dims in model crystal_symmetry (this also checks for
    #   map_model_manager as its model is used to return crystal_symmetry).
    new_mam = mam.deep_copy()
    object_list = [new_mam, new_mam.model(), new_mam.map_manager_1(),
       new_mam.map_manager().ncs_object()]
    assert all_objects_have_same_symmetry_and_shift_cart(object_list)
    assert new_mam.check_consistency(stop_on_errors=False,print_errors=True)
    cs = new_mam.crystal_symmetry()
    shift_cart = new_mam.shift_cart()
    # Now change some objects outside the map_model_manager
    params = list(cs.unit_cell().parameters()[:6])
    params[1] += delta
    new_uc = uctbx.unit_cell(params)
    new_cs = crystal.symmetry(unit_cell=new_uc, space_group=cs.space_group())
    print("\nNew CS: ",new_cs)
    new_mam.model().set_crystal_symmetry(new_cs)
    print("\nCHECK WITH DELTA",delta)
    if delta < 0.01:
      assert all_objects_have_same_symmetry_and_shift_cart(object_list)
      assert new_mam.check_consistency(stop_on_errors=False,print_errors=True)
    else:
      assert not all_objects_have_same_symmetry_and_shift_cart(object_list,
        print_errors=True)
      assert not new_mam.check_consistency(stop_on_errors=False,print_errors=True)

  for delta in [None]:
    # Check SG in map_manager
    new_mam = mam.deep_copy()
    object_list = [new_mam, new_mam.model(), new_mam.map_manager_1(),
       new_mam.map_manager().ncs_object()]
    assert all_objects_have_same_symmetry_and_shift_cart(object_list)
    new_mam.map_manager().set_crystal_symmetry_of_partial_map(
      space_group_number=6)
    print("\nCHECK WITH new SG")
    assert not all_objects_have_same_symmetry_and_shift_cart(object_list,
        print_errors=True)

  for delta in [0.00001,0.15]:
    # Check unit_cell in model unit_cell_crystal_symmetry angles
    new_mam = mam.deep_copy()
    object_list = [new_mam, new_mam.model(), new_mam.map_manager_1(),
       new_mam.map_manager().ncs_object()]
    cs = new_mam.unit_cell_crystal_symmetry()
    shift_cart = new_mam.shift_cart()
    assert all_objects_have_same_symmetry_and_shift_cart(object_list)
    params = list(cs.unit_cell().parameters()[:6])
    params[4] += delta
    new_uc = uctbx.unit_cell(params)
    new_cs = crystal.symmetry(unit_cell=new_uc, space_group=cs.space_group())
    print("\nNew CS: ",new_cs)
    new_mam.model().set_unit_cell_crystal_symmetry(new_cs)
    print("\nCHECK WITH DELTA",delta)
    if delta < 0.01:
      assert all_objects_have_same_symmetry_and_shift_cart(object_list)
    else:
      assert not all_objects_have_same_symmetry_and_shift_cart(object_list,
        print_errors=True)

  for delta in [0.00001,0.15]:
    # Check unit_cell angles in model crystal_symmetry
    new_mam = mam.deep_copy()
    object_list = [new_mam, new_mam.model(), new_mam.map_manager_1(),
       new_mam.map_manager().ncs_object()]
    cs = new_mam.crystal_symmetry()
    shift_cart = new_mam.shift_cart()
    assert all_objects_have_same_symmetry_and_shift_cart(object_list)
    params = list(cs.unit_cell().parameters()[:6])
    params[4] += delta
    new_uc = uctbx.unit_cell(params)
    new_cs = crystal.symmetry(unit_cell=new_uc, space_group=cs.space_group())
    print("\nNew CS: ",new_cs)
    new_mam.model().set_crystal_symmetry(new_cs)
    print("\nCHECK WITH DELTA",delta)
    if delta < 0.01:
      assert all_objects_have_same_symmetry_and_shift_cart(object_list)
    else:
      assert not all_objects_have_same_symmetry_and_shift_cart(object_list,
        print_errors=True)

  for delta in [0.00001,0.15]:
    # Check shift_cart
    new_mam = mam.deep_copy()
    object_list = [new_mam, new_mam.model(), new_mam.map_manager_1(),
       new_mam.map_manager().ncs_object()]
    cs = new_mam.crystal_symmetry()
    shift_cart = new_mam.shift_cart()
    assert all_objects_have_same_symmetry_and_shift_cart(object_list)
    params = list(shift_cart)
    params[1] += delta
    new_shift_cart = params
    print("\nNew SC: ",new_shift_cart)
    new_mam.model().set_shift_cart(new_shift_cart)
    print("\nCHECK SC WITH DELTA",delta)
    if delta < 0.01:
      assert all_objects_have_same_symmetry_and_shift_cart(object_list)
    else:
      assert not all_objects_have_same_symmetry_and_shift_cart(object_list,
        print_errors=True)

  for delta in [0.00001,0.15]:
    # Check shift_cart in ncs object
    new_mam = mam.deep_copy()
    object_list = [new_mam, new_mam.model(), new_mam.map_manager_1(),
       new_mam.map_manager().ncs_object()]
    cs = new_mam.crystal_symmetry()
    shift_cart = new_mam.shift_cart()
    assert all_objects_have_same_symmetry_and_shift_cart(object_list)
    params = list(shift_cart)
    params[1] += delta
    new_shift_cart = params
    print("\nNew SC: ",new_shift_cart)
    new_mam.map_manager().ncs_object().set_shift_cart(new_shift_cart)
    print("\nCHECK SC ncs WITH DELTA",delta)
    if delta < 0.01:
      assert all_objects_have_same_symmetry_and_shift_cart(object_list)
      assert new_mam.map_manager().check_consistency(print_errors=True,
         stop_on_errors=True)
    else:
      assert not all_objects_have_same_symmetry_and_shift_cart(object_list,
        print_errors=True)
      assert not new_mam.map_manager().check_consistency(print_errors=True,
         stop_on_errors=False)


  dm = DataManager(['model','miller_array', 'real_map', 'phil','ncs_spec'])
  dm.set_overwrite(True)

  # Create a model with ncs
  from iotbx.regression.ncs.tst_ncs import pdb_str_5
  file_name='tst_mam.pdb'
  f=open(file_name,'w')
  print (pdb_str_5, file = f)
  f.close()

  # Generate map data from this model (it has ncs)
  mmm=map_model_manager()
  mmm.generate_map(box_cushion=0, file_name=file_name,n_residues=500, d_min=3)
  ncs_mam=mmm.deep_copy()
  ncs_mam_copy=mmm.deep_copy()

  # Make sure this model has 126 sites (42 sites times 3-fold ncs)
  assert ncs_mam.model().get_sites_cart().size() == 126
  assert approx_equal (ncs_mam.model().get_sites_cart()[0],
    (23.560999999999996, 8.159, 10.660000000000002), 1e-4)

  # Get just unique part (42 sites)
  unique_mam=ncs_mam.extract_all_maps_around_model(select_unique_by_ncs=True)
  assert unique_mam.model().get_sites_cart().size() == 42
  assert approx_equal (unique_mam.model().get_sites_cart()[0],
    (18.740916666666664, 13.1794, 16.10544), 1e-3)

  # Make sure that the extraction did not change the original but does change
  #   the extracted part
  assert (unique_mam.model().get_sites_cart()[0] !=
     ncs_mam.model().get_sites_cart()[0])  # it was a deep copy so original stays

  # Shift back the extracted part and make sure it matches the original now
  shifted_back_unique_model=mmm.get_model_from_other(unique_mam.deep_copy())
  assert approx_equal (shifted_back_unique_model.get_sites_cart()[0],
    (23.560999999999996, 8.158999999999997, 10.66), 1e-4)

  # Change the extracted model
  sites_cart=unique_mam.model().get_sites_cart()
  sites_cart[0]=(1,1,1)
  unique_mam.model().get_hierarchy().atoms().set_xyz(sites_cart)
  # Note; setting xyz in hierarchy does not set xrs by itself. do that now:
  unique_mam.model().set_sites_cart_from_hierarchy(multiply_ncs=False)

  # Make sure we really changed it
  assert approx_equal (unique_mam.model().get_sites_cart()[0], (1,1,1))

  # Now propagate all the changes in this unique part to entire original model
  #   using NCS
  ncs_mam.propagate_model_from_other(other = unique_mam,
    model_id = 'model',
    other_model_id = 'model')
  # ...and check that copy 1 and copy 2 both change
  assert approx_equal (ncs_mam.model().get_sites_cart()[0],
     (5.820083333333333, -4.020400000000001, -4.445440000000001), 1e-3)
  assert approx_equal (ncs_mam.model().get_sites_cart()[42],
     (38.41904613024224, 17.233251085893276, 2.5547442135142524), 1e-3)

  # Find ncs from map or model
  nn=ncs_mam_copy
  nn.write_map('ncs.ccp4')
  nn.write_model('ncs.pdb')
  ncs_object=nn.get_ncs_from_model()
  dm.write_ncs_spec_file(ncs_object,'ncs.ncs_spec')
  print ("NCS from map",ncs_object)
  nn.set_ncs_object(ncs_object)
  print ("NCS now: ",nn.ncs_object())
  nn.get_ncs_from_map(ncs_object=ncs_object)
  print ("ncs cc:",nn.ncs_cc())
  assert approx_equal(nn.ncs_cc(),0.961915979834,eps=0.01)

  # Make a deep_copy
  dc=mam.deep_copy()
  new_mam=mam.deep_copy()
  assert mam.map_manager().map_data()[0]==new_mam.map_manager().map_data()[0]

  # Make a customized_copy
  new_mam=mam.customized_copy(model_dict={'model':mam.model()})
  assert new_mam.model() is mam.model()
  assert not new_mam.map_dict() is mam.map_dict()

  new_mam=mam.customized_copy(model_dict={'model':mam.model()},
    map_dict=mam.map_dict())
  assert new_mam.model() is mam.model()
  assert new_mam.map_dict() is mam.map_dict()
  print (mam)

  # Add a map
  mam = dc.deep_copy()
  print (mam.map_id_list())
  assert len(mam.map_id_list()) == 6
  mam.add_map_manager_by_id(mam.map_manager().deep_copy(),'new_map_manager')
  print (mam.map_id_list())
  assert len(mam.map_id_list()) == 7

  # duplicate a map
  mam = dc.deep_copy()
  print (mam.map_id_list())
  assert len(mam.map_id_list()) == 6
  mam.duplicate_map_manager('map_manager','new_map_manager')
  print (mam.map_id_list())
  assert len(mam.map_id_list()) == 7

  # resolution_filter a map
  mam = dc.deep_copy()
  print (mam.map_id_list())
  mam.duplicate_map_manager('map_manager','new_map_manager')
  mam.resolution_filter(map_id='new_map_manager',d_min=3.5,d_max=6)

  # Add a model
  mam = dc.deep_copy()
  print (mam.model_id_list())
  assert len(mam.model_id_list()) == 1
  mam.add_model_by_id(mam.model().deep_copy(),'new_model')
  print (mam.model_id_list())
  assert len(mam.model_id_list()) == 2

  # Initialize a map
  mam1=new_mam.deep_copy()
  mam1.initialize_maps(map_value=6)
  assert mam1.map_manager().map_data()[225] == 6

  # Create mask around density and apply to all maps
  mam1=new_mam.deep_copy()
  mam1.mask_all_maps_around_density(solvent_content=0.5,
    soft_mask=False,)
  s = (mam1.get_map_manager_by_id('mask').map_data() > 0.5)
  assert approx_equal( (s.count(True),s.size()), (1024,2048))

  # Create soft mask around density and apply to all maps
  mam1=new_mam.deep_copy()
  mam1.mask_all_maps_around_density(solvent_content=0.5,
    soft_mask=True,)
  s = (mam1.get_map_manager_by_id('mask').map_data() > 0.5)

  # Create mask around edges and apply to all maps
  mam1=new_mam.deep_copy()
  mam1.write_map('before.ccp4')
  mam1.mask_all_maps_around_edges(soft_mask_radius=8)
  mam1.write_map('after.ccp4')
  mam1.write_map(map_id = 'mask',file_name='mask.ccp4')
  s = (mam1.get_map_manager_by_id('mask').map_data() > 0.5)
  assert approx_equal( (s.count(True),s.size()), (1496, 2048))

  # Create a  mask around atoms and apply to all maps
  new_mam.mask_all_maps_around_atoms(mask_atoms_atom_radius=8,
      soft_mask=False)
  s = (new_mam.get_map_manager_by_id('mask').map_data() > 0.5)
  assert approx_equal( (s.count(True),s.size()), (138,2048))

  # Create a soft mask around atoms and apply to all maps
  new_mam.mask_all_maps_around_atoms(mask_atoms_atom_radius=8,
      soft_mask=True)
  s = (new_mam.get_map_manager_by_id('mask').map_data() > 0.5)
  assert approx_equal( (s.count(True),s.size()), (1924,2048))

  # Create a soft mask around atoms and do not do anything with it
  new_mam.create_mask_around_atoms(mask_atoms_atom_radius=8,
      soft_mask=True)
  s = (new_mam.get_map_manager_by_id('mask').map_data() > 0.5)
  assert approx_equal( (s.count(True),s.size()), (1924,2048))

  # Create a soft mask around atoms; do not do anything with it, wrapping =true
  dummy_mam=new_mam.deep_copy()
  dummy_mam.map_manager().set_wrapping(True)
  dummy_mam.create_mask_around_atoms(mask_atoms_atom_radius=8,
      soft_mask=True)
  s = (dummy_mam.get_map_manager_by_id('mask').map_data() > 0.5)
  assert approx_equal( (s.count(True),s.size()), (1924,2048))

  # Mask around edges and do not do anything with it
  mam=dc.deep_copy()
  mam.create_mask_around_edges()
  s = (mam.get_map_manager_by_id('mask').map_data() > 0.5)
  mam.write_map(map_id='mask',file_name='edge.ccp4')
  assert approx_equal( (s.count(True),s.size()), (1792,2048))

  # Mask around density and to not do anything with it
  mam=dc.deep_copy()
  mam.create_mask_around_density(soft_mask=False)
  s = (mam.get_map_manager_by_id('mask').map_data() > 0.5)
  assert approx_equal( (s.count(True),s.size()), (856,2048))

  # Apply the current mask to one map
  mam.apply_mask_to_map('map_manager')
  s = (mam.map_manager().map_data() > 0.)
  assert approx_equal( (s.count(True),s.size()), (424,2048))
  s = (mam.map_manager().map_data() != 0.)
  assert approx_equal( (s.count(True),s.size()), (856,2048))
  assert approx_equal ((mam.map_manager().map_data()[225]),-0.0418027862906)

  # Apply any mask to one map
  mam.apply_mask_to_map('map_manager',mask_id='mask')
  s = (mam.map_manager().map_data() > 0.)
  assert approx_equal( (s.count(True),s.size()), (424,2048))
  s = (mam.map_manager().map_data() != 0.)
  assert approx_equal( (s.count(True),s.size()), (856,2048))
  assert approx_equal ((mam.map_manager().map_data()[225]),-0.0418027862906)

  # Apply the mask to all maps
  mam.apply_mask_to_maps()
  s = (mam.map_manager().map_data() > 0.)
  assert approx_equal( (s.count(True),s.size()), (424,2048))
  s = (mam.map_manager().map_data() != 0.)
  assert approx_equal( (s.count(True),s.size()), (856,2048))
  assert approx_equal ((mam.map_manager().map_data()[225]),-0.0418027862906)

  # Apply the mask to all maps, setting outside value to mean inside
  mam.apply_mask_to_maps(set_outside_to_mean_inside=True)
  s = (mam.map_manager().map_data() > 0.)
  assert approx_equal( (s.count(True),s.size()), (424,2048))
  s = (mam.map_manager().map_data() != 0.)
  assert approx_equal( (s.count(True),s.size()), (2048,2048))
  assert approx_equal ((mam.map_manager().map_data()[2047]),-0.0759598612785)
  s = (mam.get_map_manager_by_id('mask').map_data() >  0).as_1d()
  inside = mam.map_manager().map_data().as_1d().select(s)
  outside = mam.map_manager().map_data().as_1d().select(~s)
  assert approx_equal ((inside.min_max_mean().max,outside.min_max_mean().max),
   (0.317014873028,-0.0159585822888))


  # Make a new map and model, get mam and box with selection
  mmm=map_model_manager()
  mmm.generate_map(box_cushion=0,wrapping=True, d_min=3)
  mam=mmm
  mam_dc=mam.deep_copy()
  new_mm_1=mam.map_manager()
  assert approx_equal( (mmm.map_data().all(),new_mm_1.map_data().all()),
     ((18, 25, 20),(18, 25, 20)))

  # Get local fsc or randomized map
  dc=mam_dc.deep_copy()
  dc.map_manager().set_wrapping(False)
  map_coeffs = dc.map_manager().map_as_fourier_coefficients(d_min=3)
  from cctbx.development.create_models_or_maps import generate_map
  new_mm_1 = generate_map(map_coeffs=map_coeffs,
    d_min=3,
    low_resolution_real_space_noise_fraction=1,
    high_resolution_real_space_noise_fraction=50,
    map_manager=dc.map_manager(),
    random_seed=124321)
  new_mm_2 = generate_map(map_coeffs=map_coeffs,
    d_min=3,
    low_resolution_real_space_noise_fraction=1,
    high_resolution_real_space_noise_fraction=50,
    map_manager=dc.map_manager(),
    random_seed=734119)
  dc.add_map_manager_by_id(new_mm_1,'map_manager_1')
  dc.add_map_manager_by_id(new_mm_2,'map_manager_2')
  cc=dc.map_map_cc()
  fsc_curve=dc.map_map_fsc()
  dc.set_log(sys.stdout)
  dc.local_fsc(n_boxes = 1)

  # Get map-map FSC
  dc=mam_dc.deep_copy()
  dc.duplicate_map_manager(map_id='map_manager',new_map_id='filtered')
  dc.resolution_filter(d_min=3.5, d_max=10, map_id='filtered')
  dc.create_mask_around_atoms()
  fsc_curve=dc.map_map_fsc(
      map_id_1='map_manager',map_id_2='filtered',mask_id='mask',
      resolution=3.5,fsc_cutoff = 0.97)
  assert approx_equal(fsc_curve.d_min, 3.44, eps=0.01)
  assert approx_equal (fsc_curve.fsc.fsc[-1],0.7097, eps = 0.01)

  # Get map-map CC
  dc=mam_dc.deep_copy()
  dc.duplicate_map_manager(map_id='map_manager',new_map_id='filtered')
  dc.resolution_filter(d_min=3.5, d_max=6, map_id='filtered')
  cc=dc.map_map_cc('map_manager','filtered')
  assert approx_equal(cc , 0.6504435255003295)

  # Get map-map CC with mask
  dc=mam_dc.deep_copy()
  dc.duplicate_map_manager(map_id='map_manager',new_map_id='filtered')
  dc.create_mask_around_density(mask_id='filtered')
  cc=dc.map_map_cc('map_manager','filtered',mask_id='mask')
  assert approx_equal(cc , 0.4515628372038732)

  # box around model
  mam=mam_dc.deep_copy()
  mam.box_all_maps_around_model_and_shift_origin(
      selection_string="resseq 221:221")
  new_mm_1=mam.map_manager()
  assert approx_equal( (mam_dc.map_data().all(),new_mm_1.map_data().all()),
    ((18, 25, 20),(24, 20, 20)))

  # box around model and add soft mask to edges
  mam=mam_dc.deep_copy()
  mam.box_all_maps_around_model_and_shift_origin(
      selection_string="resseq 221:221",
      soft_mask_around_edges = True)
  new_mm_2=mam.map_manager()
  assert approx_equal( (mam_dc.map_data().all(),new_mm_2.map_data().all()),
    ((18, 25, 20),(40,35,38)))

  # extract_around_model (get new mam)
  new_mam_dc=mam_dc.extract_all_maps_around_model(
      selection_string="resseq 221:221")
  new_mm_1a=new_mam_dc.map_manager()
  assert approx_equal( (mam_dc.map_data().all(),new_mm_1a.map_data().all()),
    ((18, 25, 20),(24, 20, 20)))
  assert approx_equal(new_mm_1.map_data(),new_mm_1a.map_data())

  # extract_around_model (get new mam) and soft_mask_around_edges
  new_mam_dc=mam_dc.extract_all_maps_around_model(
      selection_string="resseq 221:221", soft_mask_around_edges = True)
  new_mm_2a=new_mam_dc.map_manager()
  assert approx_equal( (mam_dc.map_data().all(),new_mm_2a.map_data().all()),
    ((18, 25, 20),(40,35,38)))
  assert approx_equal(new_mm_2.map_data(),new_mm_2a.map_data())

  # box around_density
  mam2=mam_dc.deep_copy()
  mam2.box_all_maps_around_density_and_shift_origin(box_cushion=0)
  new_mm_2=mam2.map_manager()
  assert approx_equal( (mam_dc.map_data().all(),new_mm_2.map_data().all()),
    ((18, 25, 20),(16, 18, 18)))

  # extract_around_density (get new mam)
  mam2=mam_dc.deep_copy()
  mam2_b=mam2.extract_all_maps_around_density(box_cushion=0)
  new_mm_2=mam2_b.map_manager()
  assert approx_equal( (mam_dc.map_data().all(),new_mm_2.map_data().all()),
    ((18, 25, 20),(16, 18, 18)))

  # Repeat as map_model_manager:
  mmm=mam_dc.as_map_model_manager().deep_copy()
  mmm.box_all_maps_around_model_and_shift_origin(
      selection_string="resseq 221:221")
  new_mm_1a=mmm.map_manager()
  assert approx_equal( (mam_dc.map_data().all(),new_mm_1a.map_data().all()),
    ((18, 25, 20),(24, 20, 20)))
  assert approx_equal(new_mm_1.map_data(),new_mm_1a.map_data())

  # box around density
  mam = mam_dc.deep_copy()
  mam.box_all_maps_around_density_and_shift_origin(box_cushion=0,soft_mask_around_edges=False)
  new_mm_1=mam.map_manager()
  assert approx_equal( (mam_dc.map_data().all(),new_mm_1.map_data().all()),
    ((18,25 , 20),(16, 18, 18)))

  # box around density and soft mask edges
  mam = mam_dc.deep_copy()
  mam.box_all_maps_around_density_and_shift_origin(box_cushion=0,
   soft_mask_around_edges=True)
  new_mm_1=mam.map_manager()
  assert approx_equal( (mam_dc.map_data().all(),new_mm_1.map_data().all()),
    ((18, 25 , 20),(18, 25, 20)))

  # extract around density (get new mam)
  mam1=mam_dc.deep_copy()
  mam1.extract_all_maps_around_density(box_cushion=0)
  new_mm_1=mam1.map_manager()
  assert approx_equal( (mmm.map_data().all(),new_mm_1.map_data().all()),
    ((24, 20, 20),(18, 25, 20)))

  # create mask around density, then box around mask (i.e., box around density)
  mam.create_mask_around_density(soft_mask=False)
  mam.box_all_maps_around_mask_and_shift_origin(box_cushion=3)
  new_mm_1=mam.map_manager()
  assert approx_equal( (mmm.map_data().all(),new_mm_1.map_data().all()),
    ((24, 20, 20),(18, 25, 20)))

  # box with bounds
  mam.box_all_maps_with_bounds_and_shift_origin(lower_bounds=(10,10,10),
     upper_bounds=(15,15,15))
  new_mm_1=mam.map_manager()
  assert approx_equal( (mmm.map_data().all(),new_mm_1.map_data().all()),
    ((24, 20, 20),(6, 6, 6)))

  # extract with bounds
  mam=mam_dc.deep_copy()
  mam_1=mam.extract_all_maps_with_bounds(lower_bounds=(10,10,10),
     upper_bounds=(15,15,15))
  new_mm_1=mam_1.map_manager()
  assert approx_equal( (mmm.map_data().all(),new_mm_1.map_data().all()),
    ((24, 20, 20),(6, 6, 6)))

  # box with unique
  mam=mam_dc.deep_copy()
  mam.box_all_maps_around_unique_and_shift_origin(
      molecular_mass=2500,resolution=3)
  new_mm_1=mam.map_manager()
  assert approx_equal( (mmm.map_data().all(),new_mm_1.map_data().all()),
    ((24, 20, 20),(18, 25, 20)))

  # extract with unique
  mam=mam_dc.deep_copy()
  mam_1=mam.extract_all_maps_around_unique(
      molecular_mass=2500,resolution=3)
  new_mm_1=mam_1.map_manager()
  assert approx_equal( (mmm.map_data().all(),new_mm_1.map_data().all()),
    ((24,20, 20),(18, 25, 20)))

  # extract a box and then restore model into same reference as current mam
  mam=mam_dc.deep_copy()
  mam.box_all_maps_with_bounds_and_shift_origin(lower_bounds=(2,2,2),
     upper_bounds=(17,17,17))
  print("mam:",mam.model().get_sites_cart()[0],mam.map_manager().origin_is_zero())
  # extract a box
  box_mam=mam.extract_all_maps_with_bounds(lower_bounds=(10,10,10),
     upper_bounds=(15,15,15))
  box_model=box_mam.model()
  matched_box_model=mam.get_model_from_other(box_mam)
  assert approx_equal(matched_box_model.get_sites_cart()[0],mam.model().get_sites_cart()[0])

  # Convert a map to fourier coefficients
  mam=mam_dc.deep_copy()
  ma=mam.map_as_fourier_coefficients(d_min=3)
  assert approx_equal(ma.d_min(),3.01655042414)


  mam.add_map_from_fourier_coefficients(ma, map_id='new_map_manager')
  cc=flex.linear_correlation(
   mam.get_map_manager_by_id('map_manager').map_data().as_1d(),
   mam.get_map_manager_by_id('new_map_manager').map_data().as_1d()).coefficient()
  assert (cc >= 0.99)

  # Get map-model CC
  dc=mam_dc.extract_all_maps_around_model(
      selection_string="(name ca or name cb or name c or name o) "+
        "and resseq 221:221", box_cushion=0)
  cc=dc.map_model_cc(resolution=3)
  assert approx_equal (cc, 0.817089390421)

  # Get map-model density
  dc=mam_dc.extract_all_maps_around_model(
      selection_string="(name ca or name cb or name c or name o) "+
        "and resseq 221:221", box_cushion=0)
  density=dc.density_at_model_sites(selection_string = 'name ca')
  assert approx_equal (density.min_max_mean().mean, 0.841152333991)


  # Remove model outside map
  dc.remove_model_outside_map(boundary=0)
  assert (mam_dc.model().get_sites_cart().size(),
     dc.model().get_sites_cart().size()) == (86, 4)

  # shift a model to match the map
  dc=mam_dc.extract_all_maps_around_model(
      selection_string="(name ca or name cb or name c or name o) "+
        "and resseq 221:221", box_cushion=0)
  actual_model=dc.model().deep_copy()
  working_model=dc.model().deep_copy()
  working_model.set_shift_cart((0,0,0))
  working_model.set_sites_cart(working_model.get_sites_cart()-actual_model.shift_cart())
  dc.shift_any_model_to_match(working_model)
  assert approx_equal (actual_model.get_sites_cart()[0],working_model.get_sites_cart()[0])


if __name__ == "__main__":
  args = sys.argv[1:]
  import libtbx.load_env
  if not args:
    file_name = libtbx.env.under_dist(
      module_name = "iotbx",
      path = "ccp4_map/tst_input.map")
    args = [file_name]
  if libtbx.env.has_module("phenix"):
    exercise(file_name = args[0])
  else:
    print("Skipped: Requires phenix module")

  print ("OK")
