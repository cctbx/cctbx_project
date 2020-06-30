from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
import os, sys
from libtbx.utils import Sorry
from libtbx.test_utils import approx_equal
from mmtbx.model import manager as model_manager
from iotbx.data_manager import DataManager

def exercise(file_name, out = sys.stdout):
  if not os.path.isfile(file_name):
    raise Sorry("Missing the file: %s" %(file_name)+"\n"+"Please run with phenix.python read_write_mrc.py my_mrcfile.mrc")

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
  mam = map_model_manager(
          map_manager =  m,
          ncs_object =  ncs_object,
          map_manager_1 =  m.deep_copy(),
          map_manager_2 =  m.deep_copy(),
          extra_map_manager_list =  [m.deep_copy(),m.deep_copy(),m.deep_copy()],
          extra_map_id_list = ["extra_1","extra_2","map_manager_mask"],
          model     = model.deep_copy(),)
  r_model=mam.as_r_model()

  print (r_model.map_manager())
  print (r_model.model())
  print (r_model.map_manager_1())
  print (r_model.map_manager_2())
  print (r_model.map_manager_mask())
  print (r_model.map_manager().ncs_object())
  all_map_names=r_model.map_id_list()
  for id in all_map_names:
    print("Map_manager %s: %s " %(id,r_model.get_map_manager_by_id(id)))

  dm = DataManager(['model','miller_array', 'real_map', 'phil'])
  dm.set_overwrite(True)


  # Make a deep_copy
  dc=r_model.deep_copy()
  new_r_model=r_model.deep_copy()
  assert r_model.map_manager().map_data()[0]==new_r_model.map_manager().map_data()[0]

  # Make a customized_copy
  new_r_model=r_model.customized_copy(model_dict={'model':r_model.model()})
  assert new_r_model.model() is r_model.model()
  assert not new_r_model.map_dict() is r_model.map_dict()

  new_r_model=r_model.customized_copy(model_dict={'model':r_model.model()},
    map_dict=r_model.map_dict())
  assert new_r_model.model() is r_model.model()
  assert new_r_model.map_dict() is r_model.map_dict()
  print (r_model)

  # Add a map
  rm = dc.deep_copy()
  print (rm.map_id_list())
  assert len(rm.map_id_list()) == 6
  rm.add_map_manager_by_id(rm.map_manager().deep_copy(),'new_map_manager')
  print (rm.map_id_list())
  assert len(rm.map_id_list()) == 7

  # duplicate a map
  rm = dc.deep_copy()
  print (rm.map_id_list())
  assert len(rm.map_id_list()) == 6
  rm.duplicate_map_manager('map_manager','new_map_manager')
  print (rm.map_id_list())
  assert len(rm.map_id_list()) == 7

  # resolution_filter a map
  rm = dc.deep_copy()
  print (rm.map_id_list())
  rm.duplicate_map_manager('map_manager','new_map_manager')
  rm.resolution_filter(map_id='new_map_manager',high_resolution=3.5,low_resolution=6)

  # Add a model
  rm = dc.deep_copy()
  print (rm.model_id_list())
  assert len(rm.model_id_list()) == 1
  rm.add_model_by_id(rm.model().deep_copy(),'new_model')
  print (rm.model_id_list())
  assert len(rm.model_id_list()) == 2

  # Initialize a map
  rm1=new_r_model.deep_copy()
  rm1.initialize_maps(map_value=6)
  assert rm1.map_manager().map_data()[225] == 6

  # Create mask around density and apply to all maps
  rm1=new_r_model.deep_copy()
  rm1.mask_all_maps_around_density(solvent_content=0.5,
    soft_mask=True,)
  s = (rm1.get_map_manager_by_id('mask').map_data() > 0.5)
  assert approx_equal( (s.count(True),s.size()), (1024,2048))

  # Create mask around edges and apply to all maps
  rm1=new_r_model.deep_copy()
  rm1.mask_all_maps_around_edges()
  s = (rm1.get_map_manager_by_id('mask').map_data() > 0.5)
  assert approx_equal( (s.count(True),s.size()), (1176,2048))

  # Create a soft mask around model and apply to all maps
  new_r_model.mask_all_maps_around_atoms(mask_atoms_atom_radius=8,
      soft_mask=True)
  s = (new_r_model.get_map_manager_by_id('mask').map_data() > 0.5)
  assert approx_equal( (s.count(True),s.size()), (339,2048))

  # Create a soft mask around model and do not do anything with it
  new_r_model.create_mask_around_atoms(mask_atoms_atom_radius=8,
      soft_mask=True)
  s = (new_r_model.get_map_manager_by_id('mask').map_data() > 0.5)
  assert approx_equal( (s.count(True),s.size()), (339,2048))

  # Create a sharp mask around model and do not do anything with it
  new_r_model.create_mask_around_atoms(soft_mask=False,
     mask_atoms_atom_radius=8)
  s = (new_r_model.get_map_manager_by_id('mask').map_data() > 0.5)
  assert approx_equal( (s.count(True),s.size()), (35,2048))

  # Mask around edges and do not do anything with it
  r_model=dc.deep_copy()
  r_model.create_mask_around_edges()
  s = (r_model.get_map_manager_by_id('mask').map_data() > 0.5)
  assert approx_equal( (s.count(True),s.size()), (1176,2048))

  # Mask around density and to not do anything with it
  r_model=dc.deep_copy()
  r_model.create_mask_around_density(soft_mask=False)
  s = (r_model.get_map_manager_by_id('mask').map_data() > 0.5)
  assert approx_equal( (s.count(True),s.size()), (856,2048))

  # Apply the current mask to one map
  r_model.apply_mask_to_map('map_manager')
  s = (r_model.map_manager().map_data() > 0.)
  assert approx_equal( (s.count(True),s.size()), (424,2048))
  s = (r_model.map_manager().map_data() != 0.)
  assert approx_equal( (s.count(True),s.size()), (856,2048))
  assert approx_equal ((r_model.map_manager().map_data()[225]),-0.0418027862906)

  # Apply any mask to one map
  r_model.apply_mask_to_map('map_manager',mask_id='mask')
  s = (r_model.map_manager().map_data() > 0.)
  assert approx_equal( (s.count(True),s.size()), (424,2048))
  s = (r_model.map_manager().map_data() != 0.)
  assert approx_equal( (s.count(True),s.size()), (856,2048))
  assert approx_equal ((r_model.map_manager().map_data()[225]),-0.0418027862906)

  # Apply the mask to all maps
  r_model.apply_mask_to_maps()
  s = (r_model.map_manager().map_data() > 0.)
  assert approx_equal( (s.count(True),s.size()), (424,2048))
  s = (r_model.map_manager().map_data() != 0.)
  assert approx_equal( (s.count(True),s.size()), (856,2048))
  assert approx_equal ((r_model.map_manager().map_data()[225]),-0.0418027862906)

  # Apply the mask to all maps, setting outside value to mean inside
  r_model.apply_mask_to_maps(set_outside_to_mean_inside=True)
  s = (r_model.map_manager().map_data() > 0.)
  assert approx_equal( (s.count(True),s.size()), (424,2048))
  s = (r_model.map_manager().map_data() != 0.)
  assert approx_equal( (s.count(True),s.size()), (2048,2048))
  assert approx_equal ((r_model.map_manager().map_data()[2047]),-0.0759598612785)
  s = (r_model.get_map_manager_by_id('mask').map_data() >  0).as_1d()
  inside = r_model.map_manager().map_data().as_1d().select(s)
  outside = r_model.map_manager().map_data().as_1d().select(~s)
  assert approx_equal ((inside.min_max_mean().max,outside.min_max_mean().max),
   (0.317014873028,-0.0159585822888))


  # Make a new map and model, get r_model and box with selection
  mmm=map_model_manager()
  mmm.generate_map(box_cushion=0)
  rm=mmm.as_r_model()
  rm_dc=rm.deep_copy()

  new_mm_1=rm.map_manager()
  assert approx_equal( (mmm.map_data().all(),new_mm_1.map_data().all()),
     ((18, 25, 20),(18, 25, 20)))


  # Get map-map FSC
  dc=rm_dc.deep_copy()
  dc.duplicate_map_manager(map_id='map_manager',new_map_id='filtered')
  dc.resolution_filter(high_resolution=3.5, low_resolution=6, map_id='filtered')
  dc.create_mask_around_atoms()
  fsc_curve=dc.map_map_fsc(
      map_id_1='map_manager',map_id_2='filtered',mask_id='mask',
      resolution=3,fsc_cutoff = 0.97)
  assert approx_equal(fsc_curve.d_min, 3.43091149578)
  assert approx_equal (fsc_curve.fsc.fsc[-1],0.942973996865)

  # Get map-map CC
  dc=rm_dc.deep_copy()
  dc.duplicate_map_manager(map_id='map_manager',new_map_id='filtered')
  dc.resolution_filter(high_resolution=3.5, low_resolution=6, map_id='filtered')
  cc=dc.map_map_cc('map_manager','filtered')
  assert approx_equal(cc , 0.955283969999)

  # Get map-map CC with mask
  dc=rm_dc.deep_copy()
  dc.duplicate_map_manager(map_id='map_manager',new_map_id='filtered')
  dc.resolution_filter(high_resolution=3.5, low_resolution=6, map_id='filtered')
  dc.create_mask_around_density(mask_id='filtered')
  cc=dc.map_map_cc('map_manager','filtered',mask_id='mask')
  assert approx_equal(cc , 0.629060115596)

  # box around model
  rm=rm_dc.deep_copy()
  rm.box_all_maps_around_model_and_shift_origin(
      selection_string="resseq 221:221")
  new_mm_1=rm.map_manager()
  assert approx_equal( (mmm.map_data().all(),new_mm_1.map_data().all()),
    ((18, 25, 20),(24, 20, 20)))

  # extract_around_model (get new r_model)
  new_rm_dc=rm_dc.extract_all_maps_around_model(
      selection_string="resseq 221:221")
  new_mm_1a=new_rm_dc.map_manager()
  assert approx_equal( (mmm.map_data().all(),new_mm_1a.map_data().all()),
    ((18, 25, 20),(24, 20, 20)))
  assert approx_equal(new_mm_1.map_data(),new_mm_1a.map_data())

  # box around_density
  rm2=rm_dc.deep_copy()
  rm2.box_all_maps_around_density_and_shift_origin(box_cushion=0)
  new_mm_2=rm2.map_manager()
  assert approx_equal( (mmm.map_data().all(),new_mm_2.map_data().all()),
    ((18, 25, 20),(16, 23, 18)))

  # extract_around_density (get new r_model)
  rm2=rm_dc.deep_copy()
  rm2_b=rm2.extract_all_maps_around_density(box_cushion=0)
  new_mm_2=rm2_b.map_manager()
  assert approx_equal( (mmm.map_data().all(),new_mm_2.map_data().all()),
    ((18, 25, 20),(16, 23, 18)))

  # Repeat as map_model_manager:
  mmm=rm_dc.as_map_model_manager().deep_copy()
  mmm.box_all_maps_around_model_and_shift_origin(
      selection_string="resseq 221:221")
  new_mm_1a=mmm.map_manager()
  assert approx_equal( (mmm.map_data().all(),new_mm_1a.map_data().all()),
    ((24, 20, 20),(24, 20, 20)))
  assert approx_equal(new_mm_1.map_data(),new_mm_1a.map_data())

  # box around density
  rm.box_all_maps_around_density_and_shift_origin(box_cushion=0)
  new_mm_1=rm.map_manager()
  assert approx_equal( (mmm.map_data().all(),new_mm_1.map_data().all()),
    ((24, 20 , 20),(22, 18, 18)))

  # extract around density (get new r_model)
  rm1=rm_dc.deep_copy()
  rm1.extract_all_maps_around_density(box_cushion=0)
  new_mm_1=rm1.map_manager()
  assert approx_equal( (mmm.map_data().all(),new_mm_1.map_data().all()),
    ((24, 20 , 20),(18, 25, 20)))

  # create mask around density, then box around mask (i.e., box around density)
  rm.create_mask_around_density(soft_mask=False)
  rm.box_all_maps_around_mask_and_shift_origin(box_cushion=3)
  new_mm_1=rm.map_manager()
  assert approx_equal( (mmm.map_data().all(),new_mm_1.map_data().all()),
    ((24, 20 , 20),(22, 18, 18)))

  # box with bounds
  rm.box_all_maps_with_bounds_and_shift_origin(lower_bounds=(10,10,10),
     upper_bounds=(15,15,15))
  new_mm_1=rm.map_manager()
  assert approx_equal( (mmm.map_data().all(),new_mm_1.map_data().all()),
    ((24, 20, 20),(6, 6, 6)))

  # extract with bounds
  rm=rm_dc.deep_copy()
  rm_1=rm.extract_all_maps_with_bounds(lower_bounds=(10,10,10),
     upper_bounds=(15,15,15))
  new_mm_1=rm_1.map_manager()
  assert approx_equal( (mmm.map_data().all(),new_mm_1.map_data().all()),
    ((24, 20, 20),(6, 6, 6)))

  # box with unique
  rm=rm_dc.deep_copy()
  rm.box_all_maps_around_unique_and_shift_origin(
      molecular_mass=2500,resolution=3)
  new_mm_1=rm.map_manager()
  assert approx_equal( (mmm.map_data().all(),new_mm_1.map_data().all()),
    ((24, 20, 20),(18, 25, 20)))

  # extract with unique
  rm=rm_dc.deep_copy()
  rm_1=rm.extract_all_maps_around_unique(
      molecular_mass=2500,resolution=3)
  new_mm_1=rm_1.map_manager()
  assert approx_equal( (mmm.map_data().all(),new_mm_1.map_data().all()),
    ((24, 20, 20),(18, 25, 20)))

  # extract a box and then restore model into same reference as current r_model
  rm=rm_dc.deep_copy()
  rm.box_all_maps_with_bounds_and_shift_origin(lower_bounds=(2,2,2),
     upper_bounds=(17,17,17))
  print("rm:",rm.model().get_sites_cart()[0],rm.map_manager().origin_is_zero())
  # extract a box
  box_rm=rm.extract_all_maps_with_bounds(lower_bounds=(10,10,10),
     upper_bounds=(15,15,15))
  box_model=box_rm.model()
  matched_box_model=rm.get_model_from_other(box_rm)
  assert approx_equal(matched_box_model.get_sites_cart()[0],rm.model().get_sites_cart()[0])

  # Convert a map to fourier coefficients
  rm=rm_dc.deep_copy()
  ma=rm.map_as_fourier_coefficients(high_resolution=3)
  assert approx_equal(ma.d_min(),3.01655042414)


  rm.add_map_from_fourier_coefficients(ma, map_id='new_map_manager')
  cc=flex.linear_correlation(
   rm.get_map_manager_by_id('map_manager').map_data().as_1d(),
   rm.get_map_manager_by_id('new_map_manager').map_data().as_1d()).coefficient()
  assert (cc >= 0.99)

  # Get map-model CC
  dc=rm_dc.extract_all_maps_around_model(
      selection_string="(name ca or name cb or name c or name o) "+
        "and resseq 221:221", box_cushion=0)
  cc=dc.map_model_cc(resolution=3)
  assert approx_equal (cc, 0.413802839326)

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
