from __future__ import absolute_import, division, print_function
import iotbx.mrcfile
from cctbx.array_family import flex
import os, sys
from libtbx.test_utils import approx_equal
from libtbx.utils import Sorry
from mmtbx.model import manager as model_manager
from iotbx.data_manager import DataManager
from iotbx import mrcfile

def exercise(file_name=None, pdb_file_name = None, map_file_name = None ,
    split_pdb_file_name = None,
    ncs_pdb_file_name = None,
    out = sys.stdout):

  # Set up source data

  if not os.path.isfile(file_name):
    raise Sorry("Missing the file: %s" %(file_name)+"\n")

  print ("Reading from %s" %(file_name))
  from iotbx.map_manager import map_manager
  m = map_manager(file_name)

  print ("Header information from %s:" %(file_name))
  m.show_summary(out = out)

  map_data = m.map_data().deep_copy()
  crystal_symmetry = m.crystal_symmetry()
  unit_cell_parameters = m.crystal_symmetry().unit_cell().parameters()

  print ("\nMap origin: %s Extent %s"  %( map_data.origin(), map_data.all()))
  print ("Original unit cell, not just unit cell of part in this file): %s" %(
     str(unit_cell_parameters)))

  grid_point = (1, 2, 3)
  if map_data.origin() !=  (0, 0, 0): # make sure it is inside
    from scitbx.matrix import col
    grid_point = tuple (col(grid_point)+col(map_data.origin()))
  print ("\nValue of map_data at grid point %s: %.3f" %(str(grid_point),
    map_data[grid_point]))
  print ("Map data is %s" %(type(map_data)))

  random_position = (10, 5, 7.9)
  point_frac = crystal_symmetry.unit_cell().fractionalize(random_position)
  value_at_point_frac = map_data.eight_point_interpolation(point_frac)
  print ("Value of map_data at coordinates %s: %.3f" %(
      str(random_position), value_at_point_frac))

  map_data_as_float = map_data.as_float()
  print ("Map data as float is %s" %(type(map_data_as_float)))


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
         crystal_symmetry = crystal_symmetry)


  # Move map and a model to place origin at (0, 0, 0)
  # map data is new copy but model is shifted in place.

  from iotbx.map_model_manager import map_model_manager
  mam = map_model_manager(
          map_manager =  m,
          model     = model.deep_copy(),
    )

  # Read in map and model and split up
  dm = DataManager()
  aa = dm.get_map_model_manager(model_file=pdb_file_name,
    map_files=map_file_name)
  cc = dm.get_map_model_manager(model_file=ncs_pdb_file_name,
    map_files=map_file_name)
  bb = dm.get_map_model_manager(model_file=split_pdb_file_name,
    map_files=map_file_name)

  # Merge by models
  a = aa.deep_copy()
  n_starting_models = len(list(a.model().get_hierarchy().models()))
  box_info = a.split_up_map_and_model_by_boxes()
  # Change the hierarchy in a box
  small_hierarchy = box_info.mmm_list[0].model().get_hierarchy()
  # delete an atom
  for m in small_hierarchy.models():
    for chain in m.chains()[:1]:
       chain.remove_residue_group(i = 0)
  box_info.mmm_list[0].model().reset_after_changing_hierarchy()  # REQUIRED
  # Put everything back together
  a.merge_split_maps_and_models(box_info = box_info,
      allow_changes_in_hierarchy = True)
  n_merged_models = len(list(a.model().get_hierarchy().models()))
  assert n_starting_models == 1
  assert n_merged_models == 7


  # Merge in various ways
  for selection_method in ['by_ncs_groups', 'by_chain', 'by_segment',
       'supplied_selections', 'boxes']:
    if selection_method == 'boxes':
      choices = [True, False]
    else:
      choices = [True]
    if selection_method == 'by_chain':
      mask_choices = [True,False]
    else:
      mask_choices = [False]
    for select_final_boxes_based_on_model in choices:
      for skip_empty_boxes in choices:
        for mask_choice in mask_choices:
          if selection_method == 'by_ncs_groups':
            a=cc.deep_copy()
          if mask_choice: # use split model
            a=bb.deep_copy()
          else: # usual
            a=aa.deep_copy()
          print ("\nRunning split_up_map_and_model with \n"+
            "select_final_boxes_based_on_model="+
           "%s   skip_empty_boxes=%s selection_method=%s" %(
            select_final_boxes_based_on_model,skip_empty_boxes,selection_method))

          if selection_method == 'by_chain':
            print ("Mask around unused atoms: %s" %(mask_choice))
            box_info = a.split_up_map_and_model_by_chain(
              mask_around_unselected_atoms=mask_choice)
          elif selection_method == 'by_segment':
            box_info = a.split_up_map_and_model_by_segment()
          elif selection_method == 'supplied_selections':
            selection = a.model().selection('all')
            box_info = a.split_up_map_and_model_by_supplied_selections(
              selection_list = [selection])
          elif selection_method == 'by_ncs_groups':
            box_info = a.split_up_map_and_model_by_ncs_groups()
          elif selection_method == 'boxes':
            box_info = a.split_up_map_and_model_by_boxes(
              skip_empty_boxes = skip_empty_boxes,
              select_final_boxes_based_on_model =
                select_final_boxes_based_on_model)
          assert box_info is not None
          print (selection_method,skip_empty_boxes,
              len(box_info.selection_list),
              box_info.selection_list[0].count(True))
          assert (selection_method,skip_empty_boxes,
              len(box_info.selection_list),
              box_info.selection_list[0].count(True)) in [
                ('by_chain',True,3,19),
                ("by_chain",True,1,86,),
                ("by_segment",True,1,86,),
                ("supplied_selections",True,1,86,),
                ("by_ncs_groups",True,1,86),
                ("boxes",True,7,9),
                ("boxes",False,12,0,),
                ("boxes",True,13,1,),
                ("boxes",False,36,0,),
                ], 'failed to find %s %s %s %s' % (
                    selection_method,
                    skip_empty_boxes,
                    len(box_info.selection_list),
                    box_info.selection_list[0].count(True))

          # Change the coordinates in one box
          small_model = box_info.mmm_list[0].model()
          small_sites_cart = small_model.get_sites_cart()
          from scitbx.matrix import col
          small_sites_cart += col((1,0,0))
          small_model.set_crystal_symmetry_and_sites_cart(
            sites_cart = small_sites_cart,
            crystal_symmetry = small_model.crystal_symmetry())
          # Put everything back together
          a.merge_split_maps_and_models(box_info = box_info)


  mam.box_all_maps_around_model_and_shift_origin()

  shifted_crystal_symmetry = mam.model().crystal_symmetry()
  shifted_model = mam.model()
  shifted_map_data = mam.map_data()

  print ("\nOriginal map origin (grid units):", map_data.origin())
  print ("Original model:\n", model.model_as_pdb())

  print ("Shifted map origin:", shifted_map_data.origin())
  print ("Shifted model:\n", shifted_model.model_as_pdb())


  # Save the map_model manager
  mam_dc=mam.deep_copy()
  print ("dc",mam)
  print ("dc mam_dc",mam_dc)

  # Mask map around atoms
  mam=mam_dc.deep_copy()
  print ("dc mam_dc dc",mam_dc)
  print (mam)
  mam.mask_all_maps_around_atoms(mask_atoms_atom_radius = 3,
     set_outside_to_mean_inside=True, soft_mask=False)
  print ("Mean before masking", mam.map_data().as_1d().min_max_mean().mean)
  assert approx_equal(mam.map_data().as_1d().min_max_mean().mean,
      -0.0585683621466)
  print ("Max before masking", mam.map_data().as_1d().min_max_mean().max)
  assert approx_equal(mam.map_data().as_1d().min_max_mean().max,
      -0.0585683621466)

  # Mask map around atoms, with soft mask
  mam=mam_dc.deep_copy()
  mam.mask_all_maps_around_atoms(mask_atoms_atom_radius = 3, soft_mask = True,
    soft_mask_radius = 5, set_outside_to_mean_inside=True)
  print ("Mean after first masking", mam.map_data().as_1d().min_max_mean().mean)
  assert approx_equal(mam.map_data().as_1d().min_max_mean().mean,
      0, eps= 0.05)
  print ("Max after first masking", mam.map_data().as_1d().min_max_mean().max)
  assert approx_equal(mam.map_data().as_1d().min_max_mean().max,
       0.10, eps = 0.05)

  # Mask map around atoms again
  mam.mask_all_maps_around_atoms(mask_atoms_atom_radius = 3,
     set_outside_to_mean_inside = True, soft_mask=False)
  print ("Mean after second masking", mam.map_data().as_1d().min_max_mean().mean)
  assert approx_equal(mam.map_data().as_1d().min_max_mean().mean,
     0, eps=0.1)
  print ("Max after second masking", mam.map_data().as_1d().min_max_mean().max)
  assert approx_equal(mam.map_data().as_1d().min_max_mean().max,
      0, eps=0.1)

  # Mask around edges
  mam=mam_dc.deep_copy()
  mam.mask_all_maps_around_edges( soft_mask_radius = 3)
  print ("Mean after masking edges", mam.map_data().as_1d().min_max_mean().mean)
  assert approx_equal(mam.map_data().as_1d().min_max_mean().mean,
      0, eps=0.05)
  print ("Max after masking edges", mam.map_data().as_1d().min_max_mean().max)
  assert approx_equal(mam.map_data().as_1d().min_max_mean().max,
      0.20, eps= 0.05)


  print ("\nWriting map_data and model in shifted position (origin at 0, 0, 0)")

  output_file_name = 'shifted_map.ccp4'
  print ("Writing to %s" %(output_file_name))
  mrcfile.write_ccp4_map(
      file_name = output_file_name,
      crystal_symmetry = shifted_crystal_symmetry,
      map_data = shifted_map_data, )

  output_file_name = 'shifted_model.pdb'
  f = open(output_file_name, 'w')
  print (shifted_model.model_as_pdb(), file=f)
  f.close()


  print ("\nWriting map_data and model in original position (origin at %s)" %(
      str(mam.map_manager().origin_shift_grid_units)))

  output_file_name = 'new_map_original_position.ccp4'
  print ("Writing to %s" %(output_file_name))
  mrcfile.write_ccp4_map(
      file_name = output_file_name,
      crystal_symmetry = shifted_crystal_symmetry,
      map_data = shifted_map_data,
      origin_shift_grid_units = mam.map_manager().origin_shift_grid_units)
  print (shifted_model.model_as_pdb())
  output_pdb_file_name = 'new_model_original_position.pdb'
  f = open(output_pdb_file_name, 'w')
  print (shifted_model.model_as_pdb(), file=f)
  f.close()

  # Write as mmcif
  output_cif_file_name = 'new_model_original_position.cif'
  f = open(output_cif_file_name, 'w')
  print (shifted_model.model_as_mmcif(),file = f)
  f.close()


  # Read the new map and model
  import iotbx.pdb
  with open(output_pdb_file_name) as f:
    lines = f.read()
  new_model =  model_manager(
     model_input = iotbx.pdb.input(
         source_info = None,
         lines = flex.split_lines(lines)),
         crystal_symmetry = crystal_symmetry)
  assert new_model.model_as_pdb() == model.model_as_pdb()

  with open(output_cif_file_name) as f:
    f.read()
  new_model_from_cif =  model_manager(
     model_input = iotbx.pdb.input(
         source_info = None,
         lines = flex.split_lines(lines)),
         crystal_symmetry = crystal_symmetry)
  assert new_model_from_cif.model_as_pdb() == model.model_as_pdb()

  # Read and box the original file again in case we modified m in any
  #   previous tests
  m = map_manager(file_name)
  mam=map_model_manager(model=model.deep_copy(),map_manager=m)
  mam.box_all_maps_around_model_and_shift_origin()

  file_name = output_file_name
  print ("Reading from %s" %(file_name))
  new_map = iotbx.mrcfile.map_reader(file_name = file_name, verbose = False)
  new_map.data = new_map.data.shift_origin()
  print ("Header information from %s:" %(file_name))
  new_map.show_summary(out = out)
  assert new_map.map_data().origin() == mam.map_manager().map_data().origin()
  assert new_map.crystal_symmetry().is_similar_symmetry(mam.map_manager().crystal_symmetry())

  # make a map_model_manager with lots of maps and model and ncs
  from mmtbx.ncs.ncs import ncs
  ncs_object=ncs()
  ncs_object.set_unit_ncs()
  mam = map_model_manager(
          map_manager =  m,
          ncs_object =  ncs_object,
          map_manager_1 =  m.deep_copy(),
          map_manager_2 =  m.deep_copy(),
          extra_model_list =  [model.deep_copy(),model.deep_copy()],
          extra_model_id_list = ["model_1","model_2"],
          extra_map_manager_list =  [m.deep_copy(),m.deep_copy()],
          extra_map_manager_id_list = ["extra_1","extra_2"],
          model     = model.deep_copy(),
    )


  # make a map_model_manager with lots of maps and model and ncs and run
  # with wrapping and ignore_symmetry_conflicts on
  from mmtbx.ncs.ncs import ncs
  ncs_object=ncs()
  ncs_object.set_unit_ncs()
  m.set_ncs_object(ncs_object.deep_copy())
  mam2 = map_model_manager(
          map_manager =  m.deep_copy(),
          ncs_object =  ncs_object.deep_copy(),
          map_manager_1 =  m.deep_copy(),
          map_manager_2 =  m.deep_copy(),
          extra_model_list =  [model.deep_copy(),model.deep_copy()],
          extra_model_id_list = ["model_1","model_2"],
          extra_map_manager_list =  [m.deep_copy(),m.deep_copy()],
          extra_map_manager_id_list = ["extra_1","extra_2"],
          model     = model.deep_copy(),
          ignore_symmetry_conflicts = True,
          wrapping = m.wrapping(),
    )
  assert mam.map_manager().is_similar(mam2.map_manager())
  assert mam.map_manager().is_similar(mam2.map_manager_1())
  for m in mam2.map_managers():
    assert mam.map_manager().is_similar(m)
  assert mam.model().shift_cart() == mam2.model().shift_cart()
  assert mam.model().shift_cart() == mam2.get_model_by_id('model_2').shift_cart()



  print ("OK")

if __name__ == "__main__":
  args = sys.argv[1:]
  if not args or len(args)<3:
    import libtbx.load_env
    file_name = libtbx.env.under_dist(
      module_name = "iotbx",
      path = "ccp4_map/tst_input.map")
    split_pdb_file_name = libtbx.env.under_dist(
      module_name = "iotbx",
      path = "regression/data/non_zero_origin_model_split.pdb")
    pdb_file_name = libtbx.env.under_dist(
      module_name = "iotbx",
      path = "regression/data/non_zero_origin_model.pdb")
    map_file_name = libtbx.env.under_dist(
      module_name = "iotbx",
      path = "regression/data/non_zero_origin_map.ccp4")
    ncs_pdb_file_name = libtbx.env.under_dist(
      module_name = "iotbx",
      path = "regression/data/non_zero_origin_ncs_model.pdb")
    args = [file_name,pdb_file_name, map_file_name,split_pdb_file_name,
        ncs_pdb_file_name]
  if libtbx.env.find_in_repositories(relative_path='chem_data') is not None:
    exercise(file_name = args[0], pdb_file_name=args[1], map_file_name=args[2],
     split_pdb_file_name=args[3],ncs_pdb_file_name=args[4])
  else:
    print('chem_data is not available, skipping')
