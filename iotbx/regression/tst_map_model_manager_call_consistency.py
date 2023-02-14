from __future__ import absolute_import, division, print_function
from iotbx.map_manager import map_manager
from iotbx.map_model_manager import map_model_manager
from cctbx.maptbx import box
from cctbx.maptbx.mask import create_mask_around_edges, \
      create_mask_around_atoms, create_mask_around_density
import cctbx
from libtbx.introspection import getfullargspec

def get_method_text(key, base_method_name):
  if key == 'init':
    return "box.%s.__init__" %(base_method_name)
  elif key == 'box_all_maps':
    return "box_all_maps_%s_and_shift_origin" %(base_method_name)
  elif key == 'extract_all_maps':
    return "extract_all_maps_%s" %(base_method_name)
  if key == 'init_mask':
    return "cctbx.maptbx.mask.create_mask_%s.__init__" %(base_method_name)
  if key == 'create_mask':
    return "map_model_manager.create_mask_%s" %(base_method_name)
  else:
    raise AssertionError("unknown key for get_method_text")

def get_method(key, base_method_name):
  if key == 'init':
    return getattr(getattr(box,base_method_name),"__init__")
  elif key == 'box_all_maps':
    method_name = "box_all_maps_%s_and_shift_origin" %(base_method_name)
    return  getattr(map_model_manager, method_name)
  elif key == 'extract_all_maps':
    method_name = "extract_all_maps_%s" %(base_method_name)
    return  getattr(map_model_manager, method_name)
  if key == 'init_mask':
    if base_method_name == 'around_edges':
      return create_mask_around_edges.__init__
    elif base_method_name == 'around_atoms':
      return create_mask_around_atoms.__init__
    elif base_method_name == 'around_density':
      return create_mask_around_density.__init__
  elif key == 'create_mask':
    method_name = "create_mask_%s" %(base_method_name)
    return  getattr(map_model_manager, method_name)
  else:
    raise AssertionError("unknown key for get_method")

def all_expected_in_found(found = None, expected = None):
  for x in expected:
    if not x in found:
      print("Expected to find",x," but did not")
      return False
  return True

def check_args(text, method,expected_args,group_text,
       allow_extra_in_found=None):
    found_args = getfullargspec(method).args
    expected_args = sorted(expected_args)
    found_args = sorted(found_args)
    print ("\n%s :\nExpected :%s \nFound   : %s" %(
       text,str(expected_args),str(found_args)))
    if expected_args == found_args:
      return
    elif allow_extra_in_found and all_expected_in_found(
      found = found_args, expected = expected_args):
      return
    else:  # give message
      error_message="""

FAIL : Args for %s do not match expected.
Expected :%s
Found    :%s

 If you change args for
    %s
Make sure that all of these are changed to match and also change expected
to match in iotbx/regression/tst_map_model_manager_call_consistency.py
""" %(
        text,str(expected_args),str(found_args),group_text)
      raise AssertionError (error_message)

def test_01():

  #  Making sure that all the convenience calls to each basic method for boxing
  #  have the same parameters as the basic method
  #  Base calls are in cctbx.maptbx.box and cctbx.maptbx.mask.
  #  Matching calls are in map_manager, map_model_manager

  # Args that should appear in all calls to all methods
  common_args = ['self',]

  # Args that should be in common for all calls to specific methods
  method_args_dict = {
    'with_bounds': ['lower_bounds', 'upper_bounds',
    'model_can_be_outside_bounds',
    'stay_inside_current_map', 'use_cubic_boxing','require_match_unit_cell_crystal_symmetry'],
    'around_model':[ 'box_cushion','model_can_be_outside_bounds',
      'stay_inside_current_map', 'use_cubic_boxing',
      'require_match_unit_cell_crystal_symmetry'],
    'around_density':[ 'box_cushion','threshold', 'get_half_height_width',
       'stay_inside_current_map', 'use_cubic_boxing',
       'model_can_be_outside_bounds',
      'require_match_unit_cell_crystal_symmetry',],
    'around_mask':[ 'box_cushion','model_can_be_outside_bounds',
      'stay_inside_current_map', 'use_cubic_boxing',
      'require_match_unit_cell_crystal_symmetry'],
    'around_unique':['box_cushion', 'target_ncs_au_model',
        'stay_inside_current_map', 'use_cubic_boxing',
        'use_symmetry_in_extract_unique', 'regions_to_keep',
         'require_match_unit_cell_crystal_symmetry',
    'residues_per_region','keep_this_region_only',
        'solvent_content', 'resolution', 'sequence', 'molecular_mass',
         'symmetry', 'chain_type', 'keep_low_density', 'soft_mask',
         'mask_expand_ratio'],
   }

  # Args that should appear in map_manager and map_model_manager calls for
  #   specific methods
  manager_method_args_dict = {
    'with_bounds': [],
    'around_model':['selection', 'selection_string', 'select_unique_by_ncs'],
    'around_density':['map_id'],
    'around_mask':['mask_id'],
    'around_unique':[],
   }

  # Args that should appear in init calls for
  #   specific methods
  init_method_args_dict = {
    'with_bounds': [],
    'around_model':[],
    'around_density':[],
    'around_mask':['mask_as_map_manager'],
    'around_unique':[],
   }

  # Args that should appear in calls that are __init__
  init_args = ['map_manager', 'model', 'wrapping', 'log']

  # Args that should appear in calls that are "extract_xxx"
  extract_args = []

  # Args that should appear in calls that are "box_all_maps_xxx"
  box_all_maps_args = ['extract_box']

  base_method_name_list = method_args_dict.keys()
  assert method_args_dict.keys() == manager_method_args_dict.keys()


  for base_method_name in base_method_name_list:
    print ("\nExpected and actual args in calls for %s" %(base_method_name))
    group_text ="""
    cctbx.mmtbx.box.%s
    iotbx.map_model_manager.extract_all_maps_%s
    iotbx.map_model_manager.box_all_maps_%s_and_shift_origin""" %(
      base_method_name,base_method_name,base_method_name)

    # Check call in cctbx.maptbx.box.xxx.__init__
    method = get_method('init',base_method_name)
    text = get_method_text('init',base_method_name)
    expected_args = common_args + method_args_dict[base_method_name] + \
      init_method_args_dict[base_method_name] + init_args
    check_args(text,method,expected_args,group_text)


    # Check call in iotbx.map_model_manager box_all_maps_xxx_and_shift_origin
    method = get_method('box_all_maps',base_method_name)
    text = get_method_text('box_all_maps',base_method_name)

    expected_args = common_args + method_args_dict[base_method_name] + \
       manager_method_args_dict[base_method_name] + box_all_maps_args
    check_args(text,method,expected_args,group_text,
       allow_extra_in_found = True)

    # Check call in iotbx.map_model_manager extract_all_maps_xxx
    method = get_method('extract_all_maps',base_method_name)
    text = get_method_text('extract_all_maps',base_method_name)

    expected_args = common_args + method_args_dict[base_method_name] + \
       manager_method_args_dict[base_method_name] + extract_args
    check_args(text,method,expected_args,group_text,
       allow_extra_in_found = True)

def test_02():
  # Make sure calls in map_model_manager to
  # create_mask_around_edges
  # create_mask_around_atoms
  # create_mask_around_density
  #  have args matching the main methods in cctbx.maptbx.mask

  # Just put what is expected here and ask programmer to change in both
  #  places if they change one

  init_arg_dict = {
     'around_atoms':['self', 'mask_atoms_atom_radius', 'model',
          'invert_mask','xray_structure', 'map_manager','n_real', 'wrapping'],
     'around_edges':['self', 'boundary_radius', 'map_manager'],
     'around_density':['self', 'map_manager','resolution','molecular_mass',
        'sequence','solvent_content'],
  }
  method_arg_dict = {
     'around_atoms':['self', 'soft_mask_radius', 'mask_atoms_atom_radius',
          'invert_mask','soft_mask', 'mask_id', 'model'],
     'around_edges':['self', 'boundary_radius', 'mask_id',],
     'around_density':['self', 'resolution', 'solvent_content', 'soft_mask',
           'soft_mask_radius', 'mask_id','map_id'],
  }
  assert init_arg_dict.keys() == method_arg_dict.keys()

  #  Check them out...
  for base_method_name in init_arg_dict.keys():
    print ("\nExpected and actual args in calls for %s" %(base_method_name))
    group_text ="""
      cctbx.mmtbx.mask.create_mask_around_%s
      iotbx.map_model_manager.create_mask_around_%s""" %(
        base_method_name,base_method_name)

    method = get_method('init_mask',base_method_name)
    text = get_method_text('init_mask',base_method_name)
    expected_args = init_arg_dict[base_method_name]
    check_args(text,method,expected_args,group_text)

    method = get_method('create_mask',base_method_name)
    text = get_method_text('create_mask',base_method_name)
    expected_args = method_arg_dict[base_method_name]
    check_args(text,method,expected_args,group_text,allow_extra_in_found=True)


if (__name__  ==  '__main__'):
  test_01()
  test_02()
  print ("OK")
