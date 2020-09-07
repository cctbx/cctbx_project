from __future__ import absolute_import, division, print_function
import sys
from libtbx.utils import Sorry
from cctbx import maptbx
from libtbx import group_args
from scitbx.array_family import flex
import iotbx.map_manager
from mmtbx.model import manager as model_manager
import mmtbx.ncs.ncs
from libtbx.utils import null_out
from libtbx.test_utils import approx_equal
from copy import deepcopy

# Reserved phil scope for MapModelManager
map_model_phil_str = '''
map_model {
  full_map = None
    .type = path
    .help = Input full map file
    .short_caption = Full map filename
  half_map = None
    .type = path
    .multiple = True
    .help = Input half map files
    .short_caption = Half map filenames
  model = None
    .type = path
    .help = Input model file
    .short_caption = Model filename
}
'''

class map_model_manager(object):

  '''
    Class for shifting origin of map(s) and model to (0, 0, 0) and keeping
    track of the shifts.

    Typical use:
    mam = map_model_manager(
      model = model,
      map_manager = map_manager,
      ncs_object = ncs_object)

    mam.box_all_maps_around_model_and_shift_origin(
        box_cushion=3)

    shifted_model = mam.model()  # at (0, 0, 0), knows about shifts
    shifted_map_manager = mam.map_manager() # also at (0, 0, 0) knows shifts
    shifted_ncs_object = mam.ncs_object() # also at (0, 0, 0) and knows shifts

    Optional after boxing:  apply soft mask to map (requires soft_mask_radius)

    The maps allowed are:
      map_dict has four special ids with interpretations:
        map_manager:  full map
        map_manager_1, map_manager_2: half-maps 1 and 2
        map_manager_mask:  a mask as a map_manager
      All other ids are any strings and are assumed to correspond to other maps

    Note:  It is permissible to call with no map_manger, but supplying
      both map_manager_1 and map_manager_2.  In this case, the working
      map_manager will be the average of map_manager_1 and map_manager_2. This
      will be created the first time map_manager is referenced.

    Note:  mam.map_manager() contains mam.ncs_object(), so it is not necessary
    to keep both.

    Note: model objects may contain internal ncs objects.  These are separate
    from those in the map_managers and separate from ncs_object in the call to
    map_model_manager.

    The ncs_object describes the NCS of the map and is a property of the map. It
    must be shared by all maps

    The model ncs objects describe the NCS of the individual model. They can
    differ between models and between models and maps.

    Note: set wrapping of all maps to match map_manager if they differ. Set
    all to be wrapping if it is set

  '''

  def __init__(self,
               model            = None,
               map_manager      = None,
               map_manager_1    = None,
               map_manager_2    = None,
               extra_model_list = None,
               extra_model_id_list = None,  # string id's for models
               extra_map_manager_list = None,
               extra_map_manager_id_list = None,  # string id's for map_managers
               ncs_object       = None,   # Overwrite map_manager ncs_objects
               ignore_symmetry_conflicts = None,  # allow mismatch of symmetry
               wrapping         = None,  # Overwrite wrapping for all maps
               absolute_angle_tolerance = 0.01,  # angle tolerance for symmetry
               absolute_length_tolerance = 0.01,  # length tolerance
               log              = None):

    # Checks
    if extra_model_list is None: extra_model_list = []
    if extra_map_manager_list is None: extra_map_manager_list = []
    for m in [model] + extra_model_list:
      assert (m is None) or isinstance(m, model_manager)
    for mm in [map_manager, map_manager_1, map_manager_2] + extra_map_manager_list:
      assert (mm is None) or isinstance(mm, iotbx.map_manager.map_manager)
    assert (ncs_object is None) or isinstance(ncs_object, mmtbx.ncs.ncs.ncs)


    # Set the log stream
    self.set_log(log = log)

    # Initialize
    self._map_dict={}
    self._model_dict = {}
    self._force_wrapping = wrapping
    self._warning_message = None




    # If no map_manager now, do not do anything and make sure there
    #    was nothing else supplied except possibly a model

    if (not map_manager) and (not map_manager_1) and (not map_manager_2):
      assert not extra_map_manager_list
      assert not ncs_object
      self._model_dict = {'model': model}
      return  # do not do anything

    # Now make sure we have map manager or half maps at least
    assert map_manager or (map_manager_1 and map_manager_2)

    # A map_manager to check against others. It will be either map_manager or
    #   map_manager_1

    if map_manager:
      any_map_manager = map_manager
      any_map_manager_is_map_manager = True
    else:
      any_map_manager = map_manager_1
      any_map_manager_is_map_manager = False

    # Overwrite wrapping if requested
    # Take wrapping from any_map_manager otherwise for all maps

    if isinstance(self._force_wrapping, bool):
      wrapping = self._force_wrapping
      if wrapping and (not any_map_manager.is_full_size()):
        raise Sorry("You cannot use wrapping=True if the map is not full size")
    else:
      wrapping = any_map_manager.wrapping()

    assert wrapping in [True, False]
    if not extra_map_manager_list:
      extra_map_manager_list=[]
    for m in [map_manager, map_manager_1, map_manager_2]+ \
       extra_map_manager_list:
      if m:
        m.set_wrapping(wrapping)

    # if ignore_symmetry_conflicts, take all symmetry information from
    #  any_map_manager and apply it to everything
    if ignore_symmetry_conflicts:
      if ncs_object:
        ncs_object.set_shift_cart(any_map_manager.shift_cart())
      if model:
        any_map_manager.set_model_symmetries_and_shift_cart_to_match_map(model)
      if extra_model_list:
        for m in extra_model_list:
          any_map_manager.set_model_symmetries_and_shift_cart_to_match_map(m)

      if map_manager_1 and (map_manager_1 is not any_map_manager):
        map_manager_1 = any_map_manager.customized_copy(
          map_data=map_manager_1.map_data())
      if map_manager_2:
        map_manager_2 = any_map_manager.customized_copy(
          map_data=map_manager_2.map_data())

      new_extra_map_manager_list = []
      for m in extra_map_manager_list:
        new_extra_map_manager_list.append(any_map_manager.customized_copy(
          map_data=m.map_data()))
      extra_map_manager_list = new_extra_map_manager_list

    # CHECKS

    # Make sure that any_map_manager is either already shifted to (0, 0, 0)
    #  or has  origin_shift_grid_unit of (0, 0, 0).
    assert any_map_manager.origin_is_zero() or \
      tuple(any_map_manager.origin_shift_grid_units) == (0, 0, 0)

    # Normally any_map_manager unit_cell_crystal_symmetry should match
    #  model original_crystal_symmetry (and also usually model.crystal_symmetry)

    # Make sure we have what is expected: optional model, mm or
    # map_manager_1 and map_manager_2 or neither,
    #   optional list of extra_map_manager_list and extra_model_list

    if extra_map_manager_list:
      if extra_map_manager_id_list:
        assert len(extra_map_manager_list) == len(extra_map_manager_id_list)
      else:
        extra_map_manager_id_list=[]
        for i in range(1,len(extra_map_manager_list)+1):
          extra_map_manager_id_list.append("extra_map_manager_%s" %(i))
    else:
      extra_map_manager_list = []
      extra_map_manager_id_list = []

    if extra_model_list:
      if extra_model_id_list:
        assert len(extra_model_list) == len(extra_model_id_list)
      else:
        extra_model_id_list=[]
        for i in range(1,len(extra_model_list)+1):
          extra_model_id_list.append("extra_model_%s" %(i))
    else:
      extra_model_list = []
      extra_model_id_list = []



    if(not [map_manager_1, map_manager_2].count(None) in [0, 2]):
      raise Sorry("None or two half-maps are required.")

    if(not any_map_manager):
      raise Sorry("A map is required.")

    # Make sure all map_managers have same gridding and symmetry
    for m in [map_manager, map_manager_1, map_manager_2]+ \
         extra_map_manager_list:
      if m and (not ignore_symmetry_conflicts):
        if not any_map_manager.is_similar(m,
           absolute_angle_tolerance = absolute_angle_tolerance,
           absolute_length_tolerance = absolute_length_tolerance,
          ):
          raise Sorry("Map manager '%s' is not similar to '%s': %s" %(
           m.file_name,any_map_manager.file_name,
            map_manager.warning_message())+
            "\nTry 'ignore_symmetry_conflicts=True'")

    # Now make sure all models match symmetry using match_map_model_ncs

    # Make a match_map_model_ncs and check unit_cell and
    #   working crystal symmetry
    #  and shift_cart for model, map, and ncs_object (if present)
    mmmn = match_map_model_ncs(
        absolute_angle_tolerance = absolute_angle_tolerance,
        absolute_length_tolerance = absolute_length_tolerance,
        ignore_symmetry_conflicts = ignore_symmetry_conflicts)
    mmmn.add_map_manager(any_map_manager)
    if model:
      mmmn.add_model(model,
        set_model_log_to_null = False,
        ) # keep the log
    if ncs_object:
      mmmn.add_ncs_object(ncs_object) # overwrites anything in any_map_manager

    # All ok here if it did not stop

    # Shift origin of model and any_map_manager and ncs_object to (0, 0, 0) with
    #    mmmn which knows about all of them

    mmmn.shift_origin()

    # any_map_manager, model, ncs_object know about shift

    any_map_manager = mmmn.map_manager()
    # Put shifted map in the right place. It is either map_manager or map_manager_1
    if any_map_manager_is_map_manager:
      map_manager = any_map_manager
    else:
      map_manager_1 = any_map_manager

    if model:
       assert mmmn.model() is not None # make sure we got it
    model = mmmn.model()  # this model knows about shift

    if model:
      # Make sure model shift manager agrees with any_map_manager shift
      assert approx_equal(model.shift_cart(), any_map_manager.shift_cart())

    # Shift origins of all maps (shifting again does nothing, but still skip if
    #    already done (it was done for any_map_manager))

    for m in [map_manager, map_manager_1, map_manager_2]+\
         extra_map_manager_list:
      if m and (not m is any_map_manager):
        m.shift_origin()

    # Shift origins of all the extra models:
    for m in extra_model_list:
      m.shift_model_and_set_crystal_symmetry(
          shift_cart=any_map_manager.shift_cart(),
          crystal_symmetry=map_manager.crystal_symmetry())
      assert approx_equal(m.shift_cart(), any_map_manager.shift_cart())

    # Transfer ncs_object to all map_managers if one is present
    if self.ncs_object():
      for m in [map_manager, map_manager_1, map_manager_2]+\
           extra_map_manager_list:
        if m:
          m.set_ncs_object(self.ncs_object())

    # Make sure all really match:
    for m in [map_manager, map_manager_1, map_manager_2]+\
        extra_map_manager_list:
      if m and not any_map_manager.is_similar(m):
          raise AssertionError(any_map_manager.warning_message())

    # Set up maps, model, as dictionaries (same as used in map_model_manager)
    self.set_up_map_dict(
      map_manager = map_manager,
      map_manager_1 = map_manager_1,
      map_manager_2 = map_manager_2,
      extra_map_manager_list = extra_map_manager_list,
      extra_map_manager_id_list = extra_map_manager_id_list)

    self.set_up_model_dict(
      model = model,
      extra_model_list = extra_model_list,
      extra_model_id_list = extra_model_id_list)

  def set_up_map_dict(self,
      map_manager = None,
      map_manager_1 = None,
      map_manager_2 = None,
      extra_map_manager_list = None,
      extra_map_manager_id_list = None):

    '''
      map_dict has four special ids with interpretations:
        map_manager:  full map
        map_manager_1, map_manager_2: half-maps 1 and 2
        map_manager_mask:  a mask in a map_manager
      All other ids are any strings and are assumed to correspond to other maps
      map_manager must be present
    '''


    assert (map_manager is not None) or (
       (map_manager_1 is not None) and (map_manager_2 is not None))
    self._map_dict={}
    self._map_dict['map_manager']=map_manager
    if map_manager_1 and map_manager_2:
      self._map_dict['map_manager_1']=map_manager_1
      self._map_dict['map_manager_2']=map_manager_2
    if extra_map_manager_id_list:
      for id, m in zip(extra_map_manager_id_list,extra_map_manager_list):
        if (id is not None) and (m is not None):
          self._map_dict[id]=m

  def set_up_model_dict(self,
      model = None,
      extra_model_list = None,
      extra_model_id_list = None):

    '''
      map_dict has one special id with interpretation:
        model:  standard model
      All other ids are any strings and are assumed to correspond to other
      models.
    '''

    self._model_dict={}
    self._model_dict['model']=model
    if extra_model_id_list:
      for id, m in zip(extra_model_id_list,extra_model_list):
        if id is not None and m is not None:
          self._model_dict[id]=m

  # prevent pickling error in Python 3 with self.log = sys.stdout
  # unpickling is limited to restoring sys.stdout
  def __getstate__(self):
    pickle_dict = self.__dict__.copy()
    if isinstance(self.log, io.TextIOWrapper):
      pickle_dict['log'] = None
    return pickle_dict

  def __setstate__(self, pickle_dict):
    self.__dict__ = pickle_dict
    if self.log is None:
      self.log = sys.stdout

  def __repr__(self):
    text = "\nMap_model_manager: \n"
    if self.model():
      text += "\n%s\n" %(str(self.model()))
    map_info = self._get_map_info()
    model_info = self._get_model_info()
    if self.map_manager():
      text += "\nmap_manager: %s\n" %(str(self.map_manager()))
    for id in map_info.other_map_id_list:
      text += "\n%s: %s\n" %(id,str(self.get_map_manager_by_id(id)))
    for id in model_info.other_model_id_list:
      text += "\n%s: %s\n" %(id,str(self.get_model_by_id(id)))
    return text

  # Methods for printing

  def set_log(self, log = sys.stdout):
    '''
       Set output log file
    '''
    if log is None:
      self.log = null_out()
    else:
      self.log = log

  def _print(self, m):
    '''
      Print to log if it is present
    '''

    if (self.log is not None) and hasattr(self.log, 'closed') and (
        not self.log.closed):
      print(m, file = self.log)

  # Methods for obtaining models, map_managers, symmetry, ncs_objects

  def crystal_symmetry(self):
    ''' Get the working crystal_symmetry'''
    return self.map_manager().crystal_symmetry()

  def unit_cell_crystal_symmetry(self):
    ''' Get the unit_cell_crystal_symmetry (full or original symmetry)'''
    return self.map_manager().unit_cell_crystal_symmetry()

  def shift_cart(self):
    ''' get the shift_cart (shift since original location)'''
    return self.map_manager().shift_cart()

  def map_dict(self):
    ''' Get the dictionary of all maps and masks as map_manager objects'''
    return self._map_dict

  def model_dict(self):
    ''' Get the dictionary of all models '''
    return self._model_dict

  def models(self):
    ''' Get all the models as a list'''
    model_list = []
    for id in self.model_id_list():
      m = self.get_model_by_id(id)
      if m is not None:
        model_list.append(m)
    return model_list

  def model(self):
    ''' Get the model '''
    return self._model_dict.get('model')

  def model_id_list(self):
    ''' Get all the names (ids) for all models'''
    mil = []
    for id in self.model_dict().keys():
      if self.get_model_by_id(id) is not None:
        mil.append(id)
    return mil

  def get_model_by_id(self, model_id):
    ''' Get a model with the name model_id'''
    return self.model_dict().get(model_id)

  def remove_model_by_id(self, model_id = 'extra'):
    '''
     Remove this model
   '''
    del self._model_dict[map_id]

  def map_managers(self):
    ''' Get all the map_managers as a list'''
    map_manager_list = []
    for id in self.map_id_list():
      mm = self.get_map_manager_by_id(id)
      if mm:
        map_manager_list.append(mm)
    return map_manager_list

  def map_manager(self):
    '''
      Get the map_manager

      If not present, calculate it from map_manager_1 and map_manager_2
      and set it.
    '''

    map_manager = self._map_dict.get('map_manager')


    if (not map_manager):
      # If map_manager_1 and map_manager_2 are supplied but no map_manager,
      #   create map_manager as average of map_manager_1 and map_manager_2

      map_manager_1 = self._map_dict.get('map_manager_1')
      map_manager_2 = self._map_dict.get('map_manager_2')
      if map_manager_1 and map_manager_2:

        map_manager = map_manager_1.customized_copy(map_data =
          0.5 * (map_manager_1.map_data() + map_manager_2.map_data()))
        self._map_dict['map_manager'] = map_manager

    return map_manager

  def map_manager_1(self):
    ''' Get half_map 1 as a map_manager object '''
    return self._map_dict.get('map_manager_1')

  def map_manager_2(self):
    ''' Get half_map 2 as a map_manager object '''
    return self._map_dict.get('map_manager_2')

  def map_manager_mask(self):
    ''' Get the mask as a map_manager object '''
    return self._map_dict.get('map_manager_mask')

  def map_id_list(self):
    ''' Get all the names (ids) for all map_managers that are present'''
    mil = []
    for id in self.map_dict().keys():
      if self.get_map_manager_by_id(id) is not None:
        mil.append(id)
    return mil

  def ncs_object(self):
    if self.map_manager():
      return self.map_manager().ncs_object()
    else:
      return None

  def experiment_type(self):
    if self.map_manager():
      return self.map_manager().experiment_type()
    else:
      return None

  def scattering_table(self):
    if self.map_manager():
      return self.map_manager().scattering_table()
    else:
      return None

  def resolution(self):
    if self.map_manager():
      return self.map_manager().resolution()
    else:
      return None

  def set_resolution(self, resolution):
    ''' Set nominal resolution '''
    # Must already have a map_manager
    assert self.map_manager() is not None
    self.map_manager().set_resolution(resolution)

  def set_scattering_table(self, scattering_table):
    ''' Set nominal scattering_table '''
    # Must already have a map_manager
    assert self.map_manager() is not None
    self.map_manager().set_scattering_table(scattering_table)

  def set_experiment_type(self, experiment_type):
    ''' Set nominal experiment_type '''
    # Must already have a map_manager
    assert self.map_manager() is not None
    self.map_manager().set_experiment_type(experiment_type)

  def _get_map_coeffs_list_from_id_list(self, id_list,
    mask_id = None):
    '''
      Get maps identified by map_id_list
      Optionally mask them with mask_id
      Return map_data from (masked) maps, converted to
        structure factors, as list
    '''
    map_data_list = self._get_map_data_list_from_id_list(id_list,
      mask_id = mask_id)
    map_coeffs_list = []
    from cctbx import miller
    for map_data in map_data_list:
      map_coeffs = miller.structure_factor_box_from_map(
        map              = map_data,
        crystal_symmetry = self.crystal_symmetry())
      map_coeffs_list.append(map_coeffs)
    return map_coeffs_list

  def _get_map_data_list_from_id_list(self, id_list,
    mask_id = None):
    '''
      Get maps identified by map_id_list
      Optionally mask them with mask_id
      Return map_data from (masked) maps as list
    '''

    map_data_list = []
    if mask_id is None: # just get the map_data
      for id in id_list:
        map_data_list.append(self.get_map_manager_by_id(id).map_data())
    else:
      assert mask_id in self.map_id_list() and \
        self.get_map_manager_by_id(mask_id).is_mask()
      # Create masked copies of all masks and get list of their id's
      new_map_id_list = self.create_masked_copies_of_maps(
         map_id_list = id_list,
         mask_id = mask_id)
      # Get their map data (masked)
      for id in new_map_id_list:
        map_data_list.append(self.get_map_manager_by_id(id).map_data())
      # Clean up dummy mask managers
      for id in new_map_id_list:
        self.remove_map_manager_by_id(id)
    return map_data_list

  def get_map_manager_by_id(self, map_id):
    '''
      Get a map_manager with the name map_id
      If map_id is 'map_manager' specifically return self.map_manager()
      so that it will create a map_manager from map_manager_1 and map_manager_2
      if map_manager is not present
    '''
    if map_id == 'map_manager':
      return self.map_manager()
    else:
      return self.map_dict().get(map_id)

  def get_map_data_by_id(self, map_id):
    ''' Get map_data from a map_manager with the name map_id'''
    map_manager = self.get_map_manager_by_id(map_id)
    if map_manager:
      return map_manager.map_data()
    else:
      return None

  def add_model_by_id(self, model, model_id,
     overwrite = True):
    '''
     Add a new model
     Must be similar to existing map_managers
     Overwrites any existing with the same id unless overwrite = False
    '''
    assert isinstance(model, mmtbx.model.manager)
    if not overwrite:
      assert not model_id in self.model_id_list() # must not duplicate
    if not self.map_manager().is_compatible_model(model): # needs shifting
      self.shift_any_model_to_match(model)
    self._model_dict[model_id] = model

  def add_map_manager_by_id(self, map_manager, map_id,
     overwrite = True):
    '''
     Add a new map_manager
     Must be similar to existing
     Overwrites any existing with the same id unless overwrite = False
     Is a mask if is_mask is set
    '''
    assert isinstance(map_manager, iotbx.map_manager.map_manager)
    assert isinstance(overwrite, bool)
    if not overwrite:
      assert not map_id in self.map_id_list() # must not duplicate
    assert map_manager.is_similar(self.map_manager())
    self._map_dict[map_id] = map_manager

  def remove_map_manager_by_id(self, map_id = 'extra'):
    '''
     Remove this map manager
     Note: you cannot remove 'map_manager' ... you can only replace it
   '''
    assert map_id != 'map_manager'
    del self._map_dict[map_id]


  def duplicate_map_manager(self,
    map_id = 'map_manager',
    new_map_id='new_map_manager'):
    '''
     Duplicate (deep_copy) map_manager
     Overwrites any existing with the new id
    '''
    map_manager = self.get_map_manager_by_id(map_id)
    assert isinstance(map_manager, iotbx.map_manager.map_manager)

    self._map_dict[new_map_id] = map_manager.deep_copy()


  # Methods for writing maps and models

  def write_map(self, file_name = None, id='map_manager'):
    if not self._map_dict.get(id):
      self._print ("No map to write out with id='%s'" %(id))
    elif not file_name:
      self._print ("Need file name to write map")
    else:
      self._map_dict.get(id).write_map(file_name = file_name)

  def write_model(self,
     file_name = None):
    if not self.model():
      self._print ("No model to write out")
    elif not file_name:
      self._print ("Need file name to write model")
    else:
      # Write out model

      f = open(file_name, 'w')
      print(self.model().model_as_pdb(), file = f)
      f.close()
      self._print("Wrote model with %s residues to %s" %(
         self.model().get_hierarchy().overall_counts().n_residues,
         file_name))

  # Methods for identifying which map_manager and model to use

  def _get_map_info(self):
    '''
      Return a group_args object specifying the map_manager and
      a list of any other maps present
    '''
    all_map_id_list=list(self._map_dict.keys())
    # We are going to need id='map_manager'   create if if missing
    assert self.map_manager() is not None # creates it
    assert all_map_id_list
    all_map_id_list.sort()
    map_id='map_manager'
    other_map_id_list=[]
    for id in all_map_id_list:
      if id != map_id:
        other_map_id_list.append(id)

    return group_args(map_id=map_id,
         other_map_id_list=other_map_id_list)

  def _get_model_info(self):
    '''
      Return a group_args object specifying the model and
      a list of any other models present
    '''
    all_model_id_list=list(self._model_dict.keys())
    assert all_model_id_list
    all_model_id_list.sort()
    model_id='model'
    other_model_id_list=[]
    for id in all_model_id_list:
      if id != model_id:
        other_model_id_list.append(id)
    if not model_id in all_model_id_list:
      model_id = None

    return group_args(model_id=model_id,
         other_model_id_list=other_model_id_list)

  # Methods for manipulation of maps

  def initialize_maps(self, map_value = 0):
    '''
      Set values of all maps to map_value
      Used to set up an empty set of maps for filling in from boxes
    '''

    for mm in self.map_managers():
      mm.initialize_map_data(map_value = map_value)

  # Methods for boxing maps (changing the dimensions of the maps)
  # box_all...methods change the contents of the current object (they do not
  #  create a new object)
  # extract_all... methods make a new object

  def extract_all_maps_with_bounds(self,
     lower_bounds,
     upper_bounds):
    '''
      Runs box_all_maps_with_bounds_and_shift_origin with extract_box=True
    '''
    return self.box_all_maps_with_bounds_and_shift_origin(
      lower_bounds = lower_bounds,
      upper_bounds = upper_bounds,
      extract_box = True)

  def box_all_maps_with_bounds_and_shift_origin(self,
     lower_bounds,
     upper_bounds,
     extract_box = False):
    '''
       Box all maps using specified bounds, shift origin of maps, model
       Replaces existing map_managers and shifts model in place

       If extract_box=True:  Creates new object with deep_copies.
       Otherwise: replaces existing map_managers and shifts model in place

       NOTE: This changes the gridding and shift_cart of the maps and model

       Can be used in map_model_manager to work with boxed maps
       and model or in map_model_manager to re-box all maps and model

       The lower_bounds and upper_bounds define the region to be boxed. These
       bounds are relative to the current map with origin at (0, 0, 0).

    '''
    assert lower_bounds is not None and upper_bounds is not None
    assert len(tuple(lower_bounds)) == 3
    assert len(tuple(upper_bounds)) == 3

    from cctbx.maptbx.box import with_bounds

    map_info=self._get_map_info()
    map_manager = self._map_dict[map_info.map_id]
    assert map_manager is not None

    model_info=self._get_model_info()
    model = self._model_dict[model_info.model_id]

    if extract_box: # make sure everything is deep_copy
      model = model.deep_copy()

    # Make box with bounds and apply it to model, first map
    box = with_bounds(
      map_manager = self._map_dict[map_info.map_id],
      lower_bounds = lower_bounds,
      upper_bounds = upper_bounds,
      model = model,
      wrapping = self._force_wrapping,
      log = self.log)
    # Now box is a copy of map_manager and model that is boxed

    # Now apply boxing to other maps and models and then insert them into
    #  either this map_model_manager object, replacing what is there (extract_box=False)
    #  or create and return a new map_model_manager object (extract_box=True)
    return self._finish_boxing(box = box, model_info = model_info,
      map_info = map_info,
      extract_box = extract_box)

  def extract_all_maps_around_model(self,
     selection_string = None,
     select_unique_by_ncs = False,
     box_cushion = 5.):
    '''
      Runs box_all_maps_around_model_and_shift_origin with extract_box=True
    '''
    return self.box_all_maps_around_model_and_shift_origin(
      selection_string = selection_string,
      box_cushion = box_cushion,
      select_unique_by_ncs = select_unique_by_ncs,
      extract_box = True)

  def box_all_maps_around_model_and_shift_origin(self,
     selection_string = None,
     box_cushion = 5.,
     select_unique_by_ncs = False,
     extract_box = False):
    '''
       Box all maps around the model, shift origin of maps, model
       If extract_box=True:  Creates new object with deep_copies.
       Otherwise: replaces existing map_managers and shifts model in place

       NOTE: This changes the gridding and shift_cart of the maps and model

       Can be used in map_model_manager to work with boxed maps
       and model or in map_model_manager to re-box all maps and model

       Requires a model

       The box_cushion defines how far away from the nearest atoms the new
       box boundaries will be placed

       The selection_string defines what part of the model to keep ('ALL' is
        default)
       If select_unique_by_ncs is set, select the unique part of the model
       automatically.  Any selection in selection_string will not be applied.
    '''
    assert isinstance(self.model(), model_manager)
    assert box_cushion is not None

    from cctbx.maptbx.box import around_model

    map_info=self._get_map_info()
    assert map_info.map_id is not None
    model_info=self._get_model_info()
    assert model_info.model_id is not None # required for box_around_model
    model = self._model_dict[model_info.model_id]

    if select_unique_by_ncs:
      model.search_for_ncs()
      sel = model.get_master_selection()
      model = model.select(sel)
    elif selection_string:
      sel = model.selection(selection_string)
      model = model.select(sel)
    elif extract_box: # make sure everything is deep_copy
      model = model.deep_copy()

    # Make box around model and apply it to model, first map
    # This step modifies model in place and creates a new map_manager
    box = around_model(
      map_manager = self._map_dict[map_info.map_id],
      model = model,
      box_cushion = box_cushion,
      wrapping = self._force_wrapping,
      log = self.log)
    # Now box is a copy of map_manager and model that is boxed

    # Now apply boxing to other maps and models and then insert them into
    #  either this map_model_manager object, replacing what is there (extract_box=False)
    #  or create and return a new map_model_manager object (extract_box=True)
    return self._finish_boxing(box = box, model_info = model_info,
      map_info = map_info,
      extract_box = extract_box)

  def extract_all_maps_around_density(self,
     box_cushion = 5.,
     threshold = 0.05,
     get_half_height_width = True,
     map_id = 'map_manager'):
    '''
      Runs box_all_maps_around_density_and_shift_origin with extract_box=True
    '''
    return self.box_all_maps_around_density_and_shift_origin(
     box_cushion = box_cushion,
     threshold = threshold,
     get_half_height_width = get_half_height_width,
     map_id = map_id,
     extract_box = True)

  def box_all_maps_around_density_and_shift_origin(self,
     box_cushion = 5.,
     threshold = 0.05,
     map_id = 'map_manager',
     get_half_height_width = True,
     extract_box = False):
    '''
       Box all maps around the density in map_id map (default is map_manager)
       shift origin of maps, model

       If extract_box=True:  Creates new object with deep_copies.
       Otherwise: replaces existing map_managers and shifts model in place

       Replaces existing map_managers and shifts model in place

       NOTE: This changes the gridding and shift_cart of the maps and model

       Can be used in map_model_manager to work with boxed maps
       and model or in map_model_manager to re-box all maps and model

       Does not require a model, but a model can be supplied.  If model is
       supplied, it is possible that the model will be outside the density
       after boxing.
       To avoid this, use box_all_maps_around_model_and_shift_origin instead.

       The box_cushion defines how far away from the nearest density the new
       box boundaries will be placed

       The threshold defines how much (relative to maximum in map)  above
       mean value of map near edges is significant and should count as density.

    '''
    assert box_cushion is not None

    from cctbx.maptbx.box import around_density

    map_info=self._get_map_info()
    assert map_info.map_id is not None
    model_info=self._get_model_info()
    model = self._model_dict[model_info.model_id]
    if extract_box: # make sure everything is deep_copy
      model = model.deep_copy()

    # Make box around model and apply it to model, first map
    box = around_density(
      map_manager = self._map_dict[map_info.map_id],
      model       = model,
      box_cushion = box_cushion,
      threshold   = threshold,
      get_half_height_width = get_half_height_width,
      wrapping    = self._force_wrapping)

    # Now box is a copy of map_manager and model that is boxed

    # Now apply boxing to other maps and models and then insert them into
    #  either this map_model_manager object, replacing what is there (extract_box=False)
    #  or create and return a new map_model_manager object (extract_box=True)
    return self._finish_boxing(box = box, model_info = model_info,
      map_info = map_info,
      extract_box = extract_box)

  def extract_all_maps_around_mask(self,
     box_cushion = 5.,
     mask_id = 'mask'):
    '''
      Runs box_all_maps_around_mask_and_shift_origin with extract_box=True
    '''
    return self.box_all_maps_around_mask_and_shift_origin(
     box_cushion = 5.,
     mask_id = mask_id,
     extract_box = True)

  def box_all_maps_around_mask_and_shift_origin(self,
     box_cushion = 5.,
     mask_id = 'mask',
     extract_box = False):
    '''
       Box all maps around specified mask, shift origin of maps, model
       Replaces existing map_managers and shifts model in place

       If extract_box=True:  Creates new object with deep_copies.
       Otherwise: replaces existing map_managers and shifts model in place

       NOTE: This changes the gridding and shift_cart of the maps and model

       Requires a mask

       The box_cushion defines how far away from the edge of the mask the new
       box boundaries will be placed

    '''
    assert isinstance(self.model(), model_manager)
    assert box_cushion is not None

    from cctbx.maptbx.box import around_mask

    map_info=self._get_map_info()
    assert map_info.map_id is not None
    map_manager = self._map_dict[map_info.map_id]

    mask_mm = self.get_map_manager_by_id(mask_id)
    assert mask_mm is not None
    assert mask_mm.is_mask()

    assert mask_mm is not map_manager  # mask and map cannot be the same

    model_info=self._get_model_info()
    model = self._model_dict[model_info.model_id]
    if extract_box: # make sure everything is deep_copy
      model = model.deep_copy()

    # Make box around mask and apply it to model, first map
    box = around_mask(
      map_manager = map_manager,
      mask_as_map_manager = mask_mm,
      model = model,
      box_cushion = box_cushion,
      wrapping = self._force_wrapping,
      log = self.log)
    # Now box is a copy of map_manager and model that is boxed

    # Now apply boxing to other maps and models and then insert them into
    #  either this map_model_manager object, replacing what is there (extract_box=False)
    #  or create and return a new map_model_manager object (extract_box=True)
    return self._finish_boxing(box = box, model_info = model_info,
      map_info = map_info,
      extract_box = extract_box)

  def extract_all_maps_around_unique(self,
     resolution = None,
     solvent_content = None,
     sequence = None,
     molecular_mass = None,
     soft_mask = True,
     chain_type = 'PROTEIN',
     box_cushion = 5,
     target_ncs_au_model = None,
     regions_to_keep = None,
     keep_low_density = True,
     symmetry = None,
     mask_expand_ratio = 1):

    '''
      Runs box_all_maps_around_mask_and_shift_origin with extract_box=True
    '''
    return self.box_all_maps_around_unique_and_shift_origin(
     resolution = resolution,
     solvent_content = solvent_content,
     sequence = sequence,
     molecular_mass = molecular_mass,
     soft_mask = soft_mask,
     chain_type = chain_type,
     box_cushion = box_cushion,
     target_ncs_au_model = target_ncs_au_model,
     regions_to_keep = regions_to_keep,
     keep_low_density = keep_low_density,
     symmetry = symmetry,
     mask_expand_ratio = mask_expand_ratio,
     extract_box = True)

  def box_all_maps_around_unique_and_shift_origin(self,
     resolution = None,
     solvent_content = None,
     sequence = None,
     molecular_mass = None,
     soft_mask = True,
     chain_type = 'PROTEIN',
     box_cushion = 5,
     target_ncs_au_model = None,
     regions_to_keep = None,
     keep_low_density = True,
     symmetry = None,
     mask_expand_ratio = 1,
     extract_box = False):
    '''
       Box all maps using bounds obtained with around_unique,
       shift origin of maps, model, and mask around unique region

       If extract_box=True:  Creates new object with deep_copies.
       Otherwise: replaces existing map_managers and shifts model in place

       Replaces existing map_managers and shifts model in place

       NOTE: This changes the gridding and shift_cart of the maps and model
       and masks the map

       Normally supply just sequence; resolution will be taken from
       map_manager resolution if present.  other options match
       all possible ways that segment_and_split_map can estimate solvent_content

       Must supply one of (sequence, solvent_content, molecular_mass)

       Symmetry is optional symmetry (i.e., D7 or C1). Used as alternative to
       ncs_object supplied in map_manager


       Additional parameters:
         mask_expand_ratio:   allows increasing masking radius beyond default at
                              final stage of masking
         solvent_content:  fraction of cell not occupied by macromolecule
         sequence:        one-letter code of sequence of unique part of molecule
         chain_type:       PROTEIN or RNA or DNA. Used with sequence to estimate
                            molecular_mass
         molecular_mass:    Molecular mass (Da) of entire molecule used to
                            estimate solvent_content
         target_ncs_au_model: model marking center of location to choose as
                              unique
         box_cushion:        buffer around unique region to be boxed
         soft_mask:  use soft mask
         keep_low_density:  keep low density regions
         regions_to_keep:   Allows choosing just highest-density contiguous
                            region (regions_to_keep=1) or a few
    '''
    from cctbx.maptbx.box import around_unique

    map_info=self._get_map_info()
    map_manager = self._map_dict[map_info.map_id]
    assert isinstance(map_manager, iotbx.map_manager.map_manager)
    if not resolution:
      resolution = self.resolution()
    assert resolution is not None
    assert (sequence, solvent_content, molecular_mass).count(None) == 2

    model_info=self._get_model_info()
    model = self._model_dict[model_info.model_id]
    if extract_box: # make sure everything is deep_copy
      model = model.deep_copy()

    # Make box with around_unique and apply it to model, first map
    box = around_unique(
      map_manager = map_manager,
      model = model,
      wrapping = self._force_wrapping,
      target_ncs_au_model = target_ncs_au_model,
      regions_to_keep = regions_to_keep,
      solvent_content = solvent_content,
      resolution = resolution,
      sequence = sequence,
      molecular_mass = molecular_mass,
      symmetry = symmetry,
      chain_type = chain_type,
      box_cushion = box_cushion,
      soft_mask = soft_mask,
      mask_expand_ratio = mask_expand_ratio,
      log = self.log)

    # Now box is a copy of map_manager and model that is boxed

    # Now apply boxing to other maps and models and then insert them into
    #  either this map_model_manager object, replacing what is there (extract_box=False)
    #  or create and return a new map_model_manager object (extract_box=True)
    other = self._finish_boxing(box = box, model_info = model_info,
      map_info = map_info,
      extract_box = extract_box)

    if not extract_box:
      other = self #  modifying this object

    # Now apply masking to all other maps (not done in _finish_boxing)
    for id in map_info.other_map_id_list:
      other._map_dict[id] = box.apply_extract_unique_mask(
        self._map_dict[id],
        resolution = resolution,
        soft_mask = soft_mask)

    if extract_box:
      return other

  def _finish_boxing(self, box, model_info, map_info,
    extract_box = False):

    '''
       If extract_box is False, modify this object in place.
       If extract_box is True , create a new object of same type and return it
    '''

    if box.warning_message():
      self._warning_message = box.warning_message()
      self._print("%s" %(box.warning_message()))

    if extract_box:
      other = self._empty_copy() # making a new object
    else:
      other = self #  modifying this object


    other._map_dict[map_info.map_id] = box.map_manager()
    other._model_dict[model_info.model_id] = box.model()

    # Apply the box to all the other maps
    for id in map_info.other_map_id_list:
      other._map_dict[id] = box.apply_to_map(self._map_dict[id])

    # Apply the box to all the other models
    for id in model_info.other_model_id_list:
      other._model_dict[id] = box.apply_to_model(self._model_dict[id])

    if extract_box:
      return other

  # Methods for masking maps ( creating masks and applying masks to maps)
  # These methods change the contents of the current object (they do not
  #  create a new object)

  def mask_all_maps_around_atoms(self,
      mask_atoms_atom_radius = 3,
      set_outside_to_mean_inside = False,
      soft_mask = False,
      soft_mask_radius = None,
      mask_id = 'mask'):
    assert mask_atoms_atom_radius is not None
    assert self.model() is not None

    '''
      Generate mask around atoms and apply to all maps.
      Overwrites values in these maps

      NOTE: Does not change the gridding or shift_cart of the maps and model

      Optionally set the value outside the mask equal to the mean inside,
        changing smoothly from actual values inside the mask to the constant
        value outside (otherwise outside everything is set to zero)

      Optional: radius around atoms for masking
      Optional: soft mask  (default = True)
        Radius will be soft_mask_radius
        (default radius is self.resolution() or resolution calculated
          from gridding)
        If soft mask is set, mask_atoms_atom_radius increased by

      Optionally use any mask specified by mask_id
    '''
    if soft_mask:
      if not soft_mask_radius:
        from cctbx.maptbx import d_min_from_map
        soft_mask_radius = d_min_from_map(
          map_data=self.map_manager().map_data(),
          unit_cell=self.map_manager().crystal_symmetry().unit_cell())
    self.create_mask_around_atoms(
         soft_mask = soft_mask,
         soft_mask_radius = soft_mask_radius,
         mask_atoms_atom_radius = mask_atoms_atom_radius)
    self.apply_mask_to_maps(mask_id = mask_id,
         set_outside_to_mean_inside = \
           set_outside_to_mean_inside)

  def mask_all_maps_around_edges(self,
      soft_mask_radius = None,
      mask_id = 'mask'):
    '''
      Apply a soft mask around edges of all maps. Overwrites values in maps
      Use 'mask' as the mask id

      NOTE: Does not change the gridding or shift_cart of the maps and model
    '''

    self.create_mask_around_edges(soft_mask_radius = soft_mask_radius,
      mask_id = mask_id)
    self.apply_mask_to_maps(mask_id = mask_id)

  def mask_all_maps_around_density(self,
     solvent_content = None,
     soft_mask = True,
     soft_mask_radius = None,
     mask_id = 'mask',
     map_id = 'map_manager'):
    '''
      Apply a soft mask around density.  Mask calculated using map_id and
      written to mask_id . Overwrites values in maps
      Default is to use 'mask' as the mask id

      NOTE: Does not change the gridding or shift_cart of the maps and model
    '''

    self.create_mask_around_density(
      solvent_content = solvent_content,
      soft_mask  = soft_mask,
      soft_mask_radius = soft_mask_radius,
      mask_id = mask_id,
      map_id = map_id)
    self.apply_mask_to_maps(mask_id = mask_id)

  def create_masked_copies_of_maps(self,
    map_id_list = None,
    mask_id = 'mask'):
   '''
    Create masked copies of all maps identified by map_id_list (default is all)
    Return list of map_id for masked versions
   '''

   new_map_id_list = []
   for id in list(map_id_list):
     new_id = self._generate_new_map_id()
     self.duplicate_map_manager(id,new_id)
     self.apply_mask_to_map(map_id=new_id, mask_id = mask_id)
     new_map_id_list.append(new_id)
   return new_map_id_list

  def apply_mask_to_map(self,
      map_id,
      mask_id = 'mask',
      set_outside_to_mean_inside = False):
    '''
      Apply the mask in 'mask' to map specified by map_id

      Optionally set the value outside the mask equal to the mean inside,
        changing smoothly from actual values inside the mask to the constant
        value outside (otherwise outside everything is set to zero)

      Optionally use any mask specified by mask_id

      NOTE: Does not change the gridding or shift_cart of the map
    '''

    self.apply_mask_to_maps(map_ids = [map_id],
      mask_id = mask_id,
      set_outside_to_mean_inside = set_outside_to_mean_inside)

  def apply_mask_to_maps(self,
      map_ids = None,
      mask_id = 'mask',
      set_outside_to_mean_inside = False):
    '''
      Apply the mask in 'mask' to maps specified by map_ids.
      If map_ids is None apply to all

      Optionally set the value outside the mask equal to the mean inside,
        changing smoothly from actual values inside the mask to the constant
        value outside (otherwise outside everything is set to zero)

      Optionally use any mask specified by mask_id

      NOTE: Does not change the gridding or shift_cart of the maps
    '''

    assert (map_ids is None) or isinstance(map_ids, list)
    assert isinstance(set_outside_to_mean_inside, bool)
    mask_mm = self.get_map_manager_by_id(mask_id)
    assert mask_mm is not None
    assert mask_mm.is_mask()

    from cctbx.maptbx.segment_and_split_map import apply_mask_to_map

    if map_ids is None:
      map_ids = list(self._map_dict.keys())
    for map_id in map_ids:
      mm=self.get_map_manager_by_id(map_id)
      if mm.is_mask(): continue  # don't apply to a mask
      if set_outside_to_mean_inside in [None, True]:
      # smoothly go from actual value inside mask to target value outside
        new_map_data = apply_mask_to_map(mask_data = mask_mm.map_data(),
          set_outside_to_mean_inside = set_outside_to_mean_inside,
          map_data = mm.map_data(),
          out = null_out())
      else: # Simple case  just multiply
        new_map_data = mask_mm.map_data() * mm.map_data()
      # Set the values in this manager
      mm.set_map_data(map_data = new_map_data)

  def create_mask_around_edges(self,
     soft_mask_radius = None,
     mask_id = 'mask' ):
    '''
      Generate new mask map_manager with soft mask around edges of mask
      Does not apply the mask to anything.
      Normally follow with apply_mask_to_map or apply_mask_to_maps

      Optional: radius around edge for masking
        (default radius is self.resolution() or resolution calculated
         from gridding)

      Generates new entry in map_manager dictionary with id of
      mask_id (default='mask') replacing any existing entry with that id
    '''

    if not soft_mask_radius:
      from cctbx.maptbx import d_min_from_map
      soft_mask_radius = d_min_from_map(
         map_data=self.map_manager().map_data(),
         unit_cell=self.map_manager().crystal_symmetry().unit_cell())

    from cctbx.maptbx.mask import create_mask_around_edges
    cm = create_mask_around_edges(map_manager = self.map_manager(),
      soft_mask_radius = soft_mask_radius)
    cm.soft_mask(soft_mask_radius = soft_mask_radius)

    # Put the mask in map_dict ided with mask_id
    self.add_map_manager_by_id(map_manager = cm.map_manager(),
      map_id = mask_id)

  def create_mask_around_atoms(self,
     mask_atoms_atom_radius = 3,
     soft_mask = False,
     soft_mask_radius = None,
     mask_id = 'mask' ):

    '''
      Generate mask based on model.  Does not apply the mask to anything.
      Normally follow with apply_mask_to_map or apply_mask_to_maps

      Optional: radius around atoms for masking
      Optional: soft mask  (default = True)
        Radius will be soft_mask_radius
        (default radius is self.resolution() or resolution calculated
           from gridding)
        If soft mask is set, mask_atoms_atom_radius increased by
          soft_mask_radius

      Generates new entry in map_manager dictionary with id of
      mask_id (default='mask') replacing any existing entry with that id
    '''

    if soft_mask:
      if not soft_mask_radius:
        from cctbx.maptbx import d_min_from_map
        soft_mask_radius = d_min_from_map(
           map_data=self.map_manager().map_data(),
           unit_cell=self.map_manager().crystal_symmetry().unit_cell())
      mask_atoms_atom_radius += soft_mask_radius

    from cctbx.maptbx.mask import create_mask_around_atoms
    cm = create_mask_around_atoms(map_manager = self.map_manager(),
      model = self.model(),
      mask_atoms_atom_radius = mask_atoms_atom_radius)

    if soft_mask: # Make the create_mask object contain a soft mask
      cm.soft_mask(soft_mask_radius = soft_mask_radius)

    # Put the mask in map_dict ided with mask_id
    self.add_map_manager_by_id(map_manager = cm.map_manager(),
      map_id = mask_id)

  def create_mask_around_density(self,
     resolution = None,
     solvent_content = None,
     soft_mask = True,
     soft_mask_radius = None,
     mask_id = 'mask',
     map_id = 'map_manager' ):

    '''
      Generate mask based on density in map_manager (map_id defines it).
      Does not apply the mask to anything.
      Normally follow with apply_mask_to_map or apply_mask_to_maps

      Optional:  supply working resolution
      Optional:  supply approximate solvent fraction

      Optional: soft mask  (default = True)
        Radius will be soft_mask_radius
        (default radius is resolution calculated from gridding)

      Generates new entry in map_manager dictionary with id of
      mask_id (default='mask') replacing any existing entry with that id
    '''

    assert solvent_content is None or isinstance(solvent_content, (int, float))
    assert soft_mask_radius is None or \
       isinstance(soft_mask_radius, (int, float))
    assert isinstance(soft_mask, bool)

    if not resolution:
      resolution = self.resolution()

    map_manager = self.get_map_manager_by_id(map_id)
    assert map_manager is not None # Need a map to create mask around density
    from cctbx.maptbx.mask import create_mask_around_density
    cm = create_mask_around_density(map_manager = map_manager,
        solvent_content = solvent_content,
        resolution = resolution)

    if soft_mask: # Make the create_mask object contain a soft mask
      if not soft_mask_radius:
        if resolution:
          soft_mask_radius = resolution
        else:
          from cctbx.maptbx import d_min_from_map
          soft_mask_radius = d_min_from_map(
           map_data=map_manager.map_data(),
           unit_cell=map_manager.crystal_symmetry().unit_cell())
      cm.soft_mask(soft_mask_radius = soft_mask_radius)

    # Put the mask in map_dict id'ed with mask_id
    self.add_map_manager_by_id(map_manager = cm.map_manager(),
      map_id = mask_id)

  # Methods for recombining models

  def propagate_model_from_other(self, other,
     model_id = 'model',
     other_model_id = 'model'):
    '''
    Import a model from other with get_model_from_other (other_model_id),
    then set coordinates of corresponding atoms in model_id

    The model in other must have been extracted from the model in this object
    or one just like it with select_unique_by_ncs=True, and no atoms can
    have been added or removed.

    '''

    if not self.model():
      return  # nothing to do

    # Get the imported hierarchy, shifted to match location of working one
    ph_imported_unique = self.get_model_from_other(other,
       other_model_id = other_model_id).get_hierarchy()

    # Get unique part of working hierarchy. Note this is not a deep_copy,
    #  so modifying it changes original hierarchy and can be propagated

    model = self.get_model_by_id(model_id)
    model.search_for_ncs()
    ph_working_unique = model.get_master_hierarchy()
    assert ph_imported_unique.is_similar_hierarchy(
       ph_working_unique) # hierarchies must match

    # Replace the coordinates in ph_working_unique with ph_imported_unique

    new_coords=ph_imported_unique.atoms().extract_xyz()
    ph_working_unique.atoms().set_xyz(new_coords)

    # And propagate these sites to rest of molecule with internal ncs
    model.set_sites_cart_from_hierarchy(multiply_ncs=True)

  def shift_any_model_to_match(self, model):
    '''
    Take any model and shift it to match the working shift_cart
    Also sets crystal_symmetry.
    Changes model in place

     Parameters:  model
    '''
    assert isinstance(model, mmtbx.model.manager)
    if not model.shift_cart():
      model.set_shift_cart((0, 0, 0))
    coordinate_shift = tuple(
      [s - o for s,o in zip(self.shift_cart(),model.shift_cart())])
    model.shift_model_and_set_crystal_symmetry(
        shift_cart = coordinate_shift,
        crystal_symmetry=self.crystal_symmetry())

  def get_model_from_other(self, other,
     other_model_id = 'model'):
    '''
     Take a model with id other_model_id from other_map_model_manager with any
     boxing and origin shifts allowed, and put it in the same reference
     frame as the current model.  Used to build up a model from pieces
     that were worked on in separate boxes.

     Changes model from other in place

     Parameters:  other:  Other map_model_manager containing a model
    '''
    assert isinstance(other, map_model_manager)
    other_model = other.get_model_by_id(other_model_id)
    assert other_model is not None # Need model for get_model_from_other
    coordinate_shift = tuple(
      [s - o for s,o in zip(self.shift_cart(),other.shift_cart())])
    other_model.shift_model_and_set_crystal_symmetry(
        shift_cart = coordinate_shift)
    matched_other_model = other_model
    return matched_other_model

  # Methods for producing Fourier coefficients and calculating maps

  def map_as_fourier_coefficients(self,
      d_min = None,
      d_max = None,
      map_id = 'map_manager'):
    '''
     Return Miller array to resolution specified based on map with id map_id

     Note that the map_manager is always zero-based (origin at (0,0,0)).
     The Fourier coefficients represent the map in this location at (0, 0, 0)
    '''

    # Checks
    map_manager = self.get_map_manager_by_id(map_id)
    assert map_manager is not None

    return map_manager.map_as_fourier_coefficients(
      d_min = d_min,
      d_max = d_max,
      )

  def add_map_from_fourier_coefficients(self,
      map_coeffs,
      map_id = 'map_from_fourier_coefficients'):
    '''
     Create map_manager from map_coeffs and add it to maps with map_id
     The map_coeffs must refer to a map with origin at (0, 0, 0) such as
     is produced by map_as_fourier_coefficients.

    '''

    # Checks
    map_manager = self.map_manager()
    assert map_manager is not None

    new_map_manager = map_manager.fourier_coefficients_as_map_manager(map_coeffs)
    self.add_map_manager_by_id(map_manager = new_map_manager,
      map_id = map_id)


  def resolution_filter(self,
      d_min = None,
      d_max = None,
      map_id = 'map_manager',
      ):
    '''
      Resolution-filter a map with range of d_min to d_max

      Typically used along with duplicate_map_manager to create a new map and
      filter it:
        rm.duplicate_map_manager(map_id='map_manager',
          new_map_id='resolution_filtered')
        rm.resolution_filter(map_id = 'resolution_filtered',)

    '''
    assert d_min is None or isinstance(d_min, (int,float))
    assert d_max is None or isinstance(d_max, (int,float))

    assert (d_min,d_max).count(None) < 2 # need some limits

    map_coeffs = self.map_as_fourier_coefficients(map_id = map_id,
      d_min = d_min,
      d_max = d_max)

    self.add_map_from_fourier_coefficients(map_coeffs,
      map_id = map_id)


  # Methods for modifying model or map

  def remove_model_outside_map(self, boundary, return_as_new_model=False):
    '''
     Remove all the atoms in the model that are well outside the map (more
     than boundary)
    '''
    assert boundary is not None

    if not self.model():
      return

    sites_frac = self.model().get_sites_frac()
    bf_a, bf_b, bf_c = self.model().crystal_symmetry().unit_cell(
        ).fractionalize((boundary, boundary, boundary))
    ub_a, ub_b, ub_c  = (1+bf_a, 1+bf_b, 1+bf_c)
    x,y,z = sites_frac.parts()
    s = (
         (x < -bf_a ) |
         (y < -bf_b) |
         (z < -bf_c) |
         (x > ub_a ) |
         (y > ub_b) |
         (z > ub_c)
         )
    if return_as_new_model:
      return self.model().select(~s)
    else:  # usual
      self.add_model_by_id( self.model().select(~s), 'model')


  # Methods for comparing maps, models and calculating FSC values

  def map_map_fsc(self,
      map_id_1 = 'map_manager',
      map_id_2 = 'map_manager',
      resolution = None,
      mask_id = None,
      mask_cutoff = 0.5,
      min_bin_width = 20,
      n_bins = 2000,
      fsc_cutoff = 0.143):
    '''
      Return the map-map FSC for these two maps, optionally masked with mask_id
      Returns fsc object which contains d_min which is d_min where fsc
        drops to fsc_cutoff, and sub-object fsc with arrays d, d_inv and
        fsc which are the FSC curve

    '''

    if not resolution:
      resolution = self.resolution()
    assert isinstance(resolution, (int, float))

    f_map_1, f_map_2 = self._get_map_coeffs_list_from_id_list(
      id_list = [map_id_1, map_id_2],
      mask_id = mask_id)

    bin_width=max(min_bin_width,int(0.5+f_map_1.size()/n_bins))

    # Get the FSC between map1 and map2
    fsc_curve = f_map_1.d_min_from_fsc(
        other = f_map_2, bin_width = bin_width, fsc_cutoff = fsc_cutoff)

    return fsc_curve


  def map_map_cc(self,
      map_id = 'map_manager_1',
      other_map_id = 'map_manager_2',
      mask_id = None,
      mask_cutoff = 0.5):

   map_map_info = self._get_map_map_info(
     map_id = map_id,
     other_map_id = other_map_id,
     mask_id = mask_id,
     mask_cutoff = mask_cutoff)
   return flex.linear_correlation(map_map_info.map_data_1d_1,
     map_map_info.map_data_1d_2).coefficient()

  def _get_map_map_info(self,
     map_id = None,
     other_map_id = None,
     mask_id = None,
     mask_cutoff = None):

   '''
     Check inputs and return selected parts of the two maps
   '''
   map1 = self.get_map_manager_by_id(map_id)
   map2 = self.get_map_manager_by_id(other_map_id)
   assert map1 and map2

   # Get the selection if any
   mask_map_manager = self.get_map_manager_by_id(mask_id)
   if mask_map_manager:
     assert mask_map_manager.is_mask()
     mask_data = mask_map_manager.map_data()
     sel = (mask_data.as_1d() > mask_cutoff)
     map_data_1d_1 = map1.map_data().as_1d().select(sel)
     map_data_1d_2 = map2.map_data().as_1d().select(sel)
   else:

     map_data_1d_1 = map1.map_data().as_1d()
     map_data_1d_2 = map2.map_data().as_1d()
   return group_args(
    map_data_1d_1 = map_data_1d_1,
    map_data_1d_2 = map_data_1d_2)


  def map_model_cc(self,
      resolution = None,
      map_id = 'map_manager',
      model_id = 'model',
      selection_string = None,
      model = None,
      ):

    if not model:
      model = self.get_model_by_id(model_id)
    if not model:
      return None
    map_manager= self.get_map_manager_by_id(map_id)
    assert model and map_manager
    if not resolution:
      resolution = self.resolution()
    assert resolution is not None

    if selection_string:
      sel = model.selection(selection_string)
      model = model.select(sel)

    import mmtbx.maps.correlation
    five_cc = mmtbx.maps.correlation.five_cc(
      map               = map_manager.map_data(),
      xray_structure    = model.get_xray_structure(),
      d_min             = resolution,
      compute_cc_mask   = True,
      compute_cc_box    = False,
      compute_cc_image  = False,
      compute_cc_volume = False,
      compute_cc_peaks  = False,)

    return five_cc.result.cc_mask

  # General methods

  def set_original_origin_grid_units(self, original_origin_grid_units = None):
    '''
     Reset (redefine) the original origin of the maps and models (apply an
      origin shift in effect).

     Procedure is: calculate shift_cart and set origin_shift_grid_units and
       shift_cart everywhere

    '''
    assert self.map_manager() is not None

    shift_cart=self.map_manager().grid_units_to_cart(
      [-i for i in original_origin_grid_units])
    for model in self.models():
      model.set_shift_cart(shift_cart)
    for map_manager in self.map_managers():
      map_manager.set_original_origin_and_gridding(
      original_origin=original_origin_grid_units)


  def _generate_new_map_id(self):
    '''
     Create a unique map_id
    '''
    used_id_list = self.map_id_list()
    i = 0
    while True:
      i += 1
      id = "temp_%s" %(i)
      if not id in used_id_list:
        return id

  def _generate_new_model_id(self):
    '''
     Create a unique model_id
    '''
    used_id_list = self.model_id_list()
    i = 0
    while (True):
      id = "temp_%s" %(i)
      if not id in used_id_list:
        return id

  def warning_message(self):
    return self._warning_message

  def show_summary(self, log = sys.stdout):
    text = self.__repr__()
    print (text, file = log)

  # Methods for accessing map_data, xrs, hierarchy directly.
  #  Perhaps remove all these

  def map_data(self):
    return self.map_manager().map_data()

  def map_data_1(self):
    if self.map_manager_1():
      return self.map_manager_1().map_data()

  def map_data_2(self):
    if self.map_manager_2():
      return self.map_manager_2().map_data()

  def map_data_list(self):
    map_data_list = []
    for mm in self.map_managers():
      map_data_list.append(mm.map_data())
    return map_data_list

  def xray_structure(self):
    if(self.model() is not None):
      return self.model().get_xray_structure()
    else:
      return None

  def hierarchy(self): return self.model().get_hierarchy()


  # Methods to be removed

  def get_counts_and_histograms(self):
    self._counts = get_map_counts(
      map_data         = self.map_data(),
      crystal_symmetry = self.crystal_symmetry())
    self._map_histograms = get_map_histograms(
        data    = self.map_data(),
        n_slots = 20,
        data_1  = self.map_data_1(),
        data_2  = self.map_data_2())

  def counts(self):
    if not hasattr(self, '_counts'):
      self.get_counts_and_histograms()
    return self._counts

  def histograms(self):
    if not hasattr(self, '_map_histograms'):
      self.get_counts_and_histograms()
    return self._map_histograms

  #  Convenience methods

  def generate_map(self,
      d_min = None,
      origin_shift_grid_units = None,
      file_name = None,
      model = None,
      n_residues = None,
      b_iso = 30,
      box_cushion = 5,
      scattering_table = 'electron',
      fractional_error = 0.0,
      gridding = None,
      wrapping = False,
     ):

    '''
      Simple interface to cctbx.development.generate_map allowing only
      a small subset of keywords. Useful for quick generation of models, map
      coefficients, and maps

      For full functionality use cctbx.development.generate_model,
      cctbx.development.generate_map_coeffs, and
      cctbx.development.generate_map

      Summary:
      --------

      Using existing model with its crystal_symmetry, if present to
      generate map.  If no existing model, use default model from library,
      box with box_cushion around it and choose n_residues to
      include default=10).

      Parameters:
      -----------

      model (model.manager object, None):    model to use (as is)
      file_name (path , None):    file containing coordinates to use (instead
                                  of default model)
      n_residues (int, 10):      Number of residues to include (from default
                                  model or file_name)
      b_iso (float, 30):         B-value (ADP) to use for all atoms
      box_cushion (float, 5):     Buffer (A) around model
      d_min (float, 3):      high_resolution limit (A)
      gridding (tuple (nx, ny, nz), None):  Gridding of map (optional)
      origin_shift_grid_units (tuple (ix, iy, iz), None):  Move location of
          origin of resulting map to (ix, iy, iz) before writing out
      wrapping:  Defines if map is to be specified as wrapped
      scattering_table (choice, 'electron'): choice of scattering table
           All choices: wk1995 it1992 n_gaussian neutron electron
      fractional_error:  resolution-dependent fractional error, ranging from
           zero at low resolution to fractional_error at d_min. Can
           be more than 1.
    '''


    if not d_min:
      d_min = self.resolution()
    if not d_min:
      d_min = 3  # default


    self._print("\nGenerating new map data\n")
    if self.map_manager():
      self._print("NOTE: replacing existing map data\n")
    if self.model() and (not model):
      self._print("NOTE: using existing model to generate map data\n")
      model = self.model()

    from cctbx.development.create_models_or_maps import generate_model, \
       generate_map_coefficients
    from cctbx.development.create_models_or_maps import generate_map \
        as generate_map_data

    if not model:
      model = generate_model(
        file_name = file_name,
        n_residues = n_residues,
        b_iso = b_iso,
        box_cushion = box_cushion,
        space_group_number = 1,
        log = self.log)
    map_coeffs = generate_map_coefficients(model = model,
        d_min = d_min,
        scattering_table = scattering_table,
        log = self.log)

    mm = generate_map_data(
      map_coeffs = map_coeffs,
      d_min = d_min,
      gridding = gridding,
      wrapping = wrapping,
      origin_shift_grid_units = origin_shift_grid_units,
      high_resolution_real_space_noise_fraction = fractional_error,
      log = self.log)

    mm.show_summary()
    self.set_up_map_dict(map_manager=mm)
    self.set_up_model_dict(model=model)

    # Set the resolution now if not already set
    if d_min and not self.resolution():
      self.set_resolution(d_min)

  def _empty_copy(self):
    '''
      Return a copy with no data
    '''
    new_mmm = map_model_manager()
    new_mmm._map_dict={}
    new_mmm._model_dict={}
    return new_mmm

  def deep_copy(self):
    '''
      Return a deep_copy of this map_manager
      Use customized copy with default map_dict and model_dict (from self)
    '''
    return self.customized_copy(map_dict = None, model_dict = None)

  def customized_copy(self, model_dict = None, map_dict = None):
    '''
      Produce a copy of this map_model object, replacing nothing,
      maps or models, or both
    '''

    # Decide what is new

    if model_dict: # take new model_dict without deep_copy
      new_model_dict = model_dict
    else:  # deep_copy existing model_dict
      new_model_dict = {}
      for id in self.model_id_list():
        new_model_dict[id]=self.get_model_by_id(id).deep_copy()

    if map_dict: # take new map_dict without deep_copy
      new_map_dict = map_dict
    else:  # deep_copy existing map_dict
      new_map_dict = {}
      for id in self.map_id_list():
        new_map_dict[id]=self.get_map_manager_by_id(id).deep_copy()

    # Build new map_manager

    new_mmm = map_model_manager()

    new_mmm._model_dict = new_model_dict
    new_mmm._map_dict = new_map_dict

    new_mmm._force_wrapping = deepcopy(self._force_wrapping)
    new_mmm._warning_message = self._warning_message

    new_mmm.set_log(self.log)
    return new_mmm


  def model_building(self,
     nproc = 1,
     soft_zero_boundary_mask = True,
     soft_zero_boundary_mask_radius = None,
     ):
    '''
     Return this object as a local_model_building object
     The model-building object has pointers to model and map_manager, not
       copies
      resolution is resolution for Fourier coefficients
      is_xray_map is True for x-ray map
      nproc is number of processors to use
    '''

    resolution = self.resolution()
    assert resolution is not None

    from phenix.model_building import local_model_building
    return local_model_building(
     map_model_manager = self, # map_model manager
     soft_zero_boundary_mask = soft_zero_boundary_mask,
     soft_zero_boundary_mask_radius = soft_zero_boundary_mask_radius,
     nproc= nproc,
    )

  def as_map_model_manager(self):
    '''
      Return this object (allows using .as_map_model_manager() on both
      map_model_manager objects and others including box.around_model() etc.

    '''
    return self


  def as_match_map_model_ncs(self):
    '''
      Return this object as a match_map_model_ncs

      Includes only the map_manager and model and ncs object, ignores all
      other maps and models (match_map_model_ncs takes only one of each).

    '''
    from iotbx.map_model_manager import match_map_model_ncs
    mmmn = match_map_model_ncs()
    if self.map_manager():
      mmmn.add_map_manager(self.map_manager())
    if self.model():
      mmmn.add_model(self.model())
    if self.ncs_object():
      mmmn.add_ncs_object(self.ncs_object())
    return mmmn


class match_map_model_ncs(object):
  '''
   match_map_model_ncs

   Use: Container to hold map, model, ncs object and check
   consistency and shift origin

   Normal usage:

     Initialize empty, then read in or add a group of model.manager,
     map_manager, and ncs objects

     Read in the models, maps, ncs objects

     Shift origin to (0, 0, 0) and save position of this (0, 0, 0) point in the
        original coordinate system so that everything can be written out
        superimposed on the original locations. This is origin_shift_grid_units
        in grid units


     NOTE: modifies the model, map_manager, and ncs objects. Call with
     deep_copy() of these if originals need to be preserved.

     Input models, maps, and ncs_object must all match in crystal_symmetry,
     original (unit_cell) crystal_symmetry, and shift_cart for maps)

     If map_manager contains an ncs_object and an ncs_object is supplied,
     the map_manager receives the supplied ncs_object

     absolute_angle_tolerance and absolute_length_tolerance are tolerances
     for crystal_symmetry.is_similar_symmetry()
  '''

  def __init__(self, log = None,
     ignore_symmetry_conflicts = None,
     absolute_angle_tolerance = 0.01,
     absolute_length_tolerance = 0.01, ):

    # Set output stream
    self.set_log(log = log)

    self._map_manager = None
    self._model = None
    self._absolute_angle_tolerance = absolute_angle_tolerance
    self._absolute_length_tolerance = absolute_length_tolerance
    self._ignore_symmetry_conflicts = ignore_symmetry_conflicts

  # prevent pickling error in Python 3 with self.log = sys.stdout
  # unpickling is limited to restoring sys.stdout
  def __getstate__(self):
    pickle_dict = self.__dict__.copy()
    if isinstance(self.log, io.TextIOWrapper):
      pickle_dict['log'] = None
    return pickle_dict

  def __setstate__(self, pickle_dict):
    self.__dict__ = pickle_dict
    if self.log is None:
      self.log = sys.stdout

  def deep_copy(self):
    new_mmmn = match_map_model_ncs()
    if self._model:
      new_mmmn.add_model(self._model.deep_copy())
    if self._map_manager:
      new_mmmn.add_map_manager(self._map_manager.deep_copy())
    return new_mmmn

  def show_summary(self):
    self._print ("Summary of maps and models")
    if self._map_manager:
      self._print("Map summary:")
      self._map_manager.show_summary(out = self.log)
    if self._model:
      self._print("Model summary:")
      self._print("Residues: %s" %(
       self._model.get_hierarchy().overall_counts().n_residues))

    if self.ncs_object():
      self._print("NCS summary:")
      self._print("Operators: %s" %(
       self.ncs_object().max_operators()))

  def set_log(self, log = sys.stdout):
    '''
       Set output log file
    '''
    if log is None:
      self.log = null_out()
    else:
      self.log = log

  def _print(self, m):
    if (self.log is not None) and hasattr(self.log, 'closed') and (
        not self.log.closed):
      self._print(m, file = self.log)

  def write_map(self, file_name = None):
    if not self._map_manager:
      self._print ("No map to write out")
    elif not file_name:
      self._print ("Need file name to write map")
    else:
      self._map_manager.write_map(file_name = file_name)

  def write_model(self,
     file_name = None):
    if not self._model:
      self._print ("No model to write out")
    elif not file_name:
      self._print ("Need file name to write model")
    else:
      # Write out model

      f = open(file_name, 'w')
      print(self._model.model_as_pdb(), file = f)
      f.close()
      self._print("Wrote model with %s residues to %s" %(
         self._model.get_hierarchy().overall_counts().n_residues,
         file_name))

  def crystal_symmetry(self):
    # Return crystal symmetry of map, or if not present, of model
    if self._map_manager:
      return self._map_manager.crystal_symmetry()
    elif self._model:
      return self._model.crystal_symmetry()
    else:
      return None

  def unit_cell_crystal_symmetry(self):
    # Return unit_cell crystal symmetry of map
    if self._map_manager:
      return self._map_manager.unit_cell_crystal_symmetry()
    else:
      return None

  def map_manager(self):
    '''
      Return the map_manager
    '''
    return self._map_manager

  def model(self):
    return self._model

  def ncs_object(self):
    if self.map_manager():
      return self.map_manager().ncs_object()
    else:
      return None


  def add_map_manager(self, map_manager):
    # Add a map and make sure its symmetry is similar to others
    assert self._map_manager is None
    self._map_manager = map_manager
    if self.model():
      self.check_model_and_set_to_match_map_if_necessary()

  def check_model_and_set_to_match_map_if_necessary(self):
    # Map, model and ncs_object all must have same symmetry and shifts at end

    if self.map_manager() and self.model():
      # Must be compatible...then set model symmetry if not set
      ok=self.map_manager().is_compatible_model(self.model(),
        absolute_angle_tolerance = self._absolute_angle_tolerance,
        absolute_length_tolerance = self._absolute_length_tolerance,
        require_match_unit_cell_crystal_symmetry=False)
      if ok or self._ignore_symmetry_conflicts:
        self.map_manager().set_model_symmetries_and_shift_cart_to_match_map(
          self.model())  # modifies self.model() in place
      else:
          raise Sorry("Model is not similar to '%s': \n%s" %(
           self.map_manager().file_name,
            self.map_manager().warning_message())+
            "\nTry 'ignore_symmetry_conflicts=True'")


  def add_model(self, model,
        set_model_log_to_null = True):
    # Add a model and make sure its symmetry is similar to others
    assert self._model is None
    # Check that model original crystal_symmetry matches full
    #    crystal_symmetry of map
    if set_model_log_to_null:
      model.set_log(null_out())
    self._model = model
    if self.map_manager():
      self.check_model_and_set_to_match_map_if_necessary()

  def add_ncs_object(self, ncs_object):
    # Add an NCS object to map_manager, overwriting any ncs object that is there
    # Must already have a map_manager. Ncs object must match shift_cart already

    assert self.map_manager() is not None
    self.map_manager().set_ncs_object(ncs_object)
    # Check to make sure its shift_cart matches
    self.check_model_and_set_to_match_map_if_necessary()

  def read_map(self, file_name):
    # Read in a map and make sure its symmetry is similar to others
    mm = map_manager(file_name)
    self.add_map_manager(mm)

  def read_model(self, file_name):
    self._print("Reading model from %s " %(file_name))
    from iotbx.pdb import input
    inp = input(file_name = file_name)
    from mmtbx.model import manager as model_manager
    model = model_manager(model_input = inp)
    self.add_model(model)


  def read_ncs_file(self, file_name):
    # Read in an NCS file and make sure its symmetry is similar to others
    from mmtbx.ncs.ncs import ncs
    ncs_object = ncs()
    ncs_object.read_ncs(file_name = file_name, log = self.log)
    if ncs_object.max_operators()<2:
       self.ncs_object.set_unit_ncs()
    self.add_ncs_object(ncs_object)

  def set_original_origin_and_gridding(self,
      original_origin = None,
      gridding = None):
    '''
     Use map_manager to reset (redefine) the original origin and gridding
     of the map.
     You can supply just the original origin in grid units, or just the
     gridding of the full unit_cell map, or both.

     Update shift_cart for model and ncs object if present.

    '''

    assert self._map_manager is not None

    self._map_manager.set_original_origin_and_gridding(
         original_origin = original_origin,
         gridding = gridding)

    # Get the current origin shift based on this new original origin
    shift_cart = self._map_manager.shift_cart()
    if self._model:
      if self._model.shift_cart() is None:
        self._model.set_unit_cell_crystal_symmetry_and_shift_cart(
          unit_cell_crystal_symmetry = \
           self._map_manager.unit_cell_crystal_symmetry())
      self._model.set_shift_cart(shift_cart)

  def shift_origin(self, desired_origin = (0, 0, 0)):
    # NOTE: desired_origin means the origin we want to achieve, not the
    #   current origin

    # shift the origin of all maps/models to desired_origin (usually (0, 0, 0))
    desired_origin = tuple(desired_origin)
    if not self._map_manager:
      self._print ("No information about origin available")
      return
    if self._map_manager.map_data().origin() == desired_origin:
      self._print("Origin is already at %s, no shifts will be applied" %(
       str(desired_origin)))

    # Figure out shift of model if incoming map and model already had a shift

    if self._model:

      # Figure out shift for model and make sure model and map agree
      shift_info = self._map_manager._get_shift_info(
         desired_origin = desired_origin)

      current_shift_cart = self._map_manager.grid_units_to_cart(
       tuple([-x for x in shift_info.current_origin_shift_grid_units]))
      expected_model_shift_cart = current_shift_cart

      shift_to_apply_cart = self._map_manager.grid_units_to_cart(
        shift_info.shift_to_apply)
      new_shift_cart = self._map_manager.grid_units_to_cart(
        tuple([-x for x in shift_info.new_origin_shift_grid_units]))
      new_full_shift_cart = new_shift_cart
      # shift_to_apply_cart is coordinate shift we are going to apply
      #  new_shift_cart is how to get to new location from original
      #   current_shift_cart is how to get to current location from original
      assert approx_equal(shift_to_apply_cart, [(a-b) for a, b in zip(
        new_shift_cart, current_shift_cart)])

      # Get shifts already applied to  model
      #    and check that they match map

      if self._model:
        existing_shift_cart = self._model.shift_cart()
        if existing_shift_cart is not None:
          assert approx_equal(existing_shift_cart, expected_model_shift_cart)

      if self._map_manager.origin_is_zero() and \
         expected_model_shift_cart == (0, 0, 0):
        pass # Need to set model shift_cart below

    # Apply shift to model, map and ncs object

    # Shift origin of map_manager
    self._map_manager.shift_origin(desired_origin = desired_origin)

    # Shift origin of model  Note this sets model shift_cart
    if self._model:
      self._model = self.shift_model_to_match_working_map(
        coordinate_shift = shift_to_apply_cart,
        new_shift_cart = new_full_shift_cart,
        model = self._model)

  def shift_ncs_to_match_working_map(self, ncs_object = None, reverse = False,
    coordinate_shift = None,
    new_shift_cart = None):

    '''
       Shift an ncs object to match the working map (based
       on self._map_manager.origin_shift_grid_units)

       The working map is the current map in its current location. Typically
       origin is at (0,0,0).

       This shifts an ncs object (typically is in its original location) to
       match this working map.

       If the ncs object was already shifted (as reflected in its shift_cart())
       it will receive the appropriate additional shift to match current map.

       If coordinate_shift is specified, it is the target final coordinate shift
       instead of the shift_cart() for the working map.
    '''

    if coordinate_shift is None:
      coordinate_shift = self.get_coordinate_shift(reverse = reverse)

    # Determine if ncs_object is already shifted
    existing_shift = ncs_object.shift_cart()

    coordinate_shift = tuple(
        [cs - es for cs, es in zip(coordinate_shift, existing_shift)])

    ncs_object = ncs_object.coordinate_offset(coordinate_shift)
    return ncs_object

  def shift_ncs_to_match_original_map(self, ncs_object = None):
    return self.shift_ncs_to_match_working_map(ncs_object = ncs_object,
      reverse = True)

  def get_coordinate_shift(self, reverse = False):
    if reverse:
       return tuple([-x for x in self._map_manager.shift_cart()])
    else:
       return self._map_manager.shift_cart()

  def shift_model_to_match_working_map(self, model = None, reverse = False,
     coordinate_shift = None,
     new_shift_cart = None):

    '''
    Shift a model based on the coordinate shift for the working map.

    Optionally specify the shift to apply (coordinate shift) and the
    new value of the shift recorded in the model (new_shift_cart)
    '''

    if coordinate_shift is None:
      coordinate_shift = self.get_coordinate_shift(
       reverse = reverse)
    if new_shift_cart is None:
      new_shift_cart = coordinate_shift

    model.shift_model_and_set_crystal_symmetry(shift_cart = coordinate_shift,
      crystal_symmetry = model.crystal_symmetry())  # keep crystal_symmetry

    # Allow specifying the final shift_cart:
    if tuple(new_shift_cart) !=  tuple(coordinate_shift):
      model.set_shift_cart(new_shift_cart)

    return model

  def shift_model_to_match_original_map(self, model = None):
    # Shift a model object to match the original map (based
    #    on -self._map_manager.origin_shift_grid_units)
    return self.shift_model_to_match_working_map(model = model, reverse = True)

  def as_map_model_manager(self):

    '''
      Return map_model_manager object with contents of this class
      (not a deepcopy)

    '''
    from iotbx.map_model_manager import map_model_manager
    mam = map_model_manager(
        map_manager = self.map_manager(),
        model = self.model(),
        )
    return mam

#   Misc methods XXX to be removed

def get_map_histograms(data, n_slots = 20, data_1 = None, data_2 = None):
  h0, h1, h2 = None, None, None
  data_min = None
  hmhcc = None
  if(data_1 is None):
    h0 = flex.histogram(data = data.as_1d(), n_slots = n_slots)
  else:
    data_min = min(flex.min(data_1), flex.min(data_2))
    data_max = max(flex.max(data_1), flex.max(data_2))
    h0 = flex.histogram(data = data.as_1d(), n_slots = n_slots)
    h1 = flex.histogram(data = data_1.as_1d(), data_min = data_min,
      data_max = data_max, n_slots = n_slots)
    h2 = flex.histogram(data = data_2.as_1d(), data_min = data_min,
      data_max = data_max, n_slots = n_slots)
    hmhcc = flex.linear_correlation(
      x = h1.slots().as_double(),
      y = h2.slots().as_double()).coefficient()
  return group_args(h_map = h0, h_half_map_1 = h1, h_half_map_2 = h2,
    _data_min = data_min, half_map_histogram_cc = hmhcc)

def get_map_counts(map_data, crystal_symmetry = None):
  a = map_data.accessor()
  map_counts = group_args(
    origin       = a.origin(),
    last         = a.last(),
    focus        = a.focus(),
    all          = a.all(),
    min_max_mean = map_data.as_1d().min_max_mean().as_tuple(),
    d_min_corner = maptbx.d_min_corner(map_data = map_data,
      unit_cell = crystal_symmetry.unit_cell()))
  return map_counts
