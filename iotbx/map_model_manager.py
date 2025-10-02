"""
High-level manager for grouping and analyzing 3D maps and models of macromolecules. This is the class to use for most map and model operations.
"""

from __future__ import absolute_import, division, print_function

import sys, os
from libtbx.utils import Sorry, Abort
from cctbx import maptbx
from cctbx import crystal
from cctbx import uctbx
from cctbx import miller
from libtbx import group_args
from cctbx.array_family import flex
from scitbx.matrix import col
from iotbx.map_manager import map_manager as MapManager
from mmtbx.model import manager as model_manager
import mmtbx.ncs.ncs
from libtbx.utils import null_out
from libtbx.test_utils import approx_equal
from copy import deepcopy
from cctbx import adptbx
from mmtbx_tls_ext import tlso, uaniso_from_tls_one_group

# Reserved phil scope for MapModelManager
map_phil_str = '''
full_map = None
  .type = path
  .help = Input full map file
  .short_caption = Map
  .style = file_type:ccp4_map input_file
half_map = None
  .type = path
  .multiple = True
  .help = Input half map files
  .short_caption = Half map
  .style = file_type:ccp4_map input_file
'''

model_phil_str = '''
model = None
  .type = path
  .help = Input model file
  .style = file_type:pdb input_file
  .short_caption = Model
'''

map_model_phil_str = '''
map_model {
  include scope iotbx.map_model_manager.map_phil_str
  include scope iotbx.map_model_manager.model_phil_str
}
'''


class map_model_manager(object):

  '''
    Class for analyzing 3D maps and models of macromolecules.

    Core functionality is the ability to shift the origin of map(s)
    and models to (0, 0, 0) and keeping track of the shifts.

    Accessory functions allow extraction of boxed maps and models,
    masking, sharpening, calculation of map-model correlations,
    and other operations.

    Uses the model manager to hold models and the map_manager to hold maps.


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

    Note:  It is permissible to call with no map_manager, but supplying
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

    If a model has no crystal_symmetry or not unit_cell, it receives the
     crystal symmetry of the map_manager and the shift_cart of the map_manager
     The same result is obtained if ignore_symmetry_conflicts is set.

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
               shift_tol = 0.001,  # shift_cart tolerance
               log              = None,
               stop_file = None,  # if present in working directory, stop
               make_cell_slightly_different_in_abc  = False,
               name = 'map_model_manager',
               verbose = False):


    # Checks
    if extra_model_list is None: extra_model_list = []
    if extra_map_manager_list is None: extra_map_manager_list = []
    for m in [model] + extra_model_list:
      assert (m is None) or isinstance(m, model_manager)
    for mm in [map_manager, map_manager_1, map_manager_2] + extra_map_manager_list:
      assert (mm is None) or isinstance(mm, MapManager)
    assert (ncs_object is None) or isinstance(ncs_object, mmtbx.ncs.ncs.ncs)


    # Set the log stream and name
    self.set_log(log = log)
    self.set_stop_file(file_name = stop_file)
    self.set_name(name)
    self.set_verbose(verbose)

    # Initialize
    self._nproc = 1
    self._multiprocessing = 'multiprocessing'
    self._queue_run_command = None
    self._resolution = None
    self._minimum_resolution = None
    self._map_dict={}
    self._model_dict = {}
    self._force_wrapping = wrapping
    self._warning_message = None
    self._scattering_table = None
    self._info = group_args(
      group_args_type = 'Information about this map_model_manager',
      )
    self._absolute_angle_tolerance = absolute_angle_tolerance
    self._absolute_length_tolerance = absolute_length_tolerance
    self._shift_tol = shift_tol

    # If no map_manager now, do not do anything and make sure there
    #    was nothing else supplied except possibly a model

    if (not map_manager) and (not map_manager_1) and (not map_manager_2):
      assert not extra_map_manager_list
      assert not ncs_object
      if model:
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

    if make_cell_slightly_different_in_abc:
      self._make_cell_slightly_different_in_abc(any_map_manager)

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

    # if the incoming model is simply shifted relative to the maps...make them
    #  match by adjusting the model shift

    if any_map_manager:
      for m in [model] + extra_model_list:
        if not m: continue
        self.add_crystal_symmetry_to_model_if_necessary(
            m, map_manager = any_map_manager)
        self.shift_any_model_to_match(m, map_manager = any_map_manager,
         set_unit_cell_crystal_symmetry = True)

    if any_map_manager and ignore_symmetry_conflicts:
      # Take all symmetry information from
      #  any_map_manager and apply it to everything

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

      if ncs_object:
        ncs_object.set_shift_cart(any_map_manager.shift_cart())

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
      if any_map_manager and m and (not ignore_symmetry_conflicts):
        if not any_map_manager.is_similar(m,
           absolute_angle_tolerance = absolute_angle_tolerance,
           absolute_length_tolerance = absolute_length_tolerance,
          ):
          raise Sorry("Map manager '%s' is not similar to '%s': %s" %(
           m.file_name,any_map_manager.file_name,
            m.warning_message())+
            "\nTry 'ignore_symmetry_conflicts=True'")

    # Now make sure all models match symmetry using _match_map_model_ncs

    # Make a _match_map_model_ncs and check unit_cell and
    #   working crystal symmetry
    #  and shift_cart for model, map, and ncs_object (if present)
    mmmn = _match_map_model_ncs(
        absolute_angle_tolerance = absolute_angle_tolerance,
        absolute_length_tolerance = absolute_length_tolerance,
        ignore_symmetry_conflicts = ignore_symmetry_conflicts)
    mmmn.set_log(self.log)
    mmmn.add_map_manager(any_map_manager)
    for m in extra_model_list: # Check all extra models
      mmmn.add_model(m.deep_copy(),
        set_model_log_to_null = False,
        ) # keep the log
      mmmn._model = None # throw it away just wanted to check it

    if model:  # Now add model for real.
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
    # Put shifted map in the right place. It is either map_manager
    #   or map_manager_1
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

    # Shift origins of all the extra models if not already done:
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

  def _make_cell_slightly_different_in_abc(self,map_manager):
    '''
    Adjust cell parameters just slightly so that gridding is not exactly the
    same in all directions.  This will make binner give uniform results
    '''
    cs=map_manager.unit_cell_crystal_symmetry()
    uc=cs.unit_cell()
    from cctbx import uctbx
    p=list(uc.parameters())
    if p[0] == p[1]:
      p[1] += 1.e-2
    if p[0] == p[2]:
      p[2] -= 1.e-2
    uc=uctbx.unit_cell(tuple(p))
    cs=cs.customized_copy(unit_cell=uc)
    map_manager.set_unit_cell_crystal_symmetry(cs)


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
    import io
    pickle_dict = self.__dict__.copy()
    if isinstance(self.log, io.TextIOWrapper):
      pickle_dict['log'] = None
    return pickle_dict

  def __setstate__(self, pickle_dict):
    self.__dict__ = pickle_dict
    if not hasattr(self, 'log') or self.log is None:
      self.log = sys.stdout

  def __repr__(self):
    text = "\nMap_model_manager '%s': \n" %(self.name)
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

  def set_name(self, name = None):
    '''
       Set name
    '''
    self.name = name


  def set_stop_file(self, file_name = None):
    '''
      Define file name that means "STOP"
    '''
    self._stop_file = file_name

  def set_verbose(self, verbose = None):
    '''
       Set verbose
    '''
    self.verbose = verbose

  # Methods to get and set info object (any information about this object)

  def set_info(self, info):
    ''' Set the information about this object'''
    self._info = info

  def info(self):
    ''' Get the information about this object'''
    return self.get_info()

  def get_info(self, item_name = None):
    ''' Get information about an item in this object'''
    if not item_name:
      return self._info
    else:
      return self._info.get(item_name)

  def add_to_info(self, item_name = None, item = None):
    ''' Add an item to information'''
    setattr(self._info,item_name, item)

  # Methods for job control

  def check_stop_file(self):
    '''Check for a stop_file'''
    if self._stop_file and os.path.isfile(self._stop_file):
      raise Abort("Stopping as the stop_file %s is present" %(
        os.path.abspath(self._stop_file)))
  # Methods for printing

  def set_log(self, log = sys.stdout):
    '''
       Set output log file
    '''
    if log is None:
      self.log = null_out()
    else:
      self.log = log

  def _print(self, m, force = False):
    '''
      Print to log if it is present
    '''

    if (self.log is not None) and hasattr(self.log, 'closed') and (
        not self.log.closed):
      print(m, file = self.log)
    elif force:
      print(m)

  # Methods for obtaining models, map_managers, symmetry, ncs_objects

  def crystal_symmetry(self):
    ''' Get the working crystal_symmetry'''
    return self.map_manager().crystal_symmetry()

  def unit_cell_crystal_symmetry(self):
    ''' Get the unit_cell_crystal_symmetry (full or original symmetry)'''
    return self.map_manager().unit_cell_crystal_symmetry()

  def shifted(self):
    ''' Determine if the maps and models in this manager are shifted relative
    to their original positions (e.g., whether the map has been boxed).
    Return True if self.map_manager() has been shifted from its
    original origin.
    '''
    return self.map_manager().shifted()

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
    del self._model_dict[model_id]

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
        # Try to make a file name
        file_name = None
        if map_manager_1.file_name and map_manager_2.file_name:
          try:
            file_name = "_and_".join([
              os.path.splitext(os.path.split(
                 map_manager_1.file_name)[-1])[0],
              os.path.splitext(os.path.split(
                 map_manager_2.file_name)[-1])[0]]) + \
              os.path.splitext(map_manager_1.file_name)[1]
          except Exception as e:
            file_name = None
        map_manager.file_name = file_name

        self._map_dict['map_manager'] = map_manager
    if self.model() and not map_manager:  # make one based on model
      crystal_symmetry = self.model().unit_cell_crystal_symmetry()
      if not crystal_symmetry:
        crystal_symmetry = self.model().crystal_symmetry()
      if not crystal_symmetry:  # make it up
        self.model().add_crystal_symmetry_if_necessary()
        crystal_symmetry = self.model().crystal_symmetry()
      if crystal_symmetry:
        from iotbx.map_manager import dummy_map_manager
        map_manager = dummy_map_manager(crystal_symmetry)
        self._map_dict['map_manager'] = map_manager
        self.info().dummy_map_manager = True # mark it
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

  def get_ncs_from_model(self):
    '''
    Return model NCS as ncs_spec object if available
    Does not set anything. If you want to save it use:
      self.set_ncs_object(self.get_ncs_from_model())
      This will set the ncs object in the map_manager (if present)
    '''
    if not self.model():
      return None
    if not self.model().get_ncs_obj():
      self.model().search_for_ncs()
    if self.model().get_ncs_obj():
      return self.model().get_ncs_obj().get_ncs_info_as_spec()
    else:
      return None

  def get_ncs_from_map(self, use_existing = True,
      include_helical_symmetry = False,
      symmetry_center = None,
      min_ncs_cc = None,
      symmetry = None,
      ncs_object = None):

    '''
    Use existing ncs object in map if present or find ncs from map
    Sets ncs_object in self.map_manager()
    Sets self._ncs_cc which can be retrieved with self.ncs_cc()
    '''
    if (not ncs_object) and use_existing:
      ncs_object = self.ncs_object()
    ncs=self.map_manager().find_map_symmetry(
        include_helical_symmetry = include_helical_symmetry,
        symmetry_center = symmetry_center,
        min_ncs_cc = min_ncs_cc,
        symmetry = symmetry,
        ncs_object = ncs_object)
    self._ncs_cc = self.map_manager().ncs_cc()
    return self.ncs_object()

  def ncs_cc(self):
    ''' Return the NCS (symmetry) correlation of the map, if present
    '''
    if hasattr(self,'_ncs_cc'):
       return self._ncs_cc

  def set_ncs_object(self, ncs_object):
    '''
    Set the ncs object of map_manager
    '''
    if not self.map_manager():
      return
    else:
      self.map_manager().set_ncs_object(ncs_object)

  def ncs_object(self):
    '''
    Get the ncs object of map_manager
    '''
    if self.map_manager():
      return self.map_manager().ncs_object()
    else:
      return None

  def experiment_type(self):
    '''
    Return the experiment_type (xray cryo_em neutron)
    '''
    if self.map_manager():
      return self.map_manager().experiment_type()
    else:
      return None

  def scattering_table(self):
    '''Return the scattering table (type of scattering)
       electron:  cryo_em
       n_gaussian x-ray (standard)
       wk1995:    x-ray (alternative)
       it1992:    x-ray (alternative)
       neutron:   neutron scattering
    '''
    if self._scattering_table:
      return self._scattering_table
    elif self.map_manager():
      return self.map_manager().scattering_table()
    else:
      return None

  def minimum_resolution(self):
    '''
     Return d_min, normally minimum available but if set, return
     value of d_min
    '''
    if not self._minimum_resolution:
      # get it and set it and return it
      self._minimum_resolution = self.map_manager().minimum_resolution()

    return self._minimum_resolution

  def nproc(self):
    ''' Return value of nproc (number of processors to use)'''
    return self._nproc

  def resolution(self,
    use_fsc_if_no_resolution_available_and_maps_available = True,
    map_id_1 = 'map_manager_1',
    map_id_2 = 'map_manager_2',
    fsc_cutoff = 0.143,
     ):
    '''Return resolution of map.  If not already calculated,
    use map-map or map-model FSC or value from map_manager'''

    if self._resolution: # have it already
      return self._resolution

    else:  # figure out resolution
      resolution = None
      if use_fsc_if_no_resolution_available_and_maps_available and \
          self.get_map_manager_by_id(map_id_1) and \
          self.get_map_manager_by_id(map_id_2):
        fsc_info = self.map_map_fsc(  # get resolution from FSC
          map_id_1 = map_id_1,
          map_id_2 = map_id_2)
        resolution = fsc_info.d_min
        if resolution is not None:
          print("\nResolution estimated from FSC of '%s' and '%s: %.3f A " %(
           map_id_1, map_id_2, resolution), file = self.log)
        elif fsc_info.fsc.fsc.min_max_mean().min > fsc_cutoff:
          print("\nResolution estimated from minimum_resolution ",
            "\nbecause FSC of '%s' and '%s is undefined" %(
           map_id_1, map_id_2), file = self.log)
          resolution = self.minimum_resolution()
        else:
          print("\nCould not obtain resolution from FSC", file = self.log)


      if (not resolution) and self.map_manager() and (
           not self.map_manager().is_dummy_map_manager()):
        # get resolution from map_manager
        resolution = self.map_manager().resolution()
        print("\nResolution obtained from map_manager: %.3f A " %(
          resolution), file = self.log)
      if resolution:
        self.set_resolution(resolution)
        return resolution
      else:
        return None

  def set_multiprocessing(self,
      nproc = None,
      multiprocessing = None,
      queue_run_command = None):
    '''  Set multiprocessing parameters'''
    if nproc:
      self._nproc = nproc

      if nproc > 1 and (multiprocessing is None):
        multiprocessing = 'multiprocessing'

    if multiprocessing:
      assert multiprocessing in ['multiprocessing','sge','lsf','pbs',
         'condor','pbspro','slurm']
      self._multiprocessing = multiprocessing

    if queue_run_command:
      self._queue_run_command = queue_run_command

  def set_resolution(self, resolution):
    ''' Set nominal resolution. Pass along to any map managers '''
    self._resolution = resolution
    for mm in self.map_managers():
      mm.set_resolution(resolution)

  def set_minimum_resolution(self, d_min):
    ''' Set minimum resolution used in calculations'''
    self._minimum_resolution = d_min

  def set_scattering_table(self, scattering_table):
    '''
     Set nominal scattering_table. Overrides anything in map_managers
       electron:  cryo_em
       n_gaussian x-ray (standard)
       wk1995:    x-ray (alternative)
       it1992:    x-ray (alternative)
       neutron:   neutron scattering
    '''
    if scattering_table is not None:
      self._scattering_table = scattering_table

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
        mm = self.get_map_manager_by_id(id)
        assert mm is not None  # map_manager_by_id must not be None
        map_data_list.append(mm.map_data())
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

  def get_any_map_manager(self):
    '''
    Return any map manager
    '''
    keys = list(self.map_dict().keys())
    if not keys:
      # Make a map dummy manager
      mm = self.map_manager() # makes a dummy one if possible
      keys = list(self.map_dict().keys())
      if keys:
        return self.map_dict()[keys[0]]
      else:
        return None
    else:
      return self.map_dict()[keys[0]]

  def get_map_data_by_id(self, map_id):
    ''' Get map_data from a map_manager with the name map_id'''
    map_manager = self.get_map_manager_by_id(map_id)
    if map_manager and (not map_manager.is_dummy_map_manager()):
      return map_manager.map_data()
    else:
      return None

  def set_model(self,model, overwrite = True):
    '''
     Overwrites existing model with id 'model'
     Allows setting model to None if overwrite is True
    '''
    self.add_model_by_id(model,'model', overwrite = overwrite)


  def add_model_by_id(self, model, model_id,
     overwrite = True):
    '''
     Add a new model
     Must be similar to existing map_managers
     Overwrites any existing with the same id unless overwrite = False
     If model is None, removes model unless overwrite is False
    '''
    if (not model) and (not overwrite):
      print("No model supplied for '%s' ... skipping addition" %(
        model_id), file = self.log)
      return

    assert (model is None) or isinstance(model, mmtbx.model.manager)
    if not overwrite:
      assert not model_id in self.model_id_list() # must not duplicate

    if model and \
      self.map_manager() and (
         not self.map_manager().is_compatible_model(model,
          require_match_unit_cell_crystal_symmetry=True)):
      # needs shifting
      self.shift_any_model_to_match(model,
         set_unit_cell_crystal_symmetry = True)
    self._model_dict[model_id] = model

  def set_map_manager(self, map_manager):
    '''
     Overwrites existing map_manager with id 'map_manager'
    '''
    self.add_map_manager_by_id(map_manager, 'map_manager')

  def add_map_manager_by_id(self, map_manager, map_id,
     overwrite = True, force = False):
    '''
     Add a new map_manager
     Must be similar to existing
     Overwrites any existing with the same id unless overwrite = False
     Is a mask if is_mask is set
    '''
    if not map_manager and not force:
      print("No map_manager supplied for '%s' ... skipping addition" %(
        map_id), file = self.log)
      return

    assert isinstance(map_manager, MapManager)
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
    assert isinstance(map_manager, MapManager)

    self._map_dict[new_map_id] = map_manager.deep_copy()


  # Methods for writing maps and models

  def write_map(self, file_name,
       map_id='map_manager'):
    '''Write out map defined by map_id (default is 'map_manager')'''

    if not self._map_dict.get(map_id):
      self._print ("No map to write out with id='%s'" %(map_id))
    elif not file_name:
      self._print ("Need file name to write map")
    else:
      self._map_dict.get(map_id).write_map(file_name = file_name)

  def write_model(self,
     file_name,
     model_id = None,
     model = None,
     data_manager = None,
     format = None,
     ):
    ''' Write a model object specified by model_id (default is 'model')
    '''
    if not model:
      if not model_id:
        model_id = 'model'
      model = self.get_model_by_id(model_id = model_id)
    if not model:
        self._print ("No model to write out")
    elif not file_name:
      self._print ("Need file name to write model")
    else:
      # Write out model
      if not data_manager:
        from iotbx.data_manager import DataManager
        data_manager = DataManager()
        data_manager.set_overwrite(True)
      file_name = data_manager.write_model_file(model,
       file_name, format = format)
      self._print("Wrote model with %s residues to %s" %(
         model.get_hierarchy().overall_counts().n_residues,
         file_name))

  # Methods for identifying which map_manager and model to use

  def _get_map_info(self):
    '''
      Return a group_args object specifying the map_manager and
      a list of any other maps present
    '''
    all_map_id_list=list(self._map_dict.keys())
    # We are going to need id='map_manager'   create if if missing
    if self.map_manager() is None: # creates it usually but if it can't ...
      return group_args(map_id=None, other_map_id_list = [])
    if not all_map_id_list:
      return group_args(map_id=None, other_map_id_list = [])
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
    if not all_model_id_list:
       return group_args(model_id=None,
         other_model_id_list=[])
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
     upper_bounds,
     boundary_to_smoothing_ratio = 2.,
     soft_mask_around_edges = None,
     soft_mask_radius = None,
     stay_inside_current_map = True,
     use_cubic_boxing = False,
     require_match_unit_cell_crystal_symmetry = None,
     model_can_be_outside_bounds = None):
    '''
      Runs box_all_maps_with_bounds_and_shift_origin with extract_box=True
    '''
    return self.box_all_maps_with_bounds_and_shift_origin(
      lower_bounds = lower_bounds,
      upper_bounds = upper_bounds,
      model_can_be_outside_bounds = model_can_be_outside_bounds,
      soft_mask_radius = soft_mask_radius,
      soft_mask_around_edges = soft_mask_around_edges,
      boundary_to_smoothing_ratio = boundary_to_smoothing_ratio,
      stay_inside_current_map = stay_inside_current_map,
      use_cubic_boxing = use_cubic_boxing,
      require_match_unit_cell_crystal_symmetry =
        require_match_unit_cell_crystal_symmetry,
      extract_box = True)

  def box_all_maps_with_bounds_and_shift_origin(self,
     lower_bounds,
     upper_bounds,
     model_can_be_outside_bounds = None,
     soft_mask_radius = None,
     soft_mask_around_edges = None,
     boundary_to_smoothing_ratio = 2.,
     stay_inside_current_map = True,
     use_cubic_boxing = False,
     require_match_unit_cell_crystal_symmetry = False,
     extract_box = False):
    '''
       Box all maps using specified bounds, shift origin of maps, model
       Replaces existing map_managers and shifts model in place

       If extract_box=True:  Creates new object with deep_copies.
       Otherwise: replaces existing map_managers and shifts model in place

       NOTE: This changes the gridding and shift_cart of the maps and model
       Also changes space group to p1

       Can be used in map_model_manager to work with boxed maps
       and model or in map_model_manager to re-box all maps and model

       The lower_bounds and upper_bounds define the region to be boxed. These
       bounds are relative to the current map with origin at (0, 0, 0).

       if soft_mask_around_edges, still uses the same bounds, but makes
         a soft mask around the edges.  Use this option if you are going
         to calculate a FT of the map or otherwise manipulate it in
         reciprocal space. Do not use this option if you are going to
         mask around atoms, density, mask or anything
         else afterwards as you should apply a mask only once.

       If use_cubic_box, make a cubic box (in grid units). If also
        stay_inside_current_map is set, keep the cubic box inside current map

       If require_match_unit_cell_crystal_symmetry is False, do not require
       unit_cell crystal symmetry to match.

    '''
    assert lower_bounds is not None and upper_bounds is not None
    assert len(tuple(lower_bounds)) == 3
    assert len(tuple(upper_bounds)) == 3

    from cctbx.maptbx.box import with_bounds

    map_info=self._get_map_info()
    map_manager = self._map_dict[map_info.map_id]
    assert map_manager is not None

    model_info=self._get_model_info()
    model = self._model_dict.get(model_info.model_id,None)

    if extract_box and model: # make sure everything is deep_copy
      model = model.deep_copy()

    if soft_mask_around_edges: # make the cushion bigger
      pass # Here we fix the bounds based on what is requested

    # Make box with bounds and apply it to model, first map
    box = with_bounds(
      map_manager = self._map_dict[map_info.map_id],
      lower_bounds = lower_bounds,
      upper_bounds = upper_bounds,
      model = model,
      wrapping = self._force_wrapping,
      model_can_be_outside_bounds = model_can_be_outside_bounds,
      stay_inside_current_map = stay_inside_current_map,
      use_cubic_boxing = use_cubic_boxing,
      require_match_unit_cell_crystal_symmetry =
         require_match_unit_cell_crystal_symmetry,
      log = self.log)
    # Now box is a copy of map_manager and model that is boxed

    # Now apply boxing to other maps and models and then insert them into
    #  either this map_model_manager object, replacing what is there
    #  (extract_box=False)
    #  or create and return a new map_model_manager object (extract_box=True)
    return self._finish_boxing(box = box, model_info = model_info,
      map_info = map_info,
      soft_mask_radius = soft_mask_radius,
      soft_mask_around_edges = soft_mask_around_edges,
      boundary_to_smoothing_ratio = boundary_to_smoothing_ratio,
      extract_box = extract_box)

  def extract_all_maps_around_model(self,
     selection_string = None,
     selection = None,
     select_unique_by_ncs = False,
     model_can_be_outside_bounds = None,
     stay_inside_current_map = None,
     box_cushion = 5.,
     boundary_to_smoothing_ratio = 2.,
     soft_mask_around_edges = None,
     soft_mask_radius = None,
     require_match_unit_cell_crystal_symmetry = None,
     use_cubic_boxing = False,
     ):
    '''
      Runs box_all_maps_around_model_and_shift_origin with extract_box=True
    '''
    return self.box_all_maps_around_model_and_shift_origin(
      selection_string = selection_string,
      selection = selection,
      box_cushion = box_cushion,
      select_unique_by_ncs = select_unique_by_ncs,
      model_can_be_outside_bounds = model_can_be_outside_bounds,
      stay_inside_current_map = stay_inside_current_map,
      soft_mask_radius = soft_mask_radius,
      soft_mask_around_edges = soft_mask_around_edges,
      boundary_to_smoothing_ratio = boundary_to_smoothing_ratio,
      require_match_unit_cell_crystal_symmetry =
        require_match_unit_cell_crystal_symmetry,
      use_cubic_boxing = use_cubic_boxing,
      extract_box = True)

  def box_all_maps_around_model_and_shift_origin(self,
     selection_string = None,
     selection = None,
     box_cushion = 5.,
     select_unique_by_ncs = False,
     model_can_be_outside_bounds = None,
     stay_inside_current_map = None,
     soft_mask_radius = None,
     soft_mask_around_edges = None,
     boundary_to_smoothing_ratio = 2.,
     use_cubic_boxing = False,
     require_match_unit_cell_crystal_symmetry = None,
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
       If selection is specified, use instead of selection_string

       If select_unique_by_ncs is set, select the unique part of the model
       automatically.  Any selection in selection_string or selection
        will not be applied.

      if soft_mask_around_edges, makes a bigger box and makes a soft mask around
       the edges.  Use this option if you are going to calculate a FT of
       the map or otherwise manipulate it in reciprocal space. Do not use this
       option if you are going to mask around atoms, density, mask or anything
       else afterwards as you should apply a mask only once.

       If use_cubic_box, make a cubic box (in grid units). If also
        stay_inside_current_map is set, keep the cubic box inside current map

       If require_match_unit_cell_crystal_symmetry is False, do not require
       unit_cell crystal symmetry to match.
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
    elif selection:
      model = model.select(selection)
    elif extract_box and model: # make sure everything is deep_copy
      model = model.deep_copy()


    if soft_mask_around_edges: # make the cushion bigger
      box_cushion += boundary_to_smoothing_ratio * self.resolution()

    # Make box around model and apply it to model, first map
    # This step modifies model in place and creates a new map_manager
    box = around_model(
      map_manager = self._map_dict[map_info.map_id],
      model = model,
      box_cushion = box_cushion,
      wrapping = self._force_wrapping,
      model_can_be_outside_bounds = model_can_be_outside_bounds,
      stay_inside_current_map = stay_inside_current_map,
      use_cubic_boxing = use_cubic_boxing,
      require_match_unit_cell_crystal_symmetry =
         require_match_unit_cell_crystal_symmetry,
      log = self.log)
    # Now box is a copy of map_manager and model that is boxed


    # Now apply boxing to other maps and models and then insert them into
    #  either this map_model_manager object, replacing what
    #  is there (extract_box=False)
    #  or create and return a new map_model_manager object (extract_box=True)
    return self._finish_boxing(box = box, model_info = model_info,
      map_info = map_info,
      soft_mask_radius = soft_mask_radius,
      soft_mask_around_edges = soft_mask_around_edges,
      boundary_to_smoothing_ratio = boundary_to_smoothing_ratio,
      extract_box = extract_box)

  def extract_all_maps_around_density(self,
     box_cushion = 5.,
     threshold = 0.05,
     get_half_height_width = True,
     model_can_be_outside_bounds = None,
     boundary_to_smoothing_ratio = 2.,
     soft_mask_around_edges = None,
     soft_mask_radius = None,
     stay_inside_current_map = True,
     require_match_unit_cell_crystal_symmetry = None,
     use_cubic_boxing = False,
     map_id = 'map_manager'):
    '''
      Runs box_all_maps_around_density_and_shift_origin with extract_box=True
    '''
    return self.box_all_maps_around_density_and_shift_origin(
      box_cushion = box_cushion,
      threshold = threshold,
      get_half_height_width = get_half_height_width,
      model_can_be_outside_bounds = model_can_be_outside_bounds,
      map_id = map_id,
      soft_mask_radius = soft_mask_radius,
      soft_mask_around_edges = soft_mask_around_edges,
      boundary_to_smoothing_ratio = boundary_to_smoothing_ratio,
      stay_inside_current_map = stay_inside_current_map,
      require_match_unit_cell_crystal_symmetry =
        require_match_unit_cell_crystal_symmetry,
      use_cubic_boxing = use_cubic_boxing,
      extract_box = True)

  def box_all_maps_around_density_and_shift_origin(self,
     box_cushion = 5.,
     threshold = 0.05,
     map_id = 'map_manager',
     get_half_height_width = True,
     model_can_be_outside_bounds = None,
     soft_mask_radius = None,
     soft_mask_around_edges = None,
     boundary_to_smoothing_ratio = 2.,
     stay_inside_current_map = None,
     use_cubic_boxing = False,
     require_match_unit_cell_crystal_symmetry = False,
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

      if soft_mask_around_edges, makes a bigger box and makes a soft mask around
       the edges.  Use this option if you are going to calculate a FT of
       the map or otherwise manipulate it in reciprocal space.

       If use_cubic_box, make a cubic box (in grid units). If also
        stay_inside_current_map is set, keep the cubic box inside current map

       If require_match_unit_cell_crystal_symmetry is False, do not require
       unit_cell crystal symmetry to match.
    '''
    assert box_cushion is not None

    from cctbx.maptbx.box import around_density

    map_info=self._get_map_info()
    assert map_info.map_id is not None
    model_info=self._get_model_info()
    model = self._model_dict[model_info.model_id]
    if extract_box and model: # make sure everything is deep_copy
      model = model.deep_copy()

    if soft_mask_around_edges: # make the cushion bigger
      box_cushion += boundary_to_smoothing_ratio * self.resolution()

    # Make box around model and apply it to model, first map
    box = around_density(
      map_manager = self._map_dict[map_info.map_id],
      model       = model,
      box_cushion = box_cushion,
      threshold   = threshold,
      get_half_height_width = get_half_height_width,
      model_can_be_outside_bounds = model_can_be_outside_bounds,
      stay_inside_current_map = stay_inside_current_map,
      require_match_unit_cell_crystal_symmetry =
         require_match_unit_cell_crystal_symmetry,
      use_cubic_boxing = use_cubic_boxing,
      wrapping    = self._force_wrapping)

    # Now box is a copy of map_manager and model that is boxed

    # Now apply boxing to other maps and models and then insert them into
    #  either this map_model_manager object, replacing what is there (extract_box=False)
    #  or create and return a new map_model_manager object (extract_box=True)
    return self._finish_boxing(box = box, model_info = model_info,
      map_info = map_info,
      soft_mask_radius = soft_mask_radius,
      soft_mask_around_edges = soft_mask_around_edges,
      boundary_to_smoothing_ratio = boundary_to_smoothing_ratio,
      extract_box = extract_box)

  def extract_all_maps_around_mask(self,
     box_cushion = 5.,
     model_can_be_outside_bounds = None,
     boundary_to_smoothing_ratio = 2.,
     soft_mask_around_edges = None,
     soft_mask_radius = None,
     stay_inside_current_map = True,
     require_match_unit_cell_crystal_symmetry = None,
     use_cubic_boxing = False,
     mask_id = 'mask'):
    '''
      Runs box_all_maps_around_mask_and_shift_origin with extract_box=True
    '''
    return self.box_all_maps_around_mask_and_shift_origin(
      box_cushion = 5.,
      mask_id = mask_id,
      model_can_be_outside_bounds = model_can_be_outside_bounds,
      soft_mask_radius = soft_mask_radius,
      soft_mask_around_edges = soft_mask_around_edges,
      boundary_to_smoothing_ratio = boundary_to_smoothing_ratio,
      stay_inside_current_map = stay_inside_current_map,
      require_match_unit_cell_crystal_symmetry =
        require_match_unit_cell_crystal_symmetry,
      use_cubic_boxing = use_cubic_boxing,
      extract_box = True)

  def box_all_maps_around_mask_and_shift_origin(self,
     box_cushion = 5.,
     mask_id = 'mask',
     model_can_be_outside_bounds = None,
     soft_mask_radius = None,
     soft_mask_around_edges = None,
     boundary_to_smoothing_ratio = 2.,
     stay_inside_current_map = True,
     use_cubic_boxing = False,
     require_match_unit_cell_crystal_symmetry = False,
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

      if soft_mask_around_edges, makes a bigger box and makes a soft mask around
       the edges.  Use this option if you are going to calculate a FT of
       the map or otherwise manipulate it in reciprocal space.

       If use_cubic_box, make a cubic box (in grid units). If also
        stay_inside_current_map is set, keep the cubic box inside current map

       If require_match_unit_cell_crystal_symmetry is False, do not require
       unit_cell crystal symmetry to match.
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
    if extract_box and model: # make sure everything is deep_copy
      model = model.deep_copy()

    if soft_mask_around_edges: # make the cushion bigger
      box_cushion += boundary_to_smoothing_ratio * self.resolution()

    # Make box around mask and apply it to model, first map
    box = around_mask(
      map_manager = map_manager,
      mask_as_map_manager = mask_mm,
      model = model,
      box_cushion = box_cushion,
      model_can_be_outside_bounds = model_can_be_outside_bounds,
      stay_inside_current_map = stay_inside_current_map,
      use_cubic_boxing = use_cubic_boxing,
      require_match_unit_cell_crystal_symmetry =
         require_match_unit_cell_crystal_symmetry,
      wrapping = self._force_wrapping,
      log = self.log)
    # Now box is a copy of map_manager and model that is boxed

    # Now apply boxing to other maps and models and then insert them into
    #  either this map_model_manager object, replacing what is there (extract_box=False)
    #  or create and return a new map_model_manager object (extract_box=True)
    return self._finish_boxing(box = box, model_info = model_info,
      map_info = map_info,
      soft_mask_radius = soft_mask_radius,
      soft_mask_around_edges = soft_mask_around_edges,
      boundary_to_smoothing_ratio = boundary_to_smoothing_ratio,
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
     require_match_unit_cell_crystal_symmetry = None,
     keep_low_density = True,
     symmetry = None,
     boundary_to_smoothing_ratio = 2.,
     soft_mask_around_edges = None,
     keep_this_region_only = None,
     residues_per_region = None,
     soft_mask_radius = None,
     mask_expand_ratio = 1,
     stay_inside_current_map = True,
     use_cubic_boxing = False,
     use_symmetry_in_extract_unique = True):

    '''
      Runs box_all_maps_around_unique_and_shift_origin with extract_box=True
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
      require_match_unit_cell_crystal_symmetry =
        require_match_unit_cell_crystal_symmetry,
      keep_low_density = keep_low_density,
      keep_this_region_only = keep_this_region_only,
      residues_per_region = residues_per_region,
      symmetry = symmetry,
      mask_expand_ratio = mask_expand_ratio,
      soft_mask_radius = soft_mask_radius,
      soft_mask_around_edges = soft_mask_around_edges,
      boundary_to_smoothing_ratio = boundary_to_smoothing_ratio,
      use_symmetry_in_extract_unique = use_symmetry_in_extract_unique,
      stay_inside_current_map = stay_inside_current_map,
      use_cubic_boxing = use_cubic_boxing,
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
     soft_mask_radius = None,
     soft_mask_around_edges = None,
     boundary_to_smoothing_ratio = 2.,
     keep_this_region_only = None,
     residues_per_region = None,
     use_symmetry_in_extract_unique = True,
     stay_inside_current_map = True,
     use_cubic_boxing = False,
     extract_box = False,
     require_match_unit_cell_crystal_symmetry = False):
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

       Use_symmetry can be set to False to ignore symmetry found in ncs_object.
         NCS object is still kept and shifted however.

      if soft_mask_around_edges, makes a bigger box and makes a soft mask around
       the edges.  Use this option if you are going to calculate a FT of
       the map or otherwise manipulate it in reciprocal space.

      If use_cubic_box, make a cubic box (in grid units). If also
        stay_inside_current_map is set, keep the cubic box inside current map

       Additional parameters:
         mask_expand_ratio:   allows increasing masking radius beyond default at
                              final stage of masking
         solvent_content:  fraction of cell not occupied by macromolecule. Can
                            be None in which case it is estimated from map
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
         residues_per_region:  Try to segment with this many residues per region
         keep_this_region_only:  Keep just this region (first one is 0 not 1)
         require_match_unit_cell_crystal_symmetry: require model unit_cell_
           crystal_symmetry to match when extracting.

    '''
    from cctbx.maptbx.box import around_unique

    map_info=self._get_map_info()
    map_manager = self._map_dict[map_info.map_id]
    assert isinstance(map_manager, MapManager)
    if not resolution:
      resolution = self.resolution()
    assert resolution is not None

    model_info=self._get_model_info()
    model = self._model_dict.get(model_info.model_id,None)
    if extract_box and model: # make sure everything is deep_copy
      model = model.deep_copy()

    if soft_mask_around_edges: # make the cushion bigger
      box_cushion += boundary_to_smoothing_ratio * self.resolution()

    # Make box with around_unique and apply it to model, first map
    box = around_unique(
      map_manager = map_manager,
      model = model,
      wrapping = self._force_wrapping,
      target_ncs_au_model = target_ncs_au_model,
      regions_to_keep = regions_to_keep,
      residues_per_region = residues_per_region,
      keep_this_region_only = keep_this_region_only,
      solvent_content = solvent_content,
      resolution = resolution,
      sequence = sequence,
      molecular_mass = molecular_mass,
      symmetry = symmetry,
      use_symmetry_in_extract_unique = use_symmetry_in_extract_unique,
      chain_type = chain_type,
      box_cushion = box_cushion,
      soft_mask = soft_mask,
      mask_expand_ratio = mask_expand_ratio,
      stay_inside_current_map = stay_inside_current_map,
      use_cubic_boxing = use_cubic_boxing,
      require_match_unit_cell_crystal_symmetry =
          require_match_unit_cell_crystal_symmetry,
      log = self.log)

    info = box.info()
    if info and hasattr(info, 'available_selected_regions'):
      self.set_info(info)  # save this information

    # Now box is a copy of map_manager and model that is boxed

    # Now apply boxing to other maps and models and then insert them into
    #  either this map_model_manager object, replacing what is there (extract_box=False)
    #  or create and return a new map_model_manager object (extract_box=True)
    other = self._finish_boxing(box = box, model_info = model_info,
      map_info = map_info,
      soft_mask_radius = soft_mask_radius,
      soft_mask_around_edges = soft_mask_around_edges,
      boundary_to_smoothing_ratio = boundary_to_smoothing_ratio,
      extract_box = extract_box)

    if not extract_box:
      other = self #  modifying this object

    # Now apply masking to all other maps (not done in _finish_boxing)
    for id in map_info.other_map_id_list:
      box.apply_around_unique_mask(
        other._map_dict[id],
        resolution = resolution,
        soft_mask = soft_mask)

    if extract_box:
      return other

  def _finish_boxing(self, box, model_info, map_info,
    soft_mask_radius = None,
    soft_mask_around_edges = None,
    boundary_to_smoothing_ratio = None,
    extract_box = False):

    '''
       Finish copying information to boxed map_model_manager

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
      other._model_dict[id] = box.apply_to_model(
          self._model_dict[id].deep_copy(),
         )

    # Copy default information over
    name = '%s_boxed' %(self.name)
    self._set_default_parameters(other, name = name)

    if soft_mask_around_edges:
      other.mask_all_maps_around_edges(
        soft_mask_radius = soft_mask_radius,
        boundary_to_smoothing_ratio = boundary_to_smoothing_ratio)

    if extract_box:
      return other

  def merge_split_maps_and_models(self,
      model_id = None,
      box_info = None,
      replace_coordinates = True,
      replace_u_aniso = False,
      allow_changes_in_hierarchy = False,
      output_model_id = None):
    '''
      Replaces coordinates in working model with those from the
        map_model_managers in box_info.  The box_info object should
        come from running split_up_map_and_model in this instance
        of the map_model_manager.

      If allow_changes_in_hierarchy is set, create a new working model where
      the hierarchy has N "models", one from each box. This allows changing
      the hierarchy structure. This will create one model in a new hierarchy
      for each box, numbered by box number.  The new hierarchy will
      be placed in a model with id  output_model_id (default is
      model_id, replacing existing model specified by model_id; usually
      this is just 'model', the default model.)
    '''
    if model_id is None:
      model_id = 'model'
    if allow_changes_in_hierarchy and output_model_id is None:
      output_model_id = 'model'

    if allow_changes_in_hierarchy:
      print(
        "\nModels from %s boxed models will be 'models' in new hierarchy" %(
        len(box_info.selection_list)), file = self.log)
      print("New model id will be: %s" %(output_model_id),
          file = self.log)
      if output_model_id in self.model_id_list():
        print("NOTE: Replacing model %s with new composite model" %(
         output_model_id), file = self.log)

      # Set up a new empty model and hierarchy
      import iotbx.pdb
      pdb_inp = iotbx.pdb.input(source_info='text', lines=[""])
      ph = pdb_inp.construct_hierarchy()
      # Make a new model and save it as output_model_id
      self.model_from_hierarchy(ph,
        model_id = output_model_id)
      # Get this hierarchy so we can add models to it:
      working_model = self.get_model_by_id(
        model_id = output_model_id)

    else:
      print(
        "\nMerging coordinates from %s boxed models into working model" %(
        len(box_info.selection_list)), file = self.log)
      print("Working model id is : %s" %(model_id),
          file = self.log)

    i = 0
    if not hasattr(box_info,'tlso_list'):
       box_info.tlso_list = len(box_info.mmm_list) * [None]

    for selection, mmm, tlso_value in zip (
      box_info.selection_list, box_info.mmm_list, box_info.tlso_list):
      i += 1
      model_to_merge = self.get_model_from_other(mmm,
        other_model_id=model_id)

      if allow_changes_in_hierarchy:  # Add a model to the hierarchy
        ph_models_in_box = 0
        for m in model_to_merge.get_hierarchy().models():
          ph_models_in_box += 1
          assert ph_models_in_box <= 1  # cannot have multiple models in box
          mm = m.detached_copy()
          mm.id = "%s" %(i) # model number is box number as string
          working_model.get_hierarchy().append_model(mm)
      else:  # replace sites and/or u_aniso values in existing model
        if replace_coordinates: # all sites
          sites_cart = self.get_model_by_id(model_id).get_sites_cart()
          #  Sites to merge from this model
          new_coords=model_to_merge.get_sites_cart()
          original_coords=sites_cart.select(selection)
          rmsd=new_coords.rms_difference(original_coords)
          print("RMSD for %s coordinates in model %s: %.3f A" %(
             original_coords.size(), i, rmsd), file = self.log)
          sites_cart.set_selected(selection, new_coords)
          self.get_model_by_id(model_id).set_crystal_symmetry_and_sites_cart(
            sites_cart = sites_cart,
            crystal_symmetry = self.get_model_by_id(
              model_id).crystal_symmetry())

        if replace_u_aniso and tlso_value: # calculate aniso U from
          print("Replacing u_cart values based on TLS info",file = self.log)
          xrs=self.get_model_by_id(model_id).get_xray_structure()
          xrs.convert_to_anisotropic()
          uc = xrs.unit_cell()
          sites_cart = xrs.sites_cart()
          u_cart=xrs.scatterers().extract_u_cart(uc)
          new_anisos= uaniso_from_tls_one_group(tlso = tlso_value,
           sites_cart = sites_cart.select(selection),
           zeroize_trace=False)
          u_cart.set_selected(selection, new_anisos)
          xrs.set_u_cart(u_cart)
          self.get_model_by_id(model_id).set_xray_structure(xrs)

    if allow_changes_in_hierarchy:
      working_model.reset_after_changing_hierarchy() # REQUIRED

  def split_up_map_and_model_by_chain(self,
    model_id = 'model',
    skip_waters = False,
    skip_hetero = False,
    box_cushion = 3,
    mask_around_unselected_atoms = None,
    mask_radius = 3,
    masked_value = -10,
    write_files = False,
    apply_box_info = True,
     ):
    '''
     Split up the map, boxing around each chain in the model.

       Returns a group_args object containing list of the map_model_manager
         objects and a list of the selection objects that define which atoms
         from the working model are in each object.

       Normally do work on each map_model_manager to create a new model with
         the same atoms, then use merge_split_maps_and_models() to replace
         coordinates in the original model with those from all the component
         models.
       Optionally carry out the step box_info = get_split_maps_and_models(...)
         separately with the keyword apply_box_info=False


       skip_waters and skip_hetero define whether waters and hetero atoms are
        ignored
       box_cushion is the padding around the model atoms when creating boxes

    '''

    return self._split_up_map_and_model(
      selection_method = 'by_chain',
      model_id = model_id,
      skip_waters = skip_waters,
      skip_hetero = skip_hetero,
      box_cushion = box_cushion,
      mask_around_unselected_atoms = mask_around_unselected_atoms,
      mask_radius = mask_radius,
      masked_value = masked_value,
      apply_box_info = apply_box_info,
      write_files = write_files)

  def split_up_map_and_model_by_segment(self,
    model_id = 'model',
    skip_waters = False,
    skip_hetero = False,
    box_cushion = 3,
    mask_around_unselected_atoms = None,
    mask_radius = 3,
    masked_value = -10,
    write_files = False,
    apply_box_info = True,
     ):
    '''
     Split up the map, boxing around each segment (each unbroken part of
      each chain) in the model

       Returns a group_args object containing list of the map_model_manager
         objects and a list of the selection objects that define which atoms
         from the working model are in each object.

       Normally do work on each map_model_manager to create a new model with
         the same atoms, then use merge_split_maps_and_models() to replace
         coordinates in the original model with those from all the component
         models.
       Optionally carry out the step box_info = get_split_maps_and_models(...)
         separately with the keyword apply_box_info=False

       skip_waters and skip_hetero define whether waters and hetero atoms are
        ignored
       box_cushion is the padding around the model atoms when creating boxes
    '''

    return self._split_up_map_and_model(
      selection_method = 'by_segment',
      model_id = model_id,
      skip_waters = skip_waters,
      skip_hetero = skip_hetero,
      box_cushion = box_cushion,
      mask_around_unselected_atoms = mask_around_unselected_atoms,
      mask_radius = mask_radius,
      masked_value = masked_value,
      apply_box_info = apply_box_info,
      write_files = write_files)

  def split_up_map_and_model_by_ncs_groups(self,
    model_id = 'model',
    box_cushion = 3,
    mask_around_unselected_atoms = None,
    mask_radius = 3,
    masked_value = -10,
    write_files = False,
    apply_box_info = True,
     ):
    '''
     Split up the map, boxing around atoms selected with each ncs group in
     ncs_groups obtained from supplied model


       Returns a group_args object containing list of the map_model_manager
         objects and a list of the selection objects that define which atoms
         from the working model are in each object.

       Normally do work on each map_model_manager to create a new model with
         the same atoms, then use merge_split_maps_and_models() to replace
         coordinates in the original model with those from all the component
         models.
       Optionally carry out the step box_info = get_split_maps_and_models(...)
         separately with the keyword apply_box_info=False

       box_cushion is the padding around the model atoms when creating boxes
    '''

    if model_id is None:
      model_id = 'model'
    model = self.get_model_by_id(model_id = model_id)
    if model is None:
      print("No model to work with", file = self.log)
      return None # no model to work with

    ncs_groups = model.get_ncs_groups()
    selection_list = []
    if ncs_groups is None or len(ncs_groups) < 1:
      selection_list =  [model.selection("all")]
    else:
      for g in ncs_groups:
        selection_list.append(g.master_iselection)

    return self._split_up_map_and_model(
      selection_method = 'supplied_selections',
      model_id = model_id,
      selection_list = selection_list,
      box_cushion = box_cushion,
      mask_around_unselected_atoms = mask_around_unselected_atoms,
      mask_radius = mask_radius,
      masked_value = masked_value,
      apply_box_info = apply_box_info,
      write_files = write_files)
  def split_up_map_and_model_by_supplied_selections(self,
    selection_list,
    model_id = 'model',
    box_cushion = 3,
    mask_around_unselected_atoms = None,
    mask_radius = 3,
    masked_value = -10,
    write_files = False,
    apply_box_info = True,
     ):
    '''
     Split up the map, boxing around atoms selected with each selection in
      selection_list
      Note: a selection can be obtained with:
        self.model().selection(selection_string)

       Returns a group_args object containing list of the map_model_manager
         objects and a list of the selection objects that define which atoms
         from the working model are in each object.

       Normally do work on each map_model_manager to create a new model with
         the same atoms, then use merge_split_maps_and_models() to replace
         coordinates in the original model with those from all the component
         models.
       Optionally carry out the step box_info = get_split_maps_and_models(...)
         separately with the keyword apply_box_info=False

       box_cushion is the padding around the model atoms when creating boxes
    '''

    return self._split_up_map_and_model(
      selection_method = 'supplied_selections',
      model_id = model_id,
      selection_list = selection_list,
      box_cushion = box_cushion,
      mask_around_unselected_atoms = mask_around_unselected_atoms,
      mask_radius = mask_radius,
      masked_value = masked_value,
      apply_box_info = apply_box_info,
      write_files = write_files)

  def split_up_map_and_model_by_boxes(self,
    model_id = 'model',
    skip_waters = False,
    skip_hetero = False,
    write_files = False,
    target_for_boxes = 24,
    select_final_boxes_based_on_model = True,
    box_cushion = 3,
    mask_around_unselected_atoms = None,
    mask_radius = 3,
    masked_value = -10,
    skip_empty_boxes = True,
    apply_box_info = True,
    mask_id = None,
    exclude_points_outside_density = None,
    minimum_boxes_inside_density = None,
     ):
    '''
     Split up the map, creating boxes that time the entire map.

     Try to get about target_for_boxes boxes. Do not go over this target

     If select_final_boxes_based_on_model then make the final boxes just go
       around the selected parts of the model with cushion defined by
       box_cushion and not tile the map.
     Otherwise select atoms inside the boxes and afterwards expand the boxes
       with box_cushion

     If skip_empty_boxes then skip boxes with no model.

     Note that this procedure just selects by atom so you can get a single atom
      in a box

       Returns a group_args object containing list of the map_model_manager
         objects and a list of the selection objects that define which atoms
         from the working model are in each object.

       Normally do work on each map_model_manager to create a new model with
         the same atoms, then use merge_split_maps_and_models() to replace
         coordinates in the original model with those from all the component
         models.
       Optionally carry out the step box_info = get_split_maps_and_models(...)
         separately with the keyword apply_box_info=False

       skip_waters and skip_hetero define whether waters and hetero atoms are
        ignored

    '''

    return self._split_up_map_and_model(
      selection_method = 'boxes',
      model_id = model_id,
      target_for_boxes = target_for_boxes,
      select_final_boxes_based_on_model = select_final_boxes_based_on_model,
      skip_empty_boxes = skip_empty_boxes,
      skip_waters = skip_waters,
      box_cushion = box_cushion,
      mask_around_unselected_atoms = mask_around_unselected_atoms,
      mask_radius = mask_radius,
      masked_value = masked_value,
      apply_box_info = apply_box_info,
      mask_id = mask_id,
      exclude_points_outside_density = exclude_points_outside_density,
      minimum_boxes_inside_density = minimum_boxes_inside_density,
      write_files = write_files)

  def _split_up_map_and_model(self,
    model_id = 'model',
    selection_method = 'by_chain',
    selection_list = None,
    skip_waters = False,
    skip_hetero = False,
    target_for_boxes = 24,
    select_final_boxes_based_on_model = True,
    skip_empty_boxes = True,
    mask_around_unselected_atoms = None,
    mask_all_maps_around_edges = None,
    mask_radius = 3,
    masked_value = -10,
    write_files = False,
    box_cushion = 3,
    apply_box_info = True,
    mask_id = None,
    exclude_points_outside_density = None,
    minimum_boxes_inside_density = None,
     ):
    '''
       Create a set of overlapping boxes and non-overlapping parts of
       the working model that cover the entire map

       Returns a group_args object containing list of the map_model_manager
         objects and a list of the selection objects that define which atoms
         from the working model are in each object.

       Normally do work on each map_model_manager to create a new model with
         the same atoms, then use merge_split_maps_and_models() to replace
         coordinates in the original model with those from all the component
         models.

       NOTE: normally you can only change the coordinates and B values in
        each overlapping box if you want to use merge_split_maps_and_models.
        If you want to change the hierarchy of the models, then when you
        split and when you merge, use the keyword:
        allow_changes_in_hierarchy=True. This will create
        one model in the hierarachy for each box, numbered by box number.

      Optionally carry out the step box_info = get_split_maps_and_models(...)
         separately with the keyword apply_box_info=False

       If selection_list (a list of selection objects matching the atoms in
         model) is supplied, use it.  Otherwise generate it using
         selection_method. Skip waters or heteroatoms or both if requested.
         If method is "boxes" then try to get about target_for_boxes boxes.

    If select_final_boxes_based_on_model and selection_method == 'boxes' then
      make the final boxes just go around the selected parts of the model and
      not tile the map.
    If skip_empty_boxes then skip anything with no model.

    If mask_around_unselected_atoms is set, then mask within each box
     around all the atoms that are not selected (including waters/hetero)
     with a mask_radius of mask_radius and set the value inside the mask to
      masked_value

     If exclude_points_outside_density and boxes method is selected,
       try to add boxes inside density (basically add the proportional
       number of boxes but put them definitely inside the density instead
       of evenly spaced.

    '''
    print ("Splitting up map and model into overlapping boxes (%s method)" %(
       selection_method), file = self.log)
    # Get selections and boxes
    box_info = get_selections_and_boxes_to_split_model(
        model_id = model_id,
        map_model_manager = self,
        selection_method = selection_method,
        selection_list = selection_list,
        skip_waters = skip_waters,
        skip_hetero = skip_hetero,
        target_for_boxes = target_for_boxes,
        select_final_boxes_based_on_model = select_final_boxes_based_on_model,
        skip_empty_boxes = skip_empty_boxes,
        box_cushion = box_cushion,
        mask_around_unselected_atoms = mask_around_unselected_atoms,
        mask_all_maps_around_edges = mask_all_maps_around_edges,
        mask_radius = mask_radius,
        masked_value = masked_value,
        mask_id = mask_id,
        exclude_points_outside_density = exclude_points_outside_density,
        minimum_boxes_inside_density = minimum_boxes_inside_density,
        log = self.log,
      )
    if mask_id and exclude_points_outside_density:
      print("Total of %s boxes considered inside density..." %(
        len(box_info.selection_list)), file = self.log)
    else:
      print("Total of %s boxes considered..." %(len(box_info.selection_list)),
        file = self.log)
    if (not apply_box_info):
      return box_info  #  run get_split_maps_and_models later

    # Get new map_model_manager for each box
    box_info = get_split_maps_and_models(
      map_model_manager = self,
      box_info = box_info)
    if write_files and box_info.mmm_list:
      from iotbx.data_manager import DataManager
      dm = DataManager()
      dm.set_overwrite(True)
      i = 0
      for mmm in box_info.mmm_list:
        i += 1
        print("Writing files for model and map: %s " %(i), file=self.log)
        model_file = "model_%s.pdb" %(i) # PDB OK
        map_file = "map_%s.ccp4" %(i)
        model_file = dm.write_model_file(mmm.model(), model_file)
        map_file = dm.write_real_map_file(mmm.map_manager(), map_file)
    return box_info

  # Methods for masking maps ( creating masks and applying masks to maps)
  # These methods change the contents of the current object (they do not
  #  create a new object)

  def mask_all_maps_around_atoms(self,
      model = None,
      mask_atoms_atom_radius = 3,
      set_outside_to_mean_inside = False,
      soft_mask = False,
      soft_mask_radius = None,
      skip_n_residues_on_ends = None,
      invert_mask = None,
      mask_id = 'mask'):
    assert mask_atoms_atom_radius is not None
    '''
      Generate mask around atoms and apply to all maps.
      Overwrites values in these maps

      NOTE: Does not change the gridding or shift_cart of the maps and model

      Optional: specify model to use
      Optionally set the value outside the mask equal to the mean inside,
        changing smoothly from actual values inside the mask to the constant
        value outside (otherwise outside everything is set to zero)

      Optional: radius around atoms for masking
      Optional: soft mask  (default = True)
        Radius will be soft_mask_radius
        (default radius is self.resolution() or resolution calculated
          from gridding)
        If soft mask is set, mask_atoms_atom_radius increased by
      Optional: invert_mask:  keep outside atoms instead of inside
      NOTE: if space group is not p1 and wrapping is True then
        atoms are expanded to P1 before calculating mask
      Optional: skip_n_residues_on_ends: any residues within
        skip_n_residues_on_ends of ends of segments
         (consecutive residue numbers or close N/P atoms) are excluded

    '''
    if not model:
      model = self.model()

    assert model is not None

    if skip_n_residues_on_ends:
      selection_string = " or ".join(get_selections_for_segments(model,
         skip_n_residues_on_ends = skip_n_residues_on_ends))
      model = model.apply_selection_string(selection_string)
      if model is None:  # nothing selected...create empty model
        from mmtbx.model import manager as model_manager
        model = model_manager.from_sites_cart(
           sites_cart = flex.vec3_double(), # empty
           crystal_symmetry = self.crystal_symmetry())
        self.map_manager(
            ).set_model_symmetries_and_shift_cart_to_match_map(model)


    if soft_mask and (not soft_mask_radius):
        soft_mask_radius = self.resolution()
    self.create_mask_around_atoms(
         model = model,
         soft_mask = soft_mask,
         soft_mask_radius = soft_mask_radius,
         mask_atoms_atom_radius = mask_atoms_atom_radius,
         invert_mask = invert_mask,
         mask_id = mask_id)
    self.apply_mask_to_maps(mask_id = mask_id,
         set_outside_to_mean_inside = \
           set_outside_to_mean_inside)

  def mask_all_maps_around_edges(self,
      soft_mask_radius = None,
      boundary_to_smoothing_ratio = 2.,
      mask_id = 'mask'):
    '''
      Apply a soft mask around edges of all maps. Overwrites values in maps
      Use 'mask' as the mask id

      NOTE: Does not change the gridding or shift_cart of the maps and model
    '''
    self.create_mask_around_edges(soft_mask_radius = soft_mask_radius,
      boundary_to_smoothing_ratio = boundary_to_smoothing_ratio,
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
     boundary_to_smoothing_ratio = 2.,
     boundary_radius = None,
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
      soft_mask_radius = self.resolution()
    if not boundary_radius:
      boundary_radius = boundary_to_smoothing_ratio * soft_mask_radius
    from cctbx.maptbx.mask import create_mask_around_edges
    cm = create_mask_around_edges(map_manager = self.map_manager(),
      boundary_radius = boundary_radius)
    cm.soft_mask(soft_mask_radius = soft_mask_radius)

    # Put the mask in map_dict ided with mask_id
    self.add_map_manager_by_id(map_manager = cm.map_manager(),
      map_id = mask_id)

  def create_mask_around_atoms(self,
     model = None,
     mask_atoms_atom_radius = 3,
     soft_mask = False,
     soft_mask_radius = None,
     invert_mask = None,
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
      Optional: invert_mask:  keep outside atoms instead of inside
      NOTE: if space group is not p1 and wrapping is True then
        atoms are expanded to P1 before calculating mask

      Generates new entry in map_manager dictionary with id of
      mask_id (default='mask') replacing any existing entry with that id
    '''

    if soft_mask:
      if not soft_mask_radius:
        soft_mask_radius = self.resolution()
      mask_atoms_atom_radius += soft_mask_radius

    if not model:
      model = self.model()

    from cctbx.maptbx.mask import create_mask_around_atoms
    cm = create_mask_around_atoms(map_manager = self.map_manager(),
      model = model,
      invert_mask = invert_mask,
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
          soft_mask_radius = self.resolution()
      cm.soft_mask(soft_mask_radius = soft_mask_radius)

    # Put the mask in map_dict id'ed with mask_id
    self.add_map_manager_by_id(map_manager = cm.map_manager(),
      map_id = mask_id)

  def create_spherical_mask(self,
     soft_mask = True,
     mask_center_cart = None,
     mask_radius = None,
     soft_mask_radius = None,
     boundary_radius = None,
     boundary_to_smoothing_ratio = 2.,
     mask_id = 'mask',
     minimum_mask_radius_ratio = 0.25 ):

    '''
      Generate spherical mask with radius mask_radius around the
      cartesian point mask_center_cart.  The value of
        mask_center_cart is relative to the
      shifted (origin at (0,0,0) ) position of the map.
      Does not apply the mask to anything.
      Normally follow with apply_mask_to_map or apply_mask_to_maps
      Default: calculate spherical soft mask centered at center of map
        soft_mask radius default is resolution()
        boundary between mask and closest edge is
            soft_mask_radius * boundary_to_smoothing_ratio

      Optional:  not soft mask.  Same as soft mask in dimensions but the
       soft_mask will not be applied

      Optional: mask_center_cart  (default is center of map)


      Generates new entry in map_manager dictionary with id of
      mask_id (default='mask') replacing any existing entry with that id
    '''

    if not soft_mask_radius:
      soft_mask_radius = self.resolution()

    if not boundary_radius:
      boundary_radius = boundary_to_smoothing_ratio * soft_mask_radius

    inner_boundary_radius = soft_mask_radius

    if not mask_center_cart:
      mask_center_cart = self.map_info(quiet=True).working_center_cart

    if not mask_radius:
      # Radius to edge of map:
      mask_radius = flex.double(mask_center_cart).min_max_mean().min

      # Back off by boundary_radius
      mask_radius -= boundary_radius
      mask_radius = max(mask_radius,
         minimum_mask_radius_ratio * boundary_radius)


    if mask_radius <= 0:
      print("\nUnable to auto-generate a mask radius for spherical mask",
        "\nwith center at (%.3f, %.3f, %.3f) A " %(tuple(mask_center_cart)),
        "\nand boundary radius of %.3f A " %( boundary_radius, ),
         file = self.log)
      return

    print("\nGenerating spherical mask ",
        "with center at (%.3f, %.3f, %.3f) A " %(tuple(mask_center_cart)),
        "\n and radius of %.1f A" %(mask_radius), file = self.log)
    if soft_mask:
      print("Mask will be a soft mask", file = self.log)
      print("Boundary radius around "+
        "mask for smoothing: %.1f A  Smoothing radius: %.1f A" %(
          boundary_radius, soft_mask_radius), file = self.log)


    sites_cart = flex.vec3_double()
    sites_cart.append(mask_center_cart)

    model = mmtbx.model.manager.from_sites_cart(
         sites_cart = sites_cart,
         crystal_symmetry = self.crystal_symmetry())

    # Get the same origin shift in model as in our maps
    self.map_manager().set_model_symmetries_and_shift_cart_to_match_map(model)
    from cctbx.maptbx.mask import create_mask_around_atoms
    cm = create_mask_around_atoms(map_manager = self.map_manager(),
      model = model,
      mask_atoms_atom_radius = mask_radius)

    if soft_mask: # Make the create_mask object contain a soft mask
      cm.soft_mask(soft_mask_radius = soft_mask_radius)

    # Put the mask in map_dict ided with mask_id
    self.add_map_manager_by_id(map_manager = cm.map_manager(),
      map_id = mask_id)

  def expand_mask(self,
     buffer_radius = 5,
     resolution = None,
     soft_mask = True,
     soft_mask_radius = None,
     mask_id = 'mask',
      ):
    """Expand the mask with id mask_id by buffer_radius"""
    assert self.get_map_manager_by_id(mask_id)


    map_manager = self.get_map_manager_by_id(mask_id)
    assert map_manager is not None # Need a map to create mask around density
    s =  (map_manager.map_data() > 0.5)
    fraction_old = s.count(True)/s.size()

    from cctbx.maptbx.mask import expand_mask
    em = expand_mask(map_manager = map_manager,
        buffer_radius = buffer_radius,
        resolution = resolution)

    if soft_mask: # Make the create_mask object contain a soft mask
      if not soft_mask_radius:
        if resolution:
          soft_mask_radius = resolution
        else:
          soft_mask_radius = self.resolution()
      em.soft_mask(soft_mask_radius = soft_mask_radius)

    # Put the mask in map_dict id'ed with mask_id
    self.add_map_manager_by_id(map_manager = em.map_manager(),
      map_id = mask_id)
    s =  (self.get_map_manager_by_id(mask_id).map_data() > 0.5)
    fraction_new= s.count(True)/s.size()
    print ( "\nExpanded mask by "+
      "%.1f A ... fraction inside changed from %.4f to %.4f" %(
     buffer_radius, fraction_old,fraction_new), file = self.log)

  # Methods for recombining and manipulating models

  def sequence(self,
    model = None,
    model_id = 'model',
    selection_string = None,
      ):
    '''
      Return sequence of model
    '''
    if not model:
      model = self.get_model_by_id(model_id = model_id)

    if selection_string:
      model = model.apply_selection_string(selection_string)

    return model.as_sequence(as_string = True)

  def rmsd_of_matching_residues(self,
      target_model_id = 'model',
      matching_model_id = None,
      target_model = None,
      matching_model = None,
      max_dist = None,
      minimum_length = None,
      chain_type = None,
      atom_name = None,
      element = None,
      ignore_element = None,
      ca_only = True,
      matching_info_list = None,
      allow_reverse = None,
      quiet = True):
    '''
    Get rmsd of all or ca(P) atoms in matching residues
    '''
    from libtbx import adopt_init_args
    kw_obj = group_args()
    adopt_init_args(kw_obj, locals())
    all_kw = kw_obj() # save calling parameters in kw as dict
    del all_kw['adopt_init_args'] # REQUIRED
    del all_kw['kw_obj']  # REQUIRED
    del all_kw['ca_only']  # REQUIRED
    del all_kw['allow_reverse']  # REQUIRED
    del all_kw['matching_info_list']  # REQUIRED

    if not matching_info_list:
      matching_info_list = self.select_matching_segments(
        max_gap = 0,
        one_to_one = True,
        **all_kw)
    overall_diffs = flex.vec3_double()
    for matching_info in matching_info_list:
      diffs = self.get_diffs_for_matching_target_and_model(
         matching_info = matching_info,
         ca_only = ca_only,
         max_dist = max_dist)
      if allow_reverse:
        reverse_diffs = self.get_diffs_for_matching_target_and_model(
         matching_info = matching_info,
         ca_only = ca_only,
         max_dist = max_dist,
         reverse = True)
        if reverse_diffs and reverse_diffs.size()>0 and \
          reverse_diffs.rms_length() < diffs.rms_length():
          diffs = reverse_diffs
      if diffs:
         overall_diffs.extend(diffs)
    all_n = self.get_info(item_name = 'matching_model_ca_size')


    if overall_diffs.size() > 0:
      self.add_to_info(item_name = 'rms_n', item = overall_diffs.size())
      self.add_to_info(item_name = 'all_n',
          item = all_n if all_n is not None else overall_diffs.size())
      self.add_to_info(item_name = 'rmsd', item = overall_diffs.rms_length())
      return overall_diffs.rms_length()
    else:
      return None

  def get_diffs_for_matching_target_and_model(self,
      matching_info = None,
      max_dist = None,
      ca_only = None,
      reverse = False,
      allow_reverse = False,
       ):
    """Get coordinate differences for matching atoms in target and model"""
    if reverse and not ca_only:
      return None # cannot do reverse for full chain
    if allow_reverse:
      diffs = self.get_diffs_for_matching_target_and_model(
         matching_info = matching_info,
         ca_only = ca_only,
         max_dist = max_dist)
      reverse_diffs = self.get_diffs_for_matching_target_and_model(
         matching_info = matching_info,
         ca_only = ca_only,
         max_dist = max_dist,
         reverse = True)
      if reverse_diffs and reverse_diffs.size()>0 and \
          reverse_diffs.rms_length() < diffs.rms_length():
          diffs = reverse_diffs
      return diffs


    target_model = matching_info.target_model
    matching_model = matching_info.matching_model
    chain_type = matching_info.chain_type
    atom_name = matching_info.atom_name
    element = matching_info.element
    ignore_element = matching_info.ignore_element

    if not atom_name or ((not ignore_element) and (not element)):
      if chain_type.upper() == "PROTEIN":
        atom_name = 'CA'
        element = 'C'
        if max_dist is None:
          max_dist = max(self.resolution(), 3.)
      else:
        atom_name = 'P'
        element = 'P'
    if ca_only:
      if ignore_element:
        target_model_ca = target_model.apply_selection_string(
        "(not hetero) and (altloc ' ' or altloc A) and name %s" %(
           atom_name))
        matching_model_ca = matching_model.apply_selection_string(
        "(not hetero) and  (altloc ' ' or altloc A) and name %s" %(
          atom_name))
      else:  # usual
        target_model_ca = target_model.apply_selection_string(
        "(not hetero) " +\
          "and (altloc ' ' or altloc A) and name %s and element %s" %(
           atom_name, element))
        matching_model_ca = matching_model.apply_selection_string(
        "(not hetero) "+\
           "and  (altloc ' ' or altloc A) and name %s and element %s" %(
           atom_name, element))
    elif (target_model.get_sites_cart().size() != \
         matching_model.get_sites_cart().size() )  or ( not
       target_model.get_hierarchy().atoms().extract_name().all_eq(
       matching_model.get_hierarchy().atoms().extract_name())):
      if self.verbose:
        for x,y in zip(
          target_model.get_hierarchy().atoms().extract_name(),
          matching_model.get_hierarchy().atoms().extract_name()):
          print(x,y, file = self.log)

        print("Target and matching model do not have the same atoms...cannot",
         "use ca_only=False in rmsd_of_matching_residues", file = self.log)
      return None

    if ca_only:
      target_sites = target_model_ca.get_sites_cart()
      matching_sites = matching_model_ca.get_sites_cart()
    else:
      target_sites = target_model.get_sites_cart()
      matching_sites = matching_model.get_sites_cart()
    if target_sites.size() != matching_sites.size():
      return None
    elif reverse:
      matching_sites = list(matching_sites)
      matching_sites.reverse()
      matching_sites = flex.vec3_double(matching_sites)
    diffs = target_sites - matching_sites
    return diffs

  def get_cb_resseq_to_skip(self,cb_resseq_list,
         far_away = None,
         far_away_n = None,
         very_far_away = None,
         very_far_away_n = None,
         ):
    '''
    Identify numbers in resseq_list that are way different than all others by
    at least max_gap
    '''
    if not far_away or not far_away_n:
      return []

    work_list = deepcopy(cb_resseq_list)
    work_list.sort()
    groups=[]
    last_i = None
    for i in work_list:
      if last_i is not None and i <= last_i + far_away:
        group.append(i)
        last_i = i
      else:
         group = [i]
         groups.append(group_args(
           group= group,
           dist_from_last = i-last_i if last_i is not None else 0)
          )
         last_i = i

    residues_to_ignore = []
    for group_info in groups:
     if len(group_info.group) <= far_away_n and \
          len(cb_resseq_list)>= 2 * len(group_info.group):
        residues_to_ignore += group_info.group
     elif group_info.dist_from_last >= very_far_away and \
          len(group_info.group) <= very_far_away_n and \
          len(cb_resseq_list)>= 2 * len(group_info.group):
        residues_to_ignore += group_info.group
    return residues_to_ignore



  def select_matching_segments(self,
      target_model_id = 'model',
      matching_model_id = None,
      target_model = None,
      matching_model = None,
      chain_type = None,
      atom_name = None,
      element = None,
      ignore_element = False,
      max_dist = None,
      max_dist_extra = None,
      minimum_length = None,
      max_gap = 5,
      far_away = 7,
      far_away_n = 2,
      very_far_away = 20,
      very_far_away_n = 1000,
      one_to_one = False,
      residue_names_must_match = False,
      minimum_match_length = 2,
      shift_without_deep_copy = False,
      quiet = True):
    """Select matching pairs of segments in two models"""

    from libtbx import adopt_init_args
    kw_obj = group_args()
    adopt_init_args(kw_obj, locals())
    all_kw = kw_obj() # save calling parameters in kw as dict
    del all_kw['adopt_init_args'] # REQUIRED
    del all_kw['kw_obj']  # REQUIRED


    '''
    Select the parts of matching_model that best match target_model
    without using matching model or target model more than once

    Only use contiguous segments of target_model (i.e., a sequence of residues
     with sequential residue numbers is a contiguous segment).  Does not
     break based on distances (assumes input numbering reflects chain breaks).

    Allow gaps of up to max_gap (and keep residues in the gaps)

    If one-to-one, then select only residues in each that match the other

    If residue_names_must_match, select only residues that match by name as
      well as position

     Return a group_args object with a list of paired sements

    If far_away is set, remove any groups of far_away_n or fewer that are
      far_away residues from all other groups

    If very_far_away and very_far_away_n are set, also remove those

    If target model or matching model have different origins from self, they
      deep-copied and shifted in place to match. The deep copy can be
      skipped if shift_without_deep_copy is set.

    If either model has multiple chains, run all pairwise combinations
    '''

    if one_to_one:
       max_gap = 0

    if target_model:
      if target_model.shift_cart() != self.shift_cart():
        if not shift_without_deep_copy:
          self.add_crystal_symmetry_to_model_if_necessary(target_model,
            map_manager = self.map_manager())
          target_model = target_model.deep_copy()
        self.shift_any_model_to_match(target_model,
         set_unit_cell_crystal_symmetry = True)
    else:
      target_model = self.get_model_by_id(target_model_id)
    if matching_model:
      if matching_model.shift_cart() != self.shift_cart():
        if not shift_without_deep_copy:
          self.add_crystal_symmetry_to_model_if_necessary(matching_model,
            map_manager = self.map_manager())
          matching_model = matching_model.deep_copy()
        self.shift_any_model_to_match(matching_model,
         set_unit_cell_crystal_symmetry = True)
    else:
      matching_model = self.get_model_by_id(matching_model_id)
    if not target_model or not matching_model:
      if not matching_model:
        print("No matching model...skipping comparison", file = self.log)
      if not target_model:
        print("No target model...skipping comparison", file = self.log)
      return []

    target_model_chain_ids = target_model.chain_ids(unique_only=True)
    matching_model_chain_ids = matching_model.chain_ids(unique_only=True)
    if len(target_model_chain_ids) > 1 or len(matching_model_chain_ids) > 1:
      target_and_matching_list = []
      for target_model_chain_id in target_model_chain_ids:
        target_model_chain = target_model.apply_selection_string(
          "chain %s" %(target_model_chain_id))
        for matching_model_chain_id in matching_model_chain_ids:
          matching_model_chain = matching_model.apply_selection_string(
            "chain %s" %(matching_model_chain_id))
          all_kw['target_model'] = target_model_chain
          all_kw['matching_model'] = matching_model_chain
          local_matching_list = self.select_matching_segments(**all_kw)
          if local_matching_list:
            target_and_matching_list += local_matching_list
      return target_and_matching_list


    if not quiet:
      print("\nFinding parts of %s that match %s " %(
       matching_model_id, target_model_id), file = self.log)

    # Get the chain type
    if chain_type is None:
      chain_type = target_model.chain_type()
      if not chain_type:
        chain_type = matching_model.chain_type()
    if not chain_type:
      print("Unable to identify chain_type of '%s' ... please set chain_type" %(
        target_model_id), file = self.log)
      return []

    res = self.resolution() if self.resolution() else 0
    if not atom_name or ((not ignore_element) and (not element)):
      if chain_type.upper() == "PROTEIN":
        atom_name = 'CA'
        element = 'C'
        if max_dist is None:
          max_dist = max(res, 3.)
      else:
        atom_name = 'P'
        element = 'P'
      if max_dist is None:
        max_dist = max(res, 7.)
    elif max_dist is None:
      if atom_name == 'CA':
        max_dist = max(res, 3.)
      else:
        max_dist = max(res, 7.)


    # Select the atoms to try and match
    if ignore_element:
      target_model_ca = target_model.apply_selection_string(
      "(not hetero) and (altloc ' ' or altloc A) and name %s" %(
         atom_name))
      matching_model_ca = matching_model.apply_selection_string(
      "(not hetero) and  (altloc ' ' or altloc A) and name %s" %(
        atom_name))
    else:  # usual
      target_model_ca = target_model.apply_selection_string(
      "(not hetero) and (altloc ' ' or altloc A) and name %s and element %s" %(
         atom_name, element))
      matching_model_ca = matching_model.apply_selection_string(
      "(not hetero) and  (altloc ' ' or altloc A) and name %s and element %s" %(
        atom_name, element))

    # Make sure we have something to work with
    if len(target_model_ca.get_sites_cart()) < 1:
      print("Target model has no sites...skipping select_matching_segments",
         file = self.log)
      return []
    if len(matching_model_ca.get_sites_cart()) < 1:
      print("Matching model has no sites...skipping select_matching_segments",
         file = self.log)
      return []

    self.add_to_info(
        item_name = 'target_model_ca_size', item = target_model_ca.size())
    self.add_to_info(item_name = 'matching_model_ca_size',
         item = matching_model_ca.size())

    if residue_names_must_match:
      ca_residue_names = target_model_ca.as_sequence(as_string = True)
      cb_residue_names = matching_model_ca.as_sequence(as_string = True)
      if ((len(ca_residue_names) != target_model_ca.get_sites_cart().size()) or
         (len(cb_residue_names) != matching_model_ca.get_sites_cart().size())):
        return [] # cannot match up due to something going wrong
    else:
      ca_residue_names = None
      cb_residue_names = None

    matching_cb_as_list_info = self.match_cb_to_ca(
      ca_sites=target_model_ca.get_sites_cart(),
      cb_as_list = list(matching_model_ca.get_sites_cart()),
      ca_residue_names = ca_residue_names,
      cb_residue_names = cb_residue_names,
      max_dist = max_dist,)
    matching_cb_as_list = matching_cb_as_list_info.cb_as_list

    cb_sites_list = list(matching_model_ca.get_sites_cart())
    cb_atoms = list (matching_model_ca.get_hierarchy().atoms())

    cb_atoms_dict = {}
    for cb_site, cb_atom in zip(cb_sites_list,cb_atoms):
     cb_atoms_dict[cb_site] = cb_atom

    ca_sites_list = list(target_model_ca.get_sites_cart())
    ca_atoms = list (target_model_ca.get_hierarchy().atoms())

    ca_atoms_dict = {}
    for ca_site, ca_atom in zip(ca_sites_list,ca_atoms):
     ca_atoms_dict[ca_site] = ca_atom

    if one_to_one:   # take only matching parts and group
      new_ca_sites_list = []
      new_matching_cb_as_list = []
      for ca_site, cb_site in zip(ca_sites_list,
          matching_cb_as_list,
          ):
        if not cb_site:  continue
        new_ca_sites_list.append(ca_site)
        new_matching_cb_as_list.append(cb_site)
      ca_sites_list = new_ca_sites_list
      matching_cb_as_list = new_matching_cb_as_list
      assert len(matching_cb_as_list) == len(ca_sites_list)

      # Select all the residues in ca_sites_list and group by segments
      ca_chain_dict = {}
      for ca_site in ca_sites_list:
        ca = ca_atoms_dict[ca_site]
        ca_chain_id = ca.parent().parent().parent().id
        ca_resseq = ca.parent().parent().resseq_as_int()
        if not ca_chain_id in ca_chain_dict: ca_chain_dict[ca_chain_id] = []
        ca_chain_dict[ca_chain_id].append(ca_resseq)
      ca_selection_list = self.get_selection_string_from_chain_dict(
        chain_dict= ca_chain_dict,
        max_gap = max_gap,
        minimum_length = minimum_length,
        return_as_group_args_list = True)

      ca_sites_groups = []
      matching_cb_as_list_groups = []

      residue_groups = []

      for segment_info in ca_selection_list:
        selection_string ="(chain %s and resseq %s:%s)" %(segment_info.chain_id,
           segment_info.first_resseq, segment_info.last_resseq)
        segment_model_ca = target_model_ca.apply_selection_string(
          selection_string)
        local_ca_sites_list = segment_model_ca.get_sites_cart()
        local_matching_cb_as_list = [] # list of cb coordinates
        for ca, ca_site_cart in zip(segment_model_ca.get_hierarchy().atoms(),
           local_ca_sites_list):
          if ca_site_cart in ca_sites_list:
            index = ca_sites_list.index(ca_site_cart)
            local_matching_cb_as_list.append(matching_cb_as_list[index])
          else:
            local_matching_cb_as_list.append(None)

        residue_groups.append(group_args(
          starting_ca_resseq = segment_info.first_resseq,
          ca_sites_list = local_ca_sites_list,
          matching_cb_as_list = local_matching_cb_as_list,))

    else:
      residue_groups = [
        group_args(
         starting_ca_resseq = None,
         ca_sites_list = ca_sites_list,
         matching_cb_as_list = matching_cb_as_list)
       ]
    sort_list = []
    for rg in residue_groups:
      sort_list.append([rg.starting_ca_resseq,rg])
    sort_list.sort()
    residue_groups = []
    for resseq,rg in sort_list:
      residue_groups.append(rg)

    # Sort the groups by starting numbers

    # Now run through matching as groups

    target_and_matching_list = []
    for residue_group in residue_groups:
      ca_sites_list = residue_group.ca_sites_list
      matching_cb_as_list = residue_group.matching_cb_as_list
      ca_chain_dict = {}
      cb_chain_dict = {}

      cb_resseq_list = []
      for ca_site, cb_site in zip(ca_sites_list,
          matching_cb_as_list, ):
        if not cb_site:  continue
        cb = cb_atoms_dict[cb_site]
        cb_chain_id = cb.parent().parent().parent().id
        cb_resseq = cb.parent().parent().resseq_as_int()
        cb_resseq_list.append(cb_resseq)
      # Are any cb_resseq way out of range of all the others?
      cb_resseq_to_skip = self.get_cb_resseq_to_skip(cb_resseq_list,
         far_away = far_away,
         far_away_n = far_away_n,
         very_far_away = very_far_away,
         very_far_away_n = very_far_away_n,
        )

      for ca_site, cb_site in zip(ca_sites_list,
          matching_cb_as_list, ):
        if not cb_site:  continue
        cb = cb_atoms_dict[cb_site]
        cb_chain_id = cb.parent().parent().parent().id
        cb_resseq = cb.parent().parent().resseq_as_int()
        if cb_resseq in cb_resseq_to_skip: continue


        ca = ca_atoms_dict[ca_site]
        ca_chain_id = ca.parent().parent().parent().id
        ca_resseq = ca.parent().parent().resseq_as_int()
        if residue_names_must_match:
          assert ca.parent().resname == cb.parent().resname

        if not cb_chain_id in cb_chain_dict: cb_chain_dict[cb_chain_id] = []
        if not ca_chain_id in ca_chain_dict: ca_chain_dict[ca_chain_id] = []
        cb_chain_dict[cb_chain_id].append(cb_resseq)
        ca_chain_dict[ca_chain_id].append(ca_resseq)

      cb_selection_string = self.get_selection_string_from_chain_dict(
       chain_dict= cb_chain_dict, max_gap = max_gap)
      local_matching_model = matching_model.apply_selection_string(
         cb_selection_string)

      if one_to_one:   # take only matching parts
        # Note they may be in different orders
        ca_selection_string = self.get_selection_string_from_chain_dict(
          chain_dict= ca_chain_dict, max_gap = max_gap)
        local_target_model = \
         target_model.apply_selection_string(ca_selection_string)
        if not local_target_model:
          continue # skip it
        target_seq = local_target_model.as_sequence(as_string = True)
        matching_seq = local_matching_model.as_sequence(as_string = True)
        if not target_seq or not matching_seq:
          continue # skip it
        target_seq=list(target_seq)
        matching_seq=list(matching_seq)
        target_seq.sort()
        matching_seq.sort()
        if residue_names_must_match:
          assert target_seq == matching_seq  # same but could be different order
        if minimum_match_length and len(target_seq) < minimum_match_length:
          continue # skip it

      else:
        local_target_model = target_model
        local_matching_model = matching_model

      t_chain_ids = local_target_model.chain_ids()
      t_chain_id = t_chain_ids[0] if t_chain_ids else None
      m_chain_ids = local_matching_model.chain_ids()
      m_chain_id = m_chain_ids[0] if m_chain_ids else None
      target_and_matching = group_args(
        group_args_type = 'target and matching residues from other',
          id = len(target_and_matching_list),
          target_model = local_target_model,
          target_model_chain_id = t_chain_id,
          target_model_start_resseq = local_target_model.first_resseq_as_int(),
          target_model_end_resseq =  local_target_model.last_resseq_as_int(),
          matching_model = local_matching_model,
          matching_model_chain_id = m_chain_id,
          matching_model_start_resseq =  local_matching_model.first_resseq_as_int(),
          matching_model_end_resseq =  local_matching_model.last_resseq_as_int(),
          chain_type = chain_type,
          atom_name = atom_name,
          element = element,
          ignore_element = ignore_element,
          )

      diffs = self.get_diffs_for_matching_target_and_model(
         matching_info = target_and_matching,
         ca_only = True,
         max_dist = max_dist)
      rms_diffs = diffs.rms_length() if diffs and diffs.size()>0 else None
      reverse_diffs = self.get_diffs_for_matching_target_and_model(
         matching_info = target_and_matching,
         ca_only = True,
         max_dist = max_dist,
         reverse = True)
      rms_reverse_diffs = \
         reverse_diffs.rms_length() if reverse_diffs and reverse_diffs.size()>0 else None
      if rms_diffs is not None and (rms_reverse_diffs is None or
          rms_diffs <= rms_reverse_diffs):
        target_and_matching.match_direction = True
        target_and_matching.rms_diffs = rms_diffs
      elif rms_diffs is not None and (rms_reverse_diffs is not None
          and rms_diffs > rms_reverse_diffs):
        target_and_matching.match_direction = False
        target_and_matching.rms_diffs = rms_reverse_diffs
      else:
        target_and_matching.match_direction = None
        target_and_matching.rms_diffs = None

      target_and_matching_list.append(target_and_matching)
    return target_and_matching_list

  def get_selection_string_from_chain_dict(self,
     chain_dict = None,
     max_gap = None,
     minimum_length = None,
     return_as_group_args_list = False):
   '''
     Return a selection string for the segments represented in chain_dict,
     allowing gaps of up to max-gap (fill them in)
     Require minimum_length if set
   '''
   selection_string_list = []
   selection_group_args_list = []
   for chain_id in chain_dict.keys():
     resseq_list = chain_dict[chain_id]
     resseq_list.sort()
     groups = []
     previous_resseq = None
     first_resseq = None
     for resseq,next_resseq in zip(resseq_list,resseq_list[1:]+[None]):
       if previous_resseq is None:
         first_resseq = resseq

       if next_resseq is None or next_resseq > resseq + max_gap+1: # break
         groups.append(
           group_args(first_resseq = first_resseq, last_resseq= resseq))
         previous_resseq = None
       else:
         previous_resseq = resseq
     for group in groups:
       if minimum_length is not None and (
          group.last_resseq-group.first_resseq+1) < minimum_length:
         continue
       selection_string_list.append("(chain %s and resseq %s:%s)" %(
          chain_id, group.first_resseq, group.last_resseq))
       selection_group_args_list.append(
         group_args( chain_id = chain_id,
             first_resseq = group.first_resseq,
             last_resseq = group.last_resseq,) )

   if return_as_group_args_list:
     return selection_group_args_list
   else:
     return " or ".join(selection_string_list)


  def choose_best_set(self,dd, max_dist = None,
      ):
    """Input is a dictionary of lists of group_args objects.  Output is
    group_args object containing a new dictionary. The new dictionary has
    only one group_args object for each key.  The one chosen is the one
    from the corresponding list that has the smallest value of dist."""

    # dd is a dict
    # values for dd are lists of group args, each has a member value of dist
    #   and a member value of id.  Choose the one that has the smallest dist

    # If there is a close alternative, also within max_dist but with residue
    #  number much closer to all the others in a group, take that one instead
    #  This is to try and go along a chain and not jump unless necessary.

    #   ca_dict[cb].append(group_args(
    #     dist=dist,
    #     id = ca,
    #     resname_id = ca_residue_names[id2], # resname of ca atom
    #     index_id = id2,  #index of ca atom
    #     other_resname_id = cb_residue_names[index_cb],  # name of cb atom
    #     other_index_id =  index_cb,  #  index of cb atob


    for key in dd.keys():
      groups = dd[key]
      best_group = None
      for group in groups:
        if (max_dist is not None) and (group.dist > max_dist):
          continue
        if best_group is None or group.dist < best_group.dist:
          best_group = group
      dd[key] = best_group

    # Now remove any duplicate id's
    target_id = self.duplicate_id(dd)
    while target_id:
      duplicate_list = []
      for key in dd.keys():
        if dd[key] and dd[key].id == target_id:
          duplicate_list.append(key)
      assert duplicate_list
      best_key = None
      for key in duplicate_list:
        if not best_key or dd[key].dist < dd[best_key].dist:
          best_key = key
      for key in duplicate_list:
        if key != best_key:
          dd[key] = None
      target_id = self.duplicate_id(dd)
    return group_args(
      dd = dd,
      )

  def duplicate_id(self,dd):
    '''Return first id in dict dd that is a duplicate'''
    id_list = []
    for key in dd.keys():
      if dd[key]:
        id = dd[key].id
        if not id in id_list:
          id_list.append(id)
        else:
          return id
    return None

  def match_cb_to_ca(self,
     ca_sites=None,
     cb_as_list = None,
     ca_residue_names = None,
     cb_residue_names = None,
     max_dist = 2.,):

    '''
     Identify cb sites that match ca sites
     If ca_residue_names and cb_residue names, require that residue names match
    '''

    cb_dict = {}
    ca_dict = {}

    if not ca_residue_names:
      ca_residue_names = ca_sites.size() * ['CA']
    if not cb_residue_names:
      cb_residue_names = len(cb_as_list) * ['CA']
    for index_cb in range(len(cb_as_list)):
      cb = cb_as_list[index_cb]
      if cb is not None and cb != (-9999, -9999, -9999):
        cb_sites = flex.vec3_double()
        cb_sites.append(cb)
        dist, id1, id2 = cb_sites.min_distance_between_any_pair_with_id(
          ca_sites)
        ca = ca_sites[id2]
        if not (ca in cb_dict.keys()):
          cb_dict[ca] = []

        if ca_residue_names[id2] != cb_residue_names[index_cb]:
          dist = 1.e+30 # do not match them
        cb_dict[ca].append(
          group_args(
            dist=dist,
            id=cb,
            resname_id = cb_residue_names[index_cb],  # name of cb atom
            index_id =  index_cb,  #  index of cb atob
            other_resname_id = ca_residue_names[id2], # resname of ca atom
            other_index_id = id2,  #index of ca atom
           ))
        if not cb in ca_dict.keys():
          ca_dict[cb] = []
        ca_dict[cb].append(group_args(
          dist=dist,
          id = ca,
          resname_id = ca_residue_names[id2], # resname of ca atom
          index_id = id2,  #index of ca atom
          other_resname_id = cb_residue_names[index_cb],  # name of cb atom
          other_index_id =  index_cb,  #  index of cb atob
          ))
    cb_dict_info = self.choose_best_set(cb_dict, max_dist = max_dist)
    ca_dict_info = self.choose_best_set(ca_dict, max_dist = max_dist)
    cb_dict = cb_dict_info.dd
    ca_dict = ca_dict_info.dd
    cb_as_list = []

    # Make a list of all the ca_dict positions and the matching (or missing)
    #  item from cb
    if not ca_residue_names:
      ca_residue_names = [None]*ca_sites.size()
    for ca,ca_resname in zip(ca_sites, ca_residue_names):

      group = cb_dict.get(ca)
      if group:
        cb_as_list.append(group.id)
        assert (not ca_resname) or (ca_resname == group.resname_id and
           ca_resname == group.other_resname_id)
      else:
        cb_as_list.append(None)

    return group_args(
       cb_as_list = cb_as_list,
       )


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

  def add_crystal_symmetry_to_model_if_necessary(self,
       model, map_manager = None):
    '''
    Take any model and add crystal symmetry if it is missing
    Changes model in place
     Parameters:  model

    Also add unit_cell_crystal_symmetry if missing
    '''
    if not model:
      return
    if not map_manager:
      map_manager = self.any_map_manager()
    assert map_manager is not None

    assert isinstance(model, mmtbx.model.manager)
    if (not model.crystal_symmetry()) or (
        not model.crystal_symmetry().unit_cell()):
      map_manager.set_model_symmetries_and_shift_cart_to_match_map(model)

  def shift_any_model_to_match(self, model, map_manager = None,
     set_unit_cell_crystal_symmetry = False):
    '''
    Take any model and shift it to match the working shift_cart
    Also sets crystal_symmetry.
    Changes model in place

     Parameters:  model
                  set_unit_cell_crystal_symmetry:  optionally set this as well
    '''
    if not model:
      return
    assert isinstance(model, mmtbx.model.manager)

    if not map_manager:
      map_manager = self.get_any_map_manager()

    if not map_manager:
      return # Do not shift model if no map_manager (no shifts known)

    if map_manager.is_dummy_map_manager():
      return # Do not shift model if no map_manager (no shifts known)

    self.add_crystal_symmetry_to_model_if_necessary(
        model, map_manager = map_manager)

    if not model.shift_cart():
      model.set_shift_cart((0, 0, 0))

    coordinate_shift = tuple(
      [s - o for s,o in zip(map_manager.shift_cart(),model.shift_cart())])

    model.shift_model_and_set_crystal_symmetry(
        shift_cart = coordinate_shift,
        crystal_symmetry=map_manager.crystal_symmetry())

    if set_unit_cell_crystal_symmetry and self.map_manager():
       self.set_model_symmetries_and_shift_cart_to_match_map(model)


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

    other_shift_cart = other_model.shift_cart()
    other_model.shift_model_back()  # removes shift cart (adds -other_shift_cart)

    other_model.set_crystal_symmetry(self.crystal_symmetry())
    other_model.set_unit_cell_crystal_symmetry(self.unit_cell_crystal_symmetry())

    other_model.shift_model_and_set_crystal_symmetry(
        shift_cart = self.shift_cart())

    return other_model

  # Methods to create a new map_model_manager with different sampling

  def as_map_model_manager_with_resampled_maps(self,
     sampling_ratio = 2):
    ''' Return a new map_model_manager with maps sampled more finely
        Parameter:  sampling_ratio  must be an integer
        Creates new maps, keeps same models
        Sampling ratio must be greater than 1
    '''

    assert sampling_ratio > 1

    n_real = tuple([
       int(sampling_ratio *n) for n in self.map_manager().map_data().all()])

    map_id = 'map_manager'
    used_map_id_list = [map_id]
    fine_mm = self.get_map_manager_by_id(map_id
         ).resample_on_different_grid(n_real)
    new_mmm = map_model_manager(map_manager = fine_mm,)

    for model_id in self.model_id_list():
      new_mmm.add_model_by_id(model = self.get_model_by_id(model_id),
        model_id = model_id)

    for map_id in self.map_id_list():
      if not map_id in used_map_id_list:
         used_map_id_list.append(map_id)
      fine_mm = self.get_map_manager_by_id(map_id
         ).resample_on_different_grid(n_real)
      new_mmm.add_map_manager_by_id(map_manager = fine_mm,
        map_id = map_id)

    return new_mmm

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
      new_map_id = 'map_manager',
      ):
    '''
      Resolution-filter a map with range of d_min to d_max and place in
      new map (can be the same)

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
      map_id = new_map_id)


  # Methods for modifying model or map

  def remove_model_outside_map(self, model = None,
        boundary = 3, return_as_new_model=False):
    '''
     Remove all the atoms in the model that are well outside the map (more
     than boundary). Boundary can be negative (remove inside box near edges)
    '''
    assert boundary is not None

    if not model:
      model = self.model()

    if not model:
      return

    sites_frac = model.get_sites_frac()
    bf_a, bf_b, bf_c = model.crystal_symmetry().unit_cell(
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
      return model.select(~s)
    else:  # usual
      self.add_model_by_id( model.select(~s), 'model')


  # Methods for sharpening and comparing maps, models and calculating FSC values

  def get_rms_f_list(self,
      map_id = 'map_manager',
      d_min = None,  # minimum resolution for calculations
      n_bins = None,
      resolution = None,  # nominal resolution
      ):
    ''' Return list of rms amplitude by bins '''
    assert d_min and n_bins
    from cctbx.maptbx.segment_and_split_map import map_coeffs_to_fp
    map_coeffs=self.get_map_manager_by_id(
      map_id).map_as_fourier_coefficients(d_min = d_min)
    f_array = get_map_coeffs_as_fp_phi(map_coeffs, n_bins = n_bins,
      d_min = d_min).f_array
    rms_f_list = flex.double()
    sthol2_list = flex.double()
    dsd = f_array.d_spacings().data()
    n_bins_use = min(n_bins,max(3,n_bins//3))
    for i_bin in f_array.binner().range_used():
      sel = f_array.binner().selection(i_bin)
      f = map_coeffs.select(sel)
      f_array_f = map_coeffs_to_fp(f)
      rms_f = f_array_f.data().norm()
      rms_f_list.append(rms_f)
      d = dsd.select(sel)
      if d.size() < 0:
        d_avg = flex.mean(d)
        sthol2 = 0.25/d_avg**2
        sthol2_list.append(sthol2)
        if i_bin-1 > n_bins_use and ((not resolution) or (d_avg >= resolution)):
          n_bins_use = i_bin - 1
      elif i_bin > 1:
        sthol2_list.append(sthol2_list[-1])
      else:
        sthol2_list.append(0)
    return group_args(
       rms_f_list = rms_f_list,
       sthol2_list = sthol2_list,
       n_bins_use = n_bins_use)


  def _update_kw_with_map_info(self, local_kw, previous_kw = None,
      text = 'overall', have_previous_scaled_data = None,
      map_id_scaled_list = None):
    """Update keywords with information from a map"""

    if have_previous_scaled_data:
      # Expect map_id_1 to be in previous_kw['map_id_to_be_scaled_list']...
      if previous_kw.get('map_id_1') and previous_kw.get('map_id_2'):
        local_kw['map_id_1'] = previous_kw['map_id_scaled_list'][
         previous_kw['map_id_to_be_scaled_list'].index(previous_kw['map_id_1'])]
        print("Map 1 to use in determining scaling: '%s' " %(
          local_kw['map_id_1']), file = self.log)
        local_kw['map_id_2'] = previous_kw['map_id_scaled_list'][
         previous_kw['map_id_to_be_scaled_list'].index(previous_kw['map_id_2'])]
        print("Map 2 to use in determining scaling: '%s' " %(
          local_kw['map_id_2']), file = self.log)

      local_kw['map_id_to_be_scaled_list'] = previous_kw['map_id_scaled_list']

      # New main map
      local_kw['map_id'] = previous_kw['map_id_scaled_list'][
          previous_kw['map_id_to_be_scaled_list'].index(previous_kw['map_id'])]
      print("New main map to use in determining scaling: '%s' " %(
          local_kw['map_id']),file = self.log)

    if map_id_scaled_list:
      local_kw['map_id_scaled_list'] = map_id_scaled_list
    else:
      local_kw['map_id_scaled_list'] = []
      for id in local_kw['map_id_to_be_scaled_list']:
        local_kw['map_id_scaled_list'].append(
            '%s_%s' %(id, text))


    print("Starting maps will come from: %s " %(
      str(local_kw['map_id_to_be_scaled_list'])),file = self.log)
    print("Sharpened maps will be in: %s " %(
       local_kw['map_id_scaled_list']), file = self.log)

    if local_kw.get('model_id') and self.get_model_by_id(local_kw['model_id']):
      cc = self.map_model_cc(map_id=local_kw['map_id'])
      print ("Current map-model CC for '%s': %.3f " %(local_kw['map_id'],cc),
         file = self.log)

    return local_kw

  def get_map_rmse(self,
      map_id = 'map_manager',
      mask_id = 'mask',
      fc_map_id = 'fc_for_model_comparison',
      fofc_map_id = 'fofc_for_model_comparison',
      model = None,
      d_min = None,
      r_free = None, r_work = None,
      match_overall_b = True,):
    """  Estimate map rms error by comparison with model in region containing
      the model. Suitable for cases where model is not refined against this
      map. If refined against this map, correct for degrees of freedom by
      estimating overfitting from ratio of Rfree/Rwork.
      If match_overall_b is True, adjust average B of model to match fall-off
      of the Fourier coefficients representing the map
      Return group_args object with:
        map_rmse = fofc_info_inside.sd  # map rms error
        map_rmse_to_sd_ratio = ratio_error_to_input_map  # rmse / sd of map
        map_sd = map_info.sd  # SD of the map (rms after setting mean to zero)

     """

    mmm = self.deep_copy() # we are going to modify things

    if r_free is not None and r_work is not None:
      overfitting_ratio = max(1., r_free/max(1.e-10, r_work))
    else:
      overfitting_ratio = 1

    if model is None:
      model = mmm.model()

    if d_min is None:
      d_min = mmm.resolution()

    if match_overall_b:   #Match B to map
      model = mmm.match_model_b_to_map(map_id = map_id, model = model,
       d_min = d_min,)

    mmm.generate_map(model = model,
         gridding=mmm.get_any_map_manager().map_data().all(),
         d_min=d_min,
         map_id = fc_map_id)

    mmm.create_mask_around_atoms(
       model = mmm.model(),
       mask_atoms_atom_radius = 2.* d_min,
       soft_mask = False,
      )

    map_info = mmm._get_mean_sd_of_map(map_id = map_id)
    fc_info= mmm._get_mean_sd_of_map(map_id = fc_map_id)
    map_info_inside = mmm._get_mean_sd_of_map(
       map_id = map_id, mask_id = mask_id)
    fc_info_inside= mmm._get_mean_sd_of_map(
       map_id = fc_map_id, mask_id = mask_id)

    # get ratio of rms inside mask to overall
    rms_inside_to_all_ratio = map_info_inside.sd / max(1.e-10, map_info.sd)

    # And offset to add to fc_map to make insides match
    offset = map_info_inside.mean - fc_info_inside.mean

    # and ratio of rms in map to fc
    ratio = map_info_inside.sd / max(1.e-10, fc_info_inside.sd)

    # Get optimal ratio to minimize fo-fc inside mask
    scale_ratio = mmm._get_scale_ratio(map_id_1 = map_id,
       map_id_2 = fc_map_id, mask_id = mask_id)

    print("Ratio of sd values: %.2f   Optimal scale_ratio: %.2f "%(
        ratio, scale_ratio), file = self.log)

    input_map = mmm.get_map_manager_by_id(map_id).map_data()

    # Set mean of input map inside mask to 0
    input_map -= map_info_inside.mean

    fc_map = mmm.get_map_manager_by_id(fc_map_id).map_data()
    # Offset and scale fc map to match input map inside mask
    fc_map += offset
    fc_map *= scale_ratio
    fofc = input_map - fc_map
    mm_fofc = mmm.get_map_manager_by_id(map_id).customized_copy(
       map_data = fofc)
    mmm.add_map_manager_by_id(map_manager = mm_fofc, map_id = fofc_map_id)


    fofc_info_inside = mmm._get_mean_sd_of_map(
       map_id = fofc_map_id, mask_id = mask_id)
    ratio_error_to_input_map_inside = fofc_info_inside.sd / max(
           1.e-10, map_info_inside.sd)
    ratio_error_to_input_map= fofc_info_inside.sd / max(
           1.e-10, map_info.sd)
    print("Mean overall for input map: %.2f  Inside mask: %.2f Diff: %.2f" %(
       map_info.mean, map_info_inside.mean,
       map_info.mean - map_info_inside.mean), file = self.log)

    print("RMS overall for input map: %.2f  Inside mask: %.2f  Ratio: %.2f" %(
       map_info.sd, map_info_inside.sd, map_info.sd /max(1.e-10,
           map_info_inside.sd) ), file = self.log)
    print(
      "RMS inside for fc map: %.2f  Ratio input map to fc map inside: %.2f" %(
       fc_info_inside.sd, ratio), file = self.log)

    print("RMS inside mask for fo-fc map: %.2f" %(fofc_info_inside.sd),
      file = self.log)

    print(
     "Ratio rms (inside mask) fo-fc to RMS (inside mask) input map: %.2f" %(
      ratio_error_to_input_map_inside), file = self.log)
    print("Ratio rms (inside mask) fo-fc to RMS (all) input map: %.2f" %(
      ratio_error_to_input_map), file = self.log)

    map_rmse = fofc_info_inside.sd * overfitting_ratio  # map rms error
    map_rmse_to_sd_ratio = ratio_error_to_input_map * overfitting_ratio
    # rmse / sd of map
    print("Estimated ratio of map rms error to map sd: %.2f" %(
       map_rmse_to_sd_ratio), file = self.log)
    print("Estimated map rms error: %.2f" %(map_rmse), file = self.log)

    return group_args(group_args_type = 'Map error estimates',
      map_rmse = map_rmse,  # map rms error in region of macromolecule
      map_rmse_to_sd_ratio = map_rmse_to_sd_ratio,  # rmse / sd of map
      map_sd = map_info.sd,  # SD of the map (rms after setting mean to zero)
      map_mean = map_info.mean, # Mean of the map
     )


  def match_model_b_to_map(self, map_id = 'map_manager', model = None,
       d_min = None,):
      """Adjust model B-value to best match a map"""
      from cctbx.maptbx.segment_and_split_map import get_b_iso
      from cctbx.maptbx.segment_and_split_map import map_coeffs_as_fp_phi
      if not d_min:
        d_min = self.d_min()
      if not model:
        model = self.model()
      map_coeffs = self.map_as_fourier_coefficients(map_id = map_id,
          d_min = d_min,)
      f,phi=map_coeffs_as_fp_phi(map_coeffs)
      b_mean,aniso_scale_and_b=get_b_iso(f,d_min=d_min,
          return_aniso_scale_and_b=True)

      mean_b_in_model = model.get_b_iso().min_max_mean().mean
      low_b_in_model = model.get_b_iso().min_max_mean().min
      offset_b = (b_mean - mean_b_in_model)
      new_b = max(low_b_in_model,model.get_b_iso() + offset_b)
      model.set_b_iso(new_b)
      print("Matching model to map by setting B average from %.2f to %.2f " %(
          mean_b_in_model,model.get_b_iso().min_max_mean().mean),
         file = self.log)
      return model

  def _get_scale_ratio(self, map_id_1 = 'map_manager',
       map_id_2 = 'fofc_for_model_comparison', mask_id = 'mask'):
    """Optimal scale to apply to map 2 to minimize rms difference from map 1"""
    map_data_1 = self.get_map_manager_by_id(map_id_1).map_data().as_1d()
    map_data_2 = self.get_map_manager_by_id(map_id_2).map_data().as_1d()
    if mask_id is not None:
      mask_data = self.get_map_manager_by_id(mask_id).map_data().as_1d()
      sel = (mask_data > 0.5)
      map_data_1 = map_data_1.select(sel)
      map_data_2 = map_data_2.select(sel)
    map_data_1 = map_data_1 - map_data_1.min_max_mean().mean
    map_data_2 = map_data_2 - map_data_2.min_max_mean().mean
    cc = flex.linear_correlation(map_data_1,map_data_2).coefficient()
    ratio = (map_data_1 * map_data_2).min_max_mean().mean / max( 1.e-10,
     flex.pow2(map_data_2).min_max_mean().mean)
    return ratio

  def _get_mean_sd_of_map(self, map_id = 'map_manager', mask_id = None):
    """
    Get mean and sd of map specified by map_id, inside mask if mask_id set
    """
    map_data = self.get_map_manager_by_id(map_id).map_data().as_1d()

    if mask_id is not None:
      mask_data = self.get_map_manager_by_id(mask_id).map_data().as_1d()
      map_data = map_data.select(mask_data > 0.5)

    mean_value = flex.mean(map_data)
    sd = map_data.sample_standard_deviation()
    return group_args(group_args_type = 'map sd',
      map_id = map_id,
      mask_id = mask_id,
      mean = mean_value,
      sd = sd,
      n=map_data.size())
  def find_k_sol_b_sol(self,
    model = None,
    d_min = None,
    model_map_id = None,
    comparison_map_id = None,
    n_bins = 5):

    ''' Routine to guess k_sol and b_sol by low-resolution Fc calculation'''

    if model_map_id is None:
      model_map_id = 'map_from_model'
    if comparison_map_id is None:
      comparison_map_id = 'map_manager'

    kb_list= [ [0,0],
                      [0.1,20], [0.1,50],
                      [0.2,20], [0.2,50],
                      [0.3,20], [0.3,50],
                      [0.15,20], [0.15,50],
                      [0.15,30], [0.15,40],
                      [0.15,10], [0.15,60],
                      [0.15,0], [0.15,5],
             ]

    from cctbx.development.create_models_or_maps import \
       generate_map_coefficients

    target_map_coeffs = self.get_map_manager_by_id(
       comparison_map_id).map_as_fourier_coefficients( d_min = d_min)
    (d_max,d_min)=target_map_coeffs.d_max_min(
       d_max_is_highest_defined_if_infinite=True)
    target_map_coeffs.setup_binner(n_bins = n_bins,
      d_max=d_max,
      d_min=d_min)
    best_kb = None
    best_cc = None
    for k_sol,b_sol in kb_list:
      map_coeffs = generate_map_coefficients(model = model,
        d_min = d_min,
        k_sol = k_sol,
        b_sol = b_sol,
        scattering_table = self.scattering_table(),
        f_obs_array = target_map_coeffs,
        log = null_out())
      sum_cc = flex.double()
      for i_bin in target_map_coeffs.binner().range_used():
        sel       = target_map_coeffs.binner().selection(i_bin)
        cc1 = map_coeffs.select(sel).map_correlation(
           target_map_coeffs.select(sel))
        cc1 = (0 if cc1 is None else cc1)
        sum_cc.append(cc1)
      cc = sum_cc.min_max_mean().mean
      if best_cc is None or cc > best_cc:
        best_cc = cc
        best_kb = [k_sol,b_sol]
    return group_args(
      k_sol = best_kb[0],
      b_sol = best_kb[1],
      cc = best_cc)


  def tls_from_map(self,
    map_id_1 = None,
    map_id_2 = None,
    map_id = None,
    model_id = None,
    mask_id = None,
    tls_by_chain = True,
    apply_tls_to_model = True,
    iterations = 1,
    skip_waters = True,
    skip_hetero = True,
    coordinate_shift_to_apply_before_tlso = None,
    core_box_size_ratio = None,
    box_cushion_ratio = None,
    exclude_points_outside_density = True,
    minimum_boxes_inside_density = True,
    d_min = None,
      **kw):
    """Estimate TLS parameters from a map"""
    if iterations:
      from libtbx import adopt_init_args
      kw_obj = group_args()
      adopt_init_args(kw_obj, locals())
      all_kw = kw_obj() # save calling parameters in kw as dict
      del all_kw['adopt_init_args'] # REQUIRED
      del all_kw['kw_obj']  # REQUIRED

      all_kw.update(kw)
      del all_kw['kw']
      all_kw['iterations'] = None
      all_kw_use = deepcopy(all_kw)
      print("\nRunning total of %s iterations of TLS from map " %(iterations),
        file = self.log)
      for iter in range(iterations-1):
        print("\nRunning iteration %s of %s of TLS from map" %(
         iter+1,iterations), file = self.log)
        result = self.tls_from_map(**all_kw_use)
      print("\nDone running extra iterations of TLS from map ",file = self.log)

    if model_id is None:
      model_id = 'model'

    # Save all keywords we want to pass on in kw
    kw['map_id_1'] = map_id_1
    kw['map_id_2'] = map_id_2
    kw['map_id'] = map_id
    kw['model_id'] = model_id
    kw['exclude_points_outside_density'] = exclude_points_outside_density
    kw['d_min'] = d_min

    # Set up list of maps to be scaled and kw
    kw = self.set_map_id_lists(kw)

    # Set keywords for tls_from_map
    kw['local_sharpen'] = True
    kw['anisotropic_sharpen'] = True
    kw['get_scale_as_aniso_u'] = True
    kw['get_tls_from_u'] = True
    kw['get_tls_info_only'] = True
    kw['replace_aniso_with_tls_equiv'] = False
    kw['overall_sharpen_before_and_after_local'] = False
    kw['coordinate_shift_to_apply_before_tlso'] =\
        coordinate_shift_to_apply_before_tlso

    print("\nRunning tls_from_map...\n",
       file = self.log)
    if kw.get('map_id_1') and kw.get('map_id_2'):
      print("\nTLS will be determined by comparison of %s and %s " %(
       kw['map_id_1'],kw['map_id_2']), file = self.log)
      method = self.half_map_sharpen
      del kw['model_id']
      del kw['map_id']
    elif kw.get('map_id') and kw.get('model_id'):
      print("\nTLS will be determined by comparison of %s and %s " %(
       kw['map_id'],kw['model_id']), file = self.log)
      method = self.model_sharpen
      del kw['map_id_1']
      del kw['map_id_2']
    else:
      raise Sorry("Need two half-maps or map and model for get_tls_from_map")

    # Run by chain if requested
    if tls_by_chain:
      print("TLS will be determined for each chain", file = self.log)
      box_info = self._split_up_map_and_model(
        model_id = model_id,
        selection_method = 'by_chain',
        skip_waters = skip_waters,
        skip_hetero = skip_hetero,
        mask_all_maps_around_edges = False,)
      tlso_list = []
      for mmm in box_info.mmm_list:
        # working shift_cart is shift from original to working xyz
        #  box shift_cart is original to box
        #  to get coords in working frame, take box xyz and subtract box shift
        #    then add working shift
        coordinate_shift = tuple(
         [working_shift - box_shift for working_shift,box_shift in zip(
           self.map_manager().shift_cart(), mmm.map_manager().shift_cart())]
         )
        kw['coordinate_shift_to_apply_before_tlso'] = coordinate_shift

        box_info =  mmm.tls_from_map(
         core_box_size_ratio = core_box_size_ratio,
         box_cushion_ratio = box_cushion_ratio,
         tls_by_chain = False,
         apply_tls_to_model = False,
         iterations = None,
          **kw)
        for tlso in box_info.tlso_list:
          tlso_list.append(tlso)
      box_info.tlso_list = tlso_list

    else:
      print("TLS will be determined for entire model as one group",
         file = self.log)
      if core_box_size_ratio and (not kw.get('core_box_size')):
        kw['core_box_size'] = core_box_size_ratio * self.resolution()  # in A,
      if box_cushion_ratio and (not kw.get('box_cushion')):
        kw['box_cushion'] = box_cushion_ratio * self.resolution()  # in  A

      tls_info = method(**kw)
          # run overall sharpening
      tlso_list = [tls_info.tlso]
      mmm_list = [self]
      box_info = group_args(
       selection_list = None,
       selection_as_text_list = None,
       tlso_list = tlso_list,
       mmm_list = [self])

    if apply_tls_to_model and model_id and \
       self.get_model_by_id(model_id = model_id):
          # set the values in the model using
      if not box_info.selection_list:
        box_info.selection_list = [
          self.get_model_by_id(model_id = model_id).selection('all')]
        box_info.selection_as_text_list = ['all']

      self.merge_split_maps_and_models(
        model_id = model_id,
        box_info = box_info,
        replace_coordinates = False,
        replace_u_aniso = True)
    return box_info

  def _sharpen_overall_local_overall(self, kw, method):
      """Run sharpening without local sharpening first
      Then set maps to scale as the scaled maps from this run"""

      assert kw.get('map_id_to_be_scaled_list') is None or (
        kw['map_id']  in kw['map_id_to_be_scaled_list']) # map_id_to_be_scaled not ok

      # Set up list of maps to be scaled
      kw['sharpen_all_maps'] = True # REQUIRED
      kw = self.set_map_id_lists(kw) # MUST COME AFTER sharpen_all_maps
      kw['overall_sharpen_before_and_after_local'] = False


      final_map_id_scaled_list = deepcopy(kw['map_id_scaled_list'])
      print("\nRunning overall sharpening, local , then overall...\n",
         file = self.log)
      if kw.get('map_id_1') and kw.get('map_id_2'):
        print("\nSharpening will be determined by comparison of %s and %s " %(
         kw['map_id_1'],kw['map_id_2']), file = self.log)
      print("Starting maps will come from: %s " %(
        str(kw['map_id_to_be_scaled_list'])),file = self.log)
      print("Final sharpened maps will be in: %s " %(final_map_id_scaled_list),
        file = self.log)

      # Run overall sharpening
      local_kw = deepcopy(kw)  # save a copy
      local_kw['local_sharpen'] = False
      local_kw['overall_sharpen_before_and_after_local'] = False
      local_kw_dc = deepcopy(local_kw)

      local_kw = deepcopy(local_kw_dc)

      print ("\n",79*"=","\nRunning overall sharpening now\n",79*"=","\n",
         file = self.log)

      local_kw = self._update_kw_with_map_info(local_kw, previous_kw = kw,
        text = 'overall')

      method( **local_kw)  # run overall sharpening

      # Now our new maps to be sharpened are in local_kw['map_id_scaled_list']

      print ("\nDone with overall sharpening\n", file = self.log)


      # Local sharpening
      print ("\n",79*"=","\nRunning local sharpening\n",79*"=","\n",
         file = self.log)

      kw = self._update_kw_with_map_info(kw, previous_kw = local_kw,
        text = 'local', have_previous_scaled_data = True)


      method( **kw)  # Run local sharpening

      print ("\n",79*"=","\nRunning final overall sharpening\n",79*"=","\n",
         file = self.log)

      local_kw = deepcopy(local_kw_dc)

      local_kw = self._update_kw_with_map_info(local_kw, previous_kw = kw,
        text = 'final_overall', have_previous_scaled_data = True,
        map_id_scaled_list = final_map_id_scaled_list)

      method( **local_kw) # Run overall sharpening

      print("Scaled maps are '%s' "%(
        str(local_kw['map_id_scaled_list'])), file = self.log)

      scaled_map_id = local_kw['map_id_scaled_list'][
         local_kw['map_id_to_be_scaled_list'].index(local_kw['map_id'])]
      print("\nFinal sharpened map is in '%s' in '%s' " %(
         scaled_map_id, self.name), file = self.log)

      if local_kw.get('model_id') and self.get_model_by_id(local_kw['model_id']):
        cc = self.map_model_cc(model_id = local_kw['model_id'],
           map_id = scaled_map_id)
        print ("Current map-model CC for '%s': %.3f " %(scaled_map_id,cc),
           file = self.log)

      print ("\n",79*"=","\nDone with local and overall sharpening\n",79*"=",
           file = self.log)

  def set_map_id_lists(self,kw):
    """Set values of keywords for maps"""
    if kw.get('overall_sharpen_before_and_after_local'):
      kw['sharpen_all_maps'] = True
    if kw.get('map_id') is None:
      kw['map_id'] = 'map_manager'
    if kw.get('map_id_to_be_scaled_list') is None:
      kw['map_id_to_be_scaled_list'] = [kw['map_id']]
      if kw.get('sharpen_all_maps') and \
            kw.get('map_id_1') and kw.get('map_id_2'): # half-map sharpening
         kw['map_id_to_be_scaled_list'].append(kw['map_id_1'])
         kw['map_id_to_be_scaled_list'].append(kw['map_id_2'])
    if kw.get('map_id_scaled_list') is None:
      kw['map_id_scaled_list'] = []
      for id in kw['map_id_to_be_scaled_list']:
         kw['map_id_scaled_list'].append("%s_scaled" %(id))
    return kw

  def external_sharpen(self,
      map_id = 'map_manager',
      map_id_external_map = 'external_map',
      map_id_to_be_scaled_list = None,
      map_id_scaled_list = None,
      exclude_points_outside_density = None,
      minimum_boxes_inside_density = None,
      resolution = None,
      d_min = None,
      k_sol = None,
      b_sol = None,
      n_bins = None,
      n_boxes = None,
      core_box_size = None,
      box_cushion = None,
      smoothing_radius = None,
      local_sharpen = None,
      anisotropic_sharpen = None,
      expected_ssqr_list = None,
      expected_ssqr_list_rms = None,
      tlso_group_info = None,
      get_tls_from_u = None,
      overall_sharpen_before_and_after_local = False,
      get_scale_as_aniso_u = None,
      use_dv_weighting = None,
      n_direction_vectors = None,
      run_analyze_anisotropy = True,
      sharpen_all_maps = False,
      nproc = None,
    ):
    '''
     Scale map_id with scale factors identified from map_id vs
      map_id_external_map
     Changes the working map_manager

     resolution is nominal resolution of map
     d_min is minimum resolution to use in calculation of Fourier coefficients
    '''

    from libtbx import adopt_init_args
    kw_obj = group_args()
    adopt_init_args(kw_obj, locals())
    kw = kw_obj() # save calling parameters in kw as dict
    del kw['adopt_init_args'] # REQUIRED
    del kw['kw_obj']  # REQUIRED

    # Checks
    assert self.get_map_manager_by_id(map_id)
    assert self.get_map_manager_by_id(map_id_external_map)

    # Allow sharpening globally before and after local sharpening
    if local_sharpen and overall_sharpen_before_and_after_local:
      print ("\nRunning external sharpening (global, local, global)\n",
       file = self.log)
      return self._sharpen_overall_local_overall(kw = kw,
        method = self.external_sharpen)

    kw = self.set_map_id_lists(kw)

    print ("\nRunning external map sharpening ", file = self.log)

    kw['map_id_2'] = map_id_external_map
    kw['is_external_based'] = True
    kw['remove_overall_anisotropy'] = False # REQUIRED
    del kw['map_id_external_map']
    kw['model_map_ids_to_leave_as_is'] = [map_id_external_map] # do not remove aniso
    self._sharpen_map(**kw)

  def half_map_sharpen(self,
      map_id = 'map_manager',
      map_id_1 = 'map_manager_1',
      map_id_2 = 'map_manager_2',
      map_id_scaled_list = None,
      map_id_to_be_scaled_list = None,
      exclude_points_outside_density = None,
      minimum_boxes_inside_density = None,
      resolution = None,
      d_min = None,
      k_sol = None,
      b_sol = None,
      n_bins = None,
      n_boxes = None,
      core_box_size = None,
      box_cushion = None,
      smoothing_radius = None,
      rmsd = None,
      local_sharpen = None,
      anisotropic_sharpen = None,
      minimum_low_res_cc = None,
      get_scale_as_aniso_u = None,
      use_dv_weighting = None,
      n_direction_vectors = None,
      run_analyze_anisotropy = True,
      spectral_scaling = True,
      expected_rms_fc_list = None,
      expected_ssqr_list = None,
      expected_ssqr_list_rms = None,
      tlso_group_info = None,
      get_tls_from_u = None,
      model_id_for_rms_fc = None,
      replace_aniso_with_tls_equiv = None,
      max_abs_b = None,
      nproc = None,
      optimize_b_eff = None,
      equalize_power = None,
      overall_sharpen_before_and_after_local = False,
      get_tls_info_only = None,
      coordinate_shift_to_apply_before_tlso = None,
      sharpen_all_maps = False,
      remove_overall_anisotropy = True,
    ):
    '''
     Scale map_id with scale factors identified from map_id_1 vs map_id_2
     Changes the working map_manager unless map_id_scaled_list is set.

     max_abs_b applies if get_scale_as_aniso_u and anisotropic_sharpen and
        local_sharpen are set. It limits range of anisotropic B.  Default is
        100 at 4 A, proportional to resolution squared

     resolution is nominal resolution of map
     d_min is minimum resolution to use in calculation of Fourier coefficients
    '''

    from libtbx import adopt_init_args
    kw_obj = group_args()
    adopt_init_args(kw_obj, locals())
    kw = kw_obj() # save calling parameters in kw as dict
    del kw['adopt_init_args'] # REQUIRED
    del kw['kw_obj']  # REQUIRED

    # Checks
    assert self.get_map_manager_by_id(map_id)

    # Set what maps are going to be sharpened and new names
    kw = self.set_map_id_lists(kw)

    # Allow sharpening globally before and after local sharpening
    if local_sharpen and overall_sharpen_before_and_after_local:
      print ("\nRunning half-map sharpening (global, local, global)\n",
       file = self.log)
      return self._sharpen_overall_local_overall(kw = kw,
        method = self.half_map_sharpen)
    print ("\nRunning half-map sharpening\n", file = self.log)
    print("Scale factors will be identified using the "+
       "maps '%s' and '%s' in map_model_manager '%s'" %(
        kw['map_id_1'],kw['map_id_2'],  self.name), file = self.log)
    print("Maps to be scaled are '%s' in map_model_manager '%s'" %(
      str(kw['map_id_to_be_scaled_list']),self.name),file = self.log)
    print("Sharpened maps after half-map sharpening will be in "+
       "'%s' in map_model_manager '%s'" %(
      str(kw['map_id_scaled_list']),self.name),file = self.log)

    if tlso_group_info:  # convert to lists
      convert_tlso_group_info_to_lists(tlso_group_info)
      if kw['get_tls_from_u'] is None:
        kw['get_tls_from_u'] = True

    # Now get scaling from comparison of the two half-maps
    #  apply the scaling to map_id_to_be_scaled

    if get_tls_info_only:
      return self._sharpen_map(**kw)
    else:
      self._sharpen_map(**kw)

  def model_sharpen(self,
      map_id = 'map_manager',
      model_id = 'model',
      map_id_scaled_list = None,
      map_id_to_be_scaled_list = None,
      exclude_points_outside_density = True,
      minimum_boxes_inside_density = True,
      resolution = None,
      d_min = None,
      k_sol = None,
      b_sol = None,
      find_k_sol_b_sol = True,
      d_min_for_k_sol_b_sol = 6.,
      n_bins = None,
      n_boxes = None,
      core_box_size = None,
      box_cushion = None,
      smoothing_radius = None,
      rmsd = None,
      local_sharpen = None,
      anisotropic_sharpen = None,
      minimum_low_res_cc = 0.20,
      get_scale_as_aniso_u = None,
      use_dv_weighting = None,
      n_direction_vectors = None,
      run_analyze_anisotropy = True,
      spectral_scaling = True,
      expected_rms_fc_list = None,
      expected_ssqr_list = None,
      expected_ssqr_list_rms = None,
      tlso_group_info = None,
      get_tls_from_u = None,
      find_tls_from_model = None,
      model_id_for_rms_fc = None,
      replace_aniso_with_tls_equiv = None,
      max_abs_b = None,
      nproc = None,
      optimize_b_eff = None,
      equalize_power = None,
      map_id_model_map = 'model_map_for_scaling',
      optimize_with_model = None,
      overall_sharpen_before_and_after_local = False,
      mask_around_model = True,
      get_tls_info_only = None,
      coordinate_shift_to_apply_before_tlso = None,
      sharpen_all_maps = False,
      remove_overall_anisotropy = True,
      save_model_map = False,
    ):
    '''
     Scale map_id with scale factors identified from map_id vs model
     Changes the working map_manager unless map_id_scaled is set.

     max_abs_b applies if get_scale_as_aniso_u and anisotropic_sharpen and
        local_sharpen are set. It limits range of anisotropic B.  Default is
        100 at 4 A, proportional to resolution squared

     resolution is nominal resolution of map
     d_min is minimum resolution to use in calculation of Fourier coefficients

    '''

    from libtbx import adopt_init_args
    kw_obj = group_args()
    adopt_init_args(kw_obj, locals())
    kw = kw_obj() # save calling parameters in kw as dict
    del kw['adopt_init_args'] # REQUIRED
    del kw['kw_obj']  # REQUIRED

    # Checks
    assert self.get_map_manager_by_id(map_id)
    assert self.get_model_by_id(model_id)

    # Set what maps are going to be sharpened and new names
    kw = self.set_map_id_lists(kw)


    print ("\nRunning model-based sharpening ", file = self.log)
    if local_sharpen:
      print("Sharpening will be local",file = self.log)
    if anisotropic_sharpen:
      if get_scale_as_aniso_u:
        if replace_aniso_with_tls_equiv:
          print("Sharpening will be anisotropic and converted to TLS",
           file = self.log)
        else:
          print("Sharpening will be anisotropic and converted to aniso U",
           file = self.log)
      else:
        print("Sharpening will be anisotropic",file = self.log)

    print("Scale factors will be identified using the "+
       "map '%s' and model '%s' in map_model_manager '%s'" %(
        kw['map_id'], kw['model_id'], self.name), file = self.log)
    print("Map to be scaled is '%s' in map_model_manager '%s'" %(
      str(kw['map_id_to_be_scaled_list']),self.name),file = self.log)
    print("Scaled map will be in '%s' in map_model_manager '%s'" %(
      str(kw['map_id_scaled_list']),self.name),file = self.log)

    if model_id_for_rms_fc is None:
      kw['model_id_for_rms_fc'] = kw['model_id']

    if tlso_group_info:  # convert to lists
      convert_tlso_group_info_to_lists(tlso_group_info)
      if get_tls_from_u is None:
        get_tls_from_u = True
    elif find_tls_from_model:
      # If we are going to use TLS groups from the model, check them here
      if not self.get_model_by_id(model_id):
        raise Sorry("Need model for find_tls_from_model")
      if get_tls_from_u is None:
        get_tls_from_u = True
      tlso_group_info = get_tlso_group_info_from_model(
         self.get_model_by_id(model_id),
         nproc = nproc,
         log = self.log)
      kw['tlso_group_info'] = tlso_group_info

    # Allow sharpening globally before and after local sharpening
    if local_sharpen and overall_sharpen_before_and_after_local:
      return self._sharpen_overall_local_overall(kw = kw,
        method = self.model_sharpen)

    del kw['find_k_sol_b_sol']  # REQUIRED
    del kw['d_min_for_k_sol_b_sol']  # REQUIRED
    del kw['mask_around_model']  # REQUIRED
    del kw['model_id']  # REQUIRED
    del kw['map_id_model_map']  # REQUIRED
    del kw['optimize_with_model']  # REQUIRED
    del kw['find_tls_from_model']  # REQUIRED
    del kw['overall_sharpen_before_and_after_local']  # REQUIRED
    del kw['save_model_map']  # REQUIRED


    # Make a copy of this map_model manager so we can modify it without
    #  changing the original
    working_mmm = self.deep_copy()
    working_mmm.set_name('working_mmm')

    # Working resolution is resolution * d_min_ratio
    if d_min is None:
      d_min = working_mmm._get_d_min_from_resolution(resolution)
    print ("High-resolution limit: "+
      "%5.2f A based on nominal resolution of %5.2f A" %(
      d_min, resolution if resolution else working_mmm.resolution()),
      file = self.log)
    map_id_to_be_scaled = kw['map_id_to_be_scaled_list'][0]
    cc = self.map_model_cc(map_id=map_id_to_be_scaled, model_id=model_id)
    print ("Map-model CC before sharpening: %.3f " %(cc), file = self.log)

    map_coeffs = working_mmm.get_map_manager_by_id(
       map_id_to_be_scaled).map_as_fourier_coefficients( d_min = d_min)

    working_n_bins =working_mmm._set_n_bins(n_bins = n_bins,
      d_min = d_min, map_coeffs = map_coeffs,
      local_sharpen = local_sharpen)

    f_array = get_map_coeffs_as_fp_phi(map_coeffs, n_bins = working_n_bins,
        d_min = d_min).f_array

    # Generate map from model using existing possibly anisotropic B
    model=working_mmm.get_model_by_id(model_id)

    if find_k_sol_b_sol and (k_sol is None) and (b_sol is None):
      # Find k_sol and b_sol
      local_mmm = working_mmm.extract_all_maps_around_model(
        stay_inside_current_map = True)
      local_mmm.mask_all_maps_around_atoms(
         mask_atoms_atom_radius = 2.* d_min,
         soft_mask =True)
      d_min_for_k_sol_b_sol = max(d_min, d_min_for_k_sol_b_sol)
      kb_info = local_mmm.find_k_sol_b_sol(local_mmm.get_model_by_id(model_id),
        d_min = d_min_for_k_sol_b_sol,
        model_map_id = map_id_model_map,
        comparison_map_id = map_id)
      if kb_info is not None:
        print("\nOptimized k_sol=%.2f  b_sol=%.1f CC=%.3f (d_min= %.2f A)" %(
         kb_info.k_sol, kb_info.b_sol, kb_info.cc, d_min_for_k_sol_b_sol),
            file = self.log)
        k_sol = kb_info.k_sol
        b_sol = kb_info.b_sol
        kw['k_sol'] = k_sol
        kw['b_sol'] = b_sol

    working_mmm.generate_map(model=model,
       gridding=working_mmm.get_any_map_manager().map_data().all(),
       d_min=d_min,
       map_id = map_id_model_map,
       k_sol = k_sol,
       b_sol = b_sol)

    # Save unmodified copy of map to be scaled
    mm_dc_list = []
    unmasked_map_id_list = []
    for id in kw['map_id_to_be_scaled_list']:
      mm_dc_list.append(working_mmm.get_map_manager_by_id(
        map_id = id).deep_copy())
      unmasked_map_id_list.append("%s_original_map" %(id))

    if mask_around_model:
      # Mask maps around the model
      print("Masking all maps in '%s' around model, saving original as '%s'" %(
        self.name,str(unmasked_map_id_list)), file = self.log)
      working_mmm.mask_all_maps_around_atoms(
       soft_mask=True, mask_atoms_atom_radius = 2.* working_mmm.resolution())

    # Put in unmasked map to be scaled
    for id, mm_dc in zip(unmasked_map_id_list, mm_dc_list):
      working_mmm.add_map_manager_by_id(map_id = id,
      map_manager = mm_dc)

    kw['map_id_2'] = map_id_model_map
    kw['is_model_based'] = True
    kw['map_id_to_be_scaled_list'] = unmasked_map_id_list
    kw['model_map_ids_to_leave_as_is'] = [map_id_model_map] # do not remove aniso

    # Now get scaling from comparison of working_map_id_to_be_scaled and
    #   map_id_model_map, and
    #  apply the scaling to  unmasked_map_id_list
    if get_tls_info_only:
      return working_mmm._sharpen_map(**kw)
    else:
      working_mmm._sharpen_map(**kw)

    # And set this map_manager
    if save_model_map:
      mm = working_mmm.get_map_manager_by_id(map_id_model_map)
      if mm:
        print("Copying map '%s' in map_manager '%s' to '%s'" %(
          map_id_model_map,working_mmm.name,self.name),
          file = self.log)
        self.add_map_manager_by_id(
          map_manager = working_mmm.get_map_manager_by_id(map_id_model_map),
          map_id = map_id_model_map)
    for id in kw['map_id_scaled_list']:
      if working_mmm.get_map_manager_by_id(id):
        print("Copying map '%s' in map_manager '%s' to '%s'" %(
          id,working_mmm.name,self.name),
          file = self.log)
        self.add_map_manager_by_id(
          map_manager = working_mmm.get_map_manager_by_id(id),
          map_id = id)
        cc = self.map_model_cc(map_id=id,
           model_id=model_id)
        print ("Map-model CC after sharpening: %.3f " %(cc), file = self.log)
      else:
        print("No sharpened map obtained "+
          "in '%s' from map_model_manager '%s'" %(
          id,working_mmm.name), file = self.log)

  def _sharpen_map(self,
      map_id = 'map_manager',
      map_id_1 = 'map_manager_1',
      map_id_2 = 'map_manager_2',
      map_id_to_be_scaled_list = None,
      map_id_scaled_list = None,
      exclude_points_outside_density = None,
      minimum_boxes_inside_density = None,
      equalize_power = None,
      resolution = None,
      d_min = None,
      k_sol = None,
      b_sol = None,
      n_bins = None,
      n_boxes = None,
      core_box_size = None,
      box_cushion = None,
      smoothing_radius = None,
      rmsd = None,
      local_sharpen = None,
      nproc = None,
      optimize_b_eff = None,
      is_model_based = False,
      is_external_based = False,
      n_bins_default = 200,
      n_bins_default_local = 20,
      anisotropic_sharpen = None,
      minimum_low_res_cc = 0.20,
      get_scale_as_aniso_u = None,
      use_dv_weighting = None,
      n_direction_vectors = None,
      run_analyze_anisotropy = None,
      spectral_scaling = None,
      expected_rms_fc_list = None,
      expected_ssqr_list = None,
      expected_ssqr_list_rms = None,
      tlso_group_info = None,
      get_tls_from_u = None,
      model_id_for_rms_fc = None,
      replace_aniso_with_tls_equiv = None,
      max_abs_b = None,
      get_tls_info_only = None,
      coordinate_shift_to_apply_before_tlso = None,
      overall_sharpen_before_and_after_local = None, # ignored
      sharpen_all_maps = False,
      model_map_ids_to_leave_as_is = None,
      remove_overall_anisotropy = None,
    ):
    '''
     Scale map_id with scale factors identified from map_id_1 and map_id_2
     Changes the working map_manager unless map_id_scaled_list
      is set in which case the map goes there.

     If spectral_scaling : multiply scale factors by expected amplitude
       vs resolution unless  map_id_2 is supplied

     if local_sharpen, use local sharpening

     If is_model_based, assume that the map_id_2 is based on a model
     If is_external_based, assume map_id_2 is external

     If anisotropic sharpening, identify resolution dependence along
      principal axes of anisotropy and apply based on position in reciprocal
      space

     max_abs_b applies if get_scale_as_aniso_u and anisotropic_sharpen and
        local_sharpen are set. It limits range of anisotropic B.  Default is
        100 at 4 A, proportional to resolution squared

     resolution is nominal resolution of map
     d_min is minimum resolution to use in calculation of Fourier coefficients
    '''
    from libtbx import adopt_init_args
    kw_obj = group_args()
    adopt_init_args(kw_obj, locals())
    kw = kw_obj() # save calling parameters in kw as dict

    # Remove things that are not actually keywords
    del kw['adopt_init_args'] # REQUIRED
    del kw['kw_obj']  # REQUIRED

    # Remove keywords that are not going to be passed on
    del kw['local_sharpen']  # REQUIRED
    del kw['overall_sharpen_before_and_after_local']  # REQUIRED
    del kw['n_bins_default']  # REQUIRED
    del kw['n_bins_default_local']  # REQUIRED
    del kw['remove_overall_anisotropy']  # REQUIRED
    del kw['model_map_ids_to_leave_as_is']  # REQUIRED

    aniso_info = None

    # Checks
    assert self.get_map_manager_by_id(map_id)
    assert (
    self.get_map_manager_by_id(map_id_1) or
        is_model_based or (
        is_external_based) and self.get_map_manager_by_id(map_id_2))
    if get_tls_from_u and n_boxes == 1:
      raise Sorry("Cannot get TLS with n_boxes==1")

    # Set get_scale_as_aniso_u if not set already
    if get_scale_as_aniso_u is None:
      if replace_aniso_with_tls_equiv or tlso_group_info:
        get_scale_as_aniso_u = True
        print("Getting scale factors as aniso U values", file = self.log)

    # remove any extra models and maps to speed up boxing and not modify orig
    map_id_list = [map_id,map_id_1,map_id_2]+kw['map_id_to_be_scaled_list']
    working_mmm = self.get_map_model_manager_with_selected(
      map_id_list = map_id_list,
      model_id_list =
          [model_id_for_rms_fc] if is_model_based and self.get_model_by_id(
            model_id_for_rms_fc) else None,
      deep_copy=True)
    working_mmm.set_log(self.log)

    # Note the map that is going to be scaled
    assert kw['map_id_to_be_scaled_list']
    assert kw['map_id_scaled_list']
    assert len(kw['map_id_to_be_scaled_list']) == len(kw['map_id_scaled_list'])

    if exclude_points_outside_density:  # Set up mask
      kw['mask_id'] = self._generate_new_map_id(prefix = 'mask_around_density')
      print("\nSetting up soft mask around density as '%s' " %(kw['mask_id']),
         file = self.log)
      working_mmm.create_mask_around_density(
        soft_mask  = True,
        mask_id = kw['mask_id'],
        map_id = map_id)
    else:
      kw['mask_id'] = None

    print("\nMaps in '%s' to be scaled will come from '%s'" %(
       self.name,str(kw['map_id_to_be_scaled_list'])),
       file = self.log)
    print("Scaled maps in '%s' will go in '%s'" %(
       self.name,str(kw['map_id_scaled_list'])), file = self.log)

    if n_bins is None:
      n_bins = n_bins_default

    # Get basic info including minimum_resolution (cutoff for map_coeffs)
    setup_info = working_mmm._get_box_setup_info(map_id_1, map_id_2,
      resolution,
      d_min,
      n_boxes = n_boxes,
      core_box_size = core_box_size,
      smoothing_radius = smoothing_radius,
      )
    if spectral_scaling and (not expected_rms_fc_list):
        from cctbx.development.approx_amplitude_vs_resolution import \
          approx_amplitude_vs_resolution
        aavr = approx_amplitude_vs_resolution(
           d_min = setup_info.minimum_resolution,
           resolution = setup_info.resolution,
           n_bins = n_bins,
           k_sol = k_sol,
           b_sol = b_sol,
           map_model_manager = working_mmm,
           model = working_mmm.get_model_by_id(model_id=model_id_for_rms_fc),
           out = self.log)
        expected_rms_fc_list = aavr.get_target_scale_factors()
        n_bins_use = aavr.n_bins_use # number of bins to resolution

    if expected_ssqr_list_rms and not expected_ssqr_list:
        from cctbx.development.approx_amplitude_vs_resolution import \
          get_expected_ssqr_list
        expected_ssqr_list = get_expected_ssqr_list(
           d_min = setup_info.minimum_resolution,
           n_bins = n_bins,
           expected_ssqr_list_rms = expected_ssqr_list_rms,
           map_model_manager = working_mmm,
           out = self.log)

    resolution = setup_info.resolution
    working_mmm.set_resolution(resolution)
    d_min = setup_info.minimum_resolution
    print ("Nominal resolution of map: %.2f A  Minimum resolution: %.2f A" %(
        resolution,d_min),
      file = self.log)
    mask_all_maps_around_edges = True
    if mask_all_maps_around_edges:
      working_mmm.mask_all_maps_around_edges(soft_mask_radius=resolution)

    map_coeffs = working_mmm.get_map_manager_by_id(
       map_id).map_as_fourier_coefficients( d_min = d_min)

    if remove_overall_anisotropy:
      b_iso = 10 * d_min  # works best if some overall B remains
      aniso_b_cart = working_mmm.remove_anisotropy(map_id = map_id,
        d_min = d_min,
        b_iso = b_iso,
        remove_from_all_maps = True,
        model_map_ids_to_leave_as_is = model_map_ids_to_leave_as_is)
    else:
      aniso_b_cart = None
      b_iso = None
    kw['b_iso'] = b_iso
    kw['aniso_b_cart'] = aniso_b_cart

    n_bins =working_mmm._set_n_bins(n_bins = n_bins,
      d_min = d_min, map_coeffs = map_coeffs,
      local_sharpen = False)

    if local_sharpen:
      # Now local sharpening
      print("\nSetting up local sharpening ...\n",file = self.log)

      if kw['n_bins'] is None: # set it here for local
        kw['n_bins'] = n_bins_default_local

      if get_tls_info_only:
        return working_mmm._local_sharpen(**kw)
      else:
        working_mmm._local_sharpen(**kw)

      for id, previous_id in zip(
          kw['map_id_scaled_list'],
          kw['map_id_to_be_scaled_list']):
        sharpened_local_mm = working_mmm.get_map_manager_by_id(id)
        if sharpened_local_mm:
          # We're done. put map in map manager and return
          print("Saving local-sharpened map '%s' from map_model_manager '%s'" %(
           id,working_mmm.name) + " as '%s' in '%s' " %(
           id,self.name), file = self.log)
          sharpened_local_mm.name = working_mmm.get_map_manager_by_id(
             previous_id).file_name

          self.add_map_manager_by_id(map_manager = sharpened_local_mm,
            map_id = id)

          # Get anisotropy (overall B) before scaling
          aniso_info = self._get_aniso_before_and_after(d_min = d_min,
            map_id = id, previous_map_id = previous_id)
          print(aniso_info.text, file = self.log)

        else:
          print("No local-sharpened map obtained "+
            "in '%s' from map_model_manager '%s'" %(
            id,working_mmm.name), file = self.log)

      return

    # Here to run overall
    print ("\nRunning overall sharpening ", file = self.log)

    if anisotropic_sharpen:
       print ("Using anisotropic sharpening ",file = self.log)
       # get scale factors in 12 directions (or 6)
       direction_vectors = working_mmm._get_aniso_direction_vectors(map_id,
         n_direction_vectors = n_direction_vectors)
    else:
       direction_vectors = [None]
    target_scale_factors_list = []

    # Mask after getting direction vectors

    if direction_vectors:
      print("\nEstimating scale factors for %s direction_vectors" %(
        len(direction_vectors)), file = self.log)
    else:
      print("Estimating scale factors ", file = self.log)
    # Analyze spectrum in map_id (will apply it to map_id_to_be_scaled)
    scaling_group_info = working_mmm._get_weights_in_shells(n_bins,
        d_min,
        map_id = map_id,
        map_id_1 = map_id_1,
        map_id_2 = map_id_2,
        rmsd = rmsd,
        resolution = resolution, # nominal resolution
        optimize_b_eff = optimize_b_eff,
        equalize_power = equalize_power,
        is_model_based = is_model_based,
        is_external_based = is_external_based,
        direction_vectors = direction_vectors,
        minimum_low_res_cc = minimum_low_res_cc,
        get_scale_as_aniso_u = get_scale_as_aniso_u,
        use_dv_weighting = use_dv_weighting,
        n_direction_vectors = n_direction_vectors,
        run_analyze_anisotropy = run_analyze_anisotropy,
        expected_rms_fc_list = expected_rms_fc_list,
        expected_ssqr_list = expected_ssqr_list,
        expected_ssqr_list_rms = expected_ssqr_list_rms,
        tlso_group_info = tlso_group_info,
        model_id_for_rms_fc = model_id_for_rms_fc,
        replace_aniso_with_tls_equiv = replace_aniso_with_tls_equiv,
        )
    """
    scaling_group_info group_args object:
      direction_vectors: direction vectors dv for anisotropy calculations
      scaling_info_list: si (scaling_info) objects, one for each dv
        each si:  si.target_scale_factors   # scale factors vs sthol2
                  si.target_sthol2 # sthol2 values  d = 0.25/sthol2**0.5
                  si.d_min_list
                  si.cc_list
                  si.low_res_cc # low-res average
      ss_b_cart_as_u_cart: anisotropic part of overall correction factor
      uu_b_cart_as_u_cart: Estimated total fall-off relative to ideal
      overall_scale: radial part of overall correction factor
    """


    #  Now apply scaling to map_id_to_be_scaled
    for id, new_id in zip(
       kw['map_id_to_be_scaled_list'],kw['map_id_scaled_list']):
      map_coeffs_to_be_scaled = self.get_map_manager_by_id(id
         ).map_as_fourier_coefficients(d_min=d_min)

      # Apply the scale factors in shells
      print("\nApplying final scale factors in shells of "+
        "resolution to map '%s' to yield '%s'" %(id, new_id), file = self.log)
      if len(direction_vectors) > 1:
        print("Using %s direction vectors" %len(direction_vectors),
          file = self.log)
        assert len(scaling_group_info.scaling_info_list) == len(
         direction_vectors)

      print ("\nApplying scale factors directly",
         file = self.log)

      new_map_manager = working_mmm._apply_scale_factors_in_shells(
          map_coeffs_to_be_scaled,
          n_bins,
          d_min,
          scaling_group_info = scaling_group_info,
          direction_vectors = direction_vectors,
          aniso_b_cart = aniso_b_cart,
          b_iso = b_iso,
          )

      if not new_map_manager:
        print("Not applying scaling to '%s'" %(id) ,file = self.log)
      else: # usual
        # All done... Set map_manager now
        print ("Scaled map is '%s' in %s" %(new_id,self.name), file = self.log)
        new_map_manager.file_name = self.get_map_manager_by_id(id).file_name
        self.add_map_manager_by_id(map_manager = new_map_manager,
          map_id = new_id)
        # Get anisotropy (overall B) before scaling
        aniso_info = self._get_aniso_before_and_after(d_min = d_min,
            map_id = new_id, previous_map_id = id)



    scale_factor_info = group_args(
       value_list=[scaling_group_info],
       xyz_list = [None],
    )

    tls_info = self._analyze_aniso(scale_factor_info,
      everything_is_inside = True,
      aniso_b_cart = aniso_b_cart,
      b_iso = b_iso,
     )
    if aniso_info:
      print(aniso_info.text, file = self.log)

    return tls_info

  def _get_aniso_before_and_after(self, d_min = None,
    map_id = None, previous_map_id = None):
    """Calculate anisotropy of map before and after sharpening"""
    prev_b_cart = self._get_aniso_of_map(d_min = d_min,
      map_id = previous_map_id)
    new_b_cart = self._get_aniso_of_map(d_min = d_min,
      map_id = map_id)
    b_sharpen = tuple(flex.double(prev_b_cart) - flex.double(new_b_cart))

    from six.moves import StringIO
    f = StringIO()
    print("\nSummary of anisotropic scaling applied for %s" %(
      self.get_map_manager_by_id(map_id).file_name), file = f)
    if prev_b_cart is not None:
      print("Original B-cart:    (%.2f, %.2f, %.2f, %.2f, %.2f, %.2f)" %(
       tuple(prev_b_cart)), file = f)
    if new_b_cart is not None:
      print("New B-cart:         (%.2f, %.2f, %.2f, %.2f, %.2f, %.2f)" %(
       tuple(new_b_cart)), file = f)
    if b_sharpen and b_sharpen is not None:
      print("Effective B-sharpen:(%.2f, %.2f, %.2f, %.2f, %.2f, %.2f)" %(
       tuple(b_sharpen)), file = f)
      print("Effective average B-sharpen: %.2f A**2" %(
        flex.double(b_sharpen[:3]).min_max_mean().mean), file = f)

    result = group_args(
     group_args_type = 'aniso_before_and_after for %s' %(previous_map_id),
     text = f.getvalue(),
     prev_b_cart = prev_b_cart,
     new_b_cart = new_b_cart,
     b_sharpen = b_sharpen,)
    return result


  def remove_anisotropy(self,
        d_min = None,
        map_coeffs = None,
        aniso_b_cart = None,
        map_id = 'map_manager',
        map_ids = None,
        remove_from_all_maps = False,
        model_map_ids_to_leave_as_is  = None,
        b_iso = None):
   '''
   Remove anisotropy from map, optionally remove anisotropy specified by
    aniso_b_cart and b_iso
   '''
   assert map_coeffs or d_min or map_id
   from cctbx.maptbx.segment_and_split_map import map_coeffs_as_fp_phi
   from cctbx.maptbx.refine_sharpening import analyze_aniso_object

   if not model_map_ids_to_leave_as_is:
     model_map_ids_to_leave_as_is = []

   if not map_coeffs:
      assert self.get_map_manager_by_id(map_id)
      map_coeffs = self.get_map_manager_by_id(map_id
        ).map_as_fourier_coefficients(d_min = d_min)
      f_array,phases=map_coeffs_as_fp_phi(map_coeffs)
   if not d_min:
     d_min = map_coeffs.d_min()


   if (not aniso_b_cart):
     aniso_b_cart = self._get_aniso_of_map(d_min = d_min, map_id = map_id)


   if remove_from_all_maps:  # remove in place from all maps
     print("Removing anisotropy from all maps", file = self.log)
     # Note: do not remove anisotropy from model-based maps..
     self._print_overall_u(aniso_b_cart,b_iso)

     if map_ids is None:
       map_ids = list(self._map_dict.keys())
     for map_id in map_ids:
       mm=self.get_map_manager_by_id(map_id)
       if mm.is_mask():
         continue # do not apply to masks
       elif map_id in model_map_ids_to_leave_as_is:
         continue # skip model maps
       map_coeffs = mm.map_as_fourier_coefficients(d_min = d_min)
       f_array,phases=map_coeffs_as_fp_phi(map_coeffs)
       analyze_aniso = analyze_aniso_object()
       analyze_aniso.set_up_aniso_correction(f_array=f_array,
         b_iso = b_iso,
         d_min = d_min,
         b_cart_to_remove = aniso_b_cart)
       scaled_f_array = analyze_aniso.apply_aniso_correction(f_array=f_array)
       new_mm = self.map_manager(
          ).fourier_coefficients_as_map_manager(
           scaled_f_array.phase_transfer(phase_source=phases,
           deg=True))
       mm.set_map_data(map_data = new_mm.map_data())
       print("Removed anisotropy from map '%s' " %(map_id), file = self.log)
     return aniso_b_cart

   else:  # apply to f_array and return
     print("Removing anisotropy and returning as map",
          "   %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f  " %(
      tuple(aniso_b_cart)),file = self.log)

     analyze_aniso = analyze_aniso_object()
     analyze_aniso.set_up_aniso_correction(f_array=f_array,
         b_iso = b_iso,
          b_cart_to_remove=aniso_b_cart)
     scaled_f_array = analyze_aniso.apply_aniso_correction(f_array=f_array)
     scaled_f_array.set_observation_type( f_array)

     return self.map_manager(
        ).fourier_coefficients_as_map_manager(
         scaled_f_array.phase_transfer(phase_source=f_array_info.phases,
         deg=True))

  def _get_aniso_of_map(self, d_min = None, map_id = 'map_manager'):
    """Get anisotropy of a map"""
    if not d_min:
      d_min = self.resolution()
    mm=self.get_map_manager_by_id(map_id)
    map_coeffs = mm.map_as_fourier_coefficients(d_min = d_min)
    from cctbx.maptbx.segment_and_split_map import map_coeffs_as_fp_phi
    f_array,phases=map_coeffs_as_fp_phi(map_coeffs)
    from cctbx.maptbx.refine_sharpening import analyze_aniso_object
    analyze_aniso = analyze_aniso_object()
    analyze_aniso.set_up_aniso_correction(f_array=f_array, d_min = d_min)
    return analyze_aniso.b_cart

  def _get_aniso_direction_vectors(self, map_id, n_direction_vectors = None,
     orient_to_axes = True):
    '''
     Find principal components of anisotropy in map and generate direction
     vectors
          Y      XY
          YZ
          Z  XZ  X

     Unique set is X Y Z XY XZ YZ
     Additional is (X,-Y,0), (X, -Z, 0), (Y, -Z, 0)
     Diagonals are (X,Y,Z), (X,Y,-Z), (X,-Y,Z), (X,-Y, -Z)
    '''
    if n_direction_vectors is not None:
      n_max = n_direction_vectors  # up to 13 minimum of 3
    else:
      n_max = 9

    ev = flex.vec3_double()
    if orient_to_axes:  # X Y Z
      ev.append((1,0,0))
      ev.append((0,1,0))
      ev.append((0,0,1))
    else:
      assert self.get_map_manager_by_id(map_id)
      map_coeffs = self.get_map_manager_by_id(map_id
        ).map_as_fourier_coefficients()
      f_array = get_map_coeffs_as_fp_phi(map_coeffs, n_bins = 1,
          d_min = self.resolution()).f_array
      from cctbx.maptbx.segment_and_split_map import get_b_iso
      b_mean,aniso_scale_and_b=get_b_iso(f_array,d_min=self.resolution(),
        return_aniso_scale_and_b=True)
      for i in range(3):
          if i >= n_max : break
          ev.append(tuple((
              aniso_scale_and_b.eigen_vectors[3*i],
              aniso_scale_and_b.eigen_vectors[3*i+1],
              aniso_scale_and_b.eigen_vectors[3*i+2])))

    if n_max >= 6:  # XY XZ YZ
      ev.append( col(ev[0]) + col(ev[1]) )
      ev.append( col(ev[0]) + col(ev[2]) )
      ev.append( col(ev[1]) + col(ev[2]) )
    if n_max >= 9:
      ev.append( col(ev[0]) - col(ev[1]) )
      ev.append( col(ev[0]) - col(ev[2]) )
      ev.append( col(ev[1]) - col(ev[2]) )
    if n_max >= 13:
      # Now add vectors between these in case that is where the variation is
      ev.append( col(ev[0]) + col(ev[1]) + col(ev[2]))
      ev.append(-col(ev[0]) + col(ev[1]) + col(ev[2]))
      ev.append( col(ev[0]) - col(ev[1]) + col(ev[2]))
      ev.append( col(ev[0]) + col(ev[1]) - col(ev[2]))
    norms = ev.norms()
    norms.set_selected((norms == 0),1)
    ev = ev/norms
    return ev

  def _set_default_parameters(self, other, name = None):
    """Set default parameters"""

    other._resolution = self._resolution
    other.set_log(self.log)
    other._nproc = self._nproc
    other._multiprocessing = self._multiprocessing
    other._queue_run_command = self._queue_run_command
    other._force_wrapping = deepcopy(self._force_wrapping)
    other._warning_message = self._warning_message
    other._stop_file = self._stop_file

    other.set_log(self.log)

    # Set up name
    if name is None:
      name = '%s_copy' %(self.name)
    other.set_name(name)

    other.set_verbose(self.verbose)
    if self._resolution:
      other.set_resolution(self._resolution)
    if self._minimum_resolution:
      other.set_minimum_resolution(self._minimum_resolution)
    if self._scattering_table:
      other.set_scattering_table(self._scattering_table)



  def get_map_model_manager_with_selected(self,
      map_id_list=None, model_id_list = None,
      deep_copy = False):
    """Create a new map_model_manager with just what we need"""
    assert map_id_list # Need maps to create map_model_manager with selected
    working_mmm = map_model_manager(
      map_manager = self.get_any_map_manager(), log = self.log,
      verbose = self.verbose)
    self._set_default_parameters(working_mmm)


    if map_id_list:
      already_copied = []
      for id in map_id_list:
        if id in already_copied: continue
        already_copied.append(id)
        if self.get_map_manager_by_id(id):
          working_mmm.add_map_manager_by_id(self.get_map_manager_by_id(id),id)
    if model_id_list:
      already_copied = []
      for id in model_id_list:
        if id in already_copied: continue
        already_copied.append(id)
        if self.get_model_by_id(id):
          working_mmm.add_model_by_id(self.get_model_by_id(id),id)
    if deep_copy:
      working_mmm = working_mmm.deep_copy()
      working_mmm.set_name("%s_with_selected_deep_copy" %(self.name))
    else:
      working_mmm.set_name("%s_with_selected" %(self.name))
    return working_mmm

  def _set_n_bins(self, n_bins = None,
      d_min = None, map_coeffs = None,
      local_sharpen = None):
    """Set number of resolution bins"""

    if n_bins is None:
      if local_sharpen:
        n_bins = 20
      else:
        n_bins = 200

    min_n_bins = 3
    original_n_bins = n_bins
    while n_bins > min_n_bins:
      f_array = get_map_coeffs_as_fp_phi(map_coeffs, n_bins = n_bins,
        d_min = d_min).f_array
      failed = False

      for i_bin in f_array.binner().range_used():
        if f_array.binner().count(i_bin)<1:
          failed = True
          break
      if failed:
        n_bins -= 1
      else: # ok
        return n_bins
    raise Sorry("Unable to set n_bins... possibly map too small?")

  def _apply_scale_factors_in_shells(self,
      map_coeffs,
      n_bins,
      d_min,
      target_scale_factors = None,
      scaling_group_info = None,
      direction_vectors= None,
      aniso_b_cart = None,
      b_iso = None,
      ):

    """Apply scale factors in shells of resolution"""
    f_array_info = get_map_coeffs_as_fp_phi(map_coeffs, n_bins = n_bins,
        d_min = d_min)

    if (not direction_vectors) or (direction_vectors[0] is None):
      target_scale_factors = scaling_group_info.scaling_info_list[0
         ].target_scale_factors
    if not target_scale_factors and (
        (not scaling_group_info) or not (scaling_group_info.scaling_info_list)
         or len(scaling_group_info.scaling_info_list)==1):
      print("Unable to scale dataset",file = self.log)
      return None

    if aniso_b_cart: # first apply aniso_b_cart to the whole array
       apply_aniso_b_cart_to_f_array_info(f_array_info,
         b_iso, d_min, aniso_b_cart)


    if target_scale_factors: # usual
      assert target_scale_factors.size() == n_bins
      assert len(list(f_array_info.f_array.binner().range_used())) == \
          target_scale_factors.size() # must be compatible binners
      scale_array=f_array_info.f_array.binner().interpolate(
        target_scale_factors, 1) # d_star_power=1
      scaled_f_array=f_array_info.f_array.customized_copy(
          data=f_array_info.f_array.data()*scale_array)

    else:  # apply anisotropic values
      assert scaling_group_info.scaling_info_list and direction_vectors
      assert len(scaling_group_info.scaling_info_list)==1 or \
         len(scaling_group_info.scaling_info_list) == direction_vectors.size()
      scale_array=flex.double(f_array_info.f_array.size(),0.)
      scale_array_weights=flex.double(f_array_info.f_array.size(),0.)
      from cctbx.maptbx.refine_sharpening import get_weights_para
      for si,direction_vector in zip(
           scaling_group_info.scaling_info_list,direction_vectors):
        if not si.target_scale_factors:
          continue  # just skip it
        assert si.target_scale_factors.size() == n_bins
        assert len(list(f_array_info.f_array.binner().range_used())) == \
          si.target_scale_factors.size() # must be compatible binners
        working_scale_array=f_array_info.f_array.binner().interpolate(
          si.target_scale_factors, 1) # d_star_power=1
        weights = get_weights_para(f_array_info.f_array, direction_vector,
          weight_by_cos = False,
          min_dot = 0.7,
          include_all_in_lowest_bin = True)
        scale_array += working_scale_array * weights
        scale_array_weights += weights
      scale_array_weights.set_selected((scale_array_weights < 1.e-10),1.e-10)
      scale_array /= scale_array_weights

      scaled_f_array=f_array_info.f_array.customized_copy(
          data=f_array_info.f_array.data()*scale_array)


    return self.map_manager(
       ).fourier_coefficients_as_map_manager(
         scaled_f_array.phase_transfer(phase_source=f_array_info.phases,
         deg=True))

  def _get_weights_in_shells(self,
     n_bins,
     d_min,
     map_id = 'map_manager',
     map_id_1 = 'map_manager_1',
     map_id_2 = 'map_manager_2',
     scale_using_last = 3,
     resolution = None,  # nominal resolution
     rmsd = None,
     cc_cut = 0.2,
     max_cc_for_rescale = 0.2,
     pseudo_likelihood = None,
     equalize_power = None,
     optimize_b_eff = None,
     is_model_based = None,
     is_external_based = None,
     minimum_low_res_cc = 0.20,
     direction_vectors = None,
     get_scale_as_aniso_u = None,
     use_dv_weighting = None,
     n_direction_vectors = None,
     run_analyze_anisotropy = None,
     expected_rms_fc_list = None,
     expected_ssqr_list = None,
     expected_ssqr_list_rms = None,
     tlso_group_info = None,
     model_id_for_rms_fc = None,
     replace_aniso_with_tls_equiv = None,
     ):
    '''
    Calculate weights in shells to yield optimal final map .
    If equalize_power, assume that perfect map has uniform power in all shells
    n_bins and d_min are required
    '''

    # Defaults:

    if not direction_vectors:
      direction_vectors = [None]

    if equalize_power is None:
      equalize_power = True

    if optimize_b_eff is None:
      if is_model_based:
        optimize_b_eff = False
      else:
        optimize_b_eff = False
    si = group_args(
      target_scale_factors = None,
      b_sharpen = 0,
      b_iso = 0,
      verbose = None,
      rmsd = rmsd,
      n_bins = n_bins,
      resolution = d_min,
      cc_cut = cc_cut,
      scale_using_last = scale_using_last,
      max_cc_for_rescale = max_cc_for_rescale,
      pseudo_likelihood = pseudo_likelihood,
      equalize_power = equalize_power,
      n_real = self.map_data().all(),
     )

    if self.get_map_manager_by_id(map_id):
      map_coeffs = self.get_map_manager_by_id(map_id
         ).map_as_fourier_coefficients(d_min=d_min)
    else:
      map_coeffs = None
    if self.get_map_manager_by_id(map_id_1):
      first_half_map_coeffs = self.get_map_manager_by_id(map_id_1
          ).map_as_fourier_coefficients(d_min=d_min)
    else:
      first_half_map_coeffs = None
    if self.get_map_manager_by_id(map_id_2):
      second_half_map_coeffs = self.get_map_manager_by_id(map_id_2
          ).map_as_fourier_coefficients(d_min=d_min)
    else:
      second_half_map_coeffs = None

    from cctbx.maptbx.refine_sharpening import calculate_fsc
    f_array = get_map_coeffs_as_fp_phi(map_coeffs, n_bins = n_bins,
        d_min = d_min).f_array
    ok_bins = True
    for i_bin in f_array.binner().range_used():
      if f_array.binner().count(i_bin)<1: # won't work...skip
        return None
    if is_external_based:
      external_map_coeffs = second_half_map_coeffs
      first_half_map_coeffs = None
      second_half_map_coeffs = None
      model_map_coeffs = None
    elif is_model_based:
      model_map_coeffs = second_half_map_coeffs
      first_half_map_coeffs = None
      second_half_map_coeffs = None
      external_map_coeffs = None
    else: # half-map
      external_map_coeffs = None
      model_map_coeffs = None

    result = calculate_fsc(
      f_array = f_array,
      map_coeffs = map_coeffs,
      first_half_map_coeffs = first_half_map_coeffs,
      second_half_map_coeffs = second_half_map_coeffs,
      model_map_coeffs=model_map_coeffs,
      external_map_coeffs=external_map_coeffs,
      si = si,
      cc_cut = si.cc_cut,
      optimize_b_eff = optimize_b_eff,
      is_model_based = is_model_based,
      scale_using_last=si.scale_using_last,
      max_cc_for_rescale=si.max_cc_for_rescale,
      pseudo_likelihood=si.pseudo_likelihood,
      equalize_power = si.equalize_power,
      direction_vectors = direction_vectors,
      smooth_fsc = False, # XXX may change
      cutoff_after_last_high_point = True,
      use_dv_weighting = use_dv_weighting,
      run_analyze_anisotropy = run_analyze_anisotropy,
      expected_rms_fc_list = expected_rms_fc_list,
      expected_ssqr_list = expected_ssqr_list,
      expected_ssqr_list_rms = expected_ssqr_list_rms,
      tlso_group_info = tlso_group_info,
      resolution = resolution, # nominal resolution
      remove_anisotropy_before_analysis = True, # work with aniso-removed
      out = self.log)
    if not hasattr(result,'scaling_info_list'):  # result is one si
      result = group_args(
        overall_si = result,
        scaling_info_list = [result],
        direction_vectors = direction_vectors,
        expected_rms_fc_list = expected_rms_fc_list,)

    # Set anything with too-low low-res CC to None for model-based run
    for si in result.scaling_info_list:
      if is_model_based and si.low_res_cc < minimum_low_res_cc:
        print("Skipping scaling with low_res_cc = %.2f" %(si.low_res_cc),
          file = self.log)
        si.target_scale_factors = None

    """
    scaling_group_info group_args object:
      direction_vectors: direction vectors dv for anisotropy calculations
      scaling_info_list: si (scaling_info) objects, one for each dv
        each si:  si.target_scale_factors   # scale factors vs sthol2
        si.target_sthol2 # sthol2 values  d = 0.25/sthol2**0.5
                  si.d_min_list
                  si.cc_list
                  si.low_res_cc # low-res average
      ss_b_cart_as_u_cart: anisotropic part of overall correction factor
      uu_b_cart_as_u_cart: Estimated total fall-off relative to ideal
      overall_scale: radial part of overall correction factor
    """
    return result

  def _remove_temp_dir(self,temp_dir):
    """Remove temporary directory"""
    if not os.path.isdir(temp_dir):
      return  # nothing to do
    else:  # remove it
     try:
       from shutil import rmtree
       rmtree(temp_dir)
     except Exception as e:
       pass # must have been removed another way

  def _create_temp_dir(self, temp_dir):
    """Create temporary directory"""
    if not os.path.isdir(temp_dir):
      os.mkdir(temp_dir)
      return temp_dir
    else:
      for i in range(1000):
        work_dir = "%s_%s" %(temp_dir,i)
        if not os.path.isdir(work_dir):
          os.mkdir(work_dir)
          return work_dir
    raise Sorry("Unable to create temporary directory", file = self.log)

  def _update_scale_factor_info_from_aniso(self, scale_factor_info,
      max_abs_b = None, get_tls_from_u = None):
    """Update scale factor information from aniso U values"""
    if max_abs_b is None:  # set default
      max_abs_b = 100 * (self.resolution()/2.5)**2  # about 100 at 2.5 A
    print("\nUpdating scale factor info from aniso U values.  \n"+
      "Maximum B correction allowed: %.2f" %(
     max_abs_b),file = self.log)
    if get_tls_from_u:
      print("\nSupplied TLS assumed to be overall fall-off with resolution",
        "\nScale factor will be reduced by 1/(1+E**2)",
           file = self.log)
    else:
      print("\nSupplied TLS assumed to be desired anisotropy of scale factors",
           file = self.log)

    map_coeffs = self.get_any_map_manager(
         ).map_as_fourier_coefficients(d_min = self.resolution())
    from cctbx.maptbx.segment_and_split_map import map_coeffs_as_fp_phi
    from cctbx.maptbx.refine_sharpening import get_nearest_lattice_points
    f_array,phases=map_coeffs_as_fp_phi(map_coeffs)
    if scale_factor_info.value_list and \
       scale_factor_info.value_list[0].overall_scale and (not f_array.binner()):
      f_array.setup_binner(
       n_bins = scale_factor_info.value_list[0].overall_scale.size(),
       d_min = scale_factor_info.d_min)

    xyz_list = scale_factor_info.xyz_list
    value_list = scale_factor_info.value_list
    for xyz, scaling_group_info in zip(xyz_list,value_list):
      direction_vectors = scaling_group_info.direction_vectors
      scaling_info_list = scaling_group_info.scaling_info_list
      if not scaling_group_info.overall_scale:
        continue # missing data
      elif get_tls_from_u:
        u_cart_to_apply = tuple(-flex.double(
          scaling_group_info.uu_b_cart_as_u_cart))
        overall_scale = None
        scale_by_ssqr_plus_one = True
      else: # usual
        # if multiply overall by interpolated scale factor
        u_cart_to_apply = tuple(flex.double(
          scaling_group_info.ss_b_cart_as_u_cart))
        overall_scale = scaling_group_info.overall_scale
        scale_by_ssqr_plus_one = False

      for si,dv in zip(scaling_info_list,direction_vectors):
        tsf = si.target_scale_factors
        target_sthol2 = si.target_sthol2
        # Recalculate target_scale_factors from ss_b_cart_as_u_cart and dv

        # NOTE: the "reflections" here are just values along a line in
        #  reciprocal space

        recip_space_vectors = flex.vec3_double()
        scale_values = flex.double()
        for sthol2 in target_sthol2:
          s = (sthol2/0.25)**0.5
          recip_space_vectors.append(col(dv) * s)
          scale_values.append(1.)

        # Create dummy array with fine spacing (will never contain much)
        #  so that there will be a lattice point very near any point in s-space
        local_f_array = create_fine_spacing_array(
          f_array.crystal_symmetry().unit_cell())

        indices = flex.miller_index(tuple(get_nearest_lattice_points(
            local_f_array.unit_cell(),recip_space_vectors)))
        scale_values_array = local_f_array.customized_copy(
              # just n_bins numbers
          data = scale_values,
          indices = indices)
        from mmtbx.scaling import absolute_scaling
        b_iso = abs(flex.double(adptbx.u_as_b(
            u_cart_to_apply))[:3].min_max_mean().mean)

        if max_abs_b and abs(b_iso) > max_abs_b:
          ratio = max_abs_b/abs(b_iso)
          u_cart_to_apply = tuple(ratio * col(
              u_cart_to_apply))
        u_star= adptbx.u_cart_as_u_star(
          scale_values_array.unit_cell(), tuple(-col(u_cart_to_apply)))
        scaled_f_array = absolute_scaling.anisotropic_correction(
          scale_values_array,0.0, u_star ,must_be_greater_than=-0.0001)

        # Now extract our values as target_scale_factors
        si.target_scale_factors = scaled_f_array.data()

        if overall_scale:  # Apply overall
          #resolution-dependent scale factors
          si.target_scale_factors *= scaling_group_info.overall_scale

        if scale_by_ssqr_plus_one:
          # Estimate best scale factor given errors
          si.target_scale_factors *= 1./(1.+si.ssqr_values)

    print("Done updating scale factors from aniso u values\n", file = self.log)

  def _analyze_aniso_replace_with_supplied(self,
     scale_factor_info,
     tlso_group_info = None,
     map_id = None,
     model_id = None,
     get_tls_from_u = None,
     aniso_b_cart = None,
     b_iso = None,
    ):
    """Analyze anisotropy and optionally replace with supplied information"""
    if model_id is None:
      model_id = 'model'

    if get_tls_from_u is None:
      get_tls_from_u = True

    # Replace aniso information with values from tlso_group_info
    #  Assume TLS is overall fall-off of data if get_tls_from_u is True
    #    set uu_b_cart_as_u_cart and ss_b_cart_as_u_cart

    #  Otherwise assume TLS is anisotropy of data
    #    set ss_b_cart_as_u_cart only

    print("Updating scale factor info from supplied aniso U values.",
      file = self.log)
    self._print_overall_u(aniso_b_cart,b_iso)


    # Get uaniso in middle of molecule
    tlso_value = tlso(
           t = tlso_group_info.T_list[0],
           l = tlso_group_info.L_list[0],
           s = tlso_group_info.S_list[0],
           origin = tlso_group_info.O_list[0],
             )
    sites_cart = flex.vec3_double()
    sites_cart.append(scale_factor_info.xyz_list.mean())

    center_uanisos= uaniso_from_tls_one_group(
         tlso = tlso_value,
         sites_cart = sites_cart,
         zeroize_trace=False)

    center_overall_u_cart = tuple(-flex.double(center_uanisos[0]) )


    # Get a standard values so we can edit it later
    # Could be better to use closest value
    average_scale_factor_info = self._average_scale_factor_info_over_xyz(
          scale_factor_info)

    std_values = average_scale_factor_info.value_list[0]
    assert std_values.scaling_info_list
    assert std_values.scaling_info_list[0].target_scale_factors

    print("\nUsing values from supplied TLS to set target scale factors",
        file = self.log)
    if get_tls_from_u:
      print ("\nSupplied TLS assumed to be overall fall-off with resolution"+
       "\nCorrection for uncertainties C = 1/(1+E**2) will be applied",
         file = self.log)
    else:
      print("\nSupplied TLS assumed to be desired anisotropy of scale factors",
         file = self.log)
    for T,L,S,origin,selection,other_shift_cart in zip(
       tlso_group_info.T_list,
       tlso_group_info.L_list,
       tlso_group_info.S_list,
       tlso_group_info.O_list,
       tlso_group_info.tlso_selection_list,
       tlso_group_info.tlso_shift_cart_list,):
      # Get mask representing this TLS
      mask_map_manager = self._create_mask_from_selection_as_string(
        map_id = map_id,
        model_id = model_id,
        selection_string=selection)

      working_scale_factor_info = self._get_scale_factor_info_inside_mask(
         scale_factor_info,
         mask_map_manager = mask_map_manager,
         inside = True)

      xyz_list = working_scale_factor_info.xyz_list
      print("\nAdding info from tls for %s" %(selection), file = self.log)

      if xyz_list.size() >= 1:
        print("Total of %s grid points inside this selection" %(
          xyz_list.size()), file = self.log)

        value_list = working_scale_factor_info.value_list
        if other_shift_cart:
          xyz_list_use = xyz_list + tuple([sc - other_sc for sc,other_sc in zip(
            self.shift_cart(), other_shift_cart)]) # XXX CHECK
        else:
          xyz_list_use = xyz_list

        print("Using TLS information from tlso_group_info", file = self.log)
        tlso_value = tlso(
             t = tlso_group_info.T_list[0],
             l = tlso_group_info.L_list[0],
             s = tlso_group_info.S_list[0],
             origin = tlso_group_info.O_list[0],
               )

        anisos_from_tls = uaniso_from_tls_one_group(
           tlso = tlso_value,
           sites_cart = xyz_list_use,
           zeroize_trace=False)

        for i in range(xyz_list.size()):

          u_cart_from_tls = tuple(flex.double(anisos_from_tls[i]))
          scaling_group_info = working_scale_factor_info.value_list[i]
          u_cart_from_tls = self._remove_overall_from_u_cart(u_cart_from_tls,
             aniso_b_cart,b_iso)

          if get_tls_from_u:
            scaling_group_info.uu_b_cart_as_u_cart = tuple(
              flex.double(u_cart_from_tls))
            tr = flex.double(u_cart_from_tls)[:3].min_max_mean().mean
            scaling_group_info.ss_b_cart_as_u_cart = tuple(
              [- (u - tr)  for u in u_cart_from_tls]) # MINUS

          else:  # usual
            scaling_group_info.ss_b_cart_as_u_cart = tuple(
              -flex.double(u_cart_from_tls)) # MINUS

          if self.verbose:
            print("\n   XYZ",xyz_list[i],
              file = self.log)
            print("   U from TLS:%s"  %str(u_cart_from_tls),
              file = self.log)
            print("   S aniso :  %s" %str(
              scaling_group_info.ss_b_cart_as_u_cart),
              file = self.log)

        if self.verbose:
          print("\nMean anisotropy as TLS:",file = self.log)
          print("T: %s" %(str(T)),file = self.log)
          print("L: %s" %(str(L)), file = self.log)
          print("S: %s" %(str(S)), file = self.log)

      else:
        print("Total of %s grid points inside this selection" %(
          xyz_list.size()), file = self.log)

      # Tack on values from points near model
      #  Make sure there are a number about equiv to the number of boxes
      #  that would be inside the model
      if self.get_model_by_id(model_id):
        new_xyz_list = self.get_model_by_id(
           model_id).apply_selection_string(selection).get_sites_cart()
        number_to_add = max(25,len(scale_factor_info.xyz_list))
        number_available = new_xyz_list.size()
        i_ratio = max(1,int((1+new_xyz_list.size())/number_to_add))
        if i_ratio > 1:
          short_xyz_list=flex.vec3_double()
          for i in range(0,number_available,i_ratio):
            short_xyz_list.append(new_xyz_list[i])
          new_xyz_list = short_xyz_list
        print("\nAdding %s u_cart values from TLS and model directly" %(
           new_xyz_list.size()), file = self.log)
        new_anisos_from_tls = uaniso_from_tls_one_group(
           tlso = tlso_value,
           sites_cart = new_xyz_list,
           zeroize_trace=False)
        for xyz,aniso in zip(new_xyz_list,new_anisos_from_tls):

          u_cart_from_tls = tuple(flex.double(aniso))
          u_cart_from_tls = self._remove_overall_from_u_cart(u_cart_from_tls,
              aniso_b_cart,b_iso)

          scaling_group_info = deepcopy(std_values)

          if get_tls_from_u:
            scaling_group_info.uu_b_cart_as_u_cart = tuple(
              flex.double(u_cart_from_tls))
            tr = flex.double(u_cart_from_tls)[:3].min_max_mean().mean
            scaling_group_info.ss_b_cart_as_u_cart = tuple(
              [- (u - tr)  for u in u_cart_from_tls]) # MINUS

          else:  # usual
            scaling_group_info.ss_b_cart_as_u_cart = tuple(
              -flex.double(u_cart_from_tls)) # MINUS

          scale_factor_info.value_list.append(scaling_group_info)
          scale_factor_info.xyz_list.append(xyz)
          if self.verbose:
            print("\n   XYZ from sites ",xyz,
              file = self.log)
            print("   U from TLS:%s"  %str(u_cart_from_tls),
              file = self.log)
            print("   S aniso :  %s" %str(
             scaling_group_info.ss_b_cart_as_u_cart),
              file = self.log)
        print("\nDone adding %s u_cart values from TLS and model directly" %(
           new_xyz_list.size()), file = self.log)

    tls_info = group_args(
         tlso = None,
         default_aniso = center_overall_u_cart
     )

    print("Done updating scale factor info from aniso U values.",
      file = self.log)
    return tls_info

  def _create_mask_from_selection_as_string(self,
      map_id = None,
      model_id = None,
      selection_string= None):
      '''
        Create a mask around density corresponding to atoms selected by
       selection_string.  If no model, just create a mask around all density
      '''
      if model_id is None:
        model_id = 'model'

      mask_id = self._generate_new_map_id(prefix = 'mask_around_density')
      if self.get_model_by_id(model_id):
        print("Creating mask around selection '%s' " %(selection_string),
          file = self.log)
        # Make a mask around the selected atoms
        self.create_mask_around_atoms(
          model = self.get_model_by_id(model_id
             ).apply_selection_string(selection_string),
          mask_atoms_atom_radius = 5,
         mask_id = mask_id)
        # multiply by density
        self.get_map_manager_by_id(map_id=mask_id).set_map_data(
          self.get_map_manager_by_id(map_id=mask_id).map_data() *
          self.get_map_manager_by_id(map_id=map_id).map_data() )
        old_mask_id = mask_id
        mask_id=old_mask_id+"_masked"
        self.create_mask_around_density(
          soft_mask  = True,
          mask_id = mask_id,
          map_id = old_mask_id)
        self.get_map_manager_by_id(map_id=mask_id).write_map('new_mask.ccp4')
      else:  # just make mask around density
        print("Creating mask around density",
          file = self.log)
        self.create_mask_around_density(
          soft_mask  = True,
          mask_id = mask_id,
          map_id = map_id)
        self.get_map_manager_by_id(map_id=mask_id).write_map('new_mask.ccp4')

      mask_map_manager = self.get_map_manager_by_id(map_id = mask_id)
      return mask_map_manager

  def _analyze_aniso(self,
     scale_factor_info,
     tlso_group_info = None,
     get_tls_from_u = None,
     map_id = None,
     mask_id = None,
     replace_inside = None,
     replace_boundary = None,
     replace_outside = None,
     coordinate_shift_to_apply_before_tlso = None,
     require_positive_definite = False,
     everything_is_inside = False,
     aniso_b_cart = None,
     b_iso = None,

    ):

    '''
      Summarize the anisotropy in the data and errors
      Optionally replace values inside mask with values from TLS analysis
      If get_tls_from_u then get tls from fo, not anisotropy correction s
      Optionally replace values with values from tlso_group_info

      Information is contained in the scale_factor_info object
      It can be location-specific or only overall
      It can be direction-specific or not.
    '''


    # Apply external values if supplied
    if tlso_group_info:
       result = self._analyze_aniso_replace_with_supplied(
         scale_factor_info,
         map_id = map_id,
         tlso_group_info = tlso_group_info,
         get_tls_from_u = get_tls_from_u,
         aniso_b_cart = aniso_b_cart,
         b_iso = b_iso,
         )
       self._summarize_scale_factor_info(scale_factor_info,
         aniso_b_cart = aniso_b_cart,
         b_iso = b_iso,)
       return result

    # Get a mask around the map if not already supplied as mask_id
    if (not everything_is_inside):
      if (not mask_id) or (not self.get_map_manager_by_id(mask_id)):
        mask_id = self._generate_new_map_id(prefix = 'mask_around_density')
        self.create_mask_around_density(
          soft_mask  = True,
          mask_id = mask_id,
          map_id = map_id)
      mask_map_manager = self.get_map_manager_by_id(mask_id)
    else:
      mask_map_manager = None

    tls_info = group_args(
       tlso = None,
       default_uaniso = None,
     )

    # Analyze scaling_group_info in relation to mask

    tlso_info_by_region={}
    inside_dict={True:'Inside mask',False:'Outside mask',None:'Edge of mask'}
    if everything_is_inside:
      inside_list = [True]
    else:
      inside_list = [True,False,None]

    for inside in inside_list:

      tlso_info_by_region[inside] = group_args(
        group_args_type = 'tlso_info_by_region',
        tlso = None,
        mean_u_cart = flex.double((0,0,0,0,0,0,)),
        mean_u_cart_n = 0,
        )

    # Now analyze by region

    for inside in inside_list:

      working_scale_factor_info = self._get_scale_factor_info_inside_mask(
         scale_factor_info,
         mask_map_manager = mask_map_manager,
         inside = inside,
         everything_is_inside = (everything_is_inside or (
           len(inside_list)==1 and len(scale_factor_info.xyz_list)==1)),
        )
      if working_scale_factor_info.xyz_list.size() < 1: continue

      scaling_group_info = working_scale_factor_info.value_list[0]
      direction_vectors = scaling_group_info.direction_vectors

      print(
       "\nLocal anisotropy and uncertainties by XYZ with inside_mask = %s:" %(
          inside), file = self.log)
      self._print_overall_u(aniso_b_cart,b_iso)

      if (inside or len(working_scale_factor_info.value_list)==1) and \
         replace_inside and \
             direction_vectors and direction_vectors != [None]:
        replace_u_cart_to_remove = True
        print("\nReplacing values inside mask with TLS-derived scale factors",
            file = self.log)
      else:
        replace_u_cart_to_remove = False
      self._analyze_scale_factor_info(
        working_scale_factor_info,
        tlso_info_by_region[inside],
        coordinate_shift_to_apply_before_tlso =
           coordinate_shift_to_apply_before_tlso,
        require_positive_definite = require_positive_definite,
        replace_u_cart_to_remove = replace_u_cart_to_remove,
        get_tls_from_u = get_tls_from_u,
        aniso_b_cart = aniso_b_cart,
        b_iso = b_iso,
        log = self.log)


    print("\nOverall average anisotropy by region:",file = self.log)
    for inside in inside_list:
      where = inside_dict[inside]
      tlso_info_by_region[inside].mean_u_cart /= max(
           1,tlso_info_by_region[inside].mean_u_cart_n)
      mean_u_cart = tuple(tlso_info_by_region[inside].mean_u_cart)
      print("%6s   (n = %4s)   (%6.2f,%6.2f,%6.2f,%6.2f,%6.2f,%6.2f) " %(
        tuple([where]+[tlso_info_by_region[inside].mean_u_cart_n]+
           list(mean_u_cart))),
        file = self.log)

    return tlso_info_by_region[True] # inside

  def _get_scale_factor_info_most_anisotropy(self,scale_factor_info,
     use_lowest = False):
    '''  Find scale_factor info with best resolution '''

    best_scaling_group_info = None
    best_xyz = None
    best_average  = None
    for scaling_group_info, xyz in zip(
       scale_factor_info.value_list,
       scale_factor_info.xyz_list):
      if scaling_group_info.get('cc_b_cart_as_u_cart'):
        average_aniso = flex.abs(flex.double(
         scaling_group_info.cc_b_cart_as_u_cart)).min_max_mean().mean

        if best_average is None or (
           ((not use_lowest) and
            average_aniso > best_average) or
           ((use_lowest) and
            average_aniso < best_average)) :
          best_average = average_aniso
          best_xyz = xyz
          best_scaling_group_info =scaling_group_info
    average_scale_factor_info = group_args(
       value_list = [best_scaling_group_info],
       xyz_list = [best_xyz])
    return average_scale_factor_info


  def _get_scale_factor_info_best_resolution(self,scale_factor_info,
     use_lowest = False):
    '''  Find scale_factor info with best resolution '''

    best_scaling_group_info = None
    best_xyz = None
    best_average  = None
    for scaling_group_info, xyz in zip(
       scale_factor_info.value_list,
       scale_factor_info.xyz_list):
      average_cc_star_list = self._get_average_cc_star_list(scaling_group_info)
      if best_average is None or (
         ((not use_lowest) and
          average_cc_star_list.min_max_mean().mean > best_average) or
         ((use_lowest) and
          average_cc_star_list.min_max_mean().mean < best_average)):
        best_average = average_cc_star_list.min_max_mean().mean
        best_xyz = xyz
        best_scaling_group_info =scaling_group_info
    average_scale_factor_info = group_args(
       value_list = [best_scaling_group_info],
       xyz_list = [best_xyz])
    return average_scale_factor_info

  def _get_average_cc_star_list(self, scaling_group_info):
    """Get average cc_star from scaling_info"""
    average_cc_star_list = None
    for si in scaling_group_info.scaling_info_list:
      if average_cc_star_list is None:
        average_cc_star_list = si.cc_list.deep_copy()
      else:
        average_cc_star_list += si.cc_list
    if average_cc_star_list:
      average_cc_star_list /= len(scaling_group_info.scaling_info_list)
    return average_cc_star_list

  def _average_scale_factor_info_over_xyz(self, scale_factor_info):
    '''
    Average scale_factor_info over xyz
    '''
    scaling_group_info_list = scale_factor_info.value_list

    # Create an average scale_factor_info object
    average_scaling_group_info = group_args(
      group_args_type = 'averaged (over xyz) scaling_info_object',
      direction_vectors = None,
      scaling_info_list = None,
      overall_si = None,
      overall_scale = None,
      aa_b_cart_as_u_cart = None,
      bb_b_cart_as_u_cart = None,
      ss_b_cart_as_u_cart = None,
      uu_b_cart_as_u_cart = None,
     )

    average_scale_factor_info = group_args(
       group_args_type = 'averaged scale_factor_info',
       value_list = [average_scaling_group_info],
       xyz_list = [None],
     )


    average_overall_scale = None

    n_used = 0
    for scaling_group_info in scaling_group_info_list:
      # Catch case with missing values and skip it
      ok = True
      for si in scaling_group_info.scaling_info_list:
        for key in ('target_scale_factors','cc_list','rms_fo_list'):
          if getattr(si,key) is None:
            ok = False
      if not ok:
        continue
      n_used += 1 # ok here

      if scaling_group_info.get('overall_scale'):
        if not average_overall_scale:
          average_overall_scale = flex.double(
            scaling_group_info.overall_scale.size(),0)
        average_overall_scale += scaling_group_info.overall_scale

      if not average_scaling_group_info.direction_vectors:
        average_scaling_group_info.direction_vectors = \
           scaling_group_info.direction_vectors

      for key in ('aa_b_cart_as_u_cart',
        'fo_b_cart_as_u_cart','uu_b_cart_as_u_cart',
        'bb_b_cart_as_u_cart', 'ss_b_cart_as_u_cart',
        ):

        if not average_scaling_group_info.get(key):
          average_scaling_group_info.add(key=key,value=flex.double(
            scaling_group_info.get(key)))

        elif scaling_group_info.get(key):
          xx = average_scaling_group_info.get(key)
          xx += flex.double(scaling_group_info.get(key))

      if not average_scaling_group_info.scaling_info_list:
        average_scaling_group_info.scaling_info_list = []
        for x in scaling_group_info.scaling_info_list:
          average_scaling_group_info.scaling_info_list.append(deepcopy(x))
      else:
        for si, average_si in zip(
            scaling_group_info.scaling_info_list,
            average_scaling_group_info.scaling_info_list):
          for key in ('target_scale_factors','cc_list','rms_fo_list'):
            setattr(average_si,key,
              getattr(average_si,key) + getattr(si,key))

      if not average_scaling_group_info.overall_si:
        average_scaling_group_info.overall_si = \
           deepcopy(scaling_group_info.overall_si)
      else:
        for key in ('target_scale_factors','cc_list','rms_fo_list'):
          setattr(average_scaling_group_info.overall_si,key,
            getattr(average_scaling_group_info.overall_si,key) +
             getattr(scaling_group_info.overall_si,key))


    if scaling_group_info_list and average_scaling_group_info.scaling_info_list:
      for key in ('target_scale_factors','cc_list','rms_fo_list'):
        for si in average_scaling_group_info.scaling_info_list:
          setattr(si,key, getattr(si,key)/max(1,n_used))
        setattr(average_scaling_group_info.overall_si,key,
           getattr(average_scaling_group_info.overall_si,key)/
           len(scaling_group_info_list))
      if average_scaling_group_info.overall_scale:
        average_scaling_group_info.overall_scale = \
          average_overall_scale/max(1, n_used)


      for key in ('aa_b_cart_as_u_cart','fo_b_cart_as_u_cart',
        'uu_b_cart_as_u_cart',
        'bb_b_cart_as_u_cart', 'ss_b_cart_as_u_cart',
        ):
        xx = average_scaling_group_info.get(key)
        if xx:
          xx /= max(1, n_used)

    avg = group_args(
      group_args_type = 'average scale_factor_info (averaged over xyz)',
      value_list = [average_scaling_group_info],
      xyz_list = [None],
     )
    return avg

  def _display_scale_values(self,
      si_list = None,
      overall_values = None,
      direction_vectors = None,
      key = 'target_scale_factors',
      text = 'Scale factors',
      extra_text = '',
      overall_text = ' ALL ',
      decimal_places = 2,
      ):
    """Display scale values"""
    assert len(list(direction_vectors))==len(si_list)

    n = len(si_list)
    sthol2_list = si_list[0].target_sthol2
    n_bins = len(sthol2_list)

    if not overall_values:
      nn = 0
      overall_values = flex.double(n_bins,0.)
      for si in si_list:
        values = si.get(key)
        if values is not None:
          overall_values += values
          nn += 1
      overall_values /= max(1,nn)

    print("\n D-min                %s by direction vector %s: " %(
       text,extra_text)+
         "\n        %s " %overall_text, file = self.log, end = "")
    for k in range(n):
      print("  %4s " %(k+1), end = "",file = self.log)
    print("\n           ", file = self.log)

    for i in range(n_bins):
      dd = 0.5/sthol2_list[i]**0.5
      if decimal_places == 1:
        print ("%6.1f  %7.1f " %(dd,overall_values[i]),
         file = self.log, end = "")
      else:
        print ("%6.2f  %7.2f " %(dd,overall_values[i]),
         file = self.log, end = "")
      for k in range(n):
        if decimal_places == 1:
          print (" %5.1f " %(si_list[k].get(key)[i] if (
            si_list[k] and si_list[k].get(key)) else 0),
            file = self.log, end= "")
        else:
          print (" %5.2f " %(si_list[k].get(key)[i] if (
            si_list[k] and si_list[k].get(key))else 0),
            file = self.log, end= "")
      print("", file = self.log)


  def _summarize_scale_factor_info(self, scale_factor_info,
     aniso_b_cart = None,
     b_iso = None,
     entry_number = 0):
    '''  Summarize scaling information'''


    if len(scale_factor_info.value_list) > 1:
      # Highest resolution
      average_scale_factor_info = self._get_scale_factor_info_best_resolution(
          scale_factor_info)
      print("\nScaling information for location with highest resolution" +
       " (Mean CC*: %.2f )" %(
           self._get_average_cc_star_list(
             average_scale_factor_info.value_list[0]).min_max_mean().mean),
          file = self.log)
      self._summarize_scale_factor_info(average_scale_factor_info,
        aniso_b_cart = aniso_b_cart,
        b_iso = b_iso
       )

      # Lowest resolution
      average_scale_factor_info = self._get_scale_factor_info_best_resolution(
          scale_factor_info, use_lowest = True)
      print("\nScaling information for location with lowest resolution" +
       " (Mean CC*: %.2f )" %(
           self._get_average_cc_star_list(
             average_scale_factor_info.value_list[0]).min_max_mean().mean),
          file = self.log)
      self._summarize_scale_factor_info(average_scale_factor_info,
        aniso_b_cart = aniso_b_cart,
        b_iso = b_iso
       )

      # Worst anisotropy
      average_scale_factor_info = self._get_scale_factor_info_most_anisotropy(
          scale_factor_info, use_lowest = True)
      if average_scale_factor_info.value_list and \
          average_scale_factor_info.value_list[0] and \
        average_scale_factor_info.value_list[0].fo_b_cart_as_u_cart:
        cc_abs = flex.abs(flex.double(
        average_scale_factor_info.value_list[0].fo_b_cart_as_u_cart)
           ).min_max_mean().mean
        print("\nScaling information for location with most anisotropy" +
         " (Mean abs(anisotropy)): %.2f )" %(cc_abs),
           file = self.log)
        self._summarize_scale_factor_info(average_scale_factor_info,
          aniso_b_cart = aniso_b_cart,
          b_iso = b_iso
         )

      # Anisotropy values vs position
      print("\nAnisotropy vs position", file = self.log)
      n_use = len(scale_factor_info.value_list)
      if n_use > 50 and not self.verbose:
        n_use = 50
        print ("First 50 listed...use verbose=True for remainder",
            file = self.log)

      for text, kw in (
          ['Estimated anisotropic fall-off of the data relative to ideal',
             'uu_b_cart_as_u_cart'],
          ['Anisotropy of the data', 'aa_b_cart_as_u_cart'],
          ['Anisotropy of the uncertainties', 'bb_b_cart_as_u_cart'],
          ['Anisotropy of the scale_factors', 'ss_b_cart_as_u_cart'],
           ):

        self._print_overall_u(aniso_b_cart,b_iso)
        self._print_aniso_by_xyz(text, kw, scale_factor_info, n_use = n_use)

      #  Overall
      print("\nSummary of overall average scaling information ",
          file = self.log)
      average_scale_factor_info = self._average_scale_factor_info_over_xyz(
          scale_factor_info)
      self._summarize_scale_factor_info(average_scale_factor_info,
          aniso_b_cart = aniso_b_cart,
          b_iso = b_iso
         )

      return


    # Summary for one location or one average

    scaling_group_info  = scale_factor_info.value_list[entry_number]

    xyz = scale_factor_info.xyz_list[entry_number]
    if xyz:
      print("\nXYZ: (%.3f, %.3f, %.3f)" %(tuple(xyz)), file = self.log)

    self._print_overall_u(aniso_b_cart,b_iso)

    self._display_scale_values(
      si_list = scaling_group_info.scaling_info_list,
      direction_vectors = scaling_group_info.direction_vectors,
      overall_values = scaling_group_info.scaling_info_list[0].rms_fc_list,
      key = 'rms_fo_list',
      text = 'RMS Fobs ',
      overall_text = 'RMS Fc',
      decimal_places = 1,
     )

    self._display_scale_values(
      si_list = scaling_group_info.scaling_info_list,
      direction_vectors = scaling_group_info.direction_vectors,
      key = 'cc_list',
      text = 'Estimated CC*',
      extra_text = " (Mean CC*: %.2f)" %(
       self._get_average_cc_star_list(scaling_group_info).min_max_mean().mean),
      decimal_places = 2,
     )

    self._display_scale_values(
      si_list = scaling_group_info.scaling_info_list,
      direction_vectors = scaling_group_info.direction_vectors,
      overall_values = scaling_group_info.overall_si.target_scale_factors,
      key = 'target_scale_factors',
      text = 'Scale factors',
      extra_text = '',
      overall_text = ' ALL ',
      decimal_places = 2,
     )

    aa_b_cart_as_u_cart = scaling_group_info.get('aa_b_cart_as_u_cart')
    bb_b_cart_as_u_cart = scaling_group_info.get('bb_b_cart_as_u_cart')
    ss_b_cart_as_u_cart = scaling_group_info.get('ss_b_cart_as_u_cart')
    uu_b_cart_as_u_cart = scaling_group_info.get('uu_b_cart_as_u_cart')
    fo_b_cart_as_u_cart = scaling_group_info.get('fo_b_cart_as_u_cart')

    if uu_b_cart_as_u_cart:
      self._print_overall_u(aniso_b_cart,b_iso)
      print("\n Estimated anisotropic fall-off of the data relative to ideal\n"+
        "(Positive means amplitudes fall off more in this direction)\n " +
       "(  X,      Y,      Z,    XY,    XZ,    YZ)\n"+
       "(%.3f, %.3f, %.3f, %.3f, %.3f, %.3f) " %(
       tuple(uu_b_cart_as_u_cart)), file = self.log)

    if aa_b_cart_as_u_cart:
      print("\n Anisotropy of the data\n"+
        "(Positive means amplitudes fall off more in this direction)\n " +
       "(  X,      Y,      Z,    XY,    XZ,    YZ)\n"+
       "(%.3f, %.3f, %.3f, %.3f, %.3f, %.3f) " %(
       tuple(aa_b_cart_as_u_cart)), file = self.log)

    if bb_b_cart_as_u_cart:
      print("\n Anisotropy of the uncertainties\n"+
        "(Positive means uncertainties decrease more in this direction)\n " +
       "(  X,      Y,      Z,    XY,    XZ,    YZ)\n"+
       "(%.3f, %.3f, %.3f, %.3f, %.3f, %.3f) " %(
       tuple(bb_b_cart_as_u_cart)), file = self.log)

    if ss_b_cart_as_u_cart:
      print("\n Anisotropy of the scale factors\n"+
       "(Positive means scale factors increase more in this direction)\n " +
       "(  X,      Y,      Z,    XY,    XZ,    YZ)\n"+
       "(%.3f, %.3f, %.3f, %.3f, %.3f, %.3f) " %(
       tuple(ss_b_cart_as_u_cart)), file = self.log)

  def _print_aniso_by_xyz(self, text, kw, scale_factor_info, n_use):
      """Summary anisotropy by region"""

      print("\n    %s by position" %(text), file = self.log)
      print("\n      Box center                       U values","\n",
           "   X      Y      Z          ",
           "  X      Y      Z     XY     XZ     YZ\n", file = self.log)

      for xyz, scaling_group_info in zip(
          scale_factor_info.xyz_list,
          scale_factor_info.value_list[:n_use]):
        u_cart = scaling_group_info.get(kw)
        if not u_cart: continue
        print( "%7.1f %7.1f %7.1f " %(xyz),
          "   %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f  " %(
            tuple(u_cart)), file = self.log)

  def _analyze_scale_factor_info(self,
        scale_factor_info,
        tlso_info_in_region,
        coordinate_shift_to_apply_before_tlso = None,
        require_positive_definite = None,
        replace_u_cart_to_remove= None,
        get_tls_from_u = None,
        aniso_b_cart = None,
        b_iso = None,
        log = sys.stdout):

      """Summarize anisotropy and scale factors"""

      self._summarize_scale_factor_info(scale_factor_info,
       aniso_b_cart = aniso_b_cart,
       b_iso = b_iso,
       )

      # Interpret as TLS  and optionally replace values

      self._summarize_and_optionally_replace_aniso_with_tls(
        scale_factor_info,
        tlso_info_in_region,
        coordinate_shift_to_apply_before_tlso =
           coordinate_shift_to_apply_before_tlso,
        require_positive_definite = require_positive_definite,
        replace_u_cart_to_remove = replace_u_cart_to_remove,
        get_tls_from_u = get_tls_from_u,
        aniso_b_cart = aniso_b_cart,
        b_iso = b_iso,
        log = self.log)


  def _summarize_and_optionally_replace_aniso_with_tls(self,
        scale_factor_info,
        tlso_info_in_region,
        coordinate_shift_to_apply_before_tlso = None,
        require_positive_definite = None,
        replace_u_cart_to_remove = None,
        get_tls_from_u = None,
        aniso_b_cart = None,
        b_iso = None,
        log = sys.stdout):

      """Summarize anisotropy and optionally replace anisotropy information
       with values calculated from TLS"""
      xyz_list = scale_factor_info.xyz_list
      scaling_group_info_list = scale_factor_info.value_list

      uanisos = flex.sym_mat3_double()
      xyz_list_use = flex.vec3_double()

      print("\nSummary of scaling as TLS",file = self.log)
      self._print_overall_u(aniso_b_cart,b_iso)

      if get_tls_from_u:
        print("\nTLS is overall fall-off with resolution",
           file = self.log)
      else:
       print("\nTLS is desired anisotropy of scale factors",
           file = self.log)

      for xyz,scaling_group_info in zip(xyz_list,scaling_group_info_list):
        if get_tls_from_u:
          u_cart = scaling_group_info.get('uu_b_cart_as_u_cart')
        else: # usual
          u_cart = scaling_group_info.get('ss_b_cart_as_u_cart')
          if u_cart:
            u_cart = tuple(-flex.double(u_cart)) # MINUS
        if not u_cart:  # missing
          continue
        # we want anisotropy, not the correction
        u_cart = self._add_overall_to_u_cart(u_cart,aniso_b_cart,b_iso)
        uanisos.append(u_cart)
        xyz_list_use.append(xyz)
        tlso_info_in_region.mean_u_cart+=flex.double(u_cart)
        tlso_info_in_region.mean_u_cart_n+=1
      if xyz_list_use.size() < 1:
        return # nothing to do
      if coordinate_shift_to_apply_before_tlso:
        xyz_list_use += coordinate_shift_to_apply_before_tlso

      from mmtbx.tls import tools
      cm = xyz_list_use.mean()
      result = tools.tls_from_uaniso_minimizer(
        uaniso         = uanisos,
        T_initial      = [0,0,0,0,0,0],
        L_initial      = [0,0,0,0,0,0],
        S_initial      = [0,0,0,0,0,0,0,0,0],
        refine_T       = True,
        refine_L       = True,
        refine_S       = True,
        origin         = cm,
        sites          = xyz_list_use,
        max_iterations = 100)
      if require_positive_definite:
        # Make sure they are positive definite
        T = adptbx.eigenvalue_filtering(result.T_min)
        L = adptbx.eigenvalue_filtering(result.L_min)
      else:
        T = result.T_min
        L = result.L_min
      S=result.S_min

      # Decide what aniso to apply. Usually the one we just calculated

      if get_tls_from_u:
        print("\nTLS analysis of anisotropic fall-off of data", file = self.log)
      else:
        print("\nTLS analysis of anisotropy of scale factors", file = self.log)

      tlso_value = tlso(t = T, l = L, s = S, origin = cm)

      tlso_info_in_region.tlso = tlso_value

      if coordinate_shift_to_apply_before_tlso:
        shift = col(coordinate_shift_to_apply_before_tlso)
      else:
        shift = col((0,0,0))

      new_anisos_as_anisotropy= uaniso_from_tls_one_group(tlso = tlso_value,
          sites_cart = xyz_list + shift, # xyz_list, not xyz_list_use
          zeroize_trace=False)

      # we want correction not aniso:
      new_anisos = []
      for u in new_anisos_as_anisotropy:
        new_anisos.append(tuple(flex.double(u)))

      for new_u_cart, xyz, scaling_group_info in zip(
          new_anisos,
          xyz_list,
          scaling_group_info_list):

        if replace_u_cart_to_remove and (not get_tls_from_u):
          scaling_group_info.ss_b_cart_as_u_cart = tuple(
             -flex.double(new_u_cart)) # MINUS
        elif replace_u_cart_to_remove and get_tls_from_u:
          scaling_group_info.uu_b_cart_as_u_cart = tuple(
             flex.double(new_u_cart))
          # ss_b_cart_as_u_cart has trace removed:
          tr = flex.double(new_u_cart[:3]).min_max_mean().mean
          tr_as_u_cart = (tr,tr,tr,0,0,0)
          scaling_group_info.ss_b_cart_as_u_cart = tuple(
            [- (u - t)  for u,t in zip(new_u_cart,tr_as_u_cart)]) # MINUS

      if get_tls_from_u:
        print("\nMean overall fall-off with resolution as TLS:",file = self.log)
      else:
        print("\nMean anisotropy as TLS:",file = self.log)
      self._print_overall_u(aniso_b_cart,b_iso)

      print("T: (%.2f, %.2f, %.2f, %.2f, %.2f, %.2f)" %(
        tuple(T)),file = self.log)
      print("L: (%.2f, %.2f, %.2f, %.2f, %.2f, %.2f)" %(
         tuple(L)), file = self.log)
      print("S: (%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f)" %(
         tuple(S)), file = self.log)

  def _remove_overall_from_u_cart(self,u_cart,aniso_b_cart,b_iso):
    """Remove overall B-iso from U-cart"""
    if not aniso_b_cart:
      return u_cart
    else:
      overall_u_value = adptbx.b_as_u(tuple(
       flex.double(aniso_b_cart)-flex.double((
          b_iso,b_iso,b_iso,0,0,0))))
      return tuple(flex.double(u_cart) - flex.double(overall_u_value))

  def _add_overall_to_u_cart(self,u_cart,aniso_b_cart,b_iso):
    """Add overall B-iso to U-cart"""
    if not aniso_b_cart:
      return u_cart
    else:
      overall_u_value = adptbx.b_as_u(tuple(
        flex.double(aniso_b_cart)-flex.double((
          b_iso,b_iso,b_iso,0,0,0))))
      return tuple(flex.double(u_cart) + flex.double(overall_u_value))

  def _print_overall_u(self,aniso_b_cart,b_iso):
    """Print overall U values"""
    if aniso_b_cart:
      u_value = adptbx.b_as_u(tuple(flex.double(aniso_b_cart)-flex.double((
          b_iso,b_iso,b_iso,0,0,0))))
      print("\nNOTE: All values apply after removal of overall U of:\n ",
          " (%.2f, %.2f, %.2f, %.2f, %.2f, %.2f)\n" %(tuple(u_value)),
          file = self.log)


  def _run_group_of_anisotropic_sharpen(self,
      map_id  = 'map_manager',
      map_id_1 = 'map_manager_1',
      map_id_2 = 'map_manager_2',
      map_id_to_be_scaled_list = None,
      map_id_scaled_list = None,
      mask_id = None,
      exclude_points_outside_density = None,
      minimum_boxes_inside_density = None,
      resolution = None,
      d_min = None,
      k_sol = None,
      b_sol = None,
      n_bins = None,
      n_boxes = None,
      core_box_size = None,
      box_cushion = None,
      smoothing_radius = None,
      rmsd = None,
      nproc = None,
      optimize_b_eff = None,
      equalize_power = None,
      is_model_based = False,
      is_external_based = False,
      get_scale_as_aniso_u = None,
      use_dv_weighting = None,
      n_direction_vectors = None,
      run_analyze_anisotropy = None,
      spectral_scaling = None,
      expected_rms_fc_list = None,
      expected_ssqr_list = None,
      expected_ssqr_list_rms = None,
      tlso_group_info = None,
      get_tls_from_u = None,
      model_id_for_rms_fc = None,
      replace_aniso_with_tls_equiv = None,
      minimum_low_res_cc = None,
      max_abs_b = None,
      get_tls_info_only = None,
      coordinate_shift_to_apply_before_tlso = None,
      sharpen_all_maps = None,
      temp_dir = 'TEMP_ANISO_LOCAL',
      aniso_b_cart = None,
      b_iso = None,
     ):
    '''
    Run local sharpening in groups with focus on reflections along
    direction vectors. Then combine results
    Summary of method:

    A map of one scale factor is the scale factor to apply in real space
       at each xyz for any contribution from an xyz in that bin.
    (1) we calculate position-dependent target_scale_factors (n_bins)
      for each direction vector (typically n=6 or 12).  Total of about
      240 bins/directions.
    (2) each resolution bin has a set of weights for all reflections w_hkl.
      These are just binner.apply_scale of (0 all other bins and 1 this bin)
    (3) each direction has a set of weights w_dv_hkl. These are just the
      dot product of the direction and the normalized (hkl). On the fly.
    (4) To sum up:
       one bin (sel), one direction vector dv, weights w_dv,
         weights_resolution_bin
       a.calculate value_map map with map_coeffs * w_dv * w_resolution_bin
       b. calculate weight map from position-dependent target_scale_factors
          for dv
       c multiply weight_map * value_map and sum over all bins, dv

    (5) To parallelize: run a group of sums, write out maps, read in and sum up.
        '''

    # Get the kw we have
    from libtbx import adopt_init_args
    kw_obj = group_args()
    adopt_init_args(kw_obj, locals())
    kw = kw_obj() # save calling parameters in kw as dict
    del kw['adopt_init_args'] # REQUIRED
    del kw['kw_obj'] # REQUIRED
    del kw['temp_dir'] # REQUIRED

    assert n_bins is not None

    print ("\nRunning anisotropic local sharpening with nproc = %s " %(
       nproc), file = self.log)

    setup_info = self._get_box_setup_info(map_id_1, map_id_2,
      resolution,
      d_min,
      smoothing_radius = smoothing_radius,
      n_boxes = n_boxes,
      core_box_size = core_box_size,
       )

    resolution = setup_info.resolution
    self.set_resolution(resolution)
    print ("Nominal resolution of map: %.2f A " %(resolution),
      file = self.log)

    if spectral_scaling and (not expected_rms_fc_list):
        from cctbx.development.approx_amplitude_vs_resolution import \
          approx_amplitude_vs_resolution
        aavr = approx_amplitude_vs_resolution(
           d_min = setup_info.minimum_resolution,
           n_bins = n_bins,
           k_sol = k_sol,
           b_sol = b_sol,
           map_model_manager = self,
           model = self.get_model_by_id(model_id=model_id_for_rms_fc))
        kw['expected_rms_fc_list'] = aavr.get_target_scale_factors()
    if expected_ssqr_list_rms and not expected_ssqr_list:
        from cctbx.development.approx_amplitude_vs_resolution import \
          get_expected_ssqr_list
        kw['expected_ssqr_list'] = get_expected_ssqr_list(
           d_min = setup_info.minimum_resolution,
           n_bins = n_bins,
           expected_ssqr_list_rms = expected_ssqr_list_rms,
           map_model_manager = self,
           out = self.log)
    # Get list of direction vectors (based on anisotropy of map)
    direction_vectors = self._get_aniso_direction_vectors(map_id,
      n_direction_vectors = n_direction_vectors)

    # Run local_fsc for each direction vector
    print("\nEstimating scale factors for %s direction_vectors" %(
        len(direction_vectors)), file = self.log)
    print("Number of resolution bins: %s  Number of processors: %s" %(
          n_bins,nproc), file = self.log)

    # Get scale factors vs resolution and location
    scale_factor_info = self.local_fsc(
        return_scale_factors = True,
        direction_vectors=direction_vectors,
         **kw)

    xyz_list = scale_factor_info.xyz_list
    if exclude_points_outside_density and (
       not scale_factor_info.exclude_points_outside_density):
      print("Points inside density are filled in so all points are now inside",
          file = self.log)
      everything_is_inside = True
    else:
      everything_is_inside = None
    replace_inside = (replace_aniso_with_tls_equiv and get_scale_as_aniso_u)
    # Summarize U vs xyz and vs inside/outside
    tls_info = self._analyze_aniso(scale_factor_info,
      tlso_group_info = tlso_group_info,
      get_tls_from_u = get_tls_from_u,
      map_id=map_id,
      mask_id=mask_id,  # can supply mask_id
      replace_inside = replace_inside,
      coordinate_shift_to_apply_before_tlso =
           coordinate_shift_to_apply_before_tlso,
      everything_is_inside = (everything_is_inside or len(xyz_list) ==1),
      aniso_b_cart = aniso_b_cart,
      b_iso = b_iso,
     )
    if get_tls_info_only:
      return tls_info

    if replace_inside:
      self._update_scale_factor_info_from_aniso(scale_factor_info,
        max_abs_b = max_abs_b,
        get_tls_from_u = get_tls_from_u)

    setup_info.kw = kw

    print("Applying interpolated scale factors (resolution and direction)",
     file = self.log)

    # Apply interpolated scale_factors (vs resolution and direction). Split
    # into groups by direction

    #  Run scaling application for each map
    for id,new_id in zip(
         kw['map_id_to_be_scaled_list'],kw['map_id_scaled_list']):
      setup_info.kw['map_id_to_be_scaled'] = id

      temp_dir = self._create_temp_dir(temp_dir)  # for returning big files
      setup_info.temp_dir = temp_dir

      # Set up to run for each direction
      index_list=[]
      for i in range(len(direction_vectors)):
        index_list.append({'i':i})

      from libtbx.easy_mp import run_parallel
      results = run_parallel(
        method = self._multiprocessing,
        qsub_command = self._queue_run_command,
        nproc = nproc,
        target_function = run_anisotropic_scaling_as_class(
           map_model_manager = self,
           direction_vectors = direction_vectors,
           scale_factor_info = scale_factor_info,
           setup_info = setup_info),
        preserve_order=False,
        kw_list = index_list)

      # Results is list of map names.  Read them in, sum up, and we're done
      from iotbx.data_manager import DataManager
      dm = DataManager()
      map_data = None
      for result in results:
        if result and result.file_name:
          mm = dm.get_real_map(result.file_name)
          mm.shift_origin()
          if map_data is None:
            map_data = mm.map_data()
          else:
            map_data += mm.map_data()
      self._remove_temp_dir(temp_dir)
      new_map_manager = self.get_any_map_manager().customized_copy(
        map_data = map_data,
        )
      new_map_manager.file_name = self.get_map_manager_by_id(id).file_name
      print("Setting map "+
       "'%s' in map_manager '%s' to local aniso-scaled map '%s'" %(
          id,self.name,new_id), file = self.log)

      self.add_map_manager_by_id(map_id = new_id,
          map_manager = new_map_manager)
    print("Done applying interpolated scale factors",
     file = self.log)

    return tls_info

  def _local_sharpen(self,
      map_id  = 'map_manager',
      map_id_1 = 'map_manager_1',
      map_id_2 = 'map_manager_2',
      map_id_to_be_scaled_list = None,
      map_id_scaled_list = None,
      mask_id = None,
      exclude_points_outside_density = None,
      minimum_boxes_inside_density = None,
      resolution = None,
      d_min = None,
      k_sol = None,
      b_sol = None,
      n_bins = None,
      n_boxes = None,
      core_box_size = None,
      box_cushion = None,
      smoothing_radius = None,
      rmsd = None,
      nproc = None,
      optimize_b_eff = None,
      equalize_power = None,
      is_model_based = False,
      is_external_based = False,
      get_scale_as_aniso_u = None,
      use_dv_weighting = None,
      n_direction_vectors = None,
      run_analyze_anisotropy = None,
      spectral_scaling = None,
      expected_rms_fc_list = None,
      expected_ssqr_list = None,
      expected_ssqr_list_rms = None,
      tlso_group_info = None,
      get_tls_from_u = None,
      model_id_for_rms_fc = None,
      replace_aniso_with_tls_equiv = None,
      anisotropic_sharpen = None,
      minimum_low_res_cc = None,
      max_abs_b = None,
      get_tls_info_only = None,
      coordinate_shift_to_apply_before_tlso = None,
      sharpen_all_maps = None,
      aniso_b_cart = None,
      b_iso = None,
     ):

    '''
     Scale map_id_to_be_scaled with local scale factors identified from map_id_1 and map_id_2
     Changes the working map_manager unless map_id_scaled is set

    '''

    # Get the kw we have
    from libtbx import adopt_init_args
    kw_obj = group_args()
    adopt_init_args(kw_obj, locals())
    kw = kw_obj() # save calling parameters in kw as dict
    del kw['adopt_init_args'] # REQUIRED
    del kw['kw_obj']  # REQUIRED
    del kw['anisotropic_sharpen']  # REQUIRED

    # Checks
    assert self.get_map_manager_by_id(map_id)
    for id in map_id_to_be_scaled_list:
      assert self.get_map_manager_by_id(id)
    assert (
    (self.get_map_manager_by_id(map_id_1) or
        is_model_based or is_external_based) and
       self.get_map_manager_by_id(map_id_2))

    assert n_bins is not None

    if nproc is None and self._nproc:
      nproc = self._nproc
    elif nproc is None:
      nproc = 1
    kw['nproc'] = nproc

    if anisotropic_sharpen:  # run N times with different direction vectors
      tls_info = self._run_group_of_anisotropic_sharpen(**kw)
      return tls_info  # tlso for region inside mask

    del kw['aniso_b_cart']  # REQUIRED
    del kw['b_iso']  # REQUIRED


    # Get scale factors vs resolution and location
    scale_factor_info = self.local_fsc(
      direction_vectors = [None],
      return_scale_factors = True, **kw)

    # scale_factor_info.value_list is a set of scaling_group_info objects.
    # scale_factor_info.xyz_list are the coordinates where these apply
    # value_list is a set of scaling_group_info objects, one per xyz.
    """
    scaling_group_info group_args object:
      direction_vectors: direction vectors dv for anisotropy calculations
      scaling_info_list: si (scaling_info) objects, one for each dv
        each si:  si.target_scale_factors   # scale factors vs sthol2
        si.target_sthol2 # sthol2 values  d = 0.25/sthol2**0.5
                  si.d_min_list
                  si.cc_list
                  si.low_res_cc # low-res average
      ss_b_cart_as_u_cart: anisotropic part of overall correction factor
      overall_scale: radial part of overall correction factor
      NOTE: total scale = total aniso scale * overall_scale[0]
    """

    xyz_list = scale_factor_info.xyz_list
    d_min = scale_factor_info.d_min
    smoothing_radius = scale_factor_info.setup_info.smoothing_radius
    assert n_bins == scale_factor_info.n_bins # must match

    self._summarize_scale_factor_info(scale_factor_info,
      aniso_b_cart = aniso_b_cart,
      b_iso = b_iso,
     )
    average_scale_factors = get_average_scale_factors(scale_factor_info)

    # Get Fourier coefficients for maps based on map_id_to_be_scaled
    for id, new_id in zip(map_id_to_be_scaled_list, map_id_scaled_list):

      map_coeffs = self.get_map_manager_by_id(id
           ).map_as_fourier_coefficients(d_min = d_min)
      f_array_info = get_map_coeffs_as_fp_phi(map_coeffs, n_bins = n_bins,
         d_min = d_min)

      if aniso_b_cart: # first apply aniso_b_cart to the whole array
        apply_aniso_b_cart_to_f_array_info(f_array_info,
         b_iso, d_min, aniso_b_cart)


      new_map_data = flex.double(flex.grid(
          self.get_map_manager_by_id(id).map_data().all()), 0.)
      # Get map for each shell of resolution

      # Get normalizations for each reflection so we can interpolate over bins
      # We are going to use each reflection in all bin calculations, but
      #   weighting them by how close they are to the center of that bin

      normalization_data = get_normalization_data_for_unit_binning(
        f_array_info.f_array)

      for i_bin in f_array_info.f_array.binner().range_used():
        # Get scale values for i_bin at all points xyz for dv 0

        scale_value_list,xyz_used_list = self._get_scale_values_for_bin(
          xyz_list=xyz_list,
          i_bin = i_bin,
          scale_factor_info = scale_factor_info,)

        # Get a map that has scale factor for this resolution vs xyz
        default_value = average_scale_factors[i_bin-1]
        weight_mm = self._create_full_size_map_manager_with_value_list(
          xyz_list = xyz_used_list,
          value_list = scale_value_list,
          smoothing_radius = smoothing_radius,
          default_value = default_value,
          n_boxes = n_boxes)

        # Multiply shell map data by weights
        weights_top_hat_shell = get_weights_for_unit_binning(
           f_array_info.f_array, i_bin) * normalization_data
        shell_map_coeffs = map_coeffs.customized_copy(
          data = map_coeffs.data() * weights_top_hat_shell)

        shell_mm = self.map_manager(
           ).fourier_coefficients_as_map_manager(shell_map_coeffs)

        new_map_data += weight_mm.map_data() * shell_mm.map_data()
      new_map_manager = self.get_any_map_manager().customized_copy(
        map_data = new_map_data)
      new_map_manager.file_name = self.get_map_manager_by_id(id).file_name
      print( "Adding locally scaled map data "+
        "as '%s' to map_model_manager '%s' as '%s' "%(
         id,self.name,new_id),file = self.log)
      self.add_map_manager_by_id(map_id = new_id,
          map_manager = new_map_manager)

  def _get_scale_factor_info_inside_mask(self,
     scale_factor_info,
     mask_map_manager = None,
     inside = True,
     everything_is_inside = None):
    """Calculate scale factors inside a mask"""
    new_xyz_list = flex.vec3_double()
    new_value_list = []
    for xyz, value in zip (scale_factor_info.xyz_list,
       scale_factor_info.value_list):
      if everything_is_inside:
          if xyz is not None:
            new_xyz_list.append(xyz)
          else:
            new_xyz_list.append((0,0,0))
          new_value_list.append(value)
      else: #check if inside
        site_frac=mask_map_manager.crystal_symmetry(
          ).unit_cell().fractionalize(xyz)
        if is_inside_mask(mask_map_manager,
            site_frac = site_frac,
            inside = inside):
          new_xyz_list.append(xyz)
          new_value_list.append(value)
    scale_factor_info_dict = scale_factor_info()
    new_scale_factor_info=group_args()
    for key in scale_factor_info_dict.keys():
      if not key in ['xyz_list', 'value_list']:
         setattr(new_scale_factor_info,key,
           deepcopy(scale_factor_info_dict[key]))
    new_scale_factor_info.value_list = new_value_list
    new_scale_factor_info.xyz_list = new_xyz_list

    return new_scale_factor_info

  def _get_scale_values_for_bin(self,
        xyz_list=None,
        i_bin = None,
        scale_factor_info = None,
        dv_id = 0):
    '''
    # Get scale values for i_bin at all points xyz for direction_vector dv_id
    Get the i_bin'th scale value for each point
    NOTE: i_bin starts at 1, not zero.
    '''
    scale_values = flex.double()
    xyz_used_list = flex.vec3_double()

    # scale_factor_info.value_list is a set of scaling_group_info objects.
    # scale_factor_info.xyz_list are the coordinates where these apply
    # scale_factor_info.n_bins is number of bins
    # value_list is a set of scaling_group_info objects, one per xyz.
    #  sgi (scaling_group_info):
    #   sgi.direction_vectors
    #   sgi.scaling_info_list: one si entry per direction
    #    si.target_scale_factors
    #    si.target_sthol2
    #    si.d_min_list
    #    si.cc_list
    #    si.low_res_cc # low-res average

    # scale_factor_info.value_list has one scaling_group_info object per xyz
    # value_list:  [ [scale_factor_info_1, scale_factor_info_2....12],[...]]

    for xyz,sgi in zip(xyz_list,scale_factor_info.value_list):
          # for one value of xyz
      # sgi.direction_vectors
      # sgi.scaling_info_list= [scaling_info_1, scaling_info_2....12]
      si = sgi.scaling_info_list[dv_id]
      if si and si.target_scale_factors:
        scale_values.append(si.target_scale_factors[i_bin-1])
        xyz_used_list.append(xyz)
      else:
        pass # failed
    return scale_values,xyz_used_list

  def local_fsc(self,
      map_id = 'map_manager',
      map_id_1 = 'map_manager_1',
      map_id_2 = 'map_manager_2',
      model_id = None,
      map_id_to_be_scaled_list = None, # NOTE: not used, just allows it in call
      map_id_scaled_list = None, # NOTE: not used, just allows it in call
      mask_id = None,
      exclude_points_outside_density = None,
      minimum_boxes_inside_density = None,
      resolution = None,
      d_min = None,
      k_sol = None,
      b_sol = None,
      max_resolution_ratio = None,
      min_bin_width = 20,
      n_bins = None,
      fsc_cutoff = 0.143,
      n_boxes = None,
      core_box_size = None,
      box_cushion = None,
      rmsd = None,
      smoothing_radius = None,
      nproc = None,
      is_model_based = None,
      optimize_b_eff = None,
      equalize_power = None,
      is_external_based = None,
      return_scale_factors = False,
      direction_vectors = None,
      minimum_low_res_cc = None,
      get_scale_as_aniso_u = None,
      use_dv_weighting = None,
      n_direction_vectors = None,
      run_analyze_anisotropy = None,
      spectral_scaling = None,
      expected_rms_fc_list = None,
      expected_ssqr_list = None,
      expected_ssqr_list_rms = None,
      tlso_group_info = None,
      get_tls_from_u = None,
      model_id_for_rms_fc = None,
      replace_aniso_with_tls_equiv = None,
      max_abs_b = None,
      get_tls_info_only = None,
      coordinate_shift_to_apply_before_tlso = None,
      sharpen_all_maps = None,
      n_bins_default = 2000,
      b_iso = None, # not used
      aniso_b_cart = None, # not used
      ):

    '''
      Calculates local Fourier Shell Correlations to estimate local resolution
      Creates map with smoothed local resolution

      Optionally estimates scale factors vs resolution at each point in map
      to apply to yield a locally-scaled map (return_scale_factors = True).

      If direction_vector is specified, weight scale factor calculation by
      dot product of reflection directions with direction_vector
    '''

    # Checks
    assert self.get_map_manager_by_id(map_id)
    if ( (is_model_based is None) and (is_external_based is None) and
       self.model() and (not
              (self.get_map_manager_by_id(map_id_1) and
              self.get_map_manager_by_id(map_id_2)) )):
         is_model_based = True # default to model-based if no info
         if model_id is not None:
           model = self.get_model_by_id(model_id)
         else:
           model = self.model()
         assert model is not None
         self.generate_map(map_id = 'model_map')
         map_id_1 = 'map_manager'
         map_id_2 = 'model_map'

    assert ( ( is_model_based or is_external_based) or
      (self.get_map_manager_by_id(map_id_1) and
       self.get_map_manager_by_id(map_id_2)))

    if n_bins is None:
      n_bins = n_bins_default
    if nproc is None and self._nproc:
      nproc = self._nproc
    elif nproc is None:
      nproc = 1

    # Get basic info including minimum_resolution (cutoff for map_coeffs)
    setup_info = self._get_box_setup_info(map_id_1, map_id_2,
      resolution,
      d_min,
      box_cushion,
      n_boxes,
      core_box_size,
      smoothing_radius = smoothing_radius)

    if spectral_scaling and (not expected_rms_fc_list):
        from cctbx.development.approx_amplitude_vs_resolution import \
          approx_amplitude_vs_resolution
        aavr = approx_amplitude_vs_resolution(
           d_min = setup_info.minimum_resolution,
           n_bins = n_bins,
           k_sol = k_sol,
           b_sol = b_sol,
           map_model_manager = self,
           model = self.get_model_by_id(model_id=model_id_for_rms_fc))
        expected_rms_fc_list = aavr.get_target_scale_factors()
    if expected_ssqr_list_rms and not expected_ssqr_list:
        from cctbx.development.approx_amplitude_vs_resolution import \
          get_expected_ssqr_list
        expected_ssqr_list = get_expected_ssqr_list(
           d_min = setup_info.minimum_resolution,
           n_bins = n_bins,
           expected_ssqr_list_rms = expected_ssqr_list_rms,
           map_model_manager = self,
           out = self.log)
    box_info = self.split_up_map_and_model_by_boxes(
      target_for_boxes = setup_info.n_boxes,
      box_cushion = setup_info.box_cushion,
      skip_empty_boxes = False,
      select_final_boxes_based_on_model = False, # required
      apply_box_info = False,
      mask_id = mask_id,
      exclude_points_outside_density = exclude_points_outside_density,
      minimum_boxes_inside_density = minimum_boxes_inside_density,
      )

    # Hold some things in box_info
    box_info.resolution = setup_info.resolution
    box_info.minimum_resolution = setup_info.minimum_resolution
    box_info.fsc_cutoff = fsc_cutoff
    box_info.n_bins = n_bins
    box_info.rmsd = rmsd
    box_info.return_scale_factors = return_scale_factors
    box_info.map_id = map_id
    box_info.map_id_1 = map_id_1
    box_info.map_id_2 = map_id_2
    box_info.is_model_based = is_model_based
    box_info.optimize_b_eff = optimize_b_eff
    box_info.equalize_power = equalize_power
    box_info.is_external_based = is_external_based
    box_info.direction_vectors = direction_vectors
    box_info.minimum_low_res_cc = minimum_low_res_cc
    box_info.get_scale_as_aniso_u = get_scale_as_aniso_u
    box_info.use_dv_weighting = use_dv_weighting
    box_info.n_direction_vectors = n_direction_vectors
    box_info.run_analyze_anisotropy = run_analyze_anisotropy
    box_info.expected_rms_fc_list = expected_rms_fc_list
    box_info.expected_ssqr_list = expected_ssqr_list
    box_info.expected_ssqr_list_rms = expected_ssqr_list_rms
    box_info.tlso_group_info = tlso_group_info
    box_info.model_id_for_rms_fc = model_id_for_rms_fc
    box_info.replace_aniso_with_tls_equiv = replace_aniso_with_tls_equiv

    log_hold = self.log
    if not self.verbose:
      self.set_log(null_out())
    results = self._run_fsc_in_boxes(
     nproc = nproc,
     box_info = box_info)
    self.set_log(log_hold)
    # results.value_list is a set of scaling_group_info objects.
    # results.xyz_list are the coordinates where these apply
    """
    scaling_group_info group_args object:
      direction_vectors: direction vectors dv for anisotropy calculations
      scaling_info_list: si (scaling_info) objects, one for each dv
        each si:  si.target_scale_factors   # scale factors vs sthol2
        si.target_sthol2 # sthol2 values  d = 0.25/sthol2**0.5
                  si.d_min_list
                  si.cc_list
                  si.low_res_cc # low-res average
      ss_b_cart_as_u_cart: anisotropic part of overall correction factor
      overall_scale: radial part of overall correction factor
    """


    results.setup_info = setup_info
    results.exclude_points_outside_density = \
        box_info.exclude_points_outside_density  # are we going to exclude them
    if return_scale_factors:
      return results

    #  Now results is a list of results. Find the good ones
    xyz_list = results.xyz_list
    d_min_list = flex.double(tuple(results.value_list))

    if xyz_list.size() == 0:
      print ("Unable to calculate local fsc map", file = self.log)
      return

    print ("D-min for overall FSC map: %.2f A " %(
      setup_info.minimum_resolution), file = self.log)
    print ("Unique values in local FSC map: %s " %(xyz_list.size()),
       file = self.log)

    x=d_min_list.min_max_mean()
    print ("Range of d_min: %.2f A to %.2f A   Mean: %.2f A " %(
      x.min, x.max, x.mean), file = self.log)

    if max_resolution_ratio is not None:
      max_value = x.min * max_resolution_ratio
      if x.max > max_value:
        s = (d_min_list > max_value)
        d_min_list.set_selected(s, max_value)
    return self._create_full_size_map_manager_with_value_list(
      xyz_list = xyz_list,
      value_list = d_min_list,
      smoothing_radius = setup_info.smoothing_radius,
      n_boxes = setup_info.n_boxes,
      small_n_real = setup_info.small_n_real,
     )

  def _create_full_size_map_manager_with_value_list(self,
      xyz_list, value_list, smoothing_radius,
      default_value = None, n_boxes = None,
      min_grid_ratio = 4,  # never bigger than 1/min_grid_ratio of full size
      small_n_real = None,
       ):
    """Create a full-size map_manager based on values at arbitrary locations
       in the map"""
    if small_n_real:
      local_n_real = small_n_real
      n_boxes = local_n_real[0]*local_n_real[1]*local_n_real[2]
    else:
      if not n_boxes:
        n_boxes = xyz_list.size()
    n_boxes = max(1,min(n_boxes,xyz_list.size()))

    # Now create a small map and fill in values
    volume_per_grid_point=self.crystal_symmetry().unit_cell(
          ).volume()/max(1,n_boxes)
    target_spacing = (volume_per_grid_point**0.33) * 0.5
    if not small_n_real:
      full_n_real = self.map_manager().map_data().all()
      local_n_real=tuple([ max(1, min(
          int((0.5+nr)/min_grid_ratio),
          int(0.5+1.5*a/target_spacing))) for
          a,nr in zip(self.crystal_symmetry().unit_cell().parameters()[:3],
           full_n_real)])
    assert value_list.size() == xyz_list.size()
    fsc_map_manager = create_map_manager_with_value_list(
       n_real = local_n_real,
       crystal_symmetry = self.crystal_symmetry(),
       value_list = value_list,
       sites_cart_list = xyz_list,
       target_spacing = target_spacing,
       default_value = default_value)

    # Get Fourier coeffs:
    map_coeffs = fsc_map_manager.map_as_fourier_coefficients()

    # Make map in full grid
    d_min_map_manager = self.map_manager(
       ).fourier_coefficients_as_map_manager(map_coeffs)
    d_min_map_manager.gaussian_filter( smoothing_radius = smoothing_radius)
    return d_min_map_manager

  def _get_d_min_from_resolution(self,resolution, d_min_ratio = 0.833):
    """Calculate value to use for d_min based on value of resolution in
    all available maps"""
    if not resolution:
      resolution = self.resolution()
    minimum_resolution = self.get_any_map_manager().resolution(
       method = 'd_min',
       set_resolution = False,
       force = True)
    return max(minimum_resolution, resolution * d_min_ratio)

  def _get_box_setup_info(self,
      map_id_1, map_id_2,
      resolution,
      d_min,
      box_cushion=None,
      n_boxes=None,
      core_box_size=None,
      smoothing_radius=None,
      box_size_ratio = 6, # full box never smaller than this ratio to resolution
      maximum_default_boxes = 2000,
      ):
    """Set up box information"""
    volume = self.crystal_symmetry().unit_cell().volume()
    if not resolution:
      resolution = self.resolution()

    if (n_boxes is not None) and (n_boxes > 0):  # n_boxes overrides
      core_box_size=int((0.5+ volume/n_boxes)**0.33)
    elif (not core_box_size):
      min_core_box_size=int((0.5+ volume/maximum_default_boxes)**0.33)
      core_box_size = max(min_core_box_size,
        int(0.5+ 3 * resolution))

    if not box_cushion:
      box_cushion = max(2.5 * resolution,
        0.5*max(0,box_size_ratio * resolution))

    core_box_size = max(
       resolution,
       box_size_ratio * resolution - 2 * box_cushion,
       core_box_size)

    n_boxes = max(1,int(0.5+volume/(core_box_size)**3))
    print ("Target core_box_size: %.2s A  Target boxes: %s Box cushion: %s" %(
        core_box_size, n_boxes, box_cushion),file = self.log)

    n_real = self.get_any_map_manager().map_data().all()
    if (not smoothing_radius):
      smoothing_radius = 0.5 * \
        self.crystal_symmetry().unit_cell().parameters()[0] * \
        min(1,core_box_size/n_real[0])
      smoothing_radius = max(2*resolution, smoothing_radius)

    small_n_real = tuple([
       max(1,min(n, int( (0.5+n)/core_box_size))) for n in n_real])
    # Working resolution is resolution * d_min_ratio
    minimum_resolution = self._get_d_min_from_resolution(resolution)
    if d_min and (minimum_resolution < d_min):
      minimum_resolution = d_min
    return group_args(
     resolution = resolution,
     box_cushion = box_cushion,
     n_boxes = n_boxes,
     small_n_real = small_n_real,
     core_box_size = core_box_size,
     smoothing_radius = smoothing_radius,
     minimum_resolution = minimum_resolution,
      )

  def _run_fsc_in_boxes(self,
     nproc = None,
     box_info = None):
    """Set up to run in each box"""

    assert box_info.n_bins is not None
    run_list=[]
    index_list=[]
    n_total = len(box_info.selection_list)
    n_in_group = int(0.5+n_total/nproc)
    for i in range(nproc):
      first_to_use = i * n_in_group + 1
      last_to_use = min(n_total,
         i * n_in_group + n_in_group )
      if i == nproc -1:
        last_to_use = n_total

      index_list.append({'i':i})
      run_list.append({'first_to_use': first_to_use,
        'last_to_use': last_to_use})

    from libtbx.easy_mp import run_parallel
    results = run_parallel(
     method = self._multiprocessing,
     qsub_command = self._queue_run_command,
     nproc = nproc,
     target_function = run_fsc_as_class(
        map_model_manager = self,
        run_list=run_list,
        box_info = box_info),
     preserve_order=False,
     kw_list = index_list)

    # Put together results
    expected_number_of_samples = len(box_info.lower_bounds_list)
    found_number_of_samples = 0
    found_number_of_samples_with_ncs = 0

    # Here each result is one scale_factor_info
    # scale_factor_info.value_list is a set of scaling_group_info objects or
    #    resolutions
    # scale_factor_info.xyz_list are the coordinates where these apply
    # scale_factor_info.n_bins is number of bins
    # value_list is a set of scaling_group_info objects, one per xyz.
    """
    scaling_group_info group_args object:
      direction_vectors: direction vectors dv for anisotropy calculations
      scaling_info_list: si (scaling_info) objects, one for each dv
        each si:  si.target_scale_factors   # scale factors vs sthol2
        si.target_sthol2 # sthol2 values  d = 0.25/sthol2**0.5
                  si.d_min_list
                  si.cc_list
                  si.low_res_cc # low-res average
      ss_b_cart_as_u_cart: anisotropic part of overall correction factor
      overall_scale: radial part of overall correction factor
    """

    all_results = None
    for scale_factor_info in results:
      if not scale_factor_info: continue
      xyz_list = scale_factor_info.xyz_list
      value_list = scale_factor_info.value_list
      found_number_of_samples += xyz_list.size()
      # Apply ncs if appropriate
      if box_info.ncs_object and box_info.ncs_object.max_operators()> 1:
        new_xyz_list = flex.vec3_double()
        new_value_list = []
        for i in range(xyz_list.size()):
          # work on one location (xyz)
          # with values: a set of scale_factor_info values, one for each
          #     direction_vector at this location
          if value_list[i] is None: continue
          new_sites,new_values = apply_ncs_to_dv_results(
            direction_vectors = box_info.direction_vectors,
            xyz = xyz_list[i],
            scaling_group_info = value_list[i],  # can be an object or a number
            ncs_object = box_info.ncs_object)
          new_xyz_list.extend(new_sites)
          new_value_list+= new_values  # n_ncs scale_factor_info objects
      else:
        new_xyz_list = xyz_list
        new_value_list = value_list
      if not all_results:
        all_results = scale_factor_info
        scale_factor_info.xyz_list = new_xyz_list
        scale_factor_info.value_list = new_value_list
      else:
        all_results.xyz_list.extend(new_xyz_list)  # vec3_double
        all_results.value_list += new_value_list   # a list
    found_number_of_samples_with_ncs = all_results.xyz_list.size()
    print ("Sampling points attempted: %s  Successful: %s  With NCS: %s" %(
      expected_number_of_samples, found_number_of_samples,
      found_number_of_samples_with_ncs), file = self.log)
    return all_results

  def local_resolution_map(self,
      map_id_1 = 'map_manager_1',
      map_id_2 = 'map_manager_2',
      map_id = 'map_manager',
      model_id = 'model',
      map_id_model_map = 'model_map',
      d_min = None,
      n_bins = 20,
      fsc_cutoff = 0.143,
      smoothing_radius = None,
      smoothing_radius_ratio = 1,
      smooth_at_end = True,
      k_sol = None,
      b_sol = None,):

    """
     Calculate local resolution map by finding resolution where local
      map correlation is fsc_cutoff at each point in the map.
     Method:  calculate local map correlation in resolution shells,
       at each point in map find the resolution where this local map
       correlation drops below fsc_cutoff.

     parameter: map_id_1:  ID of one half-map
     parameter: map_id_2:  ID of other half-map
     parameter: map_id :  ID of full map , if only one map supplied
     parameter: model_id :  ID of model, if only one map supplied
     parameter: map_id_model_map:  ID of model-map, if only one map supplied
     parameter: d_min: Finest resolution at which to calculate correlations
     parameter: n_bins: Number of resolution bins
     parameter: fsc_cutoff : value of correlation corresponding to
                             an estimated average map correlation to true map
                             of 0.5 (normally 0.143 is used)
     parameter: smoothing_radius:  radius for local correlation calculation
     parameter: smoothing_radius_ratio:  Ratio of smoothing_radius to map
                                         resolution (NOTE: to map
                                         resolution as specified in
                                         self.resolution, not to d_min).
                                         Used if smoothing_radius is None.
     parameter: smooth_at_end: smooth final local resolution map
     parameter: k_sol : k_sol for model map (if model is used)
     parameter: b_sol : b_sol for model map (if model is used)
    """

    from cctbx.maptbx.segment_and_split_map import get_smoothed_cc_map
    from cctbx.maptbx.segment_and_split_map import smooth_one_map

    hm1 = self.get_map_manager_by_id(map_id_1)
    hm2 = self.get_map_manager_by_id(map_id_2)
    full_map = self.get_map_manager_by_id(map_id)
    model = self.get_model_by_id(model_id)

    if (hm1 and hm2): # use 2 half maps as is
      pass
    elif (full_map and model):  # use full map and map from model
      hm1 = full_map
      self.generate_map(model=model,
       gridding=self.get_any_map_manager().map_data().all(),
       d_min=d_min,
       map_id = map_id_model_map,
       k_sol = k_sol,
       b_sol = b_sol)
      hm2 = self.get_map_manager_by_id(map_id_model_map)
    else:
      assert (hm1 and hm2) or (full_map and model)

    assert (hm1 and hm2)
    resolution = self.resolution()
    if d_min is None:
      d_min = self._get_d_min_from_resolution(resolution)
    if smoothing_radius is None:
      smoothing_radius = smoothing_radius_ratio * self.resolution()

    from cctbx.maptbx.segment_and_split_map import get_f_phases_from_map
    from cctbx.development.create_models_or_maps import get_map_from_map_coeffs

    crystal_symmetry = hm1.crystal_symmetry()
    mc_list = []
    for hm in [hm1, hm2]:
      mc, dummy = get_f_phases_from_map(map_data = hm.map_data(),
         crystal_symmetry = crystal_symmetry,
         d_min = d_min,
         return_as_map_coeffs = True,
         out = null_out())
      mc_list.append(mc)

    base_map_coeffs = mc_list[0]
    (d_max, d_min) = base_map_coeffs.d_max_min()
    n_bins = base_map_coeffs.safe_setup_binner(
        n_bins = n_bins, d_max = d_max, d_min = d_min)
    dsd = base_map_coeffs.d_spacings().data()
    all_bins = list(base_map_coeffs.binner().range_used())

    dsd = base_map_coeffs.d_spacings().data()

    just_above_dmin = hm1.map_data().deep_copy()
    just_above_dmin *= 0  # zero map
    just_above_cc = just_above_dmin.deep_copy()
    just_below_dmin = just_above_dmin.deep_copy()
    just_below_cc = just_above_dmin.deep_copy()
    sel       = base_map_coeffs.binner().selection(all_bins[0])
    dd         = dsd.select(sel)
    local_d_mean     = dd.min_max_mean().mean
    just_above_dmin += local_d_mean # all at high end
    just_below_dmin += d_min # all at low end
    m = hm1.deep_copy()

    n_real = hm1.map_data().all()
    for i_bin in all_bins:
      sel       = base_map_coeffs.binner().selection(i_bin)
      dd         = dsd.select(sel)
      local_d_mean     = dd.min_max_mean().mean
      map_data_list = []
      for mc in mc_list:
        mc_sel = mc.select(sel)
        map_data_list.append(get_map_from_map_coeffs(
           map_coeffs = mc_sel,
           crystal_symmetry = crystal_symmetry,
           n_real = n_real,
           apply_sigma_scaling = False))
      map_data_1, map_data_2 = map_data_list

      local_cc_map = get_smoothed_cc_map(
        map_data_1, map_data_2,
        crystal_symmetry = crystal_symmetry,
         weighting_radius = smoothing_radius)
      ss = (local_cc_map.as_1d() >= fsc_cutoff) & (
         just_above_dmin.as_1d() > local_d_mean)
      just_above_dmin.as_1d().set_selected(ss, local_d_mean)
      just_above_cc.as_1d().set_selected(ss, local_cc_map.as_1d())

      ss = (local_cc_map.as_1d() < fsc_cutoff) & (
         just_below_dmin.as_1d() < local_d_mean)
      just_below_dmin.as_1d().set_selected(ss, local_d_mean)
      just_below_cc.as_1d().set_selected(ss, local_cc_map.as_1d())

    # All set

    #  just_above_cc = smallest resolution (d-spacing) where cc >= fsc_cutoff
    # just_below_cc = largest resolution where cc < fsc_cutoff

    delta = just_above_cc - just_below_cc
    ss_neg = (delta.as_1d() <= 1.e-10)
                              # smallest d-spacing with cc >= fsc_cutoff is
                              # less than largest d_spacing with cc < fsc_cutoff
                              # happens if cc vs resolution is noisy

    delta.as_1d().set_selected(delta.as_1d() < 1.e-10, 1.e-10)
    frac = fsc_cutoff - just_below_cc
    w1 = frac / delta
    w1.as_1d().set_selected(w1.as_1d() < 0, 0)
    w1.as_1d().set_selected(w1.as_1d() > 1, 1)

    w1.as_1d().set_selected(ss_neg, 0.5) # just average if overlap
    w2 = 1 - w1
    map_data = (w1 * just_below_dmin) + (w2 * just_above_dmin)

    if smooth_at_end:
      # Finally, smooth the local resolution map
      map_data = smooth_one_map(map_data,
        crystal_symmetry = crystal_symmetry,
        smoothing_radius = smoothing_radius,
        method='exp')

    m.set_map_data(map_data)

    return m


  def map_map_fsc(self,
      map_id_1 = 'map_manager_1',
      map_id_2 = 'map_manager_2',
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
    assert n_bins is not None

    if not self.get_map_manager_by_id(map_id_1) or \
       not self.get_map_manager_by_id(map_id_2):
      return group_args(
       d_min = None,
       )
    if not resolution:
      resolution = self.resolution(
        use_fsc_if_no_resolution_available_and_maps_available=False)
    assert isinstance(resolution, (int, float))

    f_map_1, f_map_2 = self._get_map_coeffs_list_from_id_list(
      id_list = [map_id_1, map_id_2],
      mask_id = mask_id)

    bin_width=max(min_bin_width,int(0.5+f_map_1.size()/n_bins))

    # Get the FSC between map1 and map2
    fsc_curve = f_map_1.d_min_from_fsc(
        other = f_map_2, fsc_cutoff = fsc_cutoff)

    return fsc_curve

  def map_map_cc(self,
      map_id = 'map_manager_1',
      other_map_id = 'map_manager_2',
      mask_id = None,
      mask_cutoff = 0.5,
      resolution = None):

   """Calculate map-map correlation"""
   if not resolution:
     resolution = self.resolution()
   assert resolution is not None

   # resolution-filter both maps and then get cc

   map_1_id= self._generate_new_map_id()
   self.resolution_filter(
      d_min = resolution,
      map_id = map_id,
      new_map_id = map_1_id,)

   map_2_id= self._generate_new_map_id()
   self.resolution_filter(
      d_min = resolution,
      map_id = other_map_id,
      new_map_id = map_2_id,)

   map_map_info = self._get_map_map_info(
     map_id = map_1_id,
     other_map_id = map_2_id,
     mask_id = mask_id,
     mask_cutoff = mask_cutoff)
   cc = flex.linear_correlation(map_map_info.map_data_1d_1,
     map_map_info.map_data_1d_2).coefficient()
   self.remove_map_manager_by_id(map_1_id)
   self.remove_map_manager_by_id(map_2_id)
   return cc

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
   # Get mask if we set one up (excluding positions outside mask)
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


  def density_at_model_sites(self,
      map_id = 'map_manager',
      model_id = 'model',
      selection_string = None,
      model = None,
      ):
    '''
      Return density at sites in model
    '''

    if not model:
      model = self.get_model_by_id(model_id)
    else:
      model_id = '(supplied)'
    if not model:
      return None

    map_manager= self.get_map_manager_by_id(map_id)
    if not map_manager:
      raise Sorry(
     "There is no map with id='%s' available for density_at_model_sites" %(
         map_id))

    assert self.map_manager().is_compatible_model(model,
       require_match_unit_cell_crystal_symmetry=True) # model must match

    if selection_string:
      sel = model.selection(selection_string)
      model = model.select(sel)

    return map_manager.density_at_sites_cart(model.get_sites_cart())

  def map_model_cc(self,
      resolution = None,
      map_id = 'map_manager',
      model_id = 'model',
      selection_string = None,
      model = None,
      use_b_zero = True,
      ):

    """Calculate map-model correlation"""
    if not model:
      model = self.get_model_by_id(model_id)
    else:
      model_id = '(supplied)'
    if not model:
      return None
    if use_b_zero: # make a deep copy and set b values to zero
      print("Map-model CC using "+
        "map '%s' and model '%s' with B values of zero" %(map_id,model_id),
           file = self.log)
      model = model.deep_copy()
      b_iso_values = model.get_b_iso()
      model.set_b_iso(flex.double(b_iso_values.size(),0.)) # XXX should be ok
      u_cart=model.get_xray_structure().scatterers().extract_u_cart(model.get_xray_structure().crystal_symmetry().unit_cell())
      model.get_xray_structure().set_u_cart(u_cart)
      model.set_xray_structure(model.get_xray_structure())
    else:
      print("Map-model CC using "+
        "map '%s' and model '%s' with B values as is" %(map_id,model_id),
         file = self.log)

    map_manager= self.get_map_manager_by_id(map_id)
    if not map_manager:
      raise Sorry("There is no map with id='%s' available for map_model_cc" %(
         map_id))


    if not self.map_manager().is_compatible_model(model,
        require_match_unit_cell_crystal_symmetry=True): # model must match
      print("Setting model crystal symmetry to match map", file = self.log)
      model = model.deep_copy()
      model.set_crystal_symmetry(self.map_manager().crystal_symmetry())

    if not resolution:
      resolution = self.resolution()
    assert resolution is not None
    print("Resolution for map-model CC: %.3f A" %(resolution), file = self.log)

    # As five_cc is going to get the map cc without resolution cutoff on the
    #  map, do it here

    map_manager = map_manager.deep_copy()
    map_manager.resolution_filter(d_min=resolution)

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

    cc_mask = five_cc.result.cc_mask
    print("Map-model CC: %.3f" %(cc_mask), file = self.log)
    return cc_mask

  #  Methods for superposing maps

  def shift_aware_rt_to_superpose_other(self, other,
      selection_string = None):
    '''
    Identify rotation/translation to transform model from other on to
     model in this object.
    Optionally apply selection_string to both models before doing the
     mapping

    '''
    assert isinstance(other, map_model_manager)

    if selection_string:
      other_model = other.model().apply_selection_string(selection_string)
      self_model = self.model().apply_selection_string(selection_string)
    else:
      other_model = other.model()
      self_model = self.model()

    if self_model.get_sites_cart().size() == \
         other_model.get_sites_cart().size():
      # Get lsq superposition object (with r,t)
      import scitbx.math.superpose
      lsq = scitbx.math.superpose.least_squares_fit(
        reference_sites=self_model.get_sites_cart(),
        other_sites=other_model.get_sites_cart())
      other_sites_mapped = lsq.r.elems * other_model.get_sites_cart() + \
              lsq.t.elems
      starting_rmsd = self_model.get_sites_cart().rms_difference(
            other_model.get_sites_cart())
      rmsd = self_model.get_sites_cart().rms_difference(other_sites_mapped)
      print ("RMSD starting: %.3f A.  After superposition: %.3f A " %(
          starting_rmsd,rmsd), file=self.log)
    else: # use superpose_pdbs tool to try and get superposition
      try:
        from phenix.command_line import superpose_pdbs
        params = superpose_pdbs.master_params.extract()
        x = superpose_pdbs.manager(
          params,
          log = null_out(),
          write_output = False,
          save_lsq_fit_obj = True,
          pdb_hierarchy_fixed = self_model.get_hierarchy(),
          pdb_hierarchy_moving = other_model.get_hierarchy().deep_copy(),)
        lsq = x.lsq_fit_obj
        del x

      except Exception as e:
        print ("Unable to superpose other on self..", file = self.log)
        return None


    working_rt_info = group_args(
      r=lsq.r,
      t=lsq.t)

    shift_aware_rt_info = self.shift_aware_rt(
          working_rt_info=working_rt_info,
          from_obj = other,
          to_obj = self)
    return shift_aware_rt_info

  def superposed_map_manager_from_other(self,other,
     working_rt_info = None,
     absolute_rt_info = None,
     shift_aware_rt_info = None,
     selection_string = None):
    '''
    Identify rotation/translation to transform model from other on to
      model in this object.
    Optionally apply selection_string to both models before doing the
     mapping
    Then extract map from other to cover map in this object,
    Fill in with zero where undefined if wrapping is False.

    Allow specification of working_rt (applies to working coordinates in
      other and self), or absolute_rt_info (applies to absolute, original
      coordinates)

    '''

    # get the shift_aware_rt_info if not supplied
    if not shift_aware_rt_info:
      if absolute_rt_info:
        shift_aware_rt_info = self.shift_aware_rt(
          absolute_rt_info=absolute_rt_info)
      elif working_rt_info:
        shift_aware_rt_info = self.shift_aware_rt(
          working_rt_info=working_rt_info,
          from_obj = other,
          to_obj = self)
      else:
        shift_aware_rt_info = self.shift_aware_rt_to_superpose_other(other,
            selection_string = selection_string)

    rt_info = shift_aware_rt_info.working_rt_info(from_obj=other, to_obj=self)

    # Extract the other map in defined region (or all if wrapping = True)
    # Wrapping = True:  just pull from other map

    # Wrapping = False  Zero outside defined region
    #  Make a big map_model_manager for other that includes the entire
    #  region corresponding
    #  to this map.  When constructing that map, set undefined values to zero
    #  Then just pull from this big map_model_manager
    if other.map_manager().wrapping():
      other_to_use = other
    else:
      print ("Making a large version of other map where values are zero if"+
       " not defined", file = self.log)
      # other_to_use = larger_map...
      lower_bounds, upper_bounds= self._get_bounds_of_rotated_corners(
        other, rt_info)
      other_to_use=other.extract_all_maps_with_bounds(
        lower_bounds,
        upper_bounds)
      print ("Done making version of other map where values are zero if"+
       " not defined", file = self.log)

    # Ready to extract from this box with interpolation
    rt_info = shift_aware_rt_info.working_rt_info(
       from_obj=other_to_use, to_obj=self)
    r_inv = rt_info.r.inverse()
    t_inv = -r_inv*rt_info.t

    from cctbx.maptbx import superpose_maps
    superposed_map_data = superpose_maps(
      unit_cell_1        = other_to_use.crystal_symmetry().unit_cell(),
      unit_cell_2        = self.crystal_symmetry().unit_cell(),
      map_data_1         = other_to_use.map_manager().map_data(),
      n_real_2           = self.map_manager().map_data().focus(),
      rotation_matrix    = r_inv.elems,
      translation_vector = t_inv.elems,
      wrapping           = False)

    new_mm = self.map_manager().customized_copy(
      map_data = superposed_map_data)
    new_mm.set_wrapping(False) # always
    return new_mm

  def _get_bounds_of_rotated_corners(self, other, rt_info):
    '''
    Return info object with lower_bounds and upper_bounds in this map
    corresponding to the lowest and highest values of coordinates obtained
     by applying the inverse of rt_info to the corners of the map in other.
    '''

    r_inv = rt_info.r.inverse()
    t_inv = -r_inv*rt_info.t

    self_all = self.map_data().all()
    other_all = other.map_data().all()
    uc = self.crystal_symmetry().unit_cell().parameters()[:3]
    other_uc = other.crystal_symmetry().unit_cell().parameters()[:3]
    other_xyz_list = flex.vec3_double()
    for i in [0,self_all[0]]:
      x = uc[0]*i/self_all[0]
      for j in [0,self_all[1]]:
        y = uc[1]*j/self_all[1]
        for k in [0,self_all[2]]:
          z = uc[2]*k/self_all[2]
          xyz = col((x,y,z))
          other_xyz_list.append(r_inv * xyz + t_inv)
    min_xyz=other_xyz_list.min()
    max_xyz=other_xyz_list.max()
    # Bounds at least one beyond any point that could be asked for
    new_low_ijk =tuple([int(-2+xx * ii/aa) for xx, ii,aa in zip(
        min_xyz,other_all, other_uc)])
    new_high_ijk =tuple([int(2+xx * ii/aa) for xx, ii, aa in zip(
        max_xyz,other_all,other_uc)])
    return new_low_ijk,new_high_ijk


  # General methods

  def model_from_text(self,
    text,
    return_as_model = False,
    model_id = 'model_from_text'):
    '''
     Convenience method to convert txt into a model, where the
     model has symmetry and shift cart matching this manager
    '''

    from mmtbx.model import manager as model_manager
    from iotbx.pdb import input
    inp = input(source_info = 'text', lines=flex.split_lines(text))
    model = model_manager(
      model_input = inp,
      crystal_symmetry = self.crystal_symmetry(),
      log = null_out())
    self.set_model_symmetries_and_shift_cart_to_match_map(model)
    if return_as_model:
      return model
    else:
      self.add_model_by_id(model = model, model_id=model_id)

  def model_from_hierarchy(self,
    hierarchy,
    return_as_model = False,
    model_id = 'model_from_hierarchy',
    use_pdb_input = False):
    '''
     Convenience method to convert a hierarchy into a model, where the
     model has symmetry and shift cart matching this manager
    '''
    model = hierarchy.as_model_manager(
        crystal_symmetry = self.crystal_symmetry())

    self.set_model_symmetries_and_shift_cart_to_match_map(model)
    if return_as_model:
      return model
    else:
      self.add_model_by_id(model = model, model_id=model_id)


  def model_from_sites_cart(self,
    sites_cart,
    model_id = 'model_from_sites',
    return_model = None,
    **kw):
    '''
      Convenience method to use the from_sites_cart method of model manager to
      create a model.

      Model is saved with the id of model_id (default 'model_from_sites')

      Unique feature of this method is that it sets up the crystal_symmetry,
        unit_cell_crystal_symmetry, and shift_cart in the model to match
        this map_model_manager.

      Note:  sites_cart are assumed to be in the same coordinate frame as this
        map_model_manager.

      Any keywords used in from_sites_cart are allowed (i.e.,
        atom_name and scatterer or atom_name_list and scatterer_list
        occ or occ_list
        b_iso or b_iso_list
        resname or resname_list
        resseq_list
    '''


    from mmtbx.model import manager as model_manager
    model = model_manager.from_sites_cart(
       sites_cart = sites_cart,
       crystal_symmetry = self.crystal_symmetry(),
       **kw)

    self.set_model_symmetries_and_shift_cart_to_match_map(model)

    if return_model:
      return model
    else:
      self.add_model_by_id(model = model, model_id=model_id)

  def set_model_symmetries_and_shift_cart_to_match_map(self,
     model):
    '''
      Method to set the crystal_symmetry, unit_cell_crystal_symmetry and
      shift_cart of any model to match those of this map_model_manager

      This method changes the model in place.

      Assumes that the sites in this model are already in the same coordinate
      system as this map_model_manager

      This is the preferred method of making sure that a model has the
      same symmetry and shifts as a map_model_manager if the model coordinates
      already match the map_model_manager.

      If you have a model that is relative to some other coordinate system,
      first make sure that model has a complete set of symmetry and shifts
      for that coordinate system (any model coming from a map_model_manager
      will have these).  Then use either self.shift_any_model_to_match(model) to
      shift that model directly to match this map_model_manager, or use
      self.add_model_by_id(model_id=model_id, model=model) which will import
      that model into this manager and shift it to match in the process.

    '''

    self.map_manager().set_model_symmetries_and_shift_cart_to_match_map(model)



  def remove_origin_shift_and_unit_cell_crystal_symmetry(self):
    '''
       Set the origin of the current map to be (0,0,0) and reset the
       unit_cell_crystal_symmetry to match the current crystal_symmetry.

       Do not use this method if you want to use the map_model_manager or any
       of its maps and models in the usual way.  This is a method to create
       a set of maps and model that have no reference to their original
       locations.

       Typical use:
         box_mmm = mmm.extract_all_maps_around_model()
         box_mmm.remove_origin_shift_and_unit_cell_crystal_symmetry()
         model_file = data_manager.write_model_file(box_mmm.model(), model_file)
         data_manager.write_real_map_file(box_mmm.map_manager(), map_file)
       Now model_file and map_file have no origin shift and the model CRYST1
       matches the crystal_symmetry and map_file unit_cell_crystal_symmetry of
       map_file

    '''

    assert self.map_data() is not None
    # Target origin (zero)
    zero_origin = (0,0,0)
    # Gridding of full unit cell (same as current map data)
    box_gridding = self.map_manager().map_data().all()

    # Reset origin and gridding for each map_manager
    for mm in self.map_managers():
      mm.set_original_origin_and_gridding(
       original_origin = zero_origin,
       gridding = box_gridding)

    # Any map_manager
    mm = self.map_manager()
    # Reset origin for each model
    for m in self.models():
      mm.set_model_symmetries_and_shift_cart_to_match_map(m)
    # All done


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


  def _generate_new_map_id(self, prefix = 'temp'):
    '''
     Create a unique map_id
    '''
    used_id_list = self.map_id_list()
    i = 0
    while True:
      i += 1
      id = "%s_%s" %(prefix,i)
      if not id in used_id_list:
        return id

  def _generate_new_model_id(self,prefix = 'temp'):
    '''
     Create a unique model_id
    '''
    used_id_list = self.model_id_list()
    i = 0
    while (True):
      id = "%s_%s" %(prefix,i)
      if not id in used_id_list:
        return id

  def warning_message(self):
    """Return the warning message"""
    return self._warning_message

  def show_summary(self, log = sys.stdout):
    '''Show short summary of this map_model_manager'''
    text = self.__repr__()
    print (text, file = log)

  # Methods for accessing map_data, xrs, hierarchy directly.
  #  Perhaps remove all these

  def map_data(self):
    '''Return map_data from map_manager'''
    if self.map_manager() and (not self.map_manager().is_dummy_map_manager()):
      return self.map_manager().map_data()

  def map_data_1(self):
    '''Return map_data from map_manager_1'''
    if self.map_manager_1() and (
       not self.map_manager_1().is_dummy_map_manager()):
      return self.map_manager_1().map_data()

  def map_data_2(self):
    '''Return map_data from map_manager_2'''
    if self.map_manager_2() and (
       not self.map_manager_2().is_dummy_map_manager()):
      return self.map_manager_2().map_data()

  def map_data_list(self):
    '''Return list of map_data, one from from each map_manager'''
    map_data_list = []
    for mm in self.map_managers():
      map_data_list.append(mm.map_data())
    return map_data_list

  def xray_structure(self):
    '''Return xray_structure from working model'''
    if(self.model() is not None):
      return self.model().get_xray_structure()
    else:
      return None

  def hierarchy(self):
    '''Return the hierarchy from working model'''
    return self.model().get_hierarchy()


  # Methods to be removed

  def get_counts_and_histograms(self):
    '''Calculate summary information about maps and histograms of values'''
    self._counts = get_map_counts(
      map_data         = self.map_data(),
      crystal_symmetry = self.crystal_symmetry())
    self._map_histograms = get_map_histograms(
        data    = self.map_data(),
        n_slots = 20,
        data_1  = self.map_data_1(),
        data_2  = self.map_data_2())

  def counts(self):
    '''Return summary information about maps and histograms of values'''
    if not hasattr(self, '_counts'):
      self.get_counts_and_histograms()
    return self._counts

  def histograms(self):
    """Return map histograms"""
    if not hasattr(self, '_map_histograms'):
      self.get_counts_and_histograms()
    return self._map_histograms

  #  Convenience methods

  def mask_info(self,
    mask_id = 'mask',
    cutoff = 0.5,
    quiet = False,
    ):
    '''  Summarizes info about mask and returns group_args'''

    mask_mm = self.get_map_manager_by_id(map_id = mask_id)
    if not mask_mm:
      return None
    mask_info = self.map_info(quiet=True, map_id = mask_id) # basic map info

    # Change title
    mask_info.group_args_type = "Summary of mask info about mask '%s' " %(
      mask_id)

    # Add mask-specific information
    mask_info.cutoff = cutoff
    mask_info.is_mask = mask_mm.is_mask()
    mask_info.marked_points = (mask_mm.map_data() >= cutoff ).count(True)
    mask_info.fraction_marked = mask_info.marked_points/max(1.,
        mask_mm.map_data().size())

    if not quiet:
      print (mask_info, file = self.log)
    return mask_info

  def map_info(self,
    map_id = 'map_manager',
    sigma_cutoff = 1,
    quiet = False,
    ):
    '''  Summarizes info about map and returns group_args'''

    map_mm = self.get_map_manager_by_id(map_id = map_id)
    if not map_mm:
      return None

    mmm = map_mm.map_data().as_1d().min_max_mean()
    standard_deviation = map_mm.map_data(
       ).as_1d().sample_standard_deviation()
    cutoff = mmm.mean + sigma_cutoff * standard_deviation

    points_above_sigma_cutoff = (map_mm.map_data() >= cutoff ).count(True)

    map_info = group_args(
      group_args_type = "Summary of map info about map '%s' " %(map_id),
      cutoff = cutoff,
      mean = mmm.mean,
      min = mmm.min,
      max = mmm.max,
      standard_deviation = standard_deviation,
      size = map_mm.map_data().size(),
      points_above_sigma_cutoff = points_above_sigma_cutoff,
      fraction_above_sigma_cutoff= points_above_sigma_cutoff/max(
         1., map_mm.map_data().size()),
      absolute_center_cart = map_mm.absolute_center_cart(),
      working_center_cart = tuple(flex.double(map_mm.absolute_center_cart()) +
          flex.double(map_mm.shift_cart())),
      )
    if not quiet:
      print (map_info, file = self.log)
    return map_info


  def shift_aware_rt(self,
     from_obj = None,
     to_obj = None,
     working_rt_info = None,
     absolute_rt_info = None):
   '''
   Returns shift_aware_rt object

   Uses rt_info objects (group_args with members of r, t).

   Simplifies keeping track of rotation/translation between two
    objects that each may have an offset from absolute coordinates.

   absolute rt is rotation/translation when everything is in original,
      absolute cartesian coordinates.

   working_rt is rotation/translation of anything in "from_obj" object
      to anything in "to_obj" object using working coordinates in each.

   Usage:
   shift_aware_rt = self.shift_aware_rt(absolute_rt_info = rt_info)
   shift_aware_rt = self.shift_aware_rt(working_rt_info = rt_info,
      from_obj=from_obj, to_obj = to_obj)

   apply RT using working coordinates in objects
   sites_cart_to_obj = shift_aware_rt.apply_rt(sites_cart_from_obj,
      from_obj=from_obj, to_obj=to_obj)

   apply RT absolute coordinates
   sites_cart_to = shift_aware_rt.apply_rt(sites_cart_from)

   '''
   from iotbx.map_manager import shift_aware_rt

   return shift_aware_rt(
     from_obj = from_obj,
     to_obj = to_obj,
     working_rt_info = working_rt_info,
     absolute_rt_info = absolute_rt_info)

  def generate_map(self,
      d_min = None,
      origin_shift_grid_units = None,
      file_name = None,
      model = None,
      n_residues = None,
      b_iso = 30,
      k_sol = None,
      b_sol = None,
      box_cushion = 5,
      scattering_table = None,
      fractional_error = 0.0,
      gridding = None,
      wrapping = False,
      map_id = None,
      f_obs_array = None,
      resolution_factor = None,
     ):
    """
      Simple interface to cctbx.development.generate_map allowing only
      a small subset of keywords. Useful for quick generation of models, map
      coefficients, and maps

      For full functionality use cctbx.development.generate_model,
      cctbx.development.generate_map_coeffs, and
      cctbx.development.generate_map

      Summary

      If no map_manager is present, use supplied or existing model to
         generate map_manager and model.

      If map_manager is present and is not a dummy map_manager,
         use supplied or existing model as model and
         create new entry in this this map_model_manager with name map_id.
         If map_id is None, use 'model_map'

      If no existing or supplied model, use default model from library,
      box with box_cushion around it and choose n_residues to
      include (default=10).

      Parameters:

      model (model.manager object, None):    model to use (as is)
      file_name (path , None):    file containing coordinates to use (instead
                                  of default model)
      n_residues (int, 10):      Number of residues to include (from default
                                  model or file_name)
      b_iso (float, 30):         B-value (ADP) to use for all atoms if model
                                 is not supplied
      box_cushion (float, 5):     Buffer (A) around model
      d_min (float, 3):      high_resolution limit (A)
      gridding (tuple (nx, ny, nz), None):  Gridding of map (optional)
      origin_shift_grid_units (tuple (ix, iy, iz), None):  Move location of
          origin of resulting map to (ix, iy, iz) before writing out
      wrapping:  Defines if map is to be specified as wrapped
      resolution_factor:  Defines ratio of resolution to gridding if gridding is not set
      scattering_table (choice, 'electron'): choice of scattering table
           All choices: wk1995 it1992 n_gaussian neutron electron
      fractional_error:  resolution-dependent fractional error, ranging from
           zero at low resolution to fractional_error at d_min. Can
           be more than 1.
      map_id:  ID of map_manager to be created with model-map information (only
                 applies if there is an existing map_manager)
    """

    # See if we have a map_manager
    if (not self.map_manager()) or (
        self.map_manager() and self.map_manager().is_dummy_map_manager()):
      have_map_manager = False
    else:
      have_map_manager = True


    #  Choose scattering table
    if not scattering_table:
      scattering_table = self.scattering_table()

    # Set the resolution now if not already set
    if d_min is not None and have_map_manager:
      self.set_resolution(d_min)
    elif d_min is None and self.resolution():
      d_min = self.resolution()
    elif d_min is None:
      d_min = 3

    self._print("\nGenerating new map data\n")
    if self.model() and (not model):
      self._print("NOTE: using existing model to generate map data\n")
      model = self.model()

    if not map_id:
      map_id = 'model_map'
      self._print("Generated map will go in %s" %(map_id))

    if have_map_manager:
      if not gridding:
        gridding = self.map_manager().map_data().all()
        origin_shift_grid_units = self.map_manager().origin_shift_grid_units
        self._print(
          "Using existing map_manager as source of gridding and origin")
      self._print("Model map in map_model_manager "+
         "'%s' will be placed in map_manager '%s'" %(self.name,map_id))

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
        log = null_out())

    if have_map_manager:  #  make sure model matches
      if not self.map_manager().is_compatible_model(model,
          require_match_unit_cell_crystal_symmetry=True):
         self.shift_any_model_to_match(model,
         set_unit_cell_crystal_symmetry = True)

    map_coeffs = generate_map_coefficients(model = model,
        d_min = d_min,
        k_sol = k_sol,
        b_sol = b_sol,
        scattering_table = scattering_table,
        f_obs_array = f_obs_array,
        log = null_out())
    mm = generate_map_data(
      map_coeffs = map_coeffs,
      d_min = d_min,
      gridding = gridding,
      wrapping = wrapping,
      resolution_factor = resolution_factor,
      origin_shift_grid_units = origin_shift_grid_units,
      high_resolution_real_space_noise_fraction = fractional_error,
      log = null_out())

    if have_map_manager and self.get_any_map_manager():
      new_mm = self.get_any_map_manager().customized_copy(
        map_data=mm.map_data())
      self.add_map_manager_by_id(new_mm,map_id)
    else: # create map-model manager info
      self.set_up_map_dict(map_manager=mm)
      self.set_up_model_dict(model=model)
      if map_id is not None and map_id != 'map_manager':  # put it in map_id too
        new_mm = self.get_any_map_manager().customized_copy(
          map_data=mm.map_data())
        self.add_map_manager_by_id(new_mm,map_id)

    # Set the resolution if we didn't already
      self.set_resolution(d_min)

  def _empty_copy(self):
    '''
      Return a copy with no data
    '''
    new_mmm = map_model_manager(log = self.log)
    new_mmm._map_dict={}
    new_mmm._model_dict={}
    return new_mmm

  def deep_copy(self):
    '''
      Return a deep_copy of this map_manager
      Use customized copy with default map_dict and model_dict (from self)
    '''
    return self.customized_copy(map_dict = None, model_dict = None,
      name = "%s_deep_copy" %(self.name))

  def customized_copy(self, model_dict = None, map_dict = None,
      name = None):
    '''
      Produce a copy of this map_model object, replacing nothing,
      maps or models, or both
    '''

    # Decide what is new

    if model_dict is not None: # take new model_dict without deep_copy
      new_model_dict = model_dict
    else:  # deep_copy existing model_dict
      new_model_dict = {}
      for id in self.model_id_list():
        new_model_dict[id]=self.get_model_by_id(id).deep_copy()

    if map_dict is not None: # take new map_dict without deep_copy
      new_map_dict = map_dict
    else:  # deep_copy existing map_dict
      new_map_dict = {}
      for id in self.map_id_list():
        new_map_dict[id]=self.get_map_manager_by_id(id).deep_copy()

    # Build new map_manager

    new_mmm = map_model_manager()

    new_mmm._model_dict = new_model_dict
    new_mmm._map_dict = new_map_dict
    new_mmm._info = deepcopy(self._info)

    self._set_default_parameters(new_mmm, name = name)

    return new_mmm


  def shift_origin_to_match_original(self):
    '''  Shift the origin of all maps and models to match original'''

    # First shift the maps
    for mm in self.map_managers():
      mm.shift_origin_to_match_original()

    # And now models to match
    for model in self.models():
      self.get_any_map_manager().shift_model_to_match_map(model)

  def check_consistency(self, stop_on_errors = True, print_errors = True,
        absolute_angle_tolerance = None,
        absolute_length_tolerance = None,
        shift_tol = None):
    ''' Check that all component objects have the same crystal_symmetry,
      unit_cell_crystal_symmetry, and shift_cart (if these are part of the
      objects).'''

    if absolute_angle_tolerance is None:
      absolute_angle_tolerance = self._absolute_angle_tolerance
    if absolute_length_tolerance is None:
      absolute_length_tolerance = self._absolute_length_tolerance
    if shift_tol is None:
      shift_tol = self._shift_tol

    object_list = self.models() + self.map_managers()
    if self.ncs_object():
      object_list.append(self.ncs_object())

    object_list_with_self = object_list + [self]
    # Check for consistency among all objects and also include self here

    ok = all_objects_have_same_symmetry_and_shift_cart(object_list_with_self,
        absolute_angle_tolerance = absolute_angle_tolerance,
        absolute_length_tolerance = absolute_length_tolerance,
        shift_tol = shift_tol,
        print_errors = False)

    target = self.map_manager()
    if target:
      for mm in self.map_managers():
        if (not mm.is_similar(target)):
          ok = False
          text = "Map managers %s and %s are not compatible" %(mm, target)

    if (not ok):
      text = "Consistency check failure in map_model_manager '%s'" %(
           self.name)
      if print_errors:
        self._print(text)
        all_objects_have_same_symmetry_and_shift_cart(object_list_with_self,
          absolute_angle_tolerance = absolute_angle_tolerance,
          absolute_length_tolerance = absolute_length_tolerance,
          shift_tol = shift_tol,
          print_errors = True)
      if (stop_on_errors):
        raise AssertionError(text)

    # Check for consistency inside each object that has a check_consistency()
    #   method

    for object in object_list:
      if hasattr(object, 'check_consistency'):
        object.check_consistency(stop_on_errors = stop_on_errors,
          print_errors = print_errors,
          absolute_angle_tolerance = absolute_angle_tolerance,
          absolute_length_tolerance = absolute_length_tolerance,
          shift_tol = shift_tol)

    return ok

  def model_building(self,
     nproc = None,
     soft_zero_boundary_mask = True,
     soft_zero_boundary_mask_radius = None,
     model_id = 'model',
     normalize = True,
     ):
    '''
     Return this object as a local_model_building object
     The model-building object has pointers to model and map_manager, not
       copies
      resolution is resolution for Fourier coefficients
      is_xray_map is True for x-ray map
      nproc is number of processors to use

      If no map_manager is present, resolution, experiment type are not
      required
    '''

    if self.map_manager() and (not self.map_manager().is_dummy_map_manager()):
      resolution = self.resolution()
      assert resolution is not None

    if not nproc:
      nproc = self.nproc()

    from phenix.model_building import local_model_building
    mb = local_model_building(
     map_model_manager = self, # map_model manager
     soft_zero_boundary_mask = soft_zero_boundary_mask,
     soft_zero_boundary_mask_radius = soft_zero_boundary_mask_radius,
     nproc= nproc,
     model_id = model_id,
     normalize = normalize,
     log = self.log,
    )
    mb.set_defaults(debug = self.verbose)
    return mb

  def as_map_model_manager(self):
    '''
      Return this object (allows using .as_map_model_manager() on both
      map_model_manager objects and others including box.around_model() etc.

    '''
    return self


  def as__match_map_model_ncs(self):
    '''
      Return this object as a _match_map_model_ncs

      Includes only the map_manager and model and ncs object, ignores all
      other maps and models (_match_map_model_ncs takes only one of each).

      Note two underscores because the module it is returning is
      _match_map_model_ncs.

    '''
    from iotbx.map_model_manager import _match_map_model_ncs
    mmmn = _match_map_model_ncs()
    if self.map_manager():
      mmmn.add_map_manager(self.map_manager())
    if self.model():
      mmmn.add_model(self.model())
    if self.ncs_object():
      mmmn.add_ncs_object(self.ncs_object())
    return mmmn


class _match_map_model_ncs(object):
  '''
   _match_map_model_ncs

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
    import io
    pickle_dict = self.__dict__.copy()
    if isinstance(self.log, io.TextIOWrapper):
      pickle_dict['log'] = None
    return pickle_dict

  def __setstate__(self, pickle_dict):
    self.__dict__ = pickle_dict
    if self.log is None:
      self.log = sys.stdout

  def deep_copy(self):
    new_mmmn = _match_map_model_ncs()
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

  def _print(self, m, force = False):
    if (self.log is not None) and hasattr(self.log, 'closed') and (
        not self.log.closed):
      print(m, file = self.log)
    elif force:
      print(m)

  def write_map(self, file_name = None):
    if not self._map_manager:
      self._print ("No map to write out")
    elif not file_name:
      self._print ("Need file name to write map")
    else:
      self._map_manager.write_map(file_name = file_name)

  def write_model(self,
     file_name = None, data_manager = None, format = None):
    if not self._model:
      self._print ("No model to write out")
    elif not file_name:
      self._print ("Need file name to write model")
    else:
      # Write out model
      if not data_manager:
        from iotbx.data_manager import DataManager
        data_manager = DataManager()
        data_manager.set_overwrite(True)
      file_name = data_manager.write_model_file(self._model, file_name,
        format = format)
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
        model=self.model()
        self.map_manager().set_model_symmetries_and_shift_cart_to_match_map(
          self.model())  # modifies self.model() in place
        model=self.model()
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
    #  or at least be compatible

    assert self.map_manager() is not None
    self.map_manager().set_ncs_object(ncs_object) # checks for shift_cart

  def read_map(self, file_name):
    # Read in a map and make sure its symmetry is similar to others
    mm = MapManager(file_name)
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
    if self._model:
      self._map_manager.set_model_symmetries_and_shift_cart_to_match_map(
        self._model)

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

    # Shift origin of map_manager. Note if there was any problem
    self._map_manager._warning_message = ''
    self._map_manager.shift_origin(desired_origin = desired_origin)
    if self._map_manager.warning_message():
      self._print("\n%s\n" %(self._map_manager.warning_message()), force=True)

    # Shift origin of model  Note this sets model shift_cart
    if self._model:
      self._model = self.shift_model_to_match_working_map(
        coordinate_shift = shift_to_apply_cart,
        new_shift_cart = new_full_shift_cart,
        final_crystal_symmetry = self._map_manager.crystal_symmetry(),
        final_unit_cell_crystal_symmetry =
           self._map_manager.unit_cell_crystal_symmetry(),
        model = self._model)

  def shift_ncs_to_match_working_map(self, ncs_object = None,
    reverse = False,
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
     new_shift_cart = None,
    final_crystal_symmetry = None,
    final_unit_cell_crystal_symmetry = None):

    '''
    Shift a model based on the coordinate shift for the working map.

    Also match the crystal_symmetry and unit_cell_crystal_symmetry
      of the model to the map, unless specified as final_crystal_symmetry
      and final_unit_cell_crystal_symmetry.

    Optionally specify the shift to apply (coordinate shift) and the
    new value of the shift recorded in the model (new_shift_cart)
    '''

    if final_crystal_symmetry is None:
      final_crystal_symmetry = self.crystal_symmetry()
    if final_unit_cell_crystal_symmetry is None:
      final_unit_cell_crystal_symmetry = self.unit_cell_crystal_symmetry()

    if coordinate_shift is None:
      coordinate_shift = self.get_coordinate_shift(
       reverse = reverse)
    if new_shift_cart is None:
      new_shift_cart = coordinate_shift


    model.shift_model_and_set_crystal_symmetry(shift_cart = coordinate_shift,
      crystal_symmetry = final_crystal_symmetry)

    # Allow specifying the final shift_cart:
    if tuple(new_shift_cart) !=  tuple(coordinate_shift):
      model.set_shift_cart(new_shift_cart)

    return model

  def shift_model_to_match_original_map(self, model = None):
    # Shift a model object to match the original map (based
    #    on -self._map_manager.origin_shift_grid_units)
    return self.shift_model_to_match_working_map(model = model, reverse = True,
      final_crystal_symmetry = self.unit_cell_crystal_symmetry(),
      final_unit_cell_crystal_symmetry = self.unit_cell_crystal_symmetry())

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

#   Misc methods

def all_objects_have_same_symmetry_and_shift_cart(object_list,
        absolute_angle_tolerance = 0.01,
        absolute_length_tolerance = 0.01,
        shift_tol = 0.001,
        print_errors = False):

  '''  Return True if all the objects have same crystal_symmetry(),
    unit_cell_crystal_symmetry(), and shift_cart() if attributes are present.

    Note: For model, map_model_manager, and map_manager, all of these
    attributes are present; for ncs_object, only shift_cart'''

  attribute_list = [
     'crystal_symmetry',
     'unit_cell_crystal_symmetry',
     'shift_cart']

  for key in attribute_list:

    value_list = []  # all the values must match
    for object in object_list:
      if hasattr(object, key):  # only include values for attributes present
        value_list.append(getattr(object,key)())

    if None in value_list:  # all must be None or all must be not None
      if value_list.count(None) != len(value_list):
        return False  # values do not all match

    target = value_list[0]
    for other in value_list[1:]:
      if key in ['crystal_symmetry','unit_cell_crystal_symmetry']:
        if (not target.is_similar_symmetry(other,
          absolute_angle_tolerance = absolute_angle_tolerance,
          absolute_length_tolerance = absolute_length_tolerance,
            )):
          if print_errors:
            print("\nThe attribute '%s' does not match between \n%s and \n%s" %(
              key, target, other))
          return False

      elif key in ['shift_cart']:
        if (not is_same_shift_cart(target, other, tol = shift_tol)):
          if print_errors:
            print("\nThe attribute '%s' does not match between \n%s and \n%s" %(
              key, target, other))
          return False
      else:
        raise AssertionError(
         "Missing key in all_objects_have_same_symmetry_and_shift_cart")
  return True

def is_same_shift_cart(shift, other_shift, tol = 0.001):
    '''Compare shift to other_shift'''

    this_shift=flex.double(shift)
    other_shift=flex.double(other_shift)
    delta=this_shift - other_shift
    mmm=delta.min_max_mean()
    if mmm.min < -tol or mmm.max > tol: # shifts do not match
      return False
    else:
      return True


def convert_tlso_group_info_to_lists(tlso_group_info):
      #tlso_group_info.tlso_selection_list,
      # tlso_group_info.tlso_shift_cart_list,):
      """Convert from tlso group info object to lists"""
      tlso_group_info.T_list = []
      tlso_group_info.L_list = []
      tlso_group_info.S_list = []
      tlso_group_info.O_list = []
      for tlso in tlso_group_info.tlso_list:
        tlso_group_info.T_list.append(tlso.t)
        tlso_group_info.L_list.append(tlso.l)
        tlso_group_info.S_list.append(tlso.s)
        tlso_group_info.O_list.append(tlso.origin)
      tlso_group_info.tlso_list = None
      tlso_group_info.tlso_selection_list = \
         tlso_group_info.selection_as_text_list
      tlso_group_info.tlso_shift_cart_list = len(
         tlso_group_info.tlso_selection_list) *[None]

def apply_aniso_b_cart_to_f_array_info(f_array_info,
         b_iso, d_min, aniso_b_cart):
       """Apply anisotropic B factors to a structure factor array"""
       from cctbx.maptbx.refine_sharpening import analyze_aniso_object
       analyze_aniso = analyze_aniso_object()
       analyze_aniso.set_up_aniso_correction(f_array=f_array_info.f_array,
         b_iso = b_iso,
         d_min = d_min,
         b_cart_to_remove=aniso_b_cart)
       new_array  = analyze_aniso.apply_aniso_correction(
         f_array=f_array_info.f_array)
       sel = flex.bool(new_array.data().size(),True)
       f_array_info.f_array.data().set_selected(sel,new_array.data())

def get_average_scale_factors(scale_factor_info):
  """Get average values of scale factors"""
  average_scale_factors = None
  n = 0
  n_bins = None
  for sgi in scale_factor_info.value_list:
    if sgi and hasattr(sgi,'overall_scale') and sgi.overall_scale:
       if not average_scale_factors:
         average_scale_factors = flex.double(sgi.overall_scale.size(),0)
       average_scale_factors += sgi.overall_scale
       n+=1
    elif sgi and hasattr(sgi,'scaling_info_list') and sgi.scaling_info_list \
      and hasattr(sgi.scaling_info_list[0],'target_sthol2'):
     n_bins = sgi.scaling_info_list[0].target_sthol2.size()
  if n == 0:
    if n_bins:
      return flex.double(n_bins, 1.0)
    else:
      return None
  else:
    average_scale_factors = average_scale_factors/n
    return average_scale_factors

def get_tlso_group_info_from_model(model, nproc = 1, log = sys.stdout):
  ''' Extract tlso_group_info from aniso records in model'''
  print("\nExtracting tlso_group_info from model", file = log)
  model_use=model.deep_copy()  # this routine overwrites...


  xrs=model_use.get_xray_structure()
  u_iso  = xrs.extract_u_iso_or_u_equiv() # to check for positivity

  sel1 = (u_iso >= 0)
  model_use = model_use.select(sel1)

  xrs=model_use.get_xray_structure()
  uc = xrs.unit_cell()
  sites_cart = xrs.sites_cart()
  u_cart = xrs.scatterers().extract_u_cart(uc)
  u_iso  = xrs.extract_u_iso_or_u_equiv() # to check for positivity
  pdb_hierarchy = model_use.get_hierarchy()

  ok = False
  for u in u_cart:
    if tuple(u) != (0,0,0,0,0,0) and tuple(u) != (-1,-1,-1,-1,-1,-1):
      ok = True
      break
  if not ok:
    raise Sorry(
      "Aniso U values from model are all zero or missing...cannot get TLS")

  # Get the groups in this file
  from mmtbx.command_line.find_tls_groups import find_tls, master_phil
  working_phil = master_phil.fetch()
  params = working_phil.extract()
  params.nproc = nproc

  selections = find_tls(
    params,
    pdb_hierarchy,
    xrs,
    return_as_list=True,
    ignore_pdb_header_groups=False,
    out=log)

  from mmtbx.tls import tools
  T_list = []
  L_list = []
  S_list = []
  O_list = []
  tlso_selection_list= []
  tlso_shift_cart_list = []

  print("\nTLS GROUPS FOUND ", file = log)

  for selection in selections:
    sel = model_use.selection(selection)
    cm = sites_cart.select(sel).mean()

    result = tools.tls_from_uaniso_minimizer(
      uaniso         = u_cart.select(sel),
      T_initial      = [0,0,0,0,0,0],
      L_initial      = [0,0,0,0,0,0],
      S_initial      = [0,0,0,0,0,0,0,0,0],
      refine_T       = True,
      refine_L       = True,
      refine_S       = True,
      origin         = cm,
      sites          = sites_cart.select(sel),
      max_iterations = 100)

    print("Selection: %s " %selection, file = log)
    print("T: %s" %(str(result.T_min)),file = log)
    print("L: %s" %(str(result.L_min)), file = log)
    print("S: %s" %(str(result.S_min)), file = log)
    print("Origin: %s" %(str(cm)), file = log)

    rms = get_tlso_resid(
       result.T_min,result.L_min,result.S_min,cm,u_cart,sites_cart.select(sel))

    print("Final RMS fit of u_cart with TLS: %.2f " %(rms))
    T_list.append(result.T_min)
    L_list.append(result.L_min)
    S_list.append(result.S_min)
    O_list.append(cm)
    tlso_selection_list.append(selection)
    tlso_shift_cart_list.append(
       model_use.shift_cart() if model_use.shift_cart() else (0,0,0) )

  return group_args(
     T_list = T_list,
     L_list = L_list,
     S_list = S_list,
     O_list = O_list,
     tlso_selection_list = tlso_selection_list,
     tlso_shift_cart_list = tlso_shift_cart_list)

def get_tlso_resid(T,L,S,cm,u_cart,xyz):
    '''Get residual between tlso (TLS object) and target'''
    tlso_value = tlso( t = tuple(T), l = tuple(L),
         s = tuple(S), origin = tuple(cm),)

    uanisos= uaniso_from_tls_one_group(
         tlso = tlso_value,
         sites_cart = xyz,
         zeroize_trace=False)
    u_list = flex.double()
    v_list = flex.double()
    i=0
    for u,v in zip(u_cart,uanisos):
      i+=1
      u_list.extend(flex.double(u))
      v_list.extend(flex.double(v))
    diffs = u_list - v_list
    rms = diffs.rms()
    return rms

def create_fine_spacing_array(unit_cell, cell_ratio = 10):
  '''Create an array with spacing of about cell_ratio in P1 with this unit cell
  '''
  new_params= 10*flex.double(unit_cell.parameters()[:3])
  new_params.extend(flex.double(unit_cell.parameters()[3:]))
  new_params=tuple(new_params)
  from cctbx import sgtbx
  xs = crystal.symmetry(
    unit_cell = uctbx.unit_cell(new_params),
    space_group_info=sgtbx.space_group_info(symbol='p1'))
  mi = flex.miller_index()
  return miller.array(miller.set(xs,mi))

def cutoff_values(inside = True):
  '''Return a pre-defined value depending on whether inside is True
     or False'''
  inside_dict = {
      True: group_args(
        cutoff_low = 0.9,
        cutoff_high = None),
      False: group_args(
        cutoff_low = None,
        cutoff_high = 0.1),
      None: group_args(
        cutoff_low = 0.1,
        cutoff_high = 0.9),
      }
  return inside_dict[inside]

def is_inside_mask(mask_map_manager, site_frac = None,
    inside = True):
  '''Return True if site_frac (fractional coords) is inside the mask'''
  if inside not in [True, False, None]:
    return True
  cutoff_low = cutoff_values(inside).cutoff_low
  cutoff_high = cutoff_values(inside).cutoff_high
  xx = mask_map_manager.map_data().tricubic_interpolation(site_frac)
  if ((cutoff_low is None) or (xx >= cutoff_low)) and \
         ((cutoff_high is None) or (xx < cutoff_high)):
    return True
  else:
    return False

def get_weights_for_unit_binning(f_array, i_bin):
    ''' get weighting for each reflection in a scheme for
     interpolating a top-hat function over bins
    '''
    n_used = len(list(f_array.binner().range_used()))
    weights=flex.double(n_used, 0)
    weights[i_bin - 1] = 1
    scale_array=f_array.binner().interpolate( weights, 1) # d_star_power=1
    scale_array.set_selected( (scale_array < 0) , 0) # interpolation sets some neg
    return scale_array


def get_normalization_data_for_unit_binning(f_array):
    ''' Get normalizations for each reflection in a scheme for
     interpolating a top-hat function over bins
    '''
    sum_scale_data = flex.double(f_array.size(),0)
    n_used = len(list(f_array.binner().range_used()))
    for i_bin in f_array.binner().range_used(): # i_bin starts at 1 not 0
      scale_array = get_weights_for_unit_binning(f_array, i_bin)
      sum_scale_data += scale_array
    sum_scale_data.set_selected((sum_scale_data < 1.e-10), 1.e-10)
    return 1./sum_scale_data

def apply_ncs_to_dv_results(
    direction_vectors =None,
    xyz = None,
    scaling_group_info = None,
    ncs_object = None):
  '''Apply NCS (symmetry) to direction-vector results.

    scaling_group_info group_args object:
      direction_vectors: direction vectors dv for anisotropy calculations
      overall_si
      scaling_info_list: si (scaling_info) objects, one for each dv
        each si:  si.target_scale_factors   # scale factors vs sthol2
        si.target_sthol2 # sthol2 values  d = 0.25/sthol2**0.5
                  si.d_min_list
                  si.cc_list
                  si.low_res_cc # low-res average
      ss_b_cart_as_u_cart: anisotropic part of overall correction factor
      overall_scale: radial part of overall correction factor
  '''

  assert ((direction_vectors is None) or
     (list(direction_vectors) == list(scaling_group_info.direction_vectors))
     )

  # work on one location (xyz) with one value ( one scaling_group_info object)

  # Produce a set of xyz and a set of scaling_group_info objects


  # If direction vectors are None then NCS operation just multiplies all the
  #   entries without changing them
  if not direction_vectors or direction_vectors==[None]:
    new_sites = ncs_object.apply_ncs_to_sites(xyz)
    new_scaling_group_info_list = []
    for i in range(ncs_object.max_operators()):
      new_scaling_group_info = deepcopy(scaling_group_info)
      new_scaling_group_info_list.append(new_scaling_group_info)
    return new_sites, new_scaling_group_info_list


  # We want to add on ncs_n new values of xyz, each with n_dv
  #   sets of resolution-bin-values corresponding to n_dv direction vectors

  # The key is, after application of ncs operator j, what is the
  #  order of values

  new_sites = ncs_object.apply_ncs_to_sites(xyz)
  # n_ncs new sites. Now each one should get n_dv sets of values

  # Now question is mapping of which values to which new values
  pointer_to_old_dv_id_dict_list = []
  for dv in scaling_group_info.direction_vectors:
    working_dv_list = ncs_object.apply_ncs_to_sites(dv)
    pointer_to_old_dv_id_dict=get_pointer_to_old_dv_id_dict(
      working_dv_list = working_dv_list,
      dv_list = scaling_group_info.direction_vectors)
    # Now id=pointer_to_old_dv_id_dict[i] says :
    #     values for ncs operator i should come from values[id] for this dv
    pointer_to_old_dv_id_dict_list.append(pointer_to_old_dv_id_dict)

  # Add on ncs_n new values of xyz, each with n_dv
  #   sets of resolution-bin-values corresponding to n_dv direction vectors

  #   scaling_group_info.direction_vectors
  #   scaling_group_info.scaling_info_list: one si entry per direction

  new_scaling_group_info_list = []
  for i in range(ncs_object.max_operators()):
    new_scaling_info_list_by_dv = []
    # i'th ncs operator
    j = 0
    for dv in scaling_group_info.direction_vectors:
      # j'th position in direction vectors
      id = pointer_to_old_dv_id_dict_list[j][i]
      new_scaling_info_list_by_dv.append(
         scaling_group_info.scaling_info_list[id])
      j += 1
    new_scaling_group_info = deepcopy(scaling_group_info)
    new_scaling_group_info.scaling_info_list = new_scaling_info_list_by_dv
    new_scaling_group_info_list.append(new_scaling_group_info)
  # Now new_values is the rearranged version of values appropriate for
  # this xyz this direction_vector and its ncs-related points
  assert len(new_sites) == len(new_scaling_group_info_list)

  return new_sites, new_scaling_group_info_list

def get_pointer_to_old_dv_id_dict(working_dv_list = None, dv_list = None,
   very_similar = 0.95 , allow_multiple_use = True):
  '''
  For each member of working_dv_list, identify best match to member of
  dv_list. Only use each dv_list member once unless allow_multiple_use.
  ID by abs(dot product)
  allow_multiple_use is for matching any to dv_list, False is for
  rearranging only
  '''
  dot_dict={}
  pointer_to_old_dv_id_dict = {}
  n = len(working_dv_list)
  assert allow_multiple_use or (len(dv_list) == n)
  for i in range(n):
    dot_dict[i]={}
    pointer_to_old_dv_id_dict[i] = None
    for j in range(n):
      dot_dict[i][j]=0.
  for i in range(n):
    x,y,z = working_dv_list[i]
    for j in range(dv_list.size()):
      x1,y1,z1 = dv_list[j]
      dot = abs(x*x1+y*y1+z*z1)/((x**2+y**2+z**2)*(x1**2+y1**2+z1**2))**0.5
      dot_dict[i][j] = dot  # dot of working_dv_list[i] to dv_list[j]

  used_list = []
  # See if we can use original positions for any if we are matching 1:1
  if (not allow_multiple_use):
    for i in range(n):
      if dot_dict[i][i] >= very_similar:
        pointer_to_old_dv_id_dict[i] = i
        used_list.append(i)

  # Now work through best to worst
  for i_try in range(n):
    closest_i = None
    closest_j = None
    closest_dot = None
    for i in range(n):
      if pointer_to_old_dv_id_dict[i] is not None: continue
      for j in range(n):
        if (not allow_multiple_use) and  j in used_list: continue
        if not closest_dot or dot_dict[i][j] > closest_dot:
          closest_dot = dot_dict[i][j]
          closest_j = j
          closest_i = i
    if (closest_i is not None) and (closest_j is not None):
      pointer_to_old_dv_id_dict[closest_i] = closest_j
      used_list.append(closest_j)
    else:
      assert allow_multiple_use or (len(used_list) == n)
  return pointer_to_old_dv_id_dict

def get_map_coeffs_as_fp_phi(map_coeffs, d_min= None, n_bins = None):
    '''
    Get map_coeffs as fp and phi. also set up binner if n_bins is not None
    '''
    from cctbx.maptbx.segment_and_split_map import map_coeffs_as_fp_phi
    f_array,phases=map_coeffs_as_fp_phi(map_coeffs)
    if n_bins and not f_array.binner():
      f_array.setup_binner(n_bins=n_bins,d_min=d_min)
    return group_args(
      f_array = f_array,
      phases = phases,
      d_min = d_min)

def create_map_manager_with_value_list(
       n_real = None,
       crystal_symmetry = None,
       value_list = None,
       sites_cart_list = None,
       target_spacing = None,
       max_iterations = None,
       default_value = None):
    '''
      Create a map_manager with values set with a set of sites_cart and values
      Use nearest available value for each grid point, done iteratively
       with radii in shells of target_spacing/2 and up to max_iterations shells
      If default_value is set, use that for all empty locations after
      max_iterations
    '''
    if max_iterations is None:
      if default_value is None:
        max_iterations = 30  # up to 20 grid points away
      else:
        max_iterations = 10

    if default_value is None:
      default_value = 1

    fsc_map = flex.double(flex.grid(n_real),0.)
    fsc_map_manager = MapManager(
       map_data = fsc_map,
       unit_cell_grid = fsc_map.all(),
       unit_cell_crystal_symmetry = crystal_symmetry,
       wrapping = False)
    fsc_set_map_manager = fsc_map_manager.customized_copy(
      map_data = flex.double(flex.grid(n_real),0.))

    sites_frac_list=crystal_symmetry.unit_cell().fractionalize(
       sites_cart_list)
    from cctbx.maptbx import closest_grid_point
    for site_frac,value in zip(sites_frac_list,value_list):
      index = closest_grid_point(
        fsc_map_manager.map_data().accessor(), site_frac)
      fsc_map_manager.map_data()[index] = value
      fsc_set_map_manager.map_data()[index] = 1

    # find anything not set
    not_set = (fsc_map == 0)
    for k in range(max_iterations):
      radius = 0.5 * k * target_spacing
      for i in range(sites_cart_list.size()):
        set_nearby_empty_values(
          fsc_map_manager,
          fsc_set_map_manager,
          sites_cart_list[i:i+1],
          radius,
          value_list[i])
      not_set = (fsc_set_map_manager.map_data() == 0)
      if (not_set.count(True) == 0):
        break
    not_set = (fsc_set_map_manager.map_data() == 0)
    if not_set.count(True) > 0:
      fsc_map_manager.map_data().set_selected(not_set,default_value)
    return fsc_map_manager

def set_nearby_empty_values(
    map_manager,
    set_values_map_manager,
    xyz_list,
    radius,
    value):
  '''
  Set values within radii of xyz_list points to value if not already
      set
  '''
  gias = maptbx.grid_indices_around_sites(
        unit_cell=map_manager.crystal_symmetry().unit_cell(),
        fft_n_real=map_manager.map_data().all(),
        fft_m_real=map_manager.map_data().all(),
        sites_cart=xyz_list,
        site_radii=flex.double(xyz_list.size(),radius))
  for index in gias:
        if set_values_map_manager.map_data()[index] == 0:
          map_manager.map_data()[index] = value
          set_values_map_manager.map_data()[index] = 1

def get_split_maps_and_models(
      map_model_manager = None,
      box_info = None,
      first_to_use = None,
      last_to_use = None):
  '''
  Apply selections and boxing in box_info to generate a set of
  small map_model_managers

  if mask_around_unselected_atoms is set, then mask within each box
     around all the atoms that are not selected (including waters/hetero)
     with a mask_radius of mask_radius and set the value inside the mask to
      masked_value

        mask_around_unselected_atoms = mask_around_unselected_atoms,
        mask_radius = mask_radius,
        masked_value = masked_value,
  '''

  if hasattr(box_info,'tlso_group_info') and box_info.tlso_group_info:
    # cannot pickle tlso values
    box_info.tlso_group_info.tlso_list = None
  box_info = deepcopy(box_info)
  if first_to_use is not None and last_to_use is not None:
    for x in ['lower_bounds_list', 'upper_bounds_list',
       'lower_bounds_with_cushion_list','upper_bounds_with_cushion_list',
       'selection_list']:
      if getattr(box_info,x):  # select those in range
        setattr(box_info,x,getattr(box_info,x)[first_to_use-1:last_to_use])

  mmm_list = []
  if box_info.lower_bounds_with_cushion_list:
    lower_bounds_list = box_info.lower_bounds_with_cushion_list
    upper_bounds_list = box_info.upper_bounds_with_cushion_list
  else:
    lower_bounds_list = box_info.lower_bounds_list
    upper_bounds_list = box_info.upper_bounds_list

  # Make a copy of model that does not have grm and use it for
  #   masking outside model if necessary
  # Select original model atoms to use before doing boxing

  original_model = map_model_manager.model()
  if original_model and box_info.mask_around_unselected_atoms:
    model_no_grm = map_model_manager.model().deep_copy()
    model_no_grm.unset_restraints_manager()
  else:
    model_no_grm = None

  for lower_bounds, upper_bounds, selection in zip(
       lower_bounds_list,
       upper_bounds_list,
       box_info.selection_list,):
    if original_model:
      working_model = original_model.select(selection)
      # Just use the good part here, restore just below
      map_model_manager.set_model(working_model)

      if model_no_grm:
        # Add in the model_no_grm to carry along and shift
        map_model_manager.add_model_by_id(
         model_id='model_no_grm', model = model_no_grm)

    mmm=map_model_manager.extract_all_maps_with_bounds(
     lower_bounds, upper_bounds,
     model_can_be_outside_bounds = True)

    shifted_model_no_grm_not_used = None
    if original_model:
      # Restore original model
      map_model_manager.set_model(original_model)

      if model_no_grm:
        # Collect and remove the shifted model_no_grm
        shifted_model_no_grm = mmm.get_model_by_id(model_id='model_no_grm')
        map_model_manager.remove_model_by_id(model_id='model_no_grm')
        mmm.remove_model_by_id(model_id='model_no_grm')
        shifted_model_no_grm_not_used = shifted_model_no_grm.select(~selection)

    if shifted_model_no_grm_not_used and \
       box_info.mask_around_unselected_atoms:  # mask everything we didn't keep
      # NOTE: only applies mask to map_manager, not any other map_managers
      nnn=mmm.deep_copy()
      nnn.set_model(shifted_model_no_grm_not_used)
      nnn.remove_model_outside_map(boundary=box_info.mask_radius)
      if nnn.model().get_sites_cart().size() > 0: # do something
        nnn.create_mask_around_atoms(
         mask_atoms_atom_radius=box_info.mask_radius,
         mask_id = 'mask')
        mask_mm = nnn.get_map_manager_by_id(map_id = 'mask')
        s = (mask_mm.map_data() > 0.5)
        mmm.map_manager().map_data().set_selected(s,box_info.masked_value)
    elif box_info.mask_all_maps_around_edges:  # mask around edges
      mmm.mask_all_maps_around_edges()
    mmm_list.append(mmm)
  box_info.mmm_list = mmm_list

  return box_info

def get_selections_and_boxes_to_split_model(
        model_id = None,
        map_model_manager = None,
        selection_method = 'by_chain',
        selection_list = None,
        skip_waters = False,
        skip_hetero = False,
        target_for_boxes = 24,
        box_cushion = 3,
        select_final_boxes_based_on_model = None,
        skip_empty_boxes = True,
        mask_around_unselected_atoms = None,
        mask_all_maps_around_edges = None,
        mask_radius = 3,
        masked_value = -10,
        get_unique_set_for_boxes = True,
        mask_id = None,
        exclude_points_outside_density = None,
        minimum_boxes_inside_density = 25,
        log = sys.stdout,
         ):

  '''
    Split up model into pieces using selection_method
    Choices are ['by_chain', 'by_segment','all', 'boxes']
    by_chain:  each chain is a selection
    by_segment:  each unbroken part of a chain is a selection
    boxes:  map is split into target_for_boxes boxes, all atoms in
      each box selected requires map_model_manager to be present
    Skip waters or hetero atoms in selections if specified
    If select_final_boxes_based_on_model and selection_method == 'boxes' then
      make the final boxes just go around the selected parts of the model and
      not tile the map.
    If skip_empty_boxes then skip anything with no model.
    if get_unique_set_for_boxes then get a unique set for 'boxes' method
    If mask_id is set and exclude_points_outside_density , skip boxes
      outside of mask for boxes method
     If exclude_points_outside_density,
       try to add boxes inside density (basically add the proportional
       number of boxes but put them definitely inside the density instead
       of evenly spaced.

  '''

  # Checks
  assert map_model_manager is not None
  selection_method = selection_method.lower()
  assert selection_method in ['supplied_selections',
      'by_chain', 'by_segment','all', 'boxes']
  assert (selection_method != 'boxes') or (
     map_model_manager.map_manager() is not None)

  assert (selection_method != 'supplied_selections') or (
      selection_list is not None)

  # Get selection info for waters and hetero atoms
  info = get_skip_waters_and_hetero_lines(skip_waters, skip_hetero)

  if model_id is None:
    model_id = 'model'

  model = map_model_manager.get_model_by_id(model_id)
  map_manager = map_model_manager.get_any_map_manager()

  # Get the selections
  box_info = group_args(
    selection_as_text_list = [],
    selection_list = [],
    lower_bounds_list = [],
    upper_bounds_list = [],
    lower_bounds_with_cushion_list = [],
    upper_bounds_with_cushion_list = [],
    n_real = map_manager.map_data().all(),
   )

  if selection_list:
    for selection in selection_list:
      if (not skip_empty_boxes) or (selection.count(True) > 0):
        box_info.selection_list.append(selection)

  elif selection_method == 'all':
    selection = model.selection('%s' %(info.no_water_or_het))
    if (not skip_empty_boxes) or (selection.count(True) > 0):
      box_info.selection_list = [selection]
      box_info.selection_as_text_list=[info.no_water_or_het]
  elif selection_method == 'by_chain':
    for chain_id in model.chain_ids(unique_only=True):
      if chain_id.replace(" ",""):
        selection_as_text = " %s (chain %s) " %(
         info.no_water_or_het_with_and,chain_id)
      else:
        selection_as_text = " %s " %(info.no_water_or_het)
      selection = model.selection(selection_as_text)
      if (not skip_empty_boxes) or (selection.count(True) > 0):
        box_info.selection_list.append(selection)
        box_info.selection_as_text_list.append(selection_as_text)
  elif selection_method == 'by_segment':
    selection_strings= get_selections_for_segments(model,
    no_water_or_het_with_and = info.no_water_or_het_with_and)
    for selection_string in selection_strings:
      selection = model.selection(selection_string)
      if (not skip_empty_boxes) or (selection.count(True) > 0):
        box_info.selection_list.append(selection)
        box_info.selection_as_text_list.append(selection_string)

  elif selection_method == 'boxes':
    if info.no_water_or_het and info.no_water_or_het != 'all':
      overall_selection = model.selection("not (%s) " %(info.no_water_or_het))
    else:
      overall_selection = None

    # Select inside boxes without cushion and create cushion too
    if mask_id and exclude_points_outside_density:
      mask_map_manager = map_model_manager.get_map_manager_by_id(mask_id)
    else:
      mask_map_manager = None

    if exclude_points_outside_density and mask_map_manager:
      # Make sure we have some points inside the density
      inside = (mask_map_manager.map_data() > 0.5)
      fraction_inside = inside.count(True)/inside.size()
      target_n = max(minimum_boxes_inside_density,
        int(0.5+ fraction_inside * target_for_boxes))
      volume_inside = fraction_inside * mask_map_manager.crystal_symmetry(
         ).unit_cell().volume()
      dist_min = (volume_inside/target_n)**0.33  # approximate spacing
      target_xyz_center_list = mask_map_manager.trace_atoms_in_map(
         dist_min, target_n)
      if len(target_xyz_center_list) < target_n/2:  # try with smaller grid
        dist_min = dist_min/2.
        target_xyz_center_list = mask_map_manager.trace_atoms_in_map(
           dist_min, target_n)
      exclude_points_outside_density = False  # no longer need it
      print("Using %s points inside density as target_centers " %(
         target_xyz_center_list.size()),file = log)
      dist_min = max(map_model_manager.resolution(),
         dist_min*0.5) # keep those > this dist

    else:
      target_xyz_center_list = None
      dist_min = None

    # Get boxes without and with cushion (cushion may be None)
    box_info = map_manager.get_boxes_to_tile_map(
      target_for_boxes = target_for_boxes,
      do_not_go_over_target = True,
      box_cushion = box_cushion,
      get_unique_set_for_boxes = get_unique_set_for_boxes,
      dist_min = dist_min,
      target_xyz_center_list = target_xyz_center_list,
      )
    print("Ready with %s boxes to check" %(len(box_info.lower_bounds_list)),
      file = log)
    box_info = get_selections_from_boxes(
       box_info = box_info,
       model = model,
       overall_selection = overall_selection,
       skip_empty_boxes = skip_empty_boxes,
       exclude_points_outside_density = exclude_points_outside_density,
       mask_map_manager = mask_map_manager)
    print("Ready with %s ok boxes " %(len(box_info.lower_bounds_list)),
      file = log)

  if (select_final_boxes_based_on_model and model) or (
     not box_info.lower_bounds_list): # get bounds now:
    from cctbx.maptbx.box import get_bounds_around_model
    box_info.lower_bounds_list = []
    box_info.upper_bounds_list = []
    for selection in box_info.selection_list:
      model_use=model.select(selection)
      info = get_bounds_around_model(
        map_manager = map_manager,
        model = model_use,
        box_cushion = box_cushion)
      box_info.lower_bounds_list.append(info.lower_bounds)
      box_info.upper_bounds_list.append(info.upper_bounds)
    box_info.lower_bounds_with_cushion_list = [] # not using these
    box_info.upper_bounds_with_cushion_list = []
  if not box_info.get('selection_as_text_list') or (
       not box_info.selection_as_text_list):
    box_info.selection_as_text_list = [None] * len(box_info.selection_list)
  box_info.exclude_points_outside_density = exclude_points_outside_density
  box_info.mask_around_unselected_atoms = mask_around_unselected_atoms
  box_info.mask_all_maps_around_edges = mask_all_maps_around_edges
  box_info.mask_radius = mask_radius
  box_info.masked_value = masked_value
  return box_info



def get_selections_from_boxes(box_info = None,
    model = None,
    overall_selection = None,
    skip_empty_boxes = None,
    mask_map_manager = None,
    exclude_points_outside_density = None,
   ):
  '''
    Generate a list of selections that covers all the atoms in model,
     grouped by the boxes defined in box_info
  '''
  selection_list = []
  new_lower_bounds_list = []
  new_upper_bounds_list = []
  new_lower_bounds_with_cushion_list = []
  new_upper_bounds_with_cushion_list = []
  for lower_bounds, upper_bounds,lower_bounds_with_cushion, \
    upper_bounds_with_cushion in zip (
      box_info.lower_bounds_list,
      box_info.upper_bounds_list,
      box_info.lower_bounds_with_cushion_list,
      box_info.upper_bounds_with_cushion_list,
       ):
    sel = get_selection_inside_box(
     lower_bounds = lower_bounds,
     upper_bounds = upper_bounds,
     n_real = box_info.n_real,
     model = model,
     crystal_symmetry = box_info.crystal_symmetry)
    if sel and overall_selection:
      sel = (sel & overall_selection)

    # Decide if center of this box is inside mask if supplied
    inside_mask = True
    if mask_map_manager and exclude_points_outside_density:
      center_frac = get_center_of_box_frac(
        lower_bounds = lower_bounds_with_cushion,
        upper_bounds = upper_bounds_with_cushion,
        n_real = box_info.n_real,
        crystal_symmetry = box_info.crystal_symmetry)
      # Value of mask map here
      if not is_inside_mask(mask_map_manager, site_frac = center_frac):
        inside_mask = False

    if inside_mask and (
        (not sel) or (not skip_empty_boxes) or (sel.count(True) > 0)):
      selection_list.append(sel)
      new_lower_bounds_list.append(lower_bounds)
      new_upper_bounds_list.append(upper_bounds)
      new_lower_bounds_with_cushion_list.append(lower_bounds_with_cushion)
      new_upper_bounds_with_cushion_list.append(upper_bounds_with_cushion)
  return group_args(
     ncs_object = box_info.ncs_object,
     n_real = box_info.n_real,
     selection_list = selection_list,
     lower_bounds_list = new_lower_bounds_list,
     upper_bounds_list = new_upper_bounds_list,
     lower_bounds_with_cushion_list = new_lower_bounds_with_cushion_list,
     upper_bounds_with_cushion_list = new_upper_bounds_with_cushion_list,
     )

def get_center_of_box_frac(
     lower_bounds = None,
     upper_bounds = None,
     n_real = None,
     crystal_symmetry = None):
  '''
   Get center of this box in fractional coordinates
  '''

  lower_bounds_frac = tuple([lb / x for lb,x in zip(lower_bounds, n_real)])
  upper_bounds_frac = tuple([ub / x for ub,x in zip(upper_bounds, n_real)])
  average_bounds_frac = tuple ([
     0.5*(lbf+ubf) for lbf,ubf in zip (lower_bounds_frac,upper_bounds_frac)])
  return average_bounds_frac

def get_selection_inside_box(
     lower_bounds = None,
     upper_bounds = None,
     n_real = None,
     model = None,
     crystal_symmetry = None):
  '''
   Get selection for all the atoms inside this box
  '''

  if not model:
    return None
  lower_bounds_frac = tuple([lb / x for lb,x in zip(lower_bounds, n_real)])
  upper_bounds_frac = tuple([ub / x for ub,x in zip(upper_bounds, n_real)])
  sites_frac = model.get_sites_frac()
  lb_a, lb_b, lb_c = lower_bounds_frac
  ub_a, ub_b, ub_c = upper_bounds_frac
  x,y,z = sites_frac.parts()
  s = (
         (x < lb_a) |
         (y < lb_b) |
         (z < lb_c) |
         (x > ub_a) |
         (y > ub_b) |
         (z > ub_c)
         )
  return ~s

def get_skip_waters_and_hetero_lines(skip_waters = True, skip_hetero = True):
  '''Return selection string for skipping waters or hetero atoms'''
  if skip_waters and skip_hetero:
    no_water_or_het = "( (not hetero ) and (not water)) "
  elif skip_waters:
    no_water_or_het = "( (not water)) "
  elif skip_hetero:
    no_water_or_het = "( (not hetero ) ) "
  else:
    no_water_or_het = ""
  if no_water_or_het:
    no_water_or_het_with_and = " %s and " %(no_water_or_het)
  else:
    no_water_or_het_with_and = ""
    no_water_or_het = "all"
  return group_args(
     no_water_or_het = no_water_or_het,
     no_water_or_het_with_and = no_water_or_het_with_and,
     )

def get_selections_for_segments(model,
     skip_waters = True,
     skip_hetero = True,
     no_water_or_het_with_and = None,
     skip_n_residues_on_ends = None,
     minimum_length = 1,
     return_as_group_args = False):
  ''''
    Generate selections corresponding to each segment (chain or part of a chain
    that is separate from remainder of chain)

    Use skip_waters and skip_hetero to specify whether to include them
    If skip_n_residues_on_ends is set, skip residues within
      skip_n_residues_on_ends of an end
  '''

  if not model:
    return []  # nothing to do

  assert isinstance(model, mmtbx.model.manager)

  if no_water_or_het_with_and is None:
    water_het_info = get_skip_waters_and_hetero_lines(
        skip_waters = skip_waters, skip_hetero = skip_hetero)
    no_water_or_het_with_and = water_het_info.no_water_or_het_with_and

  from iotbx.pdb import resseq_encode
  selection_info_list = []
  ph = model.get_hierarchy()
  for m in ph.models()[:1]:
    for chain in m.chains():
      first_resno = None
      last_resno = None
      chain_id = chain.id
      previous_rg = None
      for rg in chain.residue_groups():
        if previous_rg and ( (not rg.link_to_previous) or (not
           residue_group_is_linked_to_previous(rg, previous_rg)) or
         (previous_rg.resseq_as_int() + 1 != rg.resseq_as_int())):
          # break here
          selection_info_list.append( group_args(
            chain_id = chain_id,
            first_resno = first_resno,
            last_resno = last_resno))
          first_resno = None
          last_resno = None
        if not first_resno:
          first_resno = rg.resseq_as_int()
        last_resno = rg.resseq_as_int()
        previous_rg = rg
      if first_resno is not None and last_resno is not None:
        selection_info_list.append( group_args(
            chain_id = chain_id,
            first_resno = first_resno,
            last_resno = last_resno))

  selection_list = []
  for si in selection_info_list:
    if skip_n_residues_on_ends:
      first_resno = si.first_resno + skip_n_residues_on_ends
      last_resno = si.last_resno - skip_n_residues_on_ends
    else:
      first_resno = si.first_resno
      last_resno = si.last_resno
    if (last_resno - first_resno) + 1 < max(1,minimum_length):
      continue
    if return_as_group_args:
      selection_list.append(si)
    else: # usual
      selection_list.append(" %s ( chain %s and resseq %s:%s ) " %(
        no_water_or_het_with_and, si.chain_id,
        resseq_encode(first_resno).strip(), resseq_encode(last_resno).strip()))
  return selection_list

def residue_group_is_linked_to_previous(rg, previous_rg):
  """Return True if this residue group is linked to the previous one"""
  from mmtbx.secondary_structure.find_ss_from_ca import is_close_to
  if is_close_to(rg,previous_rg):
    return True
  elif  rg.resseq_as_int()!=+previous_rg.resseq_as_int()+1:
    return True
  else:
    return False
def get_map_histograms(data, n_slots = 20, data_1 = None, data_2 = None):
  """Create histograms of map values"""
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
  '''Summarize information about map as group_args'''
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

class run_anisotropic_scaling_as_class:
  '''Helper class for anisotropic scaling'''
  def __init__(self, map_model_manager=None,
      direction_vectors = None,
      scale_factor_info= None,
      setup_info = None,
       ):
    self.map_model_manager = map_model_manager
    self.direction_vectors = direction_vectors
    self.scale_factor_info = scale_factor_info
    self.setup_info = setup_info

  def __call__(self,i):
    '''
     Run anisotropic scaling with direction vector i
      To sum up one partial map:
       one bin (sel), one direction vector dv, weights w_dv,
         weights_resolution_bin
       a.calculate value_map map with map_coeffs * w_dv * w_resolution_bin
       b. calculate weight map from position-dependent target_scale_factors
          for dv
       c multiply weight_map * value_map and sum over all bins, dv


    '''
    direction_vector = self.direction_vectors[i]

    # Get the partial map
    scale_factor_info = self.scale_factor_info

    # scale_factor_info.value_list is a set of scaling_group_info objects.
    # scale_factor_info.xyz_list are the coordinates where these apply
    """
    scaling_group_info group_args object:
      direction_vectors: direction vectors dv for anisotropy calculations
      scaling_info_list: si (scaling_info) objects, one for each dv
        each si:  si.target_scale_factors   # scale factors vs sthol2
        si.target_sthol2 # sthol2 values  d = 0.25/sthol2**0.5
                  si.d_min_list
                  si.cc_list
                  si.low_res_cc # low-res average
      ss_b_cart_as_u_cart: anisotropic part of overall correction factor
      overall_scale: radial part of overall correction factor
    """

    #Check for stop_file
    self.map_model_manager.check_stop_file()

    xyz_list = scale_factor_info.xyz_list
    d_min = scale_factor_info.d_min
    smoothing_radius = scale_factor_info.setup_info.smoothing_radius
    n_bins = scale_factor_info.n_bins
    map_id_to_be_scaled = self.setup_info.kw['map_id_to_be_scaled']
    map_model_manager = self.map_model_manager

    # Get Fourier coefficient for map
    map_coeffs = map_model_manager.get_map_manager_by_id(map_id_to_be_scaled
         ).map_as_fourier_coefficients(d_min = d_min)

    new_map_data = flex.double(flex.grid(
        map_model_manager.get_map_manager_by_id(map_id_to_be_scaled
        ).map_data().all()), 0.)

    # Get map for each shell of resolution, weighting by direction vector

    # direction_vector weights:
    f_array_info = get_map_coeffs_as_fp_phi(map_coeffs,
       n_bins = n_bins, d_min = d_min)

    if self.setup_info.kw['aniso_b_cart']:
      # first apply aniso_b_cart to the whole array
      apply_aniso_b_cart_to_f_array_info(f_array_info,
        self.setup_info.kw['b_iso'], d_min, self.setup_info.kw['aniso_b_cart'])

    from cctbx.maptbx.refine_sharpening import get_normalized_weights_para

    # Get weight for each reflection for this direction vector.  Normalize
    #   to make sum of weights over all direction vectors equal to one
    #   for each reflection

    current_weights = get_normalized_weights_para(f_array_info.f_array,
      self.direction_vectors, direction_vector)

    # Map coeffs weighted by alignment with this direction vector
    weighted_map_coeffs = map_coeffs.customized_copy(
      data = map_coeffs.data() * current_weights)

    # Now interpolate scale values over i_bin values
    #  (see _local_sharpen)

    normalization_data = get_normalization_data_for_unit_binning(
      f_array_info.f_array)

    average_scale_factors = get_average_scale_factors(scale_factor_info)
    if not average_scale_factors:
       raise Sorry("No scale factors obtained ... try fewer bins")

    for i_bin in f_array_info.f_array.binner().range_used():
      # Get scale values for i_bin at all points xyz for dv i

      # Get full size map with values corresponding to weighting for resolution
      #   shell i_bin at each point
      scale_value_list,xyz_used_list = \
         map_model_manager._get_scale_values_for_bin(
        xyz_list = xyz_list,
        i_bin = i_bin,
        scale_factor_info = scale_factor_info,
        dv_id = i)
      default_value = average_scale_factors[i_bin-1]
      weight_mm = \
         map_model_manager._create_full_size_map_manager_with_value_list(
        xyz_list = xyz_used_list,
        value_list = scale_value_list,
        smoothing_radius = smoothing_radius,
        default_value = default_value,
        n_boxes = self.setup_info.n_boxes,
        small_n_real = self.setup_info.small_n_real,
      )

      # Get weights on each Fourier coeff, emphasizing this bin
      weights_top_hat_shell = get_weights_for_unit_binning(
         f_array_info.f_array, i_bin) * normalization_data

      shell_map_coeffs = weighted_map_coeffs.customized_copy(
        data = weighted_map_coeffs.data() * weights_top_hat_shell)

      shell_mm = map_model_manager.map_manager(
         ).fourier_coefficients_as_map_manager(shell_map_coeffs)

      # Get value of this map at xyz:
      new_map_data += weight_mm.map_data() * shell_mm.map_data()
    mm = map_model_manager.get_map_manager_by_id(map_id_to_be_scaled
      ).customized_copy(map_data = new_map_data)

    file_name = os.path.join(
        self.setup_info.temp_dir,'partial_map_%s.ccp4' %(i))
    from iotbx.data_manager import DataManager
    dm = DataManager()
    dm.set_overwrite(True)
    dm.write_real_map_file(mm, file_name)
    result = group_args(
      file_name = file_name,
    )

    return result

class run_fsc_as_class:
  '''Helper class to run FSC calculation'''
  def __init__(self, map_model_manager=None, run_list=None,
      box_info = None):
    self.map_model_manager = map_model_manager
    self.run_list = run_list
    self.box_info = box_info

  def __call__(self,i):
    '''
     Run a group of fsc calculations with kw
     specifying which to run

    '''
    # We are going to run with the i'th set of keywords
    kw=self.run_list[i]

    # Get the method name and expected_result_names and remove them from kw
    first_to_use = kw['first_to_use']
    last_to_use = kw['last_to_use']

    xyz_list = flex.vec3_double()
    value_list = []
    # offset to map absolute on to self.map_model_manager
    offset = self.map_model_manager.get_map_manager_by_id(self.box_info.map_id
      ).shift_cart()

    expected_rms_fc_list = None

    for i in range(first_to_use, last_to_use + 1):
      self.map_model_manager.check_stop_file()
      new_box_info = get_split_maps_and_models(
        map_model_manager = self.map_model_manager,
        box_info = self.box_info,
        first_to_use = i,
        last_to_use = i)
      mmm = new_box_info.mmm_list[0]
      if mmm.verbose:
        mmm.set_log(self.map_model_manager.log)

      xyz = mmm.get_map_manager_by_id(self.box_info.map_id
         ).absolute_center_cart()

      mmm.mask_all_maps_around_edges(soft_mask_radius=self.box_info.resolution)
      # Two choices for methods to get fsc:  _get_weights_in_shells or
      #   _map_map_fsc.   The weights_in_shells method is designed for scaling
      #  and map_map_fsc is designed to get local resolution.

      if self.box_info.return_scale_factors:
        # Get scaling weights
        map_coeffs = self.map_model_manager.get_map_manager_by_id(self.box_info.map_id
         ).map_as_fourier_coefficients(d_min=self.box_info.minimum_resolution)

        scaling_group_info = mmm._get_weights_in_shells(
           map_id = self.box_info.map_id,
           map_id_1 = self.box_info.map_id_1,
           map_id_2 = self.box_info.map_id_2,
           n_bins=self.box_info.n_bins,
           is_model_based=self.box_info.is_model_based,
           optimize_b_eff=self.box_info.optimize_b_eff,
           equalize_power=self.box_info.equalize_power,
           rmsd=self.box_info.rmsd,
           is_external_based=self.box_info.is_external_based,
           resolution = self.box_info.resolution, # nominal resolution
           d_min = self.box_info.minimum_resolution,
           direction_vectors = self.box_info.direction_vectors,
           minimum_low_res_cc = self.box_info.minimum_low_res_cc,
           get_scale_as_aniso_u = self.box_info.get_scale_as_aniso_u,
           use_dv_weighting = self.box_info.use_dv_weighting,
           n_direction_vectors = self.box_info.n_direction_vectors,
           run_analyze_anisotropy = self.box_info.run_analyze_anisotropy,
           expected_rms_fc_list = self.box_info.expected_rms_fc_list,
           expected_ssqr_list = self.box_info.expected_ssqr_list,
           expected_ssqr_list_rms = self.box_info.expected_ssqr_list_rms,
           tlso_group_info = self.box_info.tlso_group_info,
           model_id_for_rms_fc = self.box_info.model_id_for_rms_fc,
           replace_aniso_with_tls_equiv =
                self.box_info.replace_aniso_with_tls_equiv)
        if scaling_group_info:
          """
    scaling_group_info group_args object:
      direction_vectors: direction vectors dv for anisotropy calculations
      scaling_info_list: si (scaling_info) objects, one for each dv
        each si:  si.target_scale_factors   # scale factors vs sthol2
        si.target_sthol2 # sthol2 values  d = 0.25/sthol2**0.5
                  si.d_min_list
                  si.cc_list
                  si.low_res_cc # low-res average
      ss_b_cart_as_u_cart: anisotropic part of overall correction factor
      overall_scale: radial part of overall correction factor
          """

          xyz_list.append(tuple(col(xyz)+col(offset) ))
          value_list.append(scaling_group_info)

      else: # Get local resolution
        d_min = mmm.map_map_fsc(fsc_cutoff = self.box_info.fsc_cutoff,
          map_id_1 = self.box_info.map_id_1,
          map_id_2 = self.box_info.map_id_2,
          n_bins=self.box_info.n_bins).d_min
        if d_min:
          d_min = max(d_min, self.box_info.minimum_resolution)
          xyz_list.append(tuple(col(xyz)+col(offset) ))
          value_list.append(d_min)
    result = group_args(
      n_bins = self.box_info.n_bins,
      d_min = self.box_info.minimum_resolution,
      xyz_list=xyz_list,
      value_list = value_list)
    return result
