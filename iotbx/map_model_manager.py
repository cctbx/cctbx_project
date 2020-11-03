from __future__ import absolute_import, division, print_function
import sys, os
from libtbx.utils import Sorry
from cctbx import maptbx
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

# Reserved phil scope for MapModelManager
map_model_phil_str = '''
map_model {
  full_map = None
    .type = path
    .help = Input full map file
    .short_caption = Full map
    .style = file_type:ccp4_map input_file

  half_map = None
    .type = path
    .multiple = True
    .help = Input half map files
    .short_caption = Half maps
    .style = file_type:ccp4_map input_file
  model = None
    .type = path
    .help = Input model file
    .style = file_type:pdb input_file
    .short_caption = Model
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
               log              = None,
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
    self.set_name(name)
    self.set_verbose(verbose)

    # Initialize
    self._nproc = None
    self._multiprocessing = None
    self._queue_run_command = None
    self._resolution = None
    self._map_dict={}
    self._model_dict = {}
    self._force_wrapping = wrapping
    self._warning_message = None




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
      if any_map_manager and m and (not ignore_symmetry_conflicts):
        if not any_map_manager.is_similar(m,
           absolute_angle_tolerance = absolute_angle_tolerance,
           absolute_length_tolerance = absolute_length_tolerance,
          ):
          raise Sorry("Map manager '%s' is not similar to '%s': %s" %(
           m.file_name,any_map_manager.file_name,
            m.warning_message())+
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
    pickle_dict = self.__dict__.copy()
    if isinstance(self.log, io.TextIOWrapper):
      pickle_dict['log'] = None
    return pickle_dict

  def __setstate__(self, pickle_dict):
    self.__dict__ = pickle_dict
    if self.log is None:
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


  def set_verbose(self, verbose = None):
    '''
       Set verbose
    '''
    self.verbose = verbose


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

  def resolution(self,
    use_fsc_if_no_resolution_available_and_maps_available = True,
    map_id_1 = 'map_manager_1',
    map_id_2 = 'map_manager_2',
     ):
    if self._resolution: # have it already
      return self._resolution

    else:  # figure out resolution
      resolution = None
      if use_fsc_if_no_resolution_available_and_maps_available and \
          self.get_map_manager_by_id(map_id_1) and \
          self.get_map_manager_by_id(map_id_2):
        resolution = self.map_map_fsc(  # get resolution from FSC
          map_id_1 = map_id_1,
          map_id_2 = map_id_2).d_min
        if resolution:
          print("\nResolution estimated from FSC of '%s' and '%s: %.3f A " %(
           map_id_1, map_id_2, resolution), file = self.log)

      if (not resolution) and self.map_manager():
        # get resolution from map_manager
        resolution = self.map_manager().resolution()
        print("\nMinimum resolution obtained from map_manager: %.3f A " %(
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

    if multiprocessing:
      assert multiprocessing in ['multiprocessing','sge','lsf','pbs',
         'condor','pbspro','slurm']
      self._multiprocessing = multiprocessing

    if queue_run_command:
      self._queue_run_command = queue_run_command

  def set_resolution(self, resolution):
    ''' Set nominal resolution '''
    self._resolution = resolution

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
    keys = self.map_dict().keys()
    if not keys:
      return
    else:
      return self.map_dict()[keys[0]]

  def get_map_data_by_id(self, map_id):
    ''' Get map_data from a map_manager with the name map_id'''
    map_manager = self.get_map_manager_by_id(map_id)
    if map_manager:
      return map_manager.map_data()
    else:
      return None

  def set_model(self,model):
    '''
     Overwrites existing model with id 'model'
    '''
    self.add_model_by_id(model,'model')


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

  def set_map_manager(self, map_manager):
    '''
     Overwrites existing map_manager with id 'map_manager'
    '''
    self.add_map_manager_by_id(map_manager, 'map_manager')

  def add_map_manager_by_id(self, map_manager, map_id,
     overwrite = True):
    '''
     Add a new map_manager
     Must be similar to existing
     Overwrites any existing with the same id unless overwrite = False
     Is a mask if is_mask is set
    '''
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
    if not self._map_dict.get(map_id):
      self._print ("No map to write out with id='%s'" %(map_id))
    elif not file_name:
      self._print ("Need file name to write map")
    else:
      self._map_dict.get(map_id).write_map(file_name = file_name)

  def write_model(self,
     file_name,
     model_id = None,
     ):
    if not model_id:
       model_id = 'model'
    if not self.get_model_by_id(model_id = model_id):
      self._print ("No model to write out")
    elif not file_name:
      self._print ("Need file name to write model")
    else:
      # Write out model

      f = open(file_name, 'w')
      print(self.get_model_by_id(model_id = model_id).model_as_pdb(), file = f)
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
     model_can_be_outside_bounds = None):
    '''
      Runs box_all_maps_with_bounds_and_shift_origin with extract_box=True
    '''
    return self.box_all_maps_with_bounds_and_shift_origin(
      lower_bounds = lower_bounds,
      upper_bounds = upper_bounds,
      model_can_be_outside_bounds = model_can_be_outside_bounds,
      extract_box = True)

  def box_all_maps_with_bounds_and_shift_origin(self,
     lower_bounds,
     upper_bounds,
     model_can_be_outside_bounds = None,
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
    model = self._model_dict.get(model_info.model_id,None)

    if extract_box and model: # make sure everything is deep_copy
      model = model.deep_copy()

    # Make box with bounds and apply it to model, first map
    box = with_bounds(
      map_manager = self._map_dict[map_info.map_id],
      lower_bounds = lower_bounds,
      upper_bounds = upper_bounds,
      model = model,
      wrapping = self._force_wrapping,
      model_can_be_outside_bounds = model_can_be_outside_bounds,
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
     selection = None,
     select_unique_by_ncs = False,
     model_can_be_outside_bounds = None,
     box_cushion = 5.):
    '''
      Runs box_all_maps_around_model_and_shift_origin with extract_box=True
      Use either selection_string or selection if present
    '''
    return self.box_all_maps_around_model_and_shift_origin(
      selection_string = selection_string,
      selection = selection,
      box_cushion = box_cushion,
      select_unique_by_ncs = select_unique_by_ncs,
      model_can_be_outside_bounds = model_can_be_outside_bounds,
      extract_box = True)

  def box_all_maps_around_model_and_shift_origin(self,
     selection_string = None,
     selection = None,
     box_cushion = 5.,
     select_unique_by_ncs = False,
     model_can_be_outside_bounds = None,
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
    elif extract_box: # make sure everything is deep_copy
      model = model.deep_copy()

    # Make box around model and apply it to model, first map
    # This step modifies model in place and creates a new map_manager
    box = around_model(
      map_manager = self._map_dict[map_info.map_id],
      model = model,
      box_cushion = box_cushion,
      wrapping = self._force_wrapping,
      model_can_be_outside_bounds = model_can_be_outside_bounds,
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
     model_can_be_outside_bounds = None,
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
     extract_box = True)

  def box_all_maps_around_density_and_shift_origin(self,
     box_cushion = 5.,
     threshold = 0.05,
     map_id = 'map_manager',
     get_half_height_width = True,
     model_can_be_outside_bounds = None,
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
      model_can_be_outside_bounds = model_can_be_outside_bounds,
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
     model_can_be_outside_bounds = None,
     mask_id = 'mask'):
    '''
      Runs box_all_maps_around_mask_and_shift_origin with extract_box=True
    '''
    return self.box_all_maps_around_mask_and_shift_origin(
     box_cushion = 5.,
     mask_id = mask_id,
      model_can_be_outside_bounds = model_can_be_outside_bounds,
     extract_box = True)

  def box_all_maps_around_mask_and_shift_origin(self,
     box_cushion = 5.,
     mask_id = 'mask',
     model_can_be_outside_bounds = None,
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
      model_can_be_outside_bounds = model_can_be_outside_bounds,
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
    assert isinstance(map_manager, MapManager)
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

  def merge_split_maps_and_models(self,
      box_info = None):
    '''
      Replaces coordinates in working model with those from the
        map_model_managers in box_info.  The box_info object should
        come from running split_up_map_and_model in this instance
        of the map_model_manager.
    '''

    print("\nMerging coordinates from %s boxed models into working model" %(
      len(box_info.selection_list)), file = self.log)

    i = 0
    for selection, mmm in zip (box_info.selection_list, box_info.mmm_list):
      i += 1
      model_to_merge = self.get_model_from_other(mmm)
      sites_cart = self.model().get_sites_cart()
      new_coords=model_to_merge.get_sites_cart()
      original_coords=sites_cart.select(selection)
      rmsd=new_coords.rms_difference(original_coords)
      print("RMSD for %s coordinates in model %s: %.3f A" %(
         original_coords.size(), i, rmsd), file = self.log)
      sites_cart.set_selected(selection, new_coords)
      self.model().set_crystal_symmetry_and_sites_cart(sites_cart = sites_cart,
        crystal_symmetry = self.model().crystal_symmetry())

  def split_up_map_and_model_by_chain(self,
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
      skip_waters = skip_waters,
      skip_hetero = skip_hetero,
      box_cushion = box_cushion,
      mask_around_unselected_atoms = mask_around_unselected_atoms,
      mask_radius = mask_radius,
      masked_value = masked_value,
      apply_box_info = apply_box_info,
      write_files = write_files)

  def split_up_map_and_model_by_segment(self,
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
      skip_waters = skip_waters,
      skip_hetero = skip_hetero,
      box_cushion = box_cushion,
      mask_around_unselected_atoms = mask_around_unselected_atoms,
      mask_radius = mask_radius,
      masked_value = masked_value,
      apply_box_info = apply_box_info,
      write_files = write_files)

  def split_up_map_and_model_by_supplied_selections(self,
    selection_list,
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
      selection_list = selection_list,
      box_cushion = box_cushion,
      mask_around_unselected_atoms = mask_around_unselected_atoms,
      mask_radius = mask_radius,
      masked_value = masked_value,
      apply_box_info = apply_box_info,
      write_files = write_files)

  def split_up_map_and_model_by_boxes(self,
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
     ):
    '''
     Split up the map, creating boxes that time the entire map.

     Try to get about target_for_boxes boxes.

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
      target_for_boxes = target_for_boxes,
      select_final_boxes_based_on_model = select_final_boxes_based_on_model,
      skip_empty_boxes = skip_empty_boxes,
      skip_waters = skip_waters,
      box_cushion = box_cushion,
      mask_around_unselected_atoms = mask_around_unselected_atoms,
      mask_radius = mask_radius,
      masked_value = masked_value,
      apply_box_info = apply_box_info,
      write_files = write_files)

  def _split_up_map_and_model(self,
    selection_method = 'by_chain',
    selection_list = None,
    skip_waters = False,
    skip_hetero = False,
    target_for_boxes = 24,
    select_final_boxes_based_on_model = True,
    skip_empty_boxes = True,
    mask_around_unselected_atoms = None,
    mask_radius = 3,
    masked_value = -10,
    write_files = False,
    box_cushion = 3,
    apply_box_info = True,
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

    if mask_around_unselected_atoms is set, then mask within each box
     around all the atoms that are not selected (including waters/hetero)
     with a mask_radius of mask_radius and set the value inside the mask to
      masked_value

    '''
    print ("Splitting up map and model into overlapping boxes (%s method)" %(
       selection_method), file = self.log)

    # Get selections and boxes
    box_info = get_selections_and_boxes_to_split_model(
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
        mask_radius = mask_radius,
        masked_value = masked_value,
      )
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
        model_file = "model_%s.pdb" %(i)
        map_file = "map_%s.ccp4" %(i)
        dm.write_model_file(mmm.model(), model_file)
        dm.write_real_map_file(mmm.map_manager(), map_file)
    return box_info

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
        soft_mask_radius = self.resolution()
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
      soft_mask_radius = self.resolution()

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
        soft_mask_radius = self.resolution()
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
          soft_mask_radius = self.resolution()
      cm.soft_mask(soft_mask_radius = soft_mask_radius)

    # Put the mask in map_dict id'ed with mask_id
    self.add_map_manager_by_id(map_manager = cm.map_manager(),
      map_id = mask_id)

  def expand_mask(self,
     buffer_radius = 5,
     resolution = None,
     soft_mask = True,
     soft_mask_radius = None,
     mask_id = 'mask',
      ):
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
     than boundary)
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

  def get_rms_f_list(self,map_id = 'map_manager',
      d_min = None,  # minimum resolution for calculations
      n_bins = None,
      resolution = None,  # nominal resolution
      ):
    ''' Return list of rms amplitude by bins '''
    assert d_min and n_bins
    from cctbx.maptbx.segment_and_split_map import map_coeffs_to_fp
    map_coeffs=self.get_map_manager_by_id(
      map_id).map_as_fourier_coefficients( d_min = d_min)
    f_array = get_map_coeffs_as_fp_phi(map_coeffs, n_bins = n_bins,
          d_min = d_min).f_array
    rms_f_list=flex.double()
    sthol2_list=flex.double()
    dsd = f_array.d_spacings().data()

    n_bins_use = min(n_bins,max(3,n_bins//3))
    for i_bin in f_array.binner().range_used():
      sel       = f_array.binner().selection(i_bin)
      f        = map_coeffs.select(sel)
      f_array_f=map_coeffs_to_fp(f)
      rms_f=f_array_f.data().norm()
      rms_f_list.append(rms_f)
      d         = dsd.select(sel)
      d_avg     = flex.mean(d)
      sthol2=0.25/d_avg**2
      sthol2_list.append(sthol2)
      if i_bin-1 > n_bins_use and (
          (not resolution) or (d_avg >= resolution)):
         n_bins_use = i_bin - 1

    return group_args(
       rms_f_list = rms_f_list,
       sthol2_list = sthol2_list,
       n_bins_use = n_bins_use)


  def _update_kw_with_map_info(self, local_kw, previous_kw = None,
      text = 'overall', have_previous_scaled_data = None,
      map_id_scaled_list = None):


    if have_previous_scaled_data:
      if previous_kw.get('map_id_1'):
        local_kw['map_id_1'] = previous_kw['map_id_scaled_list'][
         previous_kw['map_id_to_be_scaled_list'].index(previous_kw['map_id_1'])]
        print("Map 1 to use in determining scaling: '%s' " %(
          local_kw['map_id_1']), file = self.log)
      if previous_kw.get('map_id_2'):
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

    if self.model():
      cc = self.map_model_cc(map_id=local_kw['map_id'])
      print ("Current map-model CC for '%s': %.3f " %(local_kw['map_id'],cc),
         file = self.log)

    return local_kw


  def _sharpen_overall_local_overall(self, kw, method):

      assert kw.get('map_id_to_be_scaled_list') is None or (
        kw['map_id']  in kw['map_id_to_be_scaled_list']) # map_id_to_be_scaled not ok

      # Set up list of maps to be scaled
      kw = self.set_map_id_lists(kw)
      kw['overall_sharpen_before_and_after_local'] = False

      # run sharpening without local sharpening first
      #  Then set maps to scale as the scaled maps from this run

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

      if self.model():
        cc = self.map_model_cc(map_id = scaled_map_id)
        print ("Current map-model CC for '%s': %.3f " %(scaled_map_id,cc),
           file = self.log)

      print ("\n",79*"=","\nDone with local and overall sharpening\n",79*"=",
           file = self.log)

  def set_map_id_lists(self,kw):

    if kw['map_id_to_be_scaled_list'] is None:
      kw['map_id_to_be_scaled_list'] = [kw['map_id']]
      if kw.get('map_id_1') and kw.get('map_id_2'): # half-map sharpening
         kw['map_id_to_be_scaled_list'].append(kw['map_id_1'])
         kw['map_id_to_be_scaled_list'].append(kw['map_id_2'])
    if kw['map_id_scaled_list'] is None:
      kw['map_id_scaled_list'] = []
      for id in kw['map_id_to_be_scaled_list']:
         kw['map_id_scaled_list'].append("%s_scaled" %(id))
    return kw

  def external_sharpen(self,
      map_id = 'map_manager',
      map_id_external_map = 'external_map',
      map_id_to_be_scaled_list = None,
      map_id_scaled_list = None,
      resolution = None,
      n_bins = None,
      n_boxes = None,
      core_box_size = None,
      box_cushion = None,
      smoothing_radius = None,
      local_sharpen = None,
      anisotropic_sharpen = None,
      overall_sharpen_before_and_after_local = True,
      get_scale_as_aniso_u = None,
      use_dv_weighting = None,
      n_direction_vectors = None,
      run_analyze_anisotropy = None,
      nproc = None,
    ):
    '''
     Scale map_id with scale factors identified from map_id vs
      map_id_external_map
     Changes the working map_manager

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

    print ("Running external map sharpening ", file = self.log)

    kw['map_id_2'] = map_id_external_map
    kw['is_external_based'] = True
    del kw['map_id_external_map']

    self._sharpen_map( **kw)

  def half_map_sharpen(self,
      map_id = 'map_manager',
      map_id_1 = 'map_manager_1',
      map_id_2 = 'map_manager_2',
      map_id_scaled_list = None,
      map_id_to_be_scaled_list = None,
      resolution = None,
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
      run_analyze_anisotropy = None,
      spectral_scaling = True,
      expected_rms_fc_list = None,
      model_for_rms_fc = None,
      replace_aniso_with_tls_equiv = None,
      max_abs_b = None,
      nproc = None,
      optimize_b_eff = None,
      equalize_power = None,
      overall_sharpen_before_and_after_local = True,
    ):
    '''
     Scale map_id with scale factors identified from map_id_1 vs map_id_2
     Changes the working map_manager unless map_id_scaled_list is set.

     max_abs_b applies if get_scale_as_aniso_u and anisotropic_sharpen and
        local_sharpen are set. It limits range of anisotropic B.  Default is
        100 at 4 A, proportional to resolution squared

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

    # Now get scaling from comparison of the two half-maps
    #  apply the scaling to map_id_to_be_scaled

    self._sharpen_map(**kw)

  def model_sharpen(self,
      map_id = 'map_manager',
      model_id = 'model',
      map_id_scaled_list = None,
      map_id_to_be_scaled_list = None,
      resolution = None,
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
      run_analyze_anisotropy = None,
      spectral_scaling = True,
      expected_rms_fc_list = None,
      model_for_rms_fc = None,
      replace_aniso_with_tls_equiv = None,
      max_abs_b = None,
      nproc = None,
      optimize_b_eff = None,
      equalize_power = None,
      map_id_model_map = 'model_map_for_scaling',
      optimize_with_model = None,
      overall_sharpen_before_and_after_local = True,
      mask_around_model = True,
    ):
    '''
     Scale map_id with scale factors identified from map_id vs model
     Changes the working map_manager unless map_id_scaled is set.

     max_abs_b applies if get_scale_as_aniso_u and anisotropic_sharpen and
        local_sharpen are set. It limits range of anisotropic B.  Default is
        100 at 4 A, proportional to resolution squared

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

    if model_for_rms_fc is None:
      kw['model_for_rms_fc'] = self.get_model_by_id(kw['model_id'])

    # Allow sharpening globally before and after local sharpening
    if local_sharpen and overall_sharpen_before_and_after_local:
      return self._sharpen_overall_local_overall(kw = kw,
        method = self.model_sharpen)

    del kw['mask_around_model']  # REQUIRED
    del kw['model_id']  # REQUIRED
    del kw['map_id_model_map']  # REQUIRED
    del kw['optimize_with_model']  # REQUIRED
    del kw['overall_sharpen_before_and_after_local']  # REQUIRED


    # Make a copy of this map_model manager so we can modify it without
    #  changing the original
    working_mmm = self.deep_copy()
    working_mmm.set_name('working_mmm')

    # Working resolution is resolution * d_min_ratio
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

    working_mmm.generate_map(model=model,
       gridding=working_mmm.get_any_map_manager().map_data().all(),
       d_min=d_min,
       map_id = map_id_model_map)

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

    # Now get scaling from comparison of working_map_id_to_be_scaled and
    #   map_id_model_map, and
    #  apply the scaling to  unmasked_map_id_list
    working_mmm._sharpen_map(**kw)

    # And set this map_manager
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
      equalize_power = None,
      resolution = None,
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
      minimum_low_res_cc = None,
      get_scale_as_aniso_u = None,
      use_dv_weighting = None,
      n_direction_vectors = None,
      run_analyze_anisotropy = None,
      spectral_scaling = None,
      expected_rms_fc_list = None,
      model_for_rms_fc = None,
      replace_aniso_with_tls_equiv = None,
      max_abs_b = None,
      overall_sharpen_before_and_after_local = None, # ignored
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


    # Checks
    assert self.get_map_manager_by_id(map_id)
    assert (
    (self.get_map_manager_by_id(map_id_1) or
        is_model_based or is_external_based) and
       self.get_map_manager_by_id(map_id_2))

    # remove any extra models and maps to speed up boxing and not modify orig
    working_mmm = self._get_map_model_manager_with_selected(
      map_id_list=[map_id,map_id_1,map_id_2],deep_copy=True)
    working_mmm.set_log(self.log)

    # Note the map that is going to be scaled
    assert kw['map_id_to_be_scaled_list']
    assert kw['map_id_scaled_list']
    assert len(kw['map_id_to_be_scaled_list']) == len(kw['map_id_scaled_list'])


    print("\nMaps in '%s' to be scaled will come from '%s'" %(
       self.name,str(kw['map_id_to_be_scaled_list'])),
       file = self.log)
    print("Scaled maps in '%s' will go in '%s'" %(
       self.name,str(kw['map_id_scaled_list'])), file = self.log)

    if local_sharpen:
      # Now local sharpening
      print("\nSetting up local sharpening ...\n",file = self.log)
      # We need the map to be scaled here too. Deep copy it so original is not
      #  modified
      for id in kw['map_id_to_be_scaled_list']:
        working_mmm.add_map_manager_by_id(map_id=id,
          map_manager=self.get_map_manager_by_id(
            map_id=id).deep_copy())

      if n_bins is None: # set it here
        kw['n_bins'] = n_bins_default_local

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
        else:
          print("No local-sharpened map obtained "+
            "in '%s' from map_model_manager '%s'" %(
            id,working_mmm.name), file = self.log)

      return

    # Here to run overall
    print ("\nRunning overall sharpening ", file = self.log)
    if n_bins is None:
      n_bins = n_bins_default

    # Get basic info including minimum_resolution (cutoff for map_coeffs)
    setup_info = working_mmm._get_box_setup_info(map_id_1, map_id_2,
      resolution,
      smoothing_radius = smoothing_radius,
      skip_boxes = True)
    if spectral_scaling and (not expected_rms_fc_list):
        from cctbx.development.approx_amplitude_vs_resolution import \
          approx_amplitude_vs_resolution
        aavr = approx_amplitude_vs_resolution(
           d_min = setup_info.minimum_resolution,
           resolution = setup_info.resolution,
           n_bins = n_bins,
           map_model_manager = working_mmm,
           model = model_for_rms_fc,
           out = self.log)
        expected_rms_fc_list = aavr.get_target_scale_factors()
        n_bins_use = aavr.n_bins_use # number of bins to resolution

    resolution = setup_info.resolution
    working_mmm.set_resolution(resolution)
    d_min = setup_info.minimum_resolution
    print ("Nominal resolution of map: %.2f A  Minimum resolution: %.2f A" %(
        resolution,d_min),
      file = self.log)
    working_mmm.mask_all_maps_around_edges(soft_mask_radius=resolution)

    map_coeffs = working_mmm.get_map_manager_by_id(
       map_id).map_as_fourier_coefficients( d_min = d_min)

    n_bins =working_mmm._set_n_bins(n_bins = n_bins,
      d_min = d_min, map_coeffs = map_coeffs,
      local_sharpen = False)

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
        model_for_rms_fc = model_for_rms_fc,
        replace_aniso_with_tls_equiv = replace_aniso_with_tls_equiv,
        )
    # scaling_group_info group_args object has direction vectors, list of si:
    #  scaling_group_info.direction_vectors
    #  scaling_group_info.overall_u_cart_to_remove  aniso scale factor (u_cart)
    #  scaling_group_info.scaling_info_list: one si entry per direction vector
    #    si.target_scale_factors
    #    si.target_sthol2
    #    si.d_min_list
    #    si.cc_list
    #    si.low_res_cc # low-res average


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
      if get_scale_as_aniso_u and \
           hasattr(scaling_group_info,'overall_u_cart_to_remove') and \
          scaling_group_info.overall_u_cart_to_remove:
        print ("\nApplying scale factors by removing anisotropy ",
         file = self.log)

        new_map_manager = working_mmm.remove_anisotropy(
          map_coeffs = map_coeffs_to_be_scaled,
          d_min = None,
         overall_u_cart_to_remove = scaling_group_info.overall_u_cart_to_remove)

      else: # usual
        new_map_manager = working_mmm._apply_scale_factors_in_shells(
          map_coeffs_to_be_scaled,
          n_bins,
          d_min,
          scaling_group_info = scaling_group_info,
          direction_vectors = direction_vectors)
      if not new_map_manager:
        print("Not applying scaling to '%s'" %(id) ,file = self.log)
      else: # usual
        # All done... Set map_manager now
        print ("Scaled map is '%s' in %s" %(new_id,self.name), file = self.log)
        new_map_manager.file_name = self.get_map_manager_by_id(id).file_name
        self.add_map_manager_by_id(map_manager = new_map_manager,
          map_id = new_id)

  def remove_anisotropy(self,
        d_min = None,
        map_coeffs = None,
        u_cart_to_remove = None,
        u_cart_to_add = None,
        overall_u_cart_to_remove = None,
        map_id = 'map_manager'):
   '''
   Remove anisotropy from map, optionally remove anisotropy specified by
    overall_u_cart_to_remove
   '''
   assert map_coeffs or d_min

   if not map_coeffs:
      assert self.get_map_manager_by_id(map_id)
      map_coeffs = self.get_map_manager_by_id(map_id
        ).map_as_fourier_coefficients(d_min = d_min)
   if not d_min:
     d_min = map_coeffs.d_min()

   from cctbx.maptbx.refine_sharpening import analyze_aniso_object
   f_array_info = get_map_coeffs_as_fp_phi(map_coeffs, n_bins = 1,
          d_min = self.resolution())
   f_array = f_array_info.f_array
   if (not overall_u_cart_to_remove) and (not u_cart_to_remove):
     # Get anisotropy if not supplied and apply it
     analyze_aniso = analyze_aniso_object()
     analyze_aniso.set_up_aniso_correction(f_array=f_array)
     print(" Finding and removing anisotropy (positive is sharpening) "+
        "b_cart = %s " %( str(analyze_aniso.b_cart)), file = self.log)
     scaled_f_array = analyze_aniso.apply_aniso_correction(f_array=f_array)
   else:  # apply overall_u_cart_to_remove
     from mmtbx.scaling import absolute_scaling
     if overall_u_cart_to_remove:
       print(" Removing anisotropy (positive is sharpening).b_cart = %s " %(
         str(adptbx.u_as_b(tuple(col(overall_u_cart_to_remove))))),
         file = self.log)
       u_star= adptbx.u_cart_as_u_star(
         f_array.unit_cell(), overall_u_cart_to_remove)
       scaled_f_array = absolute_scaling.anisotropic_correction(
        f_array,0.0, u_star ,must_be_greater_than=-0.0001)
     else: # remove and add
       u_star= adptbx.u_cart_as_u_star(
         f_array.unit_cell(), u_cart_to_remove)
       u_star_aniso_removed = adptbx.u_cart_as_u_star(
         f_array.unit_cell(), tuple(-col(u_cart_to_add)))

       scaled_f_array = absolute_scaling.anisotropic_correction(
         f_array,0.0, u_star ,must_be_greater_than=-0.0001)

       scaled_f_array = absolute_scaling.anisotropic_correction(
         scaled_f_array,0.0,u_star_aniso_removed,must_be_greater_than=-0.0001)

     scaled_f_array.set_observation_type( f_array)

   return self.map_manager(
       ).fourier_coefficients_as_map_manager(
         scaled_f_array.phase_transfer(phase_source=f_array_info.phases,
         deg=True))

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

    other._resolution = self._resolution
    other.set_log(self.log)
    other._nproc = self._nproc
    other._multiprocessing = self._multiprocessing
    other._queue_run_command = self._queue_run_command
    other._force_wrapping = deepcopy(self._force_wrapping)
    other._warning_message = self._warning_message

    other.set_log(self.log)

    # Set up name
    if name is None:
      name = '%s_copy' %(self.name)
    other.set_name(name)

    other.set_verbose(self.verbose)
    if self._resolution:
      other.set_resolution(self._resolution)


  def _get_map_model_manager_with_selected(self,
      map_id_list=None, model_id_list = None,
      deep_copy = False):
    # Create a new map_model_manager with just what we need
    assert map_id_list # Need maps to create map_model_manager with selected
    working_mmm = map_model_manager(
      map_manager = self.get_any_map_manager(), log = self.log,
      verbose = self.verbose)
    self._set_default_parameters(working_mmm)


    if map_id_list:
      for id in map_id_list:
        if self.get_map_manager_by_id(id):
          working_mmm.add_map_manager_by_id(self.get_map_manager_by_id(id),id)
    if model_id_list:
      for id in model_id_list:
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

    if n_bins is None:
      if local_sharpen:
        n_bins = 20
      else:
        n_bins = 200

    min_n_bins = n_bins//3
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
      ):

    # scaling_group_info group_args object has direction vectors, list of si:
    #  scaling_group_info.direction_vectors
    #  scaling_group_info.scaling_info_list: one si entry per direction vector
    #    si.target_scale_factors
    #    si.target_sthol2
    #    si.d_min_list
    #    si.cc_list
    #    si.low_res_cc # low-res average

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
     model_for_rms_fc = None,
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
      get_scale_as_aniso_u = True, # always get it; use if get_scale_as_aniso_u
      use_dv_weighting = use_dv_weighting,
      run_analyze_anisotropy = run_analyze_anisotropy,
      expected_rms_fc_list = expected_rms_fc_list,
      resolution = resolution, # nominal resolution
      out = self.log)
    if not hasattr(result,'scaling_info_list'):  # result is one si
      result = group_args(
        scaling_info_list = [result],
        direction_vectors = direction_vectors,
        expected_rms_fc_list = expected_rms_fc_list,)

    # Set anything with too-low low-res CC to None for model-based run
    for si in result.scaling_info_list:
      if is_model_based and si.low_res_cc < minimum_low_res_cc:
        print("Skipping scaling with low_res_cc = %.2f" %(si.low_res_cc),
          file = self.log)
        si.target_scale_factors = None

    # result to be returned is a group_args object with
    #   direction vectors and list of si as scaling_info_list:
    #  result.direction_vectors
    #  result.scaling_info_list:
    #    si.target_scale_factors
    #    si.target_sthol2
    #    si.d_min_list
    #    si.cc_list
    #    si.low_res_cc # low-res average
    #    si.expected_rms_fc_list  # as supplied
    return result

  def _remove_temp_dir(self,temp_dir):
    if not os.path.isdir(temp_dir):
      return  # nothing to do
    else:  # remove it
     from shutil import rmtree
     rmtree(temp_dir)

  def _create_temp_dir(self, temp_dir):
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
      max_abs_b = None):

    if max_abs_b is None:  # set default
      max_abs_b = 100 * (self.resolution()/4)**2  # about 100 at 4 A
    map_coeffs = self.get_any_map_manager(
         ).map_as_fourier_coefficients(d_min = self.resolution())
    from cctbx.maptbx.segment_and_split_map import map_coeffs_as_fp_phi
    from cctbx.maptbx.refine_sharpening import get_nearest_lattice_points
    f_array,phases=map_coeffs_as_fp_phi(map_coeffs)

    # scale_factor_info.value_list is a set of scaling_group_info objects.
    # scale_factor_info.xyz_list are the coordinates where these apply
    # scale_factor_info.n_bins is number of bins
    # value_list is a set of scaling_group_info objects, one per xyz.
    #  scaling_group_info group_args object has direction vectors, list of si:
    #   scaling_group_info.direction_vectors
    #   scaling_group_info.scaling_info_list: one si entry per direction
    #    si.target_scale_factors
    #    si.target_sthol2
    #    si.d_min_list
    #    si.cc_list
    xyz_list = scale_factor_info.xyz_list
    value_list = scale_factor_info.value_list
    for xyz, scaling_group_info in zip(xyz_list,value_list):
      direction_vectors = scaling_group_info.direction_vectors
      scaling_info_list = scaling_group_info.scaling_info_list
      overall_u_cart_to_remove = scaling_group_info.overall_u_cart_to_remove
      for si,dv in zip(scaling_info_list,direction_vectors):
        if not si.target_scale_factors:
          continue
        target_sthol2 = si.target_sthol2
        # Now recalculate target_scale_factors from aniso_u_cart and dv

        recip_space_vectors = flex.vec3_double()
        scale_values = flex.double()
        for sthol2 in target_sthol2:
          s = (sthol2/0.25)**0.5
          recip_space_vectors.append(col(dv) * s)
          scale_values.append(1.)
        indices = flex.miller_index(tuple(get_nearest_lattice_points(
            f_array.unit_cell(),recip_space_vectors)))
        scale_values_array = f_array.customized_copy(
          data = scale_values,
          indices = indices)
        from mmtbx.scaling import absolute_scaling
        b_iso = abs(flex.double(adptbx.u_as_b(
            overall_u_cart_to_remove))[:3].min_max_mean().mean)

        if max_abs_b and abs(b_iso) > max_abs_b:
          ratio = max_abs_b/abs(b_iso)
          overall_u_cart_to_remove = tuple(ratio * col(
              overall_u_cart_to_remove))
        u_star= adptbx.u_cart_as_u_star(
          scale_values_array.unit_cell(),tuple(col(overall_u_cart_to_remove)))
        scaled_f_array = absolute_scaling.anisotropic_correction(
          scale_values_array,0.0, u_star ,must_be_greater_than=-0.0001)
        if 0:
          print("TARGET: BISO %.2f" %(b_iso),"BCART",str(overall_u_cart_to_remove),
            "".join(["%.2f " %(ts) for ts in si.target_scale_factors]))
          print("FITTED: ",
            "".join(["%.2f " %(sf) for sf in scaled_f_array.data()]))
        # Now extract our values as target_scale_factors
        si.target_scale_factors = scaled_f_array.data()

  def _analyze_aniso(self, scale_factor_info, map_id = None,
     replace_inside = None,
     replace_boundary = None,
     replace_outside = None):

    # Get a mask around the map
    mask_id = self._generate_new_map_id()
    self.create_mask_around_density(
      soft_mask  = True,
      mask_id = mask_id,
      map_id = map_id)
    mask_map_manager = self.get_map_manager_by_id(mask_id)

    # Analyze u_cart values in relation to mask
    mean_u_cart_dict={}
    mean_u_cart_dict_n={}
    for inside,cutoff_low,cutoff_high in zip(
      [True,False,None],
      [0.9,None,0.1],
      [None,0.1,0.9],):

      print("\nLocal B-cart by XYZ with inside_mask = %s:" %(
        inside), file = self.log)
      print("(Mask values with lower limit of %s and high limit of %s)" %(
         cutoff_low,cutoff_high), file = self.log)
      mean_u_cart_dict[inside] = col((0,0,0,0,0,0,))
      mean_u_cart_dict_n[inside] = 0
      working_scale_factor_info = self._get_scale_factor_info_inside_mask(
         scale_factor_info, mask_map_manager,
         cutoff_low = cutoff_low, cutoff_high = cutoff_high,
         skip_u_cart_of_zero = True)
      xyz_list = working_scale_factor_info.xyz_list
      if xyz_list.size() < 1: continue
      value_list = working_scale_factor_info.value_list
      uanisos = flex.sym_mat3_double()
      for i in range(xyz_list.size()):

        values = working_scale_factor_info.value_list[i]
        aniso_u_cart = values.overall_u_cart_to_remove
        uanisos.append(aniso_u_cart)
        b_cart = adptbx.u_as_b(aniso_u_cart)
        xyz=xyz_list[i]

        direction_vectors = values.direction_vectors
        scaling_info_list = values.scaling_info_list
        if 0:
          print ("XYZ (%7.1f,%7.1f,%7.1f) " %(tuple(xyz)),file=self.log)
          for dv,si in zip(direction_vectors,scaling_info_list):
           print("DV: (%7.1f,%7.1f,%7.1f)" %(dv) +
            " B (CC*): %.3f" %(si.effective_b),
            " B (FOBS): %.3f" %(si.effective_b_f_obs),file = self.log)


        mean_u_cart_dict[inside]+=col(aniso_u_cart)
        mean_u_cart_dict_n[inside]+=1
        if 0:
          print ("XYZ: "+
         "(%7.1f,%7.1f,%7.1f) B-cart:(%6.1f,%6.1f,%6.1f,%6.1f,%6.1f,%6.1f) " %(
          tuple(list(xyz)+list( b_cart))),
          file = self.log)
      from mmtbx.tls import tools
      cm = xyz_list.mean()
      result = tools.tls_from_uaniso_minimizer(
        uaniso         = uanisos,
        T_initial      = [0,0,0,0,0,0],
        L_initial      = [0,0,0,0,0,0],
        S_initial      = [0,0,0,0,0,0,0,0,0],
        refine_T       = True,
        refine_L       = True,
        refine_S       = True,
        origin         = cm,
        sites          = xyz_list,
        max_iterations = 100)

      T=result.T_min
      L=result.L_min
      S=result.S_min

      from mmtbx_tls_ext import tlso, uaniso_from_tls_one_group
      tlso_value = tlso(t = T, l = L, s = S, origin = cm)
      new_anisos= uaniso_from_tls_one_group(tlso = tlso_value,
         sites_cart = xyz_list,
         zeroize_trace=False)
      if inside and replace_inside:
        replace_u_cart_to_remove = True
        print("\nReplacing values inside mask with TLS-derived scale factors",
          file = self.log)
      else:
        replace_u_cart_to_remove = False
      for i in range(xyz_list.size()):
        values = working_scale_factor_info.value_list[i]
        aniso_u_cart = values.overall_u_cart_to_remove
        new_u_cart = new_anisos[i]
        if replace_u_cart_to_remove:
          values.overall_u_cart_to_remove = new_u_cart
        if 0:
          b_cart = adptbx.u_as_b(aniso_u_cart)
          new_b_cart = adptbx.u_as_b(new_u_cart)
          xyz=xyz_list[i]
          print ("NEW: "+
          "(%7.1f,%7.1f,%7.1f) B-cart:(%6.1f,%6.1f,%6.1f,%6.1f,%6.1f,%6.1f) " %(
           tuple(list(xyz)+list( new_b_cart))),
           file = self.log)
      print("\nMean anisotropy as TLS:",file = self.log)
      print("T: %s" %(str(T)),file = self.log)
      print("L: %s" %(str(L)), file = self.log)
      print("S: %s" %(str(S)), file = self.log)


    print("\nOverall average b_cart by region:",file = self.log)
    for inside in [True,False,None]:
      mean_u_cart_dict[inside] /= max(1,mean_u_cart_dict_n[inside])
      b_cart = adptbx.u_as_b(mean_u_cart_dict[inside])
      print("%6s   (n = %4s)   (%6.1f,%6.1f,%6.1f,%6.1f,%6.1f,%6.1f) " %(
        tuple([inside]+[mean_u_cart_dict_n[inside]]+list(b_cart))),
        file = self.log)

  def _run_group_of_anisotropic_sharpen(self,
      map_id  = 'map_manager',
      map_id_1 = 'map_manager_1',
      map_id_2 = 'map_manager_2',
      map_id_to_be_scaled_list = None,
      map_id_scaled_list = None,
      resolution = None,
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
      model_for_rms_fc = None,
      replace_aniso_with_tls_equiv = None,
      minimum_low_res_cc = None,
      max_abs_b = None,
      temp_dir = 'TEMP_ANISO_LOCAL',
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
      smoothing_radius = smoothing_radius,
      skip_boxes = True)

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
           map_model_manager = self,
           model = model_for_rms_fc)
        kw['expected_rms_fc_list'] = aavr.get_target_scale_factors()

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

    # scale_factor_info.value_list is a set of scaling_group_info objects.
    # scale_factor_info.xyz_list are the coordinates where these apply
    # scale_factor_info.n_bins is number of bins
    # value_list is a set of scaling_group_info objects, one per xyz.
    #  scaling_group_info group_args object has direction vectors, list of si:
    #   scaling_group_info.direction_vectors
    #   scaling_group_info.scaling_info_list: one si entry per direction
    #    si.target_scale_factors
    #    si.target_sthol2
    #    si.d_min_list
    #    si.cc_list
    #    si.low_res_cc # low-res average
    #    si.effective_b # B for fall-off of cc*
    #    si.effective_b_fobs  # B for fall-off of amplitudes of Fobs

    #  Have a look at scale values vs resolution along direction_vectors
    xyz_list = scale_factor_info.xyz_list
    if 0:
     for dv_id in range(direction_vectors.size()):
      print("\nScale for direction vector (%5.2f, %5.2f, %5.2f)" %(
        direction_vectors[dv_id]))
      for i in range(xyz_list.size()):
        xyz=xyz_list[i]
        print ("XYZ = (%7.1f, %7.1f, %7.1f)" %(xyz))
        values=flex.double()
        for i_bin in range(1,scale_factor_info.n_bins+1,3):
          scale_value_list,xyz_used_list = self._get_scale_values_for_bin(
            xyz_list=xyz_list,
            i_bin = i_bin,
            scale_factor_info = scale_factor_info,
            dv_id = dv_id)
          if scale_value_list.size() > 0:
            values.append(
              scale_value_list[max(0,min(
              i,scale_value_list.size()-1))]) # position i in xyz_list
          else:
            values.append(0)
        for value in values:
          print("%5.2f "  %(value), end="")
        print()
    self._analyze_aniso(scale_factor_info,map_id=map_id,
      replace_inside = (replace_aniso_with_tls_equiv and get_scale_as_aniso_u))
    if get_scale_as_aniso_u: # Summarize U vs xyz and vs inside/outside
      self._update_scale_factor_info_from_aniso(scale_factor_info,
        max_abs_b = max_abs_b)

    setup_info.kw = kw

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

  def _local_sharpen(self,
      map_id  = 'map_manager',
      map_id_1 = 'map_manager_1',
      map_id_2 = 'map_manager_2',
      map_id_to_be_scaled_list = None,
      map_id_scaled_list = None,
      resolution = None,
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
      model_for_rms_fc = None,
      replace_aniso_with_tls_equiv = None,
      anisotropic_sharpen = None,
      minimum_low_res_cc = None,
      max_abs_b = None,
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

    # NOTE: map starts out overall-sharpened.  Therefore approximate scale
    # factors in all resolution ranges are about 1.  use that as default

    if anisotropic_sharpen:  # run N times with different direction vectors
      self._run_group_of_anisotropic_sharpen(**kw)
      return

    # Get scale factors vs resolution and location
    scale_factor_info = self.local_fsc(
      direction_vectors = [None],
      return_scale_factors = True, **kw)

    # scale_factor_info.value_list is a set of scaling_group_info objects.
    # scale_factor_info.xyz_list are the coordinates where these apply
    # value_list is a set of scaling_group_info objects, one per xyz.
    #  scaling_group_info group_args object has direction vectors, list of si:
    #   scaling_group_info.direction_vectors
    #   scaling_group_info.scaling_info_list: one si entry per direction
    #    si.target_scale_factors
    #    si.target_sthol2
    #    si.d_min_list
    #    si.cc_list
    #    si.low_res_cc # low-res average

    xyz_list = scale_factor_info.xyz_list
    d_min = scale_factor_info.d_min
    smoothing_radius = scale_factor_info.setup_info.smoothing_radius
    assert n_bins == scale_factor_info.n_bins # must match

    # Get Fourier coefficients for maps based on map_id_to_be_scaled
    for id, new_id in zip(map_id_to_be_scaled_list, map_id_scaled_list):

      map_coeffs = self.get_map_manager_by_id(id
           ).map_as_fourier_coefficients(d_min = d_min)
      f_array_info = get_map_coeffs_as_fp_phi(map_coeffs, n_bins = n_bins,
         d_min = d_min)
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
        weight_mm = self._create_full_size_map_manager_with_value_list(
          xyz_list = xyz_used_list,
          value_list = scale_value_list,
          smoothing_radius = smoothing_radius,
          default_value = scale_value_list.min_max_mean().mean)

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
     scale_factor_info, mask_map_manager, cutoff_low = 0.5,
     cutoff_high = None, skip_u_cart_of_zero = None):
    new_xyz_list = flex.vec3_double()
    new_value_list = []
    for xyz, value in zip (scale_factor_info.xyz_list,
       scale_factor_info.value_list):
      if skip_u_cart_of_zero and flex.pow2(flex.double(
       tuple(value.overall_u_cart_to_remove))).min_max_mean().max ==0:
         pass
      else: #check if inside
        site_frac=mask_map_manager.crystal_symmetry(
          ).unit_cell().fractionalize(xyz)
        xx = mask_map_manager.map_data().tricubic_interpolation(site_frac)
        if ((cutoff_low is None) or (xx >= cutoff_low)) and \
         ((cutoff_high is None) or (xx < cutoff_high)):
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
      map_id_to_be_scaled_list = None, # NOTE: not used, just allows it in call
      map_id_scaled_list = None, # NOTE: not used, just allows it in call
      resolution = None,
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
      model_for_rms_fc = None,
      replace_aniso_with_tls_equiv = None,
      max_abs_b = None,
      n_bins_default = 2000):

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
    assert (
    (self.get_map_manager_by_id(map_id_1) or
        is_model_based or is_external_based) and
       self.get_map_manager_by_id(map_id_2))

    if n_bins is None:
      n_bins = n_bins_default
    if nproc is None and self._nproc:
      nproc = self._nproc
    elif nproc is None:
      nproc = 1

    # Get basic info including minimum_resolution (cutoff for map_coeffs)
    setup_info = self._get_box_setup_info(map_id_1, map_id_2,
      resolution, box_cushion,
      n_boxes,
      core_box_size,
      smoothing_radius = smoothing_radius)

    if spectral_scaling and (not expected_rms_fc_list):
        from cctbx.development.approx_amplitude_vs_resolution import \
          approx_amplitude_vs_resolution
        aavr = approx_amplitude_vs_resolution(
           d_min = setup_info.minimum_resolution,
           n_bins = n_bins,
           map_model_manager = self,
           model = model_for_rms_fc)
        expected_rms_fc_list = aavr.get_target_scale_factors()

    box_info = self.split_up_map_and_model_by_boxes(
      target_for_boxes = setup_info.n_boxes,
      box_cushion = setup_info.box_cushion,
      skip_empty_boxes = False,
      select_final_boxes_based_on_model = False, # required
      apply_box_info = False,
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
    box_info.model_for_rms_fc = model_for_rms_fc
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
    #  scaling_group_info group_args object has direction vectors, list of si:
    #   scaling_group_info.direction_vectors
    #   scaling_group_info.scaling_info_list: one si entry per direction
    #    si.target_scale_factors
    #    si.target_sthol2
    #    si.d_min_list
    #    si.cc_list
    #    si.low_res_cc # low-res average


    results.setup_info = setup_info
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
      smoothing_radius = setup_info.smoothing_radius)

  def _create_full_size_map_manager_with_value_list(self,
      xyz_list, value_list, smoothing_radius,
      default_value = None):

    # Now create a small map and fill in values
    volume_per_grid_point=self.crystal_symmetry().unit_cell(
        ).volume()/max(1,xyz_list.size())
    target_spacing = (volume_per_grid_point**0.33) * 0.5
    local_n_real=tuple([ max(1,int(0.5+1.5*a/target_spacing)) for
        a in self.crystal_symmetry().unit_cell().parameters()[:3]])
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
      box_cushion=None,
      n_boxes=None,
      core_box_size=None,
      smoothing_radius=None,
      skip_boxes = None,
      ):
    if not resolution:
      resolution = self.resolution()

    if not box_cushion:
      box_cushion = 1.5 * resolution

    if (not core_box_size):
      if n_boxes:
        volume = self.crystal_symmetry().unit_cell().volume()
        core_box_size=int(0.5+ volume/n_boxes)**0.33
      else:
        core_box_size = max(11, 3 * resolution)

    if (not skip_boxes) and (not n_boxes):
      volume = self.crystal_symmetry().unit_cell().volume()
      n_boxes = int(0.5+volume/(core_box_size)**3)
      print ("Target core_box_size: %.2s A  Target boxes: %s" %(
        core_box_size, n_boxes),file = self.log)

    if (not smoothing_radius) and (not skip_boxes):
      smoothing_radius = 0.5 * \
        self.crystal_symmetry().unit_cell().parameters()[0] * \
        min(1,core_box_size/self.get_any_map_manager().map_data().all()[0])
      smoothing_radius = max(2*resolution, smoothing_radius)

    # Working resolution is resolution * d_min_ratio
    minimum_resolution = self._get_d_min_from_resolution(resolution)

    return group_args(
     resolution = resolution,
     box_cushion = box_cushion,
     n_boxes = n_boxes,
     core_box_size = core_box_size,
     smoothing_radius = smoothing_radius,
     minimum_resolution = minimum_resolution,
      )

  def _run_fsc_in_boxes(self,
     nproc = None,
     box_info = None):

    assert box_info.n_bins is not None
    # Set up to run in each box
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
    #  scaling_group_info group_args object has direction vectors, list of si:
    #   scaling_group_info.direction_vectors
    #   scaling_group_info.scaling_info_list: one si entry per direction
    #    si.target_scale_factors
    #    si.target_sthol2
    #    si.d_min_list
    #    si.cc_list
    #    si.low_res_cc # low-res average

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
        other = f_map_2, bin_width = bin_width, fsc_cutoff = fsc_cutoff)

    return fsc_curve

  def map_map_cc(self,
      map_id = 'map_manager_1',
      other_map_id = 'map_manager_2',
      mask_id = None,
      mask_cutoff = 0.5,
      resolution = None):

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
      use_b_zero = True,
      ):

    if not model:
      model = self.get_model_by_id(model_id)
    if not model:
      return None
    if use_b_zero: # make a deep copy and set b values to zero
      print("Map-model CC using model with B values of zero", file = self.log)
      model = model.deep_copy()
      b_iso_values = model.get_b_iso()
      model.set_b_iso(flex.double(b_iso_values.size(),0.)) # XXX should be ok
      u_cart=model.get_xray_structure().scatterers().extract_u_cart(model.get_xray_structure().crystal_symmetry().unit_cell())
      model.get_xray_structure().set_u_cart(u_cart)
      model.set_xray_structure(model.get_xray_structure())
    else:
      print("Map-model CC using model with B values as is", file = self.log)

    map_manager= self.get_map_manager_by_id(map_id)
    assert model and map_manager
    if not resolution:
      resolution = self.resolution()
    assert resolution is not None
    print("Resolution for map-model CC: %.3f A" %(resolution), file = self.log)

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
    Identify rotation/translation to map model from other on to model in this
     object.
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
    Identify rotation/translation to map model from other on to model in this
     object.
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
      scattering_table = 'electron',
      fractional_error = 0.0,
      gridding = None,
      wrapping = False,
      map_id = None,
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

      If no map_manager is present, use supplied or existing model to
         generate map_manager and model.

      If map_manager is present, use supplied or existing model as model and
         create new entry in this this map_model_manager with name map_id.
         If map_id is None, use 'model_map'

      If no existing or supplied model, use default model from library,
      box with box_cushion around it and choose n_residues to
      include (default=10).

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
      map_id:  ID of map_manager to be created with model-map information (only
                 applies if there is an existing map_manager)
    '''


    # Set the resolution now if not already set
    if d_min and self.map_manager() and (not self.resolution()):
      self.set_resolution(d_min)

    # Get some value for resolution
    if not d_min:
      d_min = self.resolution()
    if not d_min:
      d_min = 3  # default


    self._print("\nGenerating new map data\n")
    if self.model() and (not model):
      self._print("NOTE: using existing model to generate map data\n")
      model = self.model()

    # See if we have a map_manager
    if self.map_manager():
      if not gridding:
        gridding = self.map_manager().map_data().all()
        origin_shift_grid_units = self.map_manager().origin_shift_grid_units
        self._print(
          "Using existing map_manager as source of gridding and origin")
      if not map_id:
        map_id = 'model_map'
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
    map_coeffs = generate_map_coefficients(model = model,
        d_min = d_min,
        k_sol = k_sol,
        b_sol = b_sol,
        scattering_table = scattering_table,
        log = null_out())

    mm = generate_map_data(
      map_coeffs = map_coeffs,
      d_min = d_min,
      gridding = gridding,
      wrapping = wrapping,
      origin_shift_grid_units = origin_shift_grid_units,
      high_resolution_real_space_noise_fraction = fractional_error,
      log = null_out())

    if self.get_any_map_manager():
      if not map_id:
        map_id = 'model_map'
      new_mm = self.get_any_map_manager().customized_copy(
        map_data=mm.map_data())
      self.add_map_manager_by_id(new_mm,map_id)
    else: # create map-model manager info
      self.set_up_map_dict(map_manager=mm)
      self.set_up_model_dict(model=model)

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
     log = self.log,
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

    # Shift origin of map_manager
    self._map_manager.shift_origin(desired_origin = desired_origin)

    # Shift origin of model  Note this sets model shift_cart
    if self._model:
      self._model = self.shift_model_to_match_working_map(
        coordinate_shift = shift_to_apply_cart,
        new_shift_cart = new_full_shift_cart,
        final_crystal_symmetry = self._map_manager.crystal_symmetry(),
        final_unit_cell_crystal_symmetry =
           self._map_manager.unit_cell_crystal_symmetry(),
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

  assert ((direction_vectors is None) or
     (list(direction_vectors) == list(scaling_group_info.direction_vectors))
     )

  # work on one location (xyz) with one value ( one scaling_group_info object)

  # Produce a set of xyz and a set of scaling_group_info objects

  #   scaling_group_info is group_args object.
  #   scaling_group_info.direction_vectors
  #   scaling_group_info.scaling_info_list: one si entry per direction
  #    si.target_scale_factors
  #    si.target_sthol2
  #    si.d_min_list
  #    si.cc_list
  #    si.low_res_cc # low-res average

  # If direction vectors are None then NCS operation just multiplies all the
  #   entries without changing them
  if not direction_vectors:
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
  #  rearranging only
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
    for j in range(n):
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
  from cctbx.maptbx import grid_indices_around_sites
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

  from iotbx.map_model_manager import map_model_manager as MapModelManager
  box_info = deepcopy(box_info)
  if first_to_use is not None and last_to_use is not None:
    for x in ['lower_bounds_with_cushion_list','upper_bounds_with_cushion_list',
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
  if not first_to_use:
    first_to_use = 1
  if not last_to_use:
    last_to_use = len(lower_bounds_list)
  for lower_bounds, upper_bounds, selection in zip(
       lower_bounds_list,
       upper_bounds_list,
       box_info.selection_list,):

    mmm=map_model_manager.extract_all_maps_with_bounds(
     lower_bounds, upper_bounds,
     model_can_be_outside_bounds = True)

    if mmm.model():
      model_to_keep = mmm.model().select(selection)
    else:
      model_to_keep = None
    if box_info.mask_around_unselected_atoms:  # mask everything we didn't keep
      # NOTE: only applies mask to map_manager, not any other map_managers
      remaining_model=mmm.model().select(~selection)
      nnn=mmm.deep_copy()
      nnn.set_model(remaining_model)
      nnn.remove_model_outside_map(boundary=box_info.mask_radius)
      if nnn.model().get_sites_cart().size() > 0: # do something
        nnn.create_mask_around_atoms(
         mask_atoms_atom_radius=box_info.mask_radius,
         mask_id = 'mask')
        mask_mm = nnn.get_map_manager_by_id(map_id = 'mask')
        s = (mask_mm.map_data() > 0.5)
        mmm.map_manager().map_data().set_selected(s,box_info.masked_value)
    if model_to_keep:
      mmm.set_model(model_to_keep)
    mmm_list.append(mmm)
  box_info.mmm_list = mmm_list
  return box_info

def get_selections_and_boxes_to_split_model(
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
        mask_radius = 3,
        masked_value = -10,
        get_unique_set_for_boxes = True,
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

  model = map_model_manager.model()
  map_manager = map_model_manager.get_any_map_manager()

  # Get the selections
  box_info = group_args(
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
  elif selection_method == 'by_chain':
    from mmtbx.secondary_structure.find_ss_from_ca import get_chain_ids
    for chain_id in get_chain_ids(model.get_hierarchy(), unique_only=True):
      if chain_id.replace(" ",""):
        selection = model.selection(" %s (chain %s) " %(
         info.no_water_or_het_with_and,chain_id))
      else:
        selection = model.selection(" %s " %(
         info.no_water_or_het))
      if (not skip_empty_boxes) or (selection.count(True) > 0):
        box_info.selection_list.append(selection)
  elif selection_method == 'by_segment':
    selection_strings= get_selections_for_segments(model,
    no_water_or_het_with_and = info.no_water_or_het_with_and)
    for selection_string in selection_strings:
      selection = model.selection(selection_string)
      if (not skip_empty_boxes) or (selection.count(True) > 0):
        box_info.selection_list.append(selection)
  elif selection_method == 'boxes':
    if info.no_water_or_het and info.no_water_or_het != 'all':
      overall_selection = model.selection("not (%s) " %(info.no_water_or_het))
    else:
      overall_selection = None

    # Get boxes without and with cushion (cushion may be None)
    box_info = map_manager.get_boxes_to_tile_map(
      target_for_boxes = target_for_boxes,
      box_cushion = box_cushion,
      get_unique_set_for_boxes = get_unique_set_for_boxes)

    # Select inside boxes without cushion and create cushion too
    box_info = get_selections_from_boxes(
       box_info = box_info,
       model = model,
       overall_selection = overall_selection,
       skip_empty_boxes = skip_empty_boxes)

  if select_final_boxes_based_on_model or (
     not box_info.lower_bounds_list): # get bounds now:
    from cctbx.maptbx.box import get_bounds_around_model
    box_info.lower_bounds_list = []
    box_info.upper_bounds_list = []
    for selection in box_info.selection_list:
      if model:
        model_use=model.select(selection)
      else:
        model_use = None
      info = get_bounds_around_model(
        map_manager = map_manager,
        model = model_use,
        box_cushion = box_cushion)
      box_info.lower_bounds_list.append(info.lower_bounds)
      box_info.upper_bounds_list.append(info.upper_bounds)
    box_info.lower_bounds_with_cushion_list = [] # not using these
    box_info.upper_bounds_with_cushion_list = []

  box_info.mask_around_unselected_atoms = mask_around_unselected_atoms
  box_info.mask_radius = mask_radius
  box_info.masked_value = masked_value
  return box_info



def get_selections_from_boxes(box_info = None,
    model = None,
    overall_selection = None,
    skip_empty_boxes = None,
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
    if (not sel) or (not skip_empty_boxes) or (sel.count(True) > 0):
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

def get_selection_inside_box(
     lower_bounds = None,
     upper_bounds = None,
     n_real = None,
     model = None,
     crystal_symmetry = None):
  '''
   get selection for all the atoms inside this box
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

def get_skip_waters_and_hetero_lines(skip_waters, skip_hetero):
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

def get_selections_for_segments(model, no_water_or_het_with_and = ''):
  '''
    Generate selections corresponding to each segment (chain or part of a chain
    that is separate from remainder of chain)
  '''
  assert isinstance(model, mmtbx.model.manager)

  from iotbx.pdb import resseq_encode
  selection_list = []
  ph = model.get_hierarchy()
  for m in ph.models()[:1]:
    for chain in m.chains():
      first_resno = None
      last_resno = None
      chain_id = chain.id
      previous_rg = None
      for rg in chain.residue_groups():
        if previous_rg and ( (not rg.link_to_previous) or (not
           residue_group_is_linked_to_previous(rg, previous_rg))):
          # break here
          selection_list.append("%s ( chain %s and resseq %s:%s ) " %(
           no_water_or_het_with_and,
            chain_id, resseq_encode(first_resno), resseq_encode(last_resno)))
          first_resno = None
          last_resno = None
        if not first_resno:
          first_resno = rg.resseq_as_int()
        last_resno = rg.resseq_as_int()
        previous_rg = rg
      if first_resno is not None and last_resno is not None:
        selection_list.append(" %s ( chain %s and resseq %s:%s ) " %(
            no_water_or_het_with_and, chain_id,
           resseq_encode(first_resno), resseq_encode(last_resno)))
  return selection_list

def residue_group_is_linked_to_previous(rg, previous_rg):
  from mmtbx.secondary_structure.find_ss_from_ca import is_close_to
  if is_close_to(rg,previous_rg):
    return True
  elif  rg.resseq_as_int()!=+previous_rg.resseq_as_int()+1:
    return True
  else:
    return False
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

class run_anisotropic_scaling_as_class:
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
    # scale_factor_info.n_bins is number of bins
    # value_list is a set of scaling_group_info objects, one per xyz.
    #  scaling_group_info group_args object has direction vectors, list of si:
    #   scaling_group_info.direction_vectors
    #   scaling_group_info.scaling_info_list: one si entry per direction
    #    si.target_scale_factors
    #    si.target_sthol2
    #    si.d_min_list
    #    si.cc_list
    #    si.low_res_cc # low-res average

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
      weight_mm = \
         map_model_manager._create_full_size_map_manager_with_value_list(
        xyz_list = xyz_used_list,
        value_list = scale_value_list,
        smoothing_radius = smoothing_radius,
        default_value = scale_value_list.min_max_mean().mean)

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
           model_for_rms_fc = self.box_info.model_for_rms_fc,
           replace_aniso_with_tls_equiv =
                self.box_info.replace_aniso_with_tls_equiv)
        if scaling_group_info:
          # scaling_group_info group_args object has direction vectors, list of si:
          #  scaling_group_info.direction_vectors
          #  scaling_group_info.scaling_info_list: one si entry per direction
          #    si.target_scale_factors
          #    si.target_sthol2
          #    si.d_min_list
          #    si.cc_list
          #    si.low_res_cc # low-res average
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
