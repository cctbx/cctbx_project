'''
Extra functions for DataManager
'''
from __future__ import absolute_import, division, print_function
from iotbx.data_manager.map_coefficients import MapCoefficientsDataManager
from iotbx.data_manager.real_map import RealMapDataManager
from iotbx.map_model_manager import map_model_manager
from libtbx import Auto
from libtbx.utils import Sorry
import mmtbx.utils
from cctbx import crystal
from iotbx import crystal_symmetry_from_any
from iotbx import extract_xtal_data
import libtbx.phil
from copy import deepcopy

# =============================================================================
# general functions for scattering type (applies to model and/or miller_array)
class scattering_table_mixins(object):
  '''
  Functions for mapping scattering table types (wk1995 it1992 n_gaussian
  neutron electron) to model and array types (x_ray, electron, neutron).
  These will handle the mapping between the X-ray scattering tables
  (wk1995 it1992 n_gaussian) to the x_ray type.
  '''
  possible_scattering_table_types = ['wk1995', 'it1992', 'n_gaussian',
    'x_ray', 'neutron',  'electron']

  def check_scattering_table_type(self, scattering_table):
    '''
    Checks that the argument is a valid scattering table type

    Parameters
    ----------
    scattering_table : str
      The scattering table type
    data_type : str, optional
      The data type, e.g. 'model' or 'array'. This is for the error message

    Returns
    -------
    scattering_table : str
      The input scattering table. If x_ray is provided, n_gaussian is returned
      as the default scattering table.
    '''
    if scattering_table not in self.possible_scattering_table_types:
      raise Sorry('Unrecognized scattering table type, "%s," possible choices are %s.' %
                  (scattering_table, ', '.join(self.possible_scattering_table_types)))
    if scattering_table == 'x_ray':
      return 'n_gaussian'
    return scattering_table

  def map_scattering_table_type(self, scattering_table):
    '''
    Returns the appropriate model/array type based on the scattering table

    Parameters
    ----------
    scattering_table : str
      The scattering table type

    Returns
    -------
    mapped_type : str
      The mapped type
    '''
    self.check_scattering_table_type(scattering_table)
    if scattering_table in ['wk1995', 'it1992', 'n_gaussian']:
      return 'x_ray'
    else:
      return scattering_table

# -----------------------------------------------------------------------------
# extra functions for model and reflections
class fmodel_mixins(object):
  '''
  Function to extract fmodel when the DataManager supports both the "model"
  and "miller_array" data types.
  '''

  def update_all_defaults(self, data_type):
    '''
    Convenience function for setting the data_type of all models and
    Miller arrays. This sets the default as well as the type for all
    currently loaded files.

    Parameters
    ----------
    data_type : str
        The type to set (e.g. 'x_ray', 'neutron', or 'electron')
    '''
    # model
    checked_type = self.map_scattering_table_type(data_type)
    model_type = [checked_type]
    self._is_valid_model_type(model_type)
    self.set_default_model_type(model_type)
    for filename in self.get_model_names():
      if 'reference' in self.get_model_type(filename):
        self.set_model_type(filename, [checked_type, 'reference'])
      else:
        self.set_model_type(filename, model_type)

    # miller_array
    array_type = [checked_type]
    self.set_default_miller_array_type(array_type)
    for filename in self.get_miller_array_names():
      for label in self.get_miller_array_all_labels(filename):
        self.set_miller_array_type(filename, label, array_type)

  def get_fmodel_params(self):
    '''
    Return the fmodel parameters as a libtbx.phil.extract object
    '''
    if self._fmodel_phil_scope is None:
      return self.export_phil_scope(as_extract=True).data_manager.fmodel
    return self._fmodel_phil_scope.extract().data_manager.fmodel

  def set_fmodel_params(self, phil_extract):
    '''
    Set the fmodel parameters. The phil_extract can be the full DataManager
    PHIL extract or the output from get_fmodel_params
    '''
    full_extract = phil_extract
    if not hasattr(full_extract, 'data_manager'):
      full_extract = self.master_phil.extract()
      full_extract.data_manager.fmodel = phil_extract
    self._fmodel_phil_scope = self.master_phil.format(python_object=full_extract)

  def _resolve_symmetry_conflicts(self, model, reflection_file_server,
                                  params=None):
    '''
    Use logic of crystal.select_crystal_symmetry to select consensus crystal
    symmetry from multiple sources.
    '''
    if(reflection_file_server is None): return
    crystal_symmetries_from_coordinate_file = []
    crystal_symmetries_from_reflection_file = []
    rfns = []
    rfns.extend(reflection_file_server.file_name_miller_arrays.keys())
    crystal_symmetry_from_any.extract_and_append(
      file_names  = rfns,
      target_list = crystal_symmetries_from_reflection_file)
    crystal_symmetries_from_coordinate_file.append(model.crystal_symmetry())
    from_parameter_file = None
    if(params is not None):
      from_parameter_file = crystal.symmetry(
          unit_cell        = params.unit_cell,
          space_group_info = params.space_group)
    crystal_symmetry = crystal.select_crystal_symmetry(
      from_parameter_file   = from_parameter_file,
      from_coordinate_files = crystal_symmetries_from_coordinate_file,
      from_reflection_files = crystal_symmetries_from_reflection_file)
    model.set_crystal_symmetry(crystal_symmetry = crystal_symmetry)
    if(reflection_file_server.crystal_symmetry is None or not
       crystal_symmetry.is_similar_symmetry(
       reflection_file_server.crystal_symmetry)):
      reflection_file_server.update_crystal_symmetry(
        crystal_symmetry = model.crystal_symmetry())

  def get_fmodel(self,
                 crystal_symmetry = None,
                 experimental_phases_params = None,# XXX Need to be part of 'parameters'
                 scattering_table = None,
                 mask_params = None,
                 sf_accuracy_params = None,
                 free_r_flags_scope = 'miller_array.labels.name',
                 model_filename = None,
                 ):
    """
    Create mmtbx.fmodel.manager object using atomic model and diffraction data.
    crystal_symmetry: comes as cctbx.crystal.symmetry or Phil scope.
    scattering_table will trigger the use of model and reflections with corrct
    array_type.
    free_r_flags_scope defines the how-to message given if no free_r flags
      are found
    """
    array_type = self.map_scattering_table_type(scattering_table)
    scattering_table = self.check_scattering_table_type(scattering_table)
    crystal_symmetry_phil = crystal_symmetry
    if(crystal_symmetry is not None):
      if(isinstance(crystal_symmetry, libtbx.phil.scope_extract)):
        crystal_symmetry = crystal.symmetry(
          unit_cell        = crystal_symmetry.unit_cell,
          space_group_info = crystal_symmetry.space_group)
      else:
        assert isinstance(crystal_symmetry, libtbx.phil.scope_extract)
    # Gather models of appropriate type
    models = []
    if model_filename:
      models.append(self.get_model(model_filename,model_type=array_type))
    else:
      for filename in self.get_model_names(model_type=array_type):
        models.append(self.get_model(filename, model_type=array_type))
    if(len(models) == 0):
      raise Sorry("No model of '%s' type found to make fmodel."%array_type)
    if(len(models) > 1):
      raise Sorry("More than one model of '%s' type found."%array_type)
    model = models[0]
    # Get reflection file server
    rfs = self.get_reflection_file_server(
      array_type       = array_type,
      crystal_symmetry = crystal_symmetry,
      ignore_intensities_if_amplitudes_present = True)
    if rfs is None:
      raise Sorry("No reflection data provided.")
    # Resolve symmetry issues (in-place)
    self._resolve_symmetry_conflicts(
      params                 = crystal_symmetry_phil,
      model                  = model,
      reflection_file_server = rfs)
    #
    fmodel_params = self.get_fmodel_params()

    if array_type == 'neutron':
      parameters = fmodel_params.neutron_data
    else:
      parameters = fmodel_params.xray_data

    #print("LOOK : parameters.r_free_flags.required", parameters.r_free_flags.required)
    #print("LOOK : parameters.force_anomalous_flag_to_be_equal_to", parameters.force_anomalous_flag_to_be_equal_to)

    #
    # XXX
    # XXX Temporary hack/work-around (REMOVE later) start
    # XXX
    tmp_p = deepcopy(parameters)
    tmp_p.__inject__("file_name", None)
    tmp_p.__inject__("labels", None)
    tmp_p.r_free_flags.__inject__("file_name", None)
    tmp_p.r_free_flags.__inject__("label", None)
    # XXX
    # XXX Temporary hack/work-around (REMOVE later) end
    # XXX
    # Get reflection data
    dpg = None
    if(model.crystal_symmetry() is not None):
      dpg = model.crystal_symmetry().space_group().build_derived_point_group()
    # if both low_resolution and high_resolution are set, check that low_resolution > high_resolution
    if parameters.low_resolution is not None and parameters.high_resolution is not None \
      and parameters.low_resolution < parameters.high_resolution:
        raise Sorry('The low_resolution parameter is less than the high_resolution parameter. Please swap those values.')
    data = extract_xtal_data.run(
      keep_going                        = not tmp_p.r_free_flags.required,
      extract_r_free_flags              = not tmp_p.r_free_flags.ignore_r_free_flags,
      reflection_file_server            = rfs,
      parameters                        = tmp_p,
      experimental_phases_params        = experimental_phases_params,
      working_point_group               = dpg,
      free_r_flags_scope                = free_r_flags_scope,
      remark_r_free_flags_md5_hexdigest = model.get_header_r_free_flags_md5_hexdigest()).result()
    #
    # Set DataManager parameters extracted from inputs
    #
    # Extract and set twin_law
    if parameters.twin_law is None or parameters.twin_law is Auto:
      parameters.twin_law = model.twin_law_from_model_input()
    # Set test flag value
    parameters.r_free_flags.test_flag_value = data.test_flag_value
    # Load all back
    self.set_fmodel_params(fmodel_params)
    #
    if(len(data.err)>0):
      raise Sorry("\n".join(data.err))
    if(data.f_obs is None):
      raise Sorry("Diffraction data are not available to make fmodel.")
    # Setup scattering table of xray_structure
    model.setup_scattering_dictionaries(
      scattering_table = scattering_table,
      d_min            = data.f_obs.d_min())
    # Create and return fmodel
    twin_law = fmodel_params.xray_data.twin_law
    fmodel = mmtbx.utils.fmodel_manager2(
      f_obs               = data.f_obs,
      r_free_flags        = data.r_free_flags,
      abcd                = data.experimental_phases,
      xray_structure      = model.get_xray_structure(),
      twin_law            = twin_law,
      mask_params         = mask_params,
      sf_accuracy_params  = sf_accuracy_params,
      ignore_r_free_flags = parameters.r_free_flags.ignore_r_free_flags,
      mtz_object          = data.mtz_object,
      data_type           = array_type)
    return fmodel

# -----------------------------------------------------------------------------
# extra functions for maps
class map_mixins(object):
  '''
  Functions that are available when the DataManager supports both the
  "real_map" and "map_coefficients" data types.
  '''
  def has_real_maps_or_map_coefficients(
    self, expected_n=1, exact_count=False, raise_sorry=False):
    '''
    Combination of has_real_maps and has_map_coefficients
    '''
    n_real_maps = len(self._get_names(RealMapDataManager.datatype))
    n_map_coefficients = len(self._get_names(MapCoefficientsDataManager.datatype))
    actual_n = n_real_maps + n_map_coefficients
    datatype = 'real_maps or map coefficient'
    return self._check_count(datatype, actual_n, expected_n, exact_count, raise_sorry)

# -----------------------------------------------------------------------------
# extra functions for real maps
class real_map_mixins(object):
  '''
  Functions that are available when the DataManager supports the
  "real_map" data type.
  '''
  def get_map_model_manager(
    self,
    model_file=None,
    map_files=None,
    from_phil=False,
    guess_files=True,
    files_to_exclude = None,
    **kwargs):
    '''
    A convenience function for constructing a map_model_manager from the
    files in the DataManager.

    Parameters
    ==========
      model_file: str
        The name of the model file
      map_files: str or list
        The name(s) of the map files. If there is only one map, the name (str)
        is sufficient. A list is expected to have only 1 (full) or 2 (two half-
        maps) or 3 (full and 2 half-maps) map files. If three map files
        are in the list,
        the first map file is assumed to be the full map.
      from_phil: bool
        If set to True, the model and map names are retrieved from the
        standard PHIL names. The model_file and map_files parameters must
        be None if this parameter is set to True.
      files_to_exclude: str or list
        Any files listed will not be used from data_manager. For example,
          this can be used to exclude a file from being considered a half map
      **kwargs: keyworded arguments
        Extra keyworded arguments for map_model_manager constructor

    Return
    ======
      map_model_manager object
    '''

    # get filenames from PHIL
    map_model = None
    if map_files and (not isinstance(map_files,list)):
      map_files = [map_files]

    if from_phil:
      if model_file is not None or map_files is not None:
        raise Sorry(
          'If from_phil is set to True, model_file and map_files must be None.')

      params = self._program.params
      if hasattr(params, 'map_model'):
        map_model = params.map_model
      elif hasattr(params, 'input_files') and hasattr(
          params.input_files, 'map_model'):
        map_model = params.input_files.map_model
      else:
        raise Sorry('Program does not have the "map_model" PHIL scope.')

      model_file = getattr(map_model,'model', None)

      map_files = []
      full_map = getattr(map_model, 'full_map', None)
      if full_map is not None:
        map_files.append(full_map)
      half_maps = getattr(map_model, 'half_map', None)

      # Catch case where full_map is also present in half_maps as the only
      #   half map
      maps_to_exclude = [full_map] if full_map else []
      if files_to_exclude:
        if isinstance(files_to_exclude,list):
          maps_to_exclude += files_to_exclude
        else:
          maps_to_exclude.append(files_to_exclude)
      for fn in (maps_to_exclude if maps_to_exclude else []):
        if fn and (len(half_maps) == 1) and (half_maps[0] == fn):
          half_maps = []
      if half_maps:
        if len(half_maps) != 2:
          raise Sorry('Please provide 2 half-maps or one full map.')
        map_files += half_maps
    # If we didn't get anything, try looking directly at the
    #  available maps and models. If there are 1, 2 or 3 maps and 1 model,
    #  take them
    if guess_files and self.supports('model'):
      if (not model_file) and self.get_model_names() and len(self.get_model_names()) == 1:
        model_file = self.get_default_model_name()
        if map_model:
          map_model.model = model_file
    if guess_files and self.supports('real_map'):
      if (not map_files) and self.get_real_map_names() and len(self.get_real_map_names()) == 1:
        map_files = [self.get_default_real_map_name()]

      elif len(self.get_real_map_names()) in [2,3]:
        map_files = self.get_real_map_names()
    if isinstance(map_files, list):
      new_map_files = []
      for fn in map_files:
        if (not files_to_exclude) or (not fn in files_to_exclude):
          new_map_files.append(fn)
      map_files = new_map_files
    # check map_files argument
    mm = None
    mm_1 = None
    mm_2 = None
    if isinstance(map_files, list):
      if len(map_files) != 0 and \
         len(map_files) != 1 and len(map_files) != 2 and len(map_files) != 3:
        msg = 'Please provide only 1 full map or 2 half maps or 1 ' +\
         'full map and 2 half maps.\n Found:\n'
        for map_file in map_files:
          msg += ('  {map_file}\n'.format(map_file=map_file))
        raise Sorry(msg)
      elif len(map_files) == 0:
        mm = None
      elif len(map_files) == 1:
        mm = self.get_real_map(map_files[0])
        if map_model:
          map_model.full_map = map_files[0]
      elif len(map_files) == 2:
        mm_1 = self.get_real_map(map_files[0])
        mm_2 = self.get_real_map(map_files[1])
        if map_model:
          map_model.half_map = map_files
      elif len(map_files) == 3:
        if map_model:
          map_model.full_map = map_files[0]
          map_model.half_map = map_files[1:3]
        mm = self.get_real_map(map_files[0])
        mm_1 = self.get_real_map(map_files[1])
        mm_2 = self.get_real_map(map_files[2])
    elif map_files:
      mm = self.get_real_map(map_files) # it is a single file name
      if map_model:
        map_model.full_map = map_files
    else:
      mm = None

    if model_file:
      model = self.get_model(model_file)
    else:
      model = None

    # Check to make sure all map_managers are similar
    managers = [mm, mm_1,mm_2]
    comp_mm = None
    for m in managers:
      if not m: continue
      if m and (not comp_mm):
        comp_mm = m
      else:
        if (not comp_mm.is_similar(m)):
          raise Sorry("Input map files need to all have the same origin, gridding, and dimensions")
    mmm = map_model_manager(model=model, map_manager=mm, map_manager_1=mm_1,
      map_manager_2=mm_2, **kwargs)

    # clean up so that another read of maps and model will read again (these
    # are shifted when map_model_manager is called)
    if isinstance(map_files, list):
      for file_name in map_files[:3]:
        if file_name and file_name in self.get_real_map_names():
          self.remove_real_map(file_name)
    elif map_files:
          self.remove_real_map(map_files)
    if model_file:
      self.remove_model(model_file)

    return mmm

# -----------------------------------------------------------------------------
# extra functions for models and real_maps
class map_model_mixins(object):
  '''
  Functions that are available when the DataManager supports both the
  "model" and "real_map" data types.
  '''
  def remove_maps_and_models(self):
    ''' Remove all existing maps and models so they are not used by default'''

    for file_name in self.get_real_map_names():   # list of previously read maps
      self.remove_real_map(file_name)   # forget previous reads
    for file_name in self.get_model_names():
      self.remove_model(file_name)   # forget previous reads
