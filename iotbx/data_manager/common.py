'''
Extra functions for DataManager
'''
from __future__ import absolute_import, division, print_function
from iotbx.data_manager.map_coefficients import MapCoefficientsDataManager
from iotbx.data_manager.real_map import RealMapDataManager
from iotbx.map_model_manager import map_model_manager
from libtbx.utils import Sorry

# =============================================================================
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

  def get_map_model_manager(
    self, model_file=None, map_files=None, from_phil=False, **kwargs):
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
      **kwargs: keyworded arguments
        Extra keyworded arguments for map_model_manager constructor

    Return
    ======
      map_model_manager object
    '''

    # get filenames from PHIL
    if from_phil:
      if model_file is not None or map_files is not None:
        raise Sorry(
          'If from_phil is set to True, model_file and map_files must be None.')

      params = self._program.params
      map_model = None
      if hasattr(params, 'map_model'):
        map_model = params.map_model
      elif hasattr(params, 'input_files') and hasattr(
          params.input_files, 'map_model'):
        map_model = params.input_files.map_model
      else:
        raise Sorry('Program does not have the "map_model" PHIL scope.')

      model_file = map_model.model

      map_files = []
      full_map = getattr(map_model, 'full_map', None)
      if full_map is not None:
        map_files.append(full_map)
      half_maps = getattr(map_model, 'half_map', None)
      if half_maps:
        if len(half_maps) != 2:
          raise Sorry('Please provide 2 half-maps.')
        map_files += half_maps

    # If we didn't get anything, try looking directly at the
    #  available maps and models. If there are 1, 2 or 3 maps and 1 model,
    #  take them
    if (not model_file) and self.get_model_names() and \
         len(self.get_model_names()) == 1:
      model_file = self.get_default_model_name()
    if not map_files and self.get_real_map_names():
      if len(self.get_real_map_names()) == 1:
        map_files = self.get_default_real_map_name()
      elif len(self.get_real_map_names()) in [2,3]:
        map_files = self.get_real_map_names()

    # check map_files argument
    mm = None
    mm_1 = None
    mm_2 = None
    if isinstance(map_files, list):
      if len(map_files) != 1 and len(map_files) != 2 and len(map_files) != 3:
        msg = 'Please provide only 1 full map or 2 half maps or 1 ' +\
         'full map and 2 half maps.\n Found:\n'
        for map_file in map_files:
          msg += ('  {map_file}\n'.format(map_file=map_file))
        raise Sorry(msg)
      if len(map_files) == 1:
        mm = self.get_real_map(map_files[0])
      elif len(map_files) == 2:
        mm_1 = self.get_real_map(map_files[0])
        mm_2 = self.get_real_map(map_files[1])
      elif len(map_files) == 3:
        mm = self.get_real_map(map_files[0])
        mm_1 = self.get_real_map(map_files[1])
        mm_2 = self.get_real_map(map_files[2])
    elif map_files:
      mm = self.get_real_map(map_files)
    else:
      mm = None

    if model_file:
      model = self.get_model(model_file)
    else:
      model = None

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

