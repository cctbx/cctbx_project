'''
Extra functions for DataManager
'''
from __future__ import absolute_import, division, print_function
from iotbx.data_manager.map_coefficients import MapCoefficientsDataManager
from iotbx.data_manager.real_map import RealMapDataManager

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
