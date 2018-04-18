from __future__ import division, print_function
'''
'''

from iotbx.data_manager.miller_array import MillerArrayDataManager
from iotbx.cif_mtz_data_labels import mtz_map_coefficient_labels, \
  cif_map_coefficient_labels
from iotbx.reflection_file_reader import any_reflection_file
from libtbx.utils import Sorry

# =============================================================================
class MapCoefficientsDataManager(MillerArrayDataManager):

  datatype = 'map_coefficients'

  # ---------------------------------------------------------------------------
  # Map coefficients

  def add_map_coefficients(self, filename, data):
    '''
    Stores the main any_file object so that all arrays are available
    '''
    return self._add(MapCoefficientsDataManager.datatype, filename, data)

  def set_default_map_coefficients(self, filename):
    return self._set_default(MapCoefficientsDataManager.datatype, filename)

  def get_map_coefficients(self, filename=None):
    '''
    Returns the main file object
    '''
    return self._get(MapCoefficientsDataManager.datatype, filename)

  def get_map_coefficients_arrays(self, filename=None):
    '''
    Returns a list of map coefficients from the file
    '''
    array_storage = '_%s_arrays' % MapCoefficientsDataManager.datatype
    if (hasattr(self, array_storage)):
      return getattr(self, array_storage)
    else:
      raise Sorry('No data file has been processed yet.')

  def get_map_coefficients_names(self):
    return self._get_names(MapCoefficientsDataManager.datatype)

  def get_default_map_coefficients_name(self):
    return self._get_default_name(MapCoefficientsDataManager.datatype)

  def remove_map_coefficients(self, filename):
    return self._remove(MapCoefficientsDataManager.datatype, filename)

  def has_map_coefficients(
      self, expected_n=1, exact_count=False, raise_sorry=False):
    return self._has_data(
      MapCoefficientsDataManager.datatype, expected_n=expected_n,
      exact_count=exact_count, raise_sorry=raise_sorry)

  def process_map_coefficients_file(self, filename):
    self.process_miller_array_file(filename)

  def filter_map_coefficients_arrays(self, filename):
    '''
    Function for checking labels in miller_arrays to determine type
    '''
    miller_arrays = self.get_miller_arrays(filename)
    self._map_coefficients_arrays = list()
    known_labels = mtz_map_coefficient_labels.union(cif_map_coefficient_labels)
    for miller_array in miller_arrays:
      labels = set(miller_array.info().labels)
      common_labels = known_labels.intersection(labels)
      if (len(common_labels) > 0):
        self._map_coefficients_arrays.append(miller_array)

    # if map coefficients exist, start tracking
    if (len(self._map_coefficients_arrays) > 1):
      self.add_map_coefficients(filename, self.get_miller_array(filename))

  def write_map_coefficients_file(
      self, filename, miller_arrays, overwrite=False):
    self.write_miller_array_file(
      filename=filename, miller_arrays=miller_arrays, overwrite=overwrite)

# =============================================================================
# end
