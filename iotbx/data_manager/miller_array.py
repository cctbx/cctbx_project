from __future__ import division, print_function
'''
'''

from iotbx.data_manager import DataManagerBase

# =============================================================================
class MillerArrayDataManager(DataManagerBase):

  datatype = 'miller_array'

  # ---------------------------------------------------------------------------
  # Miller arrays

  def add_miller_array(self, filename, data):
    return self._add(MillerArrayDataManager.datatype, filename, data)

  def set_default_miller_array(self, filename):
    return self._set_default(MillerArrayDataManager.datatype, filename)

  def get_miller_array(self, filename=None):
    '''
    Returns the main file object
    '''
    return self._get(MillerArrayDataManager.datatype, filename)

  def get_miller_arrays(self, filename=None):
    '''
    Returns a list of all arrays
    '''
    return self.get_miller_array(filename=filename).as_miller_arrays()

  def get_miller_array_names(self):
    return self._get_names(MillerArrayDataManager.datatype)

  def get_default_miller_array_name(self):
    return self._get_default_name(MillerArrayDataManager.datatype)

  def remove_miller_array(self, filename):
    return self._remove(MillerArrayDataManager.datatype, filename)

  def has_miller_arrays(self, expected_n=1, exact_count=False, raise_sorry=False):
    return self._has_data(MillerArrayDataManager.datatype, expected_n=expected_n,
                          exact_count=exact_count, raise_sorry=raise_sorry)

  def process_miller_array_file(self, filename):
    self._process_file(MillerArrayDataManager.datatype, filename)

    # filter arrays (e.g self.filter_map_coefficients_arrays)
    for datatype in self.datatypes:
      function_name = 'filter_%s_arrays' % datatype
      if (hasattr(self, function_name)):
        getattr(self, function_name)(filename)

  def write_miller_array_file(self, filename, miller_arrays, overwrite=False):
    raise NotImplementedError

# =============================================================================
# end
