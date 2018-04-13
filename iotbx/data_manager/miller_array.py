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
    Returns a list of miller arrays from the file
    '''
    return self._get(MillerArrayDataManager.datatype, filename).\
      as_miller_arrays()

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
    return self._process_file(MillerArrayDataManager.datatype, filename)

  def write_miller_array_file(self, filename, miller_array_str, overwrite=False):
    raise NotImplementedError

# =============================================================================
# end
