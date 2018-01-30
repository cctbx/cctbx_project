from __future__ import division, print_function
'''
'''

from libtbx.data_manager import DataManagerBase

# =============================================================================
class PhilDataManager(DataManagerBase):

  datatype = 'phil'

  # ---------------------------------------------------------------------------
  # Phils
  def add_phil(self, filename, data):
    return self._add(PhilDataManager.datatype, filename, data)

  def set_default_phil(self, filename):
    return self._set_default(PhilDataManager.datatype, filename)

  def get_phil(self, filename=None):
    return self._get(PhilDataManager.datatype, filename)

  def get_phil_names(self):
    return self._get_names(PhilDataManager.datatype)

  def get_default_phil_name(self):
    return self._get_default_name(PhilDataManager.datatype)

  def remove_phil(self, filename):
    return self._remove(PhilDataManager.datatype, filename)

  def has_phils(self, expected_n=1, exact_count=False, raise_sorry=False):
    return self._has_data(PhilDataManager.datatype, expected_n=expected_n,
                          exact_count=exact_count, raise_sorry=raise_sorry)

  def process_phil_file(self, filename):
    return self._process_file(PhilDataManager.datatype, filename)

  def write_phil_file(self, filename, phil_str, overwrite=False):
    # use this instead of libtbx.phil.scope.show for consistent error messages
    self._write_text(PhilDataManager.datatype, filename,
                     phil_str, overwrite=overwrite)

# =============================================================================
# end
