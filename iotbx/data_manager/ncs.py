from __future__ import division, print_function
'''
'''

from iotbx.data_manager import DataManagerBase

# =============================================================================
class NcsDataManager(DataManagerBase):

  datatype = 'ncs'

  # ---------------------------------------------------------------------------
  # Ncs spec
  def add_ncs(self, filename, data):
    return self._add(NcsDataManager.datatype, filename, data)

  def set_default_ncs(self, filename):
    return self._set_default(NcsDataManager.datatype, filename)

  def get_ncs(self, filename=None):
    return self._get(NcsDataManager.datatype, filename)

  def get_ncs_names(self):
    return self._get_names(NcsDataManager.datatype)

  def get_default_ncs_name(self):
    return self._get_default_name(NcsDataManager.datatype)

  def remove_ncs(self, filename):
    return self._remove(NcsDataManager.datatype, filename)

  def has_ncs(self, expected_n=1, exact_count=False, raise_sorry=False):
    return self._has_data(NcsDataManager.datatype, expected_n=expected_n,
                          exact_count=exact_count, raise_sorry=raise_sorry)

  def process_ncs_file(self, filename):
    return self._process_file(NcsDataManager.datatype, filename)

  def write_ncs_file(self, filename, ncs_str, overwrite=False):
    self._write_text(NcsDataManager.datatype, filename,
                     ncs_str, overwrite=overwrite)

# =============================================================================
# end
