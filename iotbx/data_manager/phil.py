from __future__ import absolute_import, division, print_function
'''
'''

from iotbx.data_manager import DataManagerBase
from libtbx import Auto

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

  def get_default_output_phil_filename(self):
    filename = self.get_default_output_filename()
    if not filename.endswith('.eff'):
      filename += '.eff'
    return filename

  def write_phil_file(self, phil_str, filename=Auto, overwrite=Auto):
    # use this instead of libtbx.phil.scope.show for consistent error messages
    if filename is Auto:
      filename = self.get_default_output_phil_filename()
    return self._write_text(PhilDataManager.datatype, phil_str,
                            filename=filename, overwrite=overwrite)

# =============================================================================
# end
