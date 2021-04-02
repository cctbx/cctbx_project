from __future__ import absolute_import, division, print_function
'''
'''
import os

from iotbx.data_manager import DataManagerBase
from libtbx import Auto
from libtbx.utils import Sorry

# =============================================================================
class RealMapDataManager(DataManagerBase):

  datatype = 'real_map'

  # ---------------------------------------------------------------------------
  # Real-space Maps
  def add_real_map(self, filename, data):
    return self._add(RealMapDataManager.datatype, filename, data)

  def set_default_real_map(self, filename):
    return self._set_default(RealMapDataManager.datatype, filename)

  def get_real_map(self, filename=None):
    return self._get(RealMapDataManager.datatype, filename)

  def get_real_map_names(self):
    return self._get_names(RealMapDataManager.datatype)

  def get_default_real_map_name(self):
    return self._get_default_name(RealMapDataManager.datatype)

  def remove_real_map(self, filename):
    return self._remove(RealMapDataManager.datatype, filename)

  def has_real_maps(self, expected_n=1, exact_count=False, raise_sorry=False):
    return self._has_data(RealMapDataManager.datatype, expected_n=expected_n,
                          exact_count=exact_count, raise_sorry=raise_sorry)

  def process_real_map_file(self, filename):
    return self._process_file(RealMapDataManager.datatype, filename)

  def get_default_output_real_map_filename(self):
    filename = self.get_default_output_filename()
    if not filename.endswith('.mrc') and not filename.endswith('.ccp4'):
      filename += '.ccp4'
    return filename

  def write_real_map_file(self, map_manager, filename=Auto, overwrite=Auto):

    # default options
    if (filename is Auto):
      filename = self.get_default_output_real_map_filename()
    if (overwrite is Auto):
      overwrite = self._overwrite

    # check arguments
    if (os.path.isfile(filename) and (not overwrite)):
      raise Sorry('%s already exists and overwrite is set to %s.' %
                  (filename, overwrite))

    try:
      map_manager.write_map(
        file_name = filename
      )
    except IOError as err:
      raise Sorry('There was an error with writing %s.\n%s' %
                  (filename, err))

    self._output_files.append(filename)
    self._output_types.append(RealMapDataManager.datatype)

    return filename

# =============================================================================
# end
