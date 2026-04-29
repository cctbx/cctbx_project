from __future__ import absolute_import, division, print_function
'''
'''

from iotbx.file_reader import any_file
from iotbx.data_manager import DataManagerBase
from libtbx import Auto
from libtbx.utils import Sorry

# =============================================================================
class RestraintDataManager(DataManagerBase):

  datatype = 'restraint'

  # ---------------------------------------------------------------------------
  # Restraints
  def add_restraint(self, filename, data):
    return self._add(RestraintDataManager.datatype, filename, data)

  def set_default_restraint(self, filename):
    return self._set_default(RestraintDataManager.datatype, filename)

  def get_restraint(self, filename=None):
    return self._get(RestraintDataManager.datatype, filename)

  def get_restraint_names(self):
    return self._get_names(RestraintDataManager.datatype)

  def get_default_restraint_name(self):
    return self._get_default_name(RestraintDataManager.datatype)

  def remove_restraint(self, filename):
    return self._remove(RestraintDataManager.datatype, filename)

  def has_restraints(self, expected_n=1, exact_count=False, raise_sorry=False):
    return self._has_data(RestraintDataManager.datatype, expected_n=expected_n,
                          exact_count=exact_count, raise_sorry=raise_sorry)

  def process_restraint_file(self, filename):
    if (filename not in self.get_restraint_names()):
      a = any_file(filename)
      if (a.file_type != 'cif'):
        raise Sorry('%s is not a recognized restraints file' % filename)
      else:
        self.add_restraint(filename, a.file_object.model())
    return filename

  def get_default_output_restraint_filename(self):
    filename = self.get_default_output_filename()
    if not filename.endswith('.cif'):
      filename += '.cif'
    return filename

  def write_restraint_file(self, restraint_str, filename=Auto, overwrite=Auto):
    if filename is Auto:
      filename = self.get_default_output_restraint_filename()
    return self._write_text(RestraintDataManager.datatype, restraint_str,
                            filename=filename, overwrite=overwrite)

# =============================================================================
# end
