'''
'''

from iotbx.data_manager import DataManagerBase
from libtbx import Auto

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

  def process_restraint_file(self, filename, cif_engine='xcif'):
    if (filename not in self.get_restraint_names()):
      from iotbx.file_io import read_file
      from iotbx.file_io.detection import _cif_block_datatypes
      result = read_file(filename, file_type='restraint', cif_engine=cif_engine)
      cif_model = result.file_object.model()
      # keep what a restraint consumer (mon_lib_srv, ener_lib) reads: drop the
      # model / reflection blocks of a combined CIF, typically the bulk of the
      # file (the atom table), unless a block is also a restraint block
      for block_name in list(cif_model.keys()):
        types = _cif_block_datatypes(block_name, cif_model[block_name])
        if ('restraint' not in types) and (types & {'model', 'miller_array'}):
          del cif_model[block_name]
      self.add_restraint(filename, cif_model)
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
