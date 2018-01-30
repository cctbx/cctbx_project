from __future__ import division, print_function
'''
'''

from libtbx.data_manager import DataManagerBase

# =============================================================================
class SequenceDataManager(DataManagerBase):

  datatype = 'sequence'

  # ---------------------------------------------------------------------------
  # Sequences
  def add_sequence(self, filename, data):
    return self._add(SequenceDataManager.datatype, filename, data)

  def set_default_sequence(self, filename):
    return self._set_default(SequenceDataManager.datatype, filename)

  def get_sequence(self, filename=None):
    return self._get(SequenceDataManager.datatype, filename)

  def get_sequence_names(self):
    return self._get_names(SequenceDataManager.datatype)

  def get_default_sequence_name(self):
    return self._get_default_name(SequenceDataManager.datatype)

  def remove_sequence(self, filename):
    return self._remove(SequenceDataManager.datatype, filename)

  def has_sequences(self, expected_n=1, exact_count=False, raise_sorry=False):
    return self._has_data(SequenceDataManager.datatype, expected_n=expected_n,
                          exact_count=exact_count, raise_sorry=raise_sorry)

  def process_sequence_file(self, filename):
    return self._process_file(SequenceDataManager.datatype, filename)

  def write_sequence_file(self, filename, sequence_str, overwrite=False):
    self._write_text(SequenceDataManager.datatype, filename,
                     sequence_str, overwrite=overwrite)

# =============================================================================
# end
