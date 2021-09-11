from __future__ import absolute_import, division, print_function
'''
'''

from iotbx.data_manager import DataManagerBase
from libtbx import Auto

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
    '''
    Returns a list of sequence objects from a file
    '''
    return self._get(SequenceDataManager.datatype, filename)

  def get_sequence_as_string(self, filename=None, width=80,
      sequence_letters_only = False):
    '''
    Same as get_sequence, but returns a single string with all the sequences
    separated by blank lines. If sequence_letters_only,remove
    all lines starting with ">" removed and all spaces and line breaks
    '''
    sequences = self.get_sequence(filename=filename)
    if sequence_letters_only:
      output = '\n\n'.join([s.sequence for s in sequences])
    else: # usual
      output = '\n\n'.join([s.format(width) for s in sequences])
    return output

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

  def get_default_output_sequence_filename(self):
    filename = self.get_default_output_filename()
    if not filename.endswith('.seq'):
      filename += '.seq'
    return filename

  def write_sequence_file(self, sequence_str, filename=Auto, overwrite=Auto):
    if filename is Auto:
      filename = self.get_default_output_sequence_filename()
    return self._write_text(SequenceDataManager.datatype, sequence_str,
                            filename=filename, overwrite=overwrite)

# =============================================================================
# end
