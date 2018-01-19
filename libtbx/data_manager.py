from __future__ import division, print_function
'''
'''

from libtbx.utils import Sorry

# =============================================================================
class DataManager(object):

  def __init__(self):
    '''
    '''
    self._datatypes = ['model', 'sequence']
    self._storage = '_%ss'
    self._default = '_default_%s'
    self._current_storage = None
    self._current_default = None

    # generate storage for each data type
    for datatype in self._datatypes:
      # data attributes for each data type
      # e.g self._models, self._default_model
      self.set_datatype(datatype)
      setattr(self, self._current_storage, dict())
      setattr(self, self._current_default, None)

  # ---------------------------------------------------------------------------
  def set_datatype(self, datatype):
    '''
    '''
    if (datatype not in self._datatypes):
      msg = '"%s" is not a recognized datatype. Only\n %s\n are recognized'
      raise Sorry( msg % (datatype, '\n  '.join(self._datatypes)))
    self._current_storage = self._storage % datatype
    self._current_default = self._default % datatype

  # ---------------------------------------------------------------------------
  # Models
  def add_model(self, filename, data):
    return self._add('model', filename, data)

  def set_default_model(self, filename):
    return self._set_default('model', filename)

  def get_model(self, filename=None):
    return self._get('model', filename)

  def get_model_names(self):
    return self._get_names('model')

  def get_default_model_name(self):
    return self._get_default_name('model')

  def remove_model(self, filename):
    return self._remove('model', filename)

  def has_models(self, expected_n=1, exact_count=False):
    return self._has_data('model', expected_n=expected_n,
                          exact_count=exact_count)

  # ---------------------------------------------------------------------------
  # Sequences
  def add_sequence(self, filename, data):
    return self._add('sequence', filename, data)

  def set_default_sequence(self, filename):
    return self._set_default('sequence', filename)

  def get_sequence(self, filename=None):
    return self._get('sequence', filename)

  def get_sequence_names(self):
    return self._get_names('sequence')

  def get_default_sequence_name(self):
    return self._get_default_name('sequence')

  def remove_sequence(self, filename):
    return self._remove('sequence', filename)

  def has_sequences(self, expected_n=1, exact_count=False):
    return self._has_data('sequence', expected_n=expected_n,
                          exact_count=exact_count)

  # ---------------------------------------------------------------------------
  # Generic functions for manipulating data
  def _get_current_storage(self):
    '''
    '''
    return getattr(self, self._current_storage)

  def _get_current_default(self):
    '''
    '''
    return getattr(self, self._current_default)

  def _add(self, datatype, filename, data):
    '''
    '''
    self.set_datatype(datatype)
    self._get_current_storage()[filename] = data
    if (self._get_current_default() is None):
      self._set_default(datatype, filename)

  def _set_default(self, datatype, filename):
    '''
    '''
    self.set_datatype(datatype)
    if (filename not in self._get_current_storage().keys()):
      raise Sorry('"%s" is not a known %s type.' % (filename, datatype))
    else:
      setattr(self, self._current_default, filename)

  def _get(self, datatype, filename=None):
    self.set_datatype(datatype)
    if (filename is None):
      default_filename = self._get_current_default()
      if (default_filename is None):
        raise Sorry('No default %s defined' % datatype)
      else:
        return self._get(datatype, filename=default_filename)
    elif (filename not in self._get_current_storage().keys()):
      raise Sorry('"%s" is not a known %s type.' % (filename, datatype))
    else:
      return self._get_current_storage()[filename]

  def _get_names(self, datatype):
    self.set_datatype(datatype)
    return self._get_current_storage().keys()

  def _get_default_name(self, datatype):
    self.set_datatype(datatype)
    return self._get_current_default()

  def _remove(self, datatype, filename):
    self.set_datatype(datatype)
    if (filename not in self._get_current_storage().keys()):
      raise Sorry('"%s" is not being managed.' % filename)
    else:
      self._get_current_storage().pop(filename)
      if (filename == self._get_current_default()):
        setattr(self, self._current_default, None)

  def _has_data(self, datatype, expected_n=1, exact_count=False):
    self.set_datatype(datatype)
    actual_n = len(self._get_names(datatype))
    v = cmp(actual_n, expected_n)
    if (exact_count):
      # exact count required
      if (v != 0):
        raise Sorry('%i %s(s) found. Expected exactly %i.' %
                    (actual_n, datatype, expected_n))
      else:
        return True
    else:
      return (v >= 0)

# =============================================================================
