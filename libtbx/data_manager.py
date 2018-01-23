from __future__ import division, print_function
'''
'''

import iotbx.pdb
import libtbx.phil
import mmtbx.model

from iotbx.file_reader import any_file
from libtbx.utils import Sorry

# mapping from DataManager datatypes to any_file file types
any_file_type = {
  'model':'pdb',
  'sequence':'seq',
  'phil':'phil'
}

# =============================================================================
class DataManager(object):

  # PHIL scope for DataManager
  master_phil_str = '''
  data_manager {

    model_files = None
      .type = path
      .multiple = True
    default_model = None
      .type = path

    sequence_files = None
      .type = path
      .multiple = True
    default_sequence = None
      .type = path

    phil_files = None
      .type = path
      .multiple = True
    default_phil = None
      .type = path

  }
  '''

  _datatypes = ['model', 'sequence', 'phil']

  def __init__(self, phil=None):
    '''
    '''

    self.master_phil = libtbx.phil.parse(DataManager.master_phil_str)

    self._storage = '_%ss'
    self._default = '_default_%s'
    self._current_storage = None
    self._current_default = None

    # generate storage for each data type
    for datatype in DataManager._datatypes:
      # data attributes for each data type
      # e.g self._models, self._default_model
      self.set_datatype(datatype)
      setattr(self, self._current_storage, dict())
      setattr(self, self._current_default, None)

    # load information from phil
    if (phil is not None):
      self.load_phil_scope(phil)

  # ---------------------------------------------------------------------------
  def export_phil_scope(self):
    '''
    Function for exporting DataManager information into a PHIL scope
    The returned PHIL scope can be used to recreate the DataManager object with
    the load_phil_scope function

    This assumes that the key names in the data structures are valid filenames.
    '''
    phil_extract = self.master_phil.extract()
    for datatype in DataManager._datatypes:
      filenames = self._get_names(datatype)
      default = self._get_default_name(datatype)
      setattr(phil_extract.data_manager, '%s_files' % datatype, filenames)
      setattr(phil_extract.data_manager, 'default_%s' % datatype, default)

    working_phil = self.master_phil.format(python_object=phil_extract)

    return working_phil

  # ---------------------------------------------------------------------------
  def load_phil_scope(self, phil):
    '''
    Function for loading information from a PHIL scope. This will append files
    to the existing ones and will NOT override the current default file
    '''
    # sanity checks
    if (type(phil) == libtbx.phil.scope):
      working_phil = self.master_phil.fetch(source=phil)
    elif (type(phil) == libtbx.phil.scope_extract):
      working_phil = self.master_phil.format(python_object=phil)
    else:
      raise Sorry('A libtbx.phil.scope or libtbx.phil.scope_extract object is required')
    phil_extract = working_phil.extract()

    if (not hasattr(phil_extract, 'data_manager')):
      raise Sorry('The phil scope does not have a DataManager scope.')

    # append files, default file is not overridden
    for datatype in DataManager._datatypes:
      filenames = getattr(phil_extract.data_manager, '%s_files' % datatype,
                          None)
      for filename in filenames:
        # call type-specific function (e.g. self.process_model())
        # checks if file is already in DataManager
        getattr(self, 'process_%s_file' % datatype)(filename)

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

  def has_models(self, expected_n=1, exact_count=False, raise_sorry=True):
    return self._has_data('model', expected_n=expected_n,
                          exact_count=exact_count, raise_sorry=raise_sorry)

  def process_model_file(self, filename):
    # unique because any_file does not return a model object
    if (filename not in self.get_model_names()):
      a = any_file(filename)
      if (a.file_type != 'pdb'):
        raise Sorry('%s is not a recognized model file' % filename)
      else:
        model_in = iotbx.pdb.input(a.file_name)
        model = mmtbx.model.manager(model_input=model_in)
        self.add_model(filename, model)

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

  def has_sequences(self, expected_n=1, exact_count=False, raise_sorry=True):
    return self._has_data('sequence', expected_n=expected_n,
                          exact_count=exact_count, raise_sorry=raise_sorry)

  def process_sequence_file(self, filename):
    return self._process_file('sequence', filename)

  # ---------------------------------------------------------------------------
  # PHIL
  def add_phil(self, filename, data):
    return self._add('phil', filename, data)

  def set_default_phil(self, filename):
    return self._set_default('phil', filename)

  def get_phil(self, filename=None):
    return self._get('phil', filename)

  def get_phil_names(self):
    return self._get_names('phil')

  def get_default_phil_name(self):
    return self._get_default_name('phil')

  def remove_phil(self, filename):
    return self._remove('phil', filename)

  def has_phils(self, expected_n=1, exact_count=False, raise_sorry=True):
    return self._has_data('phil', expected_n=expected_n,
                          exact_count=exact_count, raise_sorry=raise_sorry)

  def process_phil_file(self, filename):
    return self._process_file('phil', filename)

  # ---------------------------------------------------------------------------
  def set_datatype(self, datatype):
    '''
    '''
    if (datatype not in DataManager._datatypes):
      msg = '"%s" is not a recognized datatype. Only\n %s\n are recognized'
      raise Sorry( msg % (datatype, '\n  '.join(DataManager._datatypes)))
    self._current_storage = self._storage % datatype
    self._current_default = self._default % datatype

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

  def _has_data(self, datatype, expected_n=1, exact_count=False,
                raise_sorry=True):
    self.set_datatype(datatype)
    actual_n = len(self._get_names(datatype))
    v = cmp(actual_n, expected_n)
    if (exact_count):
      # exact count required
      if (v != 0):
        if (raise_sorry):
          raise Sorry('%i %s(s) found. Expected exactly %i.' %
                      (actual_n, datatype, expected_n))
        else:
          return False
      else:
        return True
    else:
      if (raise_sorry):
        if (v < 0):
          raise Sorry('%i %s(s) found. Expected at least %i.' %
                      (actual_n, datatype, expected_n))
      return (v >= 0)

  def _process_file(self, datatype, filename):
    if (filename not in self._get_names(datatype)):
      a = any_file(filename)
      if (a.file_type != any_file_type[datatype]):
        raise Sorry('%s is not a recognized %s file' % (filename, datatype))
      else:
        self._add(datatype, filename, a.file_object)

# =============================================================================
# end
