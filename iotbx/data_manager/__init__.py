from __future__ import division, print_function
from six.moves import range
'''
'''

import inspect, os, re, sys
import importlib

from collections import OrderedDict

import libtbx.phil

from iotbx.file_reader import any_file
from libtbx import Auto
from libtbx.utils import multi_out, Sorry

# =============================================================================
# mapping from DataManager datatypes to any_file file types
any_file_type = {
  'map_coefficients':'hkl',
  'miller_array':'hkl',
  'model':'pdb',
  'ncs_spec':'ncs',
  'phil':'phil',
  'real_map':'ccp4_map',
  'restraint':'cif',
  'sequence':'seq',
}

# reverse dictionary to map any_file types to DataManager datatypes
data_manager_type = {value:key for key,value in any_file_type.items()}
data_manager_type['hkl'] = 'miller_array'   # map hkl to parent, miller_array

# build list of supported datatypes
# datatypes have corresponding modules in iotbx/data_manager
# e.g. iotbx/data_manager/model.py
supported_datatypes = os.listdir(os.path.dirname(__file__))
re_search = re.compile('.py$')
supported_datatypes = filter(re_search.search, supported_datatypes)
supported_datatypes.remove('__init__.py')
supported_datatypes.sort()
for i in range(len(supported_datatypes)):
  supported_datatypes[i] = supported_datatypes[i].split('.')[0]

default_datatypes = ['miller_array', 'model', 'phil', 'real_map', 'restraint',
                     'sequence']

# =============================================================================
def load_datatype_modules(datatypes=None):
  '''
  Function for dynamically loading a subset of modules for the DataManager
  The default directory for modules is iotbx/data_manager
  '''

  # set default if necessary
  if (datatypes is None):
    datatypes = default_datatypes

  # check datatypes
  error_message_1 = '%s is not a supported datatype'
  error_message_2 = 'Please provide a list of supported datatypes (%s)' % \
                    ','.join(supported_datatypes)
  if (not (isinstance(datatypes, list) or isinstance(datatypes, tuple))):
    if (isinstance(datatypes, str)):
      if (datatypes in supported_datatypes):
        datatypes = [datatypes]
      else:
        raise Sorry(error_message_1 % datatypes)
    else:
      raise Sorry(error_message_2)
  for datatype in datatypes:
    if (isinstance(datatype, str)):
      if (datatype not in supported_datatypes):
        raise Sorry(error_message_1 % datatype)
    else:
      raise Sorry(error_message_2)

  # load modules
  modules = list()
  importlib.import_module('iotbx')
  for datatype in datatypes:
    modules.append(
      importlib.import_module('.' + datatype, package='iotbx.data_manager'))

  return modules

# =============================================================================
def DataManager(datatypes=None, phil=None, logger=None):
  '''
  Function for dynamically creating a DataManager instance that supports a
  specific set of datatypes.

  All DataManager modules in iotbx/data_manager follow this pattern:
    filename = <datatype>.py -> module name = iotbx.data_manager.<datatype>
    DataManagerBase subclass name = <Datatype>DataManager

  So for the models, the filename is iotbx/data_manager/model.py and in that
  file, there should be a subclass of DataManagerBase named ModelDataManager
  '''

  if (logger is None):
    logger = multi_out()

  # set default if necessary
  if (datatypes is None):
    datatypes = default_datatypes

  # get classes
  modules = load_datatype_modules(datatypes)
  manager_classes = list()
  for datatype in datatypes:
    module_name = 'iotbx.data_manager.' + datatype
    class_name =  datatype.capitalize() + 'DataManager'
    if ('_' in datatype): # real_map becomes RealMapDataManager
      class_name = ''.join(tmp_str.capitalize()
                           for tmp_str in datatype.split('_')) + 'DataManager'
    manager_classes.append(getattr(sys.modules[module_name], class_name))

  # check inheritance and add datatypes if necessary
  class_datatypes = set()
  for manager_class in manager_classes:
    if (hasattr(manager_class, 'datatype')):
      class_datatypes.add(manager_class.datatype)
    # get full inheritance order and check
    for parent_class in inspect.getmro(manager_class)[1:]:
      if (hasattr(parent_class, 'datatype')):
        class_datatypes.add(parent_class.datatype)
  datatypes = list(class_datatypes)

  # construct new class and return instance
  data_manager_class = type('DataManager', tuple(manager_classes), dict())
  return data_manager_class(datatypes=datatypes, phil=phil, logger=logger)

# =============================================================================
class DataManagerBase(object):

  def __init__(self, datatypes=list(), phil=None, logger=None):
    '''
    '''

    self.datatypes = datatypes
    if (hasattr(self,'datatype') and (self.datatype not in self.datatypes)):
      self.datatypes.append(self.datatype)

    # functions for custom PHIL
    self.add_custom_phil_str = 'add_%s_phil_str'
    self.export_custom_phil_extract = 'export_%s_phil_extract'
    self.load_custom_phil_extract = 'load_%s_phil_extract'

    # dynamically construct master PHIL string
    self.master_phil_str = 'data_manager {\n'
    self.master_phil_str += """
default_output_filename = None
  .type = str
  .help = The default filename for output files (without file extension).
overwrite = False
  .type = bool
  .help = The default setting for overwriting files with the same name.
"""
    for datatype in self.datatypes:

      # check if a datatype has a custom PHIL str
      if (hasattr(self, self.add_custom_phil_str % datatype)):
        custom_phil_str = getattr(self, self.add_custom_phil_str % datatype)()
        self.master_phil_str += custom_phil_str

      # default PHIL
      else:
        # sequence_files = None
        #   .type = path
        #   .multiple = True
        self.master_phil_str += '%s_files = None\n' % datatype
        self.master_phil_str += '.type = path\n.multiple=True\n'

        # property for wx GUI (will be removed)
        file_type = any_file_type.get(datatype, None)
        if (file_type is not None):
          self.master_phil_str += '.style = file_type:%s input_file\n' % \
                                            file_type

      # default_sequence = None
      #   .type = path
      self.master_phil_str += 'default_%s = None\n' % datatype
      self.master_phil_str += '.type = path\n'
    self.master_phil_str += '}'

    self.master_phil = libtbx.phil.parse(self.master_phil_str)

    self._storage = '_%ss'
    self._default = '_default_%s'
    self._current_storage = None
    self._current_default = None

    # generate storage for each data type
    for datatype in self.datatypes:
      # data attributes for each data type
      # e.g self._models, self._default_model
      self._set_datatype(datatype)
      setattr(self, self._current_storage, OrderedDict())
      setattr(self, self._current_default, None)

    # track output files for internal use
    self._output_files = list()
    self._output_types = list()

    # set defaults
    self._default_output_filename = 'cctbx_program'
    self._overwrite = False

    # load information from phil
    if (phil is not None):
      self.load_phil_scope(phil)

    # logger (currently used for models)
    self.logger = logger
    if (self.logger is None):
      self.logger = multi_out()

  # ---------------------------------------------------------------------------
  def export_phil_scope(self):
    '''
    Function for exporting DataManager information into a PHIL scope
    The returned PHIL scope can be used to recreate the DataManager object with
    the load_phil_scope function

    This assumes that the key names in the data structures are valid filenames.
    '''
    phil_extract = self.master_phil.extract()
    phil_extract.data_manager.default_output_filename = self._default_output_filename
    phil_extract.data_manager.overwrite = self._overwrite
    for datatype in self.datatypes:
      if (hasattr(self, self.export_custom_phil_extract % datatype)):
        setattr(phil_extract.data_manager, '%s' % datatype,
                getattr(self, self.export_custom_phil_extract % datatype)())
      else:
        filenames = self._get_names(datatype)
        setattr(phil_extract.data_manager, '%s_files' % datatype, filenames)

      default = self._get_default_name(datatype)
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
    if (isinstance(phil, libtbx.phil.scope)):
      working_phil = self.master_phil.fetch(source=phil)
    elif (isinstance(phil, libtbx.phil.scope_extract)):
      working_phil = self.master_phil.format(python_object=phil)
    else:
      raise Sorry('A libtbx.phil.scope or libtbx.phil.scope_extract object is required')
    phil_extract = working_phil.extract()

    if (not hasattr(phil_extract, 'data_manager')):
      raise Sorry('The phil scope does not have a DataManager scope.')

    # append files, default file is not overridden
    for datatype in self.datatypes:
      if (hasattr(self, self.load_custom_phil_extract % datatype)):
        getattr(self, self.load_custom_phil_extract % datatype)(phil_extract)
      else:
        filenames = getattr(phil_extract.data_manager, '%s_files' % datatype,
                            None)
        for filename in filenames:
          # call type-specific function (e.g. self.process_model_file())
          # checks if file is already in DataManager
          getattr(self, 'process_%s_file' % datatype)(filename)

    # other options
    self._default_output_filename = getattr(
      phil_extract.data_manager, 'default_output_filename', None)
    self._overwrite = getattr(
      phil_extract.data_manager, 'overwrite', False)

  # ---------------------------------------------------------------------------
  def supports(self, datatype):
    '''
    Function that checks if the DataManager instance supports a particular
    datatype
    '''
    return (datatype in self.datatypes)

  # ---------------------------------------------------------------------------
  def set_default_output_filename(self, filename):
    self._default_output_filename = filename

  def get_default_output_filename(self):
    return self._default_output_filename

  # ---------------------------------------------------------------------------
  def set_overwrite(self, overwrite):
    self._overwrite = overwrite

  def get_overwrite(self, overwrite):
    return self._overwrite

  # ---------------------------------------------------------------------------
  # Generic functions for manipulating data
  def _set_datatype(self, datatype):
    '''
    '''
    if (datatype not in self.datatypes):
      msg = '"%s" is not a recognized datatype. Only\n %s\n are recognized'
      raise Sorry( msg % (datatype, '\n  '.join(self.datatypes)))
    self._current_storage = self._storage % datatype
    self._current_default = self._default % datatype

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
    self._set_datatype(datatype)
    self._get_current_storage()[filename] = data
    if (self._get_current_default() is None):
      self._set_default(datatype, filename)

  def _set_default(self, datatype, filename):
    '''
    '''
    self._set_datatype(datatype)
    if (filename not in self._get_current_storage().keys()):
      raise Sorry('"%s" is not a known %s type.' % (filename, datatype))
    else:
      setattr(self, self._current_default, filename)

  def _get(self, datatype, filename=None):
    self._set_datatype(datatype)
    if (filename is None):
      default_filename = self._get_current_default()
      if (default_filename is None):
        raise Sorry('No default %s defined' % datatype)
      else:
        return self._get(datatype, filename=default_filename)
    elif (filename not in self._get_current_storage().keys()):
      try:
        # try to load file if not already available
        # use type-specific function call instead of _process_file because
        # process_model_file is unique
        # change to _process_file after any_file is updated
        getattr(self, 'process_%s_file' % datatype)(filename)
        return self._get_current_storage()[filename]
      except Sorry:
        raise Sorry('"%s" is not a known %s type.' % (filename, datatype))
    else:
      return self._get_current_storage()[filename]

  def _get_names(self, datatype):
    self._set_datatype(datatype)
    return self._get_current_storage().keys()

  def _get_default_name(self, datatype):
    self._set_datatype(datatype)
    return self._get_current_default()

  def _remove(self, datatype, filename):
    self._set_datatype(datatype)
    if (filename not in self._get_current_storage().keys()):
      raise Sorry('"%s" is not being managed.' % filename)
    else:
      self._get_current_storage().pop(filename)
      if (filename == self._get_current_default()):
        setattr(self, self._current_default, None)

  def _has_data(self, datatype, expected_n=1, exact_count=False,
                raise_sorry=False):
    self._set_datatype(datatype)
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

  def _write_text(self, datatype, text_str, filename=Auto, overwrite=Auto):
    '''
    Convenience function for writing text to file
    '''

    # default options
    if (filename is Auto):
      filename = self._default_output_filename
    if (overwrite is Auto):
      overwrite = self._overwrite

    # check arguments
    if (os.path.isfile(filename) and (not overwrite)):
      raise Sorry('%s already exists and overwrite is set to %s.' %
                  (filename, overwrite))
    if (not isinstance(text_str, str)):
      raise Sorry('Please provide a text string for writing.')

    try:
      with open(filename, 'w') as f:
        f.write(text_str)
    except IOError as err:
      raise Sorry('There was an error with writing %s.\n%s' %
                  (filename, err))

    self._output_files.append(filename)
    self._output_types.append(datatype)

# =============================================================================
# end
