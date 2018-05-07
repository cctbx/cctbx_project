from __future__ import division, print_function
'''
'''

import iotbx.phil

from iotbx.data_manager import DataManagerBase
from libtbx.utils import Sorry

# =============================================================================
class MillerArrayDataManager(DataManagerBase):

  datatype = 'miller_array'
  miller_array_child_datatypes = list()                     # track children

  # ---------------------------------------------------------------------------
  # Miller arrays

  def add_miller_array_phil_str(self):
    '''
    Add custom PHIL and storage for labels
    '''
    return self._add_miller_array_phil_str(MillerArrayDataManager.datatype)

  def export_miller_array_phil_extract(self):
    '''
    Export custom PHIL extract
    '''
    return self._export_miller_array_phil_extract(
      MillerArrayDataManager.datatype)

  def load_miller_array_phil_extract(self, phil_extract):
    '''
    Load custom PHIL extract
    '''
    self._load_miller_array_phil_extract(MillerArrayDataManager.datatype,
                                         phil_extract)

  def add_miller_array(self, filename, data):
    self._add(MillerArrayDataManager.datatype, filename, data)
    self._filter_miller_array_child_datatypes(filename)

  def set_default_miller_array(self, filename):
    return self._set_default(MillerArrayDataManager.datatype, filename)

  def get_miller_array(self, filename=None):
    '''
    Returns the main file object
    '''
    return self._get(MillerArrayDataManager.datatype, filename)

  def get_miller_array_labels(self, filename=None):
    '''
    Returns a list of array labels
    '''
    return self._get_array_labels(MillerArrayDataManager.datatype, filename)

  def get_miller_arrays(self, labels=None, filename=None):
    '''
    Returns a list of arrays
    '''
    return self._get_arrays(MillerArrayDataManager.datatype, labels=labels,
                            filename=filename)

  def get_miller_array_names(self):
    return self._get_names(MillerArrayDataManager.datatype)

  def get_default_miller_array_name(self):
    return self._get_default_name(MillerArrayDataManager.datatype)

  def remove_miller_array(self, filename):
    return self._remove(MillerArrayDataManager.datatype, filename)

  def has_miller_arrays(self, expected_n=1, exact_count=False, raise_sorry=False):
    return self._has_data(MillerArrayDataManager.datatype, expected_n=expected_n,
                          exact_count=exact_count, raise_sorry=raise_sorry)

  def process_miller_array_file(self, filename):
    self._process_file(MillerArrayDataManager.datatype, filename)
    self._filter_miller_array_child_datatypes(filename)

  def filter_miller_array_arrays(self, filename):
    '''
    Populate data structures with all arrays
    '''
    if (filename not in self._miller_array_arrays.keys()):
      self._miller_array_arrays[filename] = dict()
    miller_arrays = self.get_miller_array(filename).as_miller_arrays()
    labels = list()
    for array in miller_arrays:
      label = ', '.join(array.info().labels)
      labels.append(label)
      self._miller_array_arrays[filename][label] = array
    self._miller_array_labels[filename] = labels

  def write_miller_array_file(self, filename, miller_arrays, overwrite=False):
    raise NotImplementedError

  # -----------------------------------------------------------------------------
  # base functions for miller_array datatypes
  def _add_miller_array_phil_str(self, datatype):

    # set up storage
    # self._miller_array_labels = dict()      # [filename] = label list
    # self._miller_array_arrays = dict()      # [filename] = array dict
    setattr(self, '_%s_labels' % datatype, dict())
    setattr(self, '_%s_arrays' % datatype, dict())

    # custom PHIL section
    custom_phil_str = '%s\n.multiple = True\n{\n' % datatype
    custom_phil_str += 'file = None\n'
    custom_phil_str += '.type = path\n'
    custom_phil_str += 'labels = None\n'
    custom_phil_str += '.type = str\n.multiple = True\n'
    custom_phil_str += '}\n'

    # custom PHIL scope
    setattr(self, '_custom_%s_phil' % datatype,
            iotbx.phil.parse(custom_phil_str))

    # add to child datatypes
    MillerArrayDataManager.miller_array_child_datatypes.append(datatype)

    return custom_phil_str

  def _export_miller_array_phil_extract(self, datatype):
    extract = list()
    filenames = getattr(self, 'get_%s_names' % datatype)()
    for filename in filenames:
      item_extract = getattr(self, '_custom_%s_phil' % datatype).extract()
      item_extract = getattr(item_extract, '%s' % datatype)[0]
      item_extract.file = filename
      item_extract.labels = self._get_array_labels(datatype, filename=filename)
      extract.append(item_extract)
    return extract

  def _load_miller_array_phil_extract(self, datatype, phil_extract):
    extract = phil_extract.data_manager
    extract = getattr(extract, '%s' % datatype)
    for item_extract in extract:
      if ( (not hasattr(item_extract, 'file')) or
           (not hasattr(item_extract, 'labels')) ):
        raise Sorry('This PHIL is not properly defined for the %s datatype.\n There should be a parameter for the filename ("file") and labels ("labels").\n')

      # process file
      getattr(self, 'process_%s_file' % datatype)(item_extract.file)
      # update labels
      getattr(self, '_%s_labels' % datatype)[item_extract.file] = \
        item_extract.labels

  def _filter_miller_array_child_datatypes(self, filename):
    # filter arrays (e.g self.filter_map_coefficients_arrays)
    for datatype in MillerArrayDataManager.miller_array_child_datatypes:
      function_name = 'filter_%s_arrays' % datatype
      if (hasattr(self, function_name)):
        getattr(self, function_name)(filename)

  def _check_miller_array_default_filename(self, datatype, filename=None):
    '''
    Helper function for handling default filenames
    '''
    if (filename is None):
      filename = getattr(self, 'get_default_%s_name' % datatype)()
      if (filename is None):
        raise Sorry('There is no default filename for %s.' % datatype)
    return filename

  def _check_miller_array_storage_dict(self, datatype, storage_dict, filename):
    '''
    Helper function for checking filename and datatypes
    '''
    if (filename not in storage_dict.keys()):
      self.process_miller_array_file(filename)
      if (filename not in storage_dict.keys()):
        raise Sorry('There are no known %s arrays in %s' % (datatype, filename))

  def _get_array_labels(self, datatype, filename=None):
    filename = self._check_miller_array_default_filename(datatype, filename)
    storage_dict = getattr(self, '_%s_labels' % datatype)
    self._check_miller_array_storage_dict(datatype, storage_dict, filename)
    labels = storage_dict[filename]

    return labels

  def _get_arrays(self, datatype, labels=None, filename=None):
    filename = self._check_miller_array_default_filename(datatype, filename)
    if (filename not in getattr(self, 'get_%s_names' % datatype)()):
      self.process_miller_array_file(filename)
    if (labels is None):
      labels = self._get_array_labels(datatype, filename=filename)
    else:
      if (not isinstance(labels, list)):
        raise Sorry('The labels object should be a list of labels')

    storage_dict = getattr(self, '_%s_arrays' % datatype)
    self._check_miller_array_storage_dict(datatype, storage_dict, filename)
    arrays = list()
    for label in labels:
      arrays.append(self._get_array_by_label(datatype, label, filename))
    return arrays

  def _get_array_by_label(self, datatype, label, filename):
    storage_dict = getattr(self, '_%s_arrays' % datatype)
    if (label in storage_dict[filename].keys()):
      return storage_dict[filename][label]
    else:
      raise Sorry('%s does not have any arrays labeled %s' %
                  (filename, label))

  def _child_filter_arrays(self, datatype, filename, known_labels):
    '''
    Populate data structures by checking labels in miller arrays to determine
    child type
    '''
    if (datatype not in MillerArrayDataManager.miller_array_child_datatypes):
      raise Sorry('%s is not a valid child datatype of miller_array.')
    data = self.get_miller_array(filename)
    miller_arrays = data.as_miller_arrays()
    labels = list()
    for array in miller_arrays:
      label = set(array.info().labels)
      common_labels = known_labels.intersection(label)
      if (len(common_labels) > 0):
        label = ', '.join(array.info().labels)
        labels.append(label)
        datatype_dict = getattr(self, '_%s_arrays' % datatype)
        if (filename not in datatype_dict.keys()):
          datatype_dict[filename] = dict()
        datatype_dict[filename][label] = array

    # if arrays exist, start tracking
    if (len(labels) > 1):
      getattr(self, '_%s_labels' % datatype)[filename] = labels
      self._add(datatype, filename, data)

# =============================================================================
# end
