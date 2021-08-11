'''
Parent DataManager class for handling reciprocal space data
'''
from __future__ import absolute_import, division, print_function
from six.moves import zip

import os
from copy import deepcopy

import iotbx.phil

from cctbx import crystal
from iotbx.data_manager import DataManagerBase
from iotbx.reflection_file_utils import reflection_file_server
from libtbx import Auto
from libtbx.utils import Sorry

# =============================================================================
class MillerArrayDataManager(DataManagerBase):

  datatype = 'miller_array'
  miller_array_child_datatypes = []                     # track children

  # ---------------------------------------------------------------------------
  # Miller arrays

  def add_miller_array_phil_str(self):
    '''
    Add custom PHIL and storage for type and labels
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

  def set_default_miller_array_type(self, array_type=None):
    return self._set_default_miller_array_type(MillerArrayDataManager.datatype,
                                               array_type)

  def get_default_miller_array_type(self):
    return self._get_default_miller_array_type(MillerArrayDataManager.datatype)

  def set_default_miller_array(self, filename):
    return self._set_default(MillerArrayDataManager.datatype, filename)

  def get_miller_array(self, filename=None):
    '''
    Returns the main file object
    '''
    return self._get(MillerArrayDataManager.datatype, filename)

  def set_miller_array_type(self, filename=None, label=None, array_type=None):
    return self._set_miller_array_type(MillerArrayDataManager.datatype,
                                       filename, label, array_type)

  def get_miller_array_type(self, filename=None, label=None):
    return self._get_miller_array_type(MillerArrayDataManager.datatype,
                                       filename, label)

  def get_miller_array_labels(self, filename=None):
    '''
    Returns a list of array labels
    '''
    return self._get_array_labels(MillerArrayDataManager.datatype, filename)

  def get_miller_array_types(self, filename=None):
    '''
    Returns a dict of array types, keyed by label
    '''
    return self._get_array_types(MillerArrayDataManager.datatype, filename)

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
    if filename not in self._miller_array_arrays.keys():
      self._miller_array_arrays[filename] = {}
    if filename not in self._miller_array_types.keys():
      self._miller_array_types[filename] = {}
    miller_arrays = self.get_miller_array(filename).as_miller_arrays()
    labels = []
    for array in miller_arrays:
      label = array.info().label_string()
      labels.append(label)
      self._miller_array_arrays[filename][label] = array
      self._miller_array_types[filename][label] = self._default_miller_array_type
    self._miller_array_labels[filename] = labels

  def write_miller_array_file(self, mtz_object, filename=Auto, overwrite=Auto):
    '''
    Write an MTZ file

    :param mtz_object:        iotbx_mtz_ext.object
    :param filename:          str for the output filename
    :param overwrite:         bool for overwriting files

    The mtz_object can be constructed from a cctbx.miller_array object to
    make a iotbx_mtz_ext.dataset, which can then construct the mtz_object.
    '''
    # default options
    if filename is Auto:
      filename = self.get_default_output_filename() + '.mtz'
    if overwrite is Auto:
      overwrite = self._overwrite

    # check overwrite
    if os.path.isfile(filename) and (not overwrite):
      raise Sorry('%s already exists and overwrite is set to %s.' %
                  (filename, overwrite))

    try:
      mtz_object.write(file_name=filename)
    except IOError as err:
      raise Sorry('There was an error with writing %s.\n%s' %
                  (filename, err))

    self._output_files.append(filename)
    self._output_types.append(MillerArrayDataManager.datatype)

    return filename

  def get_reflection_file_server(self, filenames=None, labels=None,
                                 array_type=None,
                                 crystal_symmetry=None, force_symmetry=None,
                                 logger=None):
    '''
    Return the file server for a single miller_array file or mulitple files

    :param filenames:         list of filenames or None
    :param labels:            list of lists of labels or None
    :param array_type:        "x_ray", "neutron", "electron", or None
    :param crystal_symmetry:  cctbx.crystal.symmetry object or None
    :param force_symmetry:    bool or None
    :param logger:            defaults to self.logger (multi_out)

    The order in filenames and labels should match, e.g. labels[0] should be the
    labels for filenames[0]. The lengths of filenames and labels should be equal
    as well. If all the labels in a file are to be added, set the labels entry
    to None, e.g. labels[0] = None.

    If array_type is None, files of any type are allowed.

    None is returned if the file_server has no arrays
    '''

    # use default logger if no logger for reflection_file_server is provided
    if logger is None:
      logger = self.logger

    # use all miller_array files if no filenames are provided
    if filenames is None:
      filenames = self.get_miller_array_names()

    # set labels
    if labels is None:
      labels = [None for filename in filenames]

    # set labels if the length of labels does not match length of filenames
    assert len(filenames) >= len(labels)
    if len(filenames) > len(labels):
      labels += [None]*(len(filenames) - len(labels))
    assert len(filenames) == len(labels)

    # force crystal symmetry if a crystal symmetry is provided
    if crystal_symmetry is not None and force_symmetry is None:
      force_symmetry = True
    if force_symmetry is None:
      force_symmetry = False

    #  determine crystal symmetry from file list, if not provided
    if crystal_symmetry is None:
      try:
        from_reflection_files = []
        for filename, file_labels in zip(filenames, labels):
          miller_arrays = self.get_miller_arrays(filename=filename)
          for miller_array in miller_arrays:
            if ((file_labels is None) or
                (miller_array.info().label_string() in file_labels)):
              from_reflection_files.append(miller_array.crystal_symmetry())
        crystal_symmetry = crystal.select_crystal_symmetry(
          from_reflection_files=from_reflection_files)
      except Exception:
        raise Sorry('A unit cell and space group could not be determined from the "filenames" argument. Please make sure there are enough data arrays being selected.')

    # crystal_symmetry and force_symmetry should be set by now
    miller_arrays = []
    for filename, file_labels in zip(filenames, labels):
      file_arrays = self.get_miller_array(filename).\
        as_miller_arrays(crystal_symmetry=crystal_symmetry,
                         force_symmetry=force_symmetry)
      if file_labels is None:
        file_labels = self.get_miller_array_labels(filename)
      for miller_array in file_arrays:
        label_name = miller_array.info().label_string()
        # check array label
        if label_name in file_labels:
          # check array type
          if (array_type is None
              or array_type == self.get_miller_array_type(filename, label_name)):
            miller_arrays.append(miller_array)
    file_server = reflection_file_server(
      crystal_symmetry=crystal_symmetry,
      force_symmetry=force_symmetry,
      miller_arrays=miller_arrays,
      err=logger)
    if len(file_server.miller_arrays) == 0:
      file_server = None
    return file_server

  # -----------------------------------------------------------------------------
  # base functions for miller_array datatypes
  def _add_miller_array_phil_str(self, datatype):

    # set up storage
    # self._miller_array_types = {}       # [filename] = type dict
    # self._miller_array_labels = {}      # [filename] = label list
    # self._miller_array_arrays = {}      # [filename] = array dict
    setattr(self, '_%s_types' % datatype, {})
    setattr(self, '_default_%s_type' % datatype, 'x_ray')
    setattr(self, '_possible_%s_types' % datatype,
            ['x_ray', 'neutron', 'electron'])
    setattr(self, '_%s_labels' % datatype, {})
    setattr(self, '_%s_arrays' % datatype, {})

    # custom PHIL section
    custom_phil_str = '''
%s
  .multiple = True
{
  file = None
    .type = path
    .short_caption = MTZ file
    .style = file_type:hkl input_file
  labels
    .multiple = True
  {
    name = None
      .type = str
    type = *%s
      .type = choice(multi=False)
  }
}
''' % (datatype, ' '.join(getattr(self, '_possible_%s_types' % datatype)))

    # custom PHIL scope
    setattr(self, '_custom_%s_phil' % datatype,
            iotbx.phil.parse(custom_phil_str))

    # add to child datatypes
    MillerArrayDataManager.miller_array_child_datatypes.append(datatype)

    return custom_phil_str

  def _export_miller_array_phil_extract(self, datatype):
    extract = []
    filenames = getattr(self, 'get_%s_names' % datatype)()
    for filename in filenames:
      item_extract = getattr(self, '_custom_%s_phil' % datatype).extract()
      item_extract = deepcopy(getattr(item_extract, '%s' % datatype)[0])
      item_extract.file = filename
      labels = self._get_array_labels(datatype, filename=filename)
      types = self._get_array_types(datatype, filename=filename)
      labels_extract = []
      for label in labels:
        label_extract = deepcopy(item_extract.labels[0])
        label_extract.name = label
        label_extract.type = types[label]
        labels_extract.append(label_extract)
      item_extract.labels = labels_extract
      extract.append(item_extract)
    return extract

  def _load_miller_array_phil_extract(self, datatype, phil_extract):
    extract = phil_extract.data_manager
    extract = getattr(extract, '%s' % datatype)
    for item_extract in extract:
      if not hasattr(item_extract, 'file'):
        raise Sorry('This PHIL is not properly defined for the %s datatype.\n There should be a parameter for the filename ("file").\n')

      # process file
      getattr(self, 'process_%s_file' % datatype)(item_extract.file)

      # check labels (if available)
      if len(item_extract.labels) > 0:
        # all labels in file
        file_labels = getattr(self, '_%s_labels' % datatype)[item_extract.file]

        # labels from PHIL
        phil_labels = []
        phil_types = {}
        for label in item_extract.labels:
          if label.name not in file_labels:
            raise Sorry('The label, %s, could not be found in %s.' %
                        (label.name, item_extract.file))
          phil_labels.append(label.name)
          if label.type not in getattr(self, '_possible_%s_types' % datatype):
            raise Sorry('Unrecognized %s type, "%s," possible choices are %s.' %
                        (datatype, label.type, ', '.join(
                          getattr(self, '_possible_%s_types' % datatype))))
          phil_types[label.name] = label.type
        getattr(self, '_%s_labels' % datatype)[item_extract.file] = phil_labels
        getattr(self, '_%s_types' % datatype)[item_extract.file] = phil_types

  def _set_default_miller_array_type(self, datatype, array_type):
    if array_type not in getattr(self, '_possible_%s_types' % datatype):
      raise Sorry('Unrecognized %s type, "%s," possible choices are %s.' %
                  (datatype, array_type,
                   ', '.join(getattr(self, '_possible_%s_types' % datatype))))
    setattr(self, '_default_%s_type' % datatype, array_type)

  def _get_default_miller_array_type(self, datatype):
    return getattr(self, '_default_%s_type' % datatype)

  def _set_miller_array_type(self, datatype, filename=None, label=None,
                             array_type=None):
    if filename is None:
      filename = self._get_default_name(datatype)
    if label is None:
      label = self._get_array_labels(datatype, filename)[0]
    if array_type is None:
      array_type = getattr(self, '_default_%s_type' % datatype)
    elif array_type not in getattr(self, '_possible_%s_types' % datatype):
      raise Sorry('Unrecognized %s type, "%s," possible choices are %s.' %
                  (datatype,
                   array_type,
                   ', '.join(getattr(self, '_possible_%s_types' % datatype))))
    getattr(self, '_%s_types' % datatype)[filename][label] = array_type

  def _get_miller_array_type(self, datatype, filename=None, label=None):
    if filename is None:
      filename = self._get_default_name(datatype)
    if label is None:
      label = self._get_array_labels(datatype, filename)[0]
    types = self._get_array_types(datatype, filename)
    return types.get(label, getattr(self, '_default_%s_type' % datatype))

  def _filter_miller_array_child_datatypes(self, filename):
    # filter arrays (e.g self.filter_map_coefficients_arrays)
    for datatype in MillerArrayDataManager.miller_array_child_datatypes:
      function_name = 'filter_%s_arrays' % datatype
      if hasattr(self, function_name):
        getattr(self, function_name)(filename)

  def _check_miller_array_default_filename(self, datatype, filename=None):
    '''
    Helper function for handling default filenames
    '''
    if filename is None:
      filename = getattr(self, 'get_default_%s_name' % datatype)()
      if filename is None:
        raise Sorry('There is no default filename for %s.' % datatype)
    return filename

  def _check_miller_array_storage_dict(self, datatype, storage_dict, filename):
    '''
    Helper function for checking filename and datatypes
    '''
    if filename not in storage_dict.keys():
      self.process_miller_array_file(filename)
      if filename not in storage_dict.keys():
        raise Sorry('There are no known %s arrays in %s' % (datatype, filename))

  def _get_array_labels(self, datatype, filename=None):
    filename = self._check_miller_array_default_filename(datatype, filename)
    storage_dict = getattr(self, '_%s_labels' % datatype)
    self._check_miller_array_storage_dict(datatype, storage_dict, filename)
    labels = storage_dict[filename]
    return labels

  def _get_array_types(self, datatype, filename=None):
    filename = self._check_miller_array_default_filename(datatype, filename)
    storage_dict = getattr(self, '_%s_types' % datatype)
    self._check_miller_array_storage_dict(datatype, storage_dict, filename)
    types = storage_dict[filename]
    return types

  def _get_arrays(self, datatype, filename=None, labels=None):
    filename = self._check_miller_array_default_filename(datatype, filename)
    if filename not in getattr(self, 'get_%s_names' % datatype)():
      self.process_miller_array_file(filename)
    if labels is None:
      labels = self._get_array_labels(datatype, filename=filename)
    else:
      if not isinstance(labels, list):
        raise Sorry('The labels argument should be a list of labels')
    storage_dict = getattr(self, '_%s_arrays' % datatype)
    self._check_miller_array_storage_dict(datatype, storage_dict, filename)
    arrays = []
    for label in labels:
      arrays.append(self._get_array_by_label(datatype, filename, label))
    return arrays

  def _get_array_by_label(self, datatype, filename, label):
    storage_dict = getattr(self, '_%s_arrays' % datatype)
    if label in storage_dict[filename].keys():
      return storage_dict[filename][label]
    else:
      raise Sorry('%s does not have any arrays labeled %s' %
                  (filename, label))

  def _child_filter_arrays(self, datatype, filename, known_labels):
    '''
    Populate data structures by checking labels in miller arrays to determine
    child type
    '''
    if datatype not in MillerArrayDataManager.miller_array_child_datatypes:
      raise Sorry('%s is not a valid child datatype of miller_array.')
    data = self.get_miller_array(filename)
    miller_arrays = data.as_miller_arrays()
    labels = []
    types = {}
    datatype_dict = getattr(self, '_%s_arrays' % datatype)
    for array in miller_arrays:
      label = set(array.info().labels)
      common_labels = known_labels.intersection(label)
      if len(common_labels) > 0:
        label = array.info().label_string()
        labels.append(label)
        if filename not in datatype_dict.keys():
          datatype_dict[filename] = {}
        datatype_dict[filename][label] = array
        types[label] = getattr(self, '_default_%s_type' % datatype)

    # if arrays exist, start tracking
    if len(labels) > 1:
      getattr(self, '_%s_labels' % datatype)[filename] = labels
      getattr(self, '_%s_types' % datatype)[filename] = types
      self._add(datatype, filename, data)

# =============================================================================
# end
