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
from iotbx.reflection_file_utils import label_table, get_experimental_phases_scores, \
  get_phase_scores, get_r_free_flags_scores, reflection_file_server
from libtbx import Auto
from libtbx.utils import Sorry

# =============================================================================
class MillerArrayDataManager(DataManagerBase):

  datatype = 'miller_array'
  miller_array_child_datatypes = []                     # track children

  # template storage names
  _type_str = '_%s_types'
  _default_type_str = '_default_%s_type'
  _possible_types_str = '_possible_%s_types'
  _array_type_str = '_%s_array_types'
  _default_array_type_str = '_default_%s_array_type'
  _possible_array_types_str = '_possible_%s_array_types'
  _labels_str = '_%s_labels'
  _arrays_str = '_%s_arrays'
  _user_selected_labels_str = '_user_selected_%s_labels'
  _all_labels_str = '_all_%s_labels'

  # template error message
  _unrecognized_type_error_str = 'Unrecognized %s type, "%s," possible choices are %s.'

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

  def load_miller_array_phil_extract(self, phil_extract, process_files=True):
    '''
    Load custom PHIL extract
    '''
    self._load_miller_array_phil_extract(MillerArrayDataManager.datatype,
                                         phil_extract,
                                         process_files=process_files)

  def add_miller_array(self, filename, data):
    self._add(MillerArrayDataManager.datatype, filename, data)
    self._filter_miller_array_child_datatypes(filename)

  def set_default_miller_array_type(self, array_type=None):
    return self._set_default_miller_array_type(MillerArrayDataManager.datatype,
                                               array_type)

  def get_default_miller_array_array_type(self):
    return self._get_default_miller_array_array_type(MillerArrayDataManager.datatype)

  def set_default_miller_array_array_type(self, array_type=None):
    return self._set_default_miller_array_array_type(MillerArrayDataManager.datatype,
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

  def set_miller_array_array_type(self, filename=None, label=None, array_type=None):
    return self._set_miller_array_array_type(MillerArrayDataManager.datatype,
                                             filename, label, array_type)

  def get_miller_array_array_type(self, filename=None, label=None):
    return self._get_miller_array_array_type(MillerArrayDataManager.datatype,
                                             filename, label)

  def get_miller_array_all_labels(self, filename=None):
    return self._get_all_array_labels(MillerArrayDataManager.datatype, filename)

  def get_miller_array_labels(self, filename=None, return_all=True):
    '''
    Returns a list of array labels. If the list has no items, the
    return_all parameter controls whether or not to return all the labels

    Parameters
    ----------
    filename : str
      The filename for the labels, if None, the default filename is used
    return_all : bool, optional
      Keeps original behavior of returning all the labels

    Returns
    -------
    labels : list
      The list of labels
    '''
    labels = self._get_array_labels(MillerArrayDataManager.datatype, filename)
    if len(labels) == 0 and return_all:
      labels = self._get_all_array_labels(MillerArrayDataManager.datatype, filename)
    return labels

  def get_miller_array_user_selected_labels(self, filename=None):
    return self._get_user_selected_array_labels(MillerArrayDataManager.datatype, filename)

  def set_miller_array_user_selected_labels(self, filename=None, labels=None):
    '''
    Set the array of selected labels.

    Parameters
    ----------
    filename : str
      The filename for the labels, if None, the default filename is used
    labels : list
      The list of labels to set, if None, all labels in the file are selected
    '''
    if filename is None:
      filename = self.get_default_miller_array_name()
    all_labels = self.get_miller_array_all_labels(filename)
    if labels is None:
      labels = all_labels
    for label in labels:
      if label not in all_labels:
        raise Sorry('{} does not exist in {}.'.format(label, filename))
    getattr(self, self._user_selected_labels_str % MillerArrayDataManager.datatype)[filename] = labels

  def get_miller_array_types(self, filename=None):
    '''
    Returns a dict of array types, keyed by label
    '''
    return self._get_array_types(MillerArrayDataManager.datatype, filename)

  def get_miller_array_array_types(self, filename=None):
    '''
    Returns a dict of array array_types, keyed by label
    '''
    return self._get_array_array_types(MillerArrayDataManager.datatype, filename)

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
    if filename not in self.get_miller_array_names():
      self._process_file(MillerArrayDataManager.datatype, filename)
      self._filter_miller_array_child_datatypes(filename)
    return filename

  def _detect_miller_array_array_type(self, array):
    '''
    Convenience function for getting the array_type of a Miller array


    Parameters
    ----------
    array : cctbx.miller.array
      The Miller array object

    Returns
    -------
    array_type : str
      The array_type as determined by the is_*_array functions
    '''
    array_type = self._default_miller_array_array_type
    if array.is_xray_amplitude_array():
      array_type = 'amplitude'
    elif array.is_xray_intensity_array():
      array_type = 'intensity'
    elif array.is_complex_array():
      array_type = 'complex'
    elif array.is_hendrickson_lattman_array():
      array_type = 'hendrickson_lattman'
    elif array.is_integer_array():
      array_type = 'integer'
    elif array.is_bool_array():
      array_type = 'bool'
    elif array.is_nonsense():
      array_type = 'nonsense'
    return array_type

  def filter_miller_array_arrays(self, filename):
    '''
    Populate data structures with all arrays
    '''
    if filename not in self._miller_array_arrays.keys():
      self._miller_array_arrays[filename] = {}
    if filename not in self._miller_array_types.keys():
      self._miller_array_types[filename] = {}
    if filename not in self._miller_array_array_types.keys():
      self._miller_array_array_types[filename] = {}
    merge_equivalents = 'miller_array_skip_merge' not in self.custom_options
    miller_arrays = self.get_miller_array(filename).\
      as_miller_arrays(merge_equivalents=merge_equivalents)
    labels = []
    for array in miller_arrays:
      label = array.info().label_string()
      labels.append(label)
      self._miller_array_arrays[filename][label] = array
      self._miller_array_types[filename][label] = self._default_miller_array_type
      self._miller_array_array_types[filename][label] = self._detect_miller_array_array_type(array)
    self._all_miller_array_labels[filename] = labels
    self._miller_array_labels[filename] = []
    self._user_selected_miller_array_labels[filename] = []

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
                                 ignore_intensities_if_amplitudes_present=False,
                                 logger=None):
    '''
    Return the file server for a single miller_array file or mulitple files

    :param filenames:         list of filenames or None
    :param labels:            list of lists of labels or None
    :param array_type:        "x_ray", "neutron", "electron", or None
    :param crystal_symmetry:  cctbx.crystal.symmetry object or None
    :param force_symmetry:    bool or None
    :param ignore_intensities_if_amplitudes_present: bool
    :param logger:            defaults to self.logger (multi_out)

    The order in filenames and labels should match, e.g. labels[0] should be the
    labels for filenames[0]. The lengths of filenames and labels should be equal
    as well. If all the labels in a file are to be added, set the labels entry
    to None, e.g. labels[0] = None.

    If array_type is None, files of any type are allowed.

    If ignore_intensities_if_amplitudes_present is set to True, intensity
    arrays are not automatically added to the reflection_file_server if
    amplitude arrays have been selected (as user_selected_arrays). This
    is to avoid fmodel from prioritizing intensity arrays if both are
    present. If both intensity and amplitude arrays are selected,
    setting this parameter to True will not ignore the selected intensity
    arrays.

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

    # determine types of selected arrays across all files and decide
    # whether to ignore intensity arrays that are not explicitly selected
    selected_types = set()
    for i, filename in enumerate(filenames):
      current_selected_labels = self.get_miller_array_user_selected_labels(filename)
      if len(current_selected_labels) > 0:
        for j, label in enumerate(current_selected_labels):
          label = self._match_label(label, self.get_miller_arrays(filename=filename))
          selected_types.add(self.get_miller_array_array_types(filename)[label])
    if 'amplitude' in selected_types and ignore_intensities_if_amplitudes_present:
      selected_types.add('intensity')

    # check for user selected labels
    selected_labels = deepcopy(labels)
    for i, filename in enumerate(filenames):
      current_selected_labels = deepcopy(self.get_miller_array_user_selected_labels(filename))
      current_all_labels = labels[i]
      if labels[i] is None:
        current_all_labels = self.get_miller_array_all_labels(filename)
      if len(current_selected_labels) > 0:
        # add selected labels
        for j, label in enumerate(current_selected_labels):
          label = self._match_label(label, self.get_miller_arrays(filename=filename))
          current_selected_labels[j] = label
        # add remaining labels that are a different type
        for label in current_all_labels:
          if self.get_miller_array_array_types(filename)[label] not in selected_types \
            or self.get_miller_array_array_types(filename)[label] == 'unknown':
            current_selected_labels.append(label)
        selected_labels[i] = current_selected_labels
    labels = selected_labels

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
    miller_array_labels = []  # filename:array_label
    def get_array_label(filename, array):
      '''
      Convenience function for creating filename:array_label label
      '''
      return filename + ':' + array.info().label_string()
    merge_equivalents = 'miller_array_skip_merge' not in self.custom_options
    for filename, file_labels in zip(filenames, labels):
      file_arrays = self.get_miller_array(filename).\
        as_miller_arrays(crystal_symmetry=crystal_symmetry,
                         force_symmetry=force_symmetry,
                         merge_equivalents=merge_equivalents)
      if file_labels is None:
        file_labels = self.get_miller_array_all_labels(filename)
      for miller_array in file_arrays:
        label_name = miller_array.info().label_string()
        # check array label
        if label_name in file_labels:
          # check array type
          if (array_type is None
              or array_type in self.get_miller_array_type(filename, label_name)):
            miller_arrays.append(miller_array)
            miller_array_labels.append(get_array_label(filename, miller_array))

    # final check of selected labels with scores
    problem_arrays = set()  # arrays with the same score as selected arrays
    # -------------------------------------------------------------------------
    def filter_arrays_by_scores(filename, arrays, test_scores):
      '''
      Convenience function for using scores to select arrays for removal
      '''
      remove_set = set()
      selected_array_labels = [get_array_label(filename, selected_array) for selected_array in arrays]
      for selected_array_label in selected_array_labels:
        if selected_array_label in score_dict:
          score = score_dict[selected_array_label]
          if score > 0 and test_scores.count(score) > 1:
            for array, array_label in zip(miller_arrays, miller_array_labels):
              if score_dict[array_label] == score and array_label not in selected_array_labels:
                remove_set.add(array)
      return remove_set
    # -------------------------------------------------------------------------
    for scores in [
      get_experimental_phases_scores(miller_arrays, False),
      get_phase_scores(miller_arrays),
      get_r_free_flags_scores(miller_arrays, None).scores,
      ]:
      score_dict = {miller_array_labels[i]:scores[i] for i in range(len(miller_arrays))}
      for filename in filenames:
        selected_arrays = self.get_miller_arrays(
          labels=self.get_miller_array_user_selected_labels(filename=filename),
          filename=filename)
        problem_arrays.update(filter_arrays_by_scores(filename, selected_arrays, scores))
    new_arrays = [ma for ma in miller_arrays if ma not in problem_arrays]

    miller_arrays = new_arrays

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
    # self._miller_array_types = {}                 # [filename] = type dict
    # self._miller_array_array_types = {}           # [filename] = type dict
    # self._miller_array_all_labels = {}            # [filename] = label list
    # self._miller_array_labels = {}                # [filename] = label list
    # self._miller_array_arrays = {}                # [filename] = array dict
    # self._user_selected_miller_array_labels = {}  # [filename] = array dict
    setattr(self, self._type_str % datatype, {})
    setattr(self, self._default_type_str % datatype, ['x_ray'])
    setattr(self, self._possible_types_str % datatype,
            ['x_ray', 'neutron', 'electron'])
    setattr(self, self._array_type_str % datatype, {})
    setattr(self, self._default_array_type_str % datatype, 'unknown')
    setattr(self, self._possible_array_types_str % datatype,
            ['unknown', 'amplitude', 'bool', 'complex', 'hendrickson_lattman', 'integer',
             'intensity', 'nonsense'])
    setattr(self, self._all_labels_str % datatype, {})
    setattr(self, self._labels_str % datatype, {})
    setattr(self, self._arrays_str % datatype, {})
    setattr(self, self._user_selected_labels_str % datatype, {})

    # custom PHIL section
    # all the array labels in a file is not stored
    # user_selected_labels stores the actual use selection
    # labels stores the full selected label and metadata
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
    .help = Storage for all labels in the file
  {
    name = None
      .type = str
    type = *%s
      .type = choice(multi=True)
    array_type = *%s
      .type = choice(multi=False)
  }
  user_selected_labels = None
    .type = str
    .multiple = True
    .help = Storage for user selected labels (can be shorter)
    .style = hidden
}
''' % (datatype,
       ' '.join(getattr(self, self._possible_types_str % datatype)),
       ' '.join(getattr(self, self._possible_array_types_str % datatype)))
    # add fmodel PHIL
    if self.supports('model'):
      custom_phil_str += '''
fmodel {
  include scope iotbx.extract_xtal_data.xray_data_str_no_filenames
  include scope iotbx.extract_xtal_data.neutron_data_str_no_filenames
}
'''

    # custom PHIL scope
    setattr(self, '_custom_%s_phil' % datatype,
            iotbx.phil.parse(custom_phil_str, process_includes=True))

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
      types = self._get_array_types(datatype, filename=filename)
      array_types = self._get_array_array_types(datatype, filename=filename)
      user_selected_labels = self._get_user_selected_array_labels(datatype, filename=filename)
      arrays = self.get_miller_arrays(filename=filename)
      labels = [self._match_label(label, arrays) for label in user_selected_labels]
      if len(labels) > len(types.keys()):
        raise Sorry('Some labels do not have types.\n{}\n{}'.format(labels, list(types.keys())))
      labels_extract = []
      for label in labels:
        label_extract = deepcopy(item_extract.labels[0])
        label_extract.name = label
        if label not in types.keys():
          raise Sorry('The label, %s, is not recognized to be %s data.' % (label, datatype))
        label_extract.type = types[label]
        label_extract.array_type = array_types[label]
        labels_extract.append(label_extract)
      item_extract.labels = labels_extract
      item_extract.user_selected_labels = user_selected_labels
      extract.append(item_extract)
    return extract

  def _match_label(self, label, miller_arrays):
    '''
    Convenience function for matching partially specified labels

    Parameters
    ----------
    label : str
      label to check
    miller_arrays : list (use self.get_miller_arrays)
      arrays to check against

    Return
    ------
    matching_label: str or None
      The most closely matching label or None if no match is found
    '''
    label_match = None
    lbl_tab = label_table(miller_arrays=miller_arrays, err=os.devnull)
    match = lbl_tab.match_data_label(
      label=label, command_line_switch='miller_array.label', f=self.logger)
    if match is not None:
      label_match = match.info().label_string()
    return label_match

  def _reconcile_user_labels(self, filename):
    '''
    Ensures that user_selected_labels are consistent between miller_array
    and map_coefficient types. Also ensures that user_selected_labels
    are uniquely selecting arrays.
    '''
    miller_arrays = self.get_miller_arrays(filename=filename)
    ma_user_selected_labels = self.get_miller_array_user_selected_labels(filename)
    matched_ma_user_selected_labels = [
      self._match_label(label, miller_arrays) for label in ma_user_selected_labels]
    datatype_user_selected_labels = []
    matched_datatype_user_selected_labels = []
    if self.supports('map_coefficients') and self.has_map_coefficients():
      for datatype in self.miller_array_child_datatypes:
        datatype_user_selected_labels = self._get_user_selected_array_labels(datatype, filename)
        matched_datatype_user_selected_labels = [
          self._match_label(label, miller_arrays) for label in datatype_user_selected_labels]
        datatype_labels = self._get_array_labels(datatype, filename)
        # append to parent labels if missing
        for label in datatype_user_selected_labels:
          matched_label = self._match_label(label, miller_arrays)
          if label not in ma_user_selected_labels \
            and matched_label not in matched_ma_user_selected_labels:
            ma_user_selected_labels.append(label)
        # append to datatype labels if missing
        for label in ma_user_selected_labels:
          matched_label = self._match_label(label, miller_arrays)
          if label not in datatype_user_selected_labels \
            and matched_label not in matched_datatype_user_selected_labels \
            and matched_label in datatype_labels:
            datatype_user_selected_labels.append(label)

    # check for duplicates
    for label_pair in [(matched_ma_user_selected_labels, ma_user_selected_labels),
                       (matched_datatype_user_selected_labels, datatype_user_selected_labels)]:
      if len(set(label_pair[0])) < len(label_pair[0]):
        raise Sorry('''
There are duplicate user_selected_labels. Please only specify unique labels.
user_selected_labels: %s
      matched labels: %s
''' % (label_pair[1], label_pair[0]))

  def _load_miller_array_phil_extract(self, datatype, phil_extract, process_files=True):
    if self.supports('model'):
      self._fmodel_phil_scope = self.master_phil.format(python_object=phil_extract)
    extract = phil_extract.data_manager
    extract = getattr(extract, '%s' % datatype)
    for item_extract in extract:
      if not hasattr(item_extract, 'file'):
        raise Sorry('This PHIL is not properly defined for the %s datatype.\n There should be a parameter for the filename ("file").\n')

      # process file, or last loaded file
      if item_extract.file is None:
        filenames = self.get_miller_array_names()
        if len(filenames) == 0:
          raise Sorry('No reflection files are available to continue processing PHIL.')
        item_extract.file = filenames[-1]
      if process_files:
        getattr(self, 'process_%s_file' % datatype)(item_extract.file)

        # check labels (if available)
        if len(item_extract.labels) > 0 or len(item_extract.user_selected_labels) > 0:
          # all labels in file
          file_labels = getattr(self, self._all_labels_str % 'miller_array')[item_extract.file]

          # labels from PHIL
          phil_labels = []
          phil_types = {}
          phil_array_types = {}
          phil_user_selected_labels = []
          specified_labels = item_extract.labels + item_extract.user_selected_labels
          for label in specified_labels:
            label_name = label
            if hasattr(label_name, 'name'):
              label_name = label_name.name
            phil_user_selected_labels.append(label_name)
            if label_name not in file_labels:

              # try matching
              label_match = self._match_label(label_name, self.get_miller_arrays(filename=item_extract.file))
              if label_match is None:
                raise Sorry('The label, %s, could not be found in %s.' %
                  (label_name, item_extract.file))
              else:
                label_name = label_match
            phil_labels.append(label_name)
            if hasattr(label, 'type'):
              if not self._is_valid_miller_array_type(datatype, label.type):
                raise Sorry(self._unrecognized_type_error_str %
                            (datatype,
                            ', '.join(label.type),
                            ', '.join(getattr(self, self._possible_types_str % datatype))))
              phil_types[label_name] = label.type
            if hasattr(label, 'array_type'):
              if label.array_type not in getattr(self, self._possible_array_types_str % datatype):
                raise Sorry(self._unrecognized_type_error_str %
                            (datatype, label.array_type, ', '.join(
                              getattr(self, self._possible_array_types_str % datatype))))
              # if a non-default type is set, assume it is manually set
              # otherwise, automatically detect array_type
              array_type = label.array_type
              if array_type == self.get_default_miller_array_array_type():
                array = self.get_miller_arrays(labels=[label_name], filename=item_extract.file)[0]
                array_type = self._detect_miller_array_array_type(array)
              phil_array_types[label_name] = array_type

          # update storage
          labels_storage = getattr(self, self._labels_str % datatype)
          if item_extract.file not in labels_storage.keys():
            raise Sorry('The file, %s, does not seem to have %s data.' % (item_extract.file, datatype))
          labels_storage = labels_storage[item_extract.file]
          types_storage = getattr(self, self._type_str % datatype)[item_extract.file]
          array_types_storage = getattr(self, self._array_type_str % datatype)[item_extract.file]
          for label in phil_labels:
            if label not in labels_storage:
              labels_storage.append(label)
            if label in phil_types:
              types_storage[label] = phil_types[label]
            if label in phil_array_types:
              array_types_storage[label] = phil_array_types[label]
          user_selected_labels_storage = getattr(self, self._user_selected_labels_str % datatype)[item_extract.file]
          for label in phil_user_selected_labels:
            if label not in user_selected_labels_storage:
              user_selected_labels_storage.append(label)

        # ensure user_selected_labels are consistent
        self._reconcile_user_labels(filename=item_extract.file)

  def _is_valid_miller_array_type(self, datatype, array_type):
    """
    Convenience function for checking if the array type is valid
    This will also check that model_type is a list to conform with the
    PHIL parameter

    Parameters
    ----------
    datatype: str
      The datatype (i.e. miller_array)
    array_type: list
      The model_type(s) to check.

    Returns
    -------
    bool:
    """
    if not isinstance(array_type, list):
      raise Sorry('The type argument must be a list.')
    if len(array_type) == 0:
      return False
    valid = True
    for at in array_type:
      valid = valid and (at in getattr(self, self._possible_types_str % datatype))
    return valid

  def _set_default_miller_array_type(self, datatype, array_type):
    if not self._is_valid_miller_array_type(datatype, array_type):
      raise Sorry(self._unrecognized_type_error_str %
                  (datatype,
                   ', '.join(array_type),
                   ', '.join(getattr(self, self._possible_types_str % datatype))))
    setattr(self, self._default_type_str % datatype, array_type)

  def _get_default_miller_array_type(self, datatype):
    return getattr(self, self._default_type_str % datatype)

  def _set_default_miller_array_array_type(self, datatype, array_type):
    if array_type not in getattr(self, self._possible_array_types_str % datatype):
      raise Sorry(self._unrecognized_type_error_str %
                  (datatype, array_type,
                   ', '.join(getattr(self, self._possible_array_types_str % datatype))))
    setattr(self, self._default_array_type_str % datatype, array_type)

  def _get_default_miller_array_array_type(self, datatype):
    return getattr(self, self._default_array_type_str % datatype)

  def _set_miller_array_type(self, datatype, filename=None, label=None,
                             array_type=None):
    if filename is None:
      filename = self._get_default_name(datatype)
    if label is None:
      label = self._get_all_array_labels(datatype, filename)[0]
    if array_type is None:
      array_type = [getattr(self, self._default_type_str % datatype)]
    if not self._is_valid_miller_array_type(datatype, array_type):
      raise Sorry(self._unrecognized_type_error_str %
                  (datatype,
                  ', '.join(array_type),
                  ', '.join(getattr(self, self._possible_types_str % datatype))))
    getattr(self, self._type_str % datatype)[filename][label] = array_type

  def _get_miller_array_type(self, datatype, filename=None, label=None):
    if filename is None:
      filename = self._get_default_name(datatype)
    if label is None:
      label = self._get_all_array_labels(datatype, filename)[0]
    types = self._get_array_types(datatype, filename)
    return types.get(label, getattr(self, self._default_type_str % datatype))

  def _set_miller_array_array_type(self, datatype, filename=None, label=None,
                                   array_type=None):
    if filename is None:
      filename = self._get_default_name(datatype)
    if label is None:
      label = self._get_all_array_labels(datatype, filename)[0]
    if array_type is None:
      array_type = getattr(self, self._default_array_type_str % datatype)
    elif array_type not in getattr(self, self._possible_array_types_str % datatype):
      raise Sorry(self._unrecognized_type_error_str %
                  (datatype,
                   array_type,
                   ', '.join(getattr(self, self._possible_array_types_str % datatype))))
    getattr(self, self._array_type_str % datatype)[filename][label] = array_type

  def _get_miller_array_array_type(self, datatype, filename=None, label=None):
    if filename is None:
      filename = self._get_default_name(datatype)
    if label is None:
      label = self._get_all_array_labels(datatype, filename)[0]
    types = self._get_array_array_types(datatype, filename)
    return types.get(label, getattr(self, self._default_array_type_str % datatype))

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

  def _get_all_array_labels(self, datatype, filename=None):
    filename = self._check_miller_array_default_filename(datatype, filename)
    storage_dict = getattr(self, self._all_labels_str % datatype)
    self._check_miller_array_storage_dict(datatype, storage_dict, filename)
    labels = storage_dict[filename]
    return labels

  def _get_array_labels(self, datatype, filename=None):
    filename = self._check_miller_array_default_filename(datatype, filename)
    storage_dict = getattr(self, self._labels_str % datatype)
    self._check_miller_array_storage_dict(datatype, storage_dict, filename)
    labels = storage_dict[filename]
    return labels

  def _get_user_selected_array_labels(self, datatype, filename=None):
    filename = self._check_miller_array_default_filename(datatype, filename)
    storage_dict = getattr(self, self._user_selected_labels_str % datatype)
    self._check_miller_array_storage_dict(datatype, storage_dict, filename)
    labels = storage_dict[filename]
    return labels

  def _get_array_types(self, datatype, filename=None):
    filename = self._check_miller_array_default_filename(datatype, filename)
    storage_dict = getattr(self, self._type_str % datatype)
    self._check_miller_array_storage_dict(datatype, storage_dict, filename)
    types = storage_dict[filename]
    return types

  def _get_array_array_types(self, datatype, filename=None):
    filename = self._check_miller_array_default_filename(datatype, filename)
    storage_dict = getattr(self, self._array_type_str % datatype)
    self._check_miller_array_storage_dict(datatype, storage_dict, filename)
    types = storage_dict[filename]
    return types

  def _get_arrays(self, datatype, filename=None, labels=None):
    filename = self._check_miller_array_default_filename(datatype, filename)
    if filename not in getattr(self, 'get_%s_names' % datatype)():
      self.process_miller_array_file(filename)
    if labels is None:
      labels = self._get_all_array_labels(datatype, filename=filename)
    else:
      if not isinstance(labels, list):
        raise Sorry('The labels argument should be a list of labels')
    storage_dict = getattr(self, self._arrays_str % datatype)
    self._check_miller_array_storage_dict(datatype, storage_dict, filename)
    arrays = []
    for label in labels:
      arrays.append(self._get_array_by_label(datatype, filename, label))
    return arrays

  def _get_array_by_label(self, datatype, filename, label, try_matching=True):
    storage_dict = getattr(self, self._arrays_str % datatype)
    if label in storage_dict[filename].keys():
      return storage_dict[filename][label]
    elif label in self._get_user_selected_array_labels(datatype, filename) and try_matching:
      matched_label = self._match_label(label, self._get_arrays(datatype, filename, labels=None))
      return self._get_array_by_label(datatype, filename, matched_label, try_matching=False)
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

    user_labels = getattr(self, self._user_selected_labels_str % datatype)
    user_labels[filename] = []

    data = self.get_miller_array(filename)
    merge_equivalents = 'miller_array_skip_merge' not in self.custom_options
    miller_arrays = data.as_miller_arrays(merge_equivalents=merge_equivalents)
    labels = []
    types = {}
    # array_types = {}
    datatype_dict = getattr(self, self._arrays_str % datatype)
    for array in miller_arrays:
      label = set(array.info().labels)
      common_labels = known_labels.intersection(label)
      if len(common_labels) > 0:
        label = array.info().label_string()
        labels.append(label)
        if filename not in datatype_dict.keys():
          datatype_dict[filename] = {}
        datatype_dict[filename][label] = array
        types[label] = getattr(self, self._default_type_str % datatype)
        # array_types[label] = getattr(self, self._default_array_type_str % datatype)

    # if arrays exist, start tracking
    if len(labels) > 0:
      getattr(self, self._labels_str % datatype)[filename] = labels
      getattr(self, self._type_str % datatype)[filename] = types
      # getattr(self, self._array_type_str % datatype)[filename] = array_types
      self._add(datatype, filename, data)

# =============================================================================
# end
