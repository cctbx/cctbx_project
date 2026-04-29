'''
Child class of MillerArrayDataManager for handling map coefficients

This will eventually be deprecated.
The array_type PHIL parameter should be used instead
'''
from __future__ import absolute_import, division, print_function

from iotbx.data_manager.miller_array import MillerArrayDataManager
from iotbx.cif_mtz_data_labels import mtz_map_coefficient_labels, \
  cif_map_coefficient_labels
from libtbx import Auto

# =============================================================================
class MapCoefficientsDataManager(MillerArrayDataManager):

  datatype = 'map_coefficients'

  # ---------------------------------------------------------------------------
  # Map coefficients

  def add_map_coefficients_phil_str(self):
    '''
    Add custom PHIL and storage for labels
    '''
    return self._add_miller_array_phil_str(MapCoefficientsDataManager.datatype)

  def export_map_coefficients_phil_extract(self):
    '''
    Export custom PHIL extract
    '''
    return self._export_miller_array_phil_extract(
      MapCoefficientsDataManager.datatype)

  def load_map_coefficients_phil_extract(self, phil_extract, process_files=True):
    '''
    Load custom PHIL extract
    '''
    self._load_miller_array_phil_extract(MapCoefficientsDataManager.datatype,
                                         phil_extract,
                                         process_files=process_files)

  def add_map_coefficients(self, filename, data):
    self.add_miller_array(filename, data)

  def set_default_map_coefficients_type(self, array_type=None):
    return self._set_default_miller_array_type(
      MapCoefficientsDataManager.datatype, array_type)

  def get_default_map_coefficients_type(self):
    return self._get_default_miller_array_type(
      MapCoefficientsDataManager.datatype)

  def set_default_map_coefficients(self, filename):
    return self._set_default(MapCoefficientsDataManager.datatype, filename)

  def get_map_coefficients(self, filename=None):
    '''
    Returns the main file object
    '''
    return self._get(MapCoefficientsDataManager.datatype, filename)

  def set_map_coefficients_type(self, filename=None, label=None, array_type=None):
    return self._set_miller_array_type(MapCoefficientsDataManager.datatype,
                                       filename, label, array_type)

  def get_map_coefficients_type(self, filename=None, label=None):
    return self._get_miller_array_type(MapCoefficientsDataManager.datatype,
                                       filename, label)

  def get_map_coefficients_all_labels(self, filename=None):
    return self.get_miller_array_all_labels(filename)

  def get_map_coefficients_labels(self, filename=None):
    '''
    Returns a list of array labels
    '''
    return self._get_array_labels(MapCoefficientsDataManager.datatype, filename)

  def get_map_coefficients_user_selected_labels(self, filename=None):
    '''
    Returns a list of user selected array labels
    '''
    return self._get_user_selected_array_labels(MapCoefficientsDataManager.datatype, filename)

  def get_map_coefficients_types(self, filename=None):
    '''
    Returns a dict of array types, keyed by label
    '''
    return self._get_array_types(MapCoefficientsDataManager.datatype, filename)

  def get_map_coefficients_arrays(self, labels=None, filename=None):
    '''
    Returns a list of map coefficients from the file
    '''
    return self._get_arrays(MapCoefficientsDataManager.datatype, labels=labels,
                            filename=filename)

  def get_map_coefficients_names(self):
    return self._get_names(MapCoefficientsDataManager.datatype)

  def get_default_map_coefficients_name(self):
    return self._get_default_name(MapCoefficientsDataManager.datatype)

  def remove_map_coefficients(self, filename):
    return self._remove(MapCoefficientsDataManager.datatype, filename)

  def has_map_coefficients(
      self, expected_n=1, exact_count=False, raise_sorry=False):
    return self._has_data(
      MapCoefficientsDataManager.datatype, expected_n=expected_n,
      exact_count=exact_count, raise_sorry=raise_sorry)

  def process_map_coefficients_file(self, filename):
    self.process_miller_array_file(filename)

  def filter_map_coefficients_arrays(self, filename):
    '''
    Populate data structures by checking labels in miller_arrays to determine
    type and by setting all complex miller arrays as map coefficients
    '''
    # check for labels
    known_labels = mtz_map_coefficient_labels.union(cif_map_coefficient_labels)
    datatype = MapCoefficientsDataManager.datatype
    self._child_filter_arrays(datatype, filename, known_labels)

    # check for complex arrays
    data = self.get_miller_array(filename)
    miller_arrays = data.as_miller_arrays()
    current_labels = []
    labels = []
    types = {}
    array_types = {}
    datatype_dict = getattr(self, '_%s_arrays' % datatype)
    for array in miller_arrays:
      label = array.info().label_string()
      if array.is_complex_array() and label not in current_labels:
        labels.append(label)
        if filename not in datatype_dict.keys():
          datatype_dict[filename] = dict()
        datatype_dict[filename][label] = array
        types[label] = getattr(self, self._default_type_str % datatype)
        array_types[label] = 'complex'

    # update data structures
    if len(labels) > 0:
      current_labels = getattr(self, self._labels_str % datatype)
      if filename not in current_labels:
        current_labels[filename] = labels
      else:
        for label in labels:
          if label not in current_labels[filename]:
            current_labels[filename].append(label)
      current_types = getattr(self, self._type_str % datatype)
      if filename not in current_types:
        current_types[filename] = types
      else:
        current_types[filename].update(types)
      current_array_types = getattr(self, self._array_type_str % datatype)
      if filename not in current_array_types:
        current_array_types[filename] = array_types
      else:
        current_array_types[filename].update(array_types)
      current_user_selected_labels = getattr(self, self._user_selected_labels_str % datatype)
      current_user_selected_labels[filename] = []
      self._add(datatype, filename, data)

  def write_map_coefficients_file(
      self, mtz_object, filename=Auto, overwrite=Auto):
    return self.write_miller_array_file(
      mtz_object, filename=filename, overwrite=overwrite)

# =============================================================================
# end
