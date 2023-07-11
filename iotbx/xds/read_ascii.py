from __future__ import absolute_import, division, print_function
from cctbx import miller
from cctbx import crystal
from cctbx import uctbx
from cctbx.array_family import flex
import sys

def get_rhs(line):
  return line.split("=", 1)[1]

class reader(object):

  def __init__(self, file_handle, header_only=False, allow_unmerged=True):
    "http://www.mpimf-heidelberg.mpg.de/~kabsch/xds/"
    flds = file_handle.readline().split()
    assert flds[0] == "!FORMAT=XDS_ASCII"
    assert (allow_unmerged) or (flds[1] == "MERGE=TRUE")
    self.unmerged_data = (flds[1] == "MERGE=FALSE")
    if   (flds[2] == "FRIEDEL'S_LAW=FALSE"):
      self.anomalous_flag = True
    elif (flds[2] == "FRIEDEL'S_LAW=TRUE"):
      self.anomalous_flag = self.unmerged_data
    else:
      raise RuntimeError("Expected FRIEDEL'S_LAW=FALSE|TRUE")
    self.unit_cell = None
    self.space_group_number = None
    self.number_of_items_in_each_data_record = None
    self.miller_index_columns = [None, None, None]
    self.iobs_column = None
    self.sigma_iobs_column = None
    self.zd_column = None
    self.wavelength = None
    for line in file_handle:
      if (line.startswith("!SPACE_GROUP_NUMBER=")):
        self.space_group_number = int(get_rhs(line))
        assert 1 <= self.space_group_number <= 230
      elif (line.startswith("!UNIT_CELL_CONSTANTS=")):
        self.unit_cell = uctbx.unit_cell(get_rhs(line))
      elif (line.startswith("!NUMBER_OF_ITEMS_IN_EACH_DATA_RECORD=")):
        self.number_of_items_in_each_data_record = int(get_rhs(line))
      elif (line.startswith("!ITEM_H=")):
        self.miller_index_columns[0] = self.column_index(line)
      elif (line.startswith("!ITEM_K=")):
        self.miller_index_columns[1] = self.column_index(line)
      elif (line.startswith("!ITEM_L=")):
        self.miller_index_columns[2] = self.column_index(line)
      elif (line.startswith("!ITEM_IOBS=")):
        self.iobs_column = self.column_index(line)
      elif (line.startswith("!ITEM_SIGMA(IOBS)=")):
        self.sigma_iobs_column = self.column_index(line)
      elif (line.startswith("!ITEM_ZD=")):
        self.zd_column = self.column_index(line)
      elif (line.startswith("!X-RAY_WAVELENGTH=")):
        self.wavelength = float(get_rhs(line))
      elif (line.startswith("!END_OF_HEADER")):
        break
    assert self.unit_cell is not None
    assert self.space_group_number is not None
    assert self.column_index is not None
    assert None not in self.miller_index_columns
    assert self.iobs_column is not None
    assert self.sigma_iobs_column is not None
    if (header_only):
      self.miller_indices = None
      self.iobs = None
      self.sigma_iobs = None
      self.zd = None
    else:
      self.miller_indices = flex.miller_index()
      self.iobs = flex.double()
      self.sigma_iobs = flex.double()
      self.zd = None
      if (self.zd_column is not None):
        self.zd = flex.double()
      for line in file_handle:
        if (line.startswith("!END_OF_DATA")):
          break
        data = line.split()
        assert len(data) == self.number_of_items_in_each_data_record
        sigma_value = float(data[self.sigma_iobs_column])
        if sigma_value < 0:
          continue # negative sigma values should be ignored (outliers)
        h = [int(data[i]) for i in self.miller_index_columns]
        self.miller_indices.append(h)
        self.iobs.append(float(data[self.iobs_column]))
        self.sigma_iobs.append(sigma_value)
        if (self.zd_column is not None):
          self.zd.append(float(data[self.zd_column]))

  def column_index(self, line):
    i_column = int(get_rhs(line))-1
    assert 0 <= i_column < self.number_of_items_in_each_data_record
    return i_column

  def crystal_symmetry(self):
    return crystal.symmetry(
      unit_cell=self.unit_cell,
      space_group_symbol=self.space_group_number)

  def miller_set(self,
        crystal_symmetry=None,
        force_symmetry=False,
        anomalous=None):
    crystal_symmetry_from_file = self.crystal_symmetry()
    if anomalous is None:
      anomalous = self.anomalous_flag
    return miller.set(
        crystal_symmetry=crystal_symmetry_from_file.join_symmetry(
          other_symmetry=crystal_symmetry,
          force=force_symmetry),
        indices=self.miller_indices,
        anomalous_flag=self.anomalous_flag)

  def as_miller_array(self,
        crystal_symmetry=None,
        force_symmetry=False,
        merge_equivalents=True,
        base_array_info=None,
        anomalous=None):
    if (base_array_info is None):
      base_array_info = miller.array_info(source_type="xds_ascii")
    crystal_symmetry_from_file = self.crystal_symmetry()
    array = (miller.array(
      miller_set=self.miller_set(
        crystal_symmetry=crystal_symmetry,
        force_symmetry=force_symmetry,
        anomalous=anomalous,
      ),
      data=self.iobs,
      sigmas=self.sigma_iobs)
      .set_info(base_array_info.customized_copy(
        labels=["iobs", "sigma_iobs"],
        crystal_symmetry_from_file=crystal_symmetry_from_file,
        wavelength=self.wavelength))
      .set_observation_type_xray_intensity())
    if (merge_equivalents):
      info = array.info()
      info.merged = True
      array = array.merge_equivalents().array().set_info(info)
    return array

  def as_miller_arrays(self,
        crystal_symmetry=None,
        force_symmetry=False,
        merge_equivalents=True,
        base_array_info=None,
        anomalous=None):
    return [self.as_miller_array(
      crystal_symmetry=crystal_symmetry,
      force_symmetry=force_symmetry,
      merge_equivalents=merge_equivalents,
      base_array_info=base_array_info,
      anomalous=anomalous,
    )]

  def batch_as_miller_array(self,
        crystal_symmetry=None,
        force_symmetry=False,
        base_array_info=None):
    if (base_array_info is None):
      base_array_info = miller.array_info(source_type="xds_ascii")
    crystal_symmetry_from_file = self.crystal_symmetry()
    return miller.array(
      miller_set=self.miller_set(
          crystal_symmetry=crystal_symmetry,
          force_symmetry=force_symmetry),
      data=self.zd).set_info(
        base_array_info.customized_copy(
          labels=["ZD"],
          crystal_symmetry_from_file=crystal_symmetry_from_file,
          wavelength=self.wavelength))

if (__name__ == "__main__"):
  reader(open(sys.argv[1])).as_miller_array().show_comprehensive_summary()
