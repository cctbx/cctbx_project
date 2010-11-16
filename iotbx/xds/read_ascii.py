from cctbx import miller
from cctbx import crystal
from cctbx import uctbx
from cctbx.array_family import flex
import sys

def get_rhs(line):
  return line.split("=", 1)[1]

class reader(object):

  def __init__(self, file_handle, header_only=False):
    "http://www.mpimf-heidelberg.mpg.de/~kabsch/xds/"
    f = iter(file_handle)
    flds = f.next().split()
    assert flds[0] == "!FORMAT=XDS_ASCII"
    assert flds[1] == "MERGE=TRUE"
    if   (flds[2] == "FRIEDEL'S_LAW=FALSE"):
      self.anomalous_flag = True
    elif (flds[2] == "FRIEDEL'S_LAW=TRUE"):
      self.anomalous_flag = False
    else:
      raise RuntimeError, "Expected FRIEDEL'S_LAW=FALSE|TRUE"
    self.unit_cell = None
    self.space_group_number = None
    self.number_of_items_in_each_data_record = None
    self.miller_index_columns = [None, None, None]
    self.iobs_column = None
    self.sigma_iobs_column = None
    for line in f:
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
    else:
      self.miller_indices = flex.miller_index()
      self.iobs = flex.double()
      self.sigma_iobs = flex.double()
      for line in f:
        if (line.startswith("!END_OF_DATA")):
          break
        data = line.split()
        assert len(data) == self.number_of_items_in_each_data_record
        h = [int(data[i]) for i in self.miller_index_columns]
        self.miller_indices.append(h)
        self.iobs.append(float(data[self.iobs_column]))
        self.sigma_iobs.append(float(data[self.sigma_iobs_column]))

  def column_index(self, line):
    i_column = int(get_rhs(line))-1
    assert 0 <= i_column < self.number_of_items_in_each_data_record
    return i_column

  def crystal_symmetry(self):
    return crystal.symmetry(
      unit_cell=self.unit_cell,
      space_group_symbol=self.space_group_number)

  def as_miller_array(self,
        crystal_symmetry=None,
        force_symmetry=False,
        merge_equivalents=True,
        base_array_info=None):
    if (base_array_info is None):
      base_array_info = miller.array_info(source_type="xds_ascii")
    crystal_symmetry_from_file = self.crystal_symmetry()
    return (miller.array(
      miller_set=miller.set(
        crystal_symmetry=crystal_symmetry_from_file.join_symmetry(
          other_symmetry=crystal_symmetry,
          force=force_symmetry),
        indices=self.miller_indices,
        anomalous_flag=self.anomalous_flag),
      data=self.iobs,
      sigmas=self.sigma_iobs)
      .set_info(base_array_info.customized_copy(
        labels=["iobs", "sigma_iobs"],
        crystal_symmetry_from_file=crystal_symmetry_from_file))
      .set_observation_type_xray_intensity())

  def as_miller_arrays(self,
        crystal_symmetry=None,
        force_symmetry=False,
        merge_equivalents=True,
        base_array_info=None):
    return [self.as_miller_array(
      crystal_symmetry=crystal_symmetry,
      force_symmetry=force_symmetry,
      merge_equivalents=merge_equivalents,
      base_array_info=base_array_info)]

if (__name__ == "__main__"):
  reader(open(sys.argv[1])).as_miller_array().show_comprehensive_summary()
