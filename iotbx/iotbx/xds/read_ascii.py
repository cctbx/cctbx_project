from cctbx import miller
from cctbx import crystal
from cctbx import uctbx
from cctbx import sgtbx
from cctbx.array_family import flex
import sys

def get_rhs(line):
  return line.split("=", 1)[1]

class column_index_retriever:

  def __init__(self, line):
    self.number_of_items_in_each_data_record = int(get_rhs(line))

  def get(self, line):
    i_column = int(get_rhs(line))-1
    assert 0 <= i_column < self.number_of_items_in_each_data_record
    return i_column

def read(file_name):
  f = iter(open(file_name))
  flds = f.next().split()
  assert flds[0] == "!FORMAT=XDS_ASCII"
  assert flds[1] == "MERGE=TRUE"
  if   (flds[2] == "FRIEDEL'S_LAW=FALSE"):
    anomalous_flag = 0001
  elif (flds[2] == "FRIEDEL'S_LAW=TRUE"):
    anomalous_flag = 00000
  else:
    raise RuntimeError, "Expected FRIEDEL'S_LAW=FALSE|TRUE"
  unit_cell = None
  space_group_number = None
  column_index = None
  miller_index_columns = [None, None, None]
  iobs_column = None
  sigma_iobs_column = None
  for line in f:
    if (line.startswith("!SPACE_GROUP_NUMBER=")):
      space_group_number = int(get_rhs(line))
    elif (line.startswith("!UNIT_CELL_CONSTANTS=")):
      unit_cell = uctbx.unit_cell(get_rhs(line))
    elif (line.startswith("!NUMBER_OF_ITEMS_IN_EACH_DATA_RECORD=")):
      column_index = column_index_retriever(line)
    elif (line.startswith("!ITEM_H=")):
      miller_index_columns[0] = column_index.get(line)
    elif (line.startswith("!ITEM_K=")):
      miller_index_columns[1] = column_index.get(line)
    elif (line.startswith("!ITEM_L=")):
      miller_index_columns[2] = column_index.get(line)
    elif (line.startswith("!ITEM_IOBS=")):
      iobs_column = column_index.get(line)
    elif (line.startswith("!ITEM_SIGMA(IOBS)=")):
      sigma_iobs_column = column_index.get(line)
    elif (line.startswith("!END_OF_HEADER")):
      break
  assert unit_cell is not None
  assert space_group_number is not None
  assert column_index is not None
  assert None not in miller_index_columns
  assert iobs_column is not None
  assert sigma_iobs_column is not None
  miller_indices = flex.miller_index()
  iobs = flex.double()
  sigma_iobs = flex.double()
  for line in f:
    if (line.startswith("!END_OF_DATA")):
      break
    data = line.split()
    assert len(data) == column_index.number_of_items_in_each_data_record
    h = [int(data[i]) for i in miller_index_columns]
    miller_indices.append(h)
    iobs.append(float(data[iobs_column]))
    sigma_iobs.append(float(data[sigma_iobs_column]))
  miller_array = miller.array(
    miller_set=miller.set(
      crystal_symmetry=crystal.symmetry(
        unit_cell=unit_cell,
        space_group_symbol=space_group_number),
      indices=miller_indices,
      anomalous_flag=anomalous_flag),
    data=iobs,
    sigmas=sigma_iobs).set_observation_type_xray_intensity()
  miller_array.show_comprehensive_summary()

if (__name__ == "__main__"):
  read(sys.argv[1])
