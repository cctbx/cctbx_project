
from libtbx import easy_pickle
import wx
import re
import os

def find_integration_files (dir_name, base_name) :
  files = []
  for file_name in os.listdir(dir_name) :
    if (file_name.startswith(base_name)) :
      files.append(os.path.join(dir_name, file_name))
  return files

def load_integration_results (dir_name, base_name) :
  files = find_integration_files(dir_name, base_name)
  results = []
  summaries = []
  for file_path in files :
    result = easy_pickle.load(file_path)
    results.append(result)
    file_name = os.path.basename(file_path)
    x = int(re.sub("_.*", "", re.sub(base_name + "_", "", file_name)))
    summary = dict(
      solution=x,
      point_group=result['pointgroup'],
      beam_center=(result['xbeam'], result['ybeam']),
      distance=result['distance'],
      d_min=result['resolution'],
      mosaicity=result['mosaicity'],
      rms=result['residual'],
      bins=result['table_raw'])
    summaries.append(summary)
  r_s = list(zip(results, summaries))
  r_s_sorted = sorted(r_s, lambda x,y: cmp(y[1]['solution'], x[1]['solution']))
  return [ r for r,s in r_s_sorted ], [ s for r, s in r_s_sorted ]

class TableData (object) :
  """Base class for wx.ListCtrl data source objects in this module."""
  def __init__ (self, table) :
    assert isinstance(table, list)
    self.table = table

  def GetItemCount (self) :
    return len(self.table)

  def GetItemImage (self, item) :
    return 0

class ResultData (TableData) :
  def GetItemText (self, item, col) :
    n_items = self.GetItemCount()
    assert (item < n_items) and (0 <= col <= 6)
    result = self.table[item]
    if (col == 0) :
      return "%d" % result['solution']
    elif (col == 1) :
      return result['point_group']
    elif (col == 2) :
      return "%.2f %.2f" % result['beam_center']
    elif (col == 3) :
      return "%.2f" % result['distance']
    elif (col == 4) :
      return "%.2f" % result['d_min']
    elif (col == 5) :
      return "%.2f" % result['mosaicity']
    else :
      return "%.3f" % result['rms']

class BinData (TableData) :
  def GetItemText (self, item, col) :
    n_items = self.GetItemCount()
    assert (item < n_items) and (0 <= col <= 4)
    bin = self.table[item]
    if (col == 0) :
      return "%d" % bin.i_bin
    elif (col == 1) :
      return "%g - %g" % bin.d_max_min
    elif (col == 2) :
      return "%d / %d" % bin.completeness
    elif (col == 3) :
      return "%8.1f" % bin.mean_I
    else :
      return "%8.1f" % bin.mean_I_sigI
