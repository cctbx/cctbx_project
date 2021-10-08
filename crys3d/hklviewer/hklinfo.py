
from __future__ import absolute_import, division, print_function

import sys, time
from iotbx.data_manager import DataManager
from iotbx.gui_tools.reflections import ArrayInfo
import textwrap


def run(arrays, philparams=None, log = sys.stdout):
  column_names_selection = []
  if philparams:
    wrap_labels = philparams.wrap_labels
    for att,val in list(philparams.selected_info.__dict__.items()):
      if not att.startswith("__"):
        column_names_selection.append( (att, val) )
  else:
    column_names_selection = None
    wrap_labels = 0
  delimiter = philparams.delimiter
  array_info_format_tpl=[]
  for i,array in enumerate(arrays):
    arrayinfo = ArrayInfo(array, wrap_labels)
    array_info_format_tpl.append( arrayinfo.get_selected_info_columns_from_phil(philparams))

  print("%d Miller arrays in this dataset:" %len(arrays))
  for i,info_fmt in enumerate(array_info_format_tpl):
    headerlst, infolst, dummy, fmtlst = info_fmt
    if i==0:
      for h in headerlst:
        print(h, end=delimiter)
      print("") # print that line break
    for i,info in enumerate(infolst):
      inf = info
      if i==0 and wrap_labels>0:
        inf = textwrap.wrap(info, width=wrap_labels)
      print(fmtlst[i].format(inf), end=delimiter) # no line break
    print("") # print that line break


if (__name__ == "__main__") :
  #time.sleep(10) # enough time to attach debugger
  dm = DataManager(datatypes=["miller_array"])
  args = sys.argv[1:]
  if (len(args)) != 1:  # print how to run program and quit
    print ("Usage: 'phenix.HKLinfo <data_file>'")
  else:
    data_file = args[0]
    arrays = dm.get_miller_arrays(filename = data_file)
    run(arrays)

