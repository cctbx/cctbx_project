
from __future__ import absolute_import, division, print_function

import sys
from iotbx.gui_tools.reflections import ArrayInfo


def run(arrays, philparams=None, log = sys.stdout):
  if philparams:
    wrap_labels = philparams.wrap_labels
  else:
    wrap_labels = 0
  print("%d Miller arrays in this dataset:" %len(arrays))
  delimiter = philparams.delimiter
  array_info_format_tpl=[]
  for i,array in enumerate(arrays):
    arrayinfo = ArrayInfo(array, wrap_labels)
    info_fmt, headerstr, infostr = arrayinfo.get_selected_info_columns_from_phil(philparams)
    if i==0:
      print(headerstr)
    print(infostr)

