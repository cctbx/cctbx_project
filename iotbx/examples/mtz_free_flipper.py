"""Flip test set and working set"""
from __future__ import absolute_import, division, print_function
import iotbx.mtz
import sys, os

def run(args, label="R-free-flags"):
  assert len(args) == 1
  input_file_name = args[0]
  output_file_name = "free_flipped_"+os.path.basename(input_file_name)
  print("Reading file:", input_file_name)
  mtz_obj = iotbx.mtz.object(file_name=input_file_name)
  column = mtz_obj.get_column(label=label)
  selection_valid = column.selection_valid()
  flags = column.extract_values()
  sel_0 = (flags == 0)
  print("Number of 0:", ( sel_0 & selection_valid).count(True))
  print("Number of 1:", (~sel_0 & selection_valid).count(True))
  flags.set_selected( sel_0 & selection_valid, 1)
  flags.set_selected(~sel_0 & selection_valid, 0)
  column.set_values(values=flags, selection_valid=selection_valid)
  print("Writing file:", output_file_name)
  mtz_obj.write(file_name=output_file_name)

if (__name__ == "__main__"):
  run(sys.argv[1:])
