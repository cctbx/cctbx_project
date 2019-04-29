from __future__ import division, print_function
# LIBTBX_SET_DISPATCHER_NAME prime.change_basis

import argparse
from iotbx import reflection_file_reader
from cctbx import sgtbx
import os

def main(mtz, index_basis):
  # read in mtz
  reflection_file = reflection_file_reader.any_reflection_file(mtz)
  miller_arrays = reflection_file.as_miller_arrays()
  # apply new basis
  cb_op = sgtbx.change_of_basis_op(index_basis)
  miller_array_new = miller_arrays[0].change_basis(cb_op)
  # write out new mtz
  hklout = os.path.splitext(mtz)[0]+'_modified.mtz'
  mtz_dataset = miller_array_new.as_mtz_dataset(column_root_label="IOBS")
  mtz_dataset.mtz_object().write(file_name=hklout)
  print('Output file saved to', hklout)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(
    description='Convert an mtz file to the specified basis.'
  )
  parser.add_argument(
    'mtz',
    metavar='MTZ File',
    help='Path to an mtz file'
  )
  parser.add_argument(
    'index_basis',
    metavar='New Basis',
    help='New index basis operators'
  )
  args = parser.parse_args()
  print('Convert ', args.mtz, ' to ', args.index_basis.strip())
  main(args.mtz, args.index_basis.strip())
