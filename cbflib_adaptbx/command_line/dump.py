from __future__ import absolute_import, division, print_function
from six.moves import range
# LIBTBX_SET_DISPATCHER_NAME cbf.dump

import sys

def process(file_name, out=None):
  if (out is None): out = sys.stdout
  import pycbf
  print("File name:", file_name, file=out)
  object = pycbf.cbf_handle_struct() # FIXME
  object.read_file(file_name, pycbf.MSG_DIGEST)
  object.rewind_datablock()
  n_blocks = object.count_datablocks()
  print("Number of blocks:", n_blocks, file=out)
  for i_block in range(n_blocks):
    object.select_datablock(i_block)
    print("  Block name:", object.datablock_name(), file=out)
    object.rewind_category()
    n_categories = object.count_categories()
    print("  Number of categories:", n_categories, file=out)
    for i_category in range(n_categories):
      object.select_category(i_category)
      print("    Category name:", object.category_name(), file=out)
  print(file=out)

def run(args):
  if (len(args) == 0):
    from libtbx.utils import Usage
    import libtbx.load_env
    raise Usage("%s your.cbf ..." % libtbx.env.dispatcher_name)
  for file_name in args:
    process(file_name)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
