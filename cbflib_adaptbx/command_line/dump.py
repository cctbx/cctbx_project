# LIBTBX_SET_DISPATCHER_NAME cbf.dump

import sys

def process(file_name, out=None):
  if (out is None): out = sys.stdout
  import pycbf
  print >> out, "File name:", file_name
  object = pycbf.cbf_handle_struct() # FIXME
  object.read_file(file_name, pycbf.MSG_DIGEST)
  object.rewind_datablock()
  n_blocks = object.count_datablocks()
  print >> out, "Number of blocks:", n_blocks
  for i_block in xrange(n_blocks):
    object.select_datablock(i_block)
    print >> out, "  Block name:", object.datablock_name()
    object.rewind_category()
    n_categories = object.count_categories()
    print >> out, "  Number of categories:", n_categories
    for i_category in xrange(n_categories):
      object.select_category(i_category)
      print >> out, "    Category name:", object.category_name()
  print >> out

def run(args):
  if (len(args) == 0):
    from libtbx.utils import Usage
    import libtbx.load_env
    raise Usage("%s your.cbf ..." % libtbx.env.dispatcher_name)
  for file_name in args:
    process(file_name)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
