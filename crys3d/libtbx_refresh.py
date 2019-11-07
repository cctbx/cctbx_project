from __future__ import absolute_import, division, print_function
try:
  from PyQt4 import uic
except ImportError:
  pass
else:
  import libtbx.load_env
  ui_dir = libtbx.env.under_dist(module_name="crys3d", path="qttbx")
  print('  Processing *.ui files in "%s"' % ui_dir)
  try:
    # check in case PyQt4 looks for ucs2, our Python is built with ucs4
    uic.compileUiDir(ui_dir, recurse=True)
  except BaseException:
    pass
