try:
  from PyQt4 import uic
except ImportError:
  pass
else:
  import libtbx.load_env
  ui_dir = libtbx.env.under_dist(module_name="crys3d", path="qttbx")
  print '  Processing *.ui files in "%s"' % ui_dir
  uic.compileUiDir(ui_dir, recurse=True)
