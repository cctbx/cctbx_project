try:
  from PyQt4 import uic
  from crys3d import qttbx
  print "Processing *.ui files in qttbx"
  uic.compileUiDir(qttbx.__path__[0], recurse=True)
except ImportError:
  pass
