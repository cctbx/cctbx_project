import sys
from libtbx.command_line.file_clutter import gather

def run(args):
  flag_x = False
  flag_ni = False
  flag_dos_format = True
  paths = []
  for arg in args:
    if (arg == "-x"):
      flag_x = True
    elif (arg == "-ni"):
      flag_ni = True
    elif (arg == "-ndos"):
      flag_dos_format = False
    else:
      paths.append(arg)
  if (len(paths) == 0): paths = ["."]
  n_is_cluttered = 0
  n_has_unused_imports = 0
  first_print = True
  for info in gather(paths=paths, find_unused_imports=not flag_ni):
    if (info.is_cluttered(flag_x=flag_x)):
      n_is_cluttered += 1
      if (first_print):
        first_print = False
        print
    if (info.has_unused_imports()):
      n_has_unused_imports += 1
      if (first_print):
        first_print = False
        print
    info.show(flag_x, flag_dos_format)
  please_use = []
  if (n_is_cluttered != 0):
    please_use.append("libtbx.clean_clutter")
  if (n_has_unused_imports != 0):
    please_use.append("libtbx.find_unused_imports_crude")
  if (len(please_use) != 0):
    print
    print "*** To clean up please use: %s ***" % ", ".join(please_use)
    print

if (__name__ == "__main__"):
  run(sys.argv[1:])
