import sys
import libtbx.command_line.file_clutter

def run(args):
  flag_x = False
  flag_dos_format = True
  paths = []
  for arg in args:
    if (arg == "-x"):
      flag_x = True
    elif (arg == "-ndos"):
      flag_dos_format = False
    else:
      paths.append(arg)
  if (len(paths) == 0): paths = ["."]
  for cluttered_file in libtbx.command_line.file_clutter.gather(paths):
    cluttered_file.show(flag_x, flag_dos_format)

if (__name__ == "__main__"):
  run(sys.argv[1:])
