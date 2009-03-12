import sys
import libtbx.command_line.file_clutter

def run(args):
  flag_x = False
  paths = []
  for arg in args:
    if (arg == "-x"):
      flag_x = True
    else:
      paths.append(arg)
  if (len(paths) == 0): paths = ["."]
  for cluttered_file in libtbx.command_line.file_clutter.gather(paths):
    cluttered_file.show(flag_x)

if (__name__ == "__main__"):
  run(sys.argv[1:])
