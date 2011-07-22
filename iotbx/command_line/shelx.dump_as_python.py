from cctbx import xray

def run(file_name, exclude_hydrogens=False):
  xs = xray.structure.from_shelx(filename=file_name,
                                 strictly_shelxl=False)
  if exclude_hydrogens:
    xs = xs.select(xs.hd_selection(), negate=True)
  print xs.as_py_code()

if __name__ == '__main__':
  import sys
  run(sys.argv[1], '--exclude-hydrogens' in sys.argv[2:])
