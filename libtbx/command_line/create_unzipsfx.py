import libtbx.path
import sys

buf_size = 1000000

def copy(src, dest):
  while 1:
    buf = src.read(buf_size)
    if (buf == ""): break
    dest.write(buf)

def find_unzipsfx():
  for command in ("unzipsfx_autorun_yes.exe",
                  "unzipsfx_autorun.exe",
                  "unzipsfx.exe"):
    path_cmd = libtbx.path.full_command_path(
      command=command, search_first=["."])
    if (path_cmd is not None): return path_cmd
  return None

def create(zip_file_name, path_unzipsfx_exe=None):
  if (path_unzipsfx_exe is None):
    path_unzipsfx_exe = find_unzipsfx()
  if (path_unzipsfx_exe is None):
    raise RuntimeError("Fatal: unzipsfx executable not found.")
  assert zip_file_name.endswith(".zip")
  exe_file_name = zip_file_name[:-4] + ".exe"
  exe_file = open(exe_file_name, "wb")
  copy(open(path_unzipsfx_exe, "rb"), exe_file)
  copy(open(zip_file_name, "rb"), exe_file)
  exe_file.close()

def run(args):
  "usage: libtbx.create_unzipsfx [path_unzipsfx_exe] zip_file_name"
  if (not len(args) in (1,2) or "-h" in args or "--help" in args):
    print run.__doc__
    return
  if (len(args) == 1):
    create(zip_file_name=args[0])
  else:
    create(zip_file_name=args[1], path_unzipsfx_exe=args[0])

if (__name__ == "__main__"):
  run(sys.argv[1:])
