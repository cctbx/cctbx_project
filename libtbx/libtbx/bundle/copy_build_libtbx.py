import libtbx.bundle.utils
import libtbx.config
import libtbx.path
import shutil
import sys, os

class copy_libtbx_files:

  def __init__(self, libtbx_env):
    self.dispatcher_front_end_exe = libtbx_env.dispatcher_front_end_exe()

  def __call__(self, target_root_dir, dirname, names):
    is_libtbx_bin = False
    if (os.path.split(dirname.lower()) == (".", "bin")):
      is_libtbx_bin = True
    create_target_dir = True
    for file_name in names:
      name = file_name.lower()
      if (   name == ".sconsign"
          or name.endswith(".pyc")
          or name.endswith(".lib")
          or name.endswith(".exp")
          or name.endswith(".a")):
        continue
      src = os.path.normpath(os.path.join(dirname, file_name))
      if (os.path.isdir(src)): continue
      if (is_libtbx_bin and self.is_disposable_file(name=name, path=src)):
        continue
      dest = os.path.normpath(os.path.join(target_root_dir, src))
      if (create_target_dir):
        libtbx.path.create_target_dir(dest)
        create_target_dir = False
      shutil.copy(src, dest)

  def is_disposable_file(self, name, path):
    if (    self.dispatcher_front_end_exe is not None
        and name.endswith(".exe")):
      if (open(path, "rb").read() == self.dispatcher_front_end_exe):
        return True
    else:
      f = open(path, "r")
      lines = []
      try:
        lines.append(f.readline().strip())
        lines.append(f.readline().strip())
      except:
        pass
      f.close()
      if (len(lines) == 2
          and lines[0][0] == "#"
          and lines[1] == "# LIBTBX_DISPATCHER DO NOT EDIT"):
        return True
    return False

def copy_python_files(target_root_dir, dirname, names):
  create_target_dir = True
  for file_name in names:
    name = file_name.lower()
    if (name.endswith(".pyc")):
      continue
    src = os.path.normpath(os.path.join(dirname, file_name))
    if (os.path.isdir(src)): continue
    dest = os.path.normpath(os.path.join(target_root_dir, src))
    if (create_target_dir):
      libtbx.path.create_target_dir(dest)
      create_target_dir = False
    shutil.copy(src, dest)

def run(target_root):
  cwd = os.getcwd()
  abs_target_root = os.path.normpath(os.path.abspath(target_root))
  libtbx_env = libtbx.config.env()
  build_root = libtbx_env.LIBTBX_BUILD
  for sub_dir,visitor in (("libtbx", copy_libtbx_files(libtbx_env=libtbx_env)),
                          ("python", copy_python_files)):
    source_dir = libtbx.path.norm_join(build_root, sub_dir)
    if (os.path.isdir(source_dir)):
      target_dir = libtbx.path.norm_join(abs_target_root, sub_dir)
      os.chdir(source_dir)
      os.path.walk(".", visitor, target_dir)
  libtbx.bundle.utils.write_bundle_info(build_root, abs_target_root)
  os.chdir(cwd)

if (__name__ == "__main__"):
  assert len(sys.argv) == 2
  run(sys.argv[1])
