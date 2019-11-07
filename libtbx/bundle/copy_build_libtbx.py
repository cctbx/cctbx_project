from __future__ import absolute_import, division, print_function
import libtbx.bundle.utils
import libtbx.load_env
import libtbx.path
import shutil
import sys, os

def copy_lib_and_exe_files(target_root_dir, dirname, names):
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
    dest = os.path.normpath(os.path.join(target_root_dir, src))
    if (create_target_dir):
      libtbx.path.create_target_dir(dest)
      create_target_dir = False
    shutil.copy(src, dest)

def copy_base_files(target_root_dir, dirname, names):
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
  def copy_sub_dir(sub_dir, visitor):
    source_dir = libtbx.env.under_build(sub_dir)
    if (os.path.isdir(source_dir)):
      target_dir = libtbx.path.norm_join(abs_target_root, sub_dir)
      os.chdir(source_dir)
      for root, dirs, files in os.walk("."):
        visitor(target_dir, root, files)
  for sub_dir,visitor in (("lib", copy_lib_and_exe_files),
                          ("exe", copy_lib_and_exe_files),
                          ("base", copy_base_files)):
    copy_sub_dir(sub_dir=sub_dir, visitor=visitor)
  for module in libtbx.env.module_list:
    for name in module.names:
      copy_sub_dir(
        sub_dir=os.path.join(name, "exe"),
        visitor=copy_lib_and_exe_files)
  libtbx.bundle.utils.write_bundle_info(
    abs_target_root, write_build_options=True)
  file_name = "command_version_suffix"
  src = libtbx.env.under_build(file_name)
  if (os.path.isfile(src)):
    shutil.copy(src, libtbx.path.norm_join(abs_target_root, file_name))
  os.chdir(cwd)

if (__name__ == "__main__"):
  assert len(sys.argv) == 2
  run(sys.argv[1])
