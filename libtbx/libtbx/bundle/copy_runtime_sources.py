import libtbx.bundle.utils
import libtbx.config
import libtbx.path
import shutil
import sys, os

def copy_dist_files(dist_copy, dirname, names):
  if (dirname.lower().endswith("cvs")):
    cvs_entries = libtbx.path.norm_join(dirname, "Entries")
    if (os.path.isfile(cvs_entries)):
      return
  create_target_dir = 0001
  for file_name in names:
    name = file_name.lower()
    if (   name == "libtbx_config"
        or (name == "dispatcher_front_end.exe" and os.name == "nt")
        or name.startswith("authors")
        or name.startswith("copying")
        or name.startswith("copyright")
        or name.startswith("license")
        or name.endswith(".py")
        or name.endswith(".html")
        or name.endswith(".csh")
        or name.endswith(".sh")):
      src = libtbx.path.norm_join(dirname, file_name)
      dest = libtbx.path.norm_join(dist_copy, src)
      if (create_target_dir):
        libtbx.path.create_target_dir(dest)
        create_target_dir = 00000
      shutil.copy(src, dest)

def is_required(dist, package):
  p = libtbx.path.norm_join(dist, "libtbx_config")
  if (os.path.isfile(p)): return 0001
  p = libtbx.path.norm_join(dist, package)
  if (os.path.isdir(p)):
    p = libtbx.path.norm_join(p, "__init__.py")
    if (os.path.isfile(p)): return 0001
  return 00000

def run(target_root):
  cwd = os.getcwd()
  abs_target_root = os.path.normpath(os.path.abspath(os.path.join(target_root)))
  libtbx_env = libtbx.config.env()
  dist_root = libtbx_env.LIBTBX_DIST_ROOT
  for package in ["libtbx"] + libtbx_env.package_list:
    for package_suf in libtbx.config.package_pair(package).primary_first():
      dist = libtbx.path.norm_join(dist_root, package_suf)
      if (os.path.isdir(dist) and is_required(dist, package)):
        dist_copy = libtbx.path.norm_join(abs_target_root, package_suf)
        os.chdir(dist)
        os.path.walk(".", copy_dist_files, dist_copy)
  libtbx.bundle.utils.write_bundle_info(dist_root, abs_target_root)
  os.chdir(cwd)

if (__name__ == "__main__"):
  assert len(sys.argv) == 2
  run(sys.argv[1])
