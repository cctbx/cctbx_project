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
  create_target_dir = True
  for file_name in names:
    name = file_name.lower()
    if (   name == "libtbx_config"
        or (name == "dispatcher_front_end.exe" and os.name == "nt")
        or name.startswith("authors")
        or name.startswith("copying")
        or name.startswith("copyright")
        or name.startswith("license")
        or name == "academic_software_licence.pdf"
        or name == "symop.lib"
        or name.endswith(".py")
        or name.endswith(".params")
        or name.endswith(".html")
        or name.endswith(".csh")
        or name.endswith(".sh")):
      src = libtbx.path.norm_join(dirname, file_name)
      if (not os.path.isdir(src)):
        dest = libtbx.path.norm_join(dist_copy, src)
        if (create_target_dir):
          libtbx.path.create_target_dir(dest)
          create_target_dir = False
        shutil.copy(src, dest)

def run(target_root):
  cwd = os.getcwd()
  abs_target_root = os.path.normpath(os.path.abspath(os.path.join(
    target_root)))
  libtbx_env = libtbx.config.env()
  dist_root = libtbx_env.LIBTBX_DIST_ROOT
  for package in ["libtbx"] + libtbx_env.package_list:
    for package_suf in libtbx.config.package_pair(
                         name=package).primary_first():
      if (package_suf == "boost"):
        continue
      dist = libtbx.config.resolve_redirection(
        dist_root=dist_root,
        name=package_suf).dist_path
      if (os.path.isdir(dist)):
        dist_copy = libtbx.path.norm_join(abs_target_root, package_suf)
        os.chdir(dist)
        os.path.walk(".", copy_dist_files, dist_copy)
  libtbx.bundle.utils.write_bundle_info(dist_root, abs_target_root)
  os.chdir(cwd)

if (__name__ == "__main__"):
  assert len(sys.argv) == 2
  run(sys.argv[1])
