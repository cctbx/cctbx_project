import libtbx.bundle.utils
import libtbx.load_env
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
        or (name == "windows_dispatcher.exe" and os.name == "nt")
        or name.startswith("authors")
        or name.startswith("copying")
        or name.startswith("copyright")
        or name.startswith("license")
        or name == "academic_software_licence.pdf"
        or name == "symop.lib"
        or name.endswith(".py")
        or name.endswith(".params")
        or name.endswith(".pl")
        or name.endswith(".pm")
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
  abs_target_root = libtbx.env.abs_path_clean(target_root)
  for module in libtbx.env.module_list:
    for name,dist_path in module.name_and_dist_path_pairs():
      if (name == "boost"): continue
      dist_copy = libtbx.path.norm_join(abs_target_root, name)
      os.chdir(dist_path)
      os.path.walk(".", copy_dist_files, dist_copy)
  libtbx.bundle.utils.write_bundle_info(abs_target_root)
  os.chdir(cwd)

if (__name__ == "__main__"):
  assert len(sys.argv) == 2
  run(sys.argv[1])
