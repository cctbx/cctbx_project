import libtbx.bundle.utils
import libtbx.load_env
import libtbx.path
import shutil
import sys, os

def copy_dist_files(dist_copy, dirname, names):
  create_target_dir = True
  names_keep = []
  for file_name in names:
    name = file_name.lower()
    def name_is_sub_dir_with_file(sub_dir, file_in_subdir):
      if (sub_dir != name): return False
      path = os.path.join(dirname, name)
      if (not os.path.isdir(path)): return False
      path = os.path.join(path, file_in_subdir)
      return os.path.isfile(path)
    if (   name_is_sub_dir_with_file("cvs", "Entries")
        or name_is_sub_dir_with_file(".svn", "README.txt")
        or name_is_sub_dir_with_file(".svn", "entries")):
      continue
    names_keep.append(file_name)
    if (   name == "libtbx_config"
        or (name == "windows_dispatcher.exe" and os.name == "nt")
        or name.startswith("authors")
        or name.startswith("copying")
        or name.startswith("copyright")
        or name.startswith("license")
        or name == "cci_diffs"
        or name == "academic_software_licence.pdf"
        or name == "symop.lib"
        or name == "case_library"
        or name.endswith(".py")
        or name.endswith(".params")
        or name.endswith(".pdb")
        or name.endswith(".pl")
        or name.endswith(".pm")
        or name.endswith(".scm")
        or name.endswith(".html")
        or name.endswith(".txt")
        or name.endswith(".csh")
        or name.endswith(".sh")):
      src = libtbx.path.norm_join(dirname, file_name)
      if (not os.path.isdir(src)):
        dest = libtbx.path.norm_join(dist_copy, src)
        if (create_target_dir):
          libtbx.path.create_target_dir(dest)
          create_target_dir = False
        shutil.copy(src, dest)
  if (len(names_keep) != len(names)):
    del names[:]
    names.extend(names_keep)

def run(target_root):
  cwd = os.getcwd()
  abs_target_root = libtbx.env.abs_path_clean(target_root)
  for module in libtbx.env.module_list:
    for name,dist_path in module.name_and_dist_path_pairs():
      if (name == "boost"): continue
      dist_copy = libtbx.path.norm_join(
        abs_target_root, os.path.basename(dist_path))
      os.chdir(dist_path)
      os.path.walk(".", copy_dist_files, dist_copy)
  libtbx.bundle.utils.write_bundle_info(abs_target_root)
  os.chdir(cwd)

if (__name__ == "__main__"):
  assert len(sys.argv) == 2
  run(sys.argv[1])
