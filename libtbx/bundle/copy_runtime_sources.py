import libtbx.bundle.utils
import libtbx.load_env
import libtbx.path
import re
import shutil
import sys, os

def copy_dist_files((exclude_from_binary_bundle, dist_copy), dirname, names):
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
    if (name.startswith(".")): continue
    if (name == "sconscript"): continue
    if (name.startswith("makefile")): continue
    if (name.endswith(".exe") and os.name != "nt"): continue
    for ext in [".pyc", ".pyo", ".h", ".c", ".hpp", ".cpp", ".cc", ".f"]:
      if (name.endswith(ext)):
        break
    else:
      src = libtbx.path.norm_join(dirname, file_name)
      if (os.path.isdir(src)): continue
      for pattern in exclude_from_binary_bundle:
        if (re.search(pattern, src) is not None):
          break
      else:
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
  abs_target_root = os.path.normpath(os.path.abspath(target_root))
  for module in libtbx.env.module_list:
    for name,dist_path in module.name_and_dist_path_pairs():
      if (name == "boost"): continue
      dist_path = abs(dist_path)
      dist_copy = libtbx.path.norm_join(
        abs_target_root, os.path.basename(dist_path))
      os.chdir(dist_path)
      os.path.walk(
        ".", copy_dist_files, (module.exclude_from_binary_bundle, dist_copy))
  libtbx.bundle.utils.write_bundle_info(abs_target_root)
  os.chdir(cwd)

if (__name__ == "__main__"):
  assert len(sys.argv) == 2
  run(sys.argv[1])
