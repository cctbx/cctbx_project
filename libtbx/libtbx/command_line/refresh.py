import libtbx.env
import libtbx.config
import os
from os.path import normpath, join, isdir, isfile, splitext
norm = normpath

def create_dispatcher(target_dir, package_name, source_dir, file_name):
  source_file = norm(join(source_dir, file_name))
  if (not isfile(source_file)): return
  if (file_name.lower().startswith("__init__.py")): return
  if (file_name.lower().endswith(".pyc")): return
  if (file_name[0] == "."): return
  if (os.name == "nt"):
    ext = splitext(file_name)[1].lower()
    if (ext not in libtbx.config.windows_pathext): return
  target_file = norm(join(target_dir, package_name))
  if (file_name.lower() != "main.py"):
    target_file += "." + splitext(file_name)[0]
  libtbx.env.cache.create_dispatcher(
    source_file=source_file,
    target_file=target_file)

def create_dispatchers(target_dir, package_name, source_dir):
  if (not isdir(source_dir)): return
  print "Processing:", source_dir
  for file_name in os.listdir(source_dir):
    create_dispatcher(
      target_dir=target_dir,
      package_name=package_name,
      source_dir=source_dir,
      file_name=file_name)

def create_show_path_duplicates(target_dir):
  package_names = {}
  for file_name in os.listdir(target_dir):
    if (file_name.startswith(".")): continue
    if (file_name.startswith("libtbx.")): continue
    if (file_name == "python"): continue
    package_names[file_name.split(".")[0]] = None
  package_names = package_names.keys()
  for command in ["show_build_path", "show_dist_paths"]:
    source_file = libtbx.env.under_dist(
      "libtbx", "libtbx/command_line/"+command+".py")
    for package_name in package_names:
      target_file = os.path.join(target_dir, package_name+"."+command)
      libtbx.env.cache.create_dispatcher(
        source_file=source_file,
        target_file=target_file)

def run():
  target_dir = libtbx.env.under_build("libtbx/bin")
  if (not isdir(target_dir)):
    os.makedirs(target_dir)
  for file_name in ("libtbx.python", "python"):
    libtbx.env.cache.create_dispatcher(
      source_file=libtbx.env.cache.LIBTBX_PYTHON_EXE,
      target_file=norm(join(target_dir, file_name)))
  for dist_path in libtbx.env.cache.dist_paths.values():
    dist_root = os.path.dirname(dist_path)
    package_name = os.path.basename(dist_path)
    for dist_path_suf in libtbx.config.package_pair(
                           name=package_name,
                           dist_root=dist_root).primary_first():
      create_dispatchers(
        target_dir=target_dir,
        package_name=package_name,
        source_dir=norm(join(dist_path_suf, package_name, "command_line")))
  exe_path = norm(join(libtbx.env.cache.LIBTBX_BUILD, "exe"))
  if (os.path.isdir(exe_path)):
    print "Processing:", exe_path
    for file_name in os.listdir(exe_path):
      if (file_name[0] == "."): continue
      libtbx.env.cache.create_dispatcher(
        source_file=norm(join(exe_path, file_name)),
        target_file=norm(join(target_dir, file_name)))
  create_show_path_duplicates(target_dir)

if (__name__ == "__main__"):
  run()
