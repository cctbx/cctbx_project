import libtbx.config
import sys, os
from os.path import normpath, join, split, isdir, isfile, splitext
norm = normpath

def create_driver(target_dir, package_name, source_dir, file_name):
  source_file = norm(join(source_dir, file_name))
  if (not isfile(source_file)): return
  if (file_name.lower().startswith("__init__.py")): return
  if (file_name[0] == "."): return
  target_file = norm(join(target_dir, package_name))
  if (file_name.lower() != "main.py"):
    target_file += "." + splitext(file_name)[0]
  if (os.name == "nt"):
    if (file_name.lower().endswith(".py")):
      target_file += ".bat"
      f = open(target_file, "w")
      print >> f, "@echo off"
      print >> f, "python", source_file, "%1 %2 %3 %4 %5 %6 %7 %8 %9"
      f.close()
  else:
    f = open(target_file, "w")
    cmd = "exec"
    if (file_name.lower().endswith(".py")):
      cmd += " python"
    print >> f, "#! /bin/csh -f"
    print >> f, "set noglob"
    print >> f, cmd, source_file, "$*"
    f.close()
    os.chmod(target_file, 0755)

def create_drivers(target_dir, package_name, source_dir):
  if (not isdir(source_dir)): return
  print "Creating drivers for files in:", source_dir
  for file_name in os.listdir(source_dir):
    create_driver(target_dir, package_name, source_dir, file_name)

def run():
  libtbx_env = libtbx.config.env()
  target_dir = norm(join(libtbx_env.LIBTBX_BUILD, "libtbx/bin"))
  if (not isdir(target_dir)):
    os.makedirs(target_dir)
  for dist_path in libtbx_env.dist_paths.values():
    package_name = split(dist_path)[1]
    create_drivers(
      target_dir,
      package_name,
      source_dir=norm(join(dist_path, package_name, "command_line")))

if (__name__ == "__main__"):
  run()
