#! /usr/bin/env python

import libtbx.config
import sys, os
from os.path import normpath, join, split, isdir, isfile, splitext
norm = normpath
import shutil

def create_driver(target_dir, package_name, source_dir, file_name):
  source_file = norm(join(source_dir, file_name))
  if (not isfile(source_file)): return
  if (file_name.lower().startswith("__init__.py")): return
  if (file_name.lower().endswith(".pyc")): return
  if (file_name[0] == "."): return
  target_file = norm(join(target_dir, package_name))
  if (file_name.lower() != "main.py"):
    target_file += "." + splitext(file_name)[0]
  if (os.name == "nt"):
    if (not file_name.lower().endswith(".py")): return
    target_file += ".py"
    action = shutil.copyfile
  else:
    action = os.symlink
    try: os.chmod(source_file, 0755)
    except: pass
  if (os.path.isfile(target_file)):
    try: os.remove(target_file)
    except OSError: pass
    else: action(source_file, target_file)
  else:
    action(source_file, target_file)

def create_drivers(target_dir, package_name, source_dir):
  if (not isdir(source_dir)): return
  print "Processing:", source_dir
  for file_name in os.listdir(source_dir):
    create_driver(target_dir, package_name, source_dir, file_name)

def create_python_dispatchers(target_dir, python_exe):
  for file_name in ("libtbx.python", "python"):
    target_file = norm(join(target_dir, file_name))
    if (os.name == "nt"):
      target_file += ".exe"
      action = shutil.copyfile
    else:
      action = os.symlink
    if (os.path.isfile(target_file)):
      try: os.remove(target_file)
      except OSError: pass
      else: action(python_exe, target_file)
    else:
      action(python_exe, target_file)

def run():
  libtbx_env = libtbx.config.env()
  target_dir = norm(join(libtbx_env.LIBTBX_BUILD, "libtbx/bin"))
  if (not isdir(target_dir)):
    os.makedirs(target_dir)
  create_python_dispatchers(target_dir, libtbx_env.LIBTBX_PYTHON_EXE)
  for dist_path in libtbx_env.dist_paths.values():
    package_name = split(dist_path)[1]
    for suffix in ("", "_adaptbx"):
      create_drivers(
        target_dir,
        package_name,
        source_dir=norm(join(dist_path+suffix, package_name, "command_line")))

if (__name__ == "__main__"):
  run()
