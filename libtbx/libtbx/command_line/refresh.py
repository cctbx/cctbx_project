#! /usr/bin/env python

import libtbx.config
import sys, os
from os.path import normpath, join, split, isdir, isfile, islink, splitext
norm = normpath
import shutil

class create_bin_sh_dispatcher:

  def __init__(self, python_exe, precall_commands):
    self.python_exe = python_exe
    self.precall_commands = precall_commands

  def __call__(self, source_file, target_file):
    f = open(target_file, "w")
    print >> f, "#! /bin/sh"
    print >> f, "# LIBTBX_DISPATCHER DO NOT EDIT"
    print >> f, "unset PYTHONHOME"
    if (self.precall_commands is not None):
      for line in self.precall_commands:
        print >> f, line
    cmd = "  exec"
    if (    source_file != self.python_exe
        and source_file.lower().endswith(".py")):
      cmd += " '"+self.python_exe+"'"
    print >> f, "if [ $# -eq 0 ]; then"
    print >> f, cmd, "'"+source_file+"'"
    print >> f, "else"
    print >> f, cmd, "'"+source_file+"'", '"$@"'
    print >> f, "fi"
    f.close()
    os.chmod(target_file, 0755)

def create_python_execfile_dispatcher(source_file, target_file):
  f = open(target_file, "w")
  print >> f, 'execfile(r"%s")' % source_file
  f.close()

def create_driver(python_exe, precall_commands,
                  target_dir, package_name, source_dir, file_name):
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
    target_file += ".px"
    action = create_python_execfile_dispatcher
  else:
    action = create_bin_sh_dispatcher(
      python_exe=python_exe,
      precall_commands=precall_commands)
    try: os.chmod(source_file, 0755)
    except: pass
  if (isfile(target_file) or islink(target_file)):
    try: os.remove(target_file)
    except OSError: pass
    else: action(source_file, target_file)
  else:
    action(source_file, target_file)

def create_drivers(python_exe, precall_commands,
                   target_dir, package_name, source_dir):
  if (not isdir(source_dir)): return
  print "Processing:", source_dir
  for file_name in os.listdir(source_dir):
    create_driver(
      python_exe=python_exe,
      precall_commands=precall_commands,
      target_dir=target_dir,
      package_name=package_name,
      source_dir=source_dir,
      file_name=file_name)

def create_posix_icc_ld_preload(libtbx_env):
  path_icc = libtbx.config.full_path("icc")
  if (path_icc is None): return None
  path_lib = os.sep.join(path_icc.split(os.sep)[:-2] + ["lib"])
  if (not os.path.isdir(path_lib)): return None
  ld_preload = []
  path_libirc_a = os.path.join(path_lib, "libirc.a")
  path_libirc_so = os.path.join(path_lib, "libirc.so")
  if (os.path.isfile(path_libirc_so)):
    ld_preload.append(path_libirc_so)
  else:
    if (os.path.isfile(path_libirc_a)):
      path_libirc_so = os.path.join(
        libtbx_env.LIBTBX_BUILD, "libtbx", "libirc.so")
      cmd = "%(path_icc)s -shared -o %(path_libirc_so)s %(path_libirc_a)s" \
        % vars()
      print cmd
      sys.stdout.flush()
      os.system(cmd)
      ld_preload.append(path_libirc_so)
  path_libunwind_so = None
  best_version = None
  for file_name in os.listdir(path_lib):
    if (file_name.startswith("libunwind.so.")):
      try: version = int(file_name.split(".")[2])
      except: version = None
      if (version is not None):
        if (best_version is None or version > best_version):
          path_libunwind_so = os.path.join(path_lib, file_name)
          best_version = version
  if (path_libunwind_so is not None):
    ld_preload.append(path_libunwind_so)
  if (len(ld_preload) == 0): return None
  return [
    'LD_PRELOAD="%s"' % os.pathsep.join(ld_preload),
    'export LD_PRELOAD']

def assemble_dispatcher_precall_commands(libtbx_env):
  lines = []
  if (    libtbx_env.python_version_major_minor == (2,2)
      and sys.platform == "linux2"
      and os.path.isfile("/etc/redhat-release")):
    try: red_hat_linux_release = open("/etc/redhat-release").readline()
    except: pass
    else:
      if (    red_hat_linux_release.startswith("Red Hat Linux release")
          and red_hat_linux_release.split()[4] == "9"):
        lines.extend([
          'if [ ! -n "$LD_ASSUME_KERNEL" ]; then',
          '  LD_ASSUME_KERNEL=2.4.1',
          '  export LD_ASSUME_KERNEL',
          'fi'])
  if (os.name == "posix" and libtbx_env.compiler == "icc"):
    addl_lines = create_posix_icc_ld_preload(libtbx_env)
    if (addl_lines is None):
      raise libtbx.config.UserError("Cannot determine LD_PRELOAD for icc.")
    lines.extend(addl_lines)
  return lines

def create_python_dispatchers(libtbx_env, target_dir,
                              python_exe, precall_commands):
  for file_name in ("libtbx.python", "python"):
    target_file = norm(join(target_dir, file_name))
    if (os.name == "nt"):
      target_file += ".exe"
      action = shutil.copyfile
    else:
      action = create_bin_sh_dispatcher(
        python_exe=python_exe,
        precall_commands=precall_commands)
      try: os.chmod(source_file, 0755)
      except: pass
    if (isfile(target_file) or islink(target_file)):
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
  precall_commands = assemble_dispatcher_precall_commands(libtbx_env)
  create_python_dispatchers(
    libtbx_env=libtbx_env,
    target_dir=target_dir,
    python_exe=libtbx_env.LIBTBX_PYTHON_EXE,
    precall_commands=precall_commands)
  for dist_path in libtbx_env.dist_paths.values():
    package_name = os.path.basename(dist_path)
    for dist_path_suf in libtbx.config.package_pair(dist_path).primary_first():
      create_drivers(
        python_exe=libtbx_env.LIBTBX_PYTHON_EXE,
        precall_commands=precall_commands,
        target_dir=target_dir,
        package_name=package_name,
        source_dir=norm(join(dist_path_suf, package_name, "command_line")))

if (__name__ == "__main__"):
  run()
