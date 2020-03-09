"""Fix RPATH and ORIGIN for relocatable binaries."""
from __future__ import absolute_import, division, print_function

import optparse
import os
import re
import subprocess
import sys

import six

# This code is partially derived from a similar module I wrote
# for EMAN2, previously released under a BSD-like license.

# Binary builds: /net/cci/auto_build/phenix_installers
# export DYLD_PRINT_INITIALIZERS="files"
# For Linux, requires PatchELF.
#   http://nixos.org/patchelf.html
# For Mac, requires otool and install_name_tool,
#   which are included with XCode.

# Ugh
DRY_RUN = False

##### Helper functions #####

def find_exec(root='.'):
  """Find executables (using +x permissions)."""
  # find . -type f -perm +111 -print
  p = check_output(['find', root, '-type', 'f', '-perm', '+111'])
  # unix find may print empty lines; strip those out.
  found = filter(None, [i.strip() for i in p.split("\n")])
  # Filter by real execuable...
  found_exec = []
  for f in found:
    try:
      p = check_output(['file', f])
    except subprocess.CalledProcessError:
      p = ''
    if "Mach-O" in p:
      found_exec.append(f)
  return found_exec

def find_ext(ext='', root='.'):
  """Find files with a particular extension. Include the ".", e.g. ".txt". """
  found = []
  for root, dirs, files in os.walk(root):
    found.extend([os.path.join(root, i) for i in files if i.endswith(ext)])
  return found

def cmd(*popenargs, **kwargs):
  print("Running:", end=' ')
  print(" ".join(*popenargs))
  if DRY_RUN:
    return
  ignorefail = kwargs.pop('ignorefail', False)
  kwargs['stdout'] = subprocess.PIPE
  kwargs['stderr'] = subprocess.PIPE
  process = subprocess.Popen(*popenargs, **kwargs)
  a, b = process.communicate()
  exitcode = process.wait()
  if exitcode:
    if ignorefail:
      print("WARNING: Command returned non-zero exit code: %s"%" ".join(*popenargs))
      print(a)
      print(b)
    else:
      print(a)
      print(b)
      raise Exception("Command returned non-zero exit code")

def echo(*popenargs, **kwargs):
    print("Running:", end=' ')
    print(" ".join(*popenargs))

# "Dry-run"
# cmd = echo

def check_output(*popenargs, **kwargs):
  """Copy of subprocess.check_output()"""
  if 'stdout' in kwargs:
    raise ValueError('stdout argument not allowed, it will be overridden.')
  process = subprocess.Popen(stdout=subprocess.PIPE, *popenargs, **kwargs)
  output, unused_err = process.communicate()
  retcode = process.poll()
  if retcode:
    cmd = kwargs.get("args")
    if cmd is None:
      cmd = popenargs[0]
    raise subprocess.CalledProcessError(retcode, cmd)
  if six.PY3:
    return output.decode("latin-1")
  return output

class FixLinuxRpath(object):
  def find_deps(self, filename):
    ret = []
    p = check_output(['ldd', filename])
    for line in p.split("\n"):
      fields = line.partition("=>")
      try:
        lib, found = fields[0], fields[2]
        if "not found" in found:
          ret.append(lib.strip())
      except Exception as e:
        print("PARSER ERROR:", e)
        print(fields)
    return ret

  def run(self, root, replace=None):
    replace = replace or {}
    targets = set()
    targets |= set(find_ext('.so', root=root))
    # targets |= set(find_exec(root=root))

    # First pass: get all directories.
    origins = {}
    for target in sorted(targets):
      origins[os.path.basename(target)] = os.path.dirname(target)

    # Second pass: find linked libraries, add rpath's
    for target in sorted(targets):
      print("\n\n====", target)
      deps = self.find_deps(target)
      rpath = set()
      # Find a matching library...
      for dep in deps:
        found = [i for i in origins if i in dep]
        # print dep, found, map(origins.get, found)
        for i in found:
          relpath = os.path.relpath(origins[i], os.path.dirname(target))
          relpath = os.path.join('$ORIGIN', relpath)
          rpath.add(relpath)

      if rpath:
        cmd(['patchelf', '--set-rpath', ":".join(rpath), target])

class FixMacRpath(object):
  """Process all binary files (executables, libraries) to rename linked libraries."""

  def find_deps(self, filename):
    """Find linked libraries using otool -L."""
    p = check_output(['otool','-L',filename])
    # otool doesn't return an exit code on failure, so check..
    if "not an object file" in p:
      raise Exception("Not Mach-O binary")
    # Just get the dylib install names
    p = [i.strip().partition(" ")[0] for i in p.split("\n")[1:]]
    return p

  def id_rpath(self, filename):
    """Generate the @rpath for a file, relative to the current directory as @rpath root."""
    p = len(filename.split("/"))-1
    f = os.path.join("@loader_path", *[".."]*p)
    return f

  def run(self, root, replace=None):
    replace = replace or {}
    replace[root] = '@rpath'
    # Find all files that end in .so/.dylib
    targets = set()
    targets |= set(find_ext('.so', root=root))
    targets |= set(find_ext('.dylib', root=root))
    # Find all executable, binary (Mach-O) files.
    targets |= set(find_exec(root=root))

    print("Targets:", len(targets))
    for f in sorted(targets):
      # Get the linked libraries and
      # check if the file is a Mach-O binary
      print("\n==== Target:", f)
      try:
        libs = self.find_deps(f)
      except Exception:
        continue

      # Set the install_name id.
      install_name_id = os.path.join('@rpath', os.path.relpath(f, root))
      cmd(['install_name_tool', '-id', install_name_id, f], cwd=root, ignorefail=True)

      # Set @rpath, this is a reference to the root of the package.
      # Linked libraries will be referenced relative to this.
      rpath = os.path.join('@loader_path', os.path.relpath(root, os.path.dirname(f)))
      cmd(['install_name_tool', '-add_rpath', rpath, f], cwd=root, ignorefail=True)

      for lib in libs:
        rlib = lib
        for k,v in replace.items():
          rlib = re.sub(k, v, rlib)
        if lib != rlib:
          cmd(['install_name_tool', '-change', lib, rlib, f], cwd=root, ignorefail=True)

def run(args):
  parser = optparse.OptionParser()
  parser.add_option("--otherroot", help="Other build path")
  parser.add_option("--dry", help="Dry run", action="store_true")
  options, args = parser.parse_args(args)

  if options.dry:
    DRY_RUN = True

  # Setup args.
  root = os.path.abspath(args[-1]) # needs absolute path
  replace = {}
  replace['^lib'] = '@rpath/lib'
  if options.otherroot:
    replace[options.otherroot] = '@rpath'

  # Run the rpath fixer.
  cls = FixLinuxRpath
  if sys.platform == 'darwin':
    cls = FixMacRpath
  cls().run(root=root, replace=replace)

if __name__ == "__main__":
  run(sys.argv[1:])
