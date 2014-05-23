import os
import re
import shutil
import subprocess
import glob
import datetime
import argparse

# Binary builds: /net/cci/auto_build/phenix_installers
# export DYLD_PRINT_INITIALIZERS="files"
##### Helper functions #####

# For Linux, requires PatchELF.
#   http://nixos.org/patchelf.html
# For Mac, requires otool and install_name_tool,
#   which are included with XCode.

def find_exec(root='.'):
  """Find executables (using +x permissions)."""
  # find . -type f -perm +111 -print
  p = check_output(['find', root, '-type', 'f', '-perm', '+111'])
  # unix find may print empty lines; strip those out.
  return filter(None, [i.strip() for i in p.split("\n")])

def find_ext(ext='', root='.'):
  """Find files with a particular extension. Include the ".", e.g. ".txt". """
  found = []
  for root, dirs, files in os.walk(root):
    found.extend([os.path.join(root, i) for i in files if i.endswith(ext)])
  return found

def cmd(*popenargs, **kwargs):
  print "Running:", 
  print " ".join(*popenargs)
  kwargs['stdout'] = subprocess.PIPE
  kwargs['stderr'] = subprocess.PIPE
  process = subprocess.Popen(*popenargs, **kwargs)
  a, b = process.communicate()
  exitcode = process.wait()
  if exitcode:
    print("WARNING: Command returned non-zero exit code: %s"%" ".join(*popenargs))
    print a
    print b

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
  return output

class FixLinuxRpath(object):
  def run(self, root, replace=None):
    replace = replace or {}
    targets = set()
    targets |= set(find_ext('.so', root=root))
    targets |= set(find_ext('.dylib', root=root))
    # targets |= set(find_exec(root=root))
    for target in sorted(targets):
      if ".py" in target:
        continue
      xtarget = target.replace(root, '')
      depth = len(xtarget.split('/'))-2
      origins = ['$ORIGIN/']
      base = "".join(["../"]*depth)
      for i in ['extlib/lib']:
        origins.append('$ORIGIN/'+base+i+'/')
      try:
        cmd(['patchelf', '--set-rpath', ":".join(origins), target])
      except Exception, e:
        print "Couldnt patchelf:", e    
    
class FixMacRpath(object):
  """Process all binary files (executables, libraries) to rename linked libraries."""

  def find_deps(self, filename):
    """Find linked libraries using otool -L."""
    p = check_output(['otool','-L',filename])
    # otool doesn't return an exit code on failure, so check..
    if "not an object file" in p:
      raise Exception, "Not Mach-O binary"
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
    # Find all files that end in .so/.dylib, or are executable
    # This will include many script files, but we will ignore
    # these failures when running otool/install_name_tool
    targets = set()
    targets |= set(find_ext('.so', root=root))
    targets |= set(find_ext('.dylib', root=root))
    print "Targets:", len(targets)
    
    # targets |= set(find_exec(root=root))
    for f in sorted(targets):
      # Get the linked libraries and
      # check if the file is a Mach-O binary
      print "\n==== Target:", f
      try:
        libs = self.find_deps(f)
      except Exception, e:
        continue
        
      # Strip the absolute path down to a relative path
      frel = f.replace(root, "")[1:]

      # Set the install_name.
      install_name_id = os.path.join('@rpath', frel)
      cmd(['install_name_tool', '-id', install_name_id, f], cwd=root)

      # Set @rpath, this is a reference to the root of the package.
      # Linked libraries will be referenced relative to this.
      rpath = self.id_rpath(frel)
      try:
        cmd(['install_name_tool', '-add_rpath', rpath, f], cwd=root)
      except:
        pass

      # Process each linked library with the regexes in REPLACE.
      for lib in libs:
        olib = lib
        for k,v in replace.items():
          lib = re.sub(k, v, lib)
        if olib != lib:
          try:
            cmd(['install_name_tool', '-change', olib, lib, f], cwd=root)
          except:
            pass

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("root", help="Build path")
  args = parser.parse_args()
  
  # For Mac.
  replace = {}
  replace['^/scratch/phenix/(.+?)/'] = '@rpath/'
  replace['^lib/'] = '@rpath/build/mac-intel-osx-x86_64/lib/'
  FixMacRpath().run(root=args.root, replace=replace)
  
  