from __future__ import absolute_import, division, print_function

import collections
import errno
import glob
import os
import shutil
import subprocess
import sys
from optparse import SUPPRESS_HELP, OptionParser
from six.moves import range

import libtbx.load_env
from libtbx.auto_build.bootstrap import Toolbox

# Basically 'pip' for selected libtbx/cctbx modules.

class EpilogParser(OptionParser):
  """Simple, small OptionParser subclass to not strip epilog"""
  def format_epilog(self, formatter):
    """Don't strip newlines from this"""
    return self.epilog

def is_source_repository(path):
  return (path / '.git').isdir() or (path / '.svn').isdir()

def run(args):
  # Generate an epilog message
  possible_installs = "\nAvailable Packages:\n    " + "\n    ".join(
    x for x in sorted(warehouse.keys())
  ) + "\n"
  parser = EpilogParser(usage="libtbx.install [package]",
                        description="Installs an additional cctbx package.",
                        epilog=possible_installs)
  parser.add_option("-?", action="help", help=SUPPRESS_HELP)
  options, args = parser.parse_args(args)

  modules_directory = list(filter(lambda dir: not is_source_repository(dir), libtbx.env.repository_paths))
  if not modules_directory:
    sys.exit("No repository path candidate found. Can't install modules without an installation root.")
  if len(modules_directory) > 1:
    print("More than one repository path candidate found.")
  installation_root = modules_directory[0]
  print("Using %s as installation root" % abs(installation_root))

  packages_to_configure = set()

  errors = False
  for package in args:
    if package not in warehouse:
      print("Skipping package %s: Never heard of this before" % package)
      errors = True
      continue
    if (installation_root / package).isdir() and glob.glob(abs(installation_root / package / "*")):
      print("Skipping download of package %s: Non-empty directory already exists in installation root" % package)
      if package not in libtbx.env.module_dict and warehouse[package].get('configure', True):
        packages_to_configure.add(package)
      if warehouse[package].get('force-configure'):
        packages_to_configure.add('libtbx')
      continue

    downloaded = False
    try:
      os.makedirs(abs(installation_root / package))
    except OSError as exc:
      if exc.errno == errno.EEXIST and (installation_root / package).isdir():
        pass
      else:
        raise
    for mech, call in mechanisms.items():
      if mech in warehouse[package]:
        print("Attempting to obtain %s using %s..." % (package, mech), end='')
        sys.stdout.flush()
        if call(package=package,
                location=abs(installation_root / package),
                source=warehouse[package][mech]):
          assert (installation_root / package).isdir(), "Installation failed"
          downloaded = True
          print("success")
          break
        else:
          assert not (installation_root / package).isdir(), "Install mechanism %s did not fail cleanly" % mech
          print("failed")
    if not downloaded:
      print("Skipping package %s: Could not install" % package)
      errors = True
      continue

    if warehouse[package].get('configure', True):
      packages_to_configure.add(package)
    if warehouse[package].get('force-configure'):
      packages_to_configure.add('libtbx')

  if packages_to_configure:
    packages_to_configure = sorted(packages_to_configure)
    print("Configuring package[s] %s..." % ", ".join(packages_to_configure))

    proc = subprocess.Popen(['libtbx.configure'] + packages_to_configure,
      stdout=subprocess.PIPE, stderr=subprocess.PIPE,
      cwd=abs(libtbx.env.build_path)
      )
    out, err = proc.communicate()
    if err:
      errors = True
      print(out.decode())
      print(err.decode())
      print("Configuration failed. Run 'libtbx.configure %s' "
            "once underlying problem solved, then 'make'"
            % " ".join(packages_to_configure))
    else:
      proc = subprocess.Popen(['make'],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        cwd=abs(libtbx.env.build_path)
        )
      out, err = proc.communicate()
      if err:
        errors = True
  if errors:
    sys.exit(1)

def install_git(**kwargs):
  reference = []
  if os.name == 'posix':
    reference_repository_path = os.path.join('/dls/science/groups/scisoft/DIALS/repositories/git-reference', kwargs['package'])
    if os.path.isdir(reference_repository_path):
      reference = ['--reference', reference_repository_path]
      print("using reference repository...", end="")
      sys.stdout.flush()
  try:
    proc = subprocess.Popen(['git', 'clone', '--recursive', kwargs['source'], kwargs['location']] + reference,
      stdout=subprocess.PIPE, stderr=subprocess.PIPE,
      )
    out, err = proc.communicate()
    if err:
      if os.path.exists(kwargs['location']) and not os.listdir(kwargs['location']):
        # git-auth can leave an empty directory behind
        os.rmdir(kwargs['location'])
      return False
    if reference:
      proc = subprocess.Popen(['git', 'repack', '-a', '-d'],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        cwd=kwargs['location']
        )
      out, err = proc.communicate()
      print(out.decode())
      assert proc.returncode == 0, "Repack operation failed. Delete repository and try again."
      os.remove(os.path.join(kwargs['location'], '.git', 'objects', 'info', 'alternates'))
    Toolbox.set_git_repository_config_to_rebase(os.path.join(kwargs['location'], '.git', 'config'))
    return True
  except OSError:
    if os.path.isdir(kwargs['location']):
      shutil.rmtree(kwargs['location'])
    return False # git may not be installed

def install_tgz(**kwargs):
  location = kwargs['location']
  source = kwargs['source']
  tempfile = os.path.join(location, '.tmp.tgz')
  etagfile = os.path.join(location, '..tmp.tgz.etag')
  def cleanup():
    try: os.remove(tempfile)
    except OSError: pass
    try: os.remove(etagfile)
    except OSError: pass
    try: os.rmdir(location)
    except OSError: pass
  if Toolbox.download_to_file(source['url'], tempfile) <= 0:
    cleanup()
    return False
  import tarfile
  with tarfile.open(tempfile, 'r') as fh:
    for t in range(source.get('trim', 0)):
      location = os.path.dirname(location)
    fh.extractall(location)
  cleanup()
  return True

def install_zip(**kwargs):
  location = kwargs['location']
  source = kwargs['source']
  tempfile = os.path.join(location, '.tmp.zip')
  etagfile = os.path.join(location, '..tmp.zip.etag')
  def cleanup():
    try: os.remove(tempfile)
    except OSError: pass
    try: os.remove(etagfile)
    except OSError: pass
    try: os.rmdir(location)
    except OSError: pass
  if Toolbox.download_to_file(source['url'], tempfile) <= 0:
    cleanup()
    return False
  Toolbox.unzip(tempfile, location, trim_directory=source.get('trim', 0))
  cleanup()
  return True

def install_pip(**kwargs):
  git_installation = install_git(**kwargs)
  if not git_installation:
    return False
  proc = subprocess.Popen(['libtbx.pip', 'install', '-e', kwargs['location']],
    stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
  out, err = proc.communicate()
  print(out.decode())
  if err:
    return False
  return True

mechanisms = collections.OrderedDict((
  ('git-auth', install_git),
  ('git-anon', install_git),
  ('http-zip', install_zip),
  ('http-tgz', install_tgz),
  ('pip-auth', install_pip),
  ('pip-anon', install_pip),
))

warehouse = {
  'dials_scratch': {
    'git-auth': 'git@github.com:/dials/dials_scratch',
    'git-anon': 'https://github.com/dials/dials_scratch.git',
    'http-zip': { 'url': 'https://github.com/dials/dials_scratch/archive/master.zip', 'trim': 1 },
  },
  'dlstbx': {
    'git-auth': 'git@github.com:/DiamondLightSource/python-dlstbx.git',
  },
  'dxtbx': {
    'git-auth': 'git@github.com:cctbx/dxtbx',
    'git-anon': 'https://github.com/cctbx/dxtbx.git',
    'http-zip': { 'url': 'https://github.com/cctbx/dxtbx/archive/main.zip', 'trim': 1 },
  },
  'fast_dp': {
    'pip-auth': 'git@github.com:/DiamondLightSource/fast_dp',
    'pip-anon': 'https://github.com/DiamondLightSource/fast_dp.git',
    'configure': False,
    'force-configure': True,
  },
  "iota": {
    'pip-auth': 'git@github.com:/ssrl-px/iota',
    'pip-anon': 'https://github.com/ssrl-px/iota.git',
    'configure': False,
    'force-configure': True,
  },
  'screen19': {
    'pip-auth': 'git@github.com:/xia2/screen19',
    'pip-anon': 'https://github.com/xia2/screen19.git',
    'configure': False,
    'force-configure': True,
  },
  'msgpack': {
    'http-tgz': { 'url': 'https://gitcdn.link/repo/dials/dependencies/dials-1.13/msgpack-3.1.1.tar.gz', 'trim': 1 },
    'configure': False,
  },
}

if __name__ == '__main__':
  run(sys.argv[1:])
