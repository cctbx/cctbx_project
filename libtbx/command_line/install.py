from __future__ import absolute_import, division, print_function

import collections
import os
import shutil
import sys
from optparse import SUPPRESS_HELP, OptionParser

import procrunner

import libtbx.load_env
from libtbx.auto_build.bootstrap import Toolbox

# Basically 'pip' for selected libtbx/cctbx modules.


def is_source_repository(path):
  return (path / '.git').isdir() or (path / '.svn').isdir()

def run(args):
  parser = OptionParser(usage="libtbx.install [package]",
                        description="Installs an additional cctbx package")
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
    if (installation_root / package).isdir():
      print("Skipping download of package %s: Directory already exists in installation root" % package)
      if package not in libtbx.env.module_dict:
        packages_to_configure.add(package)
      continue
    if package not in warehouse:
      print("Skipping package %s: Never heard of this before" % package)
      errors = True
      continue

    downloaded = False
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

    packages_to_configure.add(package)

  if packages_to_configure:
    packages_to_configure = sorted(packages_to_configure)
    print("Configuring package[s] %s" % ", ".join(packages_to_configure))
    os.chdir(abs(libtbx.env.build_path))
    result = procrunner.run_process(['libtbx.configure'] + packages_to_configure,
                                    print_stdout=False, print_stderr=False)
    if result['exitcode']:
      errors = True
      print(result['stdout'])
      print(result['stderr'])
      print("Configuration failed. Run 'libtbx.configure %s' "
            "once underlying problem solved, then 'make'"
            % " ".join(packages_to_configure))
    else:
      result = procrunner.run_process(['make'])
      if result['exitcode']:
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
    result = procrunner.run_process(['git', 'clone', '--recursive', kwargs['source'], kwargs['location']] + reference,
                                    print_stdout=False)
    if result['exitcode']:
      return False
    if reference:
      oldcwd = os.getcwd()
      os.chdir(kwargs['location'])
      result = procrunner.run_process(['git', 'repack', '-a', '-d'], print_stderr=True)
      os.chdir(oldcwd)
      assert result['exitcode'] == 0, "Repack operation failed. Delete repository and try again."
      os.remove(os.path.join(kwargs['location'], '.git', 'objects', 'info', 'alternates'))
    Toolbox.set_git_repository_config_to_rebase(os.path.join(kwargs['location'], '.git', 'config'))
    return True
  except OSError:
    if os.path.isdir(kwargs['location']):
      shutil.rmtree(kwargs['location'])
    return False # git may not be installed

def install_zip(**kwargs):
  location = kwargs['location']
  source = kwargs['source']
  os.mkdir(location)
  tempfile = os.path.join(location, '.tmp.zip')
  etagfile = os.path.join(location, '..tmp.zip.etag')
  def cleanup():
    try: os.remove(tempfile)
    except OSError: pass
    try: os.remove(etagfile)
    except OSError: pass
  if Toolbox.download_to_file(source['url'], tempfile) <= 0:
    cleanup()
    os.rmdir(location)
    return False
  Toolbox.unzip(tempfile, location, trim_directory=source.get('trim', 0))
  cleanup()
  return True

mechanisms = collections.OrderedDict((
  ('git-auth', install_git),
  ('git-anon', install_git),
  ('http-zip', install_zip),
))

warehouse = {
  'dials_scratch': {
    'git-auth': 'git@github.com:/dials/dials_scratch',
    'git-anon': 'https://github.com/dials/dials_scratch.git',
    'http-zip': { 'url': 'https://github.com/dials/dials_scratch/archive/master.zip', 'trim': 1 },
  },
  'dlstbx': {
    'git-auth': 'dascgitolite@dasc-git.diamond.ac.uk:/dials/dlstbx.git',
  },
  'i19': {
    'git-auth': 'git@github.com:/xia2/i19',
    'git-anon': 'https://github.com/xia2/i19.git',
    'http-zip': { 'url': 'https://github.com/xia2/i19/archive/master.zip', 'trim': 1 },
  },
  'xia2_regression': {
    'git-auth': 'git@github.com:/xia2/xia2_regression',
    'git-anon': 'https://github.com/xia2/xia2_regression.git',
    'http-zip': { 'url': 'https://github.com/xia2/xia2_regression/archive/master.zip', 'trim': 1 },
  },
}

if __name__ == '__main__':
  run(sys.argv[1:])
