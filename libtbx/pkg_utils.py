# Python package related functions
#
# To let libtbx modules specify and satisfy their requirements and capabilities.

from __future__ import absolute_import, division, print_function

import contextlib
import itertools
import json
import os
import sys

import libtbx.load_env

try:
  import conda.cli.python_api
except ImportError:
  conda = None

try:
  import pip
  import pkg_resources

  # Don't run if pip version < 9.0.0.
  # There is no technical reason for this, the code should work still.
  # But in the interest of not upsetting legacy build systems let's be cautious.
  if not all(symbol in dir(pkg_resources) for symbol in
             ('parse_version', 'require', 'DistributionNotFound', 'VersionConflict')) \
     or pkg_resources.parse_version(pip.__version__) < pkg_resources.parse_version('9.0.0'):
    pip = None
    pkg_resources = None
  if pip:
    if pkg_resources.parse_version(pip.__version__) >= pkg_resources.parse_version('19.3.0'):
      import pip._internal.main
      pip_main = pip._internal.main.main
    elif pkg_resources.parse_version(pip.__version__) >= pkg_resources.parse_version('10.0.0'):
      import pip._internal
      pip_main = pip._internal.main
    else:
      pip_main = pip.main
except ImportError:
  pip = None
  pkg_resources = None

# Try to find packaging. This is normally(?) installed by setuptools but
# otherwise pip keeps a copy.
try:
  import packaging
  from packaging.requirements import Requirement
except ImportError:
  try:
    import pip._vendor.packaging as packaging
    from pip._vendor.packaging.requirements import Requirement
  except ImportError:
    # If all else fails then we need the symbol to check
    Requirement = None
    packaging = None

try:
  import setuptools
except ImportError:
  setuptools = None

def _notice(*lines, **context):
  print(os.linesep + "=" * 80 + os.linesep + os.linesep +
        os.linesep.join(l.format(**context) for l in lines) + os.linesep +
        os.linesep + "=" * 80 + os.linesep)

_defined_entrypoints = set()

def require(pkgname, version=None):
  '''Ensure a package requirement is met. Install or update package as required
     and print a warning message if this can't be done due to the local
     environment, or when automatic package management is disabled by setting
     the environment variable 'LIBTBX_DISABLE_UPDATES'.
     :param pkgname: A string describing the package requirement. This will
                     generally just be a package name, but package features
                     can be specified in square brackets. Features are not
                     enforced, but will be requested during installation and
                     update.
     :param version: An optional string describing version constraints. This
                     can be a minimum version, eg. '>=1.0', a maximum version,
                     eg. '<2', or both, eg. '>=4.5,<4.6'.
     :return: True when the requirement is met, False otherwise.'''

  if not pip:
    _notice("  WARNING: Can not verify python package requirements - pip/setuptools out of date",
            "  Please update pip and setuptools by running:", "",
            "    libtbx.python -m pip install pip setuptools --upgrade", "",
            "  or following the instructions at https://pip.pypa.io/en/stable/installing/")
    return False

  if not version:
    version = ''

  requirement = Requirement(pkgname+version)

  # Check if we have an environment marker in the request
  if requirement.marker and not requirement.marker.evaluate():
    # We skip dependencies that don't match our current environment
    return True
  # Erase the marker from any further output
  requirement.marker = None

  # package name without feature specification
  basepkgname = requirement.name

  requirestring = str(requirement)
  baserequirestring = requirement.name + str(requirement.specifier)
  try:
    try:
      print("requires %s, has %s" % (requirestring, pkg_resources.require(requirestring)[0].version))
      return True
    except pkg_resources.UnknownExtra:
      print("requires %s, has %s, but without features" % (requirestring, pkg_resources.require(baserequirestring)[0].version))
      return True

  except pkg_resources.DistributionNotFound:
    currentversion = '(not determined)'
    project_name = pkgname
    action = 'install'
    print("requirement %s is not currently met, package not installed" % (requirestring))

  except pkg_resources.VersionConflict:
    currentversion = pkg_resources.require(basepkgname)[0].version
    project_name = pkg_resources.require(basepkgname)[0].project_name
    action = 'update'
    print("requirement %s is not currently met, current version %s" % (requirestring, currentversion))

  # Check if package can be updated
  for path_item in sys.path:
    egg_link = os.path.join(path_item, project_name + '.egg-link')
    if os.path.isfile(egg_link):
      with open(egg_link, 'r') as fh:
        location = fh.readline().strip()
      _notice("    WARNING: Can not update package {package} automatically.", "",
              "It is installed as editable package for development purposes. The currently",
              "installed version, {currentversion}, is too old. The required version is {requirement}.",
              "Please update the package manually in its installed location:", "",
              "    {location}",
              package=pkgname, currentversion=currentversion, requirement=version, location=location)
      return False

  if not os.path.isdir(libtbx.env.under_base('.')):
    _notice("    WARNING: Can not {action} package {package} automatically.", "",
            "You are running in a base-less installation, which disables automatic package",
            "management by convention, see https://github.com/cctbx/cctbx_project/issues/151", "",
            "Please {action} the package manually.",
            package=pkgname, currentversion=currentversion, requirement=version, action=action)
    return False

  if os.getenv('LIBTBX_DISABLE_UPDATES') and os.getenv('LIBTBX_DISABLE_UPDATES').strip() not in ('0', ''):
    _notice("    WARNING: Can not {action} package {package} automatically.", "",
            "Environment variable LIBTBX_DISABLE_UPDATES is set. Please {action} manually.",
            package=pkgname, currentversion=currentversion, requirement=version, action=action)
    return False

  if libtbx.env.build_options.use_conda:
    _notice("    WARNING: Can not {action} package {package} automatically.", "",
            "You are in a conda environment. Please {action} manually.",
            package=pkgname, currentversion=currentversion, requirement=version, action=action)
    return False

  print("attempting {action} of {package}...".format(action=action, package=pkgname))
  has_req_tracker = os.environ.get('PIP_REQ_TRACKER')
  exit_code = pip_main(['install', requirestring])
  if not has_req_tracker:
    # clean up environment after pip call for next invocation
    os.environ.pop('PIP_REQ_TRACKER', None)
  if exit_code == 0:
    print("{action} successful".format(action=action))
    return True
  else:
    print("{action} failed. please check manually".format(action=action))
    return False

@contextlib.contextmanager
def _silence():
  '''Helper context which shuts up stdout.'''
  if os.name == "nt":
    # Can't silence using this method on Windows. Just leave it.
    yield
    return
  sys.stdout.flush()
  try:
    oldstdout = os.dup(sys.stdout.fileno())
    dest_file = open(os.devnull, 'w')
    os.dup2(dest_file.fileno(), sys.stdout.fileno())
    yield
  finally:
    if oldstdout is not None:
      os.dup2(oldstdout, sys.stdout.fileno())
    if dest_file is not None:
      dest_file.close()

def define_entry_points(epdict, **kwargs):
  '''A function to allow non-setuptools packages (ie. libtbx modules) to use
     the setuptools entry_points mechanism. Call this function from
     libtbx_refresh.py and pass a dictionary of entry points.'''
  # Determine the name of the calling module, and thus the internal module name
  # of the run_tests file. Use exception trick to pick up the current frame.
  try:
    raise Exception()
  except Exception:
    frame = sys.exc_info()[2].tb_frame.f_back
  # Extract the caller name
  caller = frame.f_globals['__name__']
  if caller == '__main__':
    # well that is not very informative, is it.
    caller = os.path.abspath(frame.f_code.co_filename)  # Get the full path of the libtbx_refresh.py file.
    refresh_file, _ = os.path.splitext(caller)
    if not refresh_file.endswith('libtbx_refresh'):
      raise RuntimeError('Entry points can only be defined from within libtbx_refresh.py')
    # the name of the parent directory of libtbx_refresh.py is the caller name
    caller = os.path.basename(os.path.dirname(refresh_file))
  else:
    if not caller.endswith('.libtbx_refresh'):
      raise RuntimeError('Entry points can only be defined from within libtbx_refresh.py')
    caller = caller[:-15]

  # No setuptools mechanism without setuptools.
  if not setuptools:
    raise ImportError("You must install setuptools to configure package {}. Run\n  libtbx.pip install setuptools".format(caller))

  if caller in _defined_entrypoints:
    raise RuntimeError("Entry points have already been defined for package %s. There must only be a single call to "
                       "define_entry_points() per calling package" % caller)
  _defined_entrypoints.add(caller)

  print("Updating entry points for {caller}".format(caller=caller))
  for ep in epdict:
    print("  {n} entries for entry point {ep}".format(ep=ep, n=len(epdict[ep])))

  # Temporarily change to {libtbx.env.build_path}/lib directory. This
  # is where a directory named libtbx.{caller}.egg-info will be
  # created containing the entry point info.
  try:
    curdir = os.getcwd()
    os.chdir(os.path.join(abs(libtbx.env.build_path), 'lib'))
    # Now trick setuptools into thinking it is in control here.
    try:
      argv_orig = sys.argv
      sys.argv = ['setup.py', 'egg_info']
      # And make it run quietly
      with _silence():
        setuptools.setup(
          name='libtbx.{}'.format(caller),
          description='libtbx entry point manager for {}'.format(caller),
          entry_points=epdict,
          **kwargs
        )
    finally:
      sys.argv = argv_orig
  finally:
    os.chdir(curdir)

  # With the entry points installed into build/lib check for any legacy entry
  # point definitions in both build/ and base/ and remove them where possible.
  # Code block can be removed April 2020
  try:
    legacy_path = abs(libtbx.env.build_path / 'libtbx.{}.egg-info'.format(caller))
    if os.path.exists(legacy_path):
      print("Removing legacy entry points from", legacy_path)
      import shutil
      shutil.rmtree(legacy_path)
  except Exception as e:
    print("Could not remove legacy entry points:", str(e))
  try:
    import site
    for spp in site.getsitepackages():
      legacy_path = os.path.join(spp, 'libtbx.{}.egg-link'.format(caller))
      if os.path.exists(legacy_path):
        print("Removing legacy entry points from", legacy_path)
        os.remove(legacy_path)
  except Exception as e:
    print("Could not remove legacy entry points:", str(e))


def _merge_requirements(requirements, new_req):
  # type: (List[packaging.requirements.Requirement], packaging.requirements.Requirement) -> None
  """Merge a new requirement with a list.

  If it exists in an identical form (name, marker) then the
  specifiers and extras will be merged, otherwise it will be added.

  If the environment markers are different it will assume that they
  are mutually exclusive - entries will only be merged if identical, which
  could cause problems with duplicate requirement entries if not filtered
  by pass status later.

  URL field is also not handled, as unsure how to merge these if they differ.

  Args:
    requirements (List[packaging.requirements.Requirement]): Existing.
    new_req (packaging.requirements.Requirement): New requirement
  """
  assert new_req.url is None, "URL requirement fields not handled/tested"
  matches = [
    x
    for x in requirements
    if x.name == new_req.name
    and x.marker == new_req.marker
  ]
  if len(matches) > 0:
    if len(matches) > 1:
      print("Warning: More than one match for requirement", new_req, ": ", matches)
    match = matches[0]
  else:
    match = None

  if match:
    match.specifier = match.specifier & new_req.specifier
    match.modules = match.modules | new_req.modules
    match.extras = match.extras | new_req.extras
  else:
    requirements.append(new_req)

def collate_conda_requirements(modules):
  # type: (List[libtbx.env_config.module]) -> List[pkg_resources.Requirement]
  """Combine conda requirements from a module list.

  An attempt will be made to merge any joint requirements. The requirement
  objects will have an added property 'modules', which is a set of module
  names that formed the requirement.

  Attr:
      modules (Iterable[libtbx.env_config.module]): The module list

  Returns:
      List[pkg_resources.Requirement]: The merged requirements
  """
  requirements = []
  for module, spec in itertools.chain(*[[(x.name, y) for y in x.conda_required] for x in modules if hasattr(x, "conda_required")]):
    requirement = pkg_resources.Requirement.parse(spec)
    # Track where dependencies came from
    requirement.modules = {module}
    # Attempt to merge this with any other requirements to avoid double-specifying
    _merge_requirements(requirements, requirement)
  return requirements

def collate_python_requirements(modules):
  # type: (List[libtbx.env_config.module]) -> List[packaging.requirements.Requirement]
  """Combine python requirements from a module list.

  An attempt will be made to merge any joint requirements. The requirement
  objects will have an added property 'modules', which is a set of module
  names that formed the requirement.

  Attr:
      modules (Iterable[libtbx.env_config.module]): The module list

  Returns:
      List[packaging.requirements.Requirement]: The merged requirements
  """
  requirements = []
  for module, spec in itertools.chain(*[[(x.name, y) for y in x.python_required] for x in modules if hasattr(x, "python_required")]):
    requirement = Requirement(spec)
    # Track where dependencies came from
    requirement.modules = {module}
    # Attempt to merge this with any other requirements to avoid double-specifying
    _merge_requirements(requirements, requirement)
  return requirements

def resolve_module_conda_dependencies(modules):
  # type: (List[libtbx.env_config.module]) -> None
  """Resolve all python dependencies from the list of modules"""
  # Check we can do anything here
  if not libtbx.env.build_options.use_conda:
    return

  if conda is None:
    _notice("  WARNING: Can not find conda package in your environment",
            "  You will have to keep track of dependencies yourself")
    return

  conda_list, error, return_code = conda.cli.python_api.run_command(
    conda.cli.python_api.Commands.LIST,
    "--json",
    use_exception_handler=True,
  )
  if error or return_code:
    _notice("  WARNING: Could not obtain list of conda packages in your environment",
            error)
    return
  conda_environment = {package["name"]: package["version"] for package in json.loads(conda_list)}

  requirements = collate_conda_requirements(modules)
  # Now we should have an unduplicated set of requirements
  action_list = []
  for requirement in requirements:
    # Check if package is installed in development mode
    if pkg_resources:
      try:
        currentversion = pkg_resources.require(requirement.name)[0].version
      except Exception:
        pass
      else:
        location = None
        for path_item in sys.path:
          egg_link = os.path.join(path_item, requirement.name + '.egg-link')
          if os.path.isfile(egg_link):
            with open(egg_link, 'r') as fh:
              location = fh.readline().strip()
              break
        if location and currentversion in requirement:
          print("requires conda package %s, has %s as developer installation" % (requirement, currentversion))
          continue
        elif location and currentversion not in requirement:
          _notice(
            "    WARNING: Can not update package {package} automatically.", "",
            "It is installed as editable package for development purposes. The currently",
            "installed version, {currentversion}, is too old. The required version is {requirement}.",
            "Please update the package manually in its installed location:", "",
            "    {location}",
            package=requirement.name, currentversion=currentversion,
            requirement=requirement, location=location)
          continue

    # Check if package is installed with conda
    if requirement.name in conda_environment:
      if conda_environment[requirement.name] in requirement:
        print("requires conda package %s, has %s" % (requirement, conda_environment[requirement.name]))
        continue
      print("conda requirement %s is not currently met, current version %s" % (requirement, conda_environment[requirement.name]))

    # Install/update required
    print("conda requirement %s is not currently met, package not installed" % (requirement))
    action_list.append(str(requirement))

  if not action_list:
    print("All conda requirements satisfied")
    return

  if not os.path.isdir(libtbx.env.under_base('.')):
    _notice("    WARNING: Can not update conda packages automatically.", "",
            "You are running in a base-less installation, which disables automatic package",
            "management by convention, see https://github.com/cctbx/cctbx_project/issues/151", "",
            "Please update the following packages manually:",
            "  {action_list}",
            action_list=", ".join(action_list))
    return

  if os.getenv('LIBTBX_DISABLE_UPDATES') and os.getenv('LIBTBX_DISABLE_UPDATES').strip() not in ('0', ''):
    _notice("    WARNING: Can not automatically update conda environment", "",
            "Environment variable LIBTBX_DISABLE_UPDATES is set.",
            "Please update the following packages manually:",
            "  {action_list}",
            action_list=", ".join(action_list))
    return

  print("\nUpdating conda environment for packages:" + "".join("\n - " + a for a in action_list) + "\n")
  _, _, return_code = conda.cli.python_api.run_command(
    conda.cli.python_api.Commands.INSTALL,
    *action_list,
    stdout=None,
    stderr=None,
    use_exception_handler=True
  )
  if return_code:
    _notice("    WARNING: Could not automatically update conda environment", "",
            "Please check your environment manually.")

def resolve_module_python_dependencies(modules):
  # type: (List[libtbx.env_config.module]) -> None
  """Resolve all python dependencies from the list of modules"""
  # Check we can do anything here
  if Requirement is None:
    _notice("  WARNING: Can not find package requirements tools - pip/setuptools out of date?",
            "  Please update pip and setuptools by running:", "",
            "    libtbx.python -m pip install pip setuptools --upgrade", "",
            "  or following the instructions at https://pip.pypa.io/en/stable/installing/")
    return
  requirements = collate_python_requirements(modules)
  # Now we should have an unduplicated set of requirements
  for requirement in requirements:
    # Pass everything as the requirement string rather than trying to reconstruct
    result = require(str(requirement), "")
