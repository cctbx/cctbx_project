"""
Manage conda environment for CCTBX programs

This is a minimal wrapper for conda functionality. It is meant to create
the minimum environment for building CCTBX programs and is aimed towards
developers that do not want to manage their own conda environments.
Advanced users should use conda natively so that they have full control.

This is called by bootstrap.py to set up a developer's conda
environment. If no conda installation is found, one will be downloaded.
It will be installed in the "mc3" directory at the same level as
"modules."

The minimum version for conda is 4.4. This is when conda moved towards
building packages with a set of common compilers.
"""
from __future__ import absolute_import, division, print_function

import argparse
import json
import os
import platform
import shutil
import sys
import time
import warnings

# Python 2/3 compatibility
py2 = False
py3 = True
try:
  from urllib.parse import urljoin
  from urllib.request import urlopen
except ImportError:
  from urlparse import urljoin
  from urllib2 import urlopen
  py2 = True
  py3 = False

# copied from install_base_packages.py
if __package__ is None:
  sys.path.append(os.path.dirname(os.path.abspath(__file__)))
  from installer_utils import check_output, call
else:
  from .installer_utils import check_output, call

# conda on Windows seems to need cmd and the wait() in Popen
if platform.system() == 'Windows':
  import subprocess
  import tempfile
  def check_output(command_list, *args, **kwargs):
    # check for "conda info" and "activate" commands and prepend cmd
    if 'conda.exe' in command_list[0] or 'activate' in command_list[0]:
      command_list = ['cmd', '/c'] + command_list
      output = ''
      with tempfile.TemporaryFile() as f:
        returncode = subprocess.check_call(command_list, stdout=f, *args, **kwargs)
        f.seek(0)
        output = f.read()
      return output
    # miniconda3 installation
    else:
      returncode = call(command_list, *args, **kwargs)
      return returncode

# =============================================================================
# Locations for the files defining the conda environments
# Generally, they should reside in the repository of the main program
# defined by the builder. For example, the environment file for Phenix
# will be in the Phenix source tree.
conda_platform = {
  'Darwin': 'osx-64',
  'Linux': 'linux-64',
  'Windows': 'win-64',
}

# check for Apple Silicon
if platform.system() == 'Darwin' and 'arm64' in platform.platform():
  conda_platform['Darwin'] = 'osx-arm64'

version = 'PYTHON_VERSION'
default_format = '{builder}_py{version}_{platform}.txt'
default_filename = default_format.format(builder='cctbx',
  version=version, platform=conda_platform[platform.system()])

root_dir = os.path.abspath(
  os.path.join(__file__, '..', '..', '..', '..', '..'))
default_file = os.path.join('cctbx_project', 'libtbx', 'auto_build',
                            'conda_envs', default_filename)

# =============================================================================
def download_file(url, filename):
  """
  Simple function for downloading a file from a URL
  No error checking is done since anything downloaded is necessary.

  Parameters
  ----------
  url: str
    The url pointing to the location of the file
  filename: str
    The name of the local file

  Returns
  -------
    Nothing
  """
  if py3:
    with urlopen(url) as response:
      with open(filename, 'wb') as local_file:
        shutil.copyfileobj(response, local_file)
  elif py2:
    with open(filename, 'wb') as local_file:
      response = urlopen(url).read()
      local_file.write(response)

# =============================================================================
class conda_manager(object):
  """
  Class for managing conda environments

  Attributes
  ----------
  env_locations: dict
    Dictionary for determining which environment to use for each builder.
    The keys are the builder names and must match those in bootstrap.py
    The values are the relative path from "modules" to the file. Builders
    for programs outside of cctbx_project should place their environment
    files inside their repository. For example, the Phenix environment
    file should be in the phenix repository.

  Methods
  -------
  __init__(root_dir, conda_base, conda_env, verbose, log)
    Constructor that does the basic check for a conda installation
  get_conda_exe(prefix, check_file)
    Returns the platform-dependent conda executable
  get_conda_python(conda_env, check_file)
    Returns the platform-dependent conda python
  update_environments()
    Returns a list of paths that are conda environments based on the
    environments.txt file in ${HOME}/.conda
  install_miniconda(prefix)
    Downloads and installs the latest version of miniconda3
  update_conda()
    Updates the current version of conda
  create_environment(builder, filename, copy)
    Uses the known conda installtion to create an environment
  write_conda_setpaths(prefix, build_dir, conda_env, check_file)
    Writes an additional script that activates the conda environment before
    calling the normal setpaths script.
  """

  # Currently, there is one monolithic environment
  # The key names for this dictionary must match the builder names in
  # bootstrap.py
  phenix_env = os.path.join('phenix', 'conda_envs',
    default_format.format(builder='phenix', version=version,
                          platform=conda_platform[platform.system()]))
  env_locations = {
    'cctbxlite': default_file,
    'cctbx': default_file,
    'phenix': phenix_env,
    'phenix_discamb': phenix_env + '2',
    'phenix_molstar': phenix_env,
    'phenix_voyager': phenix_env,
    'phenix_release': phenix_env,
    'xfellegacy': default_file,
    'dials-old': os.path.join('dials', '.conda-envs',
      default_format.format(builder='dials', version=version,
                            platform=conda_platform[platform.system()])),
    'dials': os.path.join('dials', '.conda-envs',
             {
                 "Darwin": "macos.txt",
                 "Linux": "linux.txt",
                 "Windows": "windows.txt",
             }[platform.system()]),
    'external': default_file,
    'molprobity': default_file,
    'qrefine': default_file,
    'phaser': default_file,
    'voyager': os.path.join('phasertng', 'conda_envs',
      default_format.format(builder='phasertng', version=version,
                            platform=conda_platform[platform.system()]))
  }
  env_locations['xfel'] = env_locations['labelit'] = env_locations['dials']
  # A set of builders where the environment files do not specify the python
  # version
  env_without_python = [ "dials","xfel","labelit"]

  # ---------------------------------------------------------------------------
  def __init__(self, root_dir=root_dir, conda_base=None, conda_env=None,
               check_file=True, max_retries=5, verbose=False, log=sys.stdout):
    """
    Constructor that performs a basic check for the conda installation.
    If an installation is not found, the latest version can be downloaded
    and installed.

    Parameters
    ----------
    root_dir: str
      Required argument for specifying the root level directory for
      CCTBX builds. This is where the "modules" and "build" directories
      reside. If a new conda installation is required, it will be placed
      here under the "mc3" directory. The conda environment for building
      will also be placed in this directory
    conda_base: str
      Required argument for specifying location of a conda installation.
      If this is None, miniconda will be installed unless conda_env is
      set to a valid path.
    conda_env: str
      Optional argument for specifying location of a conda environment.
      Since this is assumed to be a working conda environment, a new
      installation of conda will not be created. This is for developers
      that want to manage their own environments.
    check_file: bool
      Flag for checking if a file exists. A RuntimeError is raised if a
      this flag is set and a file does not exist. Used in get_conda_exe
      and get_conda_python.
    max_retries: int
      When downloading conda packages, there may be network issues that
      prevent the environment from being constructed. This parameter
      controls the number of retry attempts for constructing the conda
      environment. The first retry is attempted 1 minute after the
      initial failure. The second retry is attempted 2 minutes, etc.
    verbose: bool
      Flag for showing conda output
    log: file
      For storing log output
    """
    self.system = platform.system()

    self.root_dir = root_dir

    self.conda_base = None
    if conda_base is not None:
      self.conda_base = os.path.normpath(conda_base)

    self.conda_env = None
    if conda_env is not None:
      self.conda_env = os.path.normpath(conda_env)

    self.check_file = check_file
    self.max_retries = max_retries
    self.verbose = verbose
    self.log = log

    self.conda_exe = None
    if self.conda_base is not None:
      self.conda_exe = self.get_conda_exe(self.conda_base, check_file=False)

    # default environment file for users
    self.environment_file = os.path.join(
      os.path.expanduser('~'), '.conda', 'environments.txt')
    self.environments = self.update_environments()

    # Clean environment for external Python processes
    self.env = os.environ.copy()
    self.env['PYTHONPATH'] = ''

    # error messages
    self.conda_exe_not_found = """
The conda executable cannot be found. Please make sure the correct
directory for the base conda installation was provided. The directory
can be found by running "conda info" and looking the "base environment"
value."""

    # try to determine base environment from $PATH
    if self.conda_base is None:
      paths = os.environ.get('PATH')
      if paths is not None:
        if self.system == 'Windows':
          paths = paths.split(';')
        else:
          paths = paths.split(':')
        for path in paths:
          conda_base = os.path.abspath(os.path.join(path, '..'))
          conda_exe = self.get_conda_exe(conda_base, check_file=False)
          if os.path.isfile(conda_exe):
            self.conda_base = conda_base
            self.conda_exe = conda_exe
            break

    # try to determine base environment from .conda/environments.txt
    if self.conda_base is None:
      for environment in self.environments:
        conda_exe = self.get_conda_exe(environment, check_file=False)
        if os.path.isfile(conda_exe):
          self.conda_base = environment
          self.conda_exe = conda_exe
          break

    # -------------------------------------------------------------------------
    # helper function for searching for conda_exe
    def walk_up_check(new_conda_base, new_conda_exe):
      if new_conda_base is None:
        new_conda_base = ''
      if new_conda_exe is None:
        new_conda_exe = ''
      while not os.path.isfile(new_conda_exe):
        # move up directory and recheck
        new_conda_base = os.path.abspath(os.path.join(new_conda_base, '..'))
        new_conda_exe = self.get_conda_exe(new_conda_base, check_file=False)
        if os.path.isfile(new_conda_exe):
          return new_conda_base, new_conda_exe

        # moved to root directory
        if new_conda_base == \
          os.path.abspath(os.path.join(new_conda_base, '..')):
          return '', ''
    # -------------------------------------------------------------------------

    # try to determine base evironment from conda_env
    if (self.conda_base is None) and (self.conda_env is not None):
      conda_base, conda_exe = walk_up_check(self.conda_env, self.conda_exe)
      if os.path.isfile(conda_exe):
        self.conda_base = conda_base
        self.conda_exe = conda_exe

    # install conda if necessary
    if (self.conda_base is None) and (self.conda_env is None):
      install_dir = os.path.join(self.root_dir, 'mc3')
      if os.path.isdir(install_dir):
        print('Using default conda installation', file=self.log)
        self.conda_base = install_dir
      else:
        print('Location of conda installation not provided', file=self.log)
        print('Proceeding with a fresh installation', file=self.log)
        self.conda_base = self.install_miniconda(prefix=self.root_dir)
      self.conda_exe = self.get_conda_exe(self.conda_base, check_file=False)

    self.environments = self.update_environments()

    # verify consistency and check conda version
    if self.conda_base is not None:

      # maybe a conda evironment was provided instead of the base environment
      if not os.path.isfile(self.conda_exe):
        self.conda_base, self.conda_exe = \
          walk_up_check(self.conda_base, self.conda_exe)
        self.environments = self.update_environments()

      if not os.path.isfile(self.conda_exe):
        raise RuntimeError(self.conda_exe_not_found)

      conda_info = json.loads(check_output([self.conda_exe, 'info', '--json'],
                                           env=self.env))
      consistency_check = [self.conda_base == conda_info['root_prefix']]
      for env in self.environments:
        consistency_check.append(env in conda_info['envs'])
      if False in consistency_check:
        message = """
There is a mismatch between the conda settings in your home directory
and what "conda info" is reporting. This is not a fatal error, but if
an error is encountered, please check that your conda installation and
environments exist and are working.
"""
        warnings.warn(message, RuntimeWarning)
      split_version = conda_info['conda_version'].split('.')
      major_version = int(split_version[0])
      minor_version = int(split_version[1])
      if major_version < 4 or (major_version == 4 and minor_version < 4):
        raise RuntimeError("""
CCTBX programs require conda version 4.4 and greater to make use of the
common compilers provided by conda. Please update your version with
"conda update conda".
""")

      print('Base conda installation:\n  {base}'.format(base=self.conda_base),
            file=self.log)
      if self.verbose:
        output = check_output([self.conda_exe, 'info'], env=self.env)
        print(output, file=self.log)

    if self.conda_env is not None:
      print('Build environment:\n  {conda_env}'.format(
        conda_env=self.conda_env), file=self.log)

  # ---------------------------------------------------------------------------
  def get_conda_exe(self, prefix=None, check_file=None):
    """
    Find the conda executable. This is platform-dependent

    Parameters
    ----------
    prefix: str
      The path to the base conda environment
    check_file: bool
      Used to override the check_file attribute

    Returns
    -------
    conda_exe: str
      The path to the conda executable
    """
    if prefix is None:
      prefix = self.conda_base

    if check_file is None:
      check_file = self.check_file

    if self.system == 'Windows':
      conda_exe = os.path.join(prefix, 'Scripts', 'conda.exe')
    else:
      conda_exe = os.path.join(prefix, 'bin', 'conda')

    if check_file:
      if not os.path.isfile(conda_exe):
        raise RuntimeError(self.conda_exe_not_found)

    return conda_exe

  # ---------------------------------------------------------------------------
  def get_conda_python(self, conda_env=None, check_file=None):
    """
    Find the conda python. This is platform-dependent

    Parameters
    ----------
    conda_env: str
      The path to the conda environment
    check_file: bool
      Used to override the check_file attribute

    Returns
    -------
    conda_python: str
      The path to the conda python
    """
    if conda_env is None:
      conda_env = self.conda_env
    if conda_env is None:
      conda_env = self.conda_base
    if check_file is None:
      check_file = self.check_file

    if self.system == 'Windows':
      conda_python = os.path.join(conda_env, 'python.exe')
    elif self.system == 'Darwin':
      # use python.app for GUI applications
      conda_python = os.path.join(conda_env, 'python.app', 'Contents', 'MacOS',
                                  'python')
    else:
      conda_python = os.path.join(conda_env, 'bin', 'python')

    if check_file:
      if not os.path.isfile(conda_python):
        raise RuntimeError("""
  Python could not be located in {conda_env}. Please make sure {conda_env}
  is a valid conda environment and that Python is installed.
  """.format(conda_env=conda_env))

    return conda_python

  # ---------------------------------------------------------------------------
  def update_environments(self):
    """
    Read and check for existence of environment directories

    Returns
    -------
    environments: list
      List of paths that exist on the filesystem
    """
    environments = set()
    try:
      with open(self.environment_file) as f:
        paths = f.readlines()
      for env in paths:
        env = env.strip()
        if os.path.isdir(env):
          environments.add(os.path.normpath(env))
    except IOError:
      pass

    if self.conda_base is not None:
      env_dirs = [os.path.join(self.conda_base, 'envs'),
                  os.path.join(os.path.expanduser('~'), '.conda', 'envs')]
      for env_dir in env_dirs:
        if os.path.isdir(env_dir):
          dirs = os.listdir(env_dir)
          for _dir in dirs:
            _dir = os.path.join(env_dir, _dir)
            if os.path.isdir(_dir):
              environments.add(_dir)

    return environments

  # ---------------------------------------------------------------------------
  def install_miniconda(self, prefix=None):
    """
    Download a miniconda installer and install. The default installation
    location is at the same directory level as the "modules" directory.

    Parameters
    ----------
    prefix: str
      The installation directory for miniconda

    Returns
    -------
    conda_base: str
      The location of the "base" conda environment
    """

    if prefix is None:
      prefix = self.root_dir

    # construct Miniconda3 filename
    os_names = {
      'Darwin': 'MacOSX',
      'Linux': 'Linux',
      'Windows': 'Windows',
    }
    filename = 'Miniconda3-latest-{platform}-x86_64'.format(
      platform=os_names[self.system])
    if conda_platform[self.system] == 'osx-arm64':
      filename = 'Miniconda3-latest-{platform}-arm64'.format(
        platform=os_names[self.system])
    if self.system == 'Windows':
      filename += '.exe'
    else:
      filename += '.sh'
    url_base = 'https://repo.anaconda.com/miniconda/'
    url = urljoin(url_base, filename)
    filename = os.path.join(prefix, filename)

    # Download from public repository
    if not os.path.isfile(filename):
      print('Downloading {url}'.format(url=url), file=self.log)
      download_file(url, filename)
      print('Downloaded file to {filename}'.format(filename=filename),
        file=self.log)
    else:
      print('Using local copy at {filename}'.format(filename=filename),
        file=self.log)

    # run the installer
    install_dir = os.path.join(prefix, 'mc3')
    if self.system == 'Windows':
      flags = '/InstallationType=JustMe /RegisterPython=0 /AddToPath=0 /S /D={install_dir}'.\
        format(install_dir=install_dir)
      command_list = ['"' + filename + '"', flags]
    else:
      flags = ['-b', '-p', '{install_dir}'.format(install_dir=install_dir)]
      command_list = ['/bin/sh', filename] + flags
    print('Installing miniconda to "{install_dir}"'.format(
      install_dir=install_dir), file=self.log)
    output = check_output(command_list, env=self.env)
    if self.verbose:
      print(output, file=self.log)

    return install_dir

  # ---------------------------------------------------------------------------
  def update_conda(self):
    """
    Update the version of conda, if possible. The defaults channel is
    used because that is the default for a normal miniconda installation.

    Parameters
    ----------
      None
    """
    command_list = [self.conda_exe, 'update', '-n', 'base', '-c', 'defaults',
                    '-y', 'conda']
    try:
      output = check_output(command_list, env=self.env)
    except Exception:
      print("""
*******************************************************************************
There was a failure in updating your base conda installaion. To update
manually, try running

  conda update -n base -c defaults conda

If you are using conda from a different channel, replace "defaults" with that
channel
*******************************************************************************
""")
    else:
      print(output)

  # ---------------------------------------------------------------------------
  def _retry_command(self, command_list, text, prefix, verbose=None):
    """
    Internal convenience function for retrying creation/update of a
    conda environment.

    Parameters
    ----------
      command_list: list
        Command list for check_output
      text: str
        Text to be printed. Usually "installion into" or "update of"
      prefix: str
        Directory of conda environment
      verbose: bool
        Used to override self.verbose

    Returns
    -------
      Nothing
    """

    run_command = check_output
    if self.verbose or verbose:
      run_command = call

    for retry in range(self.max_retries):
      retry += 1
      try:
        output = run_command(command_list, env=self.env)
      except Exception:
        print("""
*******************************************************************************
There was a failure in constructing the conda environment.
Attempt {retry} of {max_retries} will start {retry} minute(s) from {t}.
*******************************************************************************
""".format(retry=retry, max_retries=self.max_retries, t=time.asctime()))
        time.sleep(retry*60)
      else:
        break
    if retry == self.max_retries:
      raise RuntimeError("""
The conda environment could not be constructed. Please check that there is a
working network connection for downloading conda packages.
""")
    print('Completed {text}:\n  {prefix}'.format(text=text, prefix=prefix),
          file=self.log)

    # check that environment file is updated
    self.environments = self.update_environments()
    if prefix not in self.environments:
      raise RuntimeError("""
The newly installed environment cannot be found in
${HOME}/.conda/environments.txt.
""")

  # ---------------------------------------------------------------------------
  def create_environment(self, builder='cctbx', filename=None, python=None,
    copy=False, offline=False):
    """
    Create the environment based on the builder and file. The
    environment name is "conda_base".

    Parameters
    ----------
    builder: str
      The builder from bootstrap.py. The default environment is defined
      by the env_locations class variable
    filename: str
      If filename is not None, the argument overrides the file defined
      in the env_locations dictionary. The filename should be a
      relative path to the "modules" directory.
    python: str
      If set, the specific Python version of the environment for the
      builder is used instead of the default. Current options are
      '27' and '36' for Python 2.7 and 3.6, respectively.
    copy: bool
      If set to True, the --copy flag is passed to conda
    offline: bool
      If set to True, the --offline flag is passed to conda
    """

    # handles check for choices in case parser is not available
    if builder not in self.env_locations:
      raise RuntimeError("""
The builder, {builder}, is not recognized. The available builders are,
{builders}
""".\
format(builder=builder, builders=', '.join(sorted(self.env_locations.keys()))))

    if self.conda_base is None:
      raise RuntimeError("""A conda installation is not available.""")

    if builder == "dials" and python in ("27", "36"):
      builder = "dials-old"

    if filename is None:
      filename = os.path.join(
        self.root_dir, 'modules', self.env_locations[builder])
      if python is not None:
        if python not in ['27', '37', '38', '39', '310', '311', '312', '313']:
          raise RuntimeError(
            """Only Python 2.7, 3.7, 3.8, 3.9, and 3.10 are currently supported.""")
        filename = filename.replace('PYTHON_VERSION', python)
    else:
      filename = os.path.abspath(filename)

    if not os.path.isfile(filename):
      raise RuntimeError("""\
The file, {filename}, is not available. Please contact the developers to make \
sure that the requested version of Python is supported for the {builder} \
builder.""".format(filename=filename, builder=builder))

    yaml_format = False
    if filename.endswith('yml') or filename.endswith('yaml'):
      yaml_format = True

    # make a new environment directory
    if self.conda_env is None:
      name = 'conda_base'
      prefix = os.path.join(self.root_dir, name)
    # or use the existing one
    else:
      prefix = os.path.abspath(self.conda_env)

    # compare time stamps of the filename and environment directory
    # only install/update if the time stamp of the filename is more recent
    file_stats = None
    env_stats = None
    if os.path.exists(filename):
      file_stats = os.stat(filename)
    if os.path.exists(prefix):
      env_stats = os.stat(prefix)

    if env_stats is not None and file_stats is not None:
      if env_stats.st_mtime > file_stats.st_mtime:
        print('The environment is newer than the environment file. Skipping update.',
              file=self.log)
        return

    # install a new environment or update and existing one
    if prefix in self.environments:
      command = 'install'
      if yaml_format:
        command = 'update'
      text_messages = ['Updating', 'update of']
    else:
      command = 'create'
      text_messages = ['Installing', 'installation into']
    command_list = [self.conda_exe, command, '--prefix', prefix,
                    '--file', filename]
    if yaml_format:
      command_list.insert(1, 'env')
    if self.system == 'Windows':
      command_list = [os.path.join(self.conda_base, 'Scripts', 'activate'),
                      'base', '&&'] + command_list
    if copy and not yaml_format:
      command_list.append('--copy')
    if offline and not yaml_format:
      command_list.append('--offline')
    if builder in ("dials", "dials-old", "xfel", "labelit") and not yaml_format:
      command_list.append("-y")
    if builder in self.env_without_python:
      python_version = tuple(int(i) for i in (python or "36"))
      python_requirement = '"conda-forge::python>=%s.%s,<%s.%s"' % (
          python_version[0],
          python_version[1],
          python_version[0],
          python_version[1] + 1,
      )
      command_list.append(python_requirement)
    # RuntimeError is raised on failure
    print('{text} {builder} environment with:\n  {filename}'.format(
          text=text_messages[0], builder=builder, filename=filename),
          file=self.log)

    self._retry_command(command_list, text_messages[1], prefix, verbose=True)

    # on Windows, also download the Visual C++ 2008 Redistributable
    # use the same version as conda-forge
    # https://github.com/conda-forge/vs2008_runtime-feedstock
    if self.system == 'Windows' and prefix.endswith('conda_base'):
      download_file(
        url='https://download.microsoft.com/download/5/D/8/5D8C65CB-C849-4025-8E95-C3966CAFD8AE/vcredist_x64.exe',
        filename=os.path.join(prefix, 'vcredist_x64.exe'))

# ---------------------------------------------------------------------------
  def create_dev_environment(self, svn=False, git=True):
    """
    Create a separate environment for development tools. Packages from
    the conda-forge channel are used.

    Package list:
      svn
      git
      git-lfs

    Parameters
    ----------
      svn: bool
        If set to true, svn is installed. This is useful for Xcode 11.4
        and later where svn is no longer available.
      git: bool
        If set to true, git is installed with git-lfs. This is needed
        for repositories that use git-lfs.

    Returns
    -------
      dev_dir: str
        The directory with the environment
    """

    package_list = []
    if svn:
      package_list.append('svn')
    if git:
      package_list.append('git')
      package_list.append('git-lfs')

    prefix = os.path.join(self.root_dir, 'dev_env')
    command = 'create'
    text_messages = ['Installing', 'installation into']
    if prefix in self.environments:
      command = 'update'
      text_messages = ['Updating', 'update of']

    command_list = [self.conda_exe, command, '-y', '-c', 'conda-forge',
                    '--prefix', prefix] + package_list

    print('-'*79, file=self.log)
    print('{text} extra development environment containing:'.format(text=text_messages[0]),
          file=self.log)
    for package in package_list:
      print('  -', package, file=self.log)

    self._retry_command(command_list, text_messages[1], prefix, verbose=True)

    print('-'*79, file=self.log)
    return prefix

  # ---------------------------------------------------------------------------
  def write_conda_setpaths(self, prefix='conda', build_dir=None,
                           conda_env=None, check_file=True):
    """
    Write a script similar to setpaths.sh/csh/bat that activates the environment
    first. This is useful for developers that want to use other software
    installed in the conda environment without having to manually activate it.

    Parameters
    ----------
      prefix: str
        The prefix for the output file. The output file will be constructed
        as <prefix>_<script name>.<extension>. The <script name> is the
        standard setpaths or unsetpaths file, and the <extension> will be
        the existing extensions for the script (sh, csh, and bat)
      build_dir: str
        The build directory where the standard setpaths script resides.
      conda_env: str
        The conda environment directory
      check_file: bool
        Used to override the check_file attribute

    Returns
    -------
      Nothing
    """

    if check_file is None:
      check_file = self.check_file

    # check that build_dir is valid
    if build_dir is None:
      build_dir = os.getenv('LIBTBX_BUILD')
    if build_dir is None:
      raise RuntimeError("""\
Please run the dispatcher version of libtbx.install_conda or provide a valid
directory as an argument to the write_conda_setpaths function.""")

    build_dir_error = """\
Please provide the directory to the setpaths script.
"""
    if build_dir is None or not os.path.isdir(build_dir):
      raise RuntimeError(build_dir_error)
    if sys.platform == 'win32':
      if not os.path.isfile(os.path.join(build_dir, 'setpaths.bat')):
        raise RuntimeError(build_dir_error)
    else:
      if not os.path.isfile(os.path.join(build_dir, 'setpaths.sh')) \
         or not os.path.isfile(os.path.join(build_dir, 'setpaths.csh')):
        raise RuntimeError(build_dir_error)

    # check conda_env
    if conda_env is None:
      conda_env = self.conda_env
    if conda_env is None:
      from libtbx.env_config import get_conda_prefix
      conda_env = get_conda_prefix()
    conda_env = os.path.abspath(conda_env)

    # -------------------------------------------------------------------------
    def do_check_file(filename):
      """
      Convenience function for checking if the script was written.
      """
      if os.path.isfile(filename):
        print('{filename} has been written successfully.'.format(
          filename=os.path.basename(filename)))
      else:
        raise RuntimeError("""
{filename} has not been written successfully.""".format(filename=filename))
    # -------------------------------------------------------------------------

    # Windows
    if sys.platform == 'win32':
      script_template = """\
@ECHO OFF
CALL {mc3_dir} {conda_env}
CALL {setpaths}
"""
      # activate
      mc3_dir = os.path.abspath(
        os.path.join(self.conda_base, 'Scripts', 'activate.bat'))
      setpaths = os.path.abspath(os.path.join(build_dir, 'setpaths.bat'))
      filename = os.path.abspath(os.path.join(build_dir, prefix + '_setpaths.bat'))
      with open(filename, 'w') as f:
        f.write(script_template.format(
          mc3_dir=mc3_dir, conda_env=conda_env, setpaths=setpaths))
      if check_file:
        do_check_file(filename)
      # deactivate
      mc3_dir = 'conda'
      conda_env = 'deactivate'
      setpaths = os.path.abspath(os.path.join(build_dir, 'unsetpaths.bat'))
      filename = os.path.abspath(
        os.path.join(build_dir, prefix + '_unsetpaths.bat'))
      with open(filename, 'w') as f:
        f.write(script_template.format(
          mc3_dir=mc3_dir, conda_env=conda_env, setpaths=setpaths))
      if check_file:
        do_check_file(filename)
    # linux and macOS
    else:
      for ext in ('sh', 'csh'):
        # activate
        script_template = """\
source {mc3_dir}
{get_old_prompt}
conda activate {conda_env}
{set_old_prompt}
{unset_old_prompt}
source {setpaths}
"""
        mc3_dir = os.path.abspath(
          os.path.join(self.conda_base, 'etc', 'profile.d', 'conda.' + ext))
        if ext == 'sh':
          get_old_prompt = 'LIBTBX_OLD_PS1=$PS1'
          set_old_prompt = 'PS1=$LIBTBX_OLD_PS1'
          unset_old_prompt = 'unset LIBTBX_OLD_PS1'
        else:
          get_old_prompt = 'set libtbx_old_prompt="$prompt"'
          set_old_prompt = 'set prompt="$libtbx_old_prompt"'
          unset_old_prompt = 'unset libtbx_old_prompt'
        setpaths = os.path.abspath(os.path.join(build_dir, 'setpaths.' + ext))
        filename = os.path.abspath(
          os.path.join(build_dir, prefix + '_setpaths.' + ext))
        with open(filename, 'w') as f:
          f.write(script_template.format(
            mc3_dir=mc3_dir, get_old_prompt=get_old_prompt, conda_env=conda_env,
            set_old_prompt=set_old_prompt, unset_old_prompt=unset_old_prompt,
            setpaths=setpaths))
        if check_file:
          do_check_file(filename)
        # deactivate
        script_template = """\
conda deactivate
source {setpaths}
"""
        setpaths = os.path.abspath(os.path.join(build_dir, 'unsetpaths.' + ext))
        filename = os.path.abspath(
          os.path.join(build_dir, prefix + '_unsetpaths.' + ext))
        with open(filename, 'w') as f:
          f.write(script_template.format(setpaths=setpaths))
        if check_file:
          do_check_file(filename)

# =============================================================================
def run():
  prog = os.environ.get('LIBTBX_DISPATCHER_NAME')
  if prog is None or prog.startswith('python') or prog.endswith('python'):
    prog = os.path.basename(sys.argv[0])

  epilog = """
Example usage:

  {prog}
    Shows this help screen

  {prog} --verbose
    Shows if a base conda installation can be found. If no installation
    is found, the latest miniconda3 will be installed.

  {prog} --install_conda --builder=<builder>
    Install conda and default environment for <builder>

  {prog} --install_conda --builder=<builder> --python=37
    Install conda and default Python 3.7 environment for <builder>

  {prog} --conda_base=<path> --builder=<builder>
    Install default environment for <builder> with known conda installation

  {prog} --conda_env=<path> --builder=<builder>
    Update conda environment for builder

""".format(prog=prog)

  parser = argparse.ArgumentParser(
    prog=prog, description=__doc__, epilog=epilog,
    formatter_class=argparse.RawDescriptionHelpFormatter)

  # CLI options
  parser.add_argument(
    '--builder', default=None, type=str, nargs='?', const='cctbx',
    choices=sorted(conda_manager.env_locations.keys()),
    help="""Install the default environment for a builder. The choices are the
      same as the ones for bootstrap.py. The default builder is "cctbx." """)
  parser.add_argument(
    '--python', default='37', type=str, nargs='?', const='37',
    choices=['27', '37', '38', '39', '310', '311', '312', '313'],
    help="""When set, a specific Python version of the environment will be used.
    This only affects environments selected with the --builder flag.""")
  parser.add_argument(
    '--install_conda', action='store_true',
    help="""When set, conda will be automatically downloaded and installed
      regardless of an existing installation.""")
  parser.add_argument(
    '--update_conda', action='store_true',
    help="""When set, conda will try to update itself to the latest version.
      This should only be used if your conda installation was installed by
      this script or if your conda is writeable and uses the "defaults"
      channel""")
  parser.add_argument(
    '--install_env', default=None, type=str, nargs='?', const='',
    metavar='ENV_FILE',
    help="""When set, the environment for the builder will be installed. The
      default environment can be overridden by providing a path to the
      environment file.""")
  parser.add_argument(
    '--install_dev_env', action='store_true',
    help="""When set, an additional environment named dev_env will be
    created that contains extra development tools (svn, git, git-lfs).""")
  parser.add_argument(
    '--conda_base', default=None, type=str, metavar='BASE_DIRECTORY',
    help="""The location of the base conda environment. This is useful for
      systems where there is a multiuser installation of conda. Providing the
      path will ensure the the correct conda is used. To determine the base
      location, run "conda info" and look for the "base environment" setting.
      This is required, otherwise, a new conda installation will be created.""")
  parser.add_argument(
    '--conda_env', default=None, type=str, metavar='ENV_DIRECTORY',
    help="""The location of the conda environment for building. This is useful
      for when the exact conda environment is known. Providing the path will
      ensure that that environment is used. Using $CONDA_PREFIX as the
      argument will use the currently active environment for building.""")
  parser.add_argument(
    '--write_setpaths', default=None, type=str, nargs='?', metavar='PREFIX',
    const='conda',
    help="""When set, another script is added that activates the conda
      environment before sourcing the setpaths script. Optionally, the prefix
      of the script can be provided. The output files will be named
      <prefix>_setpaths.<extension> and <prefix>_unsetpaths.<extension>
      where the extension will be "sh"/"csh" for linux and macOS, and
      "bat" for Windows.""")
  parser.add_argument(
    '--copy', action='store_true', default=False,
    help="""When set, the new environment has copies, not links to files. This
      should only be used when building installers.""")
  parser.add_argument(
    '--offline', action='store_true', default=False,
    help="""When set, the network will not be accessed for installing packages.
      This should only be used if the packages are already cached.""")
  parser.add_argument(
    '--verbose', action='store_true', default=False,
    help="""When set, output from conda is displayed. """)

  # show help if no arguments are provided
  if len(sys.argv) == 1:
    parser.print_help()
    parser.exit()

  namespace = parser.parse_args()

  # logic for installing conda
  conda_base = namespace.conda_base
  if namespace.install_conda:
    conda_base = None

  # logic for installing environment
  filename = None
  builder = namespace.builder
  if namespace.install_env is not None:
    tmp_filename = os.path.abspath(namespace.install_env)
    if os.path.isfile(tmp_filename):
      filename = tmp_filename
    if builder is None:
      builder = 'cctbx'

  m = conda_manager(root_dir=root_dir, conda_base=conda_base,
                    conda_env=namespace.conda_env, check_file=True,
                    verbose=namespace.verbose)

  # if --update_conda is set, try to update now
  if namespace.update_conda:
    m.update_conda()

  # if --install_dev_env is set, construct development environment
  if namespace.install_dev_env:
    m.create_dev_environment()

  # if builder is available, construct environment
  if builder is not None:
    m.create_environment(builder=builder, filename=filename,
                         python=namespace.python,
                         copy=namespace.copy, offline=namespace.offline)

  # if --write_setpaths is set, write the extra script
  if namespace.write_setpaths is not None:
    m.write_conda_setpaths(prefix=namespace.write_setpaths, check_file=True)

  return 0

# =============================================================================
if __name__ == '__main__':
  sys.exit(run())

# =============================================================================
# end
