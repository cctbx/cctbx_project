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
  create_environment(builder, filename, copy)
    Uses the known conda installtion to create an environment
  """

  # Currently, there is one monolithic environment
  # The key names for this dictionary must match the builder names in
  # bootstrap.py
  env_locations = {
    'cctbxlite': default_file,
    'cctbx': default_file,
    'phenix': os.path.join('phenix', 'conda_envs',
      default_format.format(builder='phenix', version=version,
                            platform=conda_platform[platform.system()])),
    'phenix_tng': os.path.join('phenix', 'conda_envs',
      default_format.format(builder='phenix', version=version,
                            platform=conda_platform[platform.system()])),
    'xfel': default_file,
    'labelit': default_file,
    'dials': default_file,
    'external': default_file,
    'molprobity': default_file,
    'qrefine': default_file,
    'phaser': default_file,
    'phaser_tng': default_file
  }

  # ---------------------------------------------------------------------------
  def __init__(self, root_dir=root_dir, conda_base=None, conda_env=None,
               check_file=True, verbose=False, log=sys.stdout):
    """
    Constructor that performs that basic check for the conda installation.
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
      that want to manage their own environments. If a conda environment
      is active, this will be automatically set to the CONDA_PREFIX
      environment variable.
    check_file: bool
      Flag for checking if a file exists. A RuntimeError is raised if a
      this flag is set and a file does not exist. Used in get_conda_exe
      and get_conda_python.
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
    env_variable = os.environ.get('CONDA_PREFIX')
    if env_variable is not None:
      self.conda_env = env_variable

    self.check_file = check_file
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
      if conda_info['conda_version'] < '4.4':
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
      for i, env in enumerate(paths):
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
          for dir in dirs:
            dir = os.path.join(env_dir, dir)
            if os.path.isdir(dir):
              environments.add(dir)

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
      flags = '-b -u -p {install_dir}'.format(install_dir=install_dir)
      command_list = ['/bin/sh', filename, flags]
    print('Installing miniconda to {install_dir}'.format(
      install_dir=install_dir), file=self.log)
    output = check_output(command_list, env=self.env)
    if self.verbose:
      print(output, file=self.log)

    return install_dir

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

    if filename is None:
      filename = os.path.join(
        self.root_dir, 'modules', self.env_locations[builder])
      if python is not None:
        if python not in ['27', '36']:
          raise RuntimeError(
            """Only Python 2.7 and 3.6 are currently supported.""")
        filename = filename.replace('PYTHON_VERSION', python)
    else:
      filename = os.path.abspath(filename)

    if not os.path.isfile(filename):
      raise RuntimeError("""The file, {filename}, is not available""".\
                         format(filename=filename))

    # make a new environment directory
    if self.conda_env is None:
      name = 'conda_base'
      prefix = os.path.join(self.root_dir, name)
    # or use the existing one
    else:
      prefix = os.path.abspath(self.conda_env)

    # install a new environment or update and existing one
    if prefix in self.environments:
      command = 'install'
      text_messages = ['Updating', 'update of']
    else:
      command = 'create'
      text_messages = ['Installing', 'installation into']
    command_list = [self.conda_exe, command, '--prefix', prefix,
                    '--file', filename]
    if self.system == 'Windows':
      command_list = [os.path.join(self.conda_base, 'Scripts', 'activate'),
                      'base', '&&'] + command_list
    if copy:
      command_list.append('--copy')
    if offline:
      command_list.append('--offline')
    # RuntimeError is raised on failure
    print('{text} {builder} environment with:\n  {filename}'.format(
          text=text_messages[0], builder=builder, filename=filename),
          file=self.log)
    output = check_output(command_list, env=self.env)
    if self.verbose:
      print(output, file=self.log)
    print('Completed {text}:\n  {prefix}'.format(text=text_messages[1],
          prefix=prefix), file=self.log)

    # on Windows, also download the Visual C++ 2008 Redistributable
    # use the same version as conda-forge
    # https://github.com/conda-forge/vs2008_runtime-feedstock
    if self.system == 'Windows' and prefix.endswith('conda_base'):
      download_file(
        url='https://download.microsoft.com/download/5/D/8/5D8C65CB-C849-4025-8E95-C3966CAFD8AE/vcredist_x64.exe',
        filename=os.path.join(prefix, 'vcredist_x64.exe'))

    # check that environment file is updated
    self.environments = self.update_environments()
    if prefix not in self.environments:
      raise RuntimeError("""
The newly installed environment cannot be found in
${HOME}/.conda/environments.txt.
""")

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

  {prog} --install_conda --builder=<builder> --python=36
    Install conda and default Python 3.6 environment for <builder>

  {prog} --conda_base=<path> --builder=<builder>
    Install default environment for <builder> with known conda installation

  {prog} --conda_env=<path> --builder=<builder>
    Update conda environment for builder

  In an active conda environment, the last example can just be

    {prog} --builder=<builder>

  because $CONDA_PREFIX is known.

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
    '--python', default='27', type=str, nargs='?', const='27',
    choices=['27', '36'],
    help="""When set, a specific Python version of the environment will be used.
    This only affects environments selected with the --builder flag.""")
  parser.add_argument(
    '--install_conda', action='store_true',
    help="""When set, conda will be automatically downloaded and installed
      regardless of an existing installation.""")
  parser.add_argument(
    '--install_env', default=None, type=str, nargs='?', const='',
    metavar='ENV_FILE',
    help="""When set, the environment for the builder will be installed. The
      default environment can be overridden by providing a path to the
      environment file.""")
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
      ensure that that environment is used. Also, if a conda environment
      is currently active, the CONDA_PREFIX environment variable will be
      used as this location.""")
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

  if builder is not None:
    m.create_environment(builder=builder, filename=filename,
                         python=namespace.python,
                         copy=namespace.copy, offline=namespace.offline)

# =============================================================================
if __name__ == '__main__':
  run()

# =============================================================================
# end
