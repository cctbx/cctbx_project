# -*- python -*-
from __future__ import division
import os
import sys
import subprocess
import optparse
import getpass
import urllib
import urlparse

# To download this file:
# svn export svn://svn.code.sf.net/p/cctbx/code/trunk/libtbx/auto_build/bootstrap.py

# Note: to relocate an SVN repo:
# svn relocate svn+ssh://<username>@svn.code.sf.net/p/cctbx/code/trunk

# Mock commands to run standalone, without buildbot.
class ShellCommand(object):
  def __init__(self, **kwargs):
    self.kwargs = kwargs

  def get_command(self):
    return self.kwargs['command']

  def get_workdir(self):
    return self.kwargs.get('workdir', 'build')

  def run(self):
    command = self.get_command()
    workdir = self.get_workdir()
    print "===== Running in %s:"%workdir, " ".join(command)
    if workdir:
      try:
        os.makedirs(workdir)
      except Exception, e:
        pass
    p = subprocess.Popen(
      args=command,
      cwd=workdir,
      stdout=sys.stdout,
      stderr=sys.stderr
    )
    p.wait()
    if p.returncode != 0 and self.kwargs.get('haltOnFailure'):
      raise RuntimeError, "Process failed with return code %s"%(p.returncode)

##### Modules #####

class SourceModule(object):
  _modules = {}
  module = None
  authenticated = None
  anonymous = None
  def __init__(self):
    if not self._modules:
      self.update_subclasses()

  def items(self):
    return self._modules.items()

  @classmethod
  def update_subclasses(cls):
    for i in cls.__subclasses__():
      cls._modules[i.module] = i

  def get_module(self, module):
    if module in self._modules:
      return self._modules[module]
    raise KeyError, "Unknown module: %s"%module

  def get_url(self, auth=None):
    repo = None
    try:
      repo = self.get_authenticated(auth=auth)
    except KeyError, e:
      repo = self.get_anonymous()
      if not repo:
        raise Exception('No anonymous access method defined for module: %s. Try with --%s'%(self.module, e.args[0]))
    repo = repo or self.get_anonymous()
    if not repo:
      raise Exception('No access method defined for module: %s'%self.module)
    return repo

  def get_authenticated(self, auth=None):
    auth = auth or {}
    if not self.authenticated:
      return None
    return [self.authenticated[0], self.authenticated[1]%auth]

  def get_anonymous(self):
    return self.anonymous

# Core external repositories
# The trailing slashes ARE significant.
# These must all provide anonymous access.
class ccp4io_module(SourceModule):
  module = 'ccp4io'
  anonymous = ['curl', 'http://cci.lbl.gov/repositories/ccp4io.gz']
  authenticated = ['rsync', '%(cciuser)s@cci.lbl.gov:/net/cci/auto_build/repositories/ccp4io/']

class annlib_module(SourceModule):
  module = 'annlib'
  anonymous = ['curl', 'http://cci.lbl.gov/repositories/annlib.gz']
  authenticated = ['rsync', '%(cciuser)s@cci.lbl.gov:/net/cci/auto_build/repositories/annlib/']

class scons_module(SourceModule):
  module = 'scons'
  anonymous = ['curl', 'http://cci.lbl.gov/repositories/scons.gz']
  authenticated = ['rsync', '%(cciuser)s@cci.lbl.gov:/net/cci/auto_build/repositories/scons/']

class boost_module(SourceModule):
  module = 'boost'
  anonymous = ['curl', 'http://cci.lbl.gov/repositories/boost.gz']
  authenticated = ['rsync', '%(cciuser)s@cci.lbl.gov:/net/cci/auto_build/repositories/boost_hot/']

# Core CCTBX repositories
# These must all provide anonymous access.
class cctbx_module(SourceModule):
  module = 'cctbx_project'
  anonymous = ['svn','svn://svn.code.sf.net/p/cctbx/code/trunk']
  authenticated = ['svn', 'svn+ssh://%(sfuser)s@svn.code.sf.net/p/cctbx/code/trunk']

class cbflib_module(SourceModule):
  module = 'cbflib'
  anonymous = ['svn', 'svn://svn.code.sf.net/p/cbflib/code-0/trunk/CBFlib_bleeding_edge']
  authenticated = ['svn', 'svn+ssh://%(sfuser)s@svn.code.sf.net/p/cbflib/code-0/trunk/CBFlib_bleeding_edge']

class ccp4io_adaptbx(SourceModule):
  module = 'ccp4io_adaptbx'
  anonymous = ['curl', 'http://cci.lbl.gov/repositories/ccp4io_adaptbx.gz']
  authenticated = ['svn', 'svn+ssh://%(cciuser)s@cci.lbl.gov/ccp4io_adaptbx/trunk']

class annlib_adaptbx(SourceModule):
  module = 'annlib_adaptbx'
  anonymous = ['curl', 'http://cci.lbl.gov/repositories/annlib_adaptbx.gz']
  authenticated = ['svn', 'svn+ssh://%(cciuser)s@cci.lbl.gov/annlib_adaptbx/trunk']

class tntbx_module(SourceModule):
  module = 'tntbx'
  anonymous = ['curl', 'http://cci.lbl.gov/repositories/tntbx.gz']
  authenticated = ['svn', 'svn+ssh://%(cciuser)s@cci.lbl.gov/tntbx/trunk']

class clipper_module(SourceModule):
  module = 'clipper'
  anonymous = ['curl', 'http://cci.lbl.gov/repositories/clipper.gz']
  authenticated = ['svn', 'svn+ssh://%(cciuser)s@cci.lbl.gov/clipper/trunk']

class gui_resources_module(SourceModule):
  module = 'gui_resources'
  anonymous = ['curl', 'http://cci.lbl.gov/repositories/gui_resources.gz']
  authenticated = ['svn', 'svn+ssh://%(cciuser)s@cci.lbl.gov/gui_resources/trunk']

class opt_resources_module(SourceModule):
  module = 'opt_resources'
  authenticated = ['svn', 'svn+ssh://%(cciuser)s@cci.lbl.gov/opt_resources/trunk']

# Phenix repositories
class phenix_module(SourceModule):
  module = 'phenix'
  authenticated = ['svn', 'svn+ssh://%(cciuser)s@cci.lbl.gov/phenix/trunk']

class phenix_html(SourceModule):
  module = 'phenix_html'
  authenticated = ['svn', 'svn+ssh://%(cciuser)s@cci.lbl.gov/phenix_html/trunk']

class phenix_examples(SourceModule):
  module = 'phenix_examples'
  authenticated = ['svn', 'svn+ssh://%(cciuser)s@cci.lbl.gov/phenix_examples/trunk']

class phenix_regression(SourceModule):
  module = 'phenix_regression'
  authenticated = ['svn', 'svn+ssh://%(cciuser)s@cci.lbl.gov/phenix_regression/trunk']

class plex_module(SourceModule):
  module = 'Plex'
  authenticated = ['svn', 'svn+ssh://%(cciuser)s@cci.lbl.gov/Plex/trunk']

class pyquante_module(SourceModule):
  module = 'PyQuante'
  authenticated = ['svn', 'svn+ssh://%(cciuser)s@cci.lbl.gov/PyQuante/trunk']

class chem_data_module(SourceModule):
  module = 'chem_data'
  authenticated = ['svn', 'svn+ssh://%(cciuser)s@cci.lbl.gov/chem_data/trunk']

class elbow_module(SourceModule):
  module = 'elbow'
  authenticated = ['svn', 'svn+ssh://%(cciuser)s@cci.lbl.gov/elbow/trunk']

class ksdssp_module(SourceModule):
  module = 'ksdssp'
  authenticated = ['svn', 'svn+ssh://%(cciuser)s@cci.lbl.gov/ksdssp/trunk']

class pex_module(SourceModule):
  module = 'pex'
  authenticated = ['svn', 'svn+ssh://%(cciuser)s@cci.lbl.gov/pex/trunk']

class pulchra_module(SourceModule):
  module = 'pulchra'
  authenticated = ['svn', 'svn+ssh://%(cciuser)s@cci.lbl.gov/pulchra/trunk']

class solve_resolve_module(SourceModule):
  module = 'solve_resolve'
  authenticated = ['svn', 'svn+ssh://%(cciuser)s@cci.lbl.gov/solve_resolve/trunk']

class reel_module(SourceModule):
  module = 'reel'
  authenticated = ['svn', 'svn+ssh://%(cciuser)s@cci.lbl.gov/reel/trunk']

class muscle_module(SourceModule):
  module = 'muscle'
  authenticated = ['svn', 'svn+ssh://%(cciuser)s@cci.lbl.gov/muscle/trunk']

class cxi_xdr_xes_module(SourceModule):
  module = 'cxi_xdr_xes'
  authenticated = ['svn', 'svn+ssh://%(cciuser)s@cci.lbl.gov/cxi_xdr_xes/trunk']

class buildbot_module(SourceModule):
  module = 'buildbot'
  authenticated = ['svn', 'svn+ssh://%(cciuser)s@cci.lbl.gov/buildbot/trunk']

# Phaser repositories
class phaser_module(SourceModule):
  module = 'phaser'
  authenticated = ['rsync', '%(cciuser)s@cci.lbl.gov:/net/cci/auto_build/repositories/phaser/']

class phaser_regression_module(SourceModule):
  module = 'phaser_regression'
  authenticated = ['rsync', '%(cciuser)s@cci.lbl.gov:/net/cci/auto_build/repositories/phaser_regression/']

# DIALS repositories
class labelit_module(SourceModule):
  module = 'labelit'
  authenticated = ['svn', 'svn+ssh://%(cciuser)s@cci.lbl.gov/labelit/trunk']

class labelit_regression_module(SourceModule):
  module = 'labelit_regression'
  authenticated = ['svn', 'svn+ssh://%(cciuser)s@cci.lbl.gov/labelit_regression/trunk']

class dials_module(SourceModule):
  module = 'dials'
  anonymous = ['svn', 'svn://svn.code.sf.net/p/dials/code/trunk']
  authenticated = ['svn', 'svn+ssh://%(sfuser)s@svn.code.sf.net/p/dials/code/trunk']

class dials_regression_module(SourceModule):
  module = 'dials_regression'
  authenticated = ['svn', 'svn+ssh://%(cciuser)s@cci.lbl.gov/dials_regression/trunk']

class xfel_regression_module(SourceModule):
  module = 'xfel_regression'
  authenticated = ['svn', 'svn+ssh://%(cciuser)s@cci.lbl.gov/xfel_regression/trunk']

# Duke repositories
class probe_module(SourceModule):
  module = 'probe'
  anonymous = ['svn', 'https://github.com/rlabduke/probe/trunk']

class suitename_module(SourceModule):
  module = 'suitename'
  anonymous = ['svn', 'https://github.com/rlabduke/suitename/trunk']

class reduce_module(SourceModule):
  module = 'reduce'
  anonymous = ['svn', 'https://quiddity.biochem.duke.edu/svn/reduce/trunk']

class king_module(SourceModule):
  module = 'king'
  anonymous = ['svn', 'https://quiddity.biochem.duke.edu/svn/phenix/king']

MODULES = SourceModule()

###################################
##### Base Configuration      #####
###################################

class Builder(object):
  """Create buildbot configurations for CCI and CCTBX-like software."""
  # Base packages
  BASE_PACKAGES = 'all'
  # Checkout these codebases
  CODEBASES = ['cctbx_project']
  CODEBASES_EXTRA = []
  # Copy these sources from cci.lbl.gov
  HOT = []
  HOT_EXTRA = []
  # Configure for these cctbx packages
  LIBTBX = ['cctbx']
  LIBTBX_EXTRA = []

  def __init__(
      self,
      category=None,
      platform=None,
      sep=None,
      python_system=None,
      python_base=None,
      cleanup=False,
      hot=True,
      update=True,
      base=True,
      build=True,
      tests=True,
      distribute=False,
      auth=None
    ):
    """Create and add all the steps."""
    # self.cciuser = cciuser or getpass.getuser()
    self.set_auth(auth)
    self.steps = []
    self.category = category
    self.platform = platform
    self.name = '%s-%s'%(self.category, self.platform)
    # Windows convenience hack.
    if 'windows' in self.platform:
      base = False
      sep = sep or '\\'
      python_system = python_system or ['python']
      python_base = python_base or ['..', 'base', 'Python', 'python.exe']
    # Platform configuration.
    self.sep = sep or os.sep
    self.python_system = self.opjoin(*(python_system or ['python']))
    self.python_base = self.opjoin(*(python_base or ['..', 'base', 'bin', 'python']))

    self.add_init()

    # Cleanup
    if cleanup:
      self.cleanup(['dist', 'tests', 'doc', 'tmp', 'build'])
    else:
      self.cleanup(['dist', 'tests', 'doc', 'tmp'])

    # Add 'hot' sources
    if hot:
      map(self.add_module, self.get_hot())

    # Add svn sources.
    if update:
      map(self.add_module, self.get_codebases())

    # Build base packages
    if base:
      self.add_base()

    # Configure, make
    if build:
      self.add_configure()
      self.add_make()
      self.add_install()

    # Tests, tests
    if tests:
      self.add_tests()

    # Distribute
    if distribute:
      self.add_distribute()

  def add_auth(self, account, username):
    self.auth[account] = username

  def set_auth(self, auth):
    self.auth = auth or {}

  def get_auth(self):
    return self.auth

  def shell(self, **kwargs):
    # Convenience for ShellCommand
    kwargs['haltOnFailure'] = kwargs.pop('haltOnFailure', True)
    kwargs['description'] = kwargs.get('description') or kwargs.get('name')
    kwargs['timeout'] = 60*60*2 # 2 hours
    if 'workdir' in kwargs:
      kwargs['workdir'] = self.opjoin(*kwargs['workdir'])
    return ShellCommand(**kwargs)

  def run(self):
    for i in self.steps:
      i.run()

  def opjoin(self, *args):
    return self.sep.join(args)

  def get_codebases(self):
    return self.CODEBASES + self.CODEBASES_EXTRA

  def get_hot(self):
    return self.HOT + self.HOT_EXTRA

  def get_libtbx_configure(self):
    return self.LIBTBX + self.LIBTBX_EXTRA

  def add_init(self):
    pass

  def cleanup(self, dirs=None):
    dirs = dirs or []
    self.add_step(self.shell(
      name='cleanup',
      command=['rm', '-rf'] + dirs,
      workdir=['.']
    ))

  def add_step(self, step):
    """Add a step."""
    self.steps.append(step)

  def add_module(self, module):
    method, url = MODULES.get_module(module)().get_url(auth=self.get_auth())
    if method == 'rsync':
      self._add_rsync(module, url)
    elif method == 'curl':
      self._add_curl(module, url)
    elif method == 'svn':
      self._add_svn(module, url)
    elif method == 'git':
      self._add_git(module, url)
    else:
      raise Exception('Unknown access method: %s %s'%(method, url))

  def _add_rsync(self, module, url):
    """Add packages not in source control."""
    # rsync the hot packages.
    self.add_step(self.shell(
      name='hot %s'%module,
      command=[
        'rsync',
        '-aL',
        '--delete',
        url,
        module,
      ],
      workdir=['modules']
    ))

  def _add_curl(self, module, url):
    filename = urlparse.urlparse(url).path.split('/')[-1]
    self.add_step(self.shell(
      command=['curl', url, '-o', filename],
      workdir=['modules'],
    ))
    self.add_step(self.shell(
      command=['tar', '-xvz', '-f', filename],
      workdir=['modules']
    ))

  def _add_svn(self, module, url):
    self.add_step(self.shell(
        command=['svn', 'co', url, module],
        workdir=['modules']
    ))

  def _add_git(self, module, url):
    pass

  def add_command(self, command, name=None, workdir=None, args=None, **kwargs):
    # Relative path to workdir.
    workdir = workdir or ['build']
    dots = [".."]*len(workdir)
    if workdir[0] == '.':
      dots = []
    dots.extend(['build', 'bin', command])
    self.add_step(self.shell(
      name=name or command,
      command=[self.opjoin(*dots)] + (args or []),
      workdir=workdir,
      **kwargs
    ))

  def add_test_command(self, command, name=None, workdir=None, args=None, **kwargs):
    self.add_command(
      command,
      name='test %s'%command,
      workdir=(workdir or ['tests', command]),
      haltOnFailure=False,
      **kwargs
    )

  def add_test_parallel(self, module=None):
    self.add_command(
      'libtbx.run_tests_parallel',
      name='test %s'%module,
      workdir=['tests', module],
      args=['module=%s'%module, 'nproc=auto', 'verbosity=1'],
      haltOnFailure=False
    )

  # Override these methods.
  def add_base(self):
    """Build the base dependencies, e.g. Python, HDF5, etc."""
    self.add_step(self.shell(
      name='base',
      command=[
        self.python_system,
        self.opjoin('modules', 'cctbx_project', 'libtbx', 'auto_build', 'install_base_packages.py'),
        '--python-shared',
        '--skip-if-exists',
        '--%s'%self.BASE_PACKAGES
      ],
      workdir=['.']
    ))

  def add_configure(self):
    self.add_step(self.shell(command=[
        self.python_base,
        self.opjoin('..', 'modules', 'cctbx_project', 'libtbx', 'configure.py')
        ] + self.get_libtbx_configure(),
      workdir=['build']
    ))

  def add_make(self):
    # Todo: nproc=auto
    self.add_command('libtbx.scons', args=['-j','4'])

  def add_install(self):
    """Run after compile, before tests."""
    pass

  def add_tests(self):
    """Run the unit tests."""
    pass

  def add_distribute(self):
    pass

##### Specific Configurations ######

class CCIBuilder(Builder):
  """Base class for packages that include CCTBX as a dependency."""
  # Base packages
  BASE_PACKAGES = 'all'
  # Checkout these codebases
  CODEBASES = [
    'cbflib',
    'cctbx_project',
    'gui_resources',
    'ccp4io_adaptbx',
    'annlib_adaptbx',
    'tntbx',
    'clipper'
  ]
  CODEBASES_EXTRA = []
  # Copy these sources from cci.lbl.gov
  HOT = [
    'annlib',
    'boost',
    'scons',
    'ccp4io',
  ]
  HOT_EXTRA = []
  # Configure for these cctbx packages
  LIBTBX = [
    'cctbx',
    'cbflib',
    'scitbx',
    'libtbx',
    'iotbx',
    'mmtbx',
    'smtbx',
    'dxtbx',
    'gltbx',
    'wxtbx',
  ]
  LIBTBX_EXTRA = []

##### CCTBX-derived packages #####

class CCTBXBuilder(CCIBuilder):
  BASE_PACKAGES = 'cctbx'
  CODEBASES_EXTRA = ['phenix_regression']
  LIBTBX_EXTRA = ['phenix_regression']
  def add_tests(self):
    self.add_test_command('libtbx.import_all_ext')
    self.add_test_command('libtbx.import_all_python', workdir=['modules', 'cctbx_project'])
    self.add_test_command('cctbx_regression.test_nightly')

class DIALSBuilder(CCIBuilder):
  CODEBASES_EXTRA = ['dials',]
  LIBTBX_EXTRA = ['dials',]
  def add_tests(self):
    self.add_test_parallel('dials')

class LABELITBuilder(CCIBuilder):
  CODEBASES_EXTRA = ['labelit', 'labelit_regression']
  LIBTBX_EXTRA = ['labelit', 'labelit_regression']
  def add_tests(self):
    pass

class XFELBuilder(CCIBuilder):
 CODEBASES_EXTRA = [
   'dials',
   'labelit',
   'labelit_regression',
   'xfel_regression',
   'cxi_xdr_xes'
 ]
 LIBTBX_EXTRA = [
   'dials',
   'labelit',
   'labelit_regression',
   'xfel',
   'xfel_regression',
   'cxi_xdr_xes'
 ]
 def add_tests(self):
   self.add_test_parallel('xfel_regression')

class PhenixBuilder(CCIBuilder):
  CODEBASES_EXTRA = [
    'chem_data',
    'phenix',
    'phenix_regression',
    'phenix_html',
    'phenix_examples',
    'labelit',
    'Plex',
    'PyQuante',
    'elbow',
    'ksdssp',
    'pex',
    'pulchra',
    'solve_resolve',
    'reel',
    'gui_resources',
    'opt_resources',
    'muscle',
    'labelit',
    'reduce',
    'probe',
    'king',
    'suitename',
  ]
  HOT_EXTRA = [
    'phaser',
    'phaser_regression',
  ]
  LIBTBX_EXTRA = [
    'chem_data',
    'phenix',
    'phenix_regression',
    'phenix_examples',
    'solve_resolve',
    'reel',
    'phaser',
    'phaser_regression',
    'labelit',
    'elbow',
    'reduce',
    'probe',
  ]
  def add_install(self):
    self.add_command('mmtbx.rebuild_rotarama_cache')
    self.add_command('phenix_html.rebuild_docs')

  def add_tests(self):
    # Windows convenience hack.
    if 'windows' in self.platform:
      self.add_test_command('phenix_regression.test_nightly_windows')
    else:
      self.add_test_command('phenix_regression.test_nightly')
    # Other Phenix tests.
    self.add_test_parallel(module='elbow')
    self.add_test_command('phenix_html.rebuild_docs')
    self.add_test_command('phenix_regression.run_p9_sad_benchmark')
    self.add_test_command('phenix_regression.run_hipip_refine_benchmark')

if __name__ == "__main__":
  usage = """Usage: %prog [options] [actions]
  
  You may specify one or more actions:
    hot - Update static sources (boost, scons, etc.)
    update - Update source repositories (cctbx, cbflib, etc.)
    base - Build base dependencies (python, hdf5, wxWidgets, etc.)
    build - Build
    tests - Run tests
    
  By default, all actions (hot, update, base, build, tests)
  will be run.  
    
  You can specify which package will be downloaded, configured,
  and built with "--builder". Current builders:
    cctbx, phenix, xfel, dials, labelit
    
  You can provide your SourceForge username with "--sfuser", and 
  your CCI SVN username with "--cciuser". These will checkout
  and update repositories with your credentials. Some builders,
  like phenix, require this argument for access to certain
  repositories.
  
  Example:
  
    python bootstrap.py --builder=cctbx --sfuser=ianrees hot update build tests

  """
  parser = optparse.OptionParser(usage=usage)
  parser.add_option("--builder", help="Builder: cctbx, phenix, xfel, dials, labelit", default="cctbx")
  parser.add_option("--cciuser", help="CCI SVN username.")
  parser.add_option("--sfuser", help="SourceForge SVN username.")
  options, args = parser.parse_args()

  # Check actions
  allowedargs = ['cleanup', 'hot', 'update', 'base', 'build', 'tests']
  args = args or ['hot', 'update', 'base', 'build', 'tests']
  actions = []
  for arg in args:
    if arg not in allowedargs:
      raise ValueError("Unknown action: %s"%arg)
  for arg in allowedargs:
    if arg in args:
      actions.append(arg)
  print "Performing actions:", " ".join(actions)

  # Check builder
  builders = {
    'cctbx': CCTBXBuilder,
    'phenix': PhenixBuilder,
    'xfel': XFELBuilder,
    'labelit': LABELITBuilder,
    'dials': DIALSBuilder
  }
  if options.builder not in builders:
    raise ValueError("Unknown builder: %s"%options.builder)

  auth = {}
  if options.cciuser:
    auth['cciuser'] = options.cciuser
  if options.sfuser:
    auth['sfuser'] = options.sfuser

  # Build
  builder = builders[options.builder]
  builder(
    category='cctbx',
    platform='debug',
    auth=auth,
    hot=('hot' in actions),
    update=('update' in actions),
    base=('base' in actions),
    build=('build' in actions),
    tests=('tests' in actions)
  ).run()
