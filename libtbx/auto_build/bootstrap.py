
# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
from __future__ import division
import os, os.path, posixpath
import sys
import subprocess
import optparse
#import getpass
import urlparse
import shutil, tarfile

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
    if not self.kwargs.get("quiet", False):
      print "===== Running in %s:"%workdir, " ".join(command)
    if workdir:
      try:
        os.makedirs(workdir)
      except Exception, e:
        pass

    if command[0] == 'curl':
      # XXX Ugly hack: intercept attemps to spawn external curl.
      # There is no need to depend on curl since Python has urllib2.
      Downloader().download_to_file(command[1], os.path.join(workdir, command[3]))
      return
    if command[0] == 'tar':
      try:
        # XXX use tarfile rather than unix tar command which is not platform independent
        tar = tarfile.open(os.path.join(workdir, command[2]))
        tar.extractall(path=workdir)
        if len(command) > 3 and command[3]: # rename to expected folder name, e.g. boost_hot -> boost
          module = os.path.join(workdir, command[3])
          tarfoldername = os.path.join(workdir, os.path.commonprefix(tar.getnames()).split('/')[0])
          # only rename if folder names differ
          if module != tarfoldername and os.path.exists(module):
            shutil.rmtree(module)
          os.rename(tarfoldername, module)
        tar.close()
      except Exception, e:
        print "Extracting tar archive resulted in error:"
        raise
      return
    if command[0] == 'rm' :
      # XXX use shutil rather than rm which is not platform independent
      for dir in command[2:]:
        if os.path.exists(dir):
          print 'Deleting directory : %s' % dir
          shutil.rmtree(dir)
      return
    try:
      #print "workdir, os.getcwd =", workdir, os.getcwd()
      #if not os.path.isabs(command[0]):
        # executable path isn't located relative to workdir
      #  command[0] = os.path.join(workdir, command[0])
      p = subprocess.Popen(
        args=command,
        cwd=workdir,
        stdout=sys.stdout,
        stderr=sys.stderr
      )
    except Exception, e:
      if isinstance(e, OSError):
        if e.errno == 2:
          executable = os.path.normpath(os.path.join(workdir, command[0]))
          raise RuntimeError("Could not run %s: File not found" % executable)
      if 'child_traceback' in dir(e):
        print "Calling subprocess resulted in error; ", e.child_traceback
      raise e

    p.wait()
    if p.returncode != 0 and self.kwargs.get('haltOnFailure'):
      raise RuntimeError, "Process failed with return code %s"%(p.returncode)

# Download URL to local file
class Downloader(object):
  def download_to_file(self, url, file, log=sys.stdout, status=True,
                        buffer_until_file_downloaded=True):
    """Downloads a URL to file. Returns the file size.
       Returns -1 if the downloaded file size does not match the expected file
       size
       Returns -2 if the download is skipped due to the file at the URL not
       being newer than the local copy (with matching file sizes).
    """

    import os, time, urllib2
    socket = urllib2.urlopen(url)

    file_size = int(socket.info().getheader('Content-Length'))
    # There is no guarantee that the content-length header is set

    remote_mtime = 0
    try:
      remote_mtime = time.mktime(socket.info().getdate('last-modified'))
    except:
      pass

    if (file_size > 0):
      if (remote_mtime > 0):
        # check if existing file matches remote size and timestamp
        try:
          (mode, ino, dev, nlink, uid, gid, size, atime, mtime, ctime) = os.stat(file)
          if (size == file_size) and (remote_mtime == mtime):
            log.write("local copy is current\n")
            socket.close()
            return -2
        except:
          # proceed with download if timestamp/size check fails for any reason
          pass

      hr_size = (file_size, "B")
      if (hr_size[0] > 500): hr_size = (hr_size[0] / 1024, "kB")
      if (hr_size[0] > 500): hr_size = (hr_size[0] / 1024, "MB")
      log.write("%.1f %s\n" % hr_size)
      if status:
        log.write("    [0%")
        log.flush()

    received = 0
    block_size = 8192
    progress = 1
    if not buffer_until_file_downloaded:
      # Buffering a large file is very inefficient and python is likely to crash
      # Allow for writing the file immediately so we can empty the buffer
      f2 = open(file, "wb")
    data = ""
    while 1:
      block = socket.read(block_size)
      received += len(block)
      data += block
      if not buffer_until_file_downloaded:
        f2.write(data)
        data = ""
      if status and (file_size > 0):
        while (100 * received / file_size) > progress:
          progress += 1
          if (progress % 20) == 0:
            log.write("%d%%" % progress)
          elif (progress % 2) == 0:
            log.write(".")
          log.flush()

      if not block: break
    if status and (file_size > 0):
      log.write("]\n")
    else:
      log.write("%d kB\n" % (received / 1024))
    log.flush()
    if not buffer_until_file_downloaded:
      f2.close()
    socket.close()
    # Do not write out file during the download. If a download temporarily fails we
    # may still have a clean, working (yet older) copy of the file.
    if buffer_until_file_downloaded:
      f = open(file, "wb")
      f.write(data)
      f.close()

    if (file_size > 0) and (file_size != received):
      return -1

    if remote_mtime > 0:
      # set file timestamp if timestamp information is available
      from stat import ST_ATIME
      st = os.stat(file)
      atime = st[ST_ATIME] # current access time
      os.utime(file,(atime,remote_mtime))

    return received

class cleanup_ext_class(object):
  def __init__(self, filename_ext, workdir=None):
    self.filename_ext = filename_ext
    self.workdir = workdir

  def get_command(self):
    return "delete *%s in %s" % (self.filename_ext, self.workdir).split()

  def remove_ext_files(self):
    cwd=os.getcwd()
    if self.workdir is not None:
      if os.path.exists(self.workdir):
        os.chdir(self.workdir)
      else:
        return
    print "\n  removing %s files in %s" % (self.filename_ext, os.getcwd())
    i=0
    for root, dirs, files in os.walk(".", topdown=False):
      for name in files:
        if name.endswith(self.filename_ext):
          os.remove(os.path.join(root, name))
          i+=1
    os.chdir(cwd)
    print "  removed %d files" % i

  def run(self):
    self.remove_ext_files()

##### Modules #####
class SourceModule(object):
  _modules = {}
  module = None
  authenticated = None
  authenticatedWindows = None
  authentarfile = None
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
    if self.authenticatedWindows and sys.platform == 'win32':
      return [self.authenticatedWindows[0], self.authenticatedWindows[1]%auth]
    if not self.authenticated:
      return None
    return [self.authenticated[0], self.authenticated[1]%auth]

  def get_tarauthenticated(self, auth=None):
    auth = auth or {}
    if self.authentarfile and sys.platform == 'win32':
      return [self.authentarfile[0]%auth, self.authentarfile[1], self.authentarfile[2]]
    return None, None, None

  def get_anonymous(self):
    return self.anonymous

# Core external repositories
# The trailing slashes ARE significant.
# These must all provide anonymous access.
# On Windows due to absence of rsync we use pscp from the Putty programs.
class ccp4io_module(SourceModule):
  module = 'ccp4io'
  anonymous = ['curl', 'http://cci.lbl.gov/repositories/ccp4io.gz']
  authentarfile = ['%(cciuser)s@cci.lbl.gov', 'ccp4io.tar.gz', '/net/cci/auto_build/repositories/ccp4io']
  authenticated = ['rsync', '%(cciuser)s@cci.lbl.gov:/net/cci/auto_build/repositories/ccp4io/']

class annlib_module(SourceModule):
  module = 'annlib'
  anonymous = ['curl', 'http://cci.lbl.gov/repositories/annlib.gz']
  #authenticatedWindows = anonymous
  authentarfile = ['%(cciuser)s@cci.lbl.gov', 'annlib.tar.gz', '/net/cci/auto_build/repositories/annlib']
  authenticated = ['rsync', '%(cciuser)s@cci.lbl.gov:/net/cci/auto_build/repositories/annlib/']

class scons_module(SourceModule):
  module = 'scons'
  anonymous = ['curl', 'http://cci.lbl.gov/repositories/scons.gz']
  #authenticatedWindows = anonymous
  authentarfile = ['%(cciuser)s@cci.lbl.gov', 'scons.tar.gz', '/net/cci/auto_build/repositories/scons']
  authenticated = ['rsync', '%(cciuser)s@cci.lbl.gov:/net/cci/auto_build/repositories/scons/']

class boost_module(SourceModule):
  module = 'boost'
  anonymous = ['curl', 'http://cci.lbl.gov/repositories/boost.gz']
  # Compared to rsync pscp is very slow when downloading multiple files
  # Resort to downloading the compressed archive on Windows
  #authenticatedWindows = anonymous
  authentarfile = ['%(cciuser)s@cci.lbl.gov', 'boost_hot.tar.gz', '/net/cci/auto_build/repositories/boost_hot/']
  authenticated = ['rsync', '%(cciuser)s@cci.lbl.gov:/net/cci/auto_build/repositories/boost_hot/']

class libsvm_module(SourceModule):
  module = 'libsvm'
  anonymous = ['curl', 'http://cci.lbl.gov/repositories/libsvm.gz']
  authenticated = ['rsync', '%(cciuser)s@cci.lbl.gov:/net/cci/auto_build/repositories/libsvm/']

# Core CCTBX repositories
# These must all provide anonymous access.
class cctbx_module(SourceModule):
  module = 'cctbx_project'
  anonymous = ['svn','svn://svn.code.sf.net/p/cctbx/code/trunk']
  authenticated = ['svn', '%(sfmethod)s://%(sfuser)s@svn.code.sf.net/p/cctbx/code/trunk']

class cbflib_module(SourceModule):
  module = 'cbflib'
  anonymous = ['svn', 'svn://svn.code.sf.net/p/cbflib/code-0/trunk/CBFlib_bleeding_edge']
  authenticated = ['svn', '%(sfmethod)s://%(sfuser)s@svn.code.sf.net/p/cbflib/code-0/trunk/CBFlib_bleeding_edge']

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

class amber_module(SourceModule):
  module = 'amber_adaptbx'
  authenticated = ['svn', 'svn+ssh://%(cciuser)s@cci.lbl.gov/amber_adaptbx/trunk']

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
  # Compared to rsync pscp is very slow when downloading multiple files
  # Resort to downloading a compressed archive on Windows. Must create it first
  authentarfile = ['%(cciuser)s@cci.lbl.gov', 'phaser.tar.gz', '/net/cci/auto_build/repositories/phaser']
  authenticated = ['rsync', '%(cciuser)s@cci.lbl.gov:/net/cci/auto_build/repositories/phaser/']

class phaser_regression_module(SourceModule):
  module = 'phaser_regression'
  # Compared to rsync pscp is very slow when downloading multiple files
  # Resort to downloading a compressed archive on Windows. Must create it first
  authentarfile = ['%(cciuser)s@cci.lbl.gov', 'phaser_regression.tar.gz', '/net/cci/auto_build/repositories/phaser_regression']
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
  authenticated = ['svn', '%(sfmethod)s://%(sfuser)s@svn.code.sf.net/p/dials/code/trunk']

class dials_regression_module(SourceModule):
  module = 'dials_regression'
  authenticated = ['svn', 'svn+ssh://%(cciuser)s@cci.lbl.gov/dials_regression/trunk']

class xfel_regression_module(SourceModule):
  module = 'xfel_regression'
  authenticated = ['svn', 'svn+ssh://%(cciuser)s@cci.lbl.gov/xfel_regression/trunk']

class xia2_module(SourceModule):
  module = 'xia2'
  anonymous = ['svn', 'svn://svn.code.sf.net/p/xia2/code/trunk/xia2']
  authenticated = ['svn', '%(sfmethod)s://%(sfuser)s@svn.code.sf.net/p/xia2/code/trunk/xia2']

class xia2_regression_module(SourceModule):
  module = 'xia2_regression'
  anonymous = ['svn', 'svn://svn.code.sf.net/p/xia2/code/trunk/xia2_regression']
  authenticated = ['svn', '%(sfmethod)s://%(sfuser)s@svn.code.sf.net/p/xia2/code/trunk/xia2_regression']

# Duke repositories
class probe_module(SourceModule):
  module = 'probe'
  anonymous = ['svn', 'https://github.com/rlabduke/probe.git/trunk']

class suitename_module(SourceModule):
  module = 'suitename'
  anonymous = ['svn', 'https://github.com/rlabduke/suitename.git/trunk']

class reduce_module(SourceModule):
  module = 'reduce'
  anonymous = ['svn', 'https://github.com/rlabduke/reduce.git/trunk']

class king_module(SourceModule):
  module = 'king'
  anonymous = ['svn', 'https://github.com/rlabduke/phenix_king_binaries.git/trunk']

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
      python_base=None,
      cleanup=False,
      hot=True,
      update=True,
      base=True,
      build=True,
      tests=True,
      distribute=False,
      auth=None,
      with_python=None,
      nproc=4,
      verbose=False,
      download_only=False,
      skip_base="",
      force_base_build=False,
    ):
    if nproc is None:
      self.nproc=4
    else:
      self.nproc=nproc
    """Create and add all the steps."""
    # self.cciuser = cciuser or getpass.getuser()
    self.set_auth(auth)
    self.steps = []
    self.category = category
    self.platform = platform
    self.name = '%s-%s'%(self.category, self.platform)
    # Platform configuration.
    self.python_base = self.opjoin(*['..', 'base', 'bin', 'python'])
    if 'win32'==sys.platform:
      self.python_base = self.opjoin(*[os.getcwd(), 'base', 'bin', 'python', 'python.exe'])
    self.with_python = with_python
    if self.with_python:
      self.python_base = with_python
    self.verbose = verbose
    self.download_only = download_only
    self.skip_base = skip_base
    self.force_base_build = force_base_build

    self.add_init()

    # Cleanup
    if cleanup:
      self.cleanup(['dist', 'tests', 'doc', 'tmp', 'base', 'base_tmp', 'build'])
    else:
      self.cleanup(['dist', 'tests', 'doc', 'tmp'])

    # Add 'hot' sources
    if hot:
      map(self.add_module, self.get_hot())

    # Add svn sources.
    if update:
      map(self.add_module, self.get_codebases())

    # always remove .pyc files
    self.remove_pyc()

    # Build base packages
    if base:
      self.add_base()

    # Configure, make
    if build and not self.download_only:
      self.add_configure()
      self.add_make()
      self.add_install()

    # Tests, tests
    if tests and not self.download_only:
      self.add_tests()

    # Distribute
    if distribute and not self.download_only:
      self.add_distribute()

  def add_auth(self, account, username):
    self.auth[account] = username

  def set_auth(self, auth):
    self.auth = auth or {}

  def get_auth(self):
    return self.auth

  def remove_pyc(self):
    self.add_step(cleanup_ext_class(".pyc", "modules"))

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
    #return self.sep.join(args)
    return os.path.join(*args)

  def get_codebases(self):
    return self.CODEBASES + self.CODEBASES_EXTRA

  def get_hot(self):
    return self.HOT + self.HOT_EXTRA

  def get_libtbx_configure(self):
    if sys.platform == "win32": # we can't currently compile cbflib for Windows
      return list(set(self.LIBTBX + self.LIBTBX_EXTRA) - set(['cbflib']))
    return self.LIBTBX + self.LIBTBX_EXTRA

  def add_init(self):
    pass

  def cleanup(self, dirs=None):
    dirs = dirs or []
    cmd=['rm', '-rf'] + dirs
    self.add_step(self.shell(
      name='cleanup',
      command =cmd,
      workdir=['.']
    ))

  def add_step(self, step):
    """Add a step."""
    self.steps.append(step)
    if 0:
      print "commands "*8
      for step in self.steps:
        print step
        try:    print " ".join(step.get_command())
        except: print '????'
      print "commands "*8

  def add_module(self, module):
    method, url = MODULES.get_module(module)().get_url(auth=self.get_auth())
    tarurl, arxname, dirpath = MODULES.get_module(module)().get_tarauthenticated(auth=self.get_auth())
    if sys.platform == "win32":
      if module in ["cbflib",]: # can't currently compile cbflib for Windows due to lack of HDF5 component
        return
    if method == 'rsync' and sys.platform != "win32":
      self._add_rsync(module, url)
    elif sys.platform == "win32" and method == 'pscp':
      self._add_pscp(module, url)
    elif sys.platform == "win32" and tarurl:
        self._add_remote_make_tar(module, tarurl, arxname, dirpath)
        self._add_pscp(module, tarurl + ':' + arxname)
        self._add_remote_rm_tar(module, tarurl, arxname)
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

  def _add_remote_make_tar(self, module, tarurl, arxname, dirpath):
    """Add packages not in source control."""
    # tar up hot packages for quick file transfer to windows since there's no rsync and pscp is painfully slow
    if dirpath[-1] == '/':
      dirpath = dirpath[:-1]
    basename = posixpath.basename(dirpath)
    self.add_step(self.shell( # pack directory with tar on remote system but exclude all svn files
      name='hot %s'%module,
      command=[
        'plink',
        tarurl,
        'cd',
        posixpath.split(dirpath)[0],
        '&&',
        'tar',
        'cfz',
        '~/' + arxname,
        basename + '/*',
        '--exclude',
        '*.svn'
      ],
      workdir=['modules']
    ))

  def _add_remote_rm_tar(self, module, tarurl, arxname):
    self.add_step(self.shell( # delete the tarfile on remote system
      name='hot %s'%module,
      command=[
        'plink',
        tarurl,
        'rm ',
        arxname
      ],
      workdir=['modules']
    ))
    self.add_step(self.shell( # extract the tarfile on our system
      command=['tar', 'xzf', arxname, module],
      workdir=['modules']
    ))

  def _add_pscp(self, module, url):
    """Add packages not in source control."""
    url1 = url
    if url[-1] == '/':
      url1 = url[:-1]
    self.add_step(self.shell( # copy files/directory recursively from remote system
      name='hot %s'%module,
      command=[
        'pscp',
        '-r',
        url1,
        '.',
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
      command=['tar', 'xzf', filename],
      workdir=['modules']
    ))

  def _add_svn(self, module, url):
    if os.path.exists(self.opjoin(*['modules', module, '.svn'])):
      # print "using update..."
      self.add_step(self.shell(
          command=['svn', 'update', module],
          workdir=['modules']
      ))
      self.add_step(self.shell(
          command=['svn', 'status', module],
          workdir=['modules'],
          quiet=True,
      ))
    elif os.path.exists(self.opjoin(*['modules', module])):
      print "Existing non-svn directory -- dont know what to do. skipping: %s"%module
    else:
      # print "fresh checkout..."
      self.add_step(self.shell(
          command=['svn', 'co', url, module],
          workdir=['modules']
      ))

  def _add_git(self, module, url):
    pass

  def add_command(self, command, name=None, workdir=None, args=None, **kwargs):
    if sys.platform == 'win32':
      command = command + '.bat'
    # Relative path to workdir.
    workdir = workdir or ['build']
    dots = [".."]*len(workdir)
    if workdir[0] == '.':
      dots = []
    if sys.platform == 'win32':
      dots.extend([os.getcwd(), 'build', 'bin', command])
    else:
      dots.extend(['build', 'bin', command])
    self.add_step(self.shell(
      name=name or command,
      command=[self.opjoin(*dots)] + (args or []),
      workdir=workdir,
      **kwargs
    ))

  def add_test_command(self,
                       command,
                       name=None,
                       workdir=None,
                       args=None,
                       haltOnFailure=False,
                       **kwargs
                       ):
    if name is None: name='test %s'%command
    self.add_command(
      command,
      name=name,
      workdir=(workdir or ['tests', command]),
      args=args,
      haltOnFailure=haltOnFailure,
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
  def add_base(self, extra_opts=[]):
    """Build the base dependencies, e.g. Python, HDF5, etc."""
    if self.with_python:
      extra_opts = ['--with-python', self.with_python]
    if self.verbose:
      extra_opts.append('-v')
    if self.download_only:
      extra_opts.append('--download-only')
    if self.skip_base:
      extra_opts.append('--skip-base=%s' % self.skip_base)
    if not self.force_base_build:
      extra_opts.append("--skip-if-exists")
    self.add_step(self.shell(
      name='base',
      command=[
        'python',
        self.opjoin('modules', 'cctbx_project', 'libtbx', 'auto_build', 'install_base_packages.py'),
        '--python-shared',
        '--nproc=%s' %str(self.nproc),
        '--%s'%self.BASE_PACKAGES
      ] + extra_opts,
      workdir=['.']
    ))

  def add_configure(self):
    self.add_step(self.shell(command=[
        self.python_base, # default to using our python rather than system python
        self.opjoin('..', 'modules', 'cctbx_project', 'libtbx', 'configure.py')
        ] + self.get_libtbx_configure(),
      workdir=['build'],
      description="run configure.py",
    ))
    # Prepare saving configure.py command to file should user want to manually recompile Phenix
    relpython = os.path.relpath(self.python_base, self.opjoin(os.getcwd(),'build'))
    if sys.platform != 'win32':
      relpython = os.path.relpath(self.python_base, os.getcwd())
    configcmd =[
        relpython, # default to using our python rather than system python
        self.opjoin('..', 'modules', 'cctbx_project', 'libtbx', 'configure.py')
        ] + self.get_libtbx_configure()
    fname = self.opjoin("config_modules.cmd")
    confstr = subprocess.list2cmdline(configcmd) + '\n'
    if sys.platform != 'win32':
      fname = self.opjoin("config_modules.sh")
      confstr = '#!/bin/sh\n\n' + confstr
    # klonky way of writing file later on, but it works
    self.add_step(self.shell(command=[
         'python','-c','open(r\"%s\",\"w\").write(r\"\"\"%s\"\"\")' %(fname, confstr)
         ],
      workdir=['build'],
      description="Saving configure.py command to file",
    ))


  def add_make(self):
    # Todo: nproc=auto
    self.add_command('libtbx.scons', args=['-j', str(self.nproc)])

  def add_install(self):
    """Run after compile, before tests."""
    self.add_command('mmtbx.rebuild_rotarama_cache',
                     name="rebuild rotarama",
    )

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
    #"libsvm",
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
  def add_tests(self):
#    self.add_step(cleanup_ext_class(".pyc", "modules"))
    self.add_test_command('libtbx.import_all_ext')
    self.add_test_command('libtbx.import_all_python', workdir=['modules', 'cctbx_project'])
    self.add_test_command('cctbx_regression.test_nightly')

class DIALSBuilder(CCIBuilder):
  CODEBASES_EXTRA = ['dials', 'xia2']
  LIBTBX_EXTRA = ['dials', 'xia2']
  def add_tests(self):
    self.add_test_command('libtbx.import_all_ext')
    self.add_test_command('cctbx_regression.test_nightly')
    self.add_test_parallel('dials')

  def add_base(self):
    super(DIALSBuilder, self).add_base(
      extra_opts=['--dials',
                  #'--wxpython3'
                 ])

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
   'cxi_xdr_xes',
   'prime'
 ]
 def add_tests(self):
    self.add_test_command('libtbx.import_all_ext')
    self.add_test_command('cctbx_regression.test_nightly')
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
    'amber_adaptbx',
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
    # 'king',
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
    'amber_adaptbx',
    'reduce',
    'probe',
  ]
  def add_install(self):
    Builder.add_install(self)
    self.add_command('phenix_html.rebuild_docs')

  def add_tests(self):
    # Include cctbx tests.
    self.add_test_command('libtbx.import_all_ext')
    self.add_test_command('cctbx_regression.test_nightly')
    # Windows convenience hack.
    if 'win32' == sys.platform:
      self.add_test_command('phenix_regression.test_nightly_windows')
    else:
      self.add_test_command('phenix_regression.test_nightly')
    # Other Phenix tests.
    self.add_test_parallel(module='elbow')
    self.add_test_command('phenix_html.rebuild_docs')
    self.add_test_command('phenix_regression.run_p9_sad_benchmark')
    self.add_test_command('phenix_regression.run_hipip_refine_benchmark')

def run(root=None):
  usage = """Usage: %prog [options] [actions]

  You may specify one or more actions:
    hot - Update static sources (boost, scons, etc.)
    update - Update source repositories (cctbx, cbflib, etc.)
    base - Build base dependencies (python, hdf5, wxWidgets, etc.)
    build - Build
    tests - Run tests

  The default action is to run: hot, update, base, build

  You can specify which package will be downloaded, configured,
  and built with "--builder". Current builders:
    cctbx, phenix, xfel, dials, labelit

  You can provide your SourceForge username with "--sfuser", and
  your CCI SVN username with "--cciuser". These will checkout
  and update repositories with your credentials. Some builders,
  like phenix, require this argument for access to certain
  repositories.

  You can provide the number of processes to use in compilation
  using "--nproc".
  Complete build output is shown with "-v" or "--verbose".

  Finally, you may specify a specific Python interpreter
  using "--with-python".

  Example:

    python bootstrap.py --builder=cctbx --sfuser=ianrees hot update build tests

  """
  parser = optparse.OptionParser(usage=usage)
  # parser.add_option("--root", help="Root directory; this will contain base, modules, build, etc.")
  parser.add_option("--builder", help="Builder: cctbx, phenix, xfel, dials, labelit", default="cctbx")
  parser.add_option("--cciuser", help="CCI SVN username.")
  parser.add_option("--sfuser", help="SourceForge SVN username.")
  parser.add_option("--sfmethod", help="SourceForge SVN checkout method.", default="svn+ssh")
  parser.add_option("--with-python", dest="with_python", help="Use specified Python interpreter")
  parser.add_option("--nproc", help="# processes in compile step.")
  parser.add_option("--download-only", dest="download_only", action="store_true", help="Do not build, only download prerequisites", default=False)
  parser.add_option("-v", "--verbose", dest="verbose", action="store_true", help="Verbose output", default=False)
  parser.add_option("--skip-base-packages",
                    dest="skip_base",
                    action="store",
                    default="")
  parser.add_option("--force-base-build",
                    dest="force_base_build",
                    action="store_true",
                    default=False)
  options, args = parser.parse_args()

  # Root dir
  # options.root = options.root or root

  # Check actions
  allowedargs = ['cleanup', 'hot', 'update', 'base', 'build', 'tests']
  args = args or ['hot', 'update', 'base', 'build']
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
  if options.sfmethod:
    auth['sfmethod'] = options.sfmethod

  # Build
  builder = builders[options.builder]
  builder(
    category=options.builder,
    platform='dev',
    with_python=options.with_python,
    auth=auth,
    hot=('hot' in actions),
    update=('update' in actions),
    base=('base' in actions),
    build=('build' in actions),
    tests=('tests' in actions),
    cleanup=("cleanup" in actions),
    nproc=options.nproc,
    verbose=options.verbose,
    download_only=options.download_only,
    skip_base=options.skip_base,
    force_base_build=options.force_base_build,
  ).run()
  print "\nBootstrap success: %s" % ", ".join(actions)

if __name__ == "__main__":
  run()
