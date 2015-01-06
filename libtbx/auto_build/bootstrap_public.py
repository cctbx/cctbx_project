from bootstrap import optparse,os
from bootstrap import ShellCommand, HOT, CODEBASES, ModuleManager
from bootstrap import DIALSBuilder, PHENIXBuilder, XFELBuilder, LABELITBuilder, CCTBXBuilder, Builder as CCIBuilder

HOT['annlib']='curl http://cci.lbl.gov/repositories/annlib.gz > annlib.gz'
HOT['scons']='curl http://cci.lbl.gov/repositories/scons.gz > scons.gz'
HOT['boost']='curl http://cci.lbl.gov/repositories/boost.gz > boost.gz'
HOT['ccp4io']='curl http://cci.lbl.gov/repositories/ccp4io.gz > ccp4io.gz'
HOT['docutils']='curl http://cci.lbl.gov/repositories/docutils.gz > docutils.gz'
HOT['annlib_adaptbx']='curl http://cci.lbl.gov/repositories/annlib_adaptbx.gz > annlib_adaptbx.gz'
HOT['ccp4io_adaptbx']='curl http://cci.lbl.gov/repositories/ccp4io_adaptbx.gz > ccp4io_adaptbx.gz'
HOT['clipper']='curl http://cci.lbl.gov/repositories/clipper.gz > clipper.gz'
HOT['opt_resources']='curl http://cci.lbl.gov/repositories/opt_resources.gz > opt_resources.gz'
HOT['gui_resources']='curl http://cci.lbl.gov/repositories/gui_resources.gz > gui_resources.gz'
HOT['tntbx']='curl http://cci.lbl.gov/repositories/tntbx.gz > tntbx.gz'


def add_hot(self, package):
    """Add packages not in source control."""
    # rsync the hot packages.
    print package, HOT[package]

    self.add_step(self.shell(
      name='hot %s'%package,
      command=HOT[package],
      workdir=['modules']
    ))
    # If it's a tarball, unzip it.
    self.add_step(self.shell(
        name='hot %s untar'%package,
        command=['tar', '-xzf', os.path.join("%s.gz"%package)],
        workdir=['modules']
      ))
def get_hot(self):
  return [
    'annlib',
    'boost',
    'scons',
    'ccp4io',
    'docutils',
    'gui_resources',
    'ccp4io_adaptbx',
    'annlib_adaptbx',
    'tntbx',
    'clipper',
  ]
def get_codebases(self):
  return [
    'cctbx_project',
    'dials',
    'cbflib',
  ]
def add_make(self):
  self.add_step(self.shell(command=['make'],
    workdir=['build']
  ))

def get_libtbx_configure(self):
  return [
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
    'dials'
    ]

CCIBuilder.add_hot = add_hot
CCIBuilder.get_hot = get_hot
CCIBuilder.get_codebases = get_codebases
CCIBuilder.add_make = add_make
CCIBuilder.get_libtbx_configure = get_libtbx_configure
if __name__ == "__main__":
  parser = optparse.OptionParser()
  parser.add_option("--builder", help="Builder: cctbx, phenix, xfel, dials, labelit", default="cctbx")
  parser.add_option("--cciuser", help="CCI SVN username.")
  options, args = parser.parse_args()

  # Check actions
  allowedargs = ['cleanup', 'hot', 'update', 'base', 'build', 'install', 'tests']
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
    'phenix': PHENIXBuilder,
    'xfel': XFELBuilder,
    'labelit': LABELITBuilder,
    'dials': DIALSBuilder
  }
  if options.builder not in builders:
    raise ValueError("Unknown builder: %s"%options.builder)

  # Build
  builder = builders[options.builder]
  builder(
    category='cctbx',
    platform='debug',
    cciuser=options.cciuser,
    hot=('hot' in actions),
    update=('update' in actions),
    base=('base' in actions),
    build=('build' in actions),
    install=('install' in actions),
    tests=('tests' in actions)
  ).run()
