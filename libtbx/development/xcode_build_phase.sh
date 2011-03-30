cd $cctbxroot

##############################################################################

${LIBTBX_BUILD}/bin/python <<PYSCRIPT
import sys, os, re
import libtbx.load_env
from libtbx.command_line import scons

print; print; print '-'*80

m = re.search(r"^(\d+ \. \d+ \. \d+) .*? \[\s*GCC\s* (\d+ \. \d+ \. \d+)",
              sys.version,
              re.X|re.M|re.S)
print 'Python %s (compiled with gcc %s)' % m.groups()
print

print "** %s **" % os.path.basename(libtbx.env.build_path)
libtbx.env.build_options.report()

print
print '-'*80
os.chdir(os.environ['LIBTBX_BUILD'])
sys.argv[1:] = os.environ['SCONS_OPTIONS'].split()
if os.environ.get('SCONS_LIB_TARGET'):
  libs = os.environ['SCONS_LIB_TARGET'].split()
  sys.argv.extend([ "lib/%s.so" % lib for lib in libs ])
  print "warning: only library being built: %s" % ', '.join(libs)
elif os.environ.get('SCONS_PROGRAM_TARGET'):
  programs = os.environ['SCONS_PROGRAM_TARGET'].split()
  sys.argv.extend(programs)
  print "warning: only program being built: %s" % ', '.join(programs)
scons.run()
PYSCRIPT
