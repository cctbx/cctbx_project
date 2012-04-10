cd $cctbxroot

##############################################################################

${XCODE_CCTBX_BUILD}/bin/python -u <<PYSCRIPT
import sys, os, re
import libtbx.load_env
from libtbx import easy_run
from libtbx.command_line import scons

print; print; print '-'*80

m = re.search(r"^(\d+ \. \d+ \. \d+) .*? \[\s*GCC\s* (\d+ \. \d+ \. \d+)",
              sys.version,
              re.X|re.M|re.S)
print 'Python %s (compiled with gcc %s)' % m.groups()
print

print "** %s **" % libtbx.env.build_path.basename()
libtbx.env.build_options.report()

print
print '-'*80
for compiler in ('gcc', 'clang',):
  print easy_run.fully_buffered('type %s' % compiler).stdout_lines[0]
print '-'*80

os.chdir(os.environ['XCODE_CCTBX_BUILD'])
sys.argv[1:] = os.environ['XCODE_SCONS_OPTIONS'].split()
if os.environ.get('XCODE_SCONS_LIB_TARGET'):
  libs = os.environ['XCODE_SCONS_LIB_TARGET'].split()
  sys.argv.extend([ "lib/%s.so" % lib for lib in libs ])
  print "warning: only library being built: %s" % ', '.join(libs)
elif os.environ.get('XCODE_SCONS_PROGRAM_TARGET'):
  programs = os.environ['XCODE_SCONS_PROGRAM_TARGET'].split()
  sys.argv.extend(programs)
  print "warning: only program being built: %s" % ', '.join(programs)
scons.run()
PYSCRIPT
