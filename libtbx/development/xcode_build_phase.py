from __future__ import absolute_import, division, print_function
import sys, os, re

def run_scons():
  import libtbx.load_env
  from libtbx import easy_run
  from libtbx.command_line import scons

  print(); print(); print('-'*80)

  m = re.search(r"^(\d+ \. \d+ \. \d+) .*? \[\s*GCC\s* (\d+ \. \d+ \. \d+)",
                sys.version,
                re.X|re.M|re.S)
  print('Python %s (compiled with gcc %s)' % m.groups())
  print()

  print("** %s **" % libtbx.env.build_path.basename())
  libtbx.env.build_options.report()

  print()
  print('-'*80)
  for compiler in ('gcc', 'clang',):
    print(easy_run.fully_buffered('type %s' % compiler).stdout_lines[0])
  print('-'*80)

  os.chdir(os.environ['XCODE_CCTBX_BUILD'])
  sys.argv[1:] = os.environ['XCODE_SCONS_OPTIONS'].split()
  if os.environ.get('XCODE_SCONS_LIB_TARGET'):
    libs = os.environ['XCODE_SCONS_LIB_TARGET'].split()
    sys.argv.extend([ "lib/%s.so" % lib for lib in libs ])
    print("warning: only library being built: %s" % ', '.join(libs))
  elif os.environ.get('XCODE_SCONS_PROGRAM_TARGET'):
    programs = os.environ['XCODE_SCONS_PROGRAM_TARGET'].split()
    sys.argv.extend(programs)
    print("warning: only program being built: %s" % ', '.join(programs))
  scons.run()

def run_filtered_scons():
  import libtbx.load_env
  import subprocess
  proc = subprocess.Popen([ "%s/libtbx.python" % abs(libtbx.env.bin_path),
                            '-u', __file__, 'child'],
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
  error_pat = re.compile(r'^(.*? : \d* : \d* : \s* error:)', re.X)
  prefix = abs(libtbx.env.build_path) + '/'
  for li in proc.stdout:
    m = error_pat.search(li)
    if m and not m.group(1).startswith('/'):
      sys.stdout.write(prefix)
    sys.stdout.write(li)
    sys.stdout.flush()
  proc.wait()
  return proc.returncode

if __name__ == '__main__':
  if sys.argv[1] == "child":
    run_scons()
  elif sys.argv[1] == "parent":
    sys.exit(run_filtered_scons())
  else:
    print("Usage: %s [child|parent]" % os.path.basename(__file__))
