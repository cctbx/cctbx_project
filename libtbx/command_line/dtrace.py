import sys, os
import libtbx.load_env

def run():
  if len(sys.argv) < 2: help()
  libtbx_path = libtbx.env.find_in_repositories('libtbx')
  dtrace_directory = os.path.join(libtbx_path, 'dtrace')
  dtrace_script = sys.argv[1]
  if (not os.path.isfile(dtrace_script)
      and not os.path.dirname(dtrace_script)):
    if dtrace_script[-2:] != '.d': dtrace_script += '.d'
    dtrace_script = os.path.join(dtrace_directory, dtrace_script)
  if not os.path.isfile(dtrace_script): help()
  os.execlp("dtrace", "dtrace",
            "-Z", # essential because the Python probes are dynamically
                  # created after Python has launched to run the script
            "-s", dtrace_script,
            "-c", "%s %s" % (sys.executable, " ".join(sys.argv[2:])))

def help():
  print "libtbx.dtrace dtrace-script python-script [args ...]"
  print "\tdtrace-script can be either the path of a dtrace script"
  print "\tor if there is no file at that path, a dtrace script"
  print "\tin libtbx/dtrace. In the latter case, the extension .d can"
  print "\tomitted on the command line."
  sys.exit(1)

if __name__ == '__main__':
  run()
