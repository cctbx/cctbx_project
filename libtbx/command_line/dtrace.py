import sys, os
import libtbx.load_env

def run():
  if len(sys.argv) < 2: help()
  libtbx_path = libtbx.env.find_in_repositories('libtbx')
  dtrace_directory = os.path.join(libtbx_path, 'dtrace')
  dtrace_script = sys.argv[1]
  params = []
  python_script_and_its_arguments = []
  for opt in sys.argv[2:]:
    if '=' in opt:
      params.append("-D%s" % opt)
    else:
      python_script_and_its_arguments.append(opt)
  python_script_and_its_arguments = " ".join(python_script_and_its_arguments)
  if (not os.path.isfile(dtrace_script)
      and not os.path.dirname(dtrace_script)):
    if dtrace_script[-2:] != '.d': dtrace_script += '.d'
    dtrace_script = os.path.join(dtrace_directory, dtrace_script)
  if not os.path.isfile(dtrace_script): help()
  args = [ "dtrace",
           "-Z", # essential because the Python probes are dynamically
                 # created after Python has launched to run the script
           "-C",
           "-s", dtrace_script,
         ] + params + [
            "-c", "%s %s" % (sys.executable, python_script_and_its_arguments)
         ]
  os.execvp("dtrace", args)


def help():
  print "usage: libtbx.dtrace dtrace_script [ param=value ...] ",
  print "python_script [args ...]"
  print
  print "       dtrace_script can be either the path of a dtrace script"
  print "       or if there is no file at that path, a dtrace script"
  print "       in libtbx/dtrace. In the latter case, the extension .d can"
  print "       omitted on the command line."
  sys.exit(1)

if __name__ == '__main__':
  run()
