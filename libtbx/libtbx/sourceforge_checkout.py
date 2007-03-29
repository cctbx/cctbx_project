# minimal imports to enable checkout even if libtbx sources are removed already
import sys, os

def usage():
  this = os.environ.get("LIBTBX_DISPATCHER_NAME")
  if (this is None):
    this = 'python "%s"' % sys.argv[0].replace('"', '\\"')
  print
  print "usage: %s project [--branch=name] [module]..." % this
  print
  print "example: %s cctbx libtbx" % this
  print
  sys.stdout.flush()

def show_and_run(command):
  print command
  sys.stdout.flush()
  os.system(command)

def run(args):
  branch = "trunk"
  remaining_args = []
  arg_iter = iter(args)
  for arg in arg_iter:
    if (arg in ["-h", "--help"]):
      return usage()
    if (arg.startswith("--branch")):
      if (arg[8] == "="):
        branch = arg[9:]
      elif (len(arg) > 8):
        return usage()
      else:
        try: branch = arg_iter.next()
        except StopIteration: return usage()
      branch = branch.strip()
      if (len(branch) == 0):
        return usage()
    else:
      remaining_args.append(arg)
  if (len(remaining_args) == 0):
    return usage()
  project = remaining_args[0]
  modules = remaining_args[1:]
  base_url = "https://%s.svn.sourceforge.net/svnroot/%s/%s" % (
    project, project, branch)
  if (len(modules) == 0):
    show_and_run(command="svn list %s" % base_url)
    return
  for module in modules:
    show_and_run(command="svn checkout %s/%s" % (base_url, module))
    print
    sys.stdout.flush()

if (__name__ == "__main__"):
  run(sys.argv[1:])
