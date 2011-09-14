import sys, os

def show(out):
  var_names = os.environ.keys()
  var_names.sort()
  for var_name in var_names:
    print >> out, "%s=%s" % (var_name, os.environ[var_name])

def run(args):
  assert len(args) == 0
  show(out=sys.stdout)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
