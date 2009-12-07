import libtbx.load_env
import sys

def run(args):
  input = "\n".join(sys.stdin.read().splitlines()) + "\n"
  for arg in args:
    s = "\n".join(arg.splitlines())
    if (input.find(s) >= 0): continue
    print "BEGIN OF INPUT (%s)" % libtbx.env.dispatcher_name
    print "v" * 79
    sys.stdout.write(input)
    print "^" * 79
    print "END OF INPUT (%s)" % libtbx.env.dispatcher_name
    raise AssertionError("String not in input: %s" % arg)

if (__name__ == "__main__"):
  run(sys.argv[1:])
