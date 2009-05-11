from scitbx.array_family import flex
import sys

def run(args):
  assert len(args) == 0
  print "\n".join(flex.show_sizes_int())
  print
  print "\n".join(flex.show_sizes_double())
  print

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
