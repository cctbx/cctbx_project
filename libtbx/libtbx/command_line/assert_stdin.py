import sys

def run():
  input = " ".join([line.rstrip() for line in sys.stdin.readlines()])
  pattern = " ".join(sys.argv[1:])
  if (input != pattern):
    sys.tracebacklimit = 0
    raise AssertionError('"%s" != "%s"' % (input, pattern))

if (__name__ == "__main__"):
  run()
