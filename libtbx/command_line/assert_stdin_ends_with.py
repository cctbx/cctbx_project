import sys

def run(args):
  assert len(args) == 1
  expected_last_line = args[0]
  assert sys.stdin.read().splitlines()[-1] == expected_last_line

if (__name__ == "__main__"):
  run(sys.argv[1:])
