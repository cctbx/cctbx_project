import sys

def run(args):
  assert len(args) == 2, "file1 file2"
  from libtbx.test_utils import show_diff
  texts = ["\n".join(open(arg).read().splitlines()) for arg in args]
  assert not show_diff(texts[0], texts[1])

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
