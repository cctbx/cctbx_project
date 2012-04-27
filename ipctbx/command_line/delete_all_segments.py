def run(args):
  assert len(args) == 0
  import ipctbx
  n_deleted = ipctbx.delete_all_segments()
  print "Number of shared memory segments deleted:", n_deleted

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
