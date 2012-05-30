
import libtbx.test_utils.parallel
from libtbx.utils import Sorry
import libtbx.phil
import os
import sys

master_phil = libtbx.phil.parse("""
dir_name = None
  .type = path
nproc = 1
  .type=  int
""")

def run (args) :
  cwd = os.getcwd()
  cwd_files = os.listdir(cwd)
  if (len(cwd_files) > 0) :
    raise Sorry("Please run this program in an empty directory.")
  user_phil = []
  for arg in args :
    if os.path.isdir(arg) :
      user_phil.append(libtbx.phil.parse("dir_name=%s" % arg))
    else :
      try :
        arg_phil = libtbx.phil.parse(arg)
      except RuntimeError :
        raise Sorry("Unrecognized argument '%s'" % arg)
      else :
        user_phil.append(arg_phil)
  params = master_phil.fetch(sources=user_phil).extract()
  assert (params.dir_name is not None)
  all_tests = libtbx.test_utils.parallel.find_tests(params.dir_name)
  if (len(all_tests) == 0) :
    raise Sorry("No test scripts found in %s." % params.dir_name)
  print "Running the following %d tests on %d processors:" % (len(all_tests),
    params.nproc)
  for test in all_tests :
    print "  " + test
  log = open("zlog", "wb")
  libtbx.test_utils.parallel.run_all_tests(
    files=all_tests,
    nprocs=params.nproc,
    log=log)
  log.close()

if (__name__ == "__main__") :
  run(sys.argv[1:])
