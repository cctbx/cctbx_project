def exercise(n, use_dumps=False):
  from libtbx import easy_pickle
  import time
  obj = []
  for i in xrange(n):
    obj.append([i,i])
  for dgz in ["", ".gz"]:
    t0 = time.time()
    if (use_dumps):
      print "dumps/loads"
      pickle_string = easy_pickle.dumps(obj=obj)
    else:
      file_name = "test.dat"+dgz
      print file_name
      easy_pickle.dump(file_name=file_name, obj=obj)
    print "  dump: %.2f s" % (time.time()-t0)
    del obj
    t0 = time.time()
    if (use_dumps):
      obj = easy_pickle.loads(pickle_string)
    else:
      obj = easy_pickle.load(file_name=file_name)
    print "  load buffered: %.2f s" % (time.time()-t0)
    if (use_dumps):
      break
    else:
      del obj
      t0 = time.time()
      obj = easy_pickle.load(
        file_name=file_name, faster_but_using_more_memory=False)
      print "  load direct: %.2f s" % (time.time()-t0)

def run(args):
  assert len(args) in [0,1]
  if (len(args) == 0):
    n = 100
  else:
    n = int(args[0])
    assert n >= 0
  for use_dumps in [False, True]:
    exercise(n, use_dumps)
  print "OK"

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
